/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "graph.h"
#include "recycleBin.h"
#include "passageMarker.h"
#include "graphStats.h"
#include "concatenatedGraph.h"
#include "readSet.h"

#define LONG_NODE_CUTOFF 0
#define LN2 0.693147
#define PROBABILITY_CUTOFF 5
#define MAX_READ_COUNT 100
#define MAX_READ_LENGTH 2000
#define MULTIPLICITY_CUTOFF 2

static Graph *graph = NULL;
static PassageMarker *path = NULL;
static RecycleBin *listMemory = NULL;
static boolean(*isUniqueFunction) (Node * node) = NULL;
static double expected_coverage = 1;
static TightString **sequences = NULL;

static IDnum multCounter = 0;
static IDnum dbgCounter = 0;
static IDnum nullCounter = 0;

typedef struct connection_st Connection;

struct connection_st {
	Node *node;
	IDnum multiplicity;
	PassageMarker *marker;
	Connection *next;
};

static RecycleBin *nodeListMemory = NULL;

#define BLOCKSIZE 1000

static Connection *allocateConnection()
{
	if (nodeListMemory == NULL)
		nodeListMemory =
		    newRecycleBin(sizeof(Connection), BLOCKSIZE);

	return allocatePointer(nodeListMemory);
}

static void deallocateConnection(Connection * nodeList)
{
	deallocatePointer(nodeListMemory, nodeList);
}

void setBaseCoverage(double coverage)
{
	expected_coverage = coverage;
}

boolean isUniqueBasic(Node * node)
{
	if (getNodeLength(node) < LONG_NODE_CUTOFF) {
		return false;
	}
	if (readCoverage(node) / (double) getNodeLength(node) >
	    1.5 * expected_coverage) {
		return false;
	}

	return true;
}

boolean isUniqueSolexa(Node * node)
{

	Coordinate nodeLength = getNodeLength(node);
	Coordinate nodeCoverage =
	    (getVirtualCoverage(node, 0) + getVirtualCoverage(node, 1));
	double nodeDensity, probability;

	if (nodeLength == 0) {
		return false;
	}
	if (nodeLength > LONG_NODE_CUTOFF) {
		nodeDensity = nodeCoverage / (double) nodeLength;

		probability =
		    LN2 / 2 +
		    nodeLength / (2 * expected_coverage) *
		    (expected_coverage * expected_coverage -
		     nodeDensity * nodeDensity / 2);
		return probability > PROBABILITY_CUTOFF;
	} else {
		return false;
		probability =
		    expected_coverage * nodeLength - nodeCoverage / LN2;
		return probability > 0;
	}
}

static void identifyUniqueNodes()
{
	IDnum index;
	Node *node;
	IDnum counter = 0;

	puts("Identifying unique nodes");

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (node == NULL)
			continue;

		setUniqueness(node, isUniqueFunction(node));

		if (getUniqueness(node))
			counter++;
	}

	printf("Done, %lu unique nodes counted\n", counter);
}

static void trimNodeLength(Node * node)
{
	PassageMarker *marker;
	Coordinate minOffset = getNodeLength(node);

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getFinishOffset(marker) < minOffset)
			minOffset = getFinishOffset(marker);

	if (minOffset == 0)
		return;

	clipNodeLength(node, 0, minOffset);

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		setFinishOffset(marker,
				getFinishOffset(marker) - minOffset);
}

static void simplifyNode(Node * node)
{
	Node *twin = getTwinNode(node);
	Node *destination, *twinDestination;

	if (arcCount(node) == 0)
		trimNodeLength(node);

	if (simpleArcCount(node) != 1)
		return;

	destination = getDestination(getArc(node));
	twinDestination = getTwinNode(destination);

	while (simpleArcCount(node) == 1
	       && simpleArcCount(twinDestination) == 1
	       && destination != twin && destination != node) {
		if (getNode(path) == destination
		    || getNode(path) == twinDestination) {
			return;
		}

		if (getNode(getNextInSequence(path)) == destination
		    || getNode(getNextInSequence(path)) ==
		    twinDestination) {
			return;
		}

		if (getNodeStatus(destination))
			return;

		concatenateNodes(node, destination, graph);

		if (simpleArcCount(node) != 1)
			return;

		destination = getDestination(getArc(node));
		twinDestination = getTwinNode(destination);
	}

}

static void decrementArc(Node * origin, Node * destination)
{
	Arc *arc = getArcBetweenNodes(origin, destination, graph);

	if (arc == NULL)
		return;

	changeMultiplicity(arc, -1);

	if (getMultiplicity(arc) == 0) {
		destroyArc(arc, graph);
	}
}

static boolean uniqueNodesConnect(Node * startingNode,
				  boolean(*isUnique) (Node *))
{
	Node *destination = NULL;
	PassageMarker *startMarker, *currentMarker;
	Connection *newList;
	Connection *list = NULL;
	boolean multipleHits = false;
	Coordinate debug;

	if (arcCount(startingNode) == 0)
		return false;

	if (getMarker(startingNode) == NULL)
		return false;

	dbgCounter++;

	// Checking for multiple destinations
	for (startMarker = getMarker(startingNode); startMarker != NULL;
	     startMarker = getNextInNode(startMarker)) {
		debug = 0;
		for (currentMarker = getNextInSequence(startMarker);
		     currentMarker != NULL;
		     currentMarker = getNextInSequence(currentMarker)) {
			if (!getUniqueness(getNode(currentMarker))) {
				debug +=
				    getPassageMarkerLength(currentMarker);
				continue;
			} else if (getNodeStatus(getNode(currentMarker))) {
				for (newList = list; newList != NULL;
				     newList = newList->next) {
					if (newList->node ==
					    getNode(currentMarker)) {
						newList->multiplicity++;
						break;
					}
				}
				if (newList == NULL)
					abort();
				break;
			} else {
				setSingleNodeStatus(getNode(currentMarker),
						    true);
				newList = allocateConnection();
				newList->node = getNode(currentMarker);
				newList->multiplicity = 1;
				newList->marker = startMarker;
				newList->next = list;
				list = newList;
				break;
			}
		}
	}

	while (list != NULL) {
		newList = list;
		list = newList->next;
		setSingleNodeStatus(newList->node, false);
		if (newList->multiplicity >= MULTIPLICITY_CUTOFF) {
			if (destination == NULL) {
				destination = newList->node;
				path = newList->marker;
			} else if (destination != newList->node)
				multipleHits = true;
		}
		deallocateConnection(newList);
	}

	// Aligning long reads to each other:
	// TODO 

	// Merge pairwise alignments and produce consensus
	// TODO

	if (multipleHits) {
		multCounter++;
		return false;
	}
	//DEBUG
	if (!
	    (destination != NULL && destination != startingNode
	     && destination != getTwinNode(startingNode)))
		nullCounter++;

	return (destination != NULL && destination != startingNode
		&& destination != getTwinNode(startingNode));
}

static boolean goesToNode(PassageMarker * marker, Node * node)
{
	PassageMarker *current;

	for (current = marker; current != NULL;
	     current = getNextInSequence(current))
		if (getNode(current) == node)
			return true;

	return false;
}

static void updateMembers(Node * bypass, Node * nextNode)
{
	PassageMarker *marker, *next, *tmp;
	Node *nextNextNode;
	Coordinate nextLength = getNodeLength(nextNode);

	// Remove unwanted arcs
	while (getArc(bypass) != NULL)
		destroyArc(getArc(bypass), graph);

	// Update  marker + arc info
	for (marker = getMarker(bypass); marker != NULL; marker = tmp) {
		tmp = getNextInNode(marker);

		if (!isTerminal(marker)
		    && getNode(getNextInSequence(marker)) == nextNode) {
			// Marker steps right into target
			next = getNextInSequence(marker);
			nextNextNode = getNode(getNextInSequence(next));

			createArc(bypass, nextNextNode, graph);
			decrementArc(nextNode, nextNextNode);

			disconnectNextPassageMarker(marker, graph);
			destroyPassageMarker(next);
		} else if (getUniqueness(nextNode)
			   && goesToNode(marker, nextNode)) {
			// Marker goes indirectly to target
			while (getNode(getNextInSequence(marker)) !=
			       nextNode) {
				next = getNextInSequence(marker);
				setPreviousInSequence(marker,
						      getNextInSequence
						      (next));
				disconnectNextPassageMarker(marker, graph);
				destroyPassageMarker(next);
			}

		} else if (!isTerminal(marker)
			   && getFinishOffset(marker) == 0) {
			// Marker goes somewhere else than to target
			next = getNextInSequence(marker);
			decrementArc(getNode(marker), getNode(next));

			incrementFinishOffset(marker, nextLength);
		} else {
			// Marker goes nowhere
			incrementFinishOffset(marker, nextLength);
		}
	}
}

static void admitGroupies(Node * source, Node * bypass)
{
	PassageMarker *marker, *tmpMarker;

	for (marker = getMarker(source); marker != NULL;
	     marker = tmpMarker) {
		tmpMarker = getNextInNode(marker);
		extractPassageMarker(marker);
		transposePassageMarker(marker, bypass);
		incrementFinishOffset(getTwinMarker(marker),
				      getNodeLength(source));
	}

}

static void adjustShortReads(Node * target, PassageMarker * pathMarker)
{
	ShortReadMarker *targetArray, *marker;
	IDnum targetLength, index;
	Coordinate position, nodeLength;

	if (!readStartsAreActivated(graph))
		return;

	targetArray = getNodeReads(getTwinNode(target), graph);
	targetLength = getNodeReadCount(getTwinNode(target), graph);

	nodeLength = getPassageMarkerLength(pathMarker);

	for (index = 0; index < targetLength; index++) {
		marker = getShortReadMarkerAtIndex(targetArray, index);
		position = getShortReadMarkerPosition(marker);
		position += nodeLength;
		setShortReadMarkerPosition(marker, position);
	}
}

static Node *bypass()
{
	Node *bypass = getNode(path);
	Node *next = NULL;
	Coordinate oldLength;
	PassageMarker *marker;

	for (marker = path;
	     marker != NULL && !getUniqueness(getNode(marker));
	     marker = getNextInSequence(marker))
		setNodeStatus(getNode(marker), true);

	if (marker != NULL && getUniqueness(getNode(marker)))
		setNodeStatus(getNode(marker), true);

	// Update extensive variables (length + descriptors + passage markers)
	while (!isTerminal(path)) {
		next = getNode(getNextInSequence(path));

		if (next == NULL) {
			setNodeStatus(bypass, false);
			return bypass;
		}

		setNodeStatus(next, false);

		// Overall node update 
		oldLength = getNodeLength(bypass);
		if (!getUniqueness(next)) {
			adjustShortReads(bypass, getNextInSequence(path));
			appendSequence(bypass, sequences,
				       getNextInSequence(path), graph);
		} else {
			concatenateReadStarts(bypass, next, graph);
			appendDescriptors(bypass, next);
		}

		if (!getUniqueness(bypass) && getUniqueness(next))
			setUniqueness(bypass, true);

		// Members
		updateMembers(bypass, next);

		// Clean data
		if (!isTerminal(path) && !getUniqueness(next)) {
			if (getMarker(next) == NULL) {
				destroyNode(next, graph);
			} else {
				simplifyNode(next);
				simplifyNode(getTwinNode(next));
			}
		} else {
			break;
		}
	}

	// Pathological cases
	if (next == bypass || next == getTwinNode(bypass)) {
		abort();
	}
	// Remove unique groupies from original path 
	admitGroupies(next, bypass);

	destroyNode(next, graph);

	// Update
	//if (!getNodeStatus(bypass) && isUniqueFunction(bypass))
	//      setNodeStatus(bypass, true);
	//      setNodeStatus(bypass, isUniqueFunction(bypass));

	simplifyNode(bypass);
	simplifyNode(getTwinNode(bypass));
	setNodeStatus(bypass, false);

	if (!getUniqueness(bypass) && isUniqueFunction(bypass))
		setUniqueness(bypass, true);

	return bypass;
}

static void trimLongReadTips()
{
	IDnum index;
	Node *node;
	PassageMarker *marker, *next;

	printf("Trimming read tips\n");

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getUniqueness(node))
			continue;

		for (marker = getMarker(node); marker != NULL;
		     marker = next) {
			next = getNextInNode(marker);

			if (!isInitial(marker) && !isTerminal(marker))
				continue;

			if (isTerminal(marker))
				marker = getTwinMarker(marker);

			while (!getUniqueness(getNode(marker))) {
				if (next != NULL
				    && (marker == next
					|| marker == getTwinMarker(next)))
					next = getNextInNode(next);
				if (getNextInSequence(marker) != NULL) {
					marker = getNextInSequence(marker);
					destroyPassageMarker
					    (getPreviousInSequence
					     (marker));
				} else {
					destroyPassageMarker(marker);
					break;
				}
			}
		}
	}
}

void readCoherentGraph(Graph * inGraph, boolean(*isUnique) (Node *),
		       double coverage, ReadSet * reads)
{
	IDnum nodeIndex;
	Node *node;
	IDnum previousNodeCount = 0;

	graph = inGraph;
	isUniqueFunction = isUnique;
	listMemory = newRecycleBin(sizeof(PassageMarkerList), 100000);
	expected_coverage = coverage;
	sequences = reads->tSequences;

	puts("Read coherency...");
	//reassessArcMultiplicities(graph);
	resetNodeStatus(graph);

	identifyUniqueNodes();
	trimLongReadTips();

	previousNodeCount = 0;
	while (previousNodeCount != nodeCount(graph)) {

		previousNodeCount = nodeCount(graph);

		for (nodeIndex = 1; nodeIndex <= nodeCount(graph);
		     nodeIndex++) {

			node = getNodeInGraph(graph, nodeIndex);

			if (node == NULL || !getUniqueness(node))
				continue;

			while (uniqueNodesConnect(node, isUnique))
				node = bypass();

			node = getTwinNode(node);

			while (uniqueNodesConnect(node, isUnique))
				node = bypass();

		}

		concatenateGraph(graph);
		break;
	}

	destroyRecycleBin(listMemory);

	printf("Confronted to %li multiple hits and %li null over %li\n",
	       multCounter, nullCounter, dbgCounter);

	puts("Read coherency over!");
}
