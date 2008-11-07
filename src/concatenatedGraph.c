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
#include "passageMarker.h"

void concatenateReadStarts(Node * target, Node * source, Graph * graph)
{
	IDnum sourceLength, targetLength;
	ShortReadMarker *sourceArray, *targetArray, *marker;
	IDnum index;
	Coordinate position, nodeLength;

	if (!readStartsAreActivated(graph))
		return;

	if (target == NULL || source == NULL)
		return;

	// Update Coordinates
	sourceArray = getNodeReads(source, graph);
	sourceLength = getNodeReadCount(source, graph);

	nodeLength = getNodeLength(target);
	for (index = 0; index < sourceLength; index++) {
		marker = getShortReadMarkerAtIndex(sourceArray, index);
		position = getShortReadMarkerPosition(marker);
		if (position != -1) {
			position += nodeLength;
			setShortReadMarkerPosition(marker, position);
		}
	}

	// Same but for symmetrical reads
	targetArray = getNodeReads(getTwinNode(target), graph);
	targetLength = getNodeReadCount(getTwinNode(target), graph);

	nodeLength = getNodeLength(source);
	for (index = 0; index < targetLength; index++) {
		marker = getShortReadMarkerAtIndex(targetArray, index);
		position = getShortReadMarkerPosition(marker);
		if (position != -1) {
			position += nodeLength;
			setShortReadMarkerPosition(marker, position);
		}
	}

	// Merging lists
	mergeNodeReads(target, source, graph);
	mergeNodeReads(getTwinNode(target), getTwinNode(source), graph);
}

// Replaces two consecutive nodes into a single equivalent node
// The extra memory is freed
void concatenateNodes(Node * nodeA, Node * nodeB, Graph * graph)
{
	PassageMarker *marker, *tmpMarker;
	Node *twinA = getTwinNode(nodeA);
	Node *twinB = getTwinNode(nodeB);
	Arc *arc;
	Category cat;

	// Arc management:
	// Freeing useless arcs
	while (getArc(nodeA) != NULL)
		destroyArc(getArc(nodeA), graph);

	// Correct arcs
	for (arc = getArc(nodeB); arc != NULL; arc = getNextArc(arc)) {
		if (getDestination(arc) != twinB)
			createAnalogousArc(nodeA, getDestination(arc),
					   arc, graph);
		else
			createAnalogousArc(nodeA, twinA, arc, graph);
	}

	// Passage marker management in node A:
	for (marker = getMarker(nodeA); marker != NULL;
	     marker = getNextInNode(marker))
		if (isTerminal(marker))
			incrementFinishOffset(marker,
					      getNodeLength(nodeB));

	// Swapping new born passageMarkers from B to A
	for (marker = getMarker(nodeB); marker != NULL; marker = tmpMarker) {
		tmpMarker = getNextInNode(marker);

		if (isInitial(marker)
		    || getNode(getPreviousInSequence(marker)) != nodeA) {
			extractPassageMarker(marker);
			transposePassageMarker(marker, nodeA);
			incrementFinishOffset(getTwinMarker(marker),
					      getNodeLength(nodeA));
		} else
			disconnectNextPassageMarker(getPreviousInSequence
						    (marker), graph);
	}

	// Read starts
	concatenateReadStarts(nodeA, nodeB, graph);

	// Gaps
	appendNodeGaps(nodeA, nodeB, graph);

	// Descriptor management (node)
	appendDescriptors(nodeA, nodeB);

	// Update uniqueness:
	setUniqueness(nodeA, getUniqueness(nodeA) || getUniqueness(nodeB));

	// Update virtual coverage
	for (cat = 0; cat < CATEGORIES; cat++)
		incrementVirtualCoverage(nodeA, cat,
					 getVirtualCoverage(nodeB, cat));

	// Update original virtual coverage
	for (cat = 0; cat < CATEGORIES; cat++)
		incrementOriginalVirtualCoverage(nodeA, cat,
						 getOriginalVirtualCoverage
						 (nodeB, cat));

	// Freeing gobbled node
	destroyNode(nodeB, graph);
}

// Detects sequences that could be simplified through concatentation
// Iterates till graph cannot be more simplified
// Useless nodes are freed from memory and remaining ones are renumbered
void concatenateGraph(Graph * graph)
{
	IDnum nodeIndex;
	Node *node, *twin;

	puts("Concatenation...");

	for (nodeIndex = 1; nodeIndex < nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);

		if (node == NULL)
			continue;

		twin = getTwinNode(node);
		while (simpleArcCount(node) == 1
		       &&
		       simpleArcCount(getTwinNode
				      (getDestination(getArc(node)))) ==
		       1) {
			if (getDestination(getArc(node)) == twin
			    || getDestination(getArc(node)) == node)
				break;
			concatenateNodes(node,
					 getDestination(getArc(node)),
					 graph);
		}

		while (simpleArcCount(twin) == 1
		       &&
		       simpleArcCount(getTwinNode
				      (getDestination(getArc(twin)))) ==
		       1) {
			if (getDestination(getArc(twin)) == node
			    || getDestination(getArc(twin)) == twin)
				break;
			concatenateNodes(twin,
					 getDestination(getArc(twin)),
					 graph);
		}
	}

	renumberNodes(graph);
	puts("Concatenation over!");
}
