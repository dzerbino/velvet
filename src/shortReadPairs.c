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
#include <time.h>

#include "globals.h"
#include "graph.h"
#include "concatenatedGraph.h"
#include "recycleBin.h"
#include "locallyCorrectedGraph.h"
#include "passageMarker.h"
#include "readSet.h"
#include "math.h"

#define BLOCK_SIZE  100000
#define LN2 1.4
#define BACKTRACK_CUTOFF 100

typedef struct connection_st Connection;
typedef struct miniConnection_st MiniConnection;
typedef struct readOccurence_st ReadOccurence;

struct connection_st {
	Node *destination;
	IDnum direct_count;
	IDnum paired_count;
	Coordinate distance;
	double variance;
	Connection *next;
	Connection *previous;
	Connection *twin;
};

struct miniConnection_st {
	Coordinate distance;
	double variance;
	Connection *frontReference;
	Connection *backReference;
	NodeList *nodeList;
};

struct readOccurence_st {
	IDnum nodeID;
	Coordinate position;
	Coordinate offset;
};

// Global params
static Graph *graph;
static ReadSet * reads;
static RecycleBin *connectionMemory = NULL;
static RecycleBin *nodeListMemory = NULL;
static double expected_coverage;
static IDnum UNRELIABLE_CONNECTION_CUTOFF = 10;

// Global pointers
static NodeList *markedNodes;
static Connection **scaffold = NULL;
static MiniConnection *localScaffold = NULL;

static Connection *allocateConnection()
{
	if (connectionMemory == NULL)
		connectionMemory =
		    newRecycleBin(sizeof(Connection), BLOCK_SIZE);

	return allocatePointer(connectionMemory);
}

static void deallocateConnection(Connection * connect)
{
	deallocatePointer(connectionMemory, connect);
}

static NodeList *allocateNodeList()
{
	if (nodeListMemory == NULL)
		nodeListMemory =
		    newRecycleBin(sizeof(NodeList), BLOCK_SIZE);

	return allocatePointer(nodeListMemory);
}

static void deallocateNodeList(NodeList * nodeList)
{
	deallocatePointer(nodeListMemory, nodeList);
}

static NodeList *recordNode(Node * node)
{
	NodeList *nodeList = allocateNodeList();
	nodeList->node = node;
	nodeList->next = markedNodes;
	nodeList->previous = NULL;

	if (markedNodes != NULL)
		markedNodes->previous = nodeList;

	markedNodes = nodeList;

	return nodeList;
}

static void destroyNodeList(NodeList * nodeList)
{
	//printf("Destroy NL  %p > %p > %p\n", nodeList->previous, nodeList, nodeList->next);

	if (nodeList->previous != NULL)
		nodeList->previous->next = nodeList->next;
	else
		markedNodes = nodeList->next;

	if (nodeList->next != NULL)
		nodeList->next->previous = nodeList->previous;

	nodeList->previous = nodeList->next = NULL;

	deallocateNodeList(nodeList);
}

static Node *popNodeRecord()
{
	MiniConnection *localConnect;

	NodeList *nodeList = markedNodes;
	Node *node;

	if (markedNodes == NULL)
		return NULL;

	node = nodeList->node;
	markedNodes = nodeList->next;
	if (markedNodes != NULL)
		markedNodes->previous = NULL;

	localConnect =
	    &localScaffold[getNodeID(nodeList->node) + nodeCount(graph)];
	localConnect->nodeList = NULL;

	deallocateNodeList(nodeList);
	return node;
}

static double norm(double X) {
	return 0.4 * exp(-X*X / 2);
}

static double normInt(double X, double Y) {
	return (erf(0.7 * Y) - erf(0.7 * X))/2;
}

static IDnum expectedNumberOfConnections(IDnum IDA, Connection* connect, IDnum ** counts, Category cat) {
	Node *A = getNodeInGraph(graph, IDA);
	Node* B = connect->destination;
	double left, middle, right;
	Coordinate longLength, shortLength, D;
	IDnum longCount;
	double M,N,O,P;
	Coordinate mu = getInsertLength(graph, cat);
	double sigma = sqrt(getInsertLength_var(graph, cat));	
	double result;

	if (mu <= 0)
		return 0;

	if (getNodeLength(A) < getNodeLength(B)) {
		longLength = getNodeLength(B);
		shortLength = getNodeLength(A);
		longCount = counts[cat][getNodeID(B) + nodeCount(graph)];
	} else {
		longLength = getNodeLength(A);
		shortLength = getNodeLength(B);
		longCount = counts[cat][IDA + nodeCount(graph)];
	}

	D = connect->distance - (longLength + shortLength)/2;

	M = (D - mu) / sigma;
	N = (D + shortLength - mu) / sigma;
	O = (D + longLength - mu) / sigma;
	P = (D + shortLength + longLength - mu) / sigma;	

	left = ((norm(M) - norm(N)) - M * normInt(M,N)) * sigma; 
	middle = shortLength * normInt(N,O); 
	right = ((norm(O) - norm(P)) - P * normInt(O,P)) * (-sigma);

	result = (longCount * (left + middle + right)) / longLength;

	if (result > 0)
		return result;
	else 
		return 0;
}

static void destroyConnection(Connection * connect, IDnum nodeID)
{
	Connection *previous, *next;

	//printf("Destroying connection from %li to %li\n", nodeID, getNodeID(connect->destination));

	if (connect == NULL)
		return;

	previous = connect->previous;
	next = connect->next;

	if (previous != NULL)
		previous->next = next;
	if (next != NULL)
		next->previous = previous;

	if (scaffold[nodeID + nodeCount(graph)] == connect)
		scaffold[nodeID + nodeCount(graph)] = next;

	if (connect->twin != NULL) {
		connect->twin->twin = NULL;
		destroyConnection(connect->twin,
				  getNodeID(connect->destination));
	}

	deallocateConnection(connect);
}

static boolean testConnection(IDnum IDA, Connection* connect, IDnum ** counts) {
	IDnum total = 0;
	Category cat;

	// Spare unique -> undetermined node connections
	if (!getUniqueness(connect->destination))
		return true;

	// Destroy tenuous connections
	if (connect->paired_count + connect->direct_count < UNRELIABLE_CONNECTION_CUTOFF)
		return false; 

	for (cat = 0; cat <= CATEGORIES; cat++)
		total += expectedNumberOfConnections(IDA, connect, counts, cat);

	// Remove inconsistent connections
	return connect->paired_count > total / 10; 
}

void detachImprobablePairs(ReadSet * sequences)
{
	IDnum index, nodeIndex;
	IDnum maxNodeIndex = 2 * nodeCount(graph) + 1;
	ShortReadMarker *nodeArray, *shortMarker;
	Node *node;
	IDnum nodeReadCount;
	IDnum seqID, pairID;
	IDnum *mateReads = sequences->mateReads;
	Category *cats = sequences->categories;

	for (nodeIndex = 0; nodeIndex < maxNodeIndex; nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex - nodeCount(graph));
		if (node == NULL)
			continue;

		nodeArray = getNodeReads(node, graph);
		nodeReadCount = getNodeReadCount(node, graph);

		for (index = 0; index < nodeReadCount; index++) {
			shortMarker =
			    getShortReadMarkerAtIndex(nodeArray, index);

			seqID = getShortReadMarkerID(shortMarker);
			if (mateReads[seqID] == -1)
				continue;

			if (getNodeLength(node) -
			    getShortReadMarkerPosition(shortMarker) >
			    2 * getInsertLength(graph, cats[seqID])) {
				pairID = mateReads[seqID];

				if (pairID != -1) {
					mateReads[seqID] = -1;
					mateReads[pairID] = -1;
				}
			}
		}
	}
}

static IDnum *computeReadToNodeCounts()
{
	IDnum readIndex, nodeIndex;
	IDnum maxNodeIndex = 2 * nodeCount(graph) + 1;
	IDnum maxReadIndex = sequenceCount(graph) + 1;
	IDnum *readNodeCounts = calloc(maxReadIndex, sizeof(IDnum));
	boolean *readMarker = calloc(maxReadIndex, sizeof(boolean));
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarker *marker;
	Node *node;
	IDnum nodeReadCount;

	if (readNodeCounts == NULL || readMarker == NULL) {
		puts("Calloc failure");
		exit(1);
	}

	puts("Computing read to node mapping array sizes");

	for (nodeIndex = 0; nodeIndex < maxNodeIndex; nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex - nodeCount(graph));
		if (node == NULL)
			continue;
		nodeArray = getNodeReads(node, graph);
		nodeReadCount = getNodeReadCount(node, graph);

		// Short reads
		for (readIndex = 0; readIndex < nodeReadCount; readIndex++) {
			shortMarker =
			    getShortReadMarkerAtIndex(nodeArray,
						      readIndex);
			readNodeCounts[getShortReadMarkerID
				       (shortMarker)]++;
		}

		// Long reads
		for (marker = getMarker(node); marker != NULL;
		     marker = getNextInNode(marker)) {
			readIndex = getPassageMarkerSequenceID(marker);
			if (readIndex < 0)
				continue;

			if (readMarker[readIndex])
				continue;

			readNodeCounts[readIndex]++;
			readMarker[readIndex] = true;
		}

		// Clean up marker array
		for (marker = getMarker(node); marker != NULL;
		     marker = getNextInNode(marker)) {
			readIndex = getPassageMarkerSequenceID(marker);
			if (readIndex > 0)
				readMarker[readIndex] = false;
		}
	}

	free(readMarker);
	return readNodeCounts;
}

static ReadOccurence **allocateReadToNodeTables(IDnum * readNodeCounts)
{
	IDnum readIndex;
	IDnum maxReadIndex = sequenceCount(graph) + 1;
	ReadOccurence **readNodes =
	    calloc(maxReadIndex, sizeof(ReadOccurence *));

	if (readNodes == NULL) {
		puts("Calloc failure");
		exit(1);
	}

	for (readIndex = 1; readIndex < maxReadIndex; readIndex++) {
		if (readNodeCounts[readIndex] != 0) {
			readNodes[readIndex] =
			    calloc(readNodeCounts[readIndex],
				   sizeof(ReadOccurence));
			if (readNodes[readIndex] == NULL
			    && readNodeCounts[readIndex] > 0) {
				puts("Calloc failure");
				exit(1);
			}
			readNodeCounts[readIndex] = 0;
		}
	}

	return readNodes;
}

static void computePartialReadToNodeMapping(IDnum nodeID,
					    ReadOccurence ** readNodes,
					    IDnum * readNodeCounts,
					    boolean * readMarker)
{
	ShortReadMarker *shortMarker;
	IDnum index, readIndex;
	ReadOccurence *readArray, *readOccurence;
	Node *node = getNodeInGraph(graph, nodeID);
	ShortReadMarker *nodeArray = getNodeReads(node, graph);
	IDnum nodeReadCount = getNodeReadCount(node, graph);
	PassageMarker *marker;

	for (index = 0; index < nodeReadCount; index++) {
		shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
		readIndex = getShortReadMarkerID(shortMarker);
		readArray = readNodes[readIndex];
		readOccurence = &readArray[readNodeCounts[readIndex]];
		readOccurence->nodeID = nodeID;
		readOccurence->position =
		    getShortReadMarkerPosition(shortMarker);
		readOccurence->offset =
		    getShortReadMarkerOffset(shortMarker);
		readNodeCounts[readIndex]++;
	}

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		readIndex = getPassageMarkerSequenceID(marker);
		if (readIndex < 0)
			continue;

		if (!readMarker[readIndex]) {
			readArray = readNodes[readIndex];
			readOccurence =
			    &readArray[readNodeCounts[readIndex]];
			readOccurence->nodeID = nodeID;
			readOccurence->position = getStartOffset(marker);
			readOccurence->offset =
			    getPassageMarkerStart(marker);
			readNodeCounts[readIndex]++;
			readMarker[readIndex] = true;
		} else {
			readArray = readNodes[readIndex];
			readOccurence =
			    &readArray[readNodeCounts[readIndex] - 1];
			readOccurence->position = -1;
			readOccurence->offset = -1;
		}
	}

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		readIndex = getPassageMarkerSequenceID(marker);
		if (readIndex > 0)
			readMarker[readIndex] = false;
	}
}

static ReadOccurence **computeReadToNodeMappings(IDnum * readNodeCounts)
{
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);
	ReadOccurence **readNodes =
	    allocateReadToNodeTables(readNodeCounts);
	boolean *readMarker =
	    calloc(sequenceCount(graph) + 1, sizeof(boolean));

	if (readMarker == NULL) {
		puts("Calloc failure");
		exit(1);
	}

	puts("Computing read to node mappings");

	for (nodeID = -nodes; nodeID <= nodes; nodeID++)
		if (nodeID != 0 && getNodeInGraph(graph, nodeID))
			computePartialReadToNodeMapping(nodeID, readNodes,
							readNodeCounts,
							readMarker);

	free(readMarker);
	return readNodes;
}

static Connection *findConnection(IDnum nodeID, IDnum node2ID)
{
	Node *node2 = getNodeInGraph(graph, node2ID);
	Connection *connect;

	if (node2 == NULL)
		return NULL;

	for (connect = scaffold[nodeID + nodeCount(graph)];
	     connect != NULL; connect = connect->next)
		if (connect->destination == node2)
			break;

	return connect;
}

static void createTwinConnection(IDnum nodeID, IDnum node2ID,
				 Connection * connect)
{
	Connection *newConnection = allocateConnection();
	IDnum nodeIndex = nodeID + nodeCount(graph);

	// Fill in
	newConnection->distance = connect->distance;
	newConnection->variance = connect->variance;
	newConnection->direct_count = connect->direct_count;
	newConnection->paired_count = connect->paired_count;
	newConnection->destination = getNodeInGraph(graph, node2ID);

	// Batch to twin
	newConnection->twin = connect;
	connect->twin = newConnection;

	// Insert in scaffold
	newConnection->previous = NULL;
	newConnection->next = scaffold[nodeIndex];
	if (scaffold[nodeIndex] != NULL)
		scaffold[nodeIndex]->previous = newConnection;
	scaffold[nodeIndex] = newConnection;
}

static Connection *createNewConnection(IDnum nodeID, IDnum node2ID,
				       IDnum direct_count,
				       IDnum paired_count,
				       Coordinate distance,
				       double variance)
{
	Node *destination = getNodeInGraph(graph, node2ID);
	IDnum nodeIndex = nodeID + nodeCount(graph);
	Connection *connect = allocateConnection();

	// Fill in 
	connect->destination = destination;
	connect->direct_count = direct_count;
	connect->paired_count = paired_count;
	connect->distance = distance;
	connect->variance = variance;

	// Insert in scaffold
	connect->previous = NULL;
	connect->next = scaffold[nodeIndex];
	if (scaffold[nodeIndex] != NULL)
		scaffold[nodeIndex]->previous = connect;
	scaffold[nodeIndex] = connect;

	// Event. pair up to twin
	if (getUniqueness(destination))
		createTwinConnection(node2ID, nodeID, connect);
	else
		connect->twin = NULL;

	return connect;
}

static void readjustConnection(Connection * connect, Coordinate distance,
			       double variance, IDnum direct_count, IDnum paired_count)
{
	connect->direct_count += direct_count;
	connect->paired_count += paired_count;
	connect->distance =
	    (variance * connect->distance +
	     distance * connect->variance) / (variance +
					      connect->variance);
	connect->variance =
	    (variance *
	     connect->variance) / (variance + connect->variance);

	if (connect->twin != NULL) {
		connect->twin->distance = connect->distance;
		connect->twin->variance = connect->variance;
		connect->twin->direct_count = connect->direct_count;
		connect->twin->paired_count = connect->paired_count;
	}
}

static void createConnection(IDnum nodeID, IDnum node2ID,
			     IDnum direct_count,
			     IDnum paired_count,
			     Coordinate distance, double variance)
{
	Connection *connect = findConnection(nodeID, node2ID);

	if (connect != NULL)
		readjustConnection(connect, distance, variance, direct_count, paired_count);
	else
		createNewConnection(nodeID, node2ID, direct_count, paired_count,
				    distance, variance);
}

static void projectFromSingleRead(Node * node,
				  ReadOccurence * readOccurence,
				  Coordinate position,
				  Coordinate offset, Coordinate length)
{
	Coordinate distance = 0;
	Node *target = getNodeInGraph(graph, -readOccurence->nodeID);
	double variance = 1;

	if (target == getTwinNode(node) || target == node)
		return;

	if (position < 0) {
		variance += getNodeLength(node) * getNodeLength(node) / 16;
		// distance += 0;
	} else {
		// variance += 0;
		distance += position - offset - getNodeLength(node) / 2;
	}

	if (readOccurence->position < 0) {
		variance +=
		    getNodeLength(target) * getNodeLength(target) / 16;
		//distance += 0;
	} else {
		// variance += 0;
		distance +=
		    -readOccurence->position + readOccurence->offset +
		    getNodeLength(target) / 2;
	}

	if (position < 0 || readOccurence->position < 0) {
		variance += length * length / 16;
		createConnection(getNodeID(node), getNodeID(target), 1, 0,
				 distance, variance);
		createConnection(-getNodeID(node), -getNodeID(target), 1, 0,
				 -distance, variance);
	} else if (distance > 0) {
		createConnection(getNodeID(node), getNodeID(target), 1, 0,
				 distance, variance);
	} else {
		createConnection(-getNodeID(node), -getNodeID(target), 1, 0,
				 -distance, variance);
	}
}

static void projectFromReadPair(Node * node, ReadOccurence * readOccurence,
				Coordinate position, Coordinate offset,
				Coordinate insertLength,
				double insertVariance)
{
	Coordinate distance = insertLength;
	Coordinate variance = insertVariance;
	Node *target = getNodeInGraph(graph, readOccurence->nodeID);

	if (target == getTwinNode(node) || target == node)
		return;

	if (getUniqueness(target) && getNodeID(target) < getNodeID(node))
		return;

	if (position < 0) {
		variance += getNodeLength(node) * getNodeLength(node) / 16;
		// distance += 0;
	} else {
		// variance += 0;
		distance += position - offset - getNodeLength(node) / 2;
	}

	if (readOccurence->position < 0) {
		variance +=
		    getNodeLength(target) * getNodeLength(target) / 16;
		//distance += 0;
	} else {
		// variance += 0;
		distance +=
		    readOccurence->position - readOccurence->offset -
		    getNodeLength(target) / 2;
	}

	createConnection(getNodeID(node), getNodeID(target), 0, 1, distance,
			 variance);
}

static void projectFromShortRead(Node * node,
				 ShortReadMarker * shortMarker,
				 IDnum * readPairs, Category * cats,
				 ReadOccurence ** readNodes,
				 IDnum * readNodeCounts,
				 Coordinate * lengths)
{
	IDnum index;
	IDnum readIndex = getShortReadMarkerID(shortMarker);
	ReadOccurence *readArray;
	IDnum readPairIndex;
	Category cat;
	Coordinate position = getShortReadMarkerPosition(shortMarker);
	Coordinate offset = getShortReadMarkerOffset(shortMarker);
	Coordinate length = lengths[getShortReadMarkerID(shortMarker) - 1];
	Coordinate insertLength;
	double insertVariance;

	// Going through single-read information
	if (readNodeCounts[readIndex] > 1 && position > 0) {
		readArray = readNodes[readIndex];
		for (index = 0; index < readNodeCounts[readIndex]; index++)
			projectFromSingleRead(node, &readArray[index],
					      position, offset, length);
	}
	// Going through paired read information
	if (readPairs == NULL)
		return;

	readPairIndex = readPairs[readIndex - 1] + 1;

	if (readPairIndex == 0)
		return;

	cat = cats[readIndex - 1];
	insertLength = getInsertLength(graph, cat);
	insertVariance = getInsertLength_var(graph, cat);

	readArray = readNodes[readPairIndex];
	for (index = 0; index < readNodeCounts[readPairIndex]; index++)
		projectFromReadPair(node, &readArray[index], position,
				    offset, insertLength, insertVariance);

}

static void projectFromLongRead(Node * node, PassageMarker * marker,
				IDnum * readPairs, Category * cats,
				ReadOccurence ** readNodes,
				IDnum * readNodeCounts,
				Coordinate * lengths)
{
	IDnum index;
	IDnum readIndex = getPassageMarkerSequenceID(marker);
	ReadOccurence *readArray;
	IDnum readPairIndex;
	Category cat;
	Coordinate position = getStartOffset(marker);
	Coordinate offset = getPassageMarkerStart(marker);
	Coordinate length =
	    lengths[getPassageMarkerSequenceID(marker) - 1];
	Coordinate insertLength;
	double insertVariance;

	// Going through single-read information
	if (readNodeCounts[readIndex] > 1 && position > 0) {
		readArray = readNodes[readIndex];
		for (index = 0; index < readNodeCounts[readIndex]; index++)
			projectFromSingleRead(node, &readArray[index],
					      position, offset, length);
	}
	// Going through paired read information
	if (readPairs == NULL)
		return;

	readPairIndex = readPairs[readIndex - 1] + 1;

	if (readPairIndex == 0)
		return;

	cat = cats[readIndex - 1];
	insertLength = getInsertLength(graph, cat);
	insertVariance = getInsertLength_var(graph, cat);

	readArray = readNodes[readPairIndex];
	for (index = 0; index < readNodeCounts[readPairIndex]; index++)
		projectFromReadPair(node, &readArray[index], position,
				    offset, insertLength, insertVariance);

}

static void projectFromNode(IDnum nodeID,
			    ReadOccurence ** readNodes,
			    IDnum * readNodeCounts,
			    IDnum * readPairs, Category * cats,
			    boolean * dubious, Coordinate * lengths)
{
	IDnum index;
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarker *marker;
	Node *node;
	IDnum nodeReadCount;

	node = getNodeInGraph(graph, nodeID);

	if (node == NULL || !getUniqueness(node))
		return;

	nodeArray = getNodeReads(node, graph);
	nodeReadCount = getNodeReadCount(node, graph);
	for (index = 0; index < nodeReadCount; index++) {
		shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
		if (dubious[getShortReadMarkerID(shortMarker) - 1])
			continue;
		projectFromShortRead(node, shortMarker, readPairs, cats,
				     readNodes, readNodeCounts, lengths);
	}

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		if (getPassageMarkerSequenceID(marker) > 0)
			projectFromLongRead(node, marker, readPairs, cats,
					    readNodes, readNodeCounts,
					    lengths);
	}
}

static Connection **computeNodeToNodeMappings(ReadOccurence ** readNodes,
					      IDnum * readNodeCounts,
					      IDnum * readPairs,
					      Category * cats,
					      boolean * dubious,
					      Coordinate * lengths)
{
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);
	scaffold = calloc(2 * nodes + 1, sizeof(Connection *));

	if (scaffold == NULL) {
		puts("Calloc failure");
		exit(1);
	}

	puts("Computing direct node to node mappings");

	for (nodeID = -nodes; nodeID <= nodes; nodeID++) {
		if (nodeID % 10000 == 0)
			printf("Scaffolding node %li\n", nodeID);

		projectFromNode(nodeID, readNodes, readNodeCounts,
				readPairs, cats, dubious, lengths);
	}

	return scaffold;
}

static void resetMiniConnection(Node * node, MiniConnection * localConnect,
				Coordinate distance, double variance,
				Connection * frontReference,
				Connection * backReference, boolean status)
{
	setSingleNodeStatus(node, status);
	localConnect->distance = distance;
	localConnect->variance = variance;
	localConnect->frontReference = frontReference;
	localConnect->backReference = backReference;
	localConnect->nodeList = recordNode(node);
}

static void setEmptyMiniConnection(Node * node)
{
	MiniConnection *localConnect =
	    &localScaffold[getNodeID(node) + nodeCount(graph)];
	localConnect->distance = 0;
	localConnect->variance = 1;
	localConnect->frontReference = NULL;
	localConnect->backReference = NULL;
	localConnect->nodeList = recordNode(node);
	setSingleNodeStatus(node, true);
}

static void readjustMiniConnection(Node * node,
				   MiniConnection * localConnect,
				   Coordinate distance,
				   Coordinate min_distance,
				   double variance,
				   Connection * frontReference,
				   Connection * backReference)
{

	localConnect->distance =
	    (variance * localConnect->distance +
	     distance * localConnect->variance) / (variance +
						   localConnect->variance);
	localConnect->variance =
	    (variance *
	     localConnect->variance) / (variance + localConnect->variance);

	if (frontReference != NULL)
		localConnect->frontReference = frontReference;
	if (backReference != NULL)
		localConnect->backReference = backReference;

	if (localConnect->distance > min_distance)
		setSingleNodeStatus(node, 1);
	else
		setSingleNodeStatus(node, -1);
}

static void integrateDerivativeDistances(Connection * connect,
					 Coordinate min_distance,
					 boolean direction)
{
	Node *reference = connect->destination;
	Node *destination;
	IDnum destinationID;
	Coordinate distance, baseDistance;
	double variance, baseVariance;
	Connection *connect2;
	MiniConnection *localConnect;

	// debug 
	IDnum counter = 0;

	if (!getUniqueness(reference))
		return;

	//printf("Opposite node %li length %li at %li ± %f\n", getNodeID(reference), getNodeLength(reference), connect->distance, connect->variance);

	baseDistance = connect->distance;
	baseVariance = connect->variance;

	for (connect2 =
	     scaffold[getNodeID(reference) + nodeCount(graph)];
	     connect2 != NULL; connect2 = connect2->next) {
		// Avoid null derivative
		if (connect2 == connect->twin)
			continue;

		destination = connect2->destination;

		// Beware of directionality
		if (!direction)
			destination = getTwinNode(destination);

		// Derivate values
		destinationID = getNodeID(destination);
		// Beware of directionality (bis)
		if (direction)
			distance = baseDistance - connect2->distance;
		else
			distance = connect2->distance - baseDistance;
		variance = connect2->variance + baseVariance;
		localConnect =
		    &localScaffold[destinationID + nodeCount(graph)];

		// Avoid over-projection
		if (distance < min_distance) {
			//printf("Node %li not at distance %li± %f (min %li)\n", destinationID, distance, variance, min_distance);
			continue;
		}

		counter++;

		if (getNodeStatus(destination)) {
			readjustMiniConnection(destination, localConnect,
					       distance, min_distance,
					       variance, NULL, NULL);
		} else
			resetMiniConnection(destination, localConnect,
					    distance, variance, NULL, NULL,
					    true);

		//printf("Node %li now at distance %li\n", destinationID, localConnect->distance);
	}

	//printf("%li secondary distances added\n", counter);
}

static void markInterestingNodes(Node * node)
{
	IDnum nodeID = getNodeID(node);
	IDnum nodeIndex = nodeID + nodeCount(graph);
	IDnum twinNodeIndex = -nodeID + nodeCount(graph);
	Connection *connect;
	Node *destination;
	MiniConnection *localConnect;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;

	// Mark own node
	setEmptyMiniConnection(node);

	// Loop thru primary scaffold
	for (connect = scaffold[nodeIndex]; connect != NULL;
	     connect = connect->next) {
		destination = getTwinNode(connect->destination);

		// DEBUG
		if (destination == node
		    || destination == getTwinNode(node))
			abort();

		localConnect =
		    &localScaffold[getNodeID(destination) +
				   nodeCount(graph)];

		if (getNodeStatus(destination)) {
			readjustMiniConnection(destination, localConnect,
					       connect->distance,
					       min_distance,
					       connect->variance, connect,
					       NULL);
			localConnect->backReference = NULL;
		} else {
			resetMiniConnection(destination, localConnect,
					    connect->distance,
					    connect->variance, connect,
					    NULL, true);
		}

		integrateDerivativeDistances(connect, min_distance, true);
	}

	// Loop thru twin's primary scaffold
	for (connect = scaffold[twinNodeIndex]; connect != NULL;
	     connect = connect->next) {
		destination = connect->destination;
		localConnect =
		    &localScaffold[getNodeID(destination) +
				   nodeCount(graph)];

		if (getNodeStatus(destination))
			readjustMiniConnection(destination, localConnect,
					       -connect->distance,
					       min_distance,
					       connect->variance, NULL,
					       connect);
		else
			resetMiniConnection(destination, localConnect,
					    -connect->distance,
					    connect->variance, NULL,
					    connect, -1);

		integrateDerivativeDistances(connect, min_distance, false);
	}
}

void unmarkNode(Node * node, MiniConnection * localConnect)
{
	if (localConnect->frontReference != NULL
	    || localConnect->backReference != NULL) {
		if (getNodeStatus(node) > 0)
			setSingleNodeStatus(node, 10);
		else
			setSingleNodeStatus(node, -10);
	} else {
		setSingleNodeStatus(node, false);
		destroyNodeList(localConnect->nodeList);
		localConnect->frontReference = NULL;
		localConnect->backReference = NULL;
		localConnect->nodeList = NULL;
	}
}

void handicapNode(Node * node)
{
	if (getNodeStatus(node) > 0)
		setSingleNodeStatus(node, 10);
	else
		setSingleNodeStatus(node, -10);
}

static void absorbExtension(Node * node, Node * extension)
{
	Arc *arc;

	appendNodeGaps(node, extension, graph);
	appendDescriptors(node, extension);

	// Destroy old nodes    
	while (getArc(node) != NULL)
		destroyArc(getArc(node), graph);

	// Create new
	for (arc = getArc(extension); arc != NULL; arc = getNextArc(arc))
		createAnalogousArc(node, getDestination(arc), arc, graph);
}

NodeList *getMarkedNodeList()
{
	return markedNodes;
}

static void absorbExtensionInScaffold(Node * node, Node * source)
{
	IDnum nodeID = getNodeID(node);
	IDnum sourceID = getNodeID(source);
	IDnum sourceIndex = sourceID + nodeCount(graph);
	Node *twinSource = getTwinNode(source);
	IDnum twinSourceIndex = getNodeID(twinSource) + nodeCount(graph);
	Connection *connect, *original;
	Node *destination;
	IDnum destinationID;
	Coordinate distance_shift =
	    (getNodeLength(node) - getNodeLength(source)) / 2;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;
	MiniConnection *localConnect;
	Coordinate distance;
	double variance;
	IDnum direct_count;
	IDnum paired_count;

	while ((connect = scaffold[sourceIndex])) {
		destination = getTwinNode(connect->destination);

		if (destination == getTwinNode(node)) {
			localConnect = &localScaffold[twinSourceIndex];
			localConnect->frontReference = NULL;
			unmarkNode(twinSource, localConnect);
			destroyConnection(connect, sourceID);
			continue;
		}
		if (destination == node) {
			localConnect = &localScaffold[sourceIndex];
			localConnect->backReference = NULL;
			unmarkNode(source, localConnect);
			destroyConnection(connect, sourceID);
			continue;
		}

		destinationID = getNodeID(destination);
		localConnect =
		    &localScaffold[destinationID + nodeCount(graph)];
		connect->distance += distance_shift;
		distance = connect->distance;
		variance = connect->variance;
		direct_count = connect->direct_count;
		paired_count = connect->paired_count;

		if (getNodeStatus(destination)) {
			readjustMiniConnection(destination, localConnect,
					       distance, min_distance,
					       variance, NULL, NULL);
			if ((original = localConnect->frontReference))
				readjustConnection(original, distance,
						   variance, direct_count,
						   paired_count);
			else
				localConnect->frontReference =
				    createNewConnection(nodeID,
							-destinationID,
							direct_count, 
						        paired_count, distance,
							variance);
		} else
			resetMiniConnection(destination, localConnect,
					    distance, variance,
					    createNewConnection(nodeID,
								-destinationID,
								direct_count,
								paired_count,
								distance,
								variance),
					    NULL, true);

		integrateDerivativeDistances(connect, min_distance, true);

		destroyConnection(connect, sourceID);
	}

	// Loop thru twin's primary scaffold
	while ((connect = scaffold[twinSourceIndex])) {
		destination = connect->destination;

		if (destination == node) {
			localConnect = &localScaffold[sourceIndex];
			localConnect->frontReference = NULL;
			unmarkNode(source, localConnect);
			destroyConnection(connect, -sourceID);
			continue;
		}
		if (destination == getTwinNode(node)) {
			localConnect = &localScaffold[twinSourceIndex];
			localConnect->backReference = NULL;
			unmarkNode(twinSource, localConnect);
			destroyConnection(connect, -sourceID);
			continue;
		}

		destinationID = getNodeID(destination);

		localConnect =
		    &localScaffold[destinationID + nodeCount(graph)];
		connect->distance -= distance_shift;
		distance = connect->distance;
		variance = connect->variance;
		direct_count = connect->direct_count;
		paired_count = connect->paired_count;

		if (getNodeStatus(destination) < 0) {
			readjustMiniConnection(destination, localConnect,
					       -distance, min_distance,
					       variance, NULL, NULL);
			if ((original = localConnect->backReference))
				readjustConnection(original, distance,
						   variance, direct_count,
						   paired_count);
		} else if (getNodeStatus(destination) > 0) {
			if ((original = localConnect->frontReference)) {
				destroyConnection(original, nodeID);
				localConnect->frontReference = NULL;
			}
			unmarkNode(destination, localConnect);
		} else if (distance > min_distance)
			resetMiniConnection(destination, localConnect,
					    -distance, variance, NULL,
					    createNewConnection(-nodeID,
								destinationID,
								direct_count,
								paired_count,
								distance,
								variance),
					    -1);

		integrateDerivativeDistances(connect, min_distance, true);
		destroyConnection(connect, -sourceID);
	}
}

static void recenterNode(Node * node, Coordinate oldLength)
{
	IDnum nodeID = getNodeID(node);
	IDnum nodeIndex = nodeID + nodeCount(graph);
	IDnum twinNodeIndex = -nodeID + nodeCount(graph);
	Connection *connect, *next;
	Coordinate distance_shift = (getNodeLength(node) - oldLength) / 2;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;
	MiniConnection *localConnect;

	//puts("Recentering node");

	for (connect = scaffold[nodeIndex]; connect != NULL;
	     connect = next) {
		next = connect->next;
		connect->distance -= distance_shift;

		if (connect->distance < min_distance) {
			//printf("Unrecording %li\n",
			//       -getNodeID(connect->destination));
			localConnect =
			    &localScaffold[-getNodeID(connect->destination)
					   + nodeCount(graph)];
			localConnect->frontReference = NULL;
			unmarkNode(getTwinNode(connect->destination),
				   localConnect);
			destroyConnection(connect, nodeID);
		} else if (connect->twin != NULL)
			connect->twin->distance -= distance_shift;
	}

	for (connect = scaffold[twinNodeIndex]; connect != NULL;
	     connect = next) {
		next = connect->next;
		connect->distance += distance_shift;

		if (connect->twin != NULL)
			connect->twin->distance += distance_shift;
	}
}

static void recenterLocalScaffold(Node * node, Coordinate oldLength)
{
	MiniConnection *localConnect;
	Coordinate distance_shift = (getNodeLength(node) - oldLength) / 2;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;
	NodeList *nodeList, *next;
	IDnum node2ID;
	Node *node2;

	for (nodeList = markedNodes; nodeList != NULL; nodeList = next) {
		next = nodeList->next;

		node2 = nodeList->node;

		if (node2 == node) {
			setSingleNodeStatus(node2, 1);
			continue;
		}

		node2ID = getNodeID(node2);
		localConnect = &localScaffold[node2ID + nodeCount(graph)];
		localConnect->distance -= distance_shift;

		if (localConnect->distance < min_distance
		    && localConnect->backReference == NULL
		    && localConnect->frontReference == NULL)
			unmarkNode(node2, localConnect);
		else if (getNodeStatus(node2) > 0)
			setSingleNodeStatus(node2, 1);
		else if (getNodeStatus(node2) < 0)
			setSingleNodeStatus(node2, -1);
	}
}

static void adjustShortReads(Node * target, Node * source)
{
	ShortReadMarker *targetArray, *marker;
	IDnum targetLength, index;
	Coordinate position, nodeLength;

	if (!readStartsAreActivated(graph))
		return;

	targetArray = getNodeReads(getTwinNode(target), graph);
	targetLength = getNodeReadCount(getTwinNode(target), graph);

	nodeLength = getNodeLength(source);

	for (index = 0; index < targetLength; index++) {
		marker = getShortReadMarkerAtIndex(targetArray, index);
		position = getShortReadMarkerPosition(marker);
		position += nodeLength;
		setShortReadMarkerPosition(marker, position);
	}
}

static void adjustLongReads(Node * target, Node * source)
{
	PassageMarker *marker;
	Coordinate nodeLength = getNodeLength(source);

	for (marker = getMarker(source); marker != NULL;
	     marker = getNextInNode(marker))
		incrementFinishOffset(marker, nodeLength);
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

static boolean comesFromNode(PassageMarker * marker, Node * node)
{
	Node *target = getTwinNode(node);
	PassageMarker *current;

	for (current = getTwinMarker(marker); current != NULL;
	     current = getNextInSequence(current))
		if (getNode(current) == target)
			return true;

	return false;
}

static void reconnectPassageMarker(PassageMarker * marker, Node * node,
				   PassageMarker ** ptr)
{
	PassageMarker *current;
	PassageMarker *next = getNextInSequence(marker);
	PassageMarker *tmpMarker;

	for (current = marker; getNode(current) != node;
	     current = getPreviousInSequence(current));

	setPreviousInSequence(current, next);
	concatenatePassageMarkers(current, marker);

	for (; marker != current; marker = tmpMarker) {
		tmpMarker = getPreviousInSequence(marker);
		if (*ptr == marker || *ptr == getTwinMarker(marker))
			*ptr = getNextInNode(*ptr);
		setNextInSequence(marker, NULL);
		setPreviousInSequence(NULL, marker);
		destroyPassageMarker(marker);
	}
}

static void concatenateLongReads(Node * node, Node * candidate,
				 Graph * graph)
{
	PassageMarker *marker, *tmpMarker;

	// Passage marker management in node:
	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		if (!goesToNode(marker, candidate))
			incrementFinishOffset(marker,
					      getNodeLength(candidate));
	}

	// Swapping new born passageMarkers from candidate to node
	for (marker = getMarker(candidate); marker != NULL;
	     marker = tmpMarker) {
		tmpMarker = getNextInNode(marker);

		if (!comesFromNode(marker, node)) {
			extractPassageMarker(marker);
			transposePassageMarker(marker, node);
			incrementFinishOffset(getTwinMarker(marker),
					      getNodeLength(node));
		} else {
			reconnectPassageMarker(marker, node, &tmpMarker);
		}
	}
}

static void adjustShortReadsByLength(Node * target, Coordinate nodeLength)
{
	ShortReadMarker *targetArray, *marker;
	IDnum targetLength, index;
	Coordinate position;

	if (!readStartsAreActivated(graph))
		return;

	targetArray = getNodeReads(getTwinNode(target), graph);
	targetLength = getNodeReadCount(getTwinNode(target), graph);

	for (index = 0; index < targetLength; index++) {
		marker = getShortReadMarkerAtIndex(targetArray, index);
		position = getShortReadMarkerPosition(marker);
		position += nodeLength;
		setShortReadMarkerPosition(marker, position);
	}
}

static boolean abs_bool(boolean val)
{
	return val >= 0 ? val : -val;
}

static IDnum abs_ID(IDnum val)
{
	return val >= 0 ? val : -val;
}

static NodeList *pathIsClear(Node * node, Node * oppositeNode,
			     Coordinate distance)
{
	Arc *arc;
	Node *candidate, *dest, *current;
	Coordinate extension_distance = 0;
	boolean maxRepeat = 1;
	Node *repeatEntrance = NULL;
	IDnum counter = 0;
	NodeList *path = NULL;
	NodeList *tail = path;

	setSingleNodeStatus(node, 2);

	current = node;
	while (true) {

		//////////////////////////////////
		//  Selecting destination       //
		//////////////////////////////////
		candidate = NULL;

		// First round for priority nodes
		for (arc = getArc(current); arc != NULL;
		     arc = getNextArc(arc)) {
			dest = getDestination(arc);

			if (dest == node || dest == getTwinNode(node))
				continue;

			if (getNodeStatus(dest) <= 0)
				continue;

			if (candidate == NULL
			    || getNodeStatus(candidate) >
			    getNodeStatus(dest)
			    || (getNodeStatus(candidate) ==
				getNodeStatus(dest)
				&& extension_distance >
				localScaffold[getNodeID(dest) +
					      nodeCount(graph)].
				distance - getNodeLength(dest) / 2)) {
				extension_distance =
				    localScaffold[getNodeID(dest) +
						  nodeCount(graph)].
				    distance - getNodeLength(dest) / 2;
				candidate = dest;
			}
		}

		if (candidate != NULL && repeatEntrance) {
			for (arc = getArc(node); arc != NULL;
			     arc = getNextArc(arc)) {
				dest = getDestination(arc);
				if (dest != candidate
				    && getNodeStatus(dest)) {
					break;
				}
			}
		}
		// In case of failure   
		if (candidate == NULL) {
			for (arc = getArc(current); arc != NULL;
			     arc = getNextArc(arc)) {
				dest = getDestination(arc);

				if (getNodeStatus(dest) == 0)
					continue;

				if (dest == node
				    || dest == getTwinNode(node))
					continue;

				if (candidate == NULL
				    || getNodeStatus(candidate) <
				    getNodeStatus(dest)
				    || (getNodeStatus(candidate) ==
					getNodeStatus(dest)
					&& extension_distance <
					localScaffold[getNodeID(dest) +
						      nodeCount(graph)].
					distance -
					getNodeLength(dest) / 2)) {
					extension_distance =
					    localScaffold[getNodeID(dest) +
							  nodeCount
							  (graph)].
					    distance -
					    getNodeLength(dest) / 2;
					candidate = dest;
				}
			}
		}
		if (candidate == NULL) {
			while (path) {
				tail = path->next;
				deallocateNodeList(path);
				path = tail;
			}
			return false;
		}
		// Loop detection
		if (candidate == repeatEntrance
		    && abs_bool(getNodeStatus(candidate)) ==
		    maxRepeat + 1) {
			while (path) {
				tail = path->next;
				deallocateNodeList(path);
				path = tail;
			}
			return false;
		} else if (abs_bool(getNodeStatus(candidate)) > maxRepeat) {
			maxRepeat = abs_bool(getNodeStatus(candidate));
			repeatEntrance = candidate;
		} else if (abs_bool(getNodeStatus(candidate)) == 1) {
			maxRepeat = 1;
			repeatEntrance = NULL;
		}

		if (getNodeStatus(candidate) > 0)
			setSingleNodeStatus(candidate,
					    getNodeStatus(candidate) + 1);
		else
			setSingleNodeStatus(candidate,
					    getNodeStatus(candidate) - 1);


		// DEBUG 
		if (abs_bool(getNodeStatus(candidate)) > 100
		    || counter > nodeCount(graph)) {
			while (path) {
				tail = path->next;
				deallocateNodeList(path);
				path = tail;
			}
			return false;
		}
		// DEBUG 
		if (candidate == node || candidate == getTwinNode(node))
			abort();

		// Missassembly detection
		if (getUniqueness(candidate) && oppositeNode
		    && candidate != oppositeNode
		    && extension_distance > distance) {
			while (path) {
				tail = path->next;
				deallocateNodeList(path);
				path = tail;
			}
			return false;
		}

		if (path == NULL) {
			path = allocateNodeList();
			path->next = NULL;
			path->node = candidate;
			tail = path;
		} else {
			tail->next = allocateNodeList();
			tail = tail->next;
			tail->node = candidate;
			tail->next = NULL;
		}

		if (getUniqueness(candidate))
			return path;

		current = candidate;
	}
}

static boolean pushNeighbours(Node * node, Node * oppositeNode,
			      Coordinate distance, boolean force_jumps)
{
	Node *candidate;
	Node *lastCandidate = NULL;
	Coordinate oldLength = getNodeLength(node);
	Category cat;
	MiniConnection *localConnect;
	NodeList *path, *tmp;

	if ((path = pathIsClear(node, oppositeNode, distance))) {
		while (path) {
			candidate = path->node;
			tmp = path->next;
			deallocateNodeList(path);
			path = tmp;

			///////////////////////////////////////
			//  Stepping forward to destination  //
			///////////////////////////////////////

			if (getUniqueness(candidate)) {
				concatenateReadStarts(node, candidate,
						      graph);
				concatenateLongReads(node, candidate,
						     graph);
				absorbExtension(node, candidate);

				// Scaffold changes
				recenterNode(node, oldLength);
				recenterLocalScaffold(node, oldLength);
				absorbExtensionInScaffold(node, candidate);

				// Read coverage
				for (cat = 0; cat < CATEGORIES; cat++) {
					incrementVirtualCoverage(node, cat,
								 getVirtualCoverage
								 (candidate,
								  cat));
					incrementOriginalVirtualCoverage(node, cat,
									 getOriginalVirtualCoverage
									 (candidate, cat));
				}

				if (getNodeStatus(candidate)) {
					localConnect =
					    &localScaffold[getNodeID
							   (candidate) +
							   nodeCount
							   (graph)];
					if (localConnect->frontReference) {
						destroyConnection
						    (localConnect->
						     frontReference,
						     getNodeID(node));
						localConnect->
						    frontReference = NULL;
					}
					if (localConnect->backReference) {
						destroyConnection
						    (localConnect->
						     backReference,
						     -getNodeID(node));
						localConnect->
						    backReference = NULL;
					}
					unmarkNode(candidate,
						   localConnect);
				}
				if (getNodeStatus(getTwinNode(candidate))) {
					localConnect =
					    &localScaffold[-getNodeID
							   (candidate) +
							   nodeCount
							   (graph)];
					if (localConnect->frontReference) {
						destroyConnection
						    (localConnect->
						     frontReference,
						     getNodeID(node));
						localConnect->
						    frontReference = NULL;
					}
					if (localConnect->backReference) {
						destroyConnection
						    (localConnect->
						     backReference,
						     -getNodeID(node));
						localConnect->
						    backReference = NULL;
					}
					unmarkNode(getTwinNode(candidate),
						   localConnect);
				}
				destroyNode(candidate, graph);
				return true;
			} else {
				adjustShortReads(node, candidate);
				adjustLongReads(node, candidate);
				absorbExtension(node, candidate);
				lastCandidate = candidate;
			}
		}
	}

	if (force_jumps && oppositeNode
	    && abs_ID(getNodeID(oppositeNode)) < abs_ID(getNodeID(node))) {
		distance -= getNodeLength(node) / 2;
		distance -= getNodeLength(oppositeNode) / 2;
		if (distance > 1) {
			adjustShortReadsByLength(node, distance);
			appendGap(node, distance, graph);
		} else {
			adjustShortReadsByLength(node, 1);
			appendGap(node, 1, graph);
		}

		concatenateReadStarts(node, oppositeNode, graph);
		concatenateLongReads(node, oppositeNode, graph);
		absorbExtension(node, oppositeNode);

		// Scaffold changes
		recenterNode(node, oldLength);
		recenterLocalScaffold(node, oldLength);
		absorbExtensionInScaffold(node, oppositeNode);

		// Read coverage
		for (cat = 0; cat < CATEGORIES; cat++)
			incrementVirtualCoverage(node, cat,
						 getVirtualCoverage
						 (oppositeNode, cat));

		if (getNodeStatus(oppositeNode)) {
			localConnect =
			    &localScaffold[getNodeID(oppositeNode) +
					   nodeCount(graph)];
			if (localConnect->frontReference) {
				destroyConnection(localConnect->
						  frontReference,
						  getNodeID(node));
				localConnect->frontReference = NULL;
			}
			if (localConnect->backReference) {
				destroyConnection(localConnect->
						  backReference,
						  -getNodeID(node));
				localConnect->backReference = NULL;
			}
			unmarkNode(oppositeNode, localConnect);
		}
		if (getNodeStatus(getTwinNode(oppositeNode))) {
			localConnect =
			    &localScaffold[-getNodeID(oppositeNode) +
					   nodeCount(graph)];
			if (localConnect->frontReference) {
				destroyConnection(localConnect->
						  frontReference,
						  getNodeID(node));
				localConnect->frontReference = NULL;
			}
			if (localConnect->backReference) {
				destroyConnection(localConnect->
						  backReference,
						  -getNodeID(node));
				localConnect->backReference = NULL;
			}
			unmarkNode(getTwinNode(oppositeNode),
				   localConnect);
		}

		destroyNode(oppositeNode, graph);
	}

	return false;
}

static void unmarkInterestingNodes()
{
	Node *node;
	MiniConnection *localConnect;

	while ((node = popNodeRecord())) {
		setSingleNodeStatus(node, false);
		localConnect =
		    &localScaffold[getNodeID(node) + nodeCount(graph)];
		localConnect->frontReference = NULL;
		localConnect->backReference = NULL;
		localConnect->nodeList = NULL;
	}
}

static void findOppositeNode(Node * node, Node ** oppositeNode,
			     Coordinate * distance)
{
	NodeList *nodeList;
	MiniConnection *localConnect;
	Node *node2;
	Coordinate min_distance =
	    getNodeLength(node) / 2 - BACKTRACK_CUTOFF;
	IDnum node2ID;

	*oppositeNode = NULL;
	*distance = 0;

	for (nodeList = markedNodes; nodeList != NULL;
	     nodeList = nodeList->next) {
		node2 = nodeList->node;
		node2ID = getNodeID(node2);
		localConnect = &localScaffold[node2ID + nodeCount(graph)];

		if (node2 == node)
			continue;

		if (!getUniqueness(node2))
			continue;

		if (localConnect->distance < min_distance)
			continue;

		if (*oppositeNode == NULL
		    || *distance > localConnect->distance) {
			*oppositeNode = node2;
			*distance = localConnect->distance;
		}
	}
}

static boolean expandLongNode(Node * node, boolean force_jumps)
{
	boolean hit = true;
	boolean modified = false;
	Node *oppositeNode;
	Coordinate distance = 0;

	markInterestingNodes(node);

	while (hit) {
		correctGraphLocally(node);
		findOppositeNode(node, &oppositeNode, &distance);
		hit =
		    pushNeighbours(node, oppositeNode, distance,
				   force_jumps);
		modified = modified || hit;
	}

	unmarkInterestingNodes();

	return modified;
}

static boolean expandLongNodes(boolean force_jumps)
{
	IDnum nodeID;
	Node *node;
	boolean modified = false;

	for (nodeID = 1; nodeID <= nodeCount(graph); nodeID++) {
		node = getNodeInGraph(graph, nodeID);

		if (node != NULL && getUniqueness(node)) {
			modified = expandLongNode(node, force_jumps)
			    || modified;
			modified =
			    expandLongNode(getTwinNode(node), force_jumps)
			    || modified;
		}
	}

	return modified;
}

static void cleanMemory(ReadOccurence ** readNodes, IDnum * readNodeCounts,
			Coordinate * lengths)
{
	IDnum index;

	puts("Cleaning memory");

	for (index = 1; index <= sequenceCount(graph); index++)
		free(readNodes[index]);

	free(readNodes);
	free(readNodeCounts);

	destroyRecycleBin(connectionMemory);
	connectionMemory = NULL;
	free(scaffold);

	destroyRecycleBin(nodeListMemory);
	nodeListMemory = NULL;

	free(localScaffold);

	free(lengths);
}

static IDnum ** countShortReads(Graph * graph) {
	IDnum ** counts = calloc(CATEGORIES + 1, sizeof(IDnum*)); 
	Category cat;
	IDnum nodeIndex;
	IDnum nodes = nodeCount(graph);
	Node * node;
	ShortReadMarker * array, *marker;
	IDnum readCount, readIndex, readID;

	// Allocate memory where needed
	for (cat = 0; cat <= CATEGORIES; cat++) 
		if (getInsertLength(graph, cat) > 0)
			counts[cat] = calloc(2 * nodeCount(graph) + 1, sizeof(IDnum));

	// Start fillin'
	for (nodeIndex = 0; nodeIndex < 2 * nodes + 1; nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex - nodes);	

		if (node == NULL || !getUniqueness(node))
			continue;

		array = getNodeReads(node, graph);
		readCount = getNodeReadCount(node, graph);
		for (readIndex = 0; readIndex < readCount; readIndex++) {
			marker = getShortReadMarkerAtIndex(array, readIndex);
			readID = getShortReadMarkerID(marker);
			cat = reads->categories[readID - 1];
			if (cat % 2 == 1 && counts[cat/2] != NULL)
				counts[cat/2][nodeIndex]++;
		}	
	}

	return counts;
}

static void removeUnreliableConnections()
{
	IDnum maxNodeIndex = nodeCount(graph) * 2 + 1;
	IDnum index;
	Connection *connect, *next;
	Category cat;
	IDnum ** counts = countShortReads(graph);
	IDnum nodes = nodeCount(graph);

	/*
	   Node* node;
	   puts("CONNECT IDA IDB dcount pcount dist lengthA lengthB var countA countB coordA coordB test");
	   for (index = 0; index < maxNodeIndex; index++) {
	   node = getNodeInGraph(graph, index - nodeCount(graph));
	   for (connect = scaffold[index]; connect != NULL;
	   connect = next) {
	   next = connect->next;
	   if (getUniqueness(connect->destination)) {
	   printf("CONNECT %li %li %li %li %li %li %li %f %li %li", index - nodeCount(graph), getNodeID(connect->destination), connect->direct_count, connect->paired_count, connect->distance, getNodeLength(node), getNodeLength(connect->destination), connect->variance, getNodeReadCount(node, graph), getNodeReadCount(connect->destination, graph) );
	if (markerCount(node) == 1 && markerCount(connect->destination) == 1)
		printf(" %li %li", getPassageMarkerFinish(getMarker(node)), getPassageMarkerFinish(getMarker(connect->destination)));
	else
		printf(" ? ?");
	if (testConnection(index - nodes, connect, counts))
		puts(" OK");
	else
		puts(" NG");
	   }
	   }
	   }
	*/

	for (index = 0; index < maxNodeIndex; index++) {
		for (connect = scaffold[index]; connect != NULL;
		     connect = next) {
			next = connect->next;
			if (!testConnection(index - nodes, connect, counts))
				destroyConnection(connect, index-nodes);	
		}
	}

	// Free memory
	for (cat = 0; cat <= CATEGORIES; cat++)
		if (counts[cat])
			free(counts[cat]);
	free(counts);
}

void exploitShortReadPairs(Graph * argGraph, ReadSet * argReads,
			   boolean * dubious, double exp_cov,
			   boolean force_jumps)
{
	IDnum *readPairs;
	Category *cats;
	IDnum *readNodeCounts;
	ReadOccurence **readNodes;
	boolean modified = true;
	Coordinate *lengths =
	    getSequenceLengths(argReads, getWordLength(argGraph));

	// Globals
	graph = argGraph;
	reads = argReads;
	expected_coverage = exp_cov;
	readPairs = reads->mateReads;
	cats = reads->categories;

	if (!readStartsAreActivated(graph))
		return;

	puts("Starting pebble resolution...");

	// Prepare graph
	resetNodeStatus(graph);
	prepareGraphForLocalCorrections(graph);

	// Memory allocation
	localScaffold =
	    calloc((2 * nodeCount(graph) + 1), sizeof(MiniConnection));
	if (localScaffold == NULL) {
		puts("Calloc failure");
		exit(1);
	}
	// Prepare primary scaffold
	readNodeCounts = computeReadToNodeCounts();
	readNodes = computeReadToNodeMappings(readNodeCounts);
	scaffold =
	    computeNodeToNodeMappings(readNodes, readNodeCounts,
				      readPairs, cats, dubious, lengths);
	removeUnreliableConnections();

	// Loop until convergence
	while (modified)
		modified = expandLongNodes(force_jumps);

	// Clean up memory
	cleanMemory(readNodes, readNodeCounts, lengths);
	deactivateLocalCorrectionSettings();

	sortGapMarkers(graph);

	puts("Pebble done.");
}

void setUnreliableConnectionCutoff(int val)
{
	UNRELIABLE_CONNECTION_CUTOFF = (IDnum) val;
}
