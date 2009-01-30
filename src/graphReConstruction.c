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
#include <string.h>
#include <limits.h>

#include "globals.h"
#include "graph.h"
#include "passageMarker.h"
#include "readSet.h"
#include "tightString.h"
#include "recycleBin.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

typedef struct kmerOccurence_st KmerOccurence;
typedef struct kmerOccurenceTable_st KmerOccurenceTable;
typedef struct smallNodeList_st SmallNodeList;

// Internal structure used to mark the ends of an Annotation
struct kmerOccurence_st {
	Kmer kmer;
	Coordinate position;
	IDnum nodeID;
};

struct kmerOccurenceTable_st {
	KmerOccurence *kmerTable;
	IDnum kmerTableSize;
	IDnum *accelerationTable;
	short int accelerationBits;
	short int accelerationShift;
	Kmer accelerationMask;
};

struct smallNodeList_st {
	Node *node;
	SmallNodeList *next;
};

static RecycleBin *smallNodeListMemory = NULL;
static SmallNodeList *nodePile = NULL;

#define BLOCKSIZE 1000

static SmallNodeList *allocateSmallNodeList()
{
	if (smallNodeListMemory == NULL)
		smallNodeListMemory =
		    newRecycleBin(sizeof(SmallNodeList), BLOCKSIZE);

	return allocatePointer(smallNodeListMemory);
}

static void deallocateSmallNodeList(SmallNodeList * smallNodeList)
{
	deallocatePointer(smallNodeListMemory, smallNodeList);
}

static void memorizeNode(Node * node)
{
	SmallNodeList *list = allocateSmallNodeList();
	list->node = node;
	list->next = nodePile;
	nodePile = list;
}

static void unlockMemorizedNodes()
{
	SmallNodeList *list;

	while (nodePile) {
		list = nodePile;
		nodePile = list->next;
		setSingleNodeStatus(list->node, false);
		deallocateSmallNodeList(list);
	}
}

int compareKmerOccurences(void const *A, void const *B)
{
	KmerOccurence *a = (KmerOccurence *) A;
	KmerOccurence *b = (KmerOccurence *) B;

	if (a->kmer < b->kmer)
		return -1;
	else if (a->kmer > b->kmer)
		return 1;
	else
		return 0;
}

static inline Kmer keyInAccelerationTable(Kmer kmer,
					  KmerOccurenceTable * table)
{
	Kmer key = kmer;
	key &= table->accelerationMask;
	key >>= table->accelerationShift;
	return key;
}

static KmerOccurenceTable *referenceGraphKmers(char *preGraphFilename,
					       short int accelerationBits)
{
	FILE *file = fopen(preGraphFilename, "r");
	const int maxline = MAXLINE;
	char line[MAXLINE];
	int wordLength;
	Coordinate nodeLength, lineLength, kmerCount;
	Kmer word, antiWord, wordFilter;
	KmerOccurenceTable *kmerTable = NULL;
	KmerOccurence *kmerOccurences, *kmerOccurencePtr;
	Coordinate kmerOccurenceIndex;
	IDnum index;
	IDnum nodeID = 0;
	IDnum *accelPtr = NULL;
	Kmer lastHeader = 0;
	Kmer header;

	if (file == NULL) {
		printf("Could not open %s, sorry\n", preGraphFilename);
		exit(1);
	}
	// Count kmers
	printf("Scanning pre-graph file %s for k-mers\n",
	       preGraphFilename);

	// First  line
	fgets(line, maxline, file);
	sscanf(line, "%*i\t%*i\t%i\n", &wordLength);
	wordFilter = (((Kmer) 1) << (2 * wordLength)) - 1;

	// Initialize kmer occurence table:
	kmerTable = malloc(sizeof(KmerOccurenceTable));
	if (kmerTable == NULL) {
		puts("Malloc failure");
		exit(1);
	}
	if (accelerationBits > 2 * wordLength)
		accelerationBits = 2 * wordLength;

	if (accelerationBits > 32)
		accelerationBits = 32;

	if (accelerationBits > 0) {
		kmerTable->accelerationBits = accelerationBits;
		kmerTable->accelerationMask = ((Kmer) 1);
		kmerTable->accelerationMask <<= accelerationBits;
		kmerTable->accelerationMask--;
		kmerTable->accelerationMask <<= (2 * wordLength -
						 accelerationBits);
		kmerTable->accelerationTable =
		    calloc(((size_t) 1) << accelerationBits,
			   sizeof(IDnum));
		if (kmerTable->accelerationTable == NULL) {
			puts("Calloc failure");
			exit(1);
		}
		accelPtr = kmerTable->accelerationTable;
		kmerTable->accelerationShift =
		    (short int) 2 *wordLength - accelerationBits;
	} else {
		kmerTable->accelerationBits = 0;
		kmerTable->accelerationMask = 0;
		kmerTable->accelerationTable = NULL;
		kmerTable->accelerationShift = 0;
	}

	// Read nodes
	fgets(line, maxline, file);
	kmerCount = 0;
	while (line[0] == 'N') {
		fgets(line, maxline, file);
		kmerCount += (long) strlen(line) - wordLength;
		if (fgets(line, maxline, file) == NULL)
			break;
	}
	fclose(file);

	// Create table
	printf("%li kmers found\n", kmerCount);
	kmerOccurences = calloc(kmerCount, sizeof(KmerOccurence));
	if (kmerOccurences == NULL && kmerCount > 0) {
		puts("Calloc failure");
		exit(1);
	}
	kmerOccurencePtr = kmerOccurences;
	kmerOccurenceIndex = 0;
	kmerTable->kmerTable = kmerOccurences;
	kmerTable->kmerTableSize = kmerCount;

	// Fill table
	file = fopen(preGraphFilename, "r");
	if (file == NULL) {
		printf("Could not open %s, sorry.\n", preGraphFilename);
		exit(1);
	}
	fgets(line, maxline, file);

	// Read nodes
	fgets(line, maxline, file);
	while (line[0] == 'N') {
		fgets(line, maxline, file);
		nodeID++;
		lineLength = (long) strlen(line);
		nodeLength = lineLength - wordLength;

		// Fill in the initial word : 
		word = 0;
		for (index = 0; index < wordLength - 1; index++) {
			word <<= 2;
			if (line[index] == 'A')
				word += ADENINE;
			else if (line[index] == 'C')
				word += CYTOSINE;
			else if (line[index] == 'G')
				word += GUANINE;
			else if (line[index] == 'T')
				word += THYMINE;
		}

		// Scan through node
		for (; index < lineLength - 1; index++) {
			word <<= 2;
			if (line[index] == 'A')
				word += ADENINE;
			else if (line[index] == 'C')
				word += CYTOSINE;
			else if (line[index] == 'G')
				word += GUANINE;
			else if (line[index] == 'T')
				word += THYMINE;
			word &= wordFilter;

			antiWord = reverseComplement(word, wordLength);

			if (word <= antiWord) {
				kmerOccurencePtr->kmer = word;
				kmerOccurencePtr->nodeID = nodeID;
				kmerOccurencePtr->position =
				    index - wordLength + 1;
			} else {
				kmerOccurencePtr->kmer = antiWord;
				kmerOccurencePtr->nodeID = -nodeID;
				kmerOccurencePtr->position =
				    lineLength - 2 - index;
			}

			kmerOccurencePtr++;
			kmerOccurenceIndex++;
		}

		if (fgets(line, maxline, file) == NULL)
			break;
	}

	fclose(file);

	// Sort table
	qsort(kmerOccurences, kmerCount, sizeof(KmerOccurence),
	      compareKmerOccurences);

	// Fill up acceleration table
	if (kmerTable->accelerationTable != NULL) {
		*accelPtr = (IDnum) 0;
		for (kmerOccurenceIndex = 0;
		     kmerOccurenceIndex < kmerCount;
		     kmerOccurenceIndex++) {
			header =
			    keyInAccelerationTable(kmerOccurences
						   [kmerOccurenceIndex].
						   kmer, kmerTable);
			while (lastHeader < header) {
				lastHeader++;
				accelPtr++;
				*accelPtr = kmerOccurenceIndex;
			}
		}

		while (lastHeader < ((Kmer) 1 << accelerationBits) - 1) {
			lastHeader++;
			accelPtr++;
			*accelPtr = kmerCount;
		}
	}

	return kmerTable;
}

static KmerOccurence *findKmerOccurenceInSortedTable(Kmer kmer,
						     KmerOccurenceTable *
						     table)
{
	KmerOccurence *array = table->kmerTable;
	Kmer key = keyInAccelerationTable(kmer, table);
	Coordinate leftIndex, rightIndex, middleIndex;

	if (table->accelerationTable != NULL) {
		leftIndex = table->accelerationTable[key];
		rightIndex = table->accelerationTable[key + 1];
	} else {
		leftIndex = 0;
		rightIndex = table->kmerTableSize;
	}

	while (true) {
		middleIndex = (rightIndex + leftIndex) / 2;

		if (leftIndex >= rightIndex) {
			return NULL;
		} else if ((array[middleIndex]).kmer == kmer) {
			return &(array[middleIndex]);
		} else if (leftIndex == middleIndex) {
			return NULL;
		} else if ((array[middleIndex]).kmer > kmer)
			rightIndex = middleIndex;
		else
			leftIndex = middleIndex;
	}
}

static void ghostThreadSequenceThroughGraph(TightString * tString,
					    KmerOccurenceTable *
					    kmerOccurences, Graph * graph,
					    IDnum seqID, Category category,
					    Kmer wordFilter,
					    boolean readTracking)
{
	Kmer word = 0;
	Kmer antiWord;
	Coordinate readNucleotideIndex;
	KmerOccurence *kmerOccurence;
	int wordLength = getWordLength(graph);

	Node *node;
	Node *previousNode = NULL;

	// Neglect any read which will not be short paired
	if ((!readTracking && category % 2 == 0)
	    || category / 2 >= CATEGORIES)
		return;

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < wordLength)
		return;

	// Verify that all short reads are reasonnably short
	if (getLength(tString) > USHRT_MAX) {
		printf("Short read of length %li, longer than limit %i\n",
		       getLength(tString), SHRT_MAX);
		puts("You should better declare this sequence as long, because it genuinely is!");
		exit(1);
	}
	// Allocate memory for the read pairs
	if (!readStartsAreActivated(graph))
		activateReadStarts(graph);

	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < wordLength - 1; readNucleotideIndex++) {
		word <<= 2;
		word += getNucleotide(readNucleotideIndex, tString);
	}

	// Go through sequence
	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		word <<= 2;
		word &= wordFilter;
		word += getNucleotide(readNucleotideIndex++, tString);

		// Compute reverse complement
		antiWord = reverseComplement(word, wordLength);

		// Search in table
		if (word <= antiWord
		    && (kmerOccurence =
			findKmerOccurenceInSortedTable(word,
						       kmerOccurences))) {
			node =
			    getNodeInGraph(graph, kmerOccurence->nodeID);
		} else if (word > antiWord
			   && (kmerOccurence =
			       findKmerOccurenceInSortedTable(antiWord,
							      kmerOccurences)))
		{
			node =
			    getNodeInGraph(graph, -kmerOccurence->nodeID);
		} else {
			node = NULL;
			if (previousNode)
				break;
		}

		previousNode = node;

		// Fill in graph
		if (node && !getNodeStatus(node)) {
			incrementReadStartCount(node, graph);
			setSingleNodeStatus(node, true);
			memorizeNode(node);
		}
	}

	unlockMemorizedNodes();
}

static void threadSequenceThroughGraph(TightString * tString,
				       KmerOccurenceTable * kmerOccurences,
				       Graph * graph,
				       IDnum seqID, Category category,
				       Kmer wordFilter,
				       boolean readTracking)
{
	Kmer word = 0;
	Kmer antiWord;
	Coordinate readNucleotideIndex;
	Coordinate kmerIndex;
	KmerOccurence *kmerOccurence;
	int wordLength = getWordLength(graph);

	PassageMarker *marker = NULL;
	PassageMarker *previousMarker = NULL;
	Node *node;
	Node *previousNode = NULL;
	Coordinate coord;
	Coordinate previousCoord = 0;

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < wordLength)
		return;

	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < wordLength - 1; readNucleotideIndex++) {
		word <<= 2;
		word += getNucleotide(readNucleotideIndex, tString);
	}

	// Go through sequence
	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		word <<= 2;
		word &= wordFilter;
		word += getNucleotide(readNucleotideIndex++, tString);

		// Compute reverse complement
		antiWord = reverseComplement(word, wordLength);

		// Search in table
		if (word <= antiWord
		    && (kmerOccurence =
			findKmerOccurenceInSortedTable(word,
						       kmerOccurences))) {
			node =
			    getNodeInGraph(graph, kmerOccurence->nodeID);
			coord = kmerOccurence->position;
		} else if (word > antiWord
			   && (kmerOccurence =
			       findKmerOccurenceInSortedTable(antiWord,
							      kmerOccurences)))
		{
			node =
			    getNodeInGraph(graph, -kmerOccurence->nodeID);
			coord =
			    getNodeLength(node) - kmerOccurence->position -
			    1;
		} else {
			node = NULL;
			if (previousNode) {
				break;
			}
		}

		// Fill in graph
		if (node) {
			kmerIndex = readNucleotideIndex - wordLength;

			if (previousNode == node
			    && previousCoord == coord - 1) {
				if (category / 2 >= CATEGORIES) {
					setPassageMarkerFinish(marker,
							       kmerIndex +
							       1);
					setFinishOffset(marker,
							getNodeLength(node)
							- coord - 1);
				} else {
					incrementVirtualCoverage(node,
								 category /
								 2, 1);
					incrementOriginalVirtualCoverage
					    (node, category / 2, 1);
				}

			} else {
				if (category / 2 >= CATEGORIES) {
					marker =
					    newPassageMarker(seqID,
							     kmerIndex,
							     kmerIndex + 1,
							     coord,
							     getNodeLength
							     (node) -
							     coord - 1);
					transposePassageMarker(marker,
							       node);
					connectPassageMarkers
					    (previousMarker, marker,
					     graph);
					previousMarker = marker;
				} else {
					if (readTracking) {
						if (!getNodeStatus(node)) {
							addReadStart(node,
								     seqID,
								     coord,
								     graph,
								     kmerIndex);
							setSingleNodeStatus
							    (node, true);
							memorizeNode(node);
						} else {
							blurLastShortReadMarker
							    (node, graph);
						}
					}

					incrementVirtualCoverage(node,
								 category /
								 2, 1);
					incrementOriginalVirtualCoverage
					    (node, category / 2, 1);
				}

				createArc(previousNode, node, graph);
			}

			previousNode = node;
			previousCoord = coord;
		}
	}

	unlockMemorizedNodes();
}

static void fillUpGraph(ReadSet * reads,
			KmerOccurenceTable * kmerOccurences, Graph * graph,
			boolean readTracking)
{
	IDnum readIndex;
	Kmer wordFilter = (((Kmer) 1) << (2 * getWordLength(graph))) - 1;
	Category category;

	resetNodeStatus(graph);

	for (readIndex = 0; readIndex < reads->readCount; readIndex++) {
		category = reads->categories[readIndex];
		ghostThreadSequenceThroughGraph(reads->
						tSequences[readIndex],
						kmerOccurences,
						graph, readIndex + 1,
						category, wordFilter,
						readTracking);
	}

	createNodeReadStartArrays(graph);

	for (readIndex = 0; readIndex < reads->readCount; readIndex++) {
		category = reads->categories[readIndex];

		if (readIndex % 100000 == 0)
			printf("Threading through reads %li / %li\n",
			       readIndex, reads->readCount);

		threadSequenceThroughGraph(reads->tSequences[readIndex],
					   kmerOccurences,
					   graph, readIndex + 1, category,
					   wordFilter, readTracking);
	}

	orderNodeReadStartArrays(graph);

	if (smallNodeListMemory != NULL)
		destroyRecycleBin(smallNodeListMemory);

	free(kmerOccurences->kmerTable);
	free(kmerOccurences->accelerationTable);
	free(kmerOccurences);
}

Graph *importPreGraph(char *preGraphFilename, ReadSet * reads,
		      boolean readTracking, short int accelerationBits)
{
	Graph *graph = readPreGraphFile(preGraphFilename);
	KmerOccurenceTable *kmerOccurences =
	    referenceGraphKmers(preGraphFilename, accelerationBits);
	fillUpGraph(reads, kmerOccurences, graph, readTracking);

	return graph;
}
