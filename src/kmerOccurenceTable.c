/*
Copyright 2010 Daniel Zerbino (zerbino@ebi.ac.uk)

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
#include "utility.h"
#include "kmer.h"

// Internal structure used to mark the ends of an Annotation
struct kmerOccurence_st {
	IDnum position;
	IDnum nodeID;
	Kmer kmer;
} ATTRIBUTE_PACKED;

struct kmerOccurenceTable_st {
	KmerOccurence *kmerTable;
	KmerOccurence * kmerOccurencePtr;
	IDnum *accelerationTable;
	IDnum kmerTableSize;
	IDnum kmerOccurenceIndex;
	short int accelerationShift;
	short int accelerationBits;
};

int compareKmerOccurences(void const *A, void const *B)
{
	KmerOccurence *a = (KmerOccurence *) A;
	KmerOccurence *b = (KmerOccurence *) B;
	return compareKmers(&(a->kmer), &(b->kmer));
}

static inline KmerKey keyInAccelerationTable(Kmer * kmer,
					  KmerOccurenceTable * table)
{
	return getKmerKey(kmer);
}

KmerOccurence *findKmerInKmerOccurenceTable(Kmer * kmer,
						     KmerOccurenceTable *
						     table)
{
	KmerOccurence *array = table->kmerTable;
	KmerKey key = keyInAccelerationTable(kmer, table);
	Coordinate leftIndex, rightIndex, middleIndex;
	int diff;

	if (table->accelerationTable != NULL) {
		leftIndex = table->accelerationTable[key];
		rightIndex = table->accelerationTable[key + 1];
	} else {
		leftIndex = 0;
		rightIndex = table->kmerTableSize;
	}

	while (true) {
		middleIndex = (rightIndex + leftIndex) / 2;

		if (leftIndex >= rightIndex)
			return NULL;

		diff = compareKmers(&(array[middleIndex].kmer), kmer);

		if (diff == 0) {
			while (middleIndex > 0 && compareKmers(&(array[middleIndex - 1].kmer), kmer) == 0)
				middleIndex--;
			return &(array[middleIndex]);
		} else if (leftIndex == middleIndex)
			return NULL;
		else if (diff > 0)
			rightIndex = middleIndex;
		else
			leftIndex = middleIndex;
	}
}

KmerOccurenceTable * newKmerOccurenceTable(short int accelerationBits, int wordLength) {
	KmerOccurenceTable * kmerTable = mallocOrExit(1, KmerOccurenceTable);

	if (accelerationBits > 2 * wordLength)
		accelerationBits = 2 * wordLength;

	if (accelerationBits > 32)
		accelerationBits = 32;

	if (accelerationBits > 0) {
		resetKeyFilter(accelerationBits);	
		kmerTable->accelerationBits = accelerationBits;
		kmerTable->accelerationTable =
		    callocOrExit((((size_t) 1) << accelerationBits) + 1,
			   IDnum);
		kmerTable->accelerationShift =
		    (short int) 2 *wordLength - accelerationBits;
	} else {
		kmerTable->accelerationBits = 0;
		kmerTable->accelerationTable = NULL;
		kmerTable->accelerationShift = 0;
	}

	return kmerTable;
}

void allocateKmerOccurences(IDnum kmerCount, KmerOccurenceTable * table) {
	KmerOccurence * kmerOccurences = callocOrExit(kmerCount + 1, KmerOccurence);
	kmerOccurences[kmerCount].position = -1;
	kmerOccurences[kmerCount].nodeID = 0;

	table->kmerTable = kmerOccurences;
	table->kmerTableSize = kmerCount;
	table->kmerOccurencePtr = kmerOccurences;
	table->kmerOccurenceIndex = 0;
}

void recordKmerOccurence(Kmer * kmer, IDnum nodeID, Coordinate position, KmerOccurenceTable * table) {
	if (table->kmerOccurenceIndex == table->kmerTableSize)
		abort();
	copyKmers(&(table->kmerOccurencePtr->kmer), kmer);
	table->kmerOccurencePtr->nodeID = nodeID;
	table->kmerOccurencePtr->position = position;

	table->kmerOccurencePtr++;
	table->kmerOccurenceIndex++;
}

void sortKmerOccurenceTable(KmerOccurenceTable * table) {
	KmerKey lastHeader = 0;
	KmerKey header;
	IDnum *accelPtr = NULL;
	IDnum kmerOccurenceIndex;

	velvetLog("Sorting kmer occurence table ... \n");

	qsort(table->kmerTable, table->kmerTableSize, sizeof(KmerOccurence),
	      compareKmerOccurences);

	velvetLog("Sorting done.\n");

	// Fill up acceleration table
	if (table->accelerationTable != NULL) {
		accelPtr = table->accelerationTable;
		*accelPtr = (IDnum) 0;
		for (kmerOccurenceIndex = 0;
		     kmerOccurenceIndex < table->kmerTableSize;
		     kmerOccurenceIndex++) {
			header =
			    keyInAccelerationTable(&table->kmerTable
						   [kmerOccurenceIndex].
						   kmer, table);
			while (lastHeader < header) {
				lastHeader++;
				accelPtr++;
				*accelPtr = kmerOccurenceIndex;
			}
		}

		while (lastHeader < (KmerKey) 1 << table->accelerationBits) {
			lastHeader++;
			accelPtr++;
			*accelPtr = table->kmerTableSize;
		}
	}

}

KmerOccurence * getNextKmerOccurence(KmerOccurence * current) {
	register KmerOccurence * next = current + 1;
	if (next->nodeID == 0)
		return NULL;
	if (compareKmers(&current->kmer, &next->kmer))
		return NULL;
	return next;
}

void destroyKmerOccurenceTable(KmerOccurenceTable * kmerTable) {
	if (kmerTable == NULL)
		return; 

	free(kmerTable->kmerTable);
	free(kmerTable->accelerationTable);
	free(kmerTable);
}

IDnum getKmerOccurenceNodeID(KmerOccurence * occurence) {
	return occurence->nodeID;
} 

Coordinate getKmerOccurencePosition(KmerOccurence * occurence) {
	return occurence->position;
}
