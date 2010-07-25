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
#include <time.h>
#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#endif

#include "globals.h"
#include "readSet.h"
#include "splay.h"
#include "tightString.h"
#include "crc.h"
#include "utility.h"
#include "kmer.h"
#include "kmerOccurenceTable.h"

struct splayTable_st {
	SplayTree **table;
#ifdef OPENMP
	omp_lock_t *tableLocks;
#endif
	KmerOccurenceTable *kmerOccurenceTable;
	int WORDLENGTH;
	boolean double_strand;
};

SplayTable *newSplayTable(int WORDLENGTH, boolean double_strand)
{
	SplayTable *splayTable = mallocOrExit(1, SplayTable);
	splayTable->WORDLENGTH = WORDLENGTH;
	splayTable->table = callocOrExit(CRC_HASH_BUCKETS, SplayTree *);
	splayTable->kmerOccurenceTable = NULL;
	splayTable->double_strand = double_strand;
#ifdef OPENMP
	splayTable->tableLocks = mallocOrExit(CRC_HASH_BUCKETS, omp_lock_t);
	int i;
	#pragma omp parallel for
	for (i = 0; i < CRC_HASH_BUCKETS; i++)
		omp_init_lock(splayTable->tableLocks + i);
	initSplayTreeMemory();
#endif
	return splayTable;
}

void destroySplayTable(SplayTable * splayTable)
{
	velvetLog("Destroying splay table\n");

	destroyAllSplayTrees();
	free(splayTable->table);
	destroyKmerOccurenceTable(splayTable->kmerOccurenceTable);
	free(splayTable);

	velvetLog("Splay table destroyed\n");
}

static KmerKey hash_kmer(Kmer * kmer)
{
#if KMER_LONGLONGS
	KmerKey key = kmer->longlongs[0];

#if KMER_LONGLONGS > 1
	key ^= kmer->longlongs[1];
#endif
#if KMER_LONGLONGS > 2
	key ^= kmer->longlongs[2];
#endif

	key = (~key) + (key << 21);
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8);
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4);
	key = key ^ (key >> 28);
	key = key + (key << 31);

	return key % CRC_HASH_BUCKETS;
#elif KMER_LONGS
	KmerKey key = kmer->longs[0];

	key += ~(key << 15);
	key ^= (key >> 10);
	key += (key << 3);
	key ^= (key >> 6);
	key += ~(key << 11);
	key ^= (key >> 16);

	return key % CRC_HASH_BUCKETS;

#elif KMER_INTS
	return kmer->ints % CRC_HASH_BUCKETS;
#elif KMER_CHARS
	return kmer->chars % CRC_HASH_BUCKETS;
#endif
}

static Coordinate getNearestHSPIndex(Coordinate position, IDnum * sequenceIDs, Coordinate sequenceLength) {
	Coordinate back_offset = -1;
	Coordinate front_offset = -1;

	for (back_offset = 1; position - back_offset > 0; back_offset++) 
		if (sequenceIDs[position - back_offset])
			break;

	for (front_offset = 1; position + front_offset < sequenceLength; front_offset++) 
		if (sequenceIDs[position + front_offset])
			break;

	if (back_offset == position && position + front_offset == sequenceLength) 
		return -1;
	else if (back_offset == position)
		return position + front_offset;
	else if (front_offset + position == sequenceLength)
		return position - back_offset;
	else
		return back_offset < front_offset? position - back_offset : position + front_offset;
}

static KmerOccurence * getMostAppropriateHit(Coordinate readCoord, Coordinate readLength, boolean direct, KmerOccurence * kmerOccurence, IDnum mapCount, IDnum * mapSequenceID, Coordinate * mapCoord, int wordLength) {
	KmerOccurence * current;
	KmerOccurence * best = NULL;
	Coordinate expectedPosition;
	Coordinate positionError;
	IDnum mapIndex;

	// If only one hit
	if (!getNextKmerOccurence(kmerOccurence))
		return kmerOccurence;

	// If multiple hits by unmapped read
	if (mapCount == 0)
		return NULL;

	// Compare cases
	for (current = kmerOccurence; current; current = getNextKmerOccurence(current)) {
		for (mapIndex = 0; mapIndex < mapCount; mapIndex++) {

			// If wrong sequence or unconsistent orientation	
			if ((direct && getKmerOccurenceNodeID(current) != mapSequenceID[mapIndex]) 
			    || (!direct && getKmerOccurenceNodeID(current) != -mapSequenceID[mapIndex])) 
				continue;

			// Compute where it is supposed to land on reference
			if (mapSequenceID[mapIndex] < 0)
				expectedPosition = mapCoord[mapIndex] + readLength - readCoord - 1;
			else 
				expectedPosition = mapCoord[mapIndex] + readCoord - wordLength + 1;
		
			// Compute positional error
			positionError = getKmerOccurencePosition(current) - expectedPosition;

			// If potential hit record
			if (positionError < 1 && positionError > -1) {
				if (best)
					// If competing hit, give up
					return NULL;
				else
					// Record current hit
					best = current;
			}
		}
	}

	return best;
}

static inline boolean
doFindOrInsertOccurenceInSplayTree(Kmer * kmer, IDnum * seqID,
				   Coordinate * position, SplayTable *table)
{
#ifdef OPENMP
	const KmerKey kmerHash = hash_kmer(kmer);
	boolean ret;

	omp_set_lock(table->tableLocks + kmerHash);
	ret =  findOrInsertOccurenceInSplayTree(kmer, seqID, position,
						table->table + kmerHash);
	omp_unset_lock(table->tableLocks + kmerHash);

	return ret;
#else
	return findOrInsertOccurenceInSplayTree(kmer, seqID, position,
						&table->table[hash_kmer(kmer)]);
#endif
}


static boolean findOrInsertOccurenceInSplayTable(Kmer * kmer, IDnum * seqID,
						 Coordinate * position,
						 SplayTable * table, IDnum * sequenceIDs,
						 Coordinate * coords, Coordinate readIndex, Coordinate readLength, boolean direct)
{
	KmerOccurence * hit;
	Coordinate HSPIndex;

	// Check if previous anchor
	if (sequenceIDs && sequenceIDs[readIndex]) {
		if (direct)
			*seqID = sequenceIDs[readIndex];
		else 
			*seqID = -sequenceIDs[readIndex];
		if (sequenceIDs[readIndex] > 0) 
			*position = coords[readIndex] + readIndex;
		else
			*position = coords[readIndex] - readIndex + readLength - 1;

		return true;
	}
	else if (coords && coords[readIndex]) 
		// If in buffer zone:
		return doFindOrInsertOccurenceInSplayTree(kmer, seqID, position, table);

	// Look up first in reference sequence k-mers
	if (table->kmerOccurenceTable 
	    && (hit = findKmerInKmerOccurenceTable(kmer, table->kmerOccurenceTable))) {
		if (!getNextKmerOccurence(hit)) {
			*seqID = getKmerOccurenceNodeID(hit);
			*position = getKmerOccurencePosition(hit);
			return true;
		} else if ((HSPIndex = getNearestHSPIndex(*position, sequenceIDs, readLength)) > 0) {
	    		hit = getMostAppropriateHit(readIndex, readLength, direct, hit, 1, &(sequenceIDs[HSPIndex]), &(coords[HSPIndex]), table->WORDLENGTH);
			if (hit) {
				*seqID = getKmerOccurenceNodeID(hit);
				*position = getKmerOccurencePosition(hit);
				return true;
			}

		}
	} 

	// If not, go through the novel k-mers
	return doFindOrInsertOccurenceInSplayTree(kmer, seqID, position, table);
}

/* SF TODO This will be needed somewhere else, we should probably create a
 * StringBuffer class
 */
#define BUFFER_APPEND(buffer, bufferSize, currentSize, line) \
{ \
	const int lineSize = strlen(line); \
	currentSize += lineSize; \
	while (currentSize > bufferSize) \
	{ \
		bufferSize *= 2; \
		buffer = reallocOrExit (buffer, bufferSize, char); \
	} \
	buffer = strcat(buffer, line); \
}

static void printAnnotations(IDnum *sequenceIDs, Coordinate * coords,
			     TightString * tString, SplayTable * table,
			     FILE * file, boolean second_in_pair, IDnum seqID) 
{
	Coordinate readNucleotideIndex = 0;
	Coordinate writeNucleotideIndex = 0;
	Kmer word;
	Kmer antiWord;
	boolean annotationClosed = true;
	IDnum sequenceID;
	Coordinate coord;
	boolean found;
	Coordinate position = 0;
	Coordinate start = 0;
	Coordinate finish = 0;
	IDnum referenceSequenceID = 0;
	Nucleotide nucleotide;
	char *buffer;
	char lineBuffer[128];
	int bufferSize = 1024;
	int currentSize = 1;

	clearKmer(&word);
	clearKmer(&antiWord);

	buffer = mallocOrExit(bufferSize, char);
	buffer[0] = '\0';

	sprintf(lineBuffer, "ROADMAP %d\n", seqID);
	BUFFER_APPEND(buffer, bufferSize, currentSize, lineBuffer);

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < table->WORDLENGTH) {
#ifdef OPENMP
	#pragma omp critical
#endif
		fprintf(file, "%s", buffer);
		free(buffer);
		return;
	}

	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) { 
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif

		sequenceID = seqID;
		coord = writeNucleotideIndex;

		if (table->double_strand) {
			if (compareKmers(&word, &antiWord) <= 0) {
				found =
				    findOrInsertOccurenceInSplayTable(&word,
								      &sequenceID,
								      &coord,
								      table,
								      sequenceIDs, 
								      coords,
								      readNucleotideIndex,
								      getLength(tString),
								      true);
			} else {
				sequenceID = -sequenceID;
				found =
				    findOrInsertOccurenceInSplayTable(&antiWord,
								      &sequenceID,
								      &coord,
								      table, 
								      sequenceIDs, 
								      coords,
								      readNucleotideIndex,
								      getLength(tString),
								      false);
				sequenceID = -sequenceID;
			}
		} else {
			if (!second_in_pair) {
				found =
				    findOrInsertOccurenceInSplayTable(&word,
								      &sequenceID,
								      &coord,
								      table,
								      sequenceIDs, 
								      coords,
								      readNucleotideIndex,
								      getLength(tString),
								      true);
			} else {
				sequenceID = -sequenceID;
				found =
				    findOrInsertOccurenceInSplayTable(&antiWord,
								      &sequenceID,
								      &coord,
								      table, 
								      sequenceIDs, 
								      coords,
								      readNucleotideIndex,
								      getLength(tString),
								      false);
				sequenceID = -sequenceID;
			}
		}

		if (!found) {
			writeNucleotideIndex++;
			if (!annotationClosed) {
				sprintf(lineBuffer, "%ld\t%lld\t%lld\t%lld\n",
					(long) referenceSequenceID, (long long) position,
					(long long) start, (long long) finish);
				BUFFER_APPEND(buffer, bufferSize, currentSize, lineBuffer);
			}
			annotationClosed = true;
		}
		// Otherwise create/complete annotation:
		else {
			// Forbidden k-mer
			if (sequenceID == 0) {
				break;
			}
			// Closed/inexistant annotation
			else if (annotationClosed) {
				referenceSequenceID = sequenceID;
				position = writeNucleotideIndex;
				start = finish = coord;

				if (referenceSequenceID > 0)
					finish++;
				else
					finish--;

				annotationClosed = false;
			}
			// Open annotation
			else if (sequenceID == referenceSequenceID
				 && coord == finish) {
				if (referenceSequenceID > 0)
					finish++;
				else
					finish--;
			}
			// Previous non corresponding annotation
			else {
				sprintf(lineBuffer, "%ld\t%lld\t%lld\t%lld\n",
					(long) referenceSequenceID, (long long) position,
					(long long) start, (long long) finish);
				BUFFER_APPEND(buffer, bufferSize, currentSize, lineBuffer);

				referenceSequenceID = sequenceID;
				position = writeNucleotideIndex;
				start = finish = coord;

				if (referenceSequenceID > 0)
					finish++;
				else
					finish--;
			}
		}

		readNucleotideIndex++;
	}

	if (!annotationClosed) {
		sprintf(lineBuffer, "%ld\t%lld\t%lld\t%lld\n",
			(long) referenceSequenceID, (long long) position,
			(long long) start, (long long) finish);
		BUFFER_APPEND(buffer, bufferSize, currentSize, lineBuffer);
	}
#ifdef OPENMP
	#pragma omp critical
#endif
	fprintf(file, "%s", buffer);
	free(buffer);

	return;
}

static void computeClearHSPs(TightString * tString, FILE * seqFile, boolean second_in_pair, SplayTable * table, IDnum * sequenceIDs, Coordinate * coords) {
	Coordinate readNucleotideIndex = 0;
	Kmer word;
	Kmer antiWord;
	Nucleotide nucleotide;
	KmerOccurence * hit;
	char line[MAXLINE];
	
	Coordinate mapCount = 0;
	Coordinate maxCount = 10;	
	IDnum * mapReferenceIDs = callocOrExit(maxCount, IDnum);
	Coordinate * mapCoords = callocOrExit(maxCount, Coordinate);
	long long_var;
	long long longlong_var;
	int penalty;

	// Parse file for mapping info
	while (seqFile && fgets(line, MAXLINE, seqFile)) {
		if (line[0] == '>')
			break;
	
		if (line[0] == 'M') {
			sscanf(line,"M\t%li\t%lli\n", &long_var, &longlong_var);
			mapReferenceIDs[mapCount] = (IDnum) long_var;
			mapCoords[mapCount] = (Coordinate) longlong_var;

			if (++mapCount == maxCount) {
				maxCount *= 2;
				mapReferenceIDs = reallocOrExit(mapReferenceIDs, maxCount, IDnum);
				mapCoords = reallocOrExit(mapCoords, maxCount, Coordinate);
			}
		}
	}

	// First pass for unambiguous hits
	// Fill in the initial word : 
	clearKmer(&word);
	clearKmer(&antiWord);
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) { 
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);

#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif

		if (table->double_strand) {
			if (compareKmers(&word, &antiWord) <= 0) {
				hit = findKmerInKmerOccurenceTable(&word, table->kmerOccurenceTable); 

				if (hit && (hit = getMostAppropriateHit(readNucleotideIndex, getLength(tString), true, hit, mapCount, mapReferenceIDs, mapCoords, table->WORDLENGTH))) 
					sequenceIDs[readNucleotideIndex] = getKmerOccurenceNodeID(hit);
			} else {
				hit = findKmerInKmerOccurenceTable(&antiWord, table->kmerOccurenceTable); 

				if (hit && (hit = getMostAppropriateHit(readNucleotideIndex, getLength(tString), false, hit, mapCount, mapReferenceIDs, mapCoords, table->WORDLENGTH))) 
					sequenceIDs[readNucleotideIndex] = -getKmerOccurenceNodeID(hit);
			}
		} else {
			if (!second_in_pair) {
				hit = findKmerInKmerOccurenceTable(&word, table->kmerOccurenceTable); 

				if (hit && (hit = getMostAppropriateHit(readNucleotideIndex, getLength(tString), true, hit, mapCount, mapReferenceIDs, mapCoords, table->WORDLENGTH))) 
					sequenceIDs[readNucleotideIndex] = getKmerOccurenceNodeID(hit);
			} else {
				hit = findKmerInKmerOccurenceTable(&antiWord, table->kmerOccurenceTable); 

				if (hit && (hit = getMostAppropriateHit(readNucleotideIndex, getLength(tString), false, hit, mapCount, mapReferenceIDs, mapCoords, table->WORDLENGTH))) 
					sequenceIDs[readNucleotideIndex] = -getKmerOccurenceNodeID(hit);
			}
		}

		if (sequenceIDs[readNucleotideIndex]) {
			if (sequenceIDs[readNucleotideIndex] > 0) 
				coords[readNucleotideIndex] = getKmerOccurencePosition(hit) - readNucleotideIndex;
			else
				coords[readNucleotideIndex] = getKmerOccurencePosition(hit) + readNucleotideIndex - getLength(tString) + 1;
		}
	
		// Barrier to flip-flopping
		if (sequenceIDs[readNucleotideIndex - 1] != 0
		    && (sequenceIDs[readNucleotideIndex] != sequenceIDs[readNucleotideIndex - 1]
			|| coords[readNucleotideIndex] != coords[readNucleotideIndex - 1])) {
			// Break in continuity... skip k positions 
			sequenceIDs[readNucleotideIndex] = 0;
			coords[readNucleotideIndex] = -1;
			readNucleotideIndex++;

			for (penalty = 0; penalty < table->WORDLENGTH  - 1 && readNucleotideIndex < getLength(tString); penalty++) {
				nucleotide = getNucleotide(readNucleotideIndex, tString);
				pushNucleotide(&word, nucleotide);

#ifdef COLOR
				reversePushNucleotide(&antiWord, nucleotide);
#else
				reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
				sequenceIDs[readNucleotideIndex] = 0;
				coords[readNucleotideIndex] = -1;
				readNucleotideIndex++;
			}
		} else
			readNucleotideIndex++;
		
	}

	free(mapReferenceIDs);
	free(mapCoords);
}

void inputSequenceIntoSplayTable(TightString * tString,
				 SplayTable * table,
				 FILE * file, FILE * seqFile,
				 boolean second_in_pair,
				 IDnum seqID)
{
	Coordinate length = getLength(tString);
	IDnum * sequenceIDs = NULL;
	Coordinate * coords = NULL;

	// If appropriate, get the HSPs on reference sequences
	if (table->kmerOccurenceTable) {
		length = getLength(tString);
		sequenceIDs = callocOrExit(length, IDnum);
		coords = callocOrExit(length, Coordinate);
		computeClearHSPs(tString, seqFile, second_in_pair, table, sequenceIDs, coords);
	}
	
	// Go through read, eventually with annotations
	printAnnotations(sequenceIDs, coords, tString, table, file, second_in_pair, seqID);

	// Clean up
	if (sequenceIDs) {
		free(sequenceIDs);
		free(coords);
	}
}

static
void inputReferenceIntoSplayTable(TightString * tString,
				 SplayTable * table, FILE * file, IDnum seqID)
{
	IDnum currentIndex;
	Coordinate readNucleotideIndex = 0;
	Coordinate kmerIndex = 0;
	Kmer word;
	Kmer antiWord;
	Nucleotide nucleotide;

	clearKmer(&word);
	clearKmer(&antiWord);

	currentIndex = seqID;
#ifdef OPENMP
	#pragma omp critical
#endif
	fprintf(file, "ROADMAP %d\n", currentIndex);

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < table->WORDLENGTH) {
		return;
	}

	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) { 
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
		if (table->double_strand) {
#ifdef COLOR
			reversePushNucleotide(&antiWord, nucleotide);
#else
			reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
		}
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex++, tString);
		pushNucleotide(&word, nucleotide);

		if (table->double_strand) {
#ifdef COLOR
			reversePushNucleotide(&antiWord, nucleotide);
#else
			reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
		}

		if (table->double_strand) {
			if (compareKmers(&word, &antiWord) <= 0)
				recordKmerOccurence(&word, currentIndex, 
						    kmerIndex,
						    table->kmerOccurenceTable);
			else
				recordKmerOccurence(&antiWord, -currentIndex, 
							kmerIndex,
							table->kmerOccurenceTable);
		} else {
			recordKmerOccurence(&word, currentIndex, 
					    kmerIndex,
					    table->kmerOccurenceTable);
		}
		kmerIndex++;
	}

	return;
}

static Coordinate countReferenceKmers(ReadSet * reads, int wordLength) {
	IDnum readIndex;
	Coordinate length = 0;
	

	for (readIndex = 0; readIndex < reads->readCount && reads->categories[readIndex] == REFERENCE; readIndex++)	
	{
		Coordinate tmpLength = getLength(getTightStringInArray(reads->tSequences, readIndex));
		if (tmpLength >= wordLength)
			length += tmpLength - wordLength + 1;
	}

	return length;
}

void inputSequenceArrayIntoSplayTableAndArchive(ReadSet * reads,
						SplayTable * table,
						char *filename, char* seqFilename)
{
	IDnum index;
	IDnum sequenceCount = reads->readCount;
	TightString *array;
	FILE *outfile = fopen(filename, "w");
	FILE *seqFile = NULL;
	IDnum kmerCount;
	IDnum referenceSequenceCount = 0;
	char line[MAXLINE];
	struct timeval start, end, diff;
	boolean second_in_pair;

	if (outfile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Couldn't write to file %s", filename);
	else
		velvetLog("Writing into roadmap file %s...\n", filename);

	// Count reference sequences
	for (index = 0; index < reads->readCount && reads->categories[index] == REFERENCE; index++)
		referenceSequenceCount++;

	fprintf(outfile, "%ld\t%ld\t%i\t%hi\n", (long) sequenceCount, (long) referenceSequenceCount, table->WORDLENGTH, (short) table->double_strand);

	if (reads->tSequences == NULL)
		convertSequences(reads);

	if (referenceSequenceCount && (kmerCount = countReferenceKmers(reads, table->WORDLENGTH))> 0) {
		table->kmerOccurenceTable = newKmerOccurenceTable(24 , table->WORDLENGTH);
		allocateKmerOccurences(kmerCount, table->kmerOccurenceTable);
		seqFile = fopen(seqFilename, "r");
		
		if (seqFile == NULL)
			exitErrorf(EXIT_FAILURE, true, "Couldn't write to file %s", seqFilename);
		else
			velvetLog("Reading mapping info from file %s\n", seqFilename);

		for (index = 0; index < referenceSequenceCount + 1; index++) 
			while (fgets(line, MAXLINE, seqFile))
				if (line[0] == '>')
					break;
	}

	gettimeofday(&start, NULL);
	velvetLog("Inputting sequences...\n");
	array = reads->tSequences;
#ifdef OPENMP
	#pragma omp parallel for
#endif
	for (index = 0; index < sequenceCount; index++) {
		const Category category = reads->categories[index];

		// Prorgess report on screen
		if (index % 1000000 == 0) {
			velvetLog("Inputting sequence %d / %d\n", index,
			       sequenceCount);
			fflush(stdout);
		}

		// Before starting non-reference sequences, sort kmerOccurenceTable
		if (index > 0 
		    && reads->categories[index - 1] == REFERENCE
		    && reads->categories[index] != REFERENCE)
			sortKmerOccurenceTable(table->kmerOccurenceTable);
		// Test to make sure that all the reference reads are before all the other reads
		else if (index > 0
		    && reads->categories[index - 1] != REFERENCE
		    && reads->categories[index] == REFERENCE) {
			velvetLog("Reference sequence placed after a non-reference read!\n");
			velvetLog(">> Please re-order the filenames in your command line so as to have the reference sequence files before all the others\n");
#ifdef DEBUG 
			abort();
#endif 
			exit(0);
		}

		if (category % 2)
			second_in_pair = (index - reads->categoriesOffsets[category]) % 2;
		else 
			second_in_pair = false;

		// Hashing the reads
		if (reads->categories[index] == REFERENCE)
			// Reference reads
			inputReferenceIntoSplayTable(getTightStringInArray(array, index), table, outfile, index + 1);
		else
			// Normal reads
			inputSequenceIntoSplayTable(getTightStringInArray(array, index), table, outfile, seqFile, second_in_pair, index + 1);
	}
	gettimeofday(&end, NULL);
	timersub(&end, &start, &diff);
	printf(">>> Sequences loaded in %ld.%06ld s\n", diff.tv_sec, diff.tv_usec);

	fclose(outfile);
	if (seqFile)
		fclose(seqFile);

	free(reads->tSequences);
	reads->tSequences = NULL;
	destroyReadSet(reads);
	velvetLog("Done inputting sequences\n");
}
