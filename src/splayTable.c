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
#include "readSet.h"
#include "splay.h"
#include "tightString.h"
#include "crc.h"
#include "utility.h"
#include "kmer.h"
#include "kmerOccurenceTable.h"

struct splayTable_st {
	SplayTree **table;
	KmerOccurenceTable *kmerOccurenceTable;
	IDnum lastIndex;
	int WORDLENGTH;
	boolean double_strand;
};

SplayTable *newSplayTable(int WORDLENGTH, boolean double_strand)
{
	SplayTable *splayTable = mallocOrExit(1, SplayTable);
	splayTable->WORDLENGTH = WORDLENGTH;
	splayTable->table = callocOrExit(CRC_HASH_BUCKETS, SplayTree *);
	splayTable->lastIndex = 0;
	splayTable->kmerOccurenceTable = NULL;
	splayTable->double_strand = double_strand;
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

static int hash_kmer(Kmer * kmer)
{
	return crc32_v((char *) kmer, KMER_BYTE_SIZE);
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
		return findOrInsertOccurenceInSplayTree(kmer, seqID, position,
							&table->
							table[hash_kmer(kmer)]);

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
	return findOrInsertOccurenceInSplayTree(kmer, seqID, position,
						&table->
						table[hash_kmer(kmer)]);
}

static void printAnnotations(IDnum *sequenceIDs, Coordinate * coords, TightString * tString, SplayTable * table, FILE * file, boolean second_in_pair) 
{
	IDnum currentIndex;
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

	clearKmer(&word);
	clearKmer(&antiWord);

	table->lastIndex++;

	currentIndex = table->lastIndex;
	velvetFprintf(file, "ROADMAP %li\n", (long) currentIndex);

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

		sequenceID = currentIndex;
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
				velvetFprintf(file, "%ld\t%lld\t%lld\t%lld\n",
					(long) referenceSequenceID, (long long) position,
					(long long) start, (long long) finish);
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
				velvetFprintf(file, "%ld\t%lld\t%lld\t%lld\n",
					(long) referenceSequenceID, (long long) position,
					(long long) start, (long long) finish);

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
		velvetFprintf(file, "%ld\t%lld\t%lld\t%lld\n",
			(long) referenceSequenceID, (long long) position,
			(long long) start, (long long) finish);
	}

	return;
}

static void computeClearHSPs(TightString * tString, FILE * seqFile, boolean second_in_pair, SplayTable * table, IDnum * sequenceIDs, Coordinate * coords) {
	Coordinate readNucleotideIndex = 0;
	Kmer word;
	Kmer antiWord;
	Kmer polyA;
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

	clearKmer(&polyA);

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

	// Kill silly poly-T beginnings
	while (readNucleotideIndex < getLength(tString) && (compareKmers(&antiWord, &polyA) == 0 || compareKmers(&word, &polyA) == 0)) {
		nucleotide = getNucleotide(readNucleotideIndex++, tString);
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

			for (penalty = 0; readNucleotideIndex < getLength(tString) && (penalty < table->WORDLENGTH  - 1 || compareKmers(&word, &polyA) == 0 || compareKmers(&antiWord, &polyA) == 0); penalty++) {
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
				 SplayTable * table, FILE * file, FILE * seqFile, boolean second_in_pair)
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
	printAnnotations(sequenceIDs, coords, tString, table, file, second_in_pair);

	// Clean up
	if (sequenceIDs) {
		free(sequenceIDs);
		free(coords);
	}
}

void inputReferenceIntoSplayTable(TightString * tString,
				 SplayTable * table, FILE * file)
{
	IDnum currentIndex;
	Coordinate readNucleotideIndex = 0;
	Coordinate kmerIndex = 0;
	Kmer word;
	Kmer antiWord;
	Nucleotide nucleotide;

	clearKmer(&word);
	clearKmer(&antiWord);

	table->lastIndex++;

	currentIndex = table->lastIndex;
	velvetFprintf(file, "ROADMAP %li\n", (long) currentIndex);

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
	boolean second_in_pair = false;

	if (outfile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Couldn't write to file %s", filename);
	else
		velvetLog("Writing into roadmap file %s...\n", filename);

	// Count reference sequences
	for (index = 0; index < reads->readCount && reads->categories[index] == REFERENCE; index++)
		referenceSequenceCount++;

	velvetFprintf(outfile, "%ld\t%ld\t%i\t%hi\n", (long) sequenceCount, (long) referenceSequenceCount, table->WORDLENGTH, (short) table->double_strand);

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

	velvetLog("Inputting sequences...\n");
	array = reads->tSequences;
	for (index = 0; index < sequenceCount; index++) {
		// Prorgess report on screen
		if (index % 100000 == 0) {
			velvetLog("Inputting sequence %li / %li\n", (long) index,
			       (long) sequenceCount);
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

		// Hashing the reads
		if (reads->categories[index] == REFERENCE)
			// Reference reads
			inputReferenceIntoSplayTable(getTightStringInArray(array, index), table, outfile);
		else
			// Normal reads
			inputSequenceIntoSplayTable(getTightStringInArray(array, index), table, outfile, seqFile, second_in_pair);

		if (reads->categories[index] % 2) 
			second_in_pair = (second_in_pair? false : true);
		else 
			second_in_pair = false;

	}

	fclose(outfile);
	if (seqFile)
		fclose(seqFile);

	//free(reads->tSequences);
	//reads->tSequences = NULL;
	//destroyReadSet(reads);
	velvetLog("Done inputting sequences\n");
}
