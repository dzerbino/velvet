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

struct splayTable_st {
	SplayTree **table;
	IDnum lastIndex;
	int WORDLENGTH;
	Kmer WORDFILTER;
};

SplayTable *newSplayTable(int WORDLENGTH)
{
	SplayTable *splayTable = malloc(sizeof(SplayTable));
	splayTable->WORDLENGTH = WORDLENGTH;
	splayTable->WORDFILTER = (((Kmer) 1) << (2 * WORDLENGTH)) - 1;
	splayTable->table = calloc(CRC_HASH_BUCKETS, sizeof(SplayTree *));
	splayTable->lastIndex = 0;

	if (splayTable == NULL || splayTable->table == NULL) {
		printf
		    ("could not allocate splay hash table of size %i\n",
		     CRC_HASH_BUCKETS);
		return NULL;
	}

	return splayTable;
}

void destroySplayTable(SplayTable * splayTable)
{
	puts("Destroying splay table");

	destroyAllSplayTrees();
	free(splayTable);

	puts("Splay table destroyed");
}

static int hash_kmer(Kmer kmer)
{
	return crc32_v((char *) &kmer, sizeof(Kmer));
}

static boolean findOrInsertOccurenceInSplayTable(Kmer kmer, IDnum * seqID,
						 Coordinate * position,
						 SplayTable * table)
{
	if (table == NULL) {
		puts("NULL table!");
		exit(1);
	}

	return findOrInsertOccurenceInSplayTree(kmer, seqID, position,
						&table->
						table[hash_kmer(kmer)]);
}

void inputSequenceIntoSplayTable(TightString * tString,
				 SplayTable * table, FILE * file)
{
	IDnum currentIndex;
	Coordinate readNucleotideIndex = 0;
	Coordinate writeNucleotideIndex = 0;
	Kmer word = 0;
	Kmer antiWord = 0;
	boolean annotationClosed = true;
	unsigned char nucleotide;
	IDnum sequenceID;
	Coordinate coord;
	boolean found;
	Coordinate position = 0;
	Coordinate start = 0;
	Coordinate finish = 0;
	IDnum referenceSequenceID = 0;

	table->lastIndex++;

	currentIndex = table->lastIndex;
	fprintf(file, "ROADMAP %li\n", currentIndex);

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < table->WORDLENGTH) {
		destroyTightString(tString);
		return;
	}
	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) {
		word <<= 2;
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		word += nucleotide;
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		readNucleotideIndex++;
		word <<= 2;
		word &= table->WORDFILTER;
		word += nucleotide;

		sequenceID = currentIndex;
		coord = writeNucleotideIndex;

		antiWord = reverseComplement(word, table->WORDLENGTH);

		if (word < antiWord) {
			found =
			    findOrInsertOccurenceInSplayTable(word,
							      &sequenceID,
							      &coord,
							      table);
		} else {
			sequenceID = -sequenceID;
			found =
			    findOrInsertOccurenceInSplayTable(antiWord,
							      &sequenceID,
							      &coord,
							      table);
			sequenceID = -sequenceID;
		}

		if (!found) {
			writeNucleotideIndex++;
			if (!annotationClosed)
				fprintf(file, "%li\t%li\t%li\t%li\n",
					referenceSequenceID, position,
					start, finish);
			annotationClosed = true;
		}
		// Other wise create/complete annotation:
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
				fprintf(file, "%li\t%li\t%li\t%li\n",
					referenceSequenceID, position,
					start, finish);

				referenceSequenceID = sequenceID;
				position = writeNucleotideIndex;
				start = finish = coord;

				if (referenceSequenceID > 0)
					finish++;
				else
					finish--;
			}
		}
	}

	if (!annotationClosed)
		fprintf(file, "%li\t%li\t%li\t%li\n", referenceSequenceID,
			position, start, finish);

	destroyTightString(tString);
	return;
}

void inputSequenceArrayIntoSplayTableAndArchive(ReadSet * reads,
						SplayTable * table,
						char *filename)
{
	IDnum index;
	IDnum sequenceCount = reads->readCount;
	TightString **array;
	FILE *outfile = fopen(filename, "w");

	if (outfile == NULL) {
		puts("Couldn't open file, sorry");
		exit(-1);
	} else
		printf("Writing into roadmap file %s...\n", filename);

	fprintf(outfile, "%li\t%i\n", sequenceCount, table->WORDLENGTH);

	if (reads->tSequences == NULL)
		convertSequences(reads);

	array = reads->tSequences;

	puts("Inputting sequences...");
	for (index = 0; index < sequenceCount; index++) {
		if (index % 100000 == 0) {
			printf("Inputting sequence %li / %li\n", index,
			       sequenceCount);
			fflush(stdout);
		}
		inputSequenceIntoSplayTable(array[index], table, outfile);
	}

	fclose(outfile);

	free(reads->tSequences);
	reads->tSequences = NULL;
	destroyReadSet(reads);
	puts("Done inputting sequences");
}

void inputMaskIntoSplayTable(TightString * tString, SplayTable * table)
{
	Coordinate readNucleotideIndex = 0;
	Kmer word = 0;
	Kmer antiWord = 0;
	unsigned char nucleotide;
	IDnum sequenceID = 0;
	Coordinate coord = 0;

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < table->WORDLENGTH) {
		destroyTightString(tString);
		return;
	}
	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) {
		word <<= 2;
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		word += nucleotide;
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		readNucleotideIndex++;
		word <<= 2;
		word &= table->WORDFILTER;
		word += nucleotide;

		antiWord = reverseComplement(word, table->WORDLENGTH);

		if (word < antiWord)
			findOrInsertOccurenceInSplayTable(word,
							  &sequenceID,
							  &coord, table);
		else
			findOrInsertOccurenceInSplayTable(antiWord,
							  &sequenceID,
							  &coord, table);
	}

	destroyTightString(tString);
	return;
}

void inputMaskArrayIntoSplayTable(ReadSet * reads, SplayTable * table)
{
	IDnum index;
	IDnum sequenceCount = reads->readCount;
	TightString **array;

	if (reads->tSequences == NULL)
		convertSequences(reads);

	array = reads->tSequences;

	puts("Loading masks...");
	for (index = 0; index < sequenceCount; index++)
		inputMaskIntoSplayTable(array[index], table);

	free(reads->tSequences);
	reads->tSequences = NULL;
	destroyReadSet(reads);
	puts("Done loading masks");
}
