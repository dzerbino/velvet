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
#include <math.h>
#include <time.h>

#include "globals.h"
#include "tightString.h"
#include "readSet.h"

#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#include "../third-party/zlib-1.2.3/Win32/include/zlib.h"
#else
#include "../third-party/zlib-1.2.3/zlib.h"
#endif

ReadSet *newReadSet()
{
	ReadSet *rs = calloc(1, sizeof(ReadSet));
	if (rs == NULL) {
		puts("Calloc failure");
		exit(1);
	}
	return rs;
}

ReadSet *newReadSetAroundTightStringArray(TightString ** array,
					  IDnum length)
{
	ReadSet *rs = newReadSet();
	rs->tSequences = array;
	rs->readCount = length;
	return rs;
}

void concatenateReadSets(ReadSet * A, ReadSet * B)
{
	ReadSet tmp;
	IDnum index;

	// Read count:  
	tmp.readCount = A->readCount + B->readCount;

	// Sequences
	if (A->sequences != NULL || B->sequences != NULL) {
		tmp.sequences = malloc(tmp.readCount * sizeof(char *));
		if (tmp.sequences == NULL && tmp.readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}

		if (A->sequences != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.sequences[index] = A->sequences[index];
			free(A->sequences);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.sequences[index] = NULL;

		if (B->sequences != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.sequences[A->readCount + index] =
				    B->sequences[index];
			free(B->sequences);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.sequences[A->readCount + index] = NULL;
	} else
		tmp.sequences = NULL;

	// tSequences
	if (A->tSequences != NULL || B->tSequences != NULL) {
		tmp.tSequences =
		    malloc(tmp.readCount * sizeof(TightString *));
		if (tmp.tSequences == NULL && tmp.readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}

		if (A->tSequences != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.tSequences[index] =
				    A->tSequences[index];
			free(A->tSequences);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.tSequences[index] = NULL;

		if (B->tSequences != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.tSequences[A->readCount + index] =
				    B->tSequences[index];
			free(B->tSequences);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.tSequences[A->readCount + index] =
				    NULL;
	} else
		tmp.tSequences = NULL;

	// Labels
	if (A->labels != NULL || B->labels != NULL) {
		tmp.labels = malloc(tmp.readCount * sizeof(char *));
		if (tmp.labels == NULL && tmp.readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}

		if (A->labels != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.labels[index] = A->labels[index];
			free(A->labels);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.labels[index] = NULL;

		if (B->labels != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.labels[A->readCount + index] =
				    B->labels[index];
			free(B->labels);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.labels[A->readCount + index] = NULL;
	} else
		tmp.labels = NULL;


	// Confidence scores
	if (A->confidenceScores != NULL || B->confidenceScores != NULL) {
		tmp.confidenceScores =
		    malloc(tmp.readCount * sizeof(Quality *));
		if (tmp.confidenceScores == NULL && tmp.readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}

		if (A->confidenceScores != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.confidenceScores[index] =
				    A->confidenceScores[index];
			free(A->confidenceScores);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.confidenceScores[index] = NULL;

		if (B->confidenceScores != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.confidenceScores[A->readCount +
						     index] =
				    B->confidenceScores[index];
			free(B->confidenceScores);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.confidenceScores[A->readCount +
						     index] = NULL;
	} else
		tmp.confidenceScores = NULL;

	// Kmer probabilities 
	if (A->kmerProbabilities != NULL || B->kmerProbabilities != NULL) {
		tmp.kmerProbabilities =
		    malloc(tmp.readCount * sizeof(Quality *));
		if (tmp.kmerProbabilities == NULL && tmp.readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}

		if (A->kmerProbabilities != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.kmerProbabilities[index] =
				    A->kmerProbabilities[index];
			free(A->kmerProbabilities);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.kmerProbabilities[index] = NULL;

		if (B->kmerProbabilities != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.kmerProbabilities[A->readCount +
						      index] =
				    B->kmerProbabilities[index];
			free(B->kmerProbabilities);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.kmerProbabilities[A->readCount +
						      index] = NULL;
	} else
		tmp.kmerProbabilities = NULL;

	// Mate reads 
	if (A->mateReads != NULL || B->mateReads != NULL) {
		tmp.mateReads = malloc(tmp.readCount * sizeof(IDnum));
		if (tmp.mateReads == NULL && tmp.readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}

		if (A->mateReads != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.mateReads[index] = A->mateReads[index];
			free(A->mateReads);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.mateReads[index] = 0;

		if (B->mateReads != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.mateReads[A->readCount + index] =
				    B->mateReads[index];
			free(B->mateReads);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.mateReads[A->readCount + index] = 0;
	} else
		tmp.mateReads = NULL;

	// Categories
	if (A->categories != NULL || B->categories != NULL) {
		tmp.categories = malloc(tmp.readCount * sizeof(Quality *));
		if (tmp.categories == NULL && tmp.readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}

		if (A->categories != NULL) {
			for (index = 0; index < A->readCount; index++)
				tmp.categories[index] =
				    A->categories[index];
			free(A->categories);
		} else
			for (index = 0; index < A->readCount; index++)
				tmp.categories[index] = CATEGORIES;

		if (B->categories != NULL) {
			for (index = 0; index < B->readCount; index++)
				tmp.categories[A->readCount + index] =
				    B->categories[index];
			free(B->categories);
		} else
			for (index = 0; index < B->readCount; index++)
				tmp.categories[A->readCount + index] =
				    CATEGORIES;
	} else
		tmp.categories = NULL;

	// Put everything back into A
	A->readCount = tmp.readCount;
	A->sequences = tmp.sequences;
	A->tSequences = tmp.tSequences;
	A->labels = tmp.labels;
	A->confidenceScores = tmp.confidenceScores;
	A->kmerProbabilities = tmp.kmerProbabilities;
	A->mateReads = tmp.mateReads;
	A->categories = tmp.categories;

	// Deallocate
	free(B);
}

void convertSequences(ReadSet * rs)
{
	rs->tSequences =
	    newTightStringArrayFromStringArray(rs->sequences,
					       rs->readCount);
	rs->sequences = NULL;
}

static Probability convertQualityScore(Quality score)
{
	return (Probability) 1 - pow(10, -score / ((double) 10));
}

void convertConfidenceScores(ReadSet * rs, int WORDLENGTH)
{
	Quality *baseCallerScores;
	Probability *kmerProbabilities;
	IDnum index;
	Coordinate position;
	Probability proba;

	rs->kmerProbabilities =
	    malloc(rs->readCount * sizeof(Probability *));
	if (rs->kmerProbabilities == NULL && rs->readCount > 0) {
		puts("Malloc failure");
		exit(1);
	}

	for (index = 0; index < rs->readCount; index++) {
		rs->kmerProbabilities[index] =
		    malloc((getLength(rs->tSequences[index]) - WORDLENGTH +
			    1) * sizeof(Probability));
		if (rs->kmerProbabilities[index] == NULL) {
			puts("Malloc failure");
			exit(1);
		}
		kmerProbabilities = rs->kmerProbabilities[index];
		baseCallerScores = rs->confidenceScores[index];

		proba = 1;
		for (position = 0;
		     position < getLength(rs->tSequences[index]);
		     position++) {
			proba *=
			    convertQualityScore(baseCallerScores
						[position]);
			if (position < WORDLENGTH)
				continue;

			proba /=
			    convertQualityScore(baseCallerScores
						[position - WORDLENGTH]);
			kmerProbabilities[position - WORDLENGTH + 1] =
			    proba;
		}

		rs->confidenceScores[index] = NULL;
		free(baseCallerScores);
	}

	free(rs->confidenceScores);
	rs->confidenceScores = NULL;
}

void categorizeReads(ReadSet * readSet, Category category)
{
	IDnum index;

	if (readSet->categories == NULL) {
		readSet->categories =
		    malloc(readSet->readCount * sizeof(Category));
		if (readSet->categories == NULL && readSet->readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}
	}

	for (index = 0; index < readSet->readCount; index++)
		readSet->categories[index] = category;
}

void simplifyReads(ReadSet * readSet)
{
	IDnum index;

	if (readSet->categories == NULL) {
		readSet->categories =
		    malloc(readSet->readCount * sizeof(Category));
		if (readSet->categories == NULL && readSet->readCount > 0) {
			puts("Malloc failure");
			exit(1);
		}
	}

	for (index = 0; index < readSet->readCount; index++) {
		if (readSet->categories[index] < CATEGORIES) {
			readSet->categories[index] /= 2;
			readSet->categories[index] *= 2;
		}
	}
}

void exportIDMapping(char *filename, ReadSet * reads)
{
	IDnum index;
	FILE *outfile = fopen(filename, "w");

	if (outfile == NULL) {
		printf("Couldn't open %s, sorry\n", filename);
		return;
	} else
		puts("Writing into file...");

	if (reads->labels == NULL) {
		fclose(outfile);
		return;
	}

	for (index = 0; index < reads->readCount; index++)
		if (reads->labels != NULL)
			fprintf(outfile, "s/SEQUENCE %li/%s/\n", index + 1,
				reads->labels[index]);

	fclose(outfile);

}

// Imports sequences from a fastq file 
// Memory space allocated within this function.
ReadSet *readSolexaFile(char *filename)
{
	FILE *file = fopen(filename, "r");
	IDnum lineCount = 0;
	IDnum readCount, readIndex;
	const int maxline = 500;
	char line[500];
	ReadSet *reads;

	if (file != NULL)
		printf("Reading Solexa file %s\n", filename);
	else {
		printf("Could not open %s\n", filename);
		exit(0);
	}

	reads = newReadSet();


	// Count lines:
	puts("Counting lines...");
	while (fgets(line, maxline, file) != NULL)
		if (strchr(line, '.') == NULL)
			lineCount++;

	readCount = lineCount;
	printf("%li reads found.\n", readCount);
	fclose(file);

	// Create table:
	reads->readCount = readCount;
	reads->sequences = malloc(readCount * sizeof(char *));
	if (reads->sequences == NULL && readCount > 0) {
		puts("Malloc failure");
		exit(1);
	}
	for (readIndex = 0; readIndex < readCount; readIndex++) {
		reads->sequences[readIndex] = malloc(100 * sizeof(char));
		if (reads->sequences[readIndex] == NULL) {
			puts("Malloc failure");
			exit(1);
		}
	}

	// Reopen file and memorize line:
	puts("Writing lines into string array...");
	file = fopen(filename, "r");
	readIndex = 0;
	while (fgets(line, maxline, file) != NULL)
		if (strchr(line, '.') == NULL) {
			sscanf(line, "%*i\t%*i\t%*i\t%*i\t%*c%[^\n]",
			       reads->sequences[readIndex]);
			readIndex++;
		}

	fclose(file);
	puts("Done");
	return reads;
}

ReadSet *readElandFile(char *filename)
{
	FILE *file = fopen(filename, "r");
	IDnum lineCount = 0;
	IDnum readCount, readIndex;
	const int maxline = 500;
	char line[500];
	ReadSet *reads;

	if (file != NULL)
		printf("Reading Eland file %s\n", filename);
	else {
		printf("Could not open %s\n", filename);
		exit(0);
	}

	reads = newReadSet();

	// Count lines:
	puts("Counting lines...");
	while (fgets(line, maxline, file) != NULL)
		lineCount++;

	readCount = lineCount;
	printf("%li reads found.\n", readCount);
	fclose(file);

	// Create table:
	reads->readCount = readCount;
	reads->sequences = malloc(readCount * sizeof(char *));
	if (reads->sequences == NULL && readCount > 0) {
		puts("Malloc failure");
		exit(1);
	}
	for (readIndex = 0; readIndex < readCount; readIndex++) {
		reads->sequences[readIndex] = malloc(100 * sizeof(char));
		if (reads->sequences[readIndex] == NULL) {
			puts("Malloc failure");
			exit(1);
		}
	}

	// Reopen file and memorize line:
	puts("Writing lines into string array...");
	file = fopen(filename, "r");
	readIndex = 0;
	while (fgets(line, maxline, file) != NULL) {
		sscanf(line, "%*[^\t]\t%[^\t\n]",
		       reads->sequences[readIndex]);
		readIndex++;
	}

	fclose(file);
	puts("Done");
	return reads;
}

void goToEndOfLine(char *line, FILE * file)
{
	size_t length = strlen(line);
	char c = line[length - 1];

	while (c != '\n')
		c = fgetc(file);
}

// Imports sequences from a fastq file 
// Memory space allocated within this function.
ReadSet *readFastQFile(char *filename)
{
	FILE *file = fopen(filename, "r");
	IDnum lineCount = 0;
	IDnum readCount, readIndex;
	const int maxline = 5000;
	char line[5000];
	ReadSet *reads;
	int i;

	if (file != NULL)
		printf("Reading FastQ file %s\n", filename);
	else {
		printf("Could not open %s\n", filename);
		exit(0);
	}

	reads = newReadSet();

	// Count lines:
	puts("Counting lines...");
	while (fgets(line, maxline, file) != NULL)
		lineCount++;
	readCount = lineCount / 4;
	printf("%li reads found.\n", readCount);
	fclose(file);

	// Create table:
	reads->readCount = readCount;
	reads->sequences = malloc(sizeof(char *) * (readCount));
	if (reads->sequences == NULL && readCount > 0) {
		puts("Malloc failure");
		exit(1);
	}
	// Reopen file and memorize line:
	puts("Writing lines into string array...");
	file = fopen(filename, "r");
	for (readIndex = 0; readIndex < readCount; readIndex++) {
		fgets(line, maxline, file);
		fgets(line, maxline, file);

		// newline stripping that will work on any platform
		for (i = strlen(line) - 1;
		     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
			line[i] = '\0';
		}

		reads->sequences[readIndex] = malloc(sizeof(char) * (strlen(line) + 1));	// Allocate enough space for null + line contents
		if (reads->sequences[readIndex] == NULL
		    && strlen(line) > 0) {
			puts("Malloc failure");
			exit(1);
		}

		strncpy(reads->sequences[readIndex], line, strlen(line) + 1);	// Copy line plus null terminating char

		fgets(line, maxline, file);
		fgets(line, maxline, file);
	}

	fclose(file);
	puts("Done");
	return reads;
}

// Imports sequences from a zipped rfastq file 
// Memory space allocated within this function.
ReadSet *readFastQGZFile(char *filename)
{
	gzFile file = gzopen(filename, "r");
	IDnum lineCount = 0;
	IDnum readCount, readIndex;
	const int maxline = 5000;
	char line[5000];
	ReadSet *reads;
	int i;

	if (file != NULL)
		printf("Reading zipped FastQ file %s\n", filename);
	else {
		printf("Could not open zipped file %s\n", filename);
		exit(0);
	}

	reads = newReadSet();

	// Count lines:
	puts("Counting lines...");
	while (gzgets(file, line, maxline) != NULL)
		lineCount++;
	readCount = lineCount / 4;
	printf("%li reads found.\n", readCount);
	gzclose(file);

	// Create table:
	reads->readCount = readCount;
	reads->sequences = malloc(sizeof(char *) * (readCount));
	if (reads->sequences == NULL && readCount > 0) {
		puts("Malloc failure");
		exit(1);
	}
	// Reopen file and memorize line:
	puts("Writing lines into string array...");
	file = gzopen(filename, "r");
	for (readIndex = 0; readIndex < readCount; readIndex++) {
		gzgets(file, line, maxline);
		gzgets(file, line, maxline);

		// newline stripping that will work on any platform
		for (i = strlen(line) - 1;
		     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
			line[i] = '\0';
		}

		reads->sequences[readIndex] = malloc(sizeof(char) * (strlen(line) + 1));	// Allocate enough space for null + line contents
		if (reads->sequences[readIndex] == NULL
		    && strlen(line) > 0) {
			puts("Malloc failure");
			exit(1);
		}
		strncpy(reads->sequences[readIndex], line, strlen(line) + 1);	// Copy line plus null terminating char

		gzgets(file, line, maxline);
		gzgets(file, line, maxline);
	}

	gzclose(file);
	puts("Done");
	return reads;
}

// Imports sequences from a fasta file 
// Memory is allocated within the function 
ReadSet *readFastAFile(char *filename)
{
	FILE *file = fopen(filename, "r");
	char *sequence = NULL;
	Coordinate bpCount = 0;
	const int maxline = 5000;
	char line[5000];
	IDnum sequenceCount, sequenceIndex;
	IDnum index;
	ReadSet *reads;

	if (file != NULL)
		printf("Reading FastA file %s;\n", filename);
	else {
		printf("Could not open %s\n", filename);
		exit(0);
	}

	reads = newReadSet();
	// Count number of separate sequences
	sequenceCount = 0;
	while (fgets(line, maxline, file) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	fclose(file);
	printf("%li sequences found\n", sequenceCount);

	reads->readCount = sequenceCount;
	reads->sequences = calloc(sequenceCount, sizeof(char *));
	if (reads->sequences == NULL && sequenceCount > 0) {
		puts("Malloc failure");
		exit(1);
	}
	if (sequenceCount == 0) {
		reads->sequences = NULL;
		return reads;
	}
	// Counting base pair length of each sequence:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file) != NULL) {
		if (line[0] == '>') {
			if (sequenceIndex != -1) {
				//printf("Sequence %li has length %li\n",
				//       sequenceIndex, bpCount);
				reads->sequences[sequenceIndex] =
				    calloc(sizeof(char), bpCount + 1);
				if (reads->sequences[sequenceIndex] ==
				    NULL) {
					puts("Allocation screwed up!");
					exit(1);
				}
			}
			sequenceIndex++;
			bpCount = 0;
		} else {
			bpCount += (long) strlen(line) - 1;
		}
	}

	//printf("Sequence %li has length %li\n", sequenceIndex, bpCount);
	reads->sequences[sequenceIndex] =
	    calloc(sizeof(char), bpCount + 1);
	if (reads->sequences[sequenceIndex] == NULL) {
		puts("Allocation screwed up!");
		exit(1);
	}
	fclose(file);

	// Reopen file and memorize line:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file)) {
		if (line[0] == '>') {
			if (sequenceIndex != -1) {
				sequence[bpCount] = '\0';
			}
			sequenceIndex++;
			bpCount = 0;
			//printf("Starting to read sequence %li\n",
			//       sequenceIndex);
			sequence = reads->sequences[sequenceIndex];
		} else {
			for (index = 0; index < (long) strlen(line) - 1;
			     index++)
				sequence[bpCount + index] = line[index];
			bpCount += (long) (strlen(line) - 1);
		}
	}

	if (sequenceIndex != -1) {
		sequence[bpCount] = '\0';
	}

	fclose(file);

	puts("Done");
	return reads;
}

// Imports sequences from a zipped fasta file 
// Memory is allocated within the function 
ReadSet *readFastAGZFile(char *filename)
{
	gzFile file = gzopen(filename, "r");
	char *sequence = NULL;
	Coordinate bpCount = 0;
	const int maxline = 100;
	char line[100];
	IDnum sequenceCount, sequenceIndex;
	IDnum index;
	ReadSet *reads;

	if (file != NULL)
		printf("Reading zipped FastA file %s;\n", filename);
	else {
		printf("Could not open %s\n", filename);
		exit(0);
	}

	reads = newReadSet();
	// Count number of separate sequences
	file = gzopen(filename, "r");
	sequenceCount = 0;
	while (gzgets(file, line, maxline) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	gzclose(file);
	printf("%li sequences found\n", sequenceCount);

	reads->readCount = sequenceCount;
	reads->sequences = malloc(sequenceCount * sizeof(char *));
	if (reads->sequences == NULL && sequenceCount > 0) {
		puts("Malloc failure");
		exit(1);
	}
	// Counting base pair length of each sequence:
	file = gzopen(filename, "r");
	sequenceIndex = -1;
	while (gzgets(file, line, maxline) != NULL) {
		if (line[0] == '>') {
			if (sequenceIndex != -1) {
				//printf("Sequence %li has length %li\n",
				//       sequenceIndex, bpCount);
				reads->sequences[sequenceIndex] =
				    malloc(sizeof(char) * (bpCount + 1));
				if (reads->sequences[sequenceIndex] ==
				    NULL) {
					puts("Allocation screwed up!");
					exit(1);
				}
			}
			sequenceIndex++;
			bpCount = 0;
		} else {
			bpCount += (long) strlen(line) - 1;
		}
	}

	//printf("Sequence %li has length %li\n", sequenceIndex, bpCount);
	reads->sequences[sequenceIndex] =
	    malloc(sizeof(char) * (bpCount + 1));
	if (reads->sequences[sequenceIndex] == NULL) {
		puts("Malloc failure");
		exit(1);
	}
	gzclose(file);

	// Reopen file and memorize line:
	file = gzopen(filename, "r");
	sequenceIndex = -1;
	while (gzgets(file, line, maxline)) {
		if (line[0] == '>') {
			if (sequenceIndex != -1) {
				sequence[bpCount] = '\0';
			}
			sequenceIndex++;
			bpCount = 0;
			//printf("Starting to read sequence %li\n",
			//       sequenceIndex);
			sequence = reads->sequences[sequenceIndex];
		} else {
			for (index = 0; index < (long) strlen(line) - 1;
			     index++)
				sequence[bpCount + index] = line[index];
			bpCount += (long) (strlen(line) - 1);
		}
	}

	if (sequenceIndex != -1) {
		sequence[bpCount] = '\0';
	}

	gzclose(file);

	puts("Done");
	return reads;
}

// Parser for new output
ReadSet *readMAQGZFile(char *filename)
{
	gzFile file = gzopen(filename, "r");
	const int maxline = 1000;
	char line[1000];
	IDnum sequenceCount, index;
	ReadSet *reads;

	if (file != NULL)
		printf("Reading zipped MAQ file %s;\n", filename);
	else {
		printf("Could not open %s\n", filename);
		exit(0);
	}

	reads = newReadSet();

	// Count number of separate sequences
	sequenceCount = 0;
	while (gzgets(file, line, maxline) != NULL)
		sequenceCount++;
	gzclose(file);
	printf("%li sequences found\n", sequenceCount);

	reads->readCount = sequenceCount;
	reads->sequences = malloc(sequenceCount * sizeof(char *));
	if (reads->sequences == NULL && sequenceCount > 0) {
		puts("Malloc failure");
		exit(1);
	}
	// Counting base pair length of each sequence:
	for (index = 0; index < sequenceCount; index++) {
		reads->sequences[index] = malloc(100 * sizeof(char));
		if (reads->sequences[index] == NULL) {
			puts("Malloc failure");
			exit(1);
		}
	}

	// Reopen file and memorize line:
	file = gzopen(filename, "r");
	index = 0;
	while (gzgets(file, line, maxline)) {
		sscanf(line, "%*s\t%*i\t%*i\t%*i\t%*i\t%*i\t\t%*i\t%[^\t]",
		       reads->sequences[index]);
		if (strspn(reads->sequences[index], "ATGC") < 21)
			strcpy(reads->sequences[index], "");
		index++;
	}

	gzclose(file);

	puts("Done");
	return reads;
}

#define FASTQ 1
#define FASTA 2
#define GERALD 3
#define ELAND 4
#define FASTA_GZ 5
#define FASTQ_GZ 6
#define MAQ_GZ 7

// General argument parser for most functions
// Basically a reused portion of toplevel code dumped into here
ReadSet *parseDataAndReadFiles(int argc, char **argv)
{
	int argIndex = 1;
	ReadSet *reads;
	ReadSet *allSequences = newReadSet();
	int filetype = FASTA;
	IDnum fileIndex = 0;
	Category cat = 0;

	if (argc < 2) {
		puts("Wrong number of arguments!");
		puts("Correct usage:");
		puts("run -<filetype> <list of files> [...] ");
		puts("Allowed filetypes:");
		puts("\t-fasta");
		puts("\t-fastq");
		puts("\t-solexa");
		puts("\t-eland");
		puts("If reading exclusively fasta file, the -fasta parameter is not necessary");
		exit(1);
	}

	for (argIndex = 1; argIndex < argc; argIndex++) {
		if (argv[argIndex][0] == '-') {

			if (strcmp(argv[argIndex], "-fastq") == 0)
				filetype = FASTQ;
			else if (strcmp(argv[argIndex], "-fasta") == 0)
				filetype = FASTA;
			else if (strcmp(argv[argIndex], "-gerald") == 0)
				filetype = GERALD;
			else if (strcmp(argv[argIndex], "-eland") == 0)
				filetype = ELAND;
			else if (strcmp(argv[argIndex], "-fastq.gz") == 0)
				filetype = FASTQ_GZ;
			else if (strcmp(argv[argIndex], "-fasta.gz") == 0)
				filetype = FASTA_GZ;
			else if (strcmp(argv[argIndex], "-maq.gz") == 0)
				filetype = MAQ_GZ;
			else if (strcmp(argv[argIndex], "-short") == 0)
				cat = 0;
			else if (strcmp(argv[argIndex], "-shortPaired") ==
				 0)
				cat = 1;
			else if (strncmp
				 (argv[argIndex], "-shortPaired",
				  12) == 0) {
				sscanf(argv[argIndex], "-shortPaired%hi",
				       (short int *) &cat);
				if (cat < 1 || cat > CATEGORIES) {
					printf("Unknown option: %s\n",
					       argv[argIndex]);
					exit(1);
				}
				cat--;
				cat *= 2;
				cat++;
			} else if (strncmp(argv[argIndex], "-short", 6) ==
				   0) {
				sscanf(argv[argIndex], "-short%hi",
				       (short int *) &cat);
				if (cat < 1 || cat > CATEGORIES) {
					printf("Unknown option: %s\n",
					       argv[argIndex]);
					exit(1);
				}
				cat--;
				cat *= 2;
			} else if (strcmp(argv[argIndex], "-long") == 0)
				cat = CATEGORIES * 2;
			else if (strcmp(argv[argIndex], "-longPaired") ==
				 0)
				cat = CATEGORIES * 2 + 1;


			else {
				printf("Unknown option: %s\n",
				       argv[argIndex]);
				exit(0);
			}

			continue;
		}

		if (cat == -1)
			continue;

		switch (filetype) {
		case FASTA:
			reads = readFastAFile(argv[argIndex]);
			break;
		case FASTQ:
			reads = readFastQFile(argv[argIndex]);
			break;
		case GERALD:
			reads = readSolexaFile(argv[argIndex]);
			break;
		case ELAND:
			reads = readElandFile(argv[argIndex]);
			break;
		case FASTA_GZ:
			reads = readFastAGZFile(argv[argIndex]);
			break;
		case FASTQ_GZ:
			reads = readFastQGZFile(argv[argIndex]);
			break;
		case MAQ_GZ:
			reads = readMAQGZFile(argv[argIndex]);
			break;
		default:
			puts("Screw up in parser... exiting");
			exit(0);
		}

		convertSequences(reads);
		categorizeReads(reads, cat);
		fileIndex++;
		concatenateReadSets(allSequences, reads);
	}

	return allSequences;

}

// General argument parser for most functions
// Basically a reused portion of toplevel code dumped into here
ReadSet *parseDataAndReadMaskFiles(int argc, char **argv)
{
	int argIndex = 1;
	ReadSet *reads;
	ReadSet *allSequences = newReadSet();
	int filetype = FASTA;
	boolean cat = false;

	if (argc < 2) {
		destroyReadSet(allSequences);
		return NULL;
	}

	for (argIndex = 1; argIndex < argc; argIndex++) {
		if (argv[argIndex][0] == '-') {

			if (strcmp(argv[argIndex], "-fastq") == 0)
				filetype = FASTQ;
			else if (strcmp(argv[argIndex], "-fasta") == 0)
				filetype = FASTA;
			else if (strcmp(argv[argIndex], "-gerald") == 0)
				filetype = GERALD;
			else if (strcmp(argv[argIndex], "-eland") == 0)
				filetype = ELAND;
			else if (strcmp(argv[argIndex], "-fastq.gz") == 0)
				filetype = FASTQ_GZ;
			else if (strcmp(argv[argIndex], "-fasta.gz") == 0)
				filetype = FASTA_GZ;
			else if (strcmp(argv[argIndex], "-maq.gz") == 0)
				filetype = MAQ_GZ;
			else if (strcmp(argv[argIndex], "-short") == 0)
				cat = false;
			else if (strcmp(argv[argIndex], "-shortPaired") ==
				 0)
				cat = false;
			else if (strcmp(argv[argIndex], "-short2") == 0)
				cat = false;
			else if (strcmp(argv[argIndex], "-shortPaired2") ==
				 0)
				cat = false;
			else if (strcmp(argv[argIndex], "-long") == 0)
				cat = false;
			else {
				printf("Unknown option: %s\n",
				       argv[argIndex]);
				exit(0);
			}

			continue;
		}

		if (!cat)
			continue;

		switch (filetype) {
		case FASTA:
			reads = readFastAFile(argv[argIndex]);
			break;
		case FASTQ:
			reads = readFastQFile(argv[argIndex]);
			break;
		case GERALD:
			reads = readSolexaFile(argv[argIndex]);
			break;
		case ELAND:
			reads = readElandFile(argv[argIndex]);
			break;
		case FASTA_GZ:
			reads = readFastAGZFile(argv[argIndex]);
			break;
		case FASTQ_GZ:
			reads = readFastQGZFile(argv[argIndex]);
			break;
		case MAQ_GZ:
			reads = readMAQGZFile(argv[argIndex]);
			break;
		default:
			puts("Screw up in parser... exiting");
			exit(0);
		}

		convertSequences(reads);
		concatenateReadSets(allSequences, reads);
	}

	return allSequences;

}

void importClippingData(char *filename, ReadSet * reads)
{
	FILE *file = fopen(filename, "r");
	char line[100];
	const int maxline = 5000;
	IDnum index = 0;
	Coordinate start, finish;
	TightString **sequences = reads->tSequences;

	if (file == NULL) {
		printf("Could not read %s, sorry.\n", filename);
		exit(1);
	}

	puts("Importing clip data");

	// For each other lines
	while (fgets(line, maxline, file) != NULL) {
		if (line[0] == 'F') {
			destroyTightString(sequences[index]);
			sequences[index] = NULL;
		} else {
			sscanf(line, "%*[PASFIL ]%*i  %li %li %*[^\n]",
			       &start, &finish);
			clipTightString(sequences[index], start, finish);
		}
		index++;
	}

	puts("Done");

	fclose(file);
}

void pairUpReads(ReadSet * reads, Category cat)
{
	int phase = 0;
	IDnum *mateReads = malloc(reads->readCount * sizeof(IDnum));
	IDnum index;

	if (mateReads == NULL && reads->readCount > 0) {
		puts("Malloc failure");
		exit(1);
	}

	for (index = 0; index < reads->readCount; index++) {
		if (reads->categories[index] != cat) {
			mateReads[index] = -1;
			phase = 0;
		} else if (phase == 0) {
			mateReads[index] = index + 1;
			phase = 1;
		} else {
			mateReads[index] = index - 1;
			phase = 0;
		}
	}

	free(reads->mateReads);
	reads->mateReads = mateReads;
}

void detachDubiousReads(ReadSet * reads, boolean * dubiousReads)
{
	IDnum index;
	IDnum pairID;
	IDnum sequenceCount = reads->readCount;
	IDnum *mateReads = reads->mateReads;

	if (dubiousReads == NULL || mateReads == NULL)
		return;

	for (index = 0; index < sequenceCount; index++) {
		if (!dubiousReads[index])
			continue;

		pairID = mateReads[index];

		if (pairID != -1) {
			//printf("Separating %li and %li\n", index, pairID);
			mateReads[index] = -1;
			mateReads[pairID] = -1;
		}
	}
}

static void exportRead(FILE * outfile, ReadSet * reads, IDnum index)
{
	Coordinate start, finish;
	char str[100];
	TightString *sequence = reads->tSequences[index];

	if (sequence == NULL)
		return;

	fprintf(outfile, ">SEQUENCE_%li_length_%li", index,
		getLength(sequence));

	if (reads->categories != NULL)
		fprintf(outfile, "\t%i", reads->categories[index]);

	fprintf(outfile, "\n");

	start = 0;
	while (start <= getLength(sequence)) {
		finish = start + 60;
		readTightStringFragment(sequence, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}

	fflush(outfile);
}

void exportReadSet(char *filename, ReadSet * reads)
{
	IDnum index;
	FILE *outfile = fopen(filename, "w+");

	if (outfile == NULL) {
		puts("Couldn't open file, sorry");
		return;
	} else
		printf("Writing into readset file: %s\n", filename);

	for (index = 0; index < reads->readCount; index++) {
		exportRead(outfile, reads, index);
	}

	fclose(outfile);

	puts("Done");
}

ReadSet *importReadSet(char *filename)
{
	FILE *file = fopen(filename, "r");
	char *sequence = NULL;
	Coordinate bpCount = 0;
	const int maxline = 100;
	char line[100];
	IDnum sequenceCount, sequenceIndex;
	IDnum index;
	ReadSet *reads;

	if (file != NULL)
		printf("Reading read set file %s;\n", filename);
	else {
		printf("Could not open %s\n", filename);
		exit(0);
	}

	reads = newReadSet();

	// Count number of separate sequences
	sequenceCount = 0;
	while (fgets(line, maxline, file) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	fclose(file);
	printf("%li sequences found\n", sequenceCount);

	reads->readCount = sequenceCount;
	reads->sequences = calloc(sequenceCount, sizeof(char *));
	if (reads->sequences == NULL && sequenceCount > 0) {
		puts("Calloc failure");
		exit(1);
	}
	reads->categories = calloc(sequenceCount, sizeof(Category));
	if (reads->categories == NULL && sequenceCount > 0) {
		puts("Calloc failure");
		exit(1);
	}
	// Counting base pair length of each sequence:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file) != NULL) {
		if (line[0] == '>') {

			// Reading category info
			sscanf(line, "%*[^\t]\t%hhi",
			       &(reads->categories[sequenceIndex + 1]));

			if (sequenceIndex != -1) {
				//printf("Sequence %li has length %li\n",
				//       sequenceIndex, bpCount);
				reads->sequences[sequenceIndex] =
				    malloc(sizeof(char) * (bpCount + 1));
				if (reads->sequences[sequenceIndex] ==
				    NULL) {
					puts("Allocation screwed up!");
					exit(1);
				}
			}
			sequenceIndex++;
			bpCount = 0;
		} else {
			bpCount += (long) strlen(line) - 1;
		}
	}

	//printf("Sequence %li has length %li\n", sequenceIndex, bpCount);
	reads->sequences[sequenceIndex] =
	    malloc(sizeof(char) * (bpCount + 1));
	if (reads->sequences[sequenceIndex] == NULL) {
		puts("Malloc failures");
		exit(1);
	}
	fclose(file);

	// Reopen file and memorize line:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file)) {
		if (line[0] == '>') {
			if (sequenceIndex != -1) {
				sequence[bpCount] = '\0';
			}
			sequenceIndex++;
			bpCount = 0;
			//printf("Starting to read sequence %li\n",
			//       sequenceIndex);
			sequence = reads->sequences[sequenceIndex];
		} else {
			for (index = 0; index < (long) strlen(line) - 1;
			     index++)
				sequence[bpCount + index] = line[index];
			bpCount += (long) (strlen(line) - 1);
		}
	}

	sequence[bpCount] = '\0';
	fclose(file);

	puts("Done");
	return reads;

}

void logInstructions(int argc, char **argv, char *directory)
{
	int index;
	char *logFilename =
	    malloc((strlen(directory) + 100) * sizeof(char));
	FILE *logFile;
	time_t date;
	char *string;

	if (logFilename == NULL) {
		puts("Malloc failure");
		exit(1);
	}
	time(&date);
	string = ctime(&date);

	strcpy(logFilename, directory);
	strcat(logFilename, "/Log");
	logFile = fopen(logFilename, "a");

	if (logFile == NULL) {
		printf("Could not open file %s, exiting...\n",
		       logFilename);
		exit(0);
	}

	fprintf(logFile, "%s", string);

	for (index = 0; index < argc; index++)
		fprintf(logFile, " %s", argv[index]);

	fprintf(logFile, "\n");

	fclose(logFile);
	free(logFilename);
}

void destroyReadSet(ReadSet * reads)
{
	IDnum index;

	if (reads == NULL)
		return;

	if (reads->sequences != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->sequences[index]);

	if (reads->tSequences != NULL)
		for (index = 0; index < reads->readCount; index++)
			destroyTightString(reads->tSequences[index]);

	if (reads->labels != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->labels[index]);

	if (reads->confidenceScores != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->confidenceScores[index]);

	if (reads->kmerProbabilities != NULL)
		for (index = 0; index < reads->readCount; index++)
			free(reads->kmerProbabilities[index]);

	free(reads->sequences);
	free(reads->tSequences);
	free(reads->labels);
	free(reads->confidenceScores);
	free(reads->kmerProbabilities);
	free(reads->mateReads);
	free(reads->categories);
	free(reads);
}

Coordinate *getSequenceLengths(ReadSet * reads, int wordLength)
{
	Coordinate *lengths = calloc(reads->readCount, sizeof(Coordinate));
	IDnum index;
	int lengthOffset = wordLength - 1;

	if (lengths == NULL && reads->readCount > 0) {
		puts("Calloc failure");
		exit(1);
	}

	for (index = 0; index < reads->readCount; index++)
		lengths[index] =
		    getLength(reads->tSequences[index]) - lengthOffset;

	return lengths;
}

Coordinate *getSequenceLengthsFromFile(char *filename, int wordLength)
{
	Coordinate *lengths;
	FILE *file = fopen(filename, "r");
	Coordinate bpCount = 0;
	const int maxline = 100;
	char line[100];
	IDnum sequenceCount, sequenceIndex;
	int lengthOffset = wordLength - 1;

	if (file != NULL)
		printf("Reading read set file %s;\n", filename);
	else {
		printf("Could not open %s\n", filename);
		exit(0);
	}

	// Count number of separate sequences
	sequenceCount = 0;
	while (fgets(line, maxline, file) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	fclose(file);

	lengths = calloc(sequenceCount, sizeof(Coordinate));
	if (lengths == NULL && sequenceCount > 0) {
		puts("Calloc failure");
		exit(1);
	}
	// Counting base pair length of each sequence:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file) != NULL) {
		if (line[0] == '>') {
			if (sequenceIndex != -1)
				lengths[sequenceIndex] =
				    bpCount - lengthOffset;
			sequenceIndex++;
			bpCount = 0;
		} else {
			bpCount += (long) strlen(line) - 1;
		}
	}
	fclose(file);

	return lengths;
}
