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
#include <string.h>
#include <sys/stat.h>
#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#include <uce-dirent.h>
#define Arc v_Arc
#else
#include <dirent.h>
#endif

#include "run.h"

static void printUsage()
{
	puts("Usage:");
	puts("./velveth directory hash_length {[-file_format][-read_type] filename1 [filename2 ...]} {...} [options]");
	puts("");
	puts("\tdirectory\t: directory name for output files");
	printf("\thash_length\t: odd integer (if even, it will be decremented) <= %i (if above, will be reduced)\n", MAXKMERLENGTH);
	puts("\tfilename\t: path to sequence file or - for standard input");	
	puts("");
	puts("File format options:");
	puts("\t-fasta\t-fastq\t-raw\t-fasta.gz\t-fastq.gz\t-raw.gz\t-sam\t-bam");
	puts("");
	puts("Read type options:");
	puts("\t-short\t-shortPaired");
#if CATEGORIES <= 5
	Category cat; 
	for (cat = 2; cat <= CATEGORIES; cat++)
	    printf("\t-short%i\t-shortPaired%i\n", cat, cat);
#else
	puts("\t...");
	printf("\t-short%i\t-shortPaired%i\n", CATEGORIES - 1, CATEGORIES - 1);
	printf("\t-short%i\t-shortPaired%i\n", CATEGORIES, CATEGORIES);
#endif
	puts("\t-long\t-longPaired");
	puts("");
	puts("Options:");
	puts("\t-strand_specific\t: for strand specific transcriptome sequencing data (default: off)");
	puts("");
	puts("Synopsis:");
	puts("");
	puts("- Short single end reads:");
	puts("\tvelveth Assem 29 -short -fastq s_1_sequence.txt");
	puts("");
	puts("- Paired-end short reads (remember to interleave paired reads):");
	puts("\tvelveth Assem 31 -shortPaired -fasta interleaved.fna");
	puts("");
	puts("- Two channels and some long reads:");
	puts("\tvelveth Assem 43 -short -fastq unmapped.fna -longPaired -fasta SangerReads.fasta");
	puts("");
	puts("- Three channels:");
	puts("\tvelveth Assem 35 -shortPaired -fasta pe_lib1.fasta -shortPaired2 pe_lib2.fasta -short3 se_lib1.fa");
	puts("");
	puts("Output:");
	puts("\tdirectory/Roadmaps");
	puts("\tdirectory/Sequences");
	puts("\t\t[Both files are picked up by graph, so please leave them there]");
}

int main(int argc, char **argv)
{
	ReadSet *allSequences;
	SplayTable *splayTable;
	int hashLength;
	char *directory, *filename, *seqFilename, *buf;
	boolean double_strand = true;
	DIR *dir;

	setProgramName("velveth");

	if (argc < 4) {
		puts("velveth - simple hashing program");
		printf("Version %i.%i.%2.2i\n", VERSION_NUMBER,
		       RELEASE_NUMBER, UPDATE_NUMBER);
		puts("\nCopyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)");
		puts("This is free software; see the source for copying conditions.  There is NO");
		puts("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
		puts("Compilation settings:");
		printf("CATEGORIES = %i\n", CATEGORIES);
		printf("MAXKMERLENGTH = %i\n", MAXKMERLENGTH);
		puts("");
		printUsage();
		return 0;
	}

	directory = argv[1];
	filename = mallocOrExit(strlen(directory) + 100, char);
	seqFilename = mallocOrExit(strlen(directory) + 100, char);
	buf = mallocOrExit(strlen(directory) + 100, char);

	hashLength = atoi(argv[2]);

	if (hashLength > MAXKMERLENGTH) {
		printf
		    ("Velvet can't handle k-mers as long as %i! We'll stick to %i if you don't mind.\n",
		     hashLength, MAXKMERLENGTH);
		hashLength = MAXKMERLENGTH;
	} else if (hashLength <= 0) {
		printf("Invalid hash length: %s\n", argv[2]);
		printUsage();
		return 0;
	} 

	if (hashLength % 2 == 0) {
		printf
		    ("Velvet can't work with even length k-mers, such as %i. We'll use %i instead, if you don't mind.\n",
		     hashLength, hashLength - 1);
		hashLength--;
	}
	resetWordFilter(hashLength);

	dir = opendir(directory);

	if (dir == NULL)
		mkdir(directory, 0777);
	else {
		sprintf(buf, "%s/PreGraph", directory);
		remove(buf);
		sprintf(buf, "%s/Graph", directory);
		remove(buf);
		sprintf(buf, "%s/Graph2", directory);
		remove(buf);
		sprintf(buf, "%s/Graph3", directory);
		remove(buf);
		sprintf(buf, "%s/Graph4", directory);
		remove(buf);
		sprintf(buf, "%s/Log", directory);
		remove(buf);
	}

	logInstructions(argc, argv, directory);

	strcpy(seqFilename, directory);
	strcat(seqFilename, "/Sequences");
	parseDataAndReadFiles(seqFilename, argc - 2, &(argv[2]), &double_strand);

	splayTable = newSplayTable(hashLength, double_strand);

	allSequences = importReadSet(seqFilename);
	printf("%i sequences in total.\n", allSequences->readCount);

	strcpy(filename, directory);
	strcat(filename, "/Roadmaps");
	inputSequenceArrayIntoSplayTableAndArchive(allSequences,
						   splayTable, filename, seqFilename);

	destroySplayTable(splayTable);
	if (dir)
		closedir(dir);
	free(filename);
	free(seqFilename);
	free(buf);

	return 0;
}
