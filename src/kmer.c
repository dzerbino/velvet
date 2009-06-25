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
#include "kmer.h"
#include "utility.h"

static const unsigned long long int longLongLeftFilter = (long long int) 3 << 62; 
static const unsigned long int longLeftFilter = (long int) 3 << 30; 
static const unsigned int intLeftFilter = (int) 3 << 14; 
static const unsigned char charLeftFilter = (char) 3 << 6; 

static unsigned long long int longLongWordFilter = -1; 
static unsigned long int longWordFilter = (long) (((long long int) 1) << 32) - 1; 
static unsigned int intWordFilter = (int) (((long int) 1) << 16) - 1; 
static unsigned char charWordFilter = (char) (((int) 1) << 8) - 1; 

static unsigned int longLongKmerFilterIndex = KMER_LONGLONGS;
static unsigned long long int longLongKmerFilter = -1;

void resetWordFilter(int wordLength) {
	int kmer_bit_size = wordLength * 2;
	int i;

	if (wordLength > MAXKMERLENGTH) 
		exitErrorf(EXIT_FAILURE, true, "Word length %i greater than max allowed value (%i).\nRecompile Velvet to deal with this word length.", wordLength, MAXKMERLENGTH);

#if KMER_LONGLONGS
	for (i = 0; i < KMER_LONGLONGS; i++) {
		if (kmer_bit_size > 64) {
			kmer_bit_size -= 64;
			continue;
		} else if (kmer_bit_size == 64) {
			longLongKmerFilterIndex = i;
			longLongKmerFilter = longLongWordFilter; 
			longWordFilter = 0;
			intWordFilter = 0;
			charWordFilter = 0;
		} else {
			longLongKmerFilterIndex = i;
			longLongKmerFilter = (((long long int) 1) << kmer_bit_size) - 1;	
			longWordFilter = 0;
			intWordFilter = 0;
			charWordFilter = 0;
			return;
		}
	}
#endif
#if KMER_LONGS
	if (kmer_bit_size > 32) {
		kmer_bit_size -= 32;
		;
	} else if (kmer_bit_size == 32) {
		intWordFilter = 0;
		charWordFilter = 0;
	} else {
		longWordFilter = (((long int) 1) << kmer_bit_size) - 1;	
		intWordFilter = 0;
		charWordFilter = 0;
		return;
	}
#endif
#if KMER_INTS
	if (kmer_bit_size > 16) {
		kmer_bit_size -= 16;
		;
	} else if (kmer_bit_size == 16) {
		charWordFilter = 0;
	} else {
		intWordFilter = (((int) 1) << kmer_bit_size) - 1;	
		charWordFilter = 0;
		return;
	}

#endif
#if KMER_CHARS
	if (kmer_bit_size >= 8) {
	} else {
		charWordFilter = (((char) 1) << kmer_bit_size) - 1;	
		return;
	}

#endif

}

static void shiftLeft(Kmer * kmer) {
	int i;
	unsigned long long int leftBits; 
	unsigned long long int rightBits = 0;

#if KMER_LONGLONGS
	for (i = 0; i < KMER_LONGLONGS; i++) {
		leftBits = (kmer->longlongs[i] & longLongLeftFilter);
		leftBits >>= 62;
		kmer->longlongs[i] <<= 2;
		kmer->longlongs[i] += rightBits;
		if (i < longLongKmerFilterIndex)
			kmer->longlongs[i] &= longLongWordFilter;
		else if (i == longLongKmerFilterIndex)
			kmer->longlongs[i] &= longLongKmerFilter;
		else 
			kmer->longlongs[i] = 0;
		rightBits = leftBits;
	}
#endif
#if KMER_LONGS
	leftBits = kmer->longs & longLeftFilter;
	leftBits >>= 30;
	kmer->longs <<= 2;
	kmer->longs += rightBits;
	kmer->longs &= longWordFilter;
	rightBits = leftBits;
#endif
#if KMER_INTS
	leftBits = kmer->ints & intLeftFilter;
	leftBits >>= 14;
	kmer->ints <<= 2;
	kmer->ints += rightBits;
	kmer->ints &= intWordFilter;
	rightBits = leftBits;
#endif
#if KMER_CHARS
	leftBits = kmer->chars & charLeftFilter;
	leftBits >>= 6;
	kmer->chars <<= 2;
	kmer->chars += rightBits;
	kmer->chars &= charWordFilter;
	rightBits = leftBits;
#endif
}

static void shiftRight(Kmer * kmer) {
	int i;
	unsigned long long int leftBits = 0;
	unsigned long long int rightBits;

#if KMER_CHARS
	rightBits = kmer->chars & 3;
	leftBits <<= 6;
	kmer->chars >>= 2;
	kmer->chars += (unsigned char) leftBits;
	leftBits = rightBits;
#endif
#if KMER_INTS
	rightBits = kmer->ints & 3;
	leftBits <<= 14;
	kmer->ints >>= 2;
	kmer->ints += (unsigned int) leftBits;
	leftBits = rightBits;
#endif
#if KMER_LONGS
	rightBits = kmer->longs & 3;
	leftBits <<= 30;
	kmer->longs >>= 2;
	kmer->longs += (unsigned long int) leftBits;
	rightBits = leftBits;
#endif
#if KMER_LONGLONGS
	for (i = KMER_LONGLONGS - 1; i >= 0; i--) {
		rightBits = kmer->longlongs[i] & 3;
		leftBits <<= 62;
		kmer->longlongs[i] >>= 2;
		kmer->longlongs[i] += leftBits;
		leftBits = rightBits;
	}
#endif
}

void copyKmers(Kmer* k1, Kmer* k2) {
	int i;

#if KMER_LONGLONGS
	for (i = 0; i < KMER_LONGLONGS; i++)
		k1->longlongs[i] = k2->longlongs[i];
#endif
#if KMER_LONGS
	k1->longs = k2->longs;
#endif
#if KMER_INTS
	k1->ints = k2->ints;
#endif
#if KMER_CHARS
	k1->chars = k2->chars;
#endif
}

static void incrementKmer(Kmer* kmer, char increment) {
	int i;

#if KMER_LONGLONGS
	kmer->longlongs[0] += increment;
	if (kmer->longlongs[0] >= increment)
		return;

	for (i = 1; i < KMER_LONGLONGS; i++) 
		if (++kmer->longlongs[i])
			return;
#if KMER_LONGS
	if (++kmer->longs)
		return;
#endif
#if KMER_INTS
	if (++kmer->ints)
		return;
#endif
#if KMER_CHARS
	++kmer->chars;
#endif

#else

#if KMER_LONGS
	kmer->longs += increment;
	if (kmer->longs >= increment)
		return;
#if KMER_INTS
	if (++kmer->ints)
		return;
#endif
#if KMER_CHARS
	++kmer->chars;
#endif

#else

#if KMER_INTS
	kmer->ints += increment;
	if (kmer->ints >= increment)
		return;
#if KMER_CHARS
	++kmer->chars;
#endif

#else 

#if KMER_CHARS
	kmer->chars += increment;
#endif

#endif
#endif
#endif
}

int compareKmers(Kmer* k1, Kmer* k2) {
	int i;

#if KMER_CHARS
	if (k1->chars == k2->chars)
		;
	else if (k1->chars > k2->chars)
		return 1;
	else 
		return -1;
#endif
#if KMER_INTS
	if (k1->ints == k2->ints)
		;
	else if (k1->ints > k2->ints)
		return 1;
	else 
		return -1;
#endif
#if KMER_LONGS
	if (k1->longs == k2->longs)
		;
	else if (k1->longs > k2->longs)
		return 1;
	else 
		return -1;
#endif
#if KMER_LONGLONGS
	for (i = KMER_LONGLONGS - 1; i >= 0; i--) {
		if (k1->longlongs[i] == k2->longlongs[i])
			continue;
		else if (k1->longlongs[i] > k2->longlongs[i])
			return 1;
		else 
			return -1;
	}
#endif

	return 0;
}

void clearKmer(Kmer * kmer) {
	int i;

#if KMER_LONGLONGS
	for (i = 0; i < KMER_LONGLONGS; i++)
		kmer->longlongs[i] = 0;
#endif
#if KMER_LONGS
	kmer->longs = 0;
#endif
#if KMER_INTS
	kmer->ints = 0;
#endif
#if KMER_CHARS
	kmer->chars = 0;
#endif
}

unsigned char rightMostNucleotide(Kmer * kmer) {
#if KMER_LONGLONGS
	return (unsigned char) (kmer->longlongs[0] && 3);
#endif
#if KMER_LONGS
	return (unsigned char) (kmer->longs && 3);
#endif
#if KMER_INTS
	return (unsigned char) (kmer->ints && 3);
#endif
#if KMER_CHARS
	return (unsigned char) (kmer->chars && 3);
#endif
}

void reverseComplement(Kmer* revComp, Kmer * word, int WORDLENGTH)
{
	int index;
	Kmer copy;
	Nucleotide nucleotide;
	
	clearKmer(revComp);
	copyKmers(&copy, word);

	for (index = 0; index < WORDLENGTH; index++) {
		nucleotide = popNucleotide(&copy);
#ifndef COLOR
		pushNucleotide(revComp, 3 - nucleotide);
#else
		pushNucleotide(revComp, nucleotide);
#endif
	}
}

void printKmer(Kmer * kmer) {
	int i;

#if KMER_CHARS
	printf("%hx\t", kmer->chars);
#endif
#if KMER_INTS
	printf("%x\t", kmer->ints);
#endif
#if KMER_LONGS
	printf("%lx\t", kmer->longs);
#endif
#if KMER_LONGLONGS
	for (i = KMER_LONGLONGS - 1; i >= 0; i--)
		printf("%llx\t", kmer->longlongs[i]);
#endif
	puts("");
}

void testKmers(int argc, char** argv) {
	Kmer kmer;		
	Kmer *k2;
	Kmer k3;
	Kmer k4;

	k2 = &k4;
	int i;
	
	printf("FORMATS %u %u %u %u\n", KMER_CHARS, KMER_INTS, KMER_LONGS, KMER_LONGLONGS);
	printf("FILTERS %x %x %lx %llx\n", (int) charLeftFilter, intLeftFilter, longLeftFilter, longLongLeftFilter);
	printf("FILTERS %hx %x %lx %llx\n", charWordFilter, intWordFilter, longWordFilter, longLongWordFilter);
	printKmer(&kmer);
	puts("Clear");
	clearKmer(&kmer);
	printKmer(&kmer);

	puts("Fill up");
	for (i = 0; i < MAXKMERLENGTH; i++) {
		pushNucleotide(&kmer, i % 4);
		printKmer(&kmer);
	}

	puts("Shift right");
	for (i = 0; i < MAXKMERLENGTH; i++) {
		popNucleotide(&kmer);
		printKmer(&kmer);
	}

	puts("Copy");
	copyKmers(k2, &kmer);
	printKmer(k2); 	
	printf("%i\n", compareKmers(k2, &kmer)); 

	puts("Reverse complement");
	clearKmer(k2);
	pushNucleotide(k2, 1);
	reverseComplement(&k3, k2, MAXKMERLENGTH);
	printKmer(&k3); 	
	printf("%i\n", compareKmers(k2, &kmer)); 
}

static Nucleotide getRightNucleotide(Kmer * kmer) {
#if KMER_LONGLONGS
	return kmer->longlongs[0] & 3;
#endif
#if KMER_LONGS
	return kmer->longs & 3;
#endif
#if KMER_INTS
	return kmer->ints & 3;
#endif
#if KMER_CHARS
	return kmer->chars & 3;
#endif
	
}

void pushNucleotide(Kmer * kmer, Nucleotide nucleotide) {
	shiftLeft(kmer);
	incrementKmer(kmer, nucleotide);
}

Nucleotide popNucleotide(Kmer * kmer) {
	Nucleotide nucl = getRightNucleotide(kmer);
	shiftRight(kmer);
	return nucl;
}
