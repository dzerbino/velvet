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
#include "recycleBin.h"

#define CHUNKSIZE 10000

static RecycleBin *treeMemory = NULL;

struct splayNode_st {
	Kmer kmer;
	IDnum seqID;
	Coordinate position;

	struct splayNode_st *left;
	struct splayNode_st *right;
};

typedef struct splayNode_st SplayNode;
typedef struct splayNode_st SplayTree;

static SplayNode *allocateSplayNode()
{
	if (treeMemory == NULL)
		treeMemory = newRecycleBin(sizeof(SplayNode), CHUNKSIZE);

	return (SplayNode *) allocatePointer(treeMemory);
}

static void deallocateSplayNode(SplayNode * node)
{
	deallocatePointer(treeMemory, node);
}

SplayTree *newSplayTree()
{
	return NULL;
}

void destroySplayTree(SplayTree * T)
{
	if (T == NULL)
		return;

	destroySplayTree(T->left);
	destroySplayTree(T->right);
	deallocateSplayNode(T);
}

void destroyAllSplayTrees()
{
	destroyRecycleBin(treeMemory);
	treeMemory = NULL;
}

/* This function can be called only if K2 has a left child */
/* Perform a rotate between a node (K2) and its left child */
/* Update heights, then return new root */

static SplayNode *SingleRotateWithLeft(SplayNode * K2)
{
	SplayNode *K1;

	K1 = K2->left;
	K2->left = K1->right;
	K1->right = K2;

	return K1;		/* New root */
}

/* This function can be called only if K1 has a right child */
/* Perform a rotate between a node (K1) and its right child */
/* Update heights, then return new root */

static SplayNode *SingleRotateWithRight(SplayNode * K1)
{
	SplayNode *K2;

	K2 = K1->right;
	K1->right = K2->left;
	K2->left = K1;

	return K2;		/* New root */
}

/* Top-down splay procedure, */
/* not requiring kmer to be in tree */

static SplayTree *Splay(Kmer kmer, SplayTree * T)
{
	SplayNode Header;
	SplayNode *LeftTreeMax, *RightTreeMin;

	if (T == NULL)
		return NULL;

	Header.left = Header.right = NULL;
	LeftTreeMax = RightTreeMin = &Header;

	while (kmer != T->kmer) {
		if (kmer < T->kmer) {
			if (T->left == NULL)
				break;
			if (kmer < T->left->kmer)
				T = SingleRotateWithLeft(T);
			if (T->left == NULL)
				break;
			/* Link right */
			RightTreeMin->left = T;
			RightTreeMin = T;
			T = T->left;
		} else {
			if (T->right == NULL)
				break;
			if (kmer > T->right->kmer)
				T = SingleRotateWithRight(T);
			if (T->right == NULL)
				break;
			/* Link left */
			LeftTreeMax->right = T;
			LeftTreeMax = T;
			T = T->right;
		}
	}			/* while kmer != T->kmer */

	/* Reassemble */
	LeftTreeMax->right = T->left;
	RightTreeMin->left = T->right;
	T->left = Header.right;
	T->right = Header.left;

	return T;
}

Kmer findInTree(Kmer X, SplayTree ** T)
{
	*T = Splay(X, *T);
	return (*T)->kmer;
}

void insertIntoTree(Kmer kmer, SplayTree ** T)
{
	SplayNode *newNode;

	if (*T == NULL) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->left = newNode->right = NULL;
		*T = newNode;
		return;
	}

	*T = Splay(kmer, *T);
	if (kmer < (*T)->kmer) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->left = (*T)->left;
		newNode->right = *T;
		(*T)->left = NULL;
		*T = newNode;
	} else if ((*T)->kmer < kmer) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->right = (*T)->right;
		newNode->left = *T;
		(*T)->right = NULL;
		*T = newNode;
	}
}

boolean
findOrInsertOccurenceInSplayTree(Kmer kmer, IDnum * seqID,
				 Coordinate * position, SplayTree ** T)
{
	SplayNode *newNode;

	if (*T == NULL) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->left = newNode->right = NULL;

		*T = newNode;

		return false;
	}

	*T = Splay(kmer, *T);
	if (kmer < (*T)->kmer) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->left = (*T)->left;
		newNode->right = *T;
		(*T)->left = NULL;

		*T = newNode;

		return false;
	} else if ((*T)->kmer < kmer) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->right = (*T)->right;
		newNode->left = *T;
		(*T)->right = NULL;

		*T = newNode;

		return false;
	} else {
		*seqID = (*T)->seqID;
		*position = (*T)->position;

		return true;
	}
}

boolean
placeOccurenceInSplayTree(Kmer kmer, IDnum * seqID,
			  Coordinate * position, SplayTree ** T)
{
	SplayNode *newNode;
	IDnum newID = *seqID;
	Coordinate newCoord = *position;

	if (*T == NULL) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->left = newNode->right = NULL;

		*T = newNode;

		return false;
	}

	*T = Splay(kmer, *T);
	if (kmer < (*T)->kmer) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->left = (*T)->left;
		newNode->right = *T;
		(*T)->left = NULL;

		*T = newNode;

		return false;
	} else if ((*T)->kmer < kmer) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->seqID = *seqID;
		newNode->position = *position;

		newNode->right = (*T)->right;
		newNode->left = *T;
		(*T)->right = NULL;

		*T = newNode;

		return false;
	} else {
		*seqID = (*T)->seqID;
		*position = (*T)->position;

		(*T)->seqID = newID;
		(*T)->position = newCoord;

		return true;
	}
}

void printTree(SplayTree * T)
{
	if (T == NULL)
		return;

	printTree(T->left);
	printf("%lli\n", T->kmer);
	printTree(T->right);
}



int test(int argc, char **argv)
{
	int i;
	SplayTree *T = newSplayTree();
	puts("Hello, world");

	puts("---TREE---");
	printTree(T);
	insertIntoTree(1, &T);
	puts("---TREE---");
	printTree(T);
	insertIntoTree(3, &T);
	puts("---TREE---");
	printTree(T);
	insertIntoTree(13, &T);
	puts("---TREE---");
	printTree(T);
	insertIntoTree(5, &T);
	puts("---TREE---");
	printTree(T);
	insertIntoTree(7, &T);
	puts("---TREE---");
	printTree(T);
	insertIntoTree(2, &T);
	puts("---TREE---");
	printTree(T);

	for (i = 1; i < 15; i++)
		printf("Test with %i: %lli\n", i, findInTree(i, &T));

	destroySplayTree(T);

	return 1;
}


void countOccurenceInSplayTree(Kmer kmer, SplayTree ** T, int increment)
{
	SplayNode *newNode;

	if (*T == NULL) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->position = increment;

		newNode->left = newNode->right = NULL;

		*T = newNode;

		return;
	}

	*T = Splay(kmer, *T);
	if (kmer < (*T)->kmer) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->position = increment;

		newNode->left = (*T)->left;
		newNode->right = *T;
		(*T)->left = NULL;

		*T = newNode;
	} else if ((*T)->kmer < kmer) {
		newNode = allocateSplayNode();
		newNode->kmer = kmer;
		newNode->position = increment;

		newNode->right = (*T)->right;
		newNode->left = *T;
		(*T)->right = NULL;

		*T = newNode;
	} else {
		((*T)->position) += increment;
	}

}

void filterAndExportSplayTree(FILE * file, SplayTree * T, int minCov,
			      int maxCov)
{
	if (T == NULL)
		return;

	filterAndExportSplayTree(file, T->left, minCov, maxCov);
	filterAndExportSplayTree(file, T->right, minCov, maxCov);
	if ((minCov == -1 || T->position >= minCov)
	    && (maxCov == -1 || T->position <= maxCov))
		fprintf(file, "%lli\n", T->kmer);
}

void displaySplayTreeMemory()
{
	printf("TREE MEMORY %li allocated %li free\n",
	       RecycleBin_memory_usage(treeMemory),
	       recycleBinFreeSpace(treeMemory));
}
