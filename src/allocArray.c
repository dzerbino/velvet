/*
Copyright 2009 Sylvain Foret (sylvain.foret@anu.edu.au) 

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
#include <unistd.h>
#ifdef OPENMP
#include <omp.h>
#endif

#include "allocArray.h"
#include "utility.h"

#if DEBUG
#define NB_PAGES_ALLOC    1
#define BLOCKS_ALLOC_SIZE 1
#else
#define NB_PAGES_ALLOC    128
#define BLOCKS_ALLOC_SIZE 128
#endif


struct AllocArrayFreeElement_st {
	AllocArrayFreeElement *next;
	ArrayIdx idx;
};

static void initAllocArray(AllocArray * array, size_t elementSize, char * name) 
{
	array->freeElements = NULL;
	array->elementSize = elementSize;
	array->blockSize = sysconf (_SC_PAGESIZE) * NB_PAGES_ALLOC;
	array->maxElements = array->blockSize / array->elementSize;
	array->maxBlocks = BLOCKS_ALLOC_SIZE;
	array->blocks = mallocOrExit (array->maxBlocks, void*);
	array->blocks[0] = mallocOrExit (array->blockSize, char);
	array->currentBlocks = 1;
	array->currentElements = 0;
#if DEBUG
	array->elementsRecycled = 0;
	array->elementsAllocated = 0;
	array->name = name;
#endif
#ifdef OPENMP
	array->counter = 0;
	array->shut = false;
#endif
}

AllocArray*
newAllocArray (size_t elementSize, char *name)
{
	AllocArray *array;

	if (elementSize < sizeof(AllocArrayFreeElement)) {
		velvetLog("Elements too small to create an AllocArray!\n");
		exit(-1);
	}
	array = mallocOrExit (1, AllocArray);
	initAllocArray(array, elementSize, name);
	return array;
}

void
destroyAllocArrayChunks (AllocArray *array)
{
	size_t i;
	for (i = 0; i < array->currentBlocks; i++)
		free (array->blocks[i]);
	free (array->blocks);
}

void
destroyAllocArray (AllocArray *array)
{
	if (array)
	{
		destroyAllocArrayChunks(array);
#if DEBUG
		velvetLog(">>> Allocation summary for %s\n", array->name);
		velvetLog(">>> Alloc'ed %ld bytes\n", array->blockSize * array->currentBlocks);
		velvetLog(">>> Alloc'ed %ld elements\n", array->elementsAllocated);
		velvetLog(">>> Recycled %ld elements\n", array->elementsRecycled);
#endif
		free (array);
	}
}

ArrayIdx
allocArrayAllocate (AllocArray *array)
{
	if (array->freeElements != NULL)
	{
		AllocArrayFreeElement *element;

		element = array->freeElements;
		array->freeElements = element->next;
#if DEBUG
		array->elementsRecycled++;
#endif
		return element->idx;
	}
	if (array->currentElements >= array->maxElements)
	{
		if (array->currentBlocks == array->maxBlocks)
		{
			array->maxBlocks += BLOCKS_ALLOC_SIZE;
#ifdef OPENMP
			array->shut = true;
			#pragma omp flush (array)
			while(array->counter > 0) {
				#pragma omp flush (array)
			}
#endif
			array->blocks = reallocOrExit (array->blocks, array->maxBlocks, void*);
#ifdef OPENMP
			#pragma omp flush (array)
			array->shut = false;
			#pragma omp flush (array)
#endif
		}
		array->blocks[array->currentBlocks] = mallocOrExit (array->blockSize, char);
		array->currentBlocks++;
		array->currentElements = 0;
	}
	array->currentElements++;
#if DEBUG
	if (array->maxElements * (array->currentBlocks - 1) + array->currentElements == UINT32_MAX)
	{
		velvetLog (">>> Reached maximum `%s' addressable with 32 bits\n", array->name);
		abort();
	}
	array->elementsAllocated++;
#endif
	return array->maxElements * (array->currentBlocks - 1) + array->currentElements;
}

void
allocArrayFree (AllocArray *array, ArrayIdx idx)
{
	if (idx != NULL_IDX)
	{
		AllocArrayFreeElement *freeElem;

		freeElem = allocArrayGetElement (array, idx);
		freeElem->idx = idx;
		freeElem->next = array->freeElements;
		array->freeElements = freeElem;
	}
}

void*
allocArrayGetElement (AllocArray *array, ArrayIdx idx)
{
	if (idx != NULL_IDX)
	{
		const ArrayIdx i = idx - 1;
		const ArrayIdx blockIdx = i / array->maxElements;
		const ArrayIdx elementIdx = i % array->maxElements;
		return (char *) (array->blocks[blockIdx]) + elementIdx * array->elementSize;
	}
	return NULL;
}

#ifdef OPENMP
AllocArray *newAllocArrayArray(unsigned int n,
			       size_t node_size, char * name)
{
	AllocArray *allocArray;
	int i;

	allocArray = mallocOrExit (n + 1, AllocArray);
	for (i = 0; i < n; i++)
		initAllocArray(allocArray + i, node_size, name);
	/* Last element marker */
	allocArray[n].elementSize = 0;

	return allocArray;
}

void destroyAllocArrayArray(AllocArray * allocArray)
{
	int i;

	for (i = 0; allocArray[i].elementSize != 0; i++)
		destroyAllocArrayChunks(allocArray + i);
	free(allocArray);
}

AllocArray *getAllocArrayInArray(AllocArray *allocArray,
				 int	     position)
{
	return allocArray + position;
}

ArrayIdx allocArrayArrayAllocate (AllocArray *array)
{
	int thread = omp_get_thread_num();
	size_t idx = (size_t) allocArrayAllocate(array + thread);

	if (idx == NULL_IDX)
		return NULL_IDX;

	idx *= omp_get_max_threads();
	idx += thread;

	if (idx > UINT32_MAX)
	{
		velvetLog (">>> Reached maximum addressable with 32 bits\n");
#ifdef  DEBUG
		abort();
#else
		exit(1);
#endif
	}

	if (idx < omp_get_max_threads())
		abort();

	return (ArrayIdx) idx;
}

void
allocArrayArrayFree (AllocArray *array, ArrayIdx idx)
{
	if (idx != NULL_IDX)
	{
		array += idx % omp_get_max_threads();
		idx /= omp_get_max_threads(); 
		AllocArrayFreeElement *freeElem;

		freeElem = allocArrayGetElement (array, idx);
		freeElem->idx = idx;
		freeElem->next = array->freeElements;
		array->freeElements = freeElem;
	}
}
#endif
