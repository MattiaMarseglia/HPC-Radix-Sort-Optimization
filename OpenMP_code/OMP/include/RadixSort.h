/*
 * Course: High Performance Computing 2021/2022
 *
 * Lecturer: Francesco Moscato	fmoscato@unisa.it
 *
 * Group:
 * Marseglia	Mattia		0622701697	    m.marseglia1@studenti.unisa.it
 * Spingola     Camilla		0622701698  	c.spingola@studenti.unisa.it
 * Turi		    Vito		0622701795  	v.turi3@studenti.unisa.it
 * Sica 		Ferdinando	0622701794	    f.sica24@studenti.unisa.it
 *
 * Copyright (C) 2021 - All Rights Reserved
 *
 * This file is part of Contest-OMP: RadixSort.
 *
 * Contest-OMP: RadixSort is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Contest-OMP: RadixSort is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Contest-OMP: RadixSort.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RADIXSORT_H_ /* Include guard */
#define RADIXSORT_H_

/** macros to get execution time: both macros have to be in the same scope
 *   define a double variable to use in ENDTIME before STARTTIME:
 *   double x;
 *   the variable will hold the execution time in seconds.
 */

#include <time.h>

/* Token concatenation used */
#define STARTTIME(id)                             \
    clock_t start_time_42_##id, end_time_42_##id; \
    start_time_42_##id = clock()

#define ENDTIME(id, x)          \
    end_time_42_##id = clock(); \
    x = ((double)(end_time_42_##id - start_time_42_##id)) / CLOCKS_PER_SEC

// element type, length, thread
void getMaxDigit(TYPE_OF_ELEMENTS array[], int n, TYPE_OF_ELEMENTS *array_pos, TYPE_OF_ELEMENTS *array_neg, int *max_pos, int *max_neg, int *pos, int *neg, int threads);
TYPE_OF_ELEMENTS getMaxDigit2(int n, TYPE_OF_ELEMENTS *arr, int threads);
void countingSort(TYPE_OF_ELEMENTS array[], int base, int size, int raw_index, TYPE_OF_ELEMENTS *matrix, int threads);
void brutalSort(int n, TYPE_OF_ELEMENTS *arr, int exp, int threads);
void serviceRadixsort(TYPE_OF_ELEMENTS array[], int size, int max, int threads);
void myRadixsort(TYPE_OF_ELEMENTS *array, int size, int threads);
void myRadixsort2(int n, TYPE_OF_ELEMENTS *arr, int threads);
void init_structures(TYPE_OF_ELEMENTS **array, int size, int max_digit);

#endif
