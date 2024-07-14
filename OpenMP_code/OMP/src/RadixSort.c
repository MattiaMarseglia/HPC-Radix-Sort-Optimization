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

/**
        @file RadixSort.c
*/

#include "RadixSort.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_team_num() -1
#endif

/**
 * @brief This function calculates the number of positive and negative elements in the array and puts it in the corrisponding arrays. It also calculates the maximum positive element and the minimum negative element.
 * @param array         pointer to the array with total elements to sort.
 * @param size          dimension of the array containing all elements.
 * @param array_pos     pointer to array used to store the positive elements.
 * @param array_neg     pointer to array used to store the negative elements.
 * @param max_pos       pointer to the number of digits of the maximum positive element.
 * @param max_neg       pointer to the number of digits of the minimum negative element.
 * @param pos           pointer to the number of positive elements.
 * @param neg           pointer to the number of negative elements.
 * @param threads       number of threads.
 */
void getMaxDigit(TYPE_OF_ELEMENTS array[], int size, TYPE_OF_ELEMENTS *array_pos, TYPE_OF_ELEMENTS *array_neg, int *max_pos, int *max_neg, int *pos, int *neg, int threads) {
    *max_pos = 0;
    *max_neg = 0;
    *neg = 0;
    *pos = 0;
    TYPE_OF_ELEMENTS tmp_neg = 0, tmp_pos = 0, i = 0;
    for (i = 0; i < size; i++) {
        if (array[i] < 0) {
            array_neg[*neg] = array[i];
            (*neg)++;
        }
        if (array[i] >= 0) {
            array_pos[*pos] = array[i];
            (*pos)++;
        }
    }

#pragma parallel for shared(array, size) private(i) num_threads(threads) reduction(min \                    \
                                                                                   : tmp_neg) reduction(max \
                                                                                                        : tmp_pos) schedule(runtime)
    for (i = 0; i < size; i++) {
        if (array[i] < tmp_neg)
            tmp_neg = array[i];
        if (array[i] > tmp_pos)
            tmp_pos = array[i];
    }

    while (tmp_neg < 0) {
        tmp_neg /= 10;
        (*max_neg)++;
    }
    while (tmp_pos > 0) {
        tmp_pos /= 10;
        (*max_pos)++;
    }
}

/**
 * @brief This function implements the counting sort algorithm based on significant places
 * @param array         pointer to the array with the elements to sort, that are only positive or negative.
 * @param base          the base of the elements analized.
 * @param size          dimension of the array to order.
 * @param raw_index     the index in the matrix to which frequences of digits need to be stored.
 * @param matrix        the pointer to the matrix in which for each row are stored the frequences of each positional element
 * @param threads       number of threads.
 */
void countingSort(TYPE_OF_ELEMENTS array[], int base, int size, int raw_index, TYPE_OF_ELEMENTS *matrix, int threads) {
    TYPE_OF_ELEMENTS place = 1;

    // printf("Ciao sono %d thread\n", omp_get_thread_num());
    for (int i = 0; i < raw_index; i++)
        place = place * 10;

    TYPE_OF_ELEMENTS max = (array[0] / place) % base;
    TYPE_OF_ELEMENTS min = (array[0] / place) % base;
    int i;
    for (int i = 1; i < size; i++) {
        if (((array[i] / place) % base) > max)
            max = ((array[i] / place) % base);

        if (((array[i] / place) % base) < min)
            min = ((array[i] / place) % base);
    }

    matrix[base] = min;

    TYPE_OF_ELEMENTS length = max - min + 1;
    for (int j = 0; j < size; j++)
        matrix[(((array[j] / place) % base) - min)]++;

    // Calculate cumulative count
    for (int j = 1; j < length; j++)
        matrix[j] += matrix[(j - 1)];
}

/**
 * @brief This function implements the radix sort algorithm for a only positive or only negative array, based on counting sort algorithm, storing during the computation the frequences of each element in a matrix.
 * @param array      pointer to the array with the elements to sort, that are only positive or negative.
 * @param size       dimension of the subarray containing all positive or negative elements.
 * @param max        maximum number of digits contained in the longer element of this subarray
 * @param threads    number of threads.
 */
void serviceRadixsort(TYPE_OF_ELEMENTS array[], int size, int max, int threads) {
    // Get maximum element

    int base = 10;
    TYPE_OF_ELEMENTS *vect[max];
    for (int i = 0; i < max; i++) {
        vect[i] = (TYPE_OF_ELEMENTS *)calloc((base + 1), sizeof(TYPE_OF_ELEMENTS));
    }

    TYPE_OF_ELEMENTS *output = (TYPE_OF_ELEMENTS *)calloc(size, sizeof(TYPE_OF_ELEMENTS));
    int i = 0, j = 0;

// Apply counting sort to sort elements based on place value.
#pragma omp parallel for num_threads(threads)
    for (i = 0; i < max; i++) {
        countingSort(array, base, size, i, vect[i], 1);
    }

    int place = 1;

    for (i = 0; i < max; i++) {
        for (j = size - 1; j >= 0; j--) {
            output[vect[i][(((array[j] / place) % base) - vect[i][base])] - 1] = array[j];
            vect[i][(((array[j] / place) % base) - vect[i][base])]--;
        }

        for (int j = 0; j < size; j += 2) {
            array[j] = output[j];
            array[j + 1] = output[j + 1];
        }
        place = place * 10;
    }
}

/**
 * @brief This function searches the maximum value of the array used by the sorting function to  * establish the maximum digits of the numbers.
 * @param n        number of elements in the array.
 * @param arr      pointer to the array used in the sorting function.
 * @param threads  number of threads used in the parallelized version.
 *
 */
TYPE_OF_ELEMENTS getMaxDigit2(int n, TYPE_OF_ELEMENTS *arr, int threads) {
    TYPE_OF_ELEMENTS mx = arr[0];
    int i = 1;
#pragma omp parallel for num_threads(threads) reduction(max \
                                                        : mx) schedule(runtime)
    for (i = 1; i < n; i++) {
        if (abs(arr[i]) > mx)
            mx = abs(arr[i]);
    }
    return mx;
}

/**
 * @brief This function does a brutal-force sorting algorithm on the numbers of the array, analyzing  * them digit by digit.
 * Source by: https://medium.com/geekculture/implementation-and-performance-analysis-of-parallel-and-serial-counting-sort-algorithm-using-openmp-56016f9ccb5c
 * @param n        number of elements in the array.
 * @param arr      pointer to the array used in the sorting function.
 * @param exp	   number used to establish which digit are being analyzed (units, tens, hundreds ...).
 * @param threads  number of threads used in the parallelized version..
 *
 */
void brutalSort(int n, TYPE_OF_ELEMENTS *arr, int exp, int threads) {
    int i, j, count;
    TYPE_OF_ELEMENTS *output = (TYPE_OF_ELEMENTS *)calloc(n, sizeof(TYPE_OF_ELEMENTS));
#pragma omp parallel for private(i, j, count) num_threads(threads) schedule(runtime)
    for (i = 0; i < n; i++) {
        count = 0;
        for (j = 0; j < n; j++) {
            if ((arr[i] / exp) % 10 > (arr[j] / exp) % 10)
                count++;
            else if ((arr[j] / exp) % 10 == (arr[i] / exp) % 10 && j < i)
                count++;
        }
        while (output[count] != 0)
            count++;
        output[count] = arr[i];
    }
    memcpy(arr, output, n * sizeof(output[0]));
}

/**
 * @brief This function implements the radix sort algorithm for a only positive or only negative array, based on counting sort algorithm, storing during the computation the frequences of each element in a matrix.
 * @param arr     pointer to the array with the elements to sort, that are only positive or negative.
 * @param n       dimension of the subarray containing all positive or negative elements.
 * @param threads    number of threads.
 */
void myRadixsort2(int n, TYPE_OF_ELEMENTS *arr, int threads) {
    TYPE_OF_ELEMENTS m = getMaxDigit2(n, arr, threads);
    int exp = 1;
    for (exp = 1; m / exp > 0; exp *= 10)
        brutalSort(n, arr, exp, threads);
}

/**
 * @brief This function instantiates two different arrays for store separately positive and negative values, in order to separate the operation of sorting.
 * @param array      pointer to the array with total elements to sort.
 * @param size       dimension of the array containing all elements.
 * @param threads    number of threads.
 */
void myRadixsort(TYPE_OF_ELEMENTS *array, int size, int threads) {
    TYPE_OF_ELEMENTS *array_pos = (TYPE_OF_ELEMENTS *)malloc(size * sizeof(TYPE_OF_ELEMENTS));
    if (array_pos == NULL)
        perror("Memory Allocation - a");
    TYPE_OF_ELEMENTS *array_neg = (TYPE_OF_ELEMENTS *)malloc(size * sizeof(TYPE_OF_ELEMENTS));
    if (array_neg == NULL)
        perror("Memory Allocation - a");
    int max_pos, max_neg;
    int size_pos, size_neg;

    getMaxDigit(array, size, array_pos, array_neg, &max_pos, &max_neg, &size_pos, &size_neg, threads);
    array_pos = realloc(array_pos, size_pos * sizeof(int));
    array_neg = realloc(array_neg, size_neg * sizeof(int));

    int myNumTeams = 2;
    if (threads == 1)
        myNumTeams = 1;
#pragma omp teams num_teams(myNumTeams) shared(array)
    {
        if (threads <= 1) {
            serviceRadixsort(array_neg, size_neg, max_neg, 1);
            memcpy(array, array_neg, size_neg * sizeof(int));
            serviceRadixsort(array_pos, size_pos, max_pos, 1);
            memcpy(array + size_neg, array_pos, size_pos * sizeof(int));
        } else {
            if (omp_get_team_num() == 0) {
                serviceRadixsort(array_pos, size_pos, max_pos, threads / 2);
                memcpy(array + size_neg, array_pos, size_pos * sizeof(int));
            } else {
                serviceRadixsort(array_neg, size_neg, max_neg, threads / 2);
                memcpy(array, array_neg, size_neg * sizeof(int));
            }
        }
    }
}

/**
 * @brief This function initializes the array structure needed in the program.
 * @param array         pointer to the array to sort in the sorting algorithm.
 * @param size          dimension of array, so the number of elements contained in the array.
 * @param max_digit     maximum number that will be contained in the array, so it determines the maximum element of digit of the elements contained in the array.
 */
void init_structures(TYPE_OF_ELEMENTS **array, int size, int max_digit) {
    TYPE_OF_ELEMENTS *tmp_array = (TYPE_OF_ELEMENTS *)malloc(size * sizeof(TYPE_OF_ELEMENTS));
    if (tmp_array == NULL)
        perror("Memory Allocation - tmp_array");
    srand(time(NULL));
    tmp_array[0] = (TYPE_OF_ELEMENTS)max_digit;
    for (int i = 0; i < size; i++) {
        if (i % 2)
            tmp_array[i] = -(rand() % max_digit);
        else
            tmp_array[i] = (rand() % max_digit);
    }
    *array = tmp_array;
}
