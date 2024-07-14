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
        @file main.c
*/

#include <stdio.h>
#include <stdlib.h>

#include "RadixSort.h"

int main(int argc, char const *argv[]) {
    TYPE_OF_ELEMENTS *array;
    double time_sort = 0.0, time_init = 0.0;
    int length, threads, algorithm, max_digit;
    if (argc < 5) {
        printf("ERROR! Usage: ./main rows columns threads");
        exit(1);
    }

    length = atoi(argv[1]);
    threads = atoi(argv[2]);
    algorithm = atoi(argv[3]);
    max_digit = atoi(argv[4]);

    if (algorithm == 0)
        length *= length;

    STARTTIME(1);
    init_structures(&array, length, max_digit);
    ENDTIME(1, time_init);
    if (algorithm == 0) {
        STARTTIME(2);
        myRadixsort(array, length, threads);
        ENDTIME(2, time_sort);
    } else {
        STARTTIME(2);
        myRadixsort2(length, array, threads);
        ENDTIME(2, time_sort);
    }
    printf("%d;%d;%d;%f;%f\n", algorithm, length, threads, time_init, time_sort);

    free(array);

    return 0;
}
