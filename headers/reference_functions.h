/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  funcs.h functions header file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#ifndef _REF_FUNCS_H
#define _REF_FUNCS_H

#include "common.h"

void reference_print_tally_matrix(float_t *matrix, unsigned size);
void reference_create_random_array(star_t * array, int size);
void reference_sort(star_t* array, int n);
void reference_print_stars(star_t* array, int n);

void reference_fill_matrix(star_t * array, float_t *matrix, int size);
void reference_print_matrix(float_t* matrix, int n);

hist_param_t reference_generate_histogram(float_t *matrix, int *histogram, int mat_size, int hist_size);
void reference_create_tally_matrix(float_t *in, float_t* out, unsigned N);

void reference_display_histogram(int *histogram, hist_param_t histparams);

#endif
