/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  funcs.h functions header file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#ifndef _FUNCS_H
#define _FUNCS_H

#include "common.h"

unsigned rng();
float_t rngf(float_t min, float_t max);

sym_matrix_t sym_matrix_intitialize(unsigned N);
void sym_matrix_free(sym_matrix_t sym_matrix);
// unsigned sym_matrix_coordinate(unsigned x, unsigned y, unsigned size);
// float_t sym_matrix_get(sym_matrix_t matrix, unsigned x, unsigned y, unsigned size);
// void sym_matrix_set(sym_matrix_t matrix, float_t value, unsigned x, unsigned y, unsigned size);
// void sym_matrix_set_vec(sym_matrix_t matrix, v8f value, unsigned x, unsigned y, unsigned size);

star_array_t star_array_initialize(size_t size);
void star_array_free(star_array_t star_array);
void create_random_array(star_array_t stars, unsigned size);

//void create_random_array(star_array_t star_array, int size);
void sort(star_array_t array, unsigned offset, unsigned n);
void print_stars(star_array_t array, unsigned n);

void fill_matrix(star_array_t array, float_t *matrix, unsigned size);
void print_matrix(sym_matrix_t matrix, unsigned n);

hist_param_t generate_histogram(float_t *matrix, unsigned *histogram, unsigned mat_size, unsigned hist_size);
void create_tally_matrix(float_t *in, float_t* out, unsigned N);

void display_histogram(int *histogram, hist_param_t histparams);

#endif
