/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  main.c main implementation file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#include <reference_functions.h>
#include <ref_input.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <utilities/logging.h>
#include <utilities/integer.h>
#include <utilities/benchmarking.h>
#include <utilities/vector_types.h>
#include <math.h>
#include <time.h>


int main(int argc, char **argv)
{
  struct timespec T;
  const unsigned N = 42;
  unsigned histogram[NUM_HIST_BOXES] = {0};
  star_t stars[N];
  reference_create_random_array(stars, N);
  //reference_create_ref_star_array(stars, N);
  float_t *matrix =  malloc(sizeof(float_t)*N*N);
  if(!matrix){
    report(FAIL, "Failed to allocate matrix: %s (%d)", strerror(errno), errno);
    return -1;
  }
  float_t *tally =  malloc(sizeof(float_t)*N*N);
  if(!tally){
    report(FAIL, "Failed to allocate tally: %s (%d)", strerror(errno), errno);
    return -1;
  }
  reference_sort(stars, N);
  reference_print_stars(stars, N);
  reference_fill_matrix(stars, matrix, N);
  //reference_print_matrix(matrix, N);
  reference_create_tally_matrix(matrix, tally, N);
  //reference_print_tally_matrix(tally, N-2);
  hist_param_t histparams = reference_generate_histogram(tally, histogram, N-2, NUM_HIST_BOXES);
  reference_display_histogram(histogram, histparams);
  free(tally);
  free(matrix);
}
