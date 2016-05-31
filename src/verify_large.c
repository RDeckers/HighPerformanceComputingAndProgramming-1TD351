/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  main.c main implementation file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#include <funcs.h>
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
  unsigned N = 42;
  unsigned histogram[NUM_HIST_BOXES] = {0};
  star_array_t stars = star_array_initialize(N);
  star_array_t sorted_stars = star_array_initialize(N);
  sym_matrix_t matrix = sym_matrix_intitialize(N);
  if(!matrix){
    report(FAIL, "Failed to allocate matrix: %s (%d)", strerror(errno), errno);
    return -1;
  }
  float_t *tally =  malloc(sizeof(float_t)*triag_nr(N-2));
  if(!tally){
    report(FAIL, "Failed to allocate tally: %s (%d)", strerror(errno), errno);
    return -1;
  }
  create_random_array(stars, N);
  report(INFO, "starting sort");
  optim_sort(stars, sorted_stars, N);
  report(INFO, "sorted");
  print_stars(sorted_stars, N);
  fill_matrix(sorted_stars, matrix, N);
  //print_matrix(matrix, N);
  create_tally_matrix(matrix, tally, N);
  //print_tally_matrix(tally, N-2);
  hist_param_t histparams = generate_histogram(tally, histogram, N-2, NUM_HIST_BOXES);
  display_histogram(histogram, histparams);
  star_array_free(stars);
  star_array_free(sorted_stars);
  sym_matrix_free(matrix);
  free(tally);
}
