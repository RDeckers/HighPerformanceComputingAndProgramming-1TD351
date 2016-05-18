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
#include <utilities/benchmarking.h>
#include <math.h>
#include <time.h>


int main(int argc, char **argv)
{
  int N, i;
  struct timespec T;
  int use_refference = 1;
  if(argc >= 2){
    report(WARN, "Superflous arguments given");
  }
  N = use_refference ? 7 : atoi(argv[1]);
  star_t *stars = malloc(N*sizeof(star_t));
  if(!stars){
    report(FAIL, "Failed to allocate star array: %s (%d)", strerror(errno), errno);
    return -1;
  }

  //create_random_array(stars, N);
  tick(&T);
  use_refference ? create_ref_star_array(stars, N) : create_random_array(stars, N);
  report(PASS, "Created star array in %ens", elapsed_since(&T));
  print_stars(stars, N);

  tick(&T);
  sort(stars, N);
  report(PASS, "Sorted star array in %ens", elapsed_since(&T));
  print_stars(stars, N);

  float_t *matrix = malloc(sizeof(float_t)*N*N);
  if(!matrix){
    report(FAIL, "Failed to allocate matrix: %s (%d)", strerror(errno), errno);
    return -1;
  }

  float_t *tally = malloc(sizeof(float_t)*(N-1)*(N-1));
  if(!matrix){
    report(FAIL, "Failed to allocate tally: %s (%d)", strerror(errno), errno);
    return -1;
  }

  tick(&T);
  fill_matrix(stars, matrix, N);
  report(PASS, "Filled matrix in %ens", elapsed_since(&T));

  print_matrix(matrix, N);

  tick(&T);
  create_tally_matrix(matrix, tally, N);
  report(PASS, "Created tally matrix in %ens", elapsed_since(&T));
  print_matrix(tally, N-2);
  int *histogram = calloc(NUM_HIST_BOXES,sizeof(int));
  tick(&T);
  hist_param_t histparams = generate_histogram(tally, histogram, N-2, NUM_HIST_BOXES);
  report(PASS, "generated histogram in %ens", elapsed_since(&T));

  display_histogram(histogram, histparams);
}
