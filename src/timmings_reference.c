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
  for(unsigned N = (1 << 5); N <= (1 << 15); N *= 2){
    report(INFO, "Starting on size = %u", N);
    unsigned histogram[NUM_HIST_BOXES] = {0};
    printf("%u ", N);
    double time_taken, total_time = 0;


    tick(&T);
    star_t *stars = malloc(sizeof(star_t)*N);
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
    time_taken = tock(&T);
    printf("%e ", time_taken);
    total_time += time_taken;

    tick(&T);
    reference_create_random_array(stars, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    reference_sort(stars, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    reference_fill_matrix(stars, matrix, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    reference_create_tally_matrix(matrix, tally, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    reference_generate_histogram(tally, histogram, N-2, NUM_HIST_BOXES);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);
    printf("%e\n", total_time);
    report(PASS, "Finished size %u in >%ens", N, total_time);
    free(stars);
    free(matrix);
    free(tally);
  }
}
