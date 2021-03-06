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
  for(unsigned N = (1 << 5); N <= (1 << 15); N *= 2){
    report(INFO, "Starting on size = %u", N);
    unsigned histogram[NUM_HIST_BOXES] = {0};
    printf("%u ", N);
    double time_taken, total_time = 0;


    tick(&T);
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
    time_taken = tock(&T);
    printf("%e ", time_taken);
    total_time += time_taken;

    tick(&T);
    create_random_array(stars, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    optim_sort(stars, sorted_stars, N);

    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    fill_matrix(sorted_stars, matrix, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    create_tally_matrix(matrix, tally, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    generate_histogram(tally, histogram, N-2, NUM_HIST_BOXES);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);
    printf("%e\n", total_time);
    report(PASS, "Finished size %u in >%ens", N, total_time);
    star_array_free(stars);
    star_array_free(sorted_stars);
    sym_matrix_free(matrix);
    free(tally);
  }
}
