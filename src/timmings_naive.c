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
  struct timespec T;
  for(unsigned N = (1 << 5); N <= (1 << 15); N *= 2){
    int histogram[NUM_HIST_BOXES] = {0};
    report(INFO, "Starting on size = %u", N);
    star_t *stars = malloc(N*sizeof(star_t));
    if(!stars){
      report(FAIL, "Failed to allocate star array: %s (%d)", strerror(errno), errno);
      return -1;
    }
    float_t *matrix = malloc(sizeof(float_t)*N*N);
    if(!matrix){
      report(FAIL, "Failed to allocate matrix: %s (%d)", strerror(errno), errno);
      return -1;
    }
    float_t *tally = malloc(sizeof(float_t)*(N-2)*(N-2));
    if(!matrix){
      report(FAIL, "Failed to allocate tally: %s (%d)", strerror(errno), errno);
      return -1;
    }
    printf("%u ", N);
    double time_taken, total_time = 0;
    tick(&T);
    create_random_array(stars, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    sort(stars, N);
    time_taken = tock(&T);
    total_time += time_taken;
    printf("%e ", time_taken);

    tick(&T);
    fill_matrix(stars, matrix, N);
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
    free(stars);
    free(matrix);
    free(tally);
  }
}