/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  funcs.c functions implementation file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#include <reference_functions.h>
#include <utilities/logging.h>
#include <utilities/rng.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

star_t reference_create_random_star(){
  const unsigned n_types = 9;
  const char types[9] = {'O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'T'};
  return (star_t){
    //.index = 0;//unused
    .spectralType = types[fast_rng()%n_types],           // random: O, B, A, F, G, K, M, L, T
    .subType = fast_rng() % 9,      // random: 0-9
    .magnitude = fast_rngf(-10, +20),           // random: (-10, +20)
    //.designation[9] = ""; 	       // sprintf("%c%d.%d", spectralType, subType, index)
    .position = {
      .x = fast_rngf(-1e5, 1e5),
      .y = fast_rngf(-1e5, 1e5),
      .z = fast_rngf(-3e3, 3e3)
    }
  };
}

void reference_create_random_array(star_t * stars, int size)
{
  for(unsigned u = 0; u < size; u++){
    stars[u] = reference_create_random_star();
    stars[u].index = u;
  }
}

void reference_print_tally_matrix(float_t *matrix, unsigned size){
  unsigned index = 0;
  for(unsigned y = 0; y < size; y++){
    for(unsigned x = 0; x < y; x++){
      printf("         ");
    }
    for(unsigned x = y; x < size; x++){
      printf("%1.1e ", matrix[y*size+x]);
    }
    puts("");
  }
}

void reference_print_stars(star_t* array, int n)
{
  int i;
  printf("\nprint_stars, n = %d:\n", n);
  for(i = 0; i<n; i++)
    printf("%c%d.%u ",array[i].spectralType, (int)array[i].subType, array[i].index);
    //printf("%s ",array[i].designation);
  printf("\n");
}

float_t reference_dist_3d(float_t x0, float_t y0, float_t z0, float_t x1, float_t y1, float_t z1){
  float_t dx = (x0-x1);
  float_t dy = (y0-y1);
  float_t dz = (z0-z1);
  return sqrt(dx*dx+dy*dy+dz*dz);
}

float_t reference_star_dist(star_t a, star_t b){ //symmetrical, pure
  return reference_dist_3d(a.position.x, a.position.y, a.position.z, b.position.x, b.position.y, b.position.z);
}

float_t reference_starfunc(star_t a, star_t b){ //symmetrical, pure
  unsigned short x = a.subType;
  unsigned short y = b.subType;
  return sqrt(x + y + x*y/0.6);
}


void reference_sort(star_t* array, int n){//Quicksort, modified from http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
  const star_t origin_star = {.position = {.x = 0, .y = 0, .z = 0}};
  int i, j;
  if (n < 2){
    return;
  }
  float_t pivot = reference_star_dist(array[n / 2], origin_star);
  for (i = 0, j = n - 1;; i++, j--) {
    while (reference_star_dist(array[i], origin_star) < pivot){
      i++;
    }
    while (pivot < reference_star_dist(array[j], origin_star)){
      j--;
    }
    if (i >= j){
      break;
    }
    star_t tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
  }
  reference_sort(array, i);
  reference_sort(array + i, n - i);
}

void reference_fill_matrix(star_t *array, float_t *matrix, int size)
{
  for(unsigned u = 0; u < size; u++){
    for(unsigned w = 0; w < size; w++){
      matrix[u*size+w] = reference_starfunc(array[u], array[w]) + reference_star_dist(array[u], array[w]);
    }
  }
}

void reference_print_matrix(float_t* theMatrix, int n)
{
  int i, j;
  printf("\nprint_matrix, n = %d:\n", n);
  for(i = 0 ; i < n; i++)
    {
      for(j = 0; j < i; j++){
        printf("         ");
      }
      for(j = i ; j < n ; j++)
	printf("%1.1e " , theMatrix[i*n+j]);
      putchar('\n');
    }
}

void reference_create_tally_matrix(float_t *in, float_t* out, unsigned N){
  for(unsigned u = 1; u < N-1; u++){
    for(unsigned w = 1; w < N-1; w++){
      float_t center = in[u*N+w];
      float_t up = in[(u+1)*N+w];
      float_t down = in[(u-1)*N+w];
      float_t left = in[u*N+w-1];
      float_t right = in[u*N+w+1];
      out[(u-1)*(N-2)+w-1] = 0.25*(fabs(center-up) + fabs(center-down) + fabs(center-left) + fabs(center-right));
    }
  }
}

hist_param_t reference_generate_histogram(float_t *matrix, int *histogram, int mat_size, int hist_size){
  float_t min = FLT_MAX, max = -min;
  for(unsigned u = 0; u < mat_size*mat_size; u++){
    if(matrix[u] > max){
      max = matrix[u];
    }
    if(matrix[u] < min){
      min = matrix[u];
    }
  }
  const hist_param_t param = {
    .hist_size = hist_size,
    .min = min,
    .max = max,
    .bin_size = (max-min)/hist_size
  };
  for(unsigned u = 0; u < mat_size; u++){
    for(unsigned w = 0; w < mat_size; w++){
      unsigned hist_index;
      for(hist_index = 0; hist_index < hist_size-1; hist_index++){
        if(param.min + param.bin_size*(1+hist_index) > matrix[u*mat_size+w]){
          break;
        }
      }
      histogram[hist_index] += 1;
    }
  }
  return param;
}

void reference_display_histogram(int *histogram, hist_param_t histparams)
{
  printf("\ndisplay_histogram:\n");
  int i,j;
  for(i = 0; i < histparams.hist_size && histparams.bin_size*i < histparams.max; i++)
    {
      printf("%11.3e ", histparams.bin_size*i+histparams.min);
    }
  printf("%11.3e\n", histparams.max);
  for(j = 0; j < i; j++)
    {
      printf("%11d ", histogram[j]);
    }
  printf("\n");
}
