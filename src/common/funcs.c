/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  funcs.c functions implementation file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

 #ifdef _WIN32
  #define _CRT_RAND_S
  #include <stdlib.h>
  #include <stdint.h>
  #define RAND_MAX UINT_MAX
 #else
  #include <stdlib.h>
 #endif

#include "funcs.h"
#include <utilities/logging.h>
#include <utilities/integer.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

unsigned rng(){
  #ifdef _WIN32
    unsigned r;
    rand_s(&r);//windows rng sucks.
    return r;
  #else
    return rand();
  #endif
}
float_t rngf(float_t min, float_t max){
  return min+((float_t)rng()/RAND_MAX)*(max-min);
}

sym_matrix_t sym_matrix_intitialize(unsigned N){
  return _mm_malloc(sizeof(float_t)*64*triag_nr((N+7)/8), 32);
}

void sym_matrix_free(sym_matrix_t sym_matrix){
  _mm_free(sym_matrix);
}
//
// unsigned sym_matrix_coordinate(unsigned x, unsigned y, unsigned size){ //http://stackoverflow.com/q/19143657/6019199
//   if(x < y){//TODO: can we do this branchless?
//     return sym_matrix_coordinate(y,x,size);
//   }
//   return x+triag_nr(size)-triag_nr(size-y)-y; //triag_nr can overflow for > sqrtf(UNSIGNED_MAX), 64k for 32 bit.
// }
//
// mat8x8f_t sym_matrix_get_block(sym_matrix_t matrix, unsigned x, unsigned y, unsigned size){
//   return matrix[sym_matrix_coordinate(x,y,size)];
// }
//
// void sym_matrix_set_vec(float_t *matrix, v8f value, unsigned x, unsigned y, unsigned size){
//   *(v8f*)(matrix+sym_matrix_coordinate(x/8,y/8,size/8)) = value;
// }
//
// void sym_matrix_set(float_t *matrix, float_t value, unsigned x, unsigned y, unsigned size){
//   matrix[sym_matrix_coordinate(x,y,size)] = value;
// }

unsigned sym_matrix_get_block_index(unsigned x, unsigned y, unsigned size){
  return x+triag_nr(size)-triag_nr(size-y)-y;
}

star_array_t star_array_initialize(size_t size){
  star_array_t star_array;
  star_array.spectralType = malloc(sizeof(char)*size);
  star_array.index = malloc(sizeof(unsigned)*size);
  size = (size+7) & ~7;//round up to multiple of 8.
  star_array.magnitude = _mm_malloc(sizeof(float_t)*size, 32);
  star_array.subType = _mm_malloc(sizeof(float_t)*size, 32);
  star_array.position.x = _mm_malloc(sizeof(float_t)*size, 32);
  star_array.position.y = _mm_malloc(sizeof(float_t)*size, 32);
  star_array.position.z = _mm_malloc(sizeof(float_t)*size, 32);
  return star_array;
}

void star_array_free(star_array_t star_array){
  free(star_array.spectralType);
  free(star_array.index);
  _mm_free(star_array.magnitude);
  _mm_free(star_array.subType);
  _mm_free(star_array.position.x);
  _mm_free(star_array.position.y);
  _mm_free(star_array.position.z);
}

void create_random_array(star_array_t stars, unsigned N){
  const unsigned n_types = 9;
  const char types[9] = {'O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'T'};
  for(unsigned u = 0; u < N; u++){
    stars.index[u] = u;
    stars.spectralType[u] = types[rng()%n_types];
    stars.subType[u] = rng() % 9;
    stars.magnitude[u] = rngf(-10, +20);
    stars.position.x[u] = rngf(-1e5, 1e5);
    stars.position.y[u] = rngf(-1e5, 1e5);
    stars.position.z[u] = rngf(-1e5, 1e5);
    //TODO: change print function to work with this.
    //sprintf(stars.designation[u], "%c%d.%d", stars.spectralType[u], stars.subType[u], stars.index[u]); //not safe for index >= 10^7
  }
}

void print_stars(star_array_t array, unsigned n)
{
  printf("\nprint_stars, n = %u:\n", n);
  for(unsigned i = 0; i<n; i++)
    printf("%c%d.%u ",array.spectralType[i], (int)array.subType[i], array.index[i]);
  printf("\n");
}

float_t len_3d(float_t x, float_t y, float_t z){
  return sqrtf(x*x+y*y+z*z);
}

v8f len_3d_vec(v8f x, v8f y, v8f z){
  return _mm256_sqrt_ps(x*x+y*y+z*z);
}

float_t dist_3d(float_t x0, float_t y0, float_t z0, float_t x1, float_t y1, float_t z1){
  float_t dx = (x0-x1);
  float_t dy = (y0-y1);
  float_t dz = (z0-z1);
  return len_3d(dx, dy, dz);
}

//TODO: use struct index notation.
v8f dist_3d_vec(v8f x0, v8f y0, v8f z0, v8f x1, v8f y1, v8f z1){
  v8f dx = (x0-x1);
  v8f dy = (y0-y1);
  v8f dz = (z0-z1);
  return len_3d_vec(dx, dy, dz);
}

float_t starfunc(float_t a, float_t b){ //symmetrical, pure
  //v8f x = _mm256_broadcast_ss(stars.subType[a]);
  //v8f y = *((v8d*)&stars.subType[b]);
  const float_t c = 1.0/0.6; //compile time evaluation.
  return sqrtf(b+a*(1+b*c));//2 * FMA, 1 SQRT//TODO:replace sqrt with something faster http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
}

v8f starfunc_vec(v8f a, v8f b){ //symmetrical, pure
  //v8f x = _mm256_broadcast_ss(stars.subType[a]);
  //v8f y = *((v8d*)&stars.subType[b]);
  const float_t c = 1.0/0.6; //compile time evaluation.
  return _mm256_sqrt_ps(b+a*(1+b*c));//2 * FMA, 1 SQRT//TODO:replace sqrt with something faster http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
}


//TODO: use log based enumeration sort?
void sort(star_array_t array, unsigned start, unsigned end){//Quicksort, modified from http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
  //const star_t origin_star = {.position = {.x = 0, .y = 0, .z = 0}};
  unsigned u, w;
  if (end-start < 2){
    return;
  }
  //TODO: compute and store the distances before sorting.
  float_t x, y, z;
  int pivot_loc = (start+end)/2;
  x =  array.position.x[pivot_loc];
  y =  array.position.y[pivot_loc];
  z =  array.position.z[pivot_loc];
  float_t pivot = len_3d(x, y, z);
  float_t L;
  for (u = start-1, w = end;;) {
    do{
      u++;
      x =  array.position.x[u];
      y =  array.position.y[u];
      z =  array.position.z[u];
      L = len_3d(x, y, z);
    }while (L < pivot);
    do{
      w--;
      x =  array.position.x[w];
      y =  array.position.y[w];
      z =  array.position.z[w];
      L = len_3d(x, y, z);
    }while(pivot < L);
    if (u >= w){
      break;
    }

    //swaps, big.
    char tmp_type = array.spectralType[u];
    array.spectralType[u] = array.spectralType[w];
    array.spectralType[w] = tmp_type;

    unsigned tmp_index = array.index[u];
    array.index[u] = array.index[w];
    array.index[w] = tmp_index;

    float_t tmp_magnitude = array.magnitude[u];
    array.magnitude[u] = array.magnitude[w];
    array.magnitude[w] = tmp_magnitude;

    float_t tmp_sub = array.subType[u];
    array.subType[u] = array.subType[w];
    array.subType[w] = tmp_sub;

    float_t tmp_pos = array.position.x[u];
    array.position.x[u] = array.position.x[w];
    array.position.x[w] = tmp_pos;

    tmp_pos = array.position.y[u];
    array.position.y[u] = array.position.y[w];
    array.position.y[w] = tmp_pos;

    tmp_pos = array.position.z[u];
    array.position.z[u] = array.position.z[w];
    array.position.z[w] = tmp_pos;
  }
  sort(array, start, u);
  sort(array, u, end);
}

void fill_matrix(star_array_t array, sym_matrix_t matrix, unsigned size)
{
  report(INFO, "starting to fill matrix");
  for(unsigned u = 0; u < size; u++){ //loop over each element in array
    v8f subType = _mm256_broadcast_ss(array.subType+u);
    v8f x = _mm256_broadcast_ss(array.position.x+u);
    v8f y = _mm256_broadcast_ss(array.position.y+u);
    v8f z = _mm256_broadcast_ss(array.position.z+u);
    unsigned w;
    report(INFO, "starting vector loop");
    for(w = u/8; w < size/8; w++){//set a block of the matrix
      unsigned flat_index = sym_matrix_get_block_index(w, u)
      v8f subType_b = *(v8f*)(array.subType+w);
      v8f func = starfunc_vec(subType, subType_b);
      v8f x_b = *(v8f*)(array.position.x+w);
      v8f y_b = *(v8f*)(array.position.y+w);
      v8f z_b = *(v8f*)(array.position.z+w);
      v8f dist = dist_3d_vec(
        x, y, z,
        x_b, y_b, z_b
      );
     sym_matrix_set_vec(matrix,  func+dist, w, u, size);
    }
    report(INFO, "starting scalar loop");
    for(; w < size; w++){
      float_t func = starfunc(subType[0], array.subType[w]);
      float_t dist = dist_3d(
        x[0], y[0], z[0],
        array.position.x[w], array.position.y[w], array.position.z[w]
      );
     sym_matrix_set(matrix,  func+dist, w, u, size);
    }
  }
}

void print_matrix(float_t* theMatrix, unsigned n)
{
  unsigned i, j;
  printf("\nprint_matrix, n = %u:\n", n);
  for(i = 0 ; i < n; i++)
    {
      for(j = 0 ; j < n ; j++)
	printf("%1.4e " , sym_matrix_get(theMatrix, j, i, n));
      putchar('\n');
    }
}

void create_tally_matrix(float_t *in, float_t* out, unsigned N){
  for(unsigned u = 1; u < N-1; u++){
    for(unsigned w = u; w < N-1; w++){
      float_t center = sym_matrix_get(in, w, u, N);
      float_t up = sym_matrix_get(in, w, u+1, N);
      float_t down = sym_matrix_get(in, w, u-1, N);
      float_t left = sym_matrix_get(in, w-1,u, N);
      float_t right = sym_matrix_get(in, w+1, u, N);
      float_t mean_abs_diff = 0.25*(fabs(center-up) + fabs(center-down) + fabs(center-left) + fabs(center-right));
      sym_matrix_set(out, mean_abs_diff, w-1, u-1, N-2);
    }
  }
}

unsigned get_hist_index(float_t min, float_t bin_size, unsigned hist_size, float_t val){
  unsigned hist_index;
  for(hist_index = 0; hist_index < hist_size-1; hist_index++){//TODO: direct computation
    if(min + bin_size*(1+hist_index) > val){
      break;
    }
  }
  return hist_index;
}

hist_param_t generate_histogram(float_t *matrix, unsigned *histogram, unsigned mat_size, unsigned hist_size){
  float_t min = FLT_MAX, max = -min;
  for(unsigned u = 0; u < mat_size*(1+mat_size)/2; u++){//TODO: calculate when filling matrix
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
    unsigned b = triag_nr(mat_size)-triag_nr(mat_size-u);
    unsigned index = get_hist_index(min, param.bin_size, hist_size, matrix[b]);
    //report(INFO, "[%u] adding 1x%e to index %u",b, matrix[b], index);
    histogram[index] += 1;
    for(unsigned o = 1;o < mat_size-u;o++){
      unsigned index = get_hist_index(min, param.bin_size, hist_size, matrix[b+o]);
      histogram[index] += 2;
      //report(INFO, "[%u] adding 2x%e to index %u",b+o, matrix[b+o], index);
    }
  }
  return param;
}

void display_histogram(int *histogram, hist_param_t histparams)
{
  printf("\ndisplay_histogram:\n");
  unsigned i,j;
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
