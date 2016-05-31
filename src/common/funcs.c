/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  funcs.c functions implementation file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/
#include "funcs.h"
#include <utilities/logging.h>
#include <utilities/integer.h>
#include <utilities/rng.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdint.h>

sym_matrix_t sym_matrix_intitialize(unsigned N){
  return _mm_malloc(sizeof(float_t)*64*triag_nr((N+7)/8), 32); //8x8 blocks, upper triangular part.
}

void sym_matrix_free(sym_matrix_t sym_matrix){
  _mm_free(sym_matrix);
}

star_array_t star_array_initialize(size_t size){
  star_array_t star_array;
  star_array.spectralType = malloc(sizeof(char)*size);
  star_array.index = _mm_malloc(sizeof(unsigned)*size, 32);
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
  _mm_free(star_array.index);
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
    stars.spectralType[u] = types[fast_rng()%n_types];
    stars.subType[u] = fast_rng() % 9;
    stars.magnitude[u] = fast_rngf(-10, +20);
    stars.position.x[u] = fast_rngf(-1e5, 1e5);
    stars.position.y[u] = fast_rngf(-1e5, 1e5);
    stars.position.z[u] = fast_rngf(-3e3, 3e3);
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

v8f dist_3d_vec(v8f x0, v8f y0, v8f z0, v8f x1, v8f y1, v8f z1){
  v8f dx = (x0-x1);
  v8f dy = (y0-y1);
  v8f dz = (z0-z1);
  return len_3d_vec(dx, dy, dz);
}

float_t starfunc(float_t a, float_t b){ //symmetrical, pure
  const float_t c = 1.0/0.6; //compile time evaluation.
  //TODO: lookup?
  return sqrtf(b+a*(1+b*c));//2 * FMA, 1 SQRT//TODO:replace sqrt with something faster http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
}

v8f starfunc_vec(v8f a, v8f b){ //symmetrical, pure
  const float_t c = 1.0/0.6; //compile time evaluation.
  return _mm256_sqrt_ps(b+a*(1+b*c));//2 * FMA, 1 SQRT//TODO:replace sqrt with something faster http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
}

void quicksort(uint64_t *array, unsigned size){
  unsigned u, w;
  if (size < 2){
    return;
  }
  uint64_t pivot = array[size/2];
  for (u = 0, w = size-1;; u++, w--) {
    while (array[u] < pivot){
      u++;
    }
    while(pivot < array[w]){
      w--;
    }
    if (u >= w){
      break;
    }
    uint64_t temp = array[u];
    array[u] = array[w];
    array[w] = temp;
  }
  quicksort(array, u);
  quicksort(array+u, size-u);
}

void optim_sort(star_array_t array, star_array_t out_array, unsigned size){
  union sort_tag{
    uint64_t fused;
    struct {
      uint32_t index;
      float_t distance;
    }__attribute__((packed)) seperate;
  };
  uint64_t *sortable = _mm_malloc(sizeof(uint64_t)*size,32);
  //TODO: bugfix
  // unsigned u;
  // for(u = 0; u < size/8; u++){
  //   v8f x = _mm256_load_ps(array.position.x+8*u);
  //   v8f y = _mm256_load_ps(array.position.y+8*u);
  //   v8f z = _mm256_load_ps(array.position.z+8*u);
  //   v8f lengths = len_3d_vec(x, y, z);
  //   v8f indices =  _mm256_load_ps((float_t*)(array.index+8*u));
  //   _mm256_store_ps((float_t*)(sortable+16*u), _mm256_unpackhi_ps(indices, lengths));
  //   _mm256_store_ps((float_t*)(sortable+16*u+8), _mm256_unpacklo_ps(indices, lengths));
  // }
  for(unsigned u = 0 ;u < size; u++){
    float_t x =  array.position.x[u];
    float_t y =  array.position.y[u];
    float_t z =  array.position.z[u];
    float_t L = len_3d(x, y, z);
    uint32_t index = array.index[u];
    union sort_tag unified = {.seperate = {.distance = L, .index = index}};
    sortable[u] = (unified).fused;
  }
  quicksort(sortable, size);
  for(unsigned u = 0; u < size; u++){
    unsigned w = ((union sort_tag*)sortable)[u].seperate.index;
    //swap.
    out_array.spectralType[u] = array.spectralType[w];
    out_array.index[u] = w;
    out_array.magnitude[u] = array.magnitude[w];
    out_array.subType[u] = array.subType[w];
    out_array.position.x[u] = array.position.x[w];
    out_array.position.y[u] = array.position.y[w];
    out_array.position.z[u] = array.position.z[w];
  }
  _mm_free(sortable);
}

//TODO: use log (or something) based enumeration sort first (list of varrarray?)
//TODO: create an array of struct{distance, index}(packed) and sort using uint64_t comparisons.
void sort(star_array_t array, unsigned start, unsigned end){//Quicksort, modified from http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
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

/*
  TODO: merge with tally?
*/
void fill_matrix(star_array_t array, sym_matrix_t matrix, unsigned size)
{
  const unsigned block_count = (size+7)/8;
  unsigned row_offset = 0;
  for(unsigned u = 0; u < size; u++){ //loop over each row, vertical
    v8f subType_r = _mm256_broadcast_ss(array.subType+u);
    v8f x_r = _mm256_broadcast_ss(array.position.x+u);
    v8f y_r = _mm256_broadcast_ss(array.position.y+u);
    v8f z_r = _mm256_broadcast_ss(array.position.z+u);
    for(unsigned w = 0; w < block_count-u/8; w++){//for each block in row.
      v8f subType_c = ((v8f*)array.subType)[w+u/8];
      v8f func = starfunc_vec(subType_r, subType_c);
      v8f x_c = ((v8f*)(array.position.x))[w+u/8];
      v8f y_c = ((v8f*)(array.position.y))[w+u/8];
      v8f z_c = ((v8f*)(array.position.z))[w+u/8];
      v8f dist = dist_3d_vec(
        x_r, y_r, z_r,
        x_c, y_c, z_c
      );
      ((v8f*)matrix)[row_offset+w] = func+dist;
    }
    row_offset += (block_count-u/8); //offset into start of this row.
  }
}

void print_matrix(sym_matrix_t matrix, unsigned n)//TODO: fix
{
  unsigned index =0;
  unsigned block_count = (n+7)/8;
  for(unsigned y = 0; y < n; y++){
    for(unsigned x = 0; x < y; x++){
      printf("         ");
    }
    for(unsigned x = y; x < n; x++){
      printf("%1.1e ", matrix[index+x-y+y%8]);
    }
    index += 8*(block_count-y/8);
    puts("");
  }
}
void set_min_max(float_t val, float_t *min, float_t *max){
  if(val < *min){
    *min = val;
  }
  if(val > *max){
    *max = val;
  }
}

/*TODO:
then include max/min for histogram
*/
 float_t min;
 float_t max;

v8f fabs_vec(v8f a){
  const v8f sign_mask = _mm256_set1_ps(-0.f); // -0.f = 1 << 31
  return _mm256_andnot_ps(sign_mask, a);
}

void create_tally_matrix(float_t *in, float_t* out, unsigned N){
  float_t up, center, left, down, right;
  v8f center_v, left_v, right_v, up_v;
  unsigned block_count = (N+7)/8;
  unsigned in_row_index = 8*block_count+1;//always points to the first element of the current row.
  unsigned old_row_index = 1;
  unsigned out_index = 0;
  unsigned old_out_index;
  for(unsigned y = 1; y < N-1;){
    unsigned x;
    center = in[in_row_index];
    float_t left_and_up = in[old_row_index];//symmetry on diagonal
    float_t right_and_down = in[in_row_index+1];//symmetry on diagonal
    out[out_index] = 0.5*(fabs(center-left_and_up)+fabs(center-right_and_down));
    if(y != 1){
      out[old_out_index] += 0.25*fabs(center-left_and_up);
    }
    x= 1;
    for(x = 1; x+8 < N-y-1; x+= 8){//TODO: N < 8?
      left_v = _mm256_load_ps(in+in_row_index+x-1);//in_row_index points to central element
      center_v = _mm256_loadu_ps(in+in_row_index+x);
      right_v = _mm256_loadu_ps(in+in_row_index+x+1);
      up_v = _mm256_loadu_ps(in+old_row_index+x); //old_row_index points to the start of the old row, offset by 1.
      _mm256_storeu_ps(out+out_index+x,0.25*(fabs_vec(center_v-left_v)+fabs_vec(center_v-right_v)+fabs_vec(center_v-up_v))); //left right and up
      if(y != 1){
        _mm256_storeu_ps(out+old_out_index+x,  _mm256_loadu_ps(out+old_out_index+x)+0.25*fabs_vec(center_v-up_v));//add the 'down' element to the previous row.
      }
    }
    for(;x < N-y-1;x++){
      left = in[in_row_index+x-1];
      center = in[in_row_index+x];
      right = in[in_row_index+x+1];
      up = in[old_row_index+x];
      out[out_index+x] = 0.25*(fabs(center-left)+fabs(center-right)+fabs(center-up));
      if(y != 1){
        out[old_out_index+x] += 0.25*fabs(center-up);
      }
    }
    old_out_index = out_index+1;
    out_index += (N-y-1);
    old_row_index = in_row_index+1;
    y++;
    in_row_index += 8*(block_count-y/8)+1;
    //getchar();
  }
}

unsigned get_hist_index(float_t min, float_t bin_size, unsigned hist_size, float_t val){
  unsigned hist_index = (val-min)/bin_size;
  return hist_index >= hist_size ? hist_size-1 : hist_index;
}

void print_tally_matrix(float_t *matrix, unsigned size){
  unsigned index = 0;
  for(unsigned y = 0; y < size; y++){
    for(unsigned x = 0; x < y; x++){
      printf("         ");
    }
    for(unsigned x = y; x < size; x++){
      printf("%1.1e ", matrix[index++]);
    }
    puts("");
  }
}

hist_param_t generate_histogram(float_t *matrix, unsigned *histogram, unsigned mat_size, unsigned hist_size){
  unsigned bounds = triag_nr(mat_size);
  float_t min = FLT_MAX, max = -min;
  for(unsigned u = 0; u < bounds; u++){//TODO: calculate when filling matrix
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
    .bin_size = (max-min)/(hist_size)
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
