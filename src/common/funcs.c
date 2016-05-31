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

sym_matrix_t sym_matrix_intitialize(unsigned N){
  return _mm_malloc(sizeof(float_t)*64*triag_nr((N+7)/8), 32); //8x8 blocks, upper triangular part.
}

void sym_matrix_free(sym_matrix_t sym_matrix){
  ////report(WARN, "freeing %p", sym_matrix);
  _mm_free(sym_matrix);
}
//
unsigned sym_matrix_coordinate(unsigned x, unsigned y, unsigned size){ //http://stackoverflow.com/q/19143657/6019199
  if(x < y){//TODO: can we do this branchless?
    return sym_matrix_coordinate(y,x,size); //if lower half, return upper half equivalent.
  }
  //base offset + (size of whole matrix - size of all rows below - offset into x.)
  unsigned block_count = (size+7)/8;
  unsigned base = x;
  unsigned base_offset = y;
  unsigned size_matrix = triag_nr(block_count)*64;
  unsigned size_below = triag_nr(block_count-(y/8))*64;//blocks below including current block,
  size_below -= (y%8)*(8*(block_count-y/8)); //rows into current block * size of current rows
  return base - base_offset + size_matrix - size_below; //triag_nr can overflow for > sqrt(UNSIGNED_MAX), 64k for 32 bit.
}


// void sym_matrix_set_vec(float_t *matrix, v8f value, unsigned x, unsigned y, unsigned size){
//   *(v8f*)(matrix+sym_matrix_coordinate(x/8,y/8,size/8)) = value;
// }

void sym_matrix_add(sym_matrix_t matrix, float_t value, unsigned x, unsigned y, unsigned size){
  matrix[sym_matrix_coordinate(x,y,size)] += value;
}

void sym_matrix_set(sym_matrix_t matrix, float_t value, unsigned x, unsigned y, unsigned size){
  matrix[sym_matrix_coordinate(x,y,size)] = value;
}

unsigned sym_matrix_get(sym_matrix_t matrix, unsigned x, unsigned y, unsigned size){
  return matrix[sym_matrix_coordinate(x,y,size)];
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
    stars.spectralType[u] = types[fast_rng()%n_types];
    stars.subType[u] = fast_rng() % 9;
    stars.magnitude[u] = fast_rngf(-10, +20);
    stars.position.x[u] = fast_rngf(-1e5, 1e5);
    stars.position.y[u] = fast_rngf(-1e5, 1e5);
    stars.position.z[u] = fast_rngf(-3e3, 3e3);
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
  //TODO: lookup?
  return sqrtf(b+a*(1+b*c));//2 * FMA, 1 SQRT//TODO:replace sqrt with something faster http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
}

v8f starfunc_vec(v8f a, v8f b){ //symmetrical, pure
  //v8f x = _mm256_broadcast_ss(stars.subType[a]);
  //v8f y = *((v8d*)&stars.subType[b]);
  const float_t c = 1.0/0.6; //compile time evaluation.
  return _mm256_sqrt_ps(b+a*(1+b*c));//2 * FMA, 1 SQRT//TODO:replace sqrt with something faster http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
}


//TODO: use log (or something) based enumeration sort first (list of varrarray?)
//TODO: create an array of struct{distance, index}(packed) and sort using uint64_t comparisons.
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

/*
  TODO: merge with tally?
*/
void fill_matrix(star_array_t array, sym_matrix_t matrix, unsigned size)
{
  ////report(INFO, "starting to fill matrix");
  const unsigned block_count = (size+7)/8;
  unsigned row_offset = 0;
  for(unsigned u = 0; u < size; u++){ //loop over each row, vertical
    ////report(INFO, "%u -> %u [%u]", u, row_offset, size);
    v8f subType_r = _mm256_broadcast_ss(array.subType+u);
    v8f x_r = _mm256_broadcast_ss(array.position.x+u);
    v8f y_r = _mm256_broadcast_ss(array.position.y+u);
    v8f z_r = _mm256_broadcast_ss(array.position.z+u);
    for(unsigned w = 0; w < block_count-u/8; w++){//for each block in row.
      v8f subType_c = ((v8f*)array.subType)[w];
      v8f func = starfunc_vec(subType_r, subType_c);
      v8f x_c = ((v8f*)(array.position.x))[w];
      v8f y_c = ((v8f*)(array.position.y))[w];
      v8f z_c = ((v8f*)(array.position.z))[w];
      v8f dist = dist_3d_vec(
        x_r, y_r, z_r,
        x_c, y_c, z_c
      );
      ////report(INFO, "writing to %p, bounds = [%p, %p]", ((v8f*)matrix)+row_offset+w, matrix, matrix+sym_matrix_coordinate(size-1, size-1, size));
      ((v8f*)matrix)[row_offset+w] = func+dist;
    }
    row_offset += (block_count-u/8); //offset into start of this row.
  }
}

void print_matrix(sym_matrix_t theMatrix, unsigned n)//TODO: fix
{
  for(unsigned u = 0; u < triag_nr((n+7)/8)*64; u++){
    printf("%1.4e\n", theMatrix[u]);
  }
  //  unsigned i, j;
  //  printf("\nprint_matrix, n = %u:\n", n);
  //  for(i = 0 ; i < n; i++)
  //    {
  //      for(j = 0 ; j < n ; j++)
	//  printf("%1.4e " , sym_matrix_get(theMatrix, j, i, n));
  //      putchar('\n');
  //    }
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
then vectorize
then include max/min for histogram
*/

 float_t min;
 float_t max;
void create_tally_matrix(float_t *in, float_t* out, unsigned N){
  float_t up, center, left, down;

  unsigned block_count = (N+7)/8;
  unsigned in_row_index = 8*block_count+1;//always points to the first element of the current row.
  unsigned out_index = 1;
  ////report(INFO, "row_index = %u", in_row_index);

  //top left corner
  down = in[in_row_index+8*block_count];
  center = in[in_row_index];
  left = in[in_row_index-1];
  up = in[1];
  out[0] = 0.25*(fabs(center-left)+fabs(center-up)+fabs(center-down)); //left up and down, next iteration sets right.
  //top row
  for(unsigned x = 1; x < N-2; x++){
    left = center;
    center = in[in_row_index+x];
    up = in[x+1];
    out[out_index-1] += 0.25*fabs(center-left);//right to previous element.
    out[out_index] = 0.25*(fabs(center-left)+fabs(center-up)); //left and up.
    out_index += 1;
  }
  out[out_index-1] += 0.25*(fabs(center-in[in_row_index+N-1]));//right for last element
  //report(INFO, "wrote to %u first loop", out_index);
  unsigned old_out_index = 0;
  unsigned old_row_index = in_row_index;

  min = max = out[0];

  in_row_index += 8*block_count+1;
  ////report(INFO, "row_index = %u", in_row_index);
  //rest of the rows
  unsigned row, column;
  for(row = 2; row < N-1;){
    center = in[in_row_index];
    float_t up_and_left = in[old_row_index+1];//symmetry
    float_t down = in[in_row_index+1];//symmetry
    out[out_index] = 0.5*fabs(center-up_and_left)+0.25*down;
    out[old_out_index+1] += 0.25*fabs(center-up_and_left);
    set_min_max(out[old_out_index+1], &min, &max);
    out_index += 1;
    for(column = 1; column < (N-1)-row; column ++){
      left = center;
      center = in[in_row_index+column];
      up = in[old_row_index+1+column];
      out[out_index-1] += 0.25*fabs(center-left);//right to previous element.
      out[old_out_index+1+column] += 0.25*fabs(center-up);//down to the previous element.
      set_min_max(out[old_out_index+1+column], &min, &max);
      out[out_index] = 0.25*(fabs(center-left)+fabs(center-up));//up and left of this element
      out_index += 1;
    }
    out[out_index-1] += 0.25*(fabs(center-in[in_row_index+N-1-row]));//right for last element in a row.
    set_min_max(out[out_index-1], &min, &max);
    //report(INFO, "wrote to %u loop %u", out_index, row-1);
    old_out_index = out_index - (N-1)+row;
    old_row_index = in_row_index;
    row++;
    in_row_index += 8*(block_count-row/8)+1;
    ////report(INFO, "%u %u row_index = %u", N, row, in_row_index);
  }
  ////report(PASS,"Done, sampled up to %u, wrote to %u", old_row_index+column, out_index-1);
  //float_t right = in[old_row_index+column];
  //out[out_index - 1] += 0.25*fabs(center-right);
}

unsigned get_hist_index(float_t min, float_t bin_size, unsigned hist_size, float_t val){
  unsigned hist_index = (val-min)/bin_size;
  //return hist_index;
  //report(INFO, "assigning %1.4e to index %u (%f)", val, hist_index, (val-min)/bin_size);
  return hist_index >= hist_size ? hist_size-1 : hist_index;
  // for(hist_index = 0; hist_index < hist_size-1; hist_index++){//TODO: direct computation
  //   if(min + bin_size*(1+hist_index) > val){
  //     break;
  //   }
  // }
  // return hist_index;
}

hist_param_t generate_histogram(float_t *matrix, unsigned *histogram, unsigned mat_size, unsigned hist_size){
  // float_t min = FLT_MAX, max = -min;
  unsigned bounds = triag_nr(mat_size);
  // for(unsigned u = 0; u < bounds; u++){//TODO: calculate when filling matrix
  //   if(matrix[u] > max){
  //     max = matrix[u];
  //   }
  //   if(matrix[u] < min){
  //     min = matrix[u];
  //   }
  // }
  const hist_param_t param = {
    .hist_size = hist_size,
    .min = min,
    .max = max,
    .bin_size = (max-min)/(hist_size)
  };
  //report(INFO, "bounds are: %1.4e, %1.4e", min, max);
  // unsigned loc =0;
  // for(unsigned b = 0; b < bounds; b++){
  //   unsigned index = get_hist_index(min, param.bin_size, hist_size, matrix[b]);
  //   histogram[index] += 2;
  // }
  // unsigned diagonal_index = 0;
  // for(unsigned u = 0; u < mat_size; u++){
  //   unsigned index = get_hist_index(min, param.bin_size, hist_size, matrix[diagonal_index]);
  //   histogram[index] -= 1;
  //   diagonal_index += mat_size-u;
  // }
  // for(unsigned y = 0; y < mat_size; y++){
  //   unsigned index = get_hist_index(min, param.bin_size, hist_size, matrix[loc++]);
  //   histogram[index] += 1;
  //   for(unsigned x = 1; x < mat_size-y; x++){
  //     unsigned index = get_hist_index(min, param.bin_size, hist_size, matrix[loc++]);
  //     histogram[index] += 2;
  //   }
  // }
  //report(PASS, "looked at up to %u (%1.4e)", loc, matrix[0]);

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
