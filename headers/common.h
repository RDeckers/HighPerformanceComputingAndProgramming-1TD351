/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  common.h common definitions header file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#ifndef _COMMON_H
#define _COMMON_H

typedef float float_t;
#define NUM_HIST_BOXES 10

typedef struct star{
  int index;                   // counting index
  char spectralType;           // random: O, B, A, F, G, K, M, L, T
  unsigned short subType;      // random: 0-9
  float_t magnitude;           // random: (-10, +20)
  char designation[9]; 	       // sprintf("%c%d.%d", spectralType, subType, index)
  struct {
    float_t x, y, z;          // x & y random in (-1e5, 1e5), z random in (-3e3, 3e3)
  } position;
} star_t;

typedef struct hist_params{
  int hist_size;
  float_t min, max, bin_size;
} hist_param_t;

#endif
