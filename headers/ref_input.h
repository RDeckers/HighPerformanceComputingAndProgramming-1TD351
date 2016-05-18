/*****
 *  Code for examination project in High Performance Computing and Programming
 * 
 *  ref_input.c header file for reference input that can be used to
 *  verify correctness for N=7
 *
 *  To use the reference input, simply call create_ref_star_array()
 *  instead of create_random_array() when setting up the array.
 *
 *  Author: Elias Rudberg
 *
 **/

#ifndef _REF_INPUT_H
#define _REF_INPUT_H

#include "common.h"

void create_ref_star_array(star_t * stars, int size);

#endif
