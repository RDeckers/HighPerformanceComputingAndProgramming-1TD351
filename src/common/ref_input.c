/*****
 *  Code for examination project in High Performance Computing and Programming
 *
 *  ref_input.c reference input that can be used to verify correctness
 *  for N=7
 *
 *  To use the reference input, simply call create_ref_star_array()
 *  instead of create_random_array() when setting up the array.
 *
 *  Author: Elias Rudberg
 *
 **/

#include "ref_input.h"
#include <string.h>
#include <utilities/logging.h>

static void assign(
				 star_array_t stars,
				 unsigned index,
		     char spectralType,
		     unsigned short subType,
		     float_t magnitude,
		     const char* designation,
		     float_t pos_x,
		     float_t pos_y,
		     float_t pos_z)
{
	//report(INFO, "adding reference star %u (%s)", index, designation);
  stars.index[index] = index;
	//report(PASS, "added index");
  stars.spectralType[index] = spectralType;
	//report(PASS, "added spectralType");
  stars.subType[index] = subType;
	//report(PASS, "added subType");
  stars.magnitude[index] = magnitude;
	//report(PASS, "added magnitude");
  //strcpy(star.designation, designation);
  stars.position.x[index] = pos_x;
	//report(PASS, "added pos_x");
  stars.position.y[index] = pos_y;
	//report(PASS, "added pos_y");
  stars.position.z[index] = pos_z;
	//report(PASS, "added pos_z");
}

void create_ref_star_array(star_array_t stars, int size)
{
  if(size != 7) {
    report(FAIL, "create_ref_star_array() only works for size=7 case.");
    exit(-1);
  }
  assign(stars, 0, 'K', 9,  -2.86561203, "K9.0",    -13380.09863,    -21690.53516,     -2935.21948);
  assign(stars, 1, 'F', 5,   6.27077770, "F5.1",      4111.79736,     37506.52734,       224.89719);
  assign(stars, 2, 'L', 3,  -8.00868320, "L3.2",    -78079.12500,     70167.44531,      -976.68732);
  assign(stars, 3, 'M', 1,  -4.14477539, "M1.3",    -43289.43359,     50282.46094,     -1989.34229);
  assign(stars, 4, 'K', 3,  -7.53354073, "K3.4",    -47956.98438,    -36407.84375,      -198.44850);
  assign(stars, 5, 'O', 9,  14.65140247, "O9.5",    -78079.49219,    -92333.16406,     -1720.43555);
  assign(stars, 6, 'M', 6,  -2.00733948, "M6.6",     32551.70508,     71175.64062,      2723.72803);
}
