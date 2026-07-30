/***********************************************************************
/
/  3D HILBERT CURVE
/
/  written by: John Wise
/  date:       April, 2010
/  modified1:
/
/  NOTES: Based on the Zoltan library
/         http://www.cs.sandia.gov/Zoltan/
/
/         Given a Cartesian coordinate [0..1], this routine will return
/         a 126-bit integer Hilbert key that progresses along a 3D
/         Hilbert curve.  42 octant levels resolve grid centers down to
/         2^-42 of the domain in global coordinates - deep-zoom grids
/         get distinct, correctly ordered keys with no bounding-box
/         renormalization (which would distort the curve and shift all
/         keys whenever the balanced grid set changes).  The previous
/         19-level double-valued key collapsed grids closer than 2^-19
/         onto identical keys and lost 4 more bits to the mantissa.
/         Keys are integer pairs so sorting and partitioning are exact;
/         longer keys are a strict refinement (prefix property), so any
/         two grids the old key already distinguished keep their order.
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

hilbert_key HilbertCurve3D(FLOAT *coord)
{

  const int MAXLEVEL = 42;

  static const char idata3d[] = {
    0,  7,  3,  4,  1,  6,  2,  5,
    0,  1,  3,  2,  7,  6,  4,  5,
    0,  3,  7,  4,  1,  2,  6,  5,
    2,  3,  5,  4,  1,  0,  6,  7,
    4,  5,  3,  2,  7,  6,  0,  1,
    4,  7,  3,  0,  5,  6,  2,  1,
    6,  7,  5,  4,  1,  0,  2,  3,
    0,  1,  7,  6,  3,  2,  4,  5,
    2,  1,  5,  6,  3,  0,  4,  7,
    6,  1,  5,  2,  7,  0,  4,  3,
    0,  7,  1,  6,  3,  4,  2,  5,
    2,  1,  3,  0,  5,  6,  4,  7,
    4,  7,  5,  6,  3,  0,  2,  1,
    4,  5,  7,  6,  3,  2,  0,  1,
    6,  1,  7,  0,  5,  2,  4,  3,
    0,  3,  1,  2,  7,  4,  6,  5,
    2,  3,  1,  0,  5,  4,  6,  7,
    6,  7,  1,  0,  5,  4,  2,  3,
    2,  5,  1,  6,  3,  4,  0,  7,
    4,  3,  7,  0,  5,  2,  6,  1,
    4,  3,  5,  2,  7,  0,  6,  1,
    6,  5,  1,  2,  7,  4,  0,  3,
    2,  5,  3,  4,  1,  6,  0,  7,
    6,  5,  7,  4,  1,  2,  0,  3};

  static const char istate3d[] = {
     1,  6,  3,  4,  2,  5,  0,  0,
     0,  7,  8,  1,  9,  4,  5,  1,
    15, 22, 23, 20,  0,  2, 19,  2,
     3, 23,  3, 15,  6, 20, 16, 22,
    11,  4, 12,  4, 20,  1, 22, 13,
    22, 12, 20, 11,  5,  0,  5, 19,
    17,  0,  6, 21,  3,  9,  6,  2,
    10,  1, 14, 13, 11,  7, 12,  7,
     8,  9,  8, 18, 14, 12, 10, 11,
    21,  8,  9,  9,  1,  6, 17,  7,
     7, 17, 15, 12, 16, 13, 10, 10,
    11, 14,  9,  5, 11, 22,  0,  8,
    18,  5, 12, 10, 19,  8, 12, 20,
     8, 13, 19,  7,  5, 13, 18,  4,
    23, 11,  7, 17, 14, 14,  6,  1,
     2, 18, 10, 15, 21, 19, 20, 15,
    16, 21, 17, 19, 16,  2,  3, 18,
     6, 10, 16, 14, 17, 23, 17, 15,
    18, 18, 21,  8, 17,  7, 13, 16,
     3,  4, 13, 16, 19, 19,  2,  5,
    16, 13, 20, 20,  4,  3, 15, 12,
     9, 21, 18, 21, 15, 14, 23, 10,
    22, 22,  6,  1, 23, 11,  4,  3,
    14, 23,  2,  9, 22, 23, 21,  0};

  int dim, level, idx;
  unsigned_int temp, state;
  unsigned long long c[3];
  hilbert_key key;
  double x;

  /* Error check */

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    if (coord[dim] < 0 || coord[dim] > 1)
      ENZO_VFAIL("Coordinates must be between 0 and 1.  coord[%d] = %f",
		 dim, coord[dim]);

  /* Convert xyz to 63-bit fixed point (bits 62..0).  Clamp exact 1.0
     to the largest double below it so the scaled value stays < 2^63. */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    x = (double) coord[dim];
    if (x >= 1.0)
      x = 1.0 - ldexp(1.0, -53);
    c[dim] = (unsigned long long) ldexp(x, 63);
  }

  /* Use state tables to convert nested quadrant's coordinates level
     by level */

  state = 0;
  key.hi = key.lo = 0;
  for (level = 0; level < MAXLEVEL; level++) {

    // Extract 3 bits at this level (bit 62-level of each coordinate)
    temp = (unsigned_int) (((c[0] >> (60-level)) & 4)
      | ((c[1] >> (61-level)) & 2)
      | ((c[2] >> (62-level)) & 1));
    idx = 8*state+temp;

    // Treat the key pair as a 128-bit shift register and shift in the
    // converted octant
    key.hi = (key.hi << 3) | (key.lo >> 61);
    key.lo = (key.lo << 3) | (unsigned long long) idata3d[idx];

    state = istate3d[idx];

  } // ENDFOR level

  return key;

}

/* (a - b) as a double: exact for nearby keys, approximate for distant
   ones.  Used only by the fuzzy-boundary heuristics in the balancers,
   which interpolate short spans along the curve; the exact integer key
   is what sorting and partitioning use. */

double HilbertKeyDiff(const hilbert_key &a, const hilbert_key &b)
{
  return ldexp((double) a.hi - (double) b.hi, 64) +
    ((double) a.lo - (double) b.lo);
}
