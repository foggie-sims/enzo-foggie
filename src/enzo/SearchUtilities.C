/***********************************************************************
/
/  BASIC SEARCHING UTILITIES
/
/  written by: John Wise
/  date:       August, 2010
/  modified:   
/
/  PURPOSE:
/
************************************************************************/
// std::lower_bound has too much overhead

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int search_lower_bound(int *arr, int value, int low, int high, 
		       int total)
{
  int mid, width;
  if (high < low)
    ENZO_VFAIL("high (%d) < low (%d) when searching for lower bound.",
	       high, low);
  width = high-low;
  mid = low + width/2;
  // First catch here if it's the last recursive call
//  if (debug)
//    printf("low high mid :: width value = %d (%d) %d (%d) %d (%d) :: %d %d\n", 
//	   low, arr[low], high, arr[high], mid, arr[mid], width, value);
  if (width <= 1) {
    if (mid < total-1)
      if (value >= arr[mid] && value < arr[mid+1])
	return mid;
    if (mid > 0)
      if (value >= arr[mid-1] && value < arr[mid])
	return mid-1;
  } // ENDIF width <= 1

  /* Clamp instead of recursing into an empty interval.
   *
   * A value outside [arr[0], arr[total]) reaches neither return above: the
   * first needs value >= arr[mid], the second needs mid > 0.  Control then
   * fell through to the recursion, and with mid == low that is
   * (low, low-1) -- high < low, and the ENZO_VFAIL at the top of the next
   * call aborts the run.  This killed halo5348/25Mpc_DM_256-L3-gas at
   * RD0094 (z = 2.45) with "high (-1) < low (0)", reached from
   * CommunicationTransferParticles during a level-0 RebuildHierarchy, and
   * left 25Mpc_DM_512-L0-gas restarting from DD0020 and re-crashing four
   * times over three days.
   *
   * The trigger is a particle a fraction of a cell outside the periodic
   * domain -- 0.066 of a root cell in that run -- giving a negative
   * CenterIndex.  Such particles exist because the only periodic wrap left
   * in the particle path is the one in CommunicationTransferParticles, and
   * a dump can be written between a position update and the next level-0
   * rebuild.  That is a separate bug; this function should not abort on it
   * either way.
   *
   * Clamping is what the caller wants: a particle just past the left edge
   * belongs to the first processor slab, one just past the right edge to
   * the last.  Grid_CommunicationTransferParticlesOpt already applies
   * min(GridPosition, Layout-1) to this return value; the top-level
   * CommunicationTransferParticles does not, so the upper clamp is capped
   * at total-1 here rather than returned as `high`, which would index
   * GridMap out of bounds.
   */
  if (arr[mid] > value) {
    if (mid <= low)
      return low;
    return search_lower_bound(arr, value, low, mid-1, total);
  } else if (arr[mid] < value) {
    if (mid >= high)
      return (high < total) ? high : total-1;
    return search_lower_bound(arr, value, mid+1, high, total);
  } else
    /* Exact hit.  arr[total] is the sentinel (TopGridDims), not a slab, so
     * a value landing exactly on it must clamp too -- returning `total`
     * here would index Layout[] one past the end. */
    return (mid < total) ? mid : total-1;  // found it exactly.
}

/* float version */

int search_lower_bound(float *arr, float value, int low, int high, 
		       int total)
{
  int mid, width;
  if (high < low)
    ENZO_VFAIL("high (%d) < low (%d) when searching for lower bound.",
	       high, low);
  width = high-low;
  mid = low + width/2;
  // First catch here if it's the last recursive call
//  if (debug)
//    printf("low high mid :: width value = %d (%d) %d (%d) %d (%d) :: %d %d\n", 
//	   low, arr[low], high, arr[high], mid, arr[mid], width, value);
  if (width <= 1) {
    if (mid < total-1)
      if (value >= arr[mid] && value < arr[mid+1])
	return mid;
    if (mid > 0)
      if (value >= arr[mid-1] && value < arr[mid])
	return mid-1;
  } // ENDIF width <= 1

  /* Same clamp as the int version above -- see the comment there for why.
   * Kept identical rather than factored out: the two differ only in the
   * type of arr and value, and a shared helper would need a template or a
   * macro for no benefit at this size. */
  if (arr[mid] > value) {
    if (mid <= low)
      return low;
    return search_lower_bound(arr, value, low, mid-1, total);
  } else if (arr[mid] < value) {
    if (mid >= high)
      return (high < total) ? high : total-1;
    return search_lower_bound(arr, value, mid+1, high, total);
  } else
    /* Exact hit.  arr[total] is the sentinel (TopGridDims), not a slab, so
     * a value landing exactly on it must clamp too -- returning `total`
     * here would index Layout[] one past the end. */
    return (mid < total) ? mid : total-1;  // found it exactly.
}
