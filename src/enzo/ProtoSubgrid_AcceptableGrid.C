/***********************************************************************
/
/  PROTOSUBGRID CLASS (CHECK IF A SUBGRID NEEDS TO BE SUBDIVIDED)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:  Robert Harkness
/  date:       November, 2005
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def" 
 
/* function prototypes */

/*
   These were definitions in the original version but have been
   superseded by parameters for tuning in large-scale AMR runs

   MINIMUM_SIDE_LENGTH 4
   MAXIMUM_SIZE 2000
*/
 
int ProtoSubgrid::AcceptableSubgrid()
{
 
  /* If NumberFlagged hasn't been computed yet, then compute it. */
 
  if (NumberFlagged == INT_UNDEFINED) {
    this->ComputeSignature(0);
    NumberFlagged = 0;
    for (int i = 0; i < GridDimension[0]; i++)
      NumberFlagged += Signature[0][i];
  }

  /* Compute size and efficiency. */
 
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++){
    size *= GridDimension[dim];
    if( GridDimension[dim] >= 0.5*MAX_ANY_SINGLE_DIRECTION)
        return FALSE;
  }
 
  float efficiency = float(NumberFlagged)/float(size);

  /* For size restrictions, use the number of child cells */
 
  //size *= 8;

  /* Subgrid is acceptable if it is efficient enough or small enough,
     but if we are using multiple processors, then make sure it is
     not too big (for good load balancing). */
 
  if (size <= POW(float(MinimumSubgridEdge), GridRank))
    return TRUE;
 
  if (size > MaximumSubgridSize && NumberOfProcessors > 1)
    return FALSE;

  /* Work cap (audit T2.1/T0.3 successor).  The cell cap above cannot see
     that a grid is expensive rather than merely large, which is what
     leaves a handful of grids holding several ranks' worth of a level's
     work - and a grid cannot be split across ranks.  MaximumSubgridWork
     is one rank's even share of this level's work times
     SubgridMaximumWorkFraction, and is zero when the feature is off, so
     this test is skipped entirely by default. */

  if (MaximumSubgridWork > 0 && NumberOfProcessors > 1 &&
      ParentWorkField != NULL && this->ComputeWork() > MaximumSubgridWork)
    return FALSE;
 
  if (efficiency > MinimumEfficiency)
    return TRUE;

  /* Not acceptable yet -- keep refining. */
 
  return FALSE;
}
