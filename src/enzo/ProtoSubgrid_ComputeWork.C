/***********************************************************************
/
/  PROTOSUBGRID: PREDICTED WORK OVER THIS CANDIDATE'S EXTENT
/
/  PURPOSE: Sums the parent's per-cell work field over the region this
/           ProtoSubgrid covers.  Used by the work cap in
/           AcceptableSubgrid (SubgridMaximumWorkFraction).
/
/           A ProtoSubgrid already indexes into the parent's arrays, so
/           nothing has to be copied when a candidate is split - the
/           pieces simply sum different sub-ranges of the same field.
/
/           Returns 0 when the feature is off (ParentWorkField NULL),
/           which makes every caller a no-op.
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ProtoSubgrid.h"

float ProtoSubgrid::ComputeWork()
{

  if (ParentWorkField == NULL)
    return 0.0;

  float work = 0.0;

  int kdim = (GridRank > 2) ? GridDimension[2] : 1;
  int jdim = (GridRank > 1) ? GridDimension[1] : 1;

  for (int k = 0; k < kdim; k++)
    for (int j = 0; j < jdim; j++) {
      int base = ((k + StartIndex[2]) * ParentDimension[1]
		  + (j + StartIndex[1])) * ParentDimension[0] + StartIndex[0];
      for (int i = 0; i < GridDimension[0]; i++)
	work += ParentWorkField[base + i];
    }

  return work;
}
