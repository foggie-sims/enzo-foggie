/***********************************************************************
/
/  GRID CLASS (FLAGS CELL FOR REFINEMENT DEPENDING ON ITS REGION)
/
/  written by: Anna Wright
/  date:       Jan 3, 2024
/  modified1:
/
/  PURPOSE: flags cells for refinement if they're inside of any 
/           MultiRefineRegions with MustRefine (method 12) enabled and 
/           below the minimum refinement level for any of those regions
/            
/
/  RETURNS: number of flagged cells, or -1 on failure
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* Currently only works for cubic geometry in 3D */

int grid::FlagCellsToBeRefinedByMustRefineRegion(int level)
{

  /* Return if this grid is not on this processor
     or if multi refined regions are not being used . */
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* NOTE:  This _only_ works for GridRank = 3 (i.e. 3D calculations).  It could be
     easily extended to 1 and 2D if necessary.  */
  if(GridRank < 3){
    ENZO_FAIL("FlagCellsToBeRefinedByMustRefineRegion only works in 3D!");
  }


  /* declarations */
  int i, j, k, index, dim, region, size = 1;
  FLOAT CellSize, xpos, ypos, zpos;
  int LocalMaximumRefinementLevel = 0;
  int LocalMinimumRefinementLevel = 0;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  int NumberOfFlaggedCells = 0;
  int NRegions;
  int MRind = INT_UNDEFINED;


  /* Default values */
  if (MultiRefineRegionMaximumOuterLevel == INT_UNDEFINED)
    MultiRefineRegionMaximumOuterLevel = MaximumRefinementLevel;
  if (MultiRefineRegionMinimumOuterLevel == INT_UNDEFINED)
    MultiRefineRegionMinimumOuterLevel = 0;

  /* error check */
  if (FlaggingField == NULL) 
    ENZO_FAIL("Flagging Field is undefined");

  /* loop over dimensions - I guess this is unnecessary,
   but it's handy to have shorter names */
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Start[dim] = GridStartIndex[dim];
    End[dim]   = GridEndIndex[dim];
  }

  /* compute size */
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
    
  CellSize = FLOAT(CellWidth[0][0]);

  /* Loop over cells in grid */
  for (k = Start[2]; k <= End[2]; k++){
    for (j = Start[1]; j <= End[1]; j++){
      for (i = Start[0]; i <= End[0]; i++){

        index = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];

        xpos = GridLeftEdge[0] + (FLOAT(i-Start[0])+0.5 )*CellSize;
        ypos = GridLeftEdge[1] + (FLOAT(j-Start[1])+0.5 )*CellSize;
        zpos = GridLeftEdge[2] + (FLOAT(k-Start[2])+0.5 )*CellSize;
        
        NRegions = 0;
        /* Loop over multirefinement regions */
        for (region = 0; region < NumberOfStaticMultiRefineRegions+NumberOfEnabledMultiRefineTracks; region++){
          /* Check whether this region uses MustRefine */
          for (int reftype = 0; reftype<MAX_FLAGGING_METHODS; reftype++){
            if (MultiRefineRegionFlaggingMethod[reftype] == 12){
              MRind = reftype;
            }
          }
          if (MRind != INT_UNDEFINED){ // if it does use MustRefine ...
            /* Check whether cell is within a given refinement region */
            if( (MultiRefineRegionLeftEdge[region][0] <= xpos) && (xpos <= MultiRefineRegionRightEdge[region][0]) &&
                (MultiRefineRegionLeftEdge[region][1] <= ypos) && (ypos <= MultiRefineRegionRightEdge[region][1]) &&
                (MultiRefineRegionLeftEdge[region][2] <= zpos) && (zpos <= MultiRefineRegionRightEdge[region][2]) ){
              /* Of those regions the cell is within, adopt refinement constraints of refine regions with maximum allowed refinement */
              if (LocalMaximumRefinementLevel < MultiRefineRegionMaximumLevel[region][MRind]){
                  if(debug && MyProcessorNumber == ROOT_PROCESSOR){
                    fprintf(stderr,"Maximum cell refinement level updated from %"ISYM" to %"ISYM"\n",LocalMaximumRefinementLevel,MultiRefineRegionMaximumLevel[region][MRind]);
                  }
                  LocalMaximumRefinementLevel = MultiRefineRegionMaximumLevel[region][MRind];
              }
              if (LocalMinimumRefinementLevel < MultiRefineRegionMinimumLevel[region][MRind]){
                  if(debug && MyProcessorNumber == ROOT_PROCESSOR){
                    fprintf(stderr,"Minimum cell refinement level updated from %"ISYM" to %"ISYM"\n",LocalMinimumRefinementLevel,MultiRefineRegionMinimumLevel[region][MRind]);
                  }
                  LocalMinimumRefinementLevel = MultiRefineRegionMinimumLevel[region][MRind];
              }
              NRegions ++;
            }
          }
        }
        /* Flag for refinement if cell is below minimum level allowed */
        if ((LocalMaximumRefinementLevel > 0) || (LocalMinimumRefinementLevel > 0)){ // if cell is inside of at least one refine region
          if (level<LocalMinimumRefinementLevel){
            FlaggingField[index] = 1;
          }
        else
          if (level<MultiRefineRegionMinimumOuterLevel){
            FlaggingField[index] = 1;
          }
        }
      }
    }
  }

  /* Count number of flagged Cells. */
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  return NumberOfFlaggedCells;
}
