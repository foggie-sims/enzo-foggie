/***********************************************************************
/
/  GRID CLASS (FLAGS CELL FOR REFINEMENT DEPENDING ON ITS REGION)
/
/  written by: Elizabeth Tasker
/  date:       May, 2013
/  modified1:
/
/  PURPOSE: Allows for non-cubic geometries in refined region at different levels
/
/  RETURNS: FAIL or SUCCESS
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

int grid::FlagCellsToBeRefinedByMultiRefineRegion(int level)
{

  /* Return if this grid is not on this processor
     or if multi refined regions are not being used . */
  if (MyProcessorNumber != ProcessorNumber || MultiRefineRegionGeometry[0] < 0)
    return SUCCESS;


  /* declarations */
  int i, j, k, index, dim, region, size = 1;
  float CellSize, xpos, ypos, zpos, ring_width, rad, dr;
  int LocalMaximumRefinementLevel = 0;
  int LocalMinimumRefinementLevel = 0;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION], NIter = 0;
  int NumberOfFlaggedCells = 0;
  int NRegions;


  /* Default values */
  if (MultiRefineRegionMaximumOuterLevel == INT_UNDEFINED)
    MultiRefineRegionMaximumOuterLevel = MaximumRefinementLevel;
  if (MultiRefineRegionMinimumOuterLevel == INT_UNDEFINED)
    MultiRefineRegionMinimumOuterLevel = 0;

  /* error check */
  if (FlaggingField == NULL) 
    ENZO_FAIL("Flagging Field is undefined");
    
  /* Check whether we're using evolving MultiRefine regions or not */
  if((MultiRefineRegionTimeType == 0) || (MultiRefineRegionTimeType == 1)){
    NIter = NumberOfMultiRefineTracks;
  }

  if(debug && MyProcessorNumber == ROOT_PROCESSOR){
    fprintf(stderr,"%i evolving MultiRefineRegions detected.\n",NIter);
  }

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
        for (region = 0; region < MAX_STATIC_REGIONS + NIter; region++){
          /* Check whether cell is within a given refinement region */
          if( (MultiRefineRegionLeftEdge[region][0] <= xpos) && (xpos <= MultiRefineRegionRightEdge[region][0]) &&
              (MultiRefineRegionLeftEdge[region][1] <= ypos) && (ypos <= MultiRefineRegionRightEdge[region][1]) &&
              (MultiRefineRegionLeftEdge[region][2] <= zpos) && (zpos <= MultiRefineRegionRightEdge[region][2]) ){
            /* Of those regions the cell is within, adopt refinement constraints of refine regions with maximum allowed refinement */
            if (LocalMaximumRefinementLevel < MultiRefineRegionMaximumLevel[region]){
                if(debug && MyProcessorNumber == ROOT_PROCESSOR){
                  fprintf(stderr,"Maximum cell refinement level updated from %i to %i\n",LocalMaximumRefinementLevel,MultiRefineRegionMaximumLevel[region]);
                }
                LocalMaximumRefinementLevel = MultiRefineRegionMaximumLevel[region];
            }
            if (LocalMinimumRefinementLevel < MultiRefineRegionMinimumLevel[region]){
                if(debug && MyProcessorNumber == ROOT_PROCESSOR){
                  fprintf(stderr,"Minimum cell refinement level updated from %i to %i\n",LocalMinimumRefinementLevel,MultiRefineRegionMinimumLevel[region]);
                }
                LocalMinimumRefinementLevel = MultiRefineRegionMinimumLevel[region];
            }
            NRegions ++;
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

  if (debug)
    printf("FlagCellsToBeRefinedByMultiRefineRegion: NumberOfFlaggedCells = %d (%.1f%%)\n",
	   NumberOfFlaggedCells, float(NumberOfFlaggedCells)*100.0/
	   float(size));

  return NumberOfFlaggedCells;
}
