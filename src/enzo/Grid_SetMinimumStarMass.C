/***********************************************************************
/
/  GRID CLASS (COMPUTE THE MINIMUM STAR PARTICLE MASS)
/
/  written by: Anna Wright
/  date:       January 19, 2024
/  modified1:
/
/  PURPOSE: Checks to see if grid is inside of a multirefine region and,
/           if so, sets minimum star particle mass to minimum for that 
/           region.
/  RETURNS: Success or Fail
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "Fluxes.h"
#include "typedefs.h"
#include "global_data.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int grid::SetMinimumStarMass(){
 /* Return if this grid is not on this processor or if multirefine regions with spatially 
    varying minimum stellar masses are not being used . */

  if (MyProcessorNumber != ProcessorNumber || MultiRefineRegionSpatiallyVaryingStarMass != 1)
    return SUCCESS;

  /* Declarations */
  int region, i, timestep;
  float MRRLeftEdge[MAX_DIMENSION], MRRRightEdge[MAX_DIMENSION], MRRMinimumStarMass, redshift, ctime, a, dadt;
  FILE *MRRFile;
  char MRRFilename[32];

/* If the current grid is inside of any multirefine regions with a stellar mass threshold lower than the default for this timestep, */
/* set StarMakerMinimumMass to the stellar mass threshold for the multirefine region */
  if (debug){ 
    sprintf(MRRFilename, "MRRGridPositions_%d.txt", MyProcessorNumber);
    MRRFile = fopen(MRRFilename,"a");
    fprintf(MRRFile,"Grid %"ISYM": %"PSYM",%"PSYM",%"PSYM"; %"PSYM",%"PSYM",%"PSYM" at %"FSYM".\n",
        ID,GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2],GridRightEdge[0],GridRightEdge[1],GridRightEdge[2],Time);
  }

  /* Take care of any static multirefine regions first */
  for (region = 0; region < NumberOfStaticMultiRefineRegions; region++){
    // Does this region have a set minimum star mass? Is it lower than the current minimum for this timestep?
    if (MultiRefineRegionMinimumStarMass[region]>0 && StarMakerMinimumMass>MultiRefineRegionMinimumStarMass[region]){ 
      // Is the current grid within this region?
      if ((GridLeftEdge[0] <= MultiRefineRegionRightEdge[region][0]) && (GridRightEdge[0] >= MultiRefineRegionLeftEdge[region][0]) &&
          (GridLeftEdge[1] <= MultiRefineRegionRightEdge[region][1]) && (GridRightEdge[1] >= MultiRefineRegionLeftEdge[region][1]) &&
          (GridLeftEdge[2] <= MultiRefineRegionRightEdge[region][2]) && (GridRightEdge[2] >= MultiRefineRegionLeftEdge[region][2])){
          StarMakerMinimumMass = MultiRefineRegionMinimumStarMass[region];
          if (debug){
            fprintf(stderr,"Minimum Stellar Mass updated to %"FSYM" because of overlap with Static MRR %"ISYM".\n",StarMakerMinimumMass,region);
          }
      }
    }
  }
  
  /* Now take care of any evolving multirefine regions. We must recalculate the position and stellar mass threshold for each */
  /* evolving multirefine region in order to account for the fact that the current time will be different for different levels */
  if (MultiRefineRegionTimeType == 1) {  // redshift
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    redshift = (1 + InitialRedshift)/a - 1;
    ctime = redshift;
  } else{ // code time
    ctime = Time;
  }

  if (MultiRefineRegionTimeType == 1) {  // redshift
    for(timestep=0; timestep<NumberOfMultiRefineTimeEntries; timestep++){
      if( ctime > EvolveMultiRefineRegionTime[timestep] ){
        break;
      }
    }
  } else{ // code time
    for(timestep=0; timestep<EvolveRefineRegionNtimes; timestep++){
	    if( ctime < EvolveRefineRegionTime[timestep] ){
	      break;
      }
    }
  }

  timestep -= 1;
  if (timestep < 0) return SUCCESS;

  for (region = 0; region < NumberOfMultiRefineTracks; region++){
    if (MultiRefineRegionMinimumStarMass[NumberOfStaticMultiRefineRegions+region]>0){ // Does this region have a set minimum star mass?
      if(timestep == NumberOfMultiRefineTimeEntries-1){
        MRRMinimumStarMass = EvolveMultiRefineRegionMinimumStarMass[region][timestep];
        for (i = 0; i < MAX_DIMENSION; i++){
          MRRLeftEdge[i] = EvolveMultiRefineRegionLeftEdge[region][timestep][i];
          MRRRightEdge[i] = EvolveMultiRefineRegionRightEdge[region][timestep][i];
        }
      } else {
          MRRMinimumStarMass = EvolveMultiRefineRegionMinimumStarMass[region][timestep]
          + (ctime - EvolveMultiRefineRegionTime[timestep])
          * (EvolveMultiRefineRegionMinimumStarMass[region][timestep+1]-EvolveMultiRefineRegionMinimumStarMass[region][timestep])
          / (EvolveMultiRefineRegionTime[timestep+1] - EvolveMultiRefineRegionTime[timestep]);
        for (i = 0; i < MAX_DIMENSION; i++){
          MRRLeftEdge[i] = EvolveMultiRefineRegionLeftEdge[region][timestep][i]
            + (ctime - EvolveMultiRefineRegionTime[timestep])
            * (EvolveMultiRefineRegionLeftEdge[region][timestep+1][i]-EvolveMultiRefineRegionLeftEdge[region][timestep][i])
            / (EvolveMultiRefineRegionTime[timestep+1] - EvolveMultiRefineRegionTime[timestep]);
          MRRRightEdge[i] = EvolveMultiRefineRegionRightEdge[region][timestep][i]
            + (ctime - EvolveMultiRefineRegionTime[timestep])
            * (EvolveMultiRefineRegionRightEdge[region][timestep+1][i]-EvolveMultiRefineRegionRightEdge[region][timestep][i])
            / (EvolveMultiRefineRegionTime[timestep+1] - EvolveMultiRefineRegionTime[timestep]);
        } // for (i = 0; i < MAX_DIMENSION; i++)
      } // if(timestep == NumberOfMultiRefineTimeEntries-1)
      if (debug){
        fprintf(MRRFile,"MRR %"ISYM": %"PSYM",%"PSYM",%"PSYM"; %"PSYM",%"PSYM",%"PSYM" at %"FSYM".\n",
        region,MRRLeftEdge[0],MRRLeftEdge[1],MRRLeftEdge[2],MRRRightEdge[0],MRRRightEdge[1],MRRRightEdge[2],ctime);
      }
      if ((GridLeftEdge[0] <= MRRRightEdge[0]) && (GridRightEdge[0] >= MRRLeftEdge[0]) &&
          (GridLeftEdge[1] <= MRRRightEdge[1]) && (GridRightEdge[1] >= MRRLeftEdge[1]) &&
          (GridLeftEdge[2] <= MRRRightEdge[2]) && (GridRightEdge[2] >= MRRLeftEdge[2])){
        if (StarMakerMinimumMass>MRRMinimumStarMass){
          StarMakerMinimumMass = MRRMinimumStarMass;
          if (debug){
            fprintf(stderr,"Minimum Stellar Mass updated to %"FSYM" because of overlap with Evolving MRR %"ISYM".\n",StarMakerMinimumMass,region);
          }
        }
      } // if region in MRR
    } // if region is using a specific minimum star mass
  } // for (region = 0; region < NumberOfMultiRefineTracks; region++)

  if (debug){ 
    fprintf(MRRFile,"Final SM: %"FSYM".\n",StarMakerMinimumMass);
    fclose(MRRFile);
  }

  return SUCCESS;

}