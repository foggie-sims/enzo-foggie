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
/           region if it is lower than minimum for timestep.
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
  FLOAT MRRLeftEdge[MAX_DIMENSION], MRRRightEdge[MAX_DIMENSION], redshift, ctime, a, dadt;
  float MRRMinimumStarMass;

/* If the current grid is inside of any multirefine regions with a stellar mass threshold lower than the default for this timestep,
   set StarMakerMinimumMass to the stellar mass threshold for the multirefine region */

  /* Take care of any static multirefine regions first */
  for (region = 0; region < NumberOfStaticMultiRefineRegions; region++){
    /* Does this region have a set minimum star mass? Is it lower than the current minimum for this timestep? */
    if (MultiRefineRegionMinimumStarMass[region]>0 && StarMakerMinimumMass>MultiRefineRegionMinimumStarMass[region]){ 
      /* Is the current grid within this region? */
      if ((GridLeftEdge[0] <= MultiRefineRegionRightEdge[region][0]) && (GridRightEdge[0] >= MultiRefineRegionLeftEdge[region][0]) &&
          (GridLeftEdge[1] <= MultiRefineRegionRightEdge[region][1]) && (GridRightEdge[1] >= MultiRefineRegionLeftEdge[region][1]) &&
          (GridLeftEdge[2] <= MultiRefineRegionRightEdge[region][2]) && (GridRightEdge[2] >= MultiRefineRegionLeftEdge[region][2])){
          StarMakerMinimumMass = MultiRefineRegionMinimumStarMass[region];
      }
    }
  }
  
  /* Now take care of any evolving multirefine regions. We must recalculate the position and stellar mass threshold for each
     evolving multirefine region in order to account for the fact that the current time will be different for different levels. */
     
  /* Find closest time entry with a time value earlier than the current simulation time */
  int closest_ind[NumberOfMultiRefineTracks];

  if (MultiRefineRegionTimeType == 1) {  // redshift
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    redshift = (1.0 + InitialRedshift)/a - 1.0;
    ctime = redshift;
    /* For each enabled MultiRefineRegion track, check to see if the current 
       redshift is within the bounds of the time entries. */
    for(i=0; i<NumberOfMultiRefineTracks; i++){
      closest_ind[i] = -1; // Fill with -1s for tracks that aren't enabled
      if (MRRTracks[i].Enabled == 1){
        if(ctime > MRRTracks[i].TimeEntries[0].Time || ctime < MRRTracks[i].TimeEntries[MRRTracks[i].NTimeEntries-1]){
          fprintf(stderr,"Grid_SetMinimumStarMass ERROR: current simulation redshift (%"PSYM") is outside of range of enabled MultiRefineRegion %"ISYM" (minimum:%"PSYM", maximum:%"PSYM")!",ctime,i,MRRTracks[i].TimeEntries[0].Time,MRRTracks[i].TimeEntries[MRRTracks[i].NTimeEntries].Time);
          my_exit(EXIT_FAILURE);
        }

        for(timestep=0; timestep<MRRTracks[i].NTimeEntries; timestep++){
          if( ctime > MRRTracks[i].TimeEntries[timestep] ){
            break;
          }
        }
        closest_ind[i] = timestep-1; // so that we're at the closest timestep _before_ the current time

        if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"Grid_SetMinimumStarMass: It is currently %"PSYM", which is greater than %"PSYM", so the closest timestep entry for enabled MultiRefineRegion track %"ISYM" is %"ISYM".\n",ctime,MRRTracks[0].TimeEntries[timestep],i,timestep);
        }
      } // if (MRRTracks[i].Enabled == 1)
    } // for(i=0; i<NumberOfMultiRefineTracks; i++)
  } else{ // code time
    ctime = Time;
    /* For each enabled MultiRefineRegion track, check to see if the current 
       time is within the bounds of the time entries. */
    for(i=0; i<NumberOfMultiRefineTracks; i++){
      closest_ind[i] = -1; // Fill with -1s for tracks that aren't enabled
      if (MRRTracks[i].Enabled == 1){
        if(ctime < MRRTracks[i].TimeEntries[0].Time || ctime > MRRTracks[i].TimeEntries[MRRTracks[i].NTimeEntries-1]){
          fprintf(stderr,"Grid_SetMinimumStarMass ERROR: current simulation time (%"PSYM") is outside of range of enabled MultiRefineRegion %"ISYM" (minimum:%"PSYM", maximum:%"PSYM")!",ctime,i,MRRTracks[i].TimeEntries[0].Time,MRRTracks[i].TimeEntries[MRRTracks[i].NTimeEntries].Time);
          my_exit(EXIT_FAILURE);
        }

        for(timestep=0; timestep<MRRTracks[i].NTimeEntries; timestep++){
          if( ctime < MRRTracks[i].TimeEntries[timestep] ){
            break;
          }
        }
        closest_ind[i] = timestep-1; // so that we're at the closest timestep _before_ the current time
        
        if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"Grid_SetMinimumStarMass: It is currently %"PSYM", which is greater than %"PSYM", so the closest timestep entry for enabled MultiRefineRegion track %"ISYM" is %"ISYM".\n",ctime,MRRTracks[0].TimeEntries[timestep],i,timestep);
        }
      } // if (MRRTracks[i].Enabled == 1)
    } // for(i=0; i<NumberOfMultiRefineTracks; i++)
  } // if we're using code time

  /* Interpolate the position of each enabled MultiRefineRegion and check whether or not this grid overlaps with any of them.
     for MRRs with overlap, calculate the current value of the MinimumStarMass and adopt the minimum allowed for this grid */
  for (region = 0; region < NumberOfMultiRefineTracks; region++){
    if (MRRTracks[region].Enabled == 1){ // don't need to calculate for tracks that aren't enabled
      timestep = closest_ind[region];
      /* Recalculate MRR position and minimum star particle mass */
      if(timestep == MRRTracks[region].NTimeEntries-1){ // if we're at the last specified time entry, no need to interpolate
        for (i = 0; i < MAX_DIMENSION; i++){
          MRRLeftEdge[i] = MRRTracks[region].TimeEntry[timestep].Pos[i];
          MRRRightEdge[i] = MRRTracks[region].TimeEntry[timestep].Pos[i+3];
        }  // for (i = 0; i < MAX_DIMENSION; i++)
        MRRMinimumStarMass = MRRTracks[region].TimeEntry[timestep].MinStarMass;
      } else { // if we're not at the last time, we need to interpolate
        MRRMinimumStarMass = MRRTracks[region].TimeEntry[timestep].MinStarMass
          + (time - MRRTracks[region].TimeEntry[timestep].Time)
          * (MRRTracks[region].TimeEntry[timestep+1].MinStarMass-MRRTracks[region].TimeEntry[timestep].MinStarMass)
          / (MRRTracks[region].TimeEntry[timestep+1].Time - MRRTracks[region].TimeEntry[timestep].Time);
        for (i = 0; i < MAX_DIMENSION; i++){
          MRRLeftEdge[i] = MRRTracks[region].TimeEntry[timestep].Pos[i]
            + (time - MRRTracks[region].TimeEntry[timestep].Time)
            * (MRRTracks[region].TimeEntry[timestep+1].Pos[i]-MRRTracks[region].TimeEntry[timestep].Pos[i])
            / (MRRTracks[region].TimeEntry[timestep+1].Time - MRRTracks[region].TimeEntry[timestep].Time);
          MRRRightEdge[i] = MRRTracks[region].TimeEntry[timestep].Pos[i+3]
            + (time - MRRTracks[region].TimeEntry[timestep].Time)
            * (MRRTracks[region].TimeEntry[timestep+1].Pos[i+3]-MRRTracks[region].TimeEntry[timestep].Pos[i+3])
            / (MRRTracks[region].TimeEntry[timestep+1].Time - MRRTracks[region].TimeEntry[timestep].Time);
        } // for (i = 0; i < MAX_DIMENSION; i++)
      }  // if we're not at the final timestep
      
      /* If our grid overlaps with this MRR, check to see if this MRR's minimum star mass is lower than the current one on this grid
         If so, set this grid's minimum star particle mass to that of this MRR */
      if(((GridRightEdge[0] > MRRLeftEdge[0]) && (GridLeftEdge[0] < MRRRightEdge[0]))
        && ((GridRightEdge[1] > MRRLeftEdge[1]) && (GridLeftEdge[1] < MRRRightEdge[1]))
        && ((GridRightEdge[2] > MRRLeftEdge[2]) && (GridLeftEdge[2] < MRRRightEdge[2]))){
        if (StarMakerMinimumMass>MRRMinimumStarMass){
          StarMakerMinimumMass = MRRMinimumStarMass;
        }
      }

      if (debug1 && MyProcessorNumber == ROOT_PROCESSOR){
        if (MultiRefineRegionDefaultStarMass>StarMakerMinimumMass){
          fprintf(stderr,"Grid_SetMinimumStarMass: I set StarMakerMinimumMass to %"FSYM".\n",StarMakerMinimumMass);
        }
      }
    } // if (MRRTracks[region].Enabled == 1)
  } // for (region = 0; region < NumberOfMultiRefineTracks; region++)

  return SUCCESS;

}