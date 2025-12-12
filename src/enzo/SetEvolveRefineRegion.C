/*------------------------------------------------------------------------
  SET REFINE REGION FROM EVOLVING REGION
  By John Wise

  History:
     03 May 2005 : JHW -- Created 
     26 April 2019: BWO -- updated for MustRefine and CoolingRefine regions,
                           added linear interpolation between times for all regions.
     14 December 2023: ACW - updated for MultiRefine regions
------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

void my_exit(int status);

int SetEvolveRefineRegion (FLOAT time) 
{

  int timestep, staticRegion, i, region;
  FLOAT a, dadt, redshift;
  int i=0, j=0;
  int enbctr=0;
  int NumberOfEnabledMultiRefineTracks=0;

  /* Return if we're not using evolving MultiRefineRegions */
  if (MultiRefineRegionTimeType!=0 && MultiRefineRegionTimeType!=1){
    return SUCCESS;
  }

  /* In case TimeType is redshift, calculate redshift */
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    redshift = (1.0 + InitialRedshift)/a - 1.0;
  }

  /* Find closest time entry with a time value earlier than the current simulation time */
  int closest_ind[NumberOfMultiRefineTracks];

  /* Set time=redshift if that's what we're doing. */
  if (MultiRefineRegionTimeType == 1) {  // redshift
    time = redshift;

    /* For each enabled MultiRefineRegion track, check to see if the current 
       redshift is within the bounds of the time entries. */
    for(i=0; i<NumberOfMultiRefineTracks; i++){
      closest_ind[i] = -1; // Fill with -1s for tracks that aren't enabled
      if (MRRTracks[i].Enabled == 1){
        NumberOfEnabledMultiRefineTracks++;
        if(time > MRRTracks[i].TimeEntries[0].Time || time < MRRTracks[i].TimeEntries[MRRTracks[i].NTimeEntries-1]){
          fprintf(stderr,"SetEvolveRefineRegion ERROR: current simulation redshift (%"PSYM") is outside of range of enabled MultiRefineRegion %"ISYM" (minimum:%"PSYM", maximum:%"PSYM")!",time,i,MRRTracks[i].TimeEntries[0].Time,MRRTracks[i].TimeEntries[MRRTracks[i].NTimeEntries].Time);
          my_exit(EXIT_FAILURE);
        }

        for(timestep=0; timestep<MRRTracks[i].NTimeEntries; timestep++){
          if( time > MRRTracks[i].TimeEntries[timestep] ){
            break;
          }
        }
        closest_ind[i] = timestep-1; // so that we're at the closest timestep _before_ the current time

        if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"SetEvolveRefineRegion: It is currently %"PSYM", which is greater than %"PSYM", so the closest timestep entry for enabled MultiRefineRegion track %"ISYM" is %"ISYM".\n",time,MRRTracks[0].TimeEntries[timestep],i,timestep);
        }
      } // if (MRRTracks[i].Enabled == 1)
    } // for(i=0; i<NumberOfMultiRefineTracks; i++)
  } else{ // code time

    /* For each enabled MultiRefineRegion track, check to see if the current 
       time is within the bounds of the time entries. */
    for(i=0; i<NumberOfMultiRefineTracks; i++){
      closest_ind[i] = -1; // Fill with -1s for tracks that aren't enabled
      if (MRRTracks[i].Enabled == 1){
        NumberOfEnabledMultiRefineTracks++;
        if(time < MRRTracks[i].TimeEntries[0].Time || time > MRRTracks[i].TimeEntries[MRRTracks[i].NTimeEntries-1]){
          fprintf(stderr,"SetEvolveRefineRegion ERROR: current simulation time (%"PSYM") is outside of range of enabled MultiRefineRegion %"ISYM" (minimum:%"PSYM", maximum:%"PSYM")!",time,i,MRRTracks[i].TimeEntries[0].Time,MRRTracks[i].TimeEntries[MRRTracks[i].NTimeEntries].Time);
          my_exit(EXIT_FAILURE);
        }

        for(timestep=0; timestep<MRRTracks[i].NTimeEntries; timestep++){
          if( time < MRRTracks[i].TimeEntries[timestep] ){
            break;
          }
        }
        closest_ind[i] = timestep-1; // so that we're at the closest timestep _before_ the current time
        
        if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"SetEvolveRefineRegion: It is currently %"PSYM", which is greater than %"PSYM", so the closest timestep entry for enabled MultiRefineRegion track %"ISYM" is %"ISYM".\n",time,MRRTracks[0].TimeEntries[timestep],i,timestep);
        }
      } // if (MRRTracks[i].Enabled == 1)
    } // for(i=0; i<NumberOfMultiRefineTracks; i++)
  } // if we're using code time

  /* Interpolate current position and minimum star particle mass for each MultiRefineRegion at the current time.
     Store these values, permitted refinement methods, and minimum and maximum levels for each refinement type. Note
     that minimum and maximum levels for closest time entry with a lower time than current simulation time will be
     used (i.e., these will not be interpolated) */

  if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
    fprintf(stderr,"SetEvolveRefineRegion sees %"ISYM" static MultiRefineRegions and %"ISYM" enabled evolving ones.\n",NumberOfStaticMultiRefineRegions,NumberOfEnabledMultiRefineTracks);
  }

  for (region = 0; region < NumberOfMultiRefineTracks; region++){
    if (MRRTracks[region].Enabled == 1){ // don't need to store information for tracks that aren't enabled
      timestep = closest_ind[region];
      if(timestep == MRRTracks[region].NTimeEntries-1){ // if we're at the last specified time entry, no need to interpolate
        MultiRefineRegionMinimumStarMass[enbctr+NumberOfStaticMultiRefineRegions] = MRRTracks[region].TimeEntry[timestep].MinStarMass;
        if (debug1 && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"SetEvolveRefineRegion: I set MultiRefineRegionMinimumStarMass[%"ISYM"] to %"FSYM" for timestep %"ISYM"\n.",region+NumberOfStaticMultiRefineRegions,MultiRefineRegionMinimumStarMass[region+NumberOfStaticMultiRefineRegions],timestep);
        }
        for (i = 0; i < MAX_DIMENSION; i++){
          MultiRefineRegionLeftEdge[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].TimeEntry[timestep].Pos[i];
          MultiRefineRegionRightEdge[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].TimeEntry[timestep].Pos[i+3];
        } 

        for (i=0; i<MRRTracks[region].NRefTypes; i++){
          MultiRefineRegionMinimumLevel[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].TimeEntries[timestep].MinLevels[i];
          MultiRefineRegionMaximumLevel[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].TimeEntries[timestep].MaxLevels[i];
        }
      } else { // if we're not at the last time, we need to interpolate
        MultiRefineRegionMinimumStarMass[enbctr+NumberOfStaticMultiRefineRegions] = MRRTracks[region].TimeEntry[timestep].MinStarMass
          + (time - MRRTracks[region].TimeEntry[timestep].Time)
          * (MRRTracks[region].TimeEntry[timestep+1].MinStarMass-MRRTracks[region].TimeEntry[timestep].MinStarMass)
          / (MRRTracks[region].TimeEntry[timestep+1].Time - MRRTracks[region].TimeEntry[timestep].Time);
        if (debug1 && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"SetEvolveRefineRegion: I set MultiRefineRegionMinimumStarMass[%"ISYM"] to %"FSYM" for inbtwn timestep %"ISYM".\n",region+NumberOfStaticMultiRefineRegions,MultiRefineRegionMinimumStarMass[region+NumberOfStaticMultiRefineRegions],timestep);
        }
        for (i = 0; i < MAX_DIMENSION; i++){
          MultiRefineRegionLeftEdge[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].TimeEntry[timestep].Pos[i]
            + (time - MRRTracks[region].TimeEntry[timestep].Time)
            * (MRRTracks[region].TimeEntry[timestep+1].Pos[i]-MRRTracks[region].TimeEntry[timestep].Pos[i])
            / (MRRTracks[region].TimeEntry[timestep+1].Time - MRRTracks[region].TimeEntry[timestep].Time);

          MultiRefineRegionRightEdge[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].TimeEntry[timestep].Pos[i+3]
            + (time - MRRTracks[region].TimeEntry[timestep].Time)
            * (MRRTracks[region].TimeEntry[timestep+1].Pos[i+3]-MRRTracks[region].TimeEntry[timestep].Pos[i+3])
            / (MRRTracks[region].TimeEntry[timestep+1].Time - MRRTracks[region].TimeEntry[timestep].Time);
        } // for (i = 0; i < MAX_DIMENSION; i++)

        for (i=0; i<MRRTracks[region].NRefTypes; i++){
          MultiRefineRegionMinimumLevel[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].TimeEntries[timestep].MinLevels[i];
          MultiRefineRegionMaximumLevel[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].TimeEntries[timestep].MaxLevels[i];
          MultiRefineRegionFlaggingMethod[enbctr+NumberOfStaticMultiRefineRegions][i] = MRRTracks[region].RefTypes[i];
        }
      }  // if we're not at the final timestep
      
      if (debug1 && MyProcessorNumber == ROOT_PROCESSOR){
        fprintf(stdout, "SetEvolveRefineRegion: EvolveMultiRefineRegion: %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM" %"ISYM" %"FSYM"\n",
          MultiRefineRegionLeftEdge[enbctr+NumberOfStaticMultiRefineRegions][0], MultiRefineRegionLeftEdge[enbctr+NumberOfStaticMultiRefineRegions][1],
          MultiRefineRegionLeftEdge[enbctr+NumberOfStaticMultiRefineRegions][2], MultiRefineRegionRightEdge[enbctr+NumberOfStaticMultiRefineRegions][0],
          MultiRefineRegionRightEdge[enbctr+NumberOfStaticMultiRefineRegions][1], MultiRefineRegionRightEdge[enbctr+NumberOfStaticMultiRefineRegions][2],
          MultiRefineRegionMinimumStarMass[enbctr+NumberOfStaticMultiRefineRegions]);
        for (i=0; i<MRRTracks[region].NRefTypes; i++)
          fprintf(stdout, "                                                %"ISYM": %"ISYM" %"ISYM"\n", 
          MultiRefineRegionFlaggingMethod[enbctr+NumberOfStaticMultiRefineRegions][i],MultiRefineRegionMinimumLevel[enbctr+NumberOfStaticMultiRefineRegions][i],
          MultiRefineRegionMaximumLevel[enbctr+NumberOfStaticMultiRefineRegions]);
      }
    
      enbctr++;
    } // if (MRRTracks[region].Enabled == 1)
  } // for (region = 0; region < NumberOfMultiRefineTracks; region++)

  if (enbctr != NumberOfEnabledMultiRefineTracks){
    fprintf(stderr,"SetEvolveRefineRegion lost track of an enabled MRR somewhere. There are %"ISYM" but I only wrote %"ISYM"\n.",NumberOfEnabledMultiRefineTracks,enbctr);
    my_exit(EXIT_FAILURE);
  }
  return SUCCESS;

}
