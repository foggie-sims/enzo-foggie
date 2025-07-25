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

  /* Return if not used */

  if (EvolveRefineRegionTime[0] <= FLOAT_UNDEFINED
      && EvolveMustRefineRegionTime[0] <= FLOAT_UNDEFINED
      && EvolveCoolingRefineRegionTime[0] <= FLOAT_UNDEFINED
      && EvolveMultiRefineRegionTime[0] <= FLOAT_UNDEFINED)
    return SUCCESS;

  if(debug1){
    printf("SetEvolveRefineRegion: setting region.\n");
    printf("SetEvolveRefineRegion: EvolveRefineRegionTime[0]:        %"PSYM"\n",EvolveRefineRegionTime[0]);
    printf("SetEvolveRefineRegion: EvolveMustRefineRegionTime[0]:    %"PSYM"\n",EvolveMustRefineRegionTime[0]);
    printf("SetEvolveRefineRegion: EvolveCoolingRefineRegionTime[0]: %"PSYM"\n",EvolveCoolingRefineRegionTime[0]);
    printf("SetEvolveRefineRegion: EvolveMultiRefineRegionTime[0]:   %"PSYM"\n",EvolveMultiRefineRegionTime[0]);
  }
  
  /* If TimeType is redshift, calculate redshift */

  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    redshift = (1.0 + InitialRedshift)/a - 1.0;
  }

  /* does time evolution for standard refinement regions */
  if((RefineRegionTimeType == 1) || (RefineRegionTimeType == 0) ){
  
    /* Find closest time step with <time */
    /* Set time=redshift if that's what we're doing. */
    if (RefineRegionTimeType == 1) {  // redshift

      time = redshift;

      /* Check to see if the current redshift is within the bounds of the given time file.  This assumes
	 that EvolveRefineRegionTime[0] is the highest redshift and the last output is the lowest (for the 
	 given track file).  */
      if(time > EvolveRefineRegionTime[0] || time < EvolveRefineRegionTime[EvolveRefineRegionNtimes-1]){
	fprintf(stderr,"SetEvolveRefineRegion ERROR: current simulation redshift is outside of range of track file redshifts!");
	my_exit(EXIT_FAILURE);
      }

      for(timestep=0; timestep<EvolveRefineRegionNtimes; timestep++){
	if( time > EvolveRefineRegionTime[timestep] ){
	  break;
	}
      }
    }else{  // code time

      /* Check to see if the current time is within the bounds of the given time file.  This assumes
	 that EvolveRefineRegionTime[0] is the earliest time and the last output is the latest time (for
	 the given track file). */
      if(time < EvolveRefineRegionTime[0] || time > EvolveRefineRegionTime[EvolveRefineRegionNtimes-1]){
	fprintf(stderr,"SetEvolveRefineRegion ERROR: current simulation time is outside of range of track file times!");
	my_exit(EXIT_FAILURE);
      }

      for(timestep=0; timestep<EvolveRefineRegionNtimes; timestep++){
	if( time < EvolveRefineRegionTime[timestep] ){
	  break;
	}
      }
    }
    timestep -= 1;
    if (timestep < 0) return SUCCESS;

    /* Set RefineRegion to EvolveRefineRegion or innermost StaticRefineRegion */

    staticRegion = 0;
    while (StaticRefineRegionLevel[staticRegion] != INT_UNDEFINED)
      staticRegion++;
    staticRegion -= 1;

    for (i = 0; i < MAX_DIMENSION; i++)
      if (staticRegion < 0) {

	/* If we're at the last timestep in our EvolveRefineRegion track, just use that;
	   otherwise, linearly interpolate between this time and the next time to avoid
	   the refinement region jumping around. */
	if(timestep == EvolveRefineRegionNtimes-1){

	  RefineRegionLeftEdge[i] = EvolveRefineRegionLeftEdge[timestep][i];
	  RefineRegionRightEdge[i] = EvolveRefineRegionRightEdge[timestep][i];

	} else {

	  RefineRegionLeftEdge[i] = EvolveRefineRegionLeftEdge[timestep][i] +
	    (time - EvolveRefineRegionTime[timestep])
	    *(EvolveRefineRegionLeftEdge[timestep+1][i]-EvolveRefineRegionLeftEdge[timestep][i])
	    / (EvolveRefineRegionTime[timestep+1] - EvolveRefineRegionTime[timestep]);

	  RefineRegionRightEdge[i] = EvolveRefineRegionRightEdge[timestep][i] +
	    (time - EvolveRefineRegionTime[timestep])
	    *(EvolveRefineRegionRightEdge[timestep+1][i]-EvolveRefineRegionRightEdge[timestep][i])
	    / (EvolveRefineRegionTime[timestep+1] - EvolveRefineRegionTime[timestep]);

	} // if(timestep == EvolveRefineRegionNtimes-1)

      } else {
	RefineRegionLeftEdge[i] = max(EvolveRefineRegionLeftEdge[timestep][i],
				      StaticRefineRegionLeftEdge[staticRegion][i]);
	RefineRegionRightEdge[i] = 
	  min(EvolveRefineRegionRightEdge[timestep][i],
	      StaticRefineRegionRightEdge[staticRegion][i]);
      } // if (staticRegion < 0) {

    if (debug1)
      fprintf(stdout, "SetEvolveRefineRegion: EvolveRegion: %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM"\n",
	      RefineRegionLeftEdge[0], RefineRegionLeftEdge[1], 
	      RefineRegionLeftEdge[2], RefineRegionRightEdge[0],
	      RefineRegionRightEdge[1], RefineRegionRightEdge[2]);

    for (int flag_method = 0; flag_method < MAX_FLAGGING_METHODS; flag_method++) {
      if( CellFlaggingMethod[flag_method] == 12 ){
	for (int dim = 0; dim < MAX_DIMENSION; dim++){
	  MustRefineRegionLeftEdge[dim]=EvolveRefineRegionLeftEdge[timestep][dim];
	  MustRefineRegionRightEdge[dim]=EvolveRefineRegionRightEdge[timestep][dim];
	}
      }
    }

  } // if((RefineRegionTimeType == 1) || (RefineRegionTimeType == 0) )



  /* does time evolution for MustRefineRegion ONLY! */
  if((MustRefineRegionTimeType == 1) || (MustRefineRegionTimeType == 0) ){
  
    /* Find closest time step with <time */
    /* Set time=redshift if that's what we're doing. */
    if (MustRefineRegionTimeType == 1) {  // redshift
      time = redshift;

      /* Check to see if the current redshift is within the bounds of the given time file.  This assumes
	 that EvolveMustRefineRegionTime[0] is the highest redshift and the last output is the lowest (for the 
	 given track file).  */
      if(time > EvolveMustRefineRegionTime[0] || time < EvolveMustRefineRegionTime[EvolveMustRefineRegionNtimes-1]){
	fprintf(stderr,"SetEvolveRefineRegion ERROR for EvolveMustRefineRegions: current simulation redshift is outside of range of track file redshifts!");
	my_exit(EXIT_FAILURE);
      }

      
      for(timestep=0; timestep<EvolveMustRefineRegionNtimes; timestep++){
	if( time > EvolveMustRefineRegionTime[timestep] ){
	  break;
	}
      }
    }else{ // code time

      /* Check to see if the current time is within the bounds of the given time file.  This assumes
	 that EvolveMustRefineRegionTime[0] is the earliest time and the last output is the latest time (for
	 the given track file). */
      if(time < EvolveMustRefineRegionTime[0] || time > EvolveMustRefineRegionTime[EvolveMustRefineRegionNtimes-1]){
	fprintf(stderr,"SetEvolveRefineRegion ERROR for EvolveMustRefineRegion: current simulation time is outside of range of track file times!");
	my_exit(EXIT_FAILURE);
      }

      for(timestep=0; timestep<EvolveMustRefineRegionNtimes; timestep++){
	if( time < EvolveMustRefineRegionTime[timestep] ){
	  break;
	}
      }
    }
    timestep -= 1;
    if (timestep < 0) return SUCCESS;

    /* Set MustRefineRegion to EvolveMustRefineRegion */

    for (i = 0; i < MAX_DIMENSION; i++){

      if(timestep == EvolveMustRefineRegionNtimes-1){
	MustRefineRegionLeftEdge[i] = EvolveMustRefineRegionLeftEdge[timestep][i];
	MustRefineRegionRightEdge[i] = EvolveMustRefineRegionRightEdge[timestep][i];
      } else {

	  MustRefineRegionLeftEdge[i] = EvolveMustRefineRegionLeftEdge[timestep][i] +
	    (time - EvolveMustRefineRegionTime[timestep])
	    *(EvolveMustRefineRegionLeftEdge[timestep+1][i]-EvolveMustRefineRegionLeftEdge[timestep][i])
	    / (EvolveMustRefineRegionTime[timestep+1] - EvolveMustRefineRegionTime[timestep]);

	  MustRefineRegionRightEdge[i] = EvolveMustRefineRegionRightEdge[timestep][i] +
	    (time - EvolveMustRefineRegionTime[timestep])
	    *(EvolveMustRefineRegionRightEdge[timestep+1][i]-EvolveMustRefineRegionRightEdge[timestep][i])
	    / (EvolveMustRefineRegionTime[timestep+1] - EvolveMustRefineRegionTime[timestep]);

      }
    } // for (i = 0; i < MAX_DIMENSION; i++)

    MustRefineRegionMinRefinementLevel = EvolveMustRefineRegionMinLevel[timestep];

    if (debug1)
      fprintf(stdout, "SetEvolveRefineRegion: EvolveMustRefineRegion: %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM"\n",
	      MustRefineRegionLeftEdge[0], MustRefineRegionLeftEdge[1], 
	      MustRefineRegionLeftEdge[2], MustRefineRegionRightEdge[0],
	      MustRefineRegionRightEdge[1], MustRefineRegionRightEdge[2],
	      MustRefineRegionMinRefinementLevel);

  } // if((MustRefineRegionTimeType == 1) || (MustRefineRegionTimeType == 0) )


  /* does time evolution for CoolingRefineRegion ONLY! */
  
  if((CoolingRefineRegionTimeType == 1) || (CoolingRefineRegionTimeType == 0) ){
  
    /* Find closest time step with <time */
    /* Set time=redshift if that's what we're doing. */
    if (CoolingRefineRegionTimeType == 1) {  // redshift
      time = redshift;

      /* Check to see if the current redshift is within the bounds of the given time file.  This assumes
	 that EvolveCoolingRefineRegionTime[0] is the highest redshift and the last output is the lowest (for the 
	 given track file).  */
      if(time > EvolveCoolingRefineRegionTime[0] || time < EvolveCoolingRefineRegionTime[EvolveCoolingRefineRegionNtimes-1]){
	fprintf(stderr,"SetEvolveRefineRegion ERROR: current simulation redshift is outside of range of track file redshifts!  (Cooling time box)");
	my_exit(EXIT_FAILURE);
      }

      for(timestep=0; timestep<EvolveCoolingRefineRegionNtimes; timestep++){
	if( time > EvolveCoolingRefineRegionTime[timestep] ){
	  break;
	}
      }
    }else{  // code time

      /* Check to see if the current time is within the bounds of the given time file.  This assumes
	 that EvolveRefineRegionTime[0] is the earliest time and the last output is the latest time (for
	 the given track file). */
      if(time < EvolveCoolingRefineRegionTime[0] || time > EvolveCoolingRefineRegionTime[EvolveCoolingRefineRegionNtimes-1]){
	fprintf(stderr,"SetEvolveRefineRegion ERROR: current simulation time is outside of range of track file times! (Cooling time box)");
	my_exit(EXIT_FAILURE);
      }

      for(timestep=0; timestep<EvolveCoolingRefineRegionNtimes; timestep++){
	if( time < EvolveCoolingRefineRegionTime[timestep] ){
	  break;
	}
      }
    }
    timestep -= 1;
    if (timestep < 0) return SUCCESS;

    /* Set CoolingRefineRegion to EvolveCoolingRefineRegion */

    for (i = 0; i < MAX_DIMENSION; i++){

      CoolingRefineRegionLeftEdge[i] = EvolveCoolingRefineRegionLeftEdge[timestep][i];
      CoolingRefineRegionRightEdge[i] = EvolveCoolingRefineRegionRightEdge[timestep][i];

      /* If we're at the last timestep in our EvolveCoolingRefineRegion track, just use that;
	 otherwise, linearly interpolate between this time and the next time to avoid
	 the refinement region jumping around. */
      if(timestep == EvolveCoolingRefineRegionNtimes-1){

	CoolingRefineRegionLeftEdge[i] = EvolveCoolingRefineRegionLeftEdge[timestep][i];
	CoolingRefineRegionRightEdge[i] = EvolveCoolingRefineRegionRightEdge[timestep][i];

      } else {

	CoolingRefineRegionLeftEdge[i] = EvolveCoolingRefineRegionLeftEdge[timestep][i] +
	  (time - EvolveCoolingRefineRegionTime[timestep])
	  *(EvolveCoolingRefineRegionLeftEdge[timestep+1][i]-EvolveCoolingRefineRegionLeftEdge[timestep][i])
	  / (EvolveCoolingRefineRegionTime[timestep+1] - EvolveCoolingRefineRegionTime[timestep]);

	CoolingRefineRegionRightEdge[i] = EvolveCoolingRefineRegionRightEdge[timestep][i] +
	  (time - EvolveCoolingRefineRegionTime[timestep])
	  *(EvolveCoolingRefineRegionRightEdge[timestep+1][i]-EvolveCoolingRefineRegionRightEdge[timestep][i])
	  / (EvolveCoolingRefineRegionTime[timestep+1] - EvolveCoolingRefineRegionTime[timestep]);

      } // if(timestep == EvolveCoolingRefineRegionNtimes-1)

    } // for (i = 0; i < MAX_DIMENSION; i++){

    if (debug1)
      fprintf(stdout, "SetEvolveRefineRegion: EvolveCoolingRefineRegion: %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM"\n",
	      CoolingRefineRegionLeftEdge[0], CoolingRefineRegionLeftEdge[1], 
	      CoolingRefineRegionLeftEdge[2], CoolingRefineRegionRightEdge[0],
	      CoolingRefineRegionRightEdge[1], CoolingRefineRegionRightEdge[2]);

  } // if((CoolingRefineRegionTimeType == 1) || (CoolingRefineRegionTimeType == 0) )
    
  /* does time evolution for MultiRefineRegion ONLY! */
  if((MultiRefineRegionTimeType == 1) || (MultiRefineRegionTimeType == 0) ){

    if (debug1){
      fprintf(stdout, "SetEvolveRefineRegion: Original MultiRefineRegion: %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM" %"ISYM" %"FSYM"\n",
      EvolveMultiRefineRegionLeftEdge[region][0][0], EvolveMultiRefineRegionLeftEdge[region][0][1],
      EvolveMultiRefineRegionLeftEdge[region][0][2], EvolveMultiRefineRegionRightEdge[region][0][0],
      EvolveMultiRefineRegionRightEdge[region][0][1], EvolveMultiRefineRegionRightEdge[region][0][2],
      EvolveMultiRefineRegionMinimumLevel[region], EvolveMultiRefineRegionMaximumLevel[region],
      EvolveMultiRefineRegionMinimumStarMass[region][0]);
    }
    
    /* Find closest time step with <time */
    /* Set time=redshift if that's what we're doing. */
    if (MultiRefineRegionTimeType == 1) {  // redshift
      time = redshift;

      /* Check to see if the current redshift is within the bounds of the given time file.  This assumes
       that EvolveMultiRefineRegionTime[0] is the highest redshift and the last output is the lowest (for the
       given track file).  */
      if(time > EvolveMultiRefineRegionTime[0] || time < EvolveMultiRefineRegionTime[NumberOfMultiRefineTimeEntries-1]){
        fprintf(stderr,"SetEvolveRefineRegion ERROR for EvolveMultiRefineRegions: current simulation redshift (%"PSYM") is outside of range of track file redshifts (minimum:%"PSYM", maximum:%"PSYM")!",time,EvolveMultiRefineRegionTime[0],EvolveMultiRefineRegionTime[NumberOfMultiRefineTimeEntries-1]);
        my_exit(EXIT_FAILURE);
      }

      for(timestep=0; timestep<NumberOfMultiRefineTimeEntries; timestep++){
        if( time > EvolveMultiRefineRegionTime[timestep] ){
          break;
        }
      }
    }else{ // code time

      /* Check to see if the current time is within the bounds of the given time file.  This assumes
       that EvolveMultiRefineRegionTime[0] is the earliest time and the last output is the latest time (for
       the given track file). */
      if(time < EvolveMultiRefineRegionTime[0] || time > EvolveMultiRefineRegionTime[NumberOfMultiRefineTimeEntries-1]){
        fprintf(stderr,"SetEvolveRefineRegion ERROR for EvolveMultiRefineRegion: current simulation time (%"PSYM") is outside of range of track file times (minimum:%"PSYM", maximum:%"PSYM")!",time,EvolveMultiRefineRegionTime[0],EvolveMultiRefineRegionTime[NumberOfMultiRefineTimeEntries-1]);
        my_exit(EXIT_FAILURE);
      }



      for(timestep=0; timestep<NumberOfMultiRefineTimeEntries; timestep++){
        if( time < EvolveMultiRefineRegionTime[timestep] ){
          break;
        }
      }
    }
    timestep -= 1;
    if (timestep < 0) return SUCCESS;
    if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
      fprintf(stderr,"SetEvolveRefineRegion: It is currently %"PSYM", which is greater than %"PSYM", so the closest timestep entry is %"ISYM".\n",time,EvolveMultiRefineRegionTime[timestep],timestep);
    }

    /* Set MultiRefineRegion to EvolveMultiRefineRegion */
    if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
      fprintf(stderr,"SetEvolveRefineRegion sees %"ISYM" static MultiRefineRegions and %"ISYM" evolving ones.\n",NumberOfStaticMultiRefineRegions,NumberOfMultiRefineTracks);
      fprintf(stderr,"SetEvolveRefineRegion says you have %"ISYM" time entries per track.\n",NumberOfMultiRefineTimeEntries);
    }

    for (region = 0; region < NumberOfMultiRefineTracks; region++){
      for (i=0; i<NumberOfMultiRefineTimeEntries; i++){
        if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"SetEvolveRefineRegion: Evolving MRR %"ISYM" has SM %"FSYM" at time entry %"ISYM".\n",region,EvolveMultiRefineRegionMinimumStarMass[region][i],i);
        }
      }
    }
    
    for (region = 0; region < NumberOfMultiRefineTracks; region++){
      if(timestep == NumberOfMultiRefineTimeEntries-1){
        MultiRefineRegionMinimumStarMass[region+NumberOfStaticMultiRefineRegions] = EvolveMultiRefineRegionMinimumStarMass[region][timestep];
        if (debug1 && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"SetEvolveRefineRegion: I set MultiRefineRegionMinimumStarMass[%"ISYM"] to %"FSYM" for timestep %"ISYM"\n.",region+NumberOfStaticMultiRefineRegions,MultiRefineRegionMinimumStarMass[region+NumberOfStaticMultiRefineRegions],timestep);
        }
        for (i = 0; i < MAX_DIMENSION; i++){
          MultiRefineRegionLeftEdge[region+NumberOfStaticMultiRefineRegions][i] = EvolveMultiRefineRegionLeftEdge[region][timestep][i];
          MultiRefineRegionRightEdge[region+NumberOfStaticMultiRefineRegions][i] = EvolveMultiRefineRegionRightEdge[region][timestep][i];
        } // for (i = 0; i < MAX_DIMENSION; i++)
      } else {
        MultiRefineRegionMinimumStarMass[region+NumberOfStaticMultiRefineRegions] = EvolveMultiRefineRegionMinimumStarMass[region][timestep]
          + (time - EvolveMultiRefineRegionTime[timestep])
          * (EvolveMultiRefineRegionMinimumStarMass[region][timestep+1]-EvolveMultiRefineRegionMinimumStarMass[region][timestep])
          / (EvolveMultiRefineRegionTime[timestep+1] - EvolveMultiRefineRegionTime[timestep]);
        if (debug && MyProcessorNumber == ROOT_PROCESSOR){
          fprintf(stderr,"SetEvolveRefineRegion: I set MultiRefineRegionMinimumStarMass[%"ISYM"] to %"FSYM" for inbtwn timestep %"ISYM".\n",region+NumberOfStaticMultiRefineRegions,MultiRefineRegionMinimumStarMass[region+NumberOfStaticMultiRefineRegions],timestep);
        }
        for (i = 0; i < MAX_DIMENSION; i++){
          MultiRefineRegionLeftEdge[region+NumberOfStaticMultiRefineRegions][i] = EvolveMultiRefineRegionLeftEdge[region][timestep][i]
            + (time - EvolveMultiRefineRegionTime[timestep])
            * (EvolveMultiRefineRegionLeftEdge[region][timestep+1][i]-EvolveMultiRefineRegionLeftEdge[region][timestep][i])
            / (EvolveMultiRefineRegionTime[timestep+1] - EvolveMultiRefineRegionTime[timestep]);

          MultiRefineRegionRightEdge[region+NumberOfStaticMultiRefineRegions][i] = EvolveMultiRefineRegionRightEdge[region][timestep][i]
            + (time - EvolveMultiRefineRegionTime[timestep])
            * (EvolveMultiRefineRegionRightEdge[region][timestep+1][i]-EvolveMultiRefineRegionRightEdge[region][timestep][i])
            / (EvolveMultiRefineRegionTime[timestep+1] - EvolveMultiRefineRegionTime[timestep]);
        } // for (i = 0; i < MAX_DIMENSION; i++)
      }  // else
        
      MultiRefineRegionMinimumLevel[region+NumberOfStaticMultiRefineRegions] = EvolveMultiRefineRegionMinimumLevel[region];
      MultiRefineRegionMaximumLevel[region+NumberOfStaticMultiRefineRegions] = EvolveMultiRefineRegionMaximumLevel[region];
      if (debug1 && MyProcessorNumber == ROOT_PROCESSOR){
        fprintf(stdout, "SetEvolveRefineRegion: EvolveMultiRefineRegion: %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM" %"ISYM" %"FSYM"\n",
          MultiRefineRegionLeftEdge[region+NumberOfStaticMultiRefineRegions][0], MultiRefineRegionLeftEdge[region+NumberOfStaticMultiRefineRegions][1],
          MultiRefineRegionLeftEdge[region+NumberOfStaticMultiRefineRegions][2], MultiRefineRegionRightEdge[region+NumberOfStaticMultiRefineRegions][0],
          MultiRefineRegionRightEdge[region+NumberOfStaticMultiRefineRegions][1], MultiRefineRegionRightEdge[region+NumberOfStaticMultiRefineRegions][2],
          MultiRefineRegionMinimumLevel[region+NumberOfStaticMultiRefineRegions], MultiRefineRegionMaximumLevel[region+NumberOfStaticMultiRefineRegions],
          MultiRefineRegionMinimumStarMass[region+NumberOfStaticMultiRefineRegions]);
      }
        
    } // for (region = 0; region < NumberOfMultiRefineTracks; region++)

  } // if((MustRefineRegionTimeType == 1) || (MustRefineRegionTimeType == 0) )

  return SUCCESS;

}
