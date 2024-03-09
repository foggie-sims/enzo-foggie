/*------------------------------------------------------------------------

  READ EVOLVING REFINE REGION FILE
  By John Wise

  File format: 
    EvolveRefineRegion:        (time or redshift) x_left y_left z_left x_right y_right z_right
    EvolveMustRefineRegion:    (time or redshift) x_left y_left z_left x_right y_right z_right min_level
    EvolveCoolingRefineRegion: (time or redshift) x_left y_left z_left x_right y_right z_right min_level

  Notes:  A "RefineRegion" adjusts the boundaries of the rectangular solid area within which 
          refinement based on any refinement criteria is allowed to occur.  A "MustRefineRegion"
	  forces refinement to min_level, and within that region additional refinement can occur,
	  up to MaximumRefinementLevel, based on any refinement criteria.  A CoolingRefineRegion
	  restricts refinement based on the cooling time to only a subvolume of the simulation, 
	  and is useful because the cooling time criterion can result in very aggressive refinement
	  in some circumstances.  **IN PRINCIPLE ALL OF THESE CAN BE USED SIMULTANEOUSLY**

  History:
     03 May 2005 : JHW -- Created
     26 April 2019: BWO -- Added evolving MustRefine and CoolingRefine regions.

------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int ReadEvolveRefineFile(void)
{

  FILE *fptr;

  /* Read in data file for an evolving RefineRegion */
  if((RefineRegionTimeType == 0) || (RefineRegionTimeType == 1)){
    
    if ((fptr = fopen(RefineRegionFile, "r")) == NULL) {
      fprintf(stderr, "Error opening refine region file %s.\n", RefineRegionFile);
      return FAIL;
    }

    char line[MAX_LINE_LENGTH];
    int nret, i=0;
    EvolveRefineRegionNtimes=0;
    while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL)){
      nret = sscanf(line, "%"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM,
		    &(EvolveRefineRegionTime[i]),
		    &(EvolveRefineRegionLeftEdge[i][0]),
		    &(EvolveRefineRegionLeftEdge[i][1]),
		    &(EvolveRefineRegionLeftEdge[i][2]),
		    &(EvolveRefineRegionRightEdge[i][0]),
		    &(EvolveRefineRegionRightEdge[i][1]),
		    &(EvolveRefineRegionRightEdge[i][2])); 
      if( nret != 7 ){
	fprintf(stderr,"WARNING: ReadEvolveRefineFile cannot interpret line %s",line);
	continue;
      }
      i++;
      EvolveRefineRegionNtimes++;
    }
    
    fclose(fptr);

    // Error check - are there too many input times?
    if(EvolveRefineRegionNtimes > MAX_REFINE_REGIONS){
      fprintf(stderr, "Too many EvolveRefineRegion times in your file!\nIncrease MAX_REFINE_REGIONS in macros_and_parameters.h!\n");
      return FAIL;
    }
    
  }

  /* Read in data file for an evolving MustRefineRegion  */
  if((MustRefineRegionTimeType == 0) || (MustRefineRegionTimeType == 1)){

    if ((fptr = fopen(MustRefineRegionFile, "r")) == NULL) {
      fprintf(stderr, "Error opening MustRefine region file %s.\n", MustRefineRegionFile);
      return FAIL;
    }

    char line[MAX_LINE_LENGTH];
    int nret, i=0;
    EvolveMustRefineRegionNtimes=0;
    while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL)){
      nret = sscanf(line, "%"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM,
		    &(EvolveMustRefineRegionTime[i]),
		    &(EvolveMustRefineRegionLeftEdge[i][0]),
		    &(EvolveMustRefineRegionLeftEdge[i][1]),
		    &(EvolveMustRefineRegionLeftEdge[i][2]),
		    &(EvolveMustRefineRegionRightEdge[i][0]),
		    &(EvolveMustRefineRegionRightEdge[i][1]),
		    &(EvolveMustRefineRegionRightEdge[i][2]),
      		    &(EvolveMustRefineRegionMinLevel[i]));
      if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
         fprintf(stderr,"Here is the line (MustRefineRegion): %s \n",line);
         fprintf(stderr,". . . and here is the value (MustRefineRegion): %i \n",EvolveMustRefineRegionMinLevel[i]);
         } 
      if( nret != 8 ){
	fprintf(stderr,"WARNING: ReadEvolveRefineFile (MustRefineRegion) cannot interpret line %s",line);
	continue;
      }
      i++;
      EvolveMustRefineRegionNtimes++;
    }

    fclose(fptr);

    // Error check - are there too many input times?
    if(EvolveMustRefineRegionNtimes > MAX_REFINE_REGIONS){
      fprintf(stderr, "Too many EvolveMustRefineRegion times in your file!\nIncrease MAX_REFINE_REGIONS in macros_and_parameters.h!\n");
      return FAIL;
    }

    /* print out debugging information for evolving MustRefine region */
    if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){

      printf("ReadEvolveMustRefineFile: I have a MustRefineRegion with TimeType %"ISYM" \n",
	     MustRefineRegionTimeType);
      
      printf("ReadEvolveRefineFile: And here is what I think my times, edges, and minimum levels are:\n");

      for(int i=0; i<EvolveMustRefineRegionNtimes; i++){
	printf("ReadEvolveRefineFile (MustRefineRegion): %"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM"\n",
	       EvolveMustRefineRegionTime[i],
	       EvolveMustRefineRegionLeftEdge[i][0],
	       EvolveMustRefineRegionLeftEdge[i][1],
	       EvolveMustRefineRegionLeftEdge[i][2],
	       EvolveMustRefineRegionRightEdge[i][0],
	       EvolveMustRefineRegionRightEdge[i][1],
	       EvolveMustRefineRegionRightEdge[i][2],
	       EvolveMustRefineRegionMinLevel[i]);
      } // for loop

      fflush(stdout);

    } // if(debug1 && MyProcessorNumber == ROOT_PROCESSOR)

  } // if((MustRefineRegionTimeType == 0) || (MustRefineRegionTimeType == 1))


  /* Read in CoolingRefineRegion information.  Note that this requires a file that is 
     EXACTLY like the MustRefineRegion file, which includes a level as the last entry in 
     each row of the file.  This is NOT USED but must be there, and is done because we often use
     the same file for both criteria.  */
  if((CoolingRefineRegionTimeType == 0) || (CoolingRefineRegionTimeType == 1)){

    if ((fptr = fopen(CoolingRefineRegionFile, "r")) == NULL) {
      fprintf(stderr, "Error opening CoolingRefine region file %s.\n", CoolingRefineRegionFile);
      return FAIL;
    }

    char line[MAX_LINE_LENGTH];
    int nret, i=0, dummy;
    EvolveCoolingRefineRegionNtimes=0;
    while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL)){
      nret = sscanf(line, "%"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM,
		    &(EvolveCoolingRefineRegionTime[i]),
		    &(EvolveCoolingRefineRegionLeftEdge[i][0]),
		    &(EvolveCoolingRefineRegionLeftEdge[i][1]),
		    &(EvolveCoolingRefineRegionLeftEdge[i][2]),
		    &(EvolveCoolingRefineRegionRightEdge[i][0]),
		    &(EvolveCoolingRefineRegionRightEdge[i][1]),
		    &(EvolveCoolingRefineRegionRightEdge[i][2]),
      		    &dummy);
      if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
         fprintf(stderr,"Here is the line (CoolingRefineRegion): %s \n",line);
         fprintf(stderr,". . . and here is the value (CoolingRefineRegion): %i \n",dummy);
         } 
      if( nret != 8 ){
	fprintf(stderr,"WARNING: ReadEvolveRefineFile cannot interpret line %s",line);
	continue;
      }
      i++;
      EvolveCoolingRefineRegionNtimes++;
    }

    fclose(fptr);

    // Error check - are there too many input times?
    if(EvolveCoolingRefineRegionNtimes > MAX_REFINE_REGIONS){
      fprintf(stderr, "Too many EvolveCoolingRefineRegion times in your file!\nIncrease MAX_REFINE_REGIONS in macros_and_parameters.h!\n");
      return FAIL;
    }

    /* print out debugging information for evolving MustRefine region */
    if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){

      printf("ReadEvolveRefineFile: I have a CoolingRefineRegion with TimeType %"ISYM" \n",
	     CoolingRefineRegionTimeType);
      
      printf("ReadEvolveRefineFile: And here is what I think my times, edges, and minimum levels are:\n");

      for(int i=0; i<EvolveCoolingRefineRegionNtimes; i++){
	printf("ReadEvolveRefineFile (CoolingRefineRegion): %"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM"\n",
	       EvolveCoolingRefineRegionTime[i],
	       EvolveCoolingRefineRegionLeftEdge[i][0],
	       EvolveCoolingRefineRegionLeftEdge[i][1],
	       EvolveCoolingRefineRegionLeftEdge[i][2],
	       EvolveCoolingRefineRegionRightEdge[i][0],
	       EvolveCoolingRefineRegionRightEdge[i][1],
	       EvolveCoolingRefineRegionRightEdge[i][2]);
      } // for loop

      fflush(stdout);

    } // if(debug1 && MyProcessorNumber == ROOT_PROCESSOR)

  } // if((CoolingRefineRegionTimeType == 0) || (CoolingRefineRegionTimeType == 1))
    
  /* Read in data file for an evolving MultiRefineRegion  */
  if((MultiRefineRegionTimeType == 0) || (MultiRefineRegionTimeType == 1)){

      if ((fptr = fopen(MultiRefineRegionFile, "r")) == NULL) {
        fprintf(stderr, "Error opening MultiRefine region file %s.\n", MultiRefineRegionFile);
        return FAIL;
      }

      char line[MAX_LINE_LENGTH];
      int nret, i=0, j=0;
      NumberOfMultiRefineTracks = 0;
      NumberOfMultiRefineTimeEntries = 0;
      int TrackInd, TimeInd;
      
      /* Read in header and verify that values are reasonable */
      fgets(line, MAX_LINE_LENGTH, fptr);
      nret = sscanf(line, "%"ISYM, &NumberOfMultiRefineTracks);
      if (nret != 1){
        fprintf(stderr, "WARNING: ReadEvolveRefineFile (MultiRefineRegion) cannot interpret the number of tracks in your track file.");
      }
      fgets(line, MAX_LINE_LENGTH, fptr);
      nret = sscanf(line, "%"ISYM, &NumberOfMultiRefineTimeEntries);
      if (nret != 1){
        fprintf(stderr, "WARNING: ReadEvolveRefineFile (MultiRefineRegion) cannot interpret the number of time entries per track in your track file.");
      }
      
      if(NumberOfMultiRefineTracks > MAX_TRACKS){
        fprintf(stderr, "Too many EvolveMultiRefineRegion tracks in your file!\nIncrease MAX_TRACKS in macros_and_parameters.h!\n");
        return FAIL;
      }
      
      if(NumberOfMultiRefineTimeEntries > MAX_TIME_ENTRIES){
        fprintf(stderr, "Too many EvolveMultiRefineRegion times per track in your file!\nIncrease MAX_TIME_ENTRIES in macros_and_parameters.h!\n");
        return FAIL;
      }
      fclose(fptr);
 
      /* Re-open file and read in the rest of the file */
      if ((fptr = fopen(MultiRefineRegionFile, "r")) == NULL) {
        fprintf(stderr, "Error opening MultiRefine region file %s.\n", MultiRefineRegionFile);
        return FAIL;
      }
      
      fgets(line, MAX_LINE_LENGTH, fptr);  // Read in header again
      fgets(line, MAX_LINE_LENGTH, fptr);
      i = 2; // Number of lines in header
      int trackID[NumberOfMultiRefineTracks];
      while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL)){
        TimeInd = (i-2)%NumberOfMultiRefineTimeEntries;
        TrackInd = int(float(i-2)/float(NumberOfMultiRefineTimeEntries));
        nret = sscanf(line, "%"ISYM "%"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM "%"ISYM "%"FSYM,
              &(trackID[TrackInd]),
              &(EvolveMultiRefineRegionTime[TimeInd]),
              &(EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][0]),
              &(EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][1]),
              &(EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][2]),
              &(EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][0]),
              &(EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][1]),
              &(EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][2]),
              &(EvolveMultiRefineRegionMinimumLevel[TrackInd]),
              &(EvolveMultiRefineRegionMaximumLevel[TrackInd]),
              &(EvolveMultiRefineRegionMinimumStarMass[TrackInd][TimeInd]));
        /* Make sure your arrays correspond to the tracks they should */
        if(TrackInd != trackID[TrackInd]){
          fprintf(stderr, "ReadEvolveRefineFile (MultiRefineRegion) says your track IDs do not match up!\n Calculated: %i; Actual: %i\n",TrackInd,trackID[TrackInd]);
          return FAIL;
        }
        /* Make sure your refine regions are within the simulation volume */
        if( (EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][0] < 0) || (EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][0] > 1) ||
            (EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][1] < 0) || (EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][1] > 1) ||
            (EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][2] < 0) || (EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][2] > 1) ||
            (EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][0] < 0) || (EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][0] > 1) ||
            (EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][1] < 0) || (EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][1] > 1) ||
            (EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][2] < 0) || (EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][2] > 1)){
            fprintf(stderr, "ReadEvolveRefineFile (MultiRefineRegion) says the position of the refine region on line %i of your track file is out of bounds\n", i);
          return FAIL;
        }
        if(debug && MyProcessorNumber == ROOT_PROCESSOR){
           fprintf(stderr,"Here is the line (MultiRefineRegion): %s \n",line);
           fprintf(stderr,". . . and here is the minimum value (MultiRefineRegion): %i \n",EvolveMultiRefineRegionMinimumLevel[TrackInd]);
           fprintf(stderr,". . . and here is the maximum value (MultiRefineRegion): %i \n",EvolveMultiRefineRegionMaximumLevel[TrackInd]);
           fprintf(stderr,". . . and here is my initial minimum stellar mass (MultiRefineRegion): %f \n",EvolveMultiRefineRegionMinimumStarMass[TrackInd][0])
        }
        if( nret != 11 ){
          fprintf(stderr,"WARNING: ReadEvolveRefineFile (MultiRefineRegion) cannot interpret line %s",line);
          continue;
        }
        i++;
      }

      fclose(fptr);

      /* Make sure that all time values are greater than 0 and increase */
          for (i=0; i<NumberOfMultiRefineTimeEntries; i++){
        if (EvolveMultiRefineRegionTime[i] < 0){
          fprintf(stderr, "ReadEvolveRefineRegion (MultiRefineRegion) has found a negative time in your track file.\n");
          return FAIL;
        }
      }
      
      if (MultiRefineRegionTimeType == 0){ // If we're using code time
        for(i=1; i<=NumberOfMultiRefineTimeEntries-1; i++){
          if (EvolveMultiRefineRegionTime[i]-EvolveMultiRefineRegionTime[i-1]<0){
            fprintf(stderr, "ReadEvolveRefineRegion (MultiRefineRegion) has found that the times in your track box decrease.\n Set MultiRefineRegionTimeType=1 if using redshift.\n");
            return FAIL;
          }
        }
      }
      if (MultiRefineRegionTimeType == 1){ // If we're using redshift
        for(i=1; i<=NumberOfMultiRefineTimeEntries-1; i++){
          if (EvolveMultiRefineRegionTime[i]-EvolveMultiRefineRegionTime[i-1]>0){
            fprintf(stderr, "ReadEvolveRefineRegion (MultiRefineRegion) has found that the redshifts in your track box increase.\n Set MultiRefineRegionTimeType=0 if using code time.\n");
            return FAIL;
          }
        }
      }
          
      /* Make sure that minimum and maximum refinement levels are reasonable */
      for (i=0; i<NumberOfMultiRefineTracks; i++){
        if( (EvolveMultiRefineRegionMinimumLevel[i] < 0) || (EvolveMultiRefineRegionMinimumLevel[i] > MaximumRefinementLevel) ||
            (EvolveMultiRefineRegionMaximumLevel[i] < 0) || (EvolveMultiRefineRegionMaximumLevel[i] > MaximumRefinementLevel) ||
            (EvolveMultiRefineRegionMaximumLevel[i] < EvolveMultiRefineRegionMinimumLevel[i])){
          fprintf(stderr, "ReadEvolveRefineRegion (MultiRefineRegion) has found unreasonable refinement levels requested in your track file for track %i.\n",i);
          return FAIL;
        }
      }
          
      /* Make sure that minimum stellar mass for multirefine regions make sense */
      for (i=0; i<NumberOfMultiRefineTracks; i++){
        for (j=0; j<NumberOfMultiRefineTimeEntries; j++){
          if(EvolveMultiRefineRegionMinimumStarMass[i][j]<0){
            fprintf(stderr, "ReadEvolveRefineRegion (MultiRefineRegion) has found a negative minimum stellar mass requested in your track file for track %i.\n",i);
            return FAIL;
          }
          if(EvolveMultiRefineRegionMinimumStarMass[i][j]>pow(10,20)){
            fprintf(stderr, "ReadEvolveRefineRegion (MultiRefineRegion) has found an unreasonably high minimum stellar mass requested in your track file for track %i.\n",i);
            return FAIL;
          }
        }
      }

      /* print out debugging information for evolving MustRefine region */
      if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){

        printf("ReadEvolveMultiRefineFile: I have a MultiRefineRegion with TimeType %"ISYM" \n",
           MultiRefineRegionTimeType);
        
        printf("ReadEvolveRefineFile: And here is what I think my times, edges, minimum and maximum levels, and minimum stellar masses are:\n");

        for(int i=2; i<(NumberOfMultiRefineTimeEntries*NumberOfMultiRefineTracks)+2; i++){
          TimeInd = (i-2)%NumberOfMultiRefineTimeEntries;
          TrackInd = int(float(i-2)/float(NumberOfMultiRefineTimeEntries));
          printf("ReadEvolveRefineFile (MustRefineRegion): %"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM " %"ISYM " %"FSYM"\n",
                 EvolveMultiRefineRegionTime[TimeInd],
                 EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][0],
                 EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][1],
                 EvolveMultiRefineRegionLeftEdge[TrackInd][TimeInd][2],
                 EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][0],
                 EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][1],
                 EvolveMultiRefineRegionRightEdge[TrackInd][TimeInd][2],
                 EvolveMultiRefineRegionMinimumLevel[TrackInd],
                 EvolveMultiRefineRegionMaximumLevel[TrackInd],
                 EvolveMultiRefineRegionMinimumStarMass[TrackInd][TimeInd]);
        } // for loop

        fflush(stdout);

      } // if(debug1 && MyProcessorNumber == ROOT_PROCESSOR)

    } // if((MustRefineRegionTimeType == 0) || (MustRefineRegionTimeType == 1))



  return SUCCESS;
}
