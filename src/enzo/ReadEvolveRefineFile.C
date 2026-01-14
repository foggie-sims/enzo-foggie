/*------------------------------------------------------------------------

  READ EVOLVING REFINE REGION FILE
  By John Wise

  File format: 
    EvolveRefineRegion:        (time or redshift) x_left y_left z_left x_right y_right z_right
    EvolveMustRefineRegion:    (time or redshift) x_left y_left z_left x_right y_right z_right min_level
    EvolveCoolingRefineRegion: (time or redshift) x_left y_left z_left x_right y_right z_right min_level
    EvolveMultiRefineRegion:   track_index (time or redshift) x_left y_left z_left x_right y_right z_right min_level max_level min_star_mass

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
     14 December 2023: ACW -- Added evolving MultiRefine regions

------------------------------------------------------------------------*/

#include <cstdio>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "hdf5.h"
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

void WriteListOfInts(FILE *fptr, int N, int nums[]);

int ReadEvolveRefineFile(void)
{

  hid_t  file_id, grp_id, dset_id, attr_id, trk_id, tme_id;
  herr_t status;
  herr_t h5_error = -1;
  int i=0, j=0, k=0;
  NumberOfMultiRefineTracks = 0;
  int TotalNumberOfFlaggingMethods=0;
  int CMind;
  int ListOfFlaggingMethodsInUse[MAX_FLAGGING_METHODS];

  /* Read in data file for an evolving MultiRefineRegion */
  if((MultiRefineRegionTimeType == 0) || (MultiRefineRegionTimeType == 1)){
    if (MyProcessorNumber == ROOT_PROCESSOR) {

      /* Count up number of global cell flagging methods in use */
      for (i=0; i<MAX_FLAGGING_METHODS; i++){
        if (CellFlaggingMethod[i]!=INT_UNDEFINED){
          ListOfFlaggingMethodsInUse[TotalNumberOfFlaggingMethods] = CellFlaggingMethod[i];
          TotalNumberOfFlaggingMethods++;
        }
      }

      /* Add in any cell flagging methods being used in static MultiRefineRegions and
         make sure we're not going above the maximum for all MRRs+global */
      for (i=0; i<NumberOfStaticMultiRefineRegions; i++){
        for (j=0; j<MAX_FLAGGING_METHODS; j++){
          if (MultiRefineRegionFlaggingMethod[i][j] != INT_UNDEFINED){
            CMind = INT_UNDEFINED;
            for (k=0; k<MAX_FLAGGING_METHODS; k++){
              if (MultiRefineRegionFlaggingMethod[i][j] == CellFlaggingMethod[k]){
                CMind = j;
              }
            }
            if (CMind == INT_UNDEFINED){
              if (TotalNumberOfFlaggingMethods >= MAX_FLAGGING_METHODS){
                fprintf(stderr, 'ReadEvolveRefineFile: Too many cell flagging methods requested. Increase MAX_FLAGGING_METHODS in macros_and_parameters.h.')
              }
              else{
                ListOfFlaggingMethods[TotalNumberOfFlaggingMethods] = MultiRefineRegionFlaggingMethod[i][j];
                TotalNumberOfFlaggingMethods++;
              }
            }
          }
        }
      }

      if (debug){
        fprintf(stderr,"Reading MultiRefineRegion data from %s.\n",MultiRefineRegionFile);
      } 

      /* Read in group containing data for all tracks */
      file_id = H5Fopen(MultiRefineRegionFile, H5F_ACC_RDONLY, H5P_DEFAULT);
      grp_id = H5Gopen(file_id, "AllTracks");
      if (grp_id == h5_error) {
        fprintf(stderr, "Can't open AllTracks group in %s.\n",MultiRefineRegionFile);
      }

      /* Read in number of tracks and verify that it is reasonable */
      attr_id = H5Aopen_name(grp_id, "NTracks");
      if (attr_id == h5_error) {
        fprintf(stderr,"Failed to open NTracks attribute in %s.\n",MultiRefineRegionFile);
        return FAIL;
      }
      status = H5Aread(attr_id, HDF5_I8, NumberOfMultiRefineTracks);
      if (attr_id == h5_error) {
        fprintf(stderr,"Failed to read NTracks in AllTracks group of %s.\n",MultiRefineRegionFile);
        return FAIL;
      }
      status = H5Aclose(attr_id);
      if (attr_id == h5_error) {
        fprintf(stderr,"Failed to close NTracks in AllTracks group of %s.\n",MultiRefineRegionFile);
        return FAIL;
      }

      if(NumberOfMultiRefineTracks > MAX_TRACKS){
        fprintf(stderr, "Too many EvolveMultiRefineRegion tracks in %s!\nIncrease MAX_TRACKS in macros_and_parameters.h!\n",MultiRefineRegionFile);
        return FAIL;
      }
    } // if (MyProcessorNumber == ROOT_PROCESSOR)

    /* Broadcast number of tracks and create struct array of tracks*/
    #ifdef MPI
      MPI_Bcast(NumberOfMultiRefineTracks, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
    #endif
    MRRTracks = new MultiRefineRegionTrackType[NumberOfMultiRefineTracks];

    /* For each track, check whether it is enabled and, if so, read in number of time entries, 
      number of refinement types allowed, and list of refinement types */
    for (i=0; i<NumberOfMultiRefineTracks; i++){
      MRRTracks[i].Enabled = 0;

      if (MyProcessorNumber == ROOT_PROCESSOR) {
        char trkname[11]; // note this assumes MAX_TRACKS<1e4
        sprintf(trkname,"Track_%04d",i);

        trk_id = H5Gopen(grp_id, trkname);
        if (trk_id == h5_error) {
          fprintf(stderr, "Can't open %s group in %s.\n",trkname,MultiRefineRegionFile);
          return FAIL;
        }

        /* Is this track enabled? */
        attr_id = H5Aopen_name(trk_id, "Enabled");
        if (attr_id == h5_error) {
          fprintf(stderr,"Failed to open Enabled attribute in %s for %s.\n",MultiRefineRegionFile,trkname);
          return FAIL;
        }
        status = H5Aread(attr_id, HDF5_I8, MRRTracks[i].Enabled);
        if (attr_id == h5_error) {
          fprintf(stderr,"Failed to read Enabled in %s group of %s.\n",trkname,MultiRefineRegionFile);
          return FAIL;
        }
        status = H5Aclose(attr_id);
        if (attr_id == h5_error) {
          fprintf(stderr,"Failed to close Enabled in %s group of %s.\n",trkname,MultiRefineRegionFile);
          return FAIL;
        } 
      } // if (MyProcessorNumber == ROOT_PROCESSOR)

      #ifdef MPI
        MPI_Bcast(MRRTracks[i].Enabled, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
      #endif
     
      /* If it is enabled, read in everything else */
      if (MRRTracks[i].Enabled==1){

        if (MyProcessorNumber == ROOT_PROCESSOR) {
          /* Read in number of refinement types being used and how many time entries there are */
          attr_id = H5Aopen_name(trk_id, "NRef");
          if (attr_id == h5_error) {
            fprintf(stderr,"Failed to open NRef attribute in %s for %s.\n",MultiRefineRegionFile,trkname);
            return FAIL;
          }
          status = H5Aread(attr_id, HDF5_I8, MRRTracks[i].NRefTypes);
          if (attr_id == h5_error) {
            fprintf(stderr,"Failed to read NRef in %s group of %s.\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }
          status = H5Aclose(attr_id);
          if (attr_id == h5_error) {
            fprintf(stderr,"Failed to close NRef in %s group of %s.\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }

          if(MRRTracks[i].NRefTypes>MAX_FLAGGING_METHODS){
            fprintf(stderr;"Too many flagging methods requested for %s in %s.\nIncrease MAX_FLAGGING_METHODS in macros_and_parameters.h.\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }

          attr_id = H5Aopen_name(trk_id, "NTimes");
          if (attr_id == h5_error) {
            fprintf(stderr,"Failed to open NTimes attribute in %s for %s.\n",MultiRefineRegionFile,trkname);
            return FAIL;
          }
          status = H5Aread(attr_id, HDF5_I8, MRRTracks[i].NTimeEntries);
          if (attr_id == h5_error) {
            fprintf(stderr,"Failed to read NTimes in %s group of %s.\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }
          status = H5Aclose(attr_id);
          if (attr_id == h5_error) {
            fprintf(stderr,"Failed to close NTimes in %s group of %s.\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }

          if(MRRTracks[i].NTimeEntries>MAX_TIME_ENTRIES){
            fprintf(stderr;"Too many time entries requested for %s in %s.\nIncrease MAX_TIME_ENTRIES in macros_and_parameters.h!\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }
        } // if (MyProcessorNumber == ROOT_PROCESSOR)

        /* Broadcast values we need to make new arrays */
        #ifdef MPI
          MPI_Bcast(MRRTracks[i].NRefTypes, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
          MPI_Bcast(MRRTracks[i].NTimeEntries, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
        #endif

        MRRTracks[i].RefTypes = new int[MRRTracks[i].NRefTypes];
        MRRTracks[i].TimeEntries = new MultiRefineRegionTimeEntry[MRRTracks[i].NTimeEntries];

        if (MyProcessorNumber == ROOT_PROCESSOR){

          /* Read in refinement types and make sure a reasonable number are requested */
          dset_id = H5Dopen(trk_id, "RefTypes");
          if (dset_id == h5_error) {
            fprintf(stderr,"Can't open RefTypes in %s group of %s.\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }
          status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                            H5S_ALL, H5P_DEFAULT, MRRTracks[i].RefTypes);
          if (status == h5_error) {
            fprintf(stderr, "Failed to read RefTypes in %s group of %s.\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }
          status = H5Dclose(dset_id);
          if (status == h5_error) {
            fprintf(stderr,"Failed to close RefTypes in %s group of %s.\n",trkname,MultiRefineRegionFile);
            return FAIL;
          }

          /* Check to see if any new refinement methods are being requested for this region. 
             If so, make sure we're not going above the maximum for all MRRs+global */
          for (j=0; j<MAX_FLAGGING_METHODS; j++){
            if (MRRTracks[i].RefTypes[j] != INT_UNDEFINED){
              CMind = INT_UNDEFINED;
              for (k=0; k<MAX_FLAGGING_METHODS; k++){
                if (MRRTracks[i].RefTypes[j] == CellFlaggingMethod[k]){
                  CMind = j;
                }
              }
              if (CMind == INT_UNDEFINED){
                if (TotalNumberOfFlaggingMethods >= MAX_FLAGGING_METHODS){
                  fprintf(stderr, 'ReadEvolveRefineFile: Too many cell flagging methods requested; %'ISYM' in addition to: ',MRRTracks[i].RefTypes[j]);
                  WriteListOfInts(stderr,MAX_FLAGGING_METHODS,ListOfFlaggingMethods);
                  fprintf(stderr, 'Please increase MAX_FLAGGING_METHODS in macros_and_parameters.h.\n');
                }
                else{
                  ListOfFlaggingMethods[TotalNumberOfFlaggingMethods] = MRRTracks[i].RefTypes[j];
                  TotalNumberOfFlaggingMethods++;
                }
              }
            }
          }
      


        } // if (MyProcessorNumber == ROOT_PROCESSOR)

        /* Broadcast refinement type array */
        #ifdef MPI
          MPI_Bcast(MRRTracks[i].RefTypes, MRRTracks[i].NRefTypes, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
        #endif

        /* For each time entry, read in position of track, minimum star particle mass,
           and minimum and maximum levels for each type of refinement */
        for(j=0; j<MRRTracks[i].NTimeEntries; j++){

          /* create arrays to hold min and max levels for each refinement type */
          MRRTracks[i].TimeEntries[j].MaxLevels = new int[MRRTracks[i].NRefTypes];
          MRRTracks[i].TimeEntries[j].MinLevels = new int[MRRTracks[i].NRefTypes];

          if (MyProcessorNumber == ROOT_PROCESSOR){

            char tmename[10]; // note this assumes MAX_TIME_ENTRIES<1e4
            sprintf(trkname,"Time_%04d",j);

            tme_id = H5Gopen(trk_id, tmename);
            if (tme_id == h5_error) {
              fprintf(stderr, "Can't open %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            /* Read in time associated with this entry and check that it's sane */
            attr_id = H5Aopen_name(tme_id, "TimeValue");
            if (attr_id == h5_error) {
              fprintf(stderr,"Failed to open TimeValue attribute in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Aread(attr_id, HDF5_R8, MRRTracks[i].TimeEntries[j].Time);
            if (attr_id == h5_error) {
              fprintf(stderr,"Failed to read TimeValue in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Aclose(attr_id);
            if (attr_id == h5_error) {
              fprintf(stderr,"Failed to close TimeValue in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            if (MRRTracks[i].TimeEntries[j].Time < 0.0){
              fprintf(stderr, "ReadEvolveRefineRegion has found a negative time in %s in %s.\n",trkname,MultiRefineRegionFile);
              return FAIL;
            }

            if (j>0){
              if (MultiRefineRegionTimeType == 0){ // if we're using code time...
                if (MRRTracks[i].TimeEntries[j].Time-MRRTracks[i].TimeEntries[j-1].Time < 0.0){
                  fprintf(stderr, "ReadEvolveRefineRegion has found that the times in %s of %s decrease.\n Set MultiRefineRegionTimeType=1 if using redshift.\n",trkname,MultiRefineRegionFile);
                  return FAIL;
                }
              } 
              if (MultiRefineRegionTimeType == 1){ // If we're using redshift
                  if (MRRTracks[i].TimeEntries[j].Time-MRRTracks[i].TimeEntries[j-1].Time > 0.0){
                  fprintf(stderr, "ReadEvolveRefineRegion has found that the redshifts in %s of %s increase.\n Set MultiRefineRegionTimeType=0 if using code time.\n",trkname,MultiRefineRegionFile);
                  return FAIL;
                }             
              }
            } // if (j>0)

            /* Read in minimum star particle mass at this entry and check that it's sane */
            attr_id = H5Aopen_name(tme_id, "MinStarMass");
            if (attr_id == h5_error) {
              fprintf(stderr,"Failed to open MinStarMass attribute in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Aread(attr_id, HDF5_R8, MRRTracks[i].TimeEntries[j].MinStarMass);
            if (attr_id == h5_error) {
              fprintf(stderr,"Failed to read MinStarMass in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Aclose(attr_id);
            if (attr_id == h5_error) {
              fprintf(stderr,"Failed to close MinStarMass in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            if(MRRTracks[i].TimeEntries[j].MinStarMass>1.0e+20){
              fprintf(stderr;"Unreasonably high minimum star particle mass requested for %s in %s in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            if(MRRTracks[i].TimeEntries[j].MinStarMass<0.0){
              fprintf(stderr;"Negative minimum star particle mass requested for %s in %s in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            /* Read in position at this entry and check that it's sane */
            dset_id = H5Dopen(tme_id, "Position");
            if (dset_id == h5_error) {
              fprintf(stderr,"Can't open Position in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                              H5S_ALL, H5P_DEFAULT, MRRTracks[i].TimeEntries[j].Pos);
            if (status == h5_error) {
              fprintf(stderr, "Failed to read Position in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Dclose(dset_id);
            if (status == h5_error) {
              fprintf(stderr,"Failed to close Position in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            if( (MRRTracks[i].TimeEntries[j].Pos[0] < 0.0) || (MRRTracks[i].TimeEntries[j].Pos[0] > 1.0) ||
                (MRRTracks[i].TimeEntries[j].Pos[1] < 0.0) || (MRRTracks[i].TimeEntries[j].Pos[1] > 1.0) ||
                (MRRTracks[i].TimeEntries[j].Pos[2] < 0.0) || (MRRTracks[i].TimeEntries[j].Pos[2] > 1.0) ||
                (MRRTracks[i].TimeEntries[j].Pos[3] < 0.0) || (MRRTracks[i].TimeEntries[j].Pos[3] > 1.0) ||
                (MRRTracks[i].TimeEntries[j].Pos[4] < 0.0) || (MRRTracks[i].TimeEntries[j].Pos[4] > 1.0) ||
                (MRRTracks[i].TimeEntries[j].Pos[5] < 0.0) || (MRRTracks[i].TimeEntries[j].Pos[5] > 1.0)){
              fprintf(stderr, "ReadEvolveRefineFile has found that the position of %s in %s in %s is out of bounds\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            /* Read in minimum and maximum refinement levels for each refinement type at this entry */

            dset_id = H5Dopen(tme_id, "MaxLevels");
            if (dset_id == h5_error) {
              fprintf(stderr,"Can't open MaxLevels in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Dread(dset_id, HDF5_I8, H5S_ALL, 
                              H5S_ALL, H5P_DEFAULT, MRRTracks[i].TimeEntries[j].MaxLevels);
            if (status == h5_error) {
              fprintf(stderr, "Failed to read MaxLevels in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Dclose(dset_id);
            if (status == h5_error) {
              fprintf(stderr,"Failed to close MaxLevels in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            dset_id = H5Dopen(tme_id, "MinLevels");
            if (dset_id == h5_error) {
              fprintf(stderr,"Can't open MinLevels in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Dread(dset_id, HDF5_I8, H5S_ALL, 
                              H5S_ALL, H5P_DEFAULT, MRRTracks[i].TimeEntries[j].MinLevels);
            if (status == h5_error) {
              fprintf(stderr, "Failed to read MinLevels in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }
            status = H5Dclose(dset_id);
            if (status == h5_error) {
              fprintf(stderr,"Failed to close MinLevels in %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

            for (k=0; k<MRRTracks[i].NRefTypes; k++){
              if( (MRRTracks[i].TimeEntries[j].MinLevels[k] < 0) || (MRRTracks[i].TimeEntries[j].MinLevels[k] > MaximumRefinementLevel) ||
                  (MRRTracks[i].TimeEntries[j].MaxLevels[k] < 0) || (MRRTracks[i].TimeEntries[j].MinLevels[k] > MaximumRefinementLevel) ||
                  (MRRTracks[i].TimeEntries[j].MaxLevels[k] < MRRTracks[i].TimeEntries[j].MinLevels[k])){
                fprintf(stderr, "ReadEvolveRefineRegion has found unreasonable refinement levels requested for %s in %s in %s.\n",tmename,trkname,MultiRefineRegionFile);
                return FAIL;
              }
            }

            status = H5Gclose(tme_id);
            if (status == h5_error) {
              fprintf(stderr, "Failed to close %s group in %s group in %s.\n",tmename,trkname,MultiRefineRegionFile);
              return FAIL;
            }

          } // if (MyProcessorNumber == ROOT_PROCESSOR)

          /* Broadcast all TimeEntry variables */
          #ifdef MPI
            MPI_Bcast(MRRTracks[i].TimeEntries[j].Pos, 6, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
            MPI_Bcast(MRRTracks[i].TimeEntries[j].MinLevels, MRRTracks[i].NRefTypes, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
            MPI_Bcast(MRRTracks[i].TimeEntries[j].MaxLevels, MRRTracks[i].NRefTypes, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
            MPI_Bcast(MRRTracks[i].TimeEntries[j].Time, 1, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
            MPI_Bcast(MRRTracks[i].TimeEntries[j].MinStarMass, 1, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
          #endif
        } // for(j=0; j<MRRTracks[i].NTimeEntries; j++)


      } // if (MRRTracks[i].Enabled==1)
      /* If the track is not enabled, we don't need to do anything else*/
      else{ 
        if (debug){
          fprintf(stderr,"%s is not enabled.\n",trkname);
        }
      } 

      /* Close the open groups and the file */
      if (MyProcessorNumber == ROOT_PROCESSOR){
        status = H5Gclose(trk_id);
        if (status == h5_error) {
          fprintf(stderr, "Failed to close %s group in %s.\n",trkname,MultiRefineRegionFile);
          return FAIL;
        }
      } // if (MyProcessorNumber == ROOT_PROCESSOR)
    } // for (i=0; i<NumberOfMultiRefineTracks; i++)

    if (MyProcessorNumber == ROOT_PROCESSOR){
      status = H5Gclose(grp_id);
      if (status == h5_error) {
        fprintf(stderr, "Failed to close AllTracks group in %s.\n",MultiRefineRegionFile);
        return FAIL;
      }

      status = H5Fclose(file_id);
      if (status == h5_error) {
        fprintf(stderr, "Failed to close %s.\n",MultiRefineRegionFile);
        return FAIL;
      }
    } // if (MyProcessorNumber == ROOT_PROCESSOR)
  } // if((MultiRefineRegionTimeType == 0) || (MultiRefineRegionTimeType == 1))
  return SUCCESS;
}