/***********************************************************************
/
/
************************************************************************/

#include <cstdio>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <math.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/**************************** Functions Prototypes ******************************/

int ReadPreSNFeedbackTable(char *name)
{

  hid_t  file_id, grp_id, dset_id, dspace_id, attr_id; 
  herr_t status;
  herr_t h5_error = -1;

  long_int *num_met = new long_int[1];
  long_int *num_age = new long_int[1];

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    
    /* Open indexer group whose data will help us navigate the tables */

    if (debug) fprintf(stderr,"Reading from %s.\n",name);
    file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
    grp_id = H5Gopen(file_id, "indexer");
    if (grp_id == h5_error) {
      fprintf(stderr, "Can't open indexer group in %s.\n",name);
    }

    status = H5Gclose(grp_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close indexer group in %s.\n",name);
      return FAIL;
    }

    /* Read indexer arrays (initial metal frac & population age) */
    dset_id = H5Dopen(file_id, "/indexer/initial_metal_fraction");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/initial_metal_fraction in %s.\n", name);
      return FAIL;
    }
    dspace_id = H5Dget_space(dset_id);
    if (dspace_id == h5_error) {
      fprintf(stderr, "Can't get data space for /indexer/initial_metal_fraction in %s.\n", name);
      return FAIL;
    }
    *num_met = H5Sget_simple_extent_npoints(dspace_id);
    if (*num_met == h5_error) {
      fprintf(stderr, "Unable to get size of /indexer/initial_metal_fraction in %s.",name);
    }
    status = H5Sclose(dspace_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close indexer/initial_metal_fraction data space in %s.",name);
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close indexer/initial_metal_fraction data set in %s.",name);
    }

    dset_id = H5Dopen(file_id, "/indexer/population_age");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/population_age in %s.\n", name);
      return FAIL;
    }
    dspace_id = H5Dget_space(dset_id);
    if (dspace_id == h5_error) {
      fprintf(stderr, "Can't get data space for /indexer/population_age in %s.\n", name);
      return FAIL;
    }
    *num_age = H5Sget_simple_extent_npoints(dspace_id);
    if (*num_age == h5_error) {
      fprintf(stderr, "Unable to get size of /indexer/population_age in %s.",name);
    }
    status = H5Sclose(dspace_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close indexer/population_age data space in %s.",name);
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close /indexer/population_age data set in %s.",name);
    }

  } // end root

  /* Store array sizes for later */

#ifdef USE_MPI
  MPI_Bcast(num_met, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(num_age, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  pSNFBTable.n_met = *num_met;
  pSNFBTable.n_age = *num_age;
  if (debug) 
    fprintf(stderr, "Pre-SN Feedback table has %d initial metal fractions & %d ages.\n",
            pSNFBTable.n_met, pSNFBTable.n_age);
  delete [] num_met;
  delete [] num_age;

  /* get and broadcast the rest of the data */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    /* Read initial metal fractions */
    pSNFBTable.ini_met = new double[pSNFBTable.n_met];
    dset_id = H5Dopen(file_id, "/indexer/initial_metal_fraction");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/initial_metal_fraction in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, pSNFBTable.ini_met);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /indexer/initial_metal_fraction in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /indexer/initial_metal_fraction in %s.\n",name);
      return FAIL;
    }

    /* Read population ages */
    pSNFBTable.pop_age = new double[pSNFBTable.n_age];
    dset_id = H5Dopen(file_id, "/indexer/population_age");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/population_age in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, pSNFBTable.pop_age);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /indexer/population_age in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /indexer/population_age in %s.\n",name);
      return FAIL;
    }

    /* Read mass yields */
    pSNFBTable.mass_yield = new double[pSNFBTable.n_met*pSNFBTable.n_age];
    dset_id = H5Dopen(file_id, "/SB99_models/wind_mass_rate");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /SB99_models/wind_mass_rate in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, pSNFBTable.mass_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /SB99_models/wind_mass_rate in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /SB99_models/wind_mass_rate in %s.\n",name);
      return FAIL;
    }

    /* Read metal mass yields */
    pSNFBTable.metm_yield = new double[pSNFBTable.n_met*pSNFBTable.n_age];
    dset_id = H5Dopen(file_id, "/SB99_models/wind_metal_mass_rate");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /SB99_models/wind_metal_mass_rate in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, pSNFBTable.metm_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /SB99_models/wind_metal_mass_rate in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /SB99_models/wind_metal_mass_rate in %s.\n",name);
      return FAIL;
    }

    /* Read wind momentum */
    pSNFBTable.mom_rate = new double[pSNFBTable.n_met*pSNFBTable.n_age];
    dset_id = H5Dopen(file_id, "/SB99_models/wind_and_Lbol_momentum");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /SB99_models/wind_and_Lbol_momentum in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, pSNFBTable.mom_rate);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /SB99_models/wind_and_Lbol_momentum in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /SB99_models/wind_and_Lbol_momentum in %s.\n",name);
      return FAIL;
    }

    /* Close file */
    status = H5Fclose (file_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close file %s",name);
    }

  } else { // not root processor

    pSNFBTable.ini_met = new double[pSNFBTable.n_met];
    pSNFBTable.pop_age = new double[pSNFBTable.n_age];
    pSNFBTable.mass_yield = new double[pSNFBTable.n_met*pSNFBTable.n_age];
    pSNFBTable.metm_yield = new double[pSNFBTable.n_met*pSNFBTable.n_age];
    pSNFBTable.mom_rate = new double[pSNFBTable.n_met*pSNFBTable.n_age];

  } // end not root

  // broadcast
#ifdef USE_MPI
  MPI_Bcast(pSNFBTable.ini_met, pSNFBTable.n_met, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(pSNFBTable.pop_age, pSNFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(pSNFBTable.mass_yield, pSNFBTable.n_met*pSNFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(pSNFBTable.metm_yield, pSNFBTable.n_met*pSNFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(pSNFBTable.mom_rate, pSNFBTable.n_met*pSNFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif
  
  return SUCCESS;
}

