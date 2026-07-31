/***********************************************************************
/
/  Read Chemical Feedback Table
/
/  Populates ChemFeedbackTableType struct with element-specific yields
/
/  Written following the pattern of ReadFeedbackTable.C
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

int ReadChemFeedbackTable(char *name)
{

  hid_t  file_id, grp_id, dset_id, dspace_id, attr_id;
  herr_t status;
  herr_t h5_error = -1;

  long_int C_index;
  long_int O_index;
  long_int Mg_index;
  long_int Si_index;
  long_int Fe_index;

  long_int num_met;
  long_int num_age;

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    
    /* Open indexer group whose data will help us navigate the tables */

    if (debug) fprintf(stderr,"Reading chemical feedback from %s.\n",name);
    file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
    grp_id = H5Gopen(file_id, "indexer");
    if (grp_id == h5_error) {
      fprintf(stderr, "Can't open data group in %s.\n",name);
    }

    /* Get the indices of each chemical element */
    attr_id = H5Aopen_name(grp_id, "C_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open C_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, &C_index);
    if (status == h5_error) {
      fprintf(stderr,"Failed to read C_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close C_index in indexer group of %s.\n",name);
      return FAIL;
    }

    attr_id = H5Aopen_name(grp_id, "O_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open O_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, &O_index);
    if (status == h5_error) {
      fprintf(stderr,"Failed to read O_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close O_index in indexer group of %s.\n",name);
      return FAIL;
    }

    attr_id = H5Aopen_name(grp_id, "Mg_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open Mg_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, &Mg_index);
    if (status == h5_error) {
      fprintf(stderr,"Failed to read Mg_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close Mg_index in indexer group of %s.\n",name);
      return FAIL;
    }

    attr_id = H5Aopen_name(grp_id, "Si_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open Si_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, &Si_index);
    if (status == h5_error) {
      fprintf(stderr,"Failed to read Si_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close Si_index in indexer group of %s.\n",name);
      return FAIL;
    }

    attr_id = H5Aopen_name(grp_id, "Fe_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open Fe_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, &Fe_index);
    if (status == h5_error) {
      fprintf(stderr,"Failed to read Fe_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close Fe_index in indexer group of %s.\n",name);
      return FAIL;
    }

    /* finished reading element indexes */

    status = H5Gclose(grp_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close indexer group in %s.\n",name);
      return FAIL;
    }

    /* Check index labels against internal enum (see typedefs.h) */
    if ((C_index != TabC) || (O_index != TabO) ||
        (Mg_index != TabMg) || (Si_index != TabSi) ||
        (Fe_index != TabFe)){
          fprintf(stderr, "Element indexes in %s don't follow expected order.\n",name);
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
    num_met = H5Sget_simple_extent_npoints(dspace_id);
    if (num_met == h5_error) {
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
    num_age = H5Sget_simple_extent_npoints(dspace_id);
    if (num_age == h5_error) {
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

    status = H5Gclose(grp_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close indexer group in %s.\n",name);
      return FAIL;
    }

  } // end root

  /* Store array sizes for later */

#ifdef USE_MPI
  MPI_Bcast(&num_met, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(&num_age, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  ChemFBTable.n_met = num_met;
  ChemFBTable.n_age = num_age;
  if (debug) 
    fprintf(stderr, "Chemical feedback table has %d initial metal fractions & %d ages.\n",
            ChemFBTable.n_met, ChemFBTable.n_age);
  delete [] num_met;
  delete [] num_age;

  /* get and broadcast the rest of the data */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    /* Read initial metal fractions */
    ChemFBTable.ini_met = new double[ChemFBTable.n_met];
    dset_id = H5Dopen(file_id, "/indexer/initial_metal_fraction");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/initial_metal_fraction in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, ChemFBTable.ini_met);
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
    ChemFBTable.pop_age = new double[ChemFBTable.n_age];
    dset_id = H5Dopen(file_id, "/indexer/population_age");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/population_age in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, ChemFBTable.pop_age);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /indexer/population_age in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /indexer/population_age in %s.\n",name);
      return FAIL;
    }

    /* Read Carbon yields */
    ChemFBTable.C_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    dset_id = H5Dopen(file_id, "/sygma_models/C_yield");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /sygma_models/C_yield in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, ChemFBTable.C_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /sygma_models/C_yield in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /sygma_models/C_yield in %s.\n",name);
      return FAIL;
    }

    /* Read Oxygen yields */
    ChemFBTable.O_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    dset_id = H5Dopen(file_id, "/sygma_models/O_yield");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /sygma_models/O_yield in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, ChemFBTable.O_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /sygma_models/O_yield in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /sygma_models/O_yield in %s.\n",name);
      return FAIL;
    }

    /* Read Magnesium yields */
    ChemFBTable.Mg_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    dset_id = H5Dopen(file_id, "/sygma_models/Mg_yield");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /sygma_models/Mg_yield in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, ChemFBTable.Mg_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /sygma_models/Mg_yield in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /sygma_models/Mg_yield in %s.\n",name);
      return FAIL;
    }

    /* Read Silicon yields */
    ChemFBTable.Si_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    dset_id = H5Dopen(file_id, "/sygma_models/Si_yield");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /sygma_models/Si_yield in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, ChemFBTable.Si_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /sygma_models/Si_yield in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /sygma_models/Si_yield in %s.\n",name);
      return FAIL;
    }

    /* Read Iron yields */
    ChemFBTable.Fe_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    dset_id = H5Dopen(file_id, "/sygma_models/Fe_yield");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /sygma_models/Fe_yield in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, ChemFBTable.Fe_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /sygma_models/Fe_yield in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /sygma_models/Fe_yield in %s.\n",name);
      return FAIL;
    }

    /* Close file */
    status = H5Fclose (file_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close file %s",name);
    }

  } else { // not root processor

    ChemFBTable.ini_met = new double[ChemFBTable.n_met];
    ChemFBTable.pop_age = new double[ChemFBTable.n_age];
    ChemFBTable.C_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    ChemFBTable.O_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    ChemFBTable.Mg_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    ChemFBTable.Si_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];
    ChemFBTable.Fe_yield = new double[ChemFBTable.n_met*ChemFBTable.n_age];

  } // end not root

  // broadcast
#ifdef USE_MPI
  MPI_Bcast(ChemFBTable.ini_met, ChemFBTable.n_met, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(ChemFBTable.pop_age, ChemFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(ChemFBTable.C_yield, ChemFBTable.n_met*ChemFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(ChemFBTable.O_yield, ChemFBTable.n_met*ChemFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(ChemFBTable.Mg_yield, ChemFBTable.n_met*ChemFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(ChemFBTable.Si_yield, ChemFBTable.n_met*ChemFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(ChemFBTable.Fe_yield, ChemFBTable.n_met*ChemFBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif
  
  return SUCCESS;
}
