/***********************************************************************
/
/  GRID CLASS (COMPUTE FLOORS AND CEILINGS FOR SELECT BARYON FIELDS)
/
/  written by: Brian O'Shea
/  date:       May 2025
/  modified1:
/
/  PURPOSE:  Apply bounds (floors and ceilings) to select baryon fields.
/
/  RETURNS:  Modified BaryonField quantities.  
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"

/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int grid::ApplyBoundsToBaryonFields()
{

  // I need to fix this!
  if (HydroMethod == Zeus_Hydro){
    ENZO_FAIL("Grid::ApplyBoundsToBaryonFields: does not support ZEUS hydro method (HydroMethod=2) yet!\n");  
  }

  // I need to fix this!
  if (GridRank != 3){
    ENZO_FAIL("Grid::ApplyBoundsToBaryonFields: does not support non-3D simulations yet!  (GridRank != 3)\n");  
  }

  /* declarations */
  float orig_ge_from_te, new_ge_from_te;
  int i, size = 1, velocity_changed = 0;

  /* more declarations.  TODO: THESE ARE MEANT TO BE REPLACED BY MORE
     SENSIBLE BOUNDS (AND PROBABLY INPUT PARAMETERS) AND ARE JUST
     HERE AS AN EXAMPLE.  */
  float rho_floor = 1.0e-5, rho_ceiling = 1.0e+15,
    vel_min = -1.0e+5, vel_max = 1.0e+5,  
    ge_floor = 1.0e-20, ge_ceiling = 1.0e+10;

  /* Compute the size of the baryon fields given grid dimensions. */
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, internal and total energy, velocity, magnetic fields*/
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
    Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {     
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    } 
 
// Loop over the baryon fields, including ghost zones
for (i = 0; i < size; i++) {

  // In loop, check for NaN/inf and CRASH if you have those
  // (we could in principle replace these with a min or max value or the value of a nearby cell, 
  // but a NaN/inf is usually a sign of a bigger problem)
  if (isfinite(BaryonField[DensNum][i]) == FALSE || isfinite(BaryonField[TENum][i]) == FALSE){
    ENZO_FAIL("Error in grid->ApplyBoundsToBaryonFields. Density or Total Energy is NaN or inf.\n");
  }

  if(DualEnergyFormalism)
    if(isfinite(BaryonField[GENum][i]) == FALSE){
      ENZO_FAIL("Error in grid->ApplyBoundsToBaryonFields. Internal energy is NaN or inf.\n");
    }

  // TODO: if we apply density floors/ceilings we're going to have to rescale multispecies fields too!

  // Check for density floor
  if(BaryonField[DensNum][i] < rho_floor){
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Density %"FSYM" replaced with lower bound %"FSYM, BaryonField[DensNum][i], rho_floor);
    }
    BaryonField[DensNum][i] = rho_floor;
  }

  // Check for density ceiling
  if(BaryonField[DensNum][i] > rho_ceiling){
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Density %"FSYM" replaced with upper bound %"FSYM, BaryonField[DensNum][i], rho_ceiling);
    }
    BaryonField[DensNum][i] = rho_ceiling;
  }

  // set these to physically unreasonable values so we can check against it later.
  orig_ge_from_te = new_ge_from_te = FLOAT_UNDEFINED;
  // this will be our cue that velocity is modified
  velocity_changed = 0;

  // calculate gas energy from total energy (both are specific energy numbers)
  // TODO: This is just for hydro now, will have to update for MHD
  orig_ge_from_te = BaryonField[TENum][i] - 0.5*( POW(BaryonField[Vel1Num][i],2.0) + POW(BaryonField[Vel2Num][i],2.0) + POW(BaryonField[Vel3Num][i],2.0) );

  // TODO -- MAKE OUTPUT MORE VERBOSE
  // check for internal energy floor
  if(orig_ge_from_te < ge_floor){
    new_ge_from_te = ge_floor;
  }

  // check for internal energy ceiling
  if(orig_ge_from_te > ge_ceiling){
    new_ge_from_te = ge_ceiling;
  }

  // TODO -- MAKE OUTPUT MORE VERBOSE
  // x-velocity - make sure it's within bounds.
  if(BaryonField[Vel1Num][i] < vel_min){
    BaryonField[Vel1Num][i] = vel_min;
    velocity_changed++;
  }
  if(BaryonField[Vel1Num][i] > vel_max){
    BaryonField[Vel1Num][i] = vel_max;
    velocity_changed++;
  }

  // y-velocity - make sure it's within bounds.
  if(BaryonField[Vel2Num][i] < vel_min){
    BaryonField[Vel2Num][i] = vel_min;
    velocity_changed++;
  }
  if(BaryonField[Vel2Num][i] > vel_max){
    BaryonField[Vel2Num][i] = vel_max;
    velocity_changed++;
  }

  // z-velocity - make sure it's within bounds.
  if(BaryonField[Vel3Num][i] < vel_min){
    BaryonField[Vel1Num][i] = vel_min;
    velocity_changed++;
  }
  if(BaryonField[Vel3Num][i] > vel_max){
    BaryonField[Vel1Num][i] = vel_max;
    velocity_changed++;
  }

  // if internal specific energy or a velocity component is changed, we need to update total energy field
  if(new_ge_from_te > FLOAT_UNDEFINED || velocity_changed > 0){
    BaryonField[TENum][i] = new_ge_from_te + 0.5*( POW(BaryonField[Vel1Num][i],2.0) + POW(BaryonField[Vel2Num][i],2.0) + POW(BaryonField[Vel3Num][i],2.0) );
  }

  // if we have a gas energy field and it has been modified, we need to reset it.
  if(DualEnergyFormalism && new_ge_from_te > FLOAT_UNDEFINED){
    BaryonField[TENum][i] = new_ge_from_te;
  }

} // for (i = 0; i < size; i++)

  return SUCCESS;
}
