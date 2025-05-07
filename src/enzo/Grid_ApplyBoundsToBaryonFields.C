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
  float *VelocityUnits, float *MassUnits, FLOAT Time);
 
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

  /* Get Units. */
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1;
  double MassUnits = 1;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* declarations */
  float orig_ge_from_te, new_ge_from_te, vel_mag, vel_ratio;
  int i, size = 1;

  /* more declarations.  TODO: THESE ARE MEANT TO BE REPLACED BY INPUT PARAMETERS.
     Density cap: 1e15 g/cm^3
     Velocity cap: 3000 km/s
     GE cap corresponds to a temperature of 10^9 K, assuming mean molecular weight is 0.5
     No floors set  */
  float rho_floor = tiny_number, rho_ceiling = 1.0e+15/DensityUnits,
    vel_min = tiny_number, vel_max = 3.0e+8/(VelocityUnits),  
    ge_floor = tiny_number, ge_ceiling = 2.476e17/POW(VelocityUnits,2.0);

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
  if (std::isfinite(BaryonField[DensNum][i]) == FALSE || std::isfinite(BaryonField[TENum][i]) == FALSE){
    ENZO_FAIL("Error in grid->ApplyBoundsToBaryonFields. Density or Total Energy is NaN or inf.\n");
  }

  if(DualEnergyFormalism)
    if(std::isfinite(BaryonField[GENum][i]) == FALSE){
      ENZO_FAIL("Error in grid->ApplyBoundsToBaryonFields. Internal energy is NaN or inf.\n");
    }

  // TODO: if we apply density floors/ceilings we're going to have to rescale multispecies fields too!

  // Check for density floor
  if(RestrictDensity && BaryonField[DensNum][i] < rho_floor){
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Density %"ESYM" g/cm^3 replaced with lower bound %"ESYM" g/cm^3\n", BaryonField[DensNum][i]*DensityUnits, rho_floor*DensityUnits);
    }
    BaryonField[DensNum][i] = rho_floor;
  }

  // Check for density ceiling
  if(RestrictDensity && BaryonField[DensNum][i] > rho_ceiling){
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Density %"ESYM" g/cm^3 replaced with upper bound %"ESYM" g/cm^3\n", BaryonField[DensNum][i]*DensityUnits, rho_ceiling*DensityUnits);
    }
    BaryonField[DensNum][i] = rho_ceiling;
  }

  // set these to physically unreasonable values so we can check against it later.
  orig_ge_from_te = new_ge_from_te = FLOAT_UNDEFINED;
  // this will be our cue that velocity is modified
  vel_ratio = FLOAT_UNDEFINED;

  // calculate gas energy from total energy (both are specific energy numbers)
  // TODO: This is just for hydro now, will have to update for MHD
  orig_ge_from_te = BaryonField[TENum][i] - 0.5*( POW(BaryonField[Vel1Num][i],2.0) + POW(BaryonField[Vel2Num][i],2.0) + POW(BaryonField[Vel3Num][i],2.0) );

  // check for internal energy floor - manually turned this off for now
/*
  if(orig_ge_from_te < ge_floor){
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Internal Energy %"FSYM" replaced with lower bound %"FSYM, orig_ge_from_te, ge_floor);
    }
    new_ge_from_te = ge_floor;
  }
*/
  // check for internal energy ceiling
  if(RestrictTemperature && orig_ge_from_te > ge_ceiling){
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Internal Energy %"ESYM" erg/g replaced with upper bound %"ESYM" erg/g\n", orig_ge_from_te*POW(VelocityUnits, 2.0), ge_ceiling*POW(VelocityUnits, 2.0));
      printf("Grid::ApplyBoundsToBaryonFields: Corresponding temperature %"ESYM" K replaced with upper bound %"ESYM" K\n", orig_ge_from_te*POW(VelocityUnits, 2.0)*4.04e-9, ge_ceiling*POW(VelocityUnits, 2.0)*4.04e-9);
    }
    new_ge_from_te = ge_ceiling;
  }

  // velocity magnitude - make sure it's within bounds.
  vel_mag = POW( POW(BaryonField[Vel1Num][i],2.0) + POW(BaryonField[Vel2Num][i],2.0) + POW(BaryonField[Vel3Num][i],2.0), 0.5);
/* Manually turning off velocity floor for now
  if(vel_mag < vel_min){
    vel_ratio = vel_min/vel_mag;
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Velocity magnitude %"FSYM" replaced with lower bound %"FSYM, vel_mag, vel_min);
    }
    BaryonField[Vel1Num][i] = BaryonField[Vel1Num][i]*vel_ratio;
    BaryonField[Vel2Num][i] = BaryonField[Vel2Num][i]*vel_ratio;
    BaryonField[Vel3Num][i] = BaryonField[Vel3Num][i]*vel_ratio;
  }
*/
  if(RestrictVelocity && vel_mag > vel_max){
    vel_ratio = vel_max/vel_mag;
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Velocity magnitude %"FSYM" km/s replaced with upper bound %"FSYM" km/s\n", vel_mag*VelocityUnits/1.0e+5, vel_max*VelocityUnits/1.0e+5);
    }
    BaryonField[Vel1Num][i] = BaryonField[Vel1Num][i]*vel_ratio;
    BaryonField[Vel2Num][i] = BaryonField[Vel2Num][i]*vel_ratio;
    BaryonField[Vel3Num][i] = BaryonField[Vel3Num][i]*vel_ratio;
  }

  // if internal specific energy or a velocity component is changed, we need to update total energy field
  if(new_ge_from_te > FLOAT_UNDEFINED || vel_ratio > FLOAT_UNDEFINED){
    if(new_ge_from_te > FLOAT_UNDEFINED){
      BaryonField[TENum][i] = new_ge_from_te + 0.5*( POW(BaryonField[Vel1Num][i],2.0) + POW(BaryonField[Vel2Num][i],2.0) + POW(BaryonField[Vel3Num][i],2.0) );
    }
    // If velocity changed but not internal specific energy, use original internal specific energy
    else{
      BaryonField[TENum][i] = orig_ge_from_te + 0.5*( POW(BaryonField[Vel1Num][i],2.0) + POW(BaryonField[Vel2Num][i],2.0) + POW(BaryonField[Vel3Num][i],2.0) );
    }
  }

  // if we have a gas energy field and it has been modified, we need to reset it.
  if(DualEnergyFormalism && new_ge_from_te > FLOAT_UNDEFINED){
    BaryonField[GENum][i] = new_ge_from_te;
  }

} // for (i = 0; i < size; i++)

  return SUCCESS;
}
