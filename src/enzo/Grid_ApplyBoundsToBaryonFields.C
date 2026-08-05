/***********************************************************************
/
/  GRID CLASS (COMPUTE FLOORS AND CEILINGS FOR SELECT BARYON FIELDS)
/
/  written by: Brian O'Shea and Cassi Lochhaas
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
  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Currently this method only works with PPM hydro. If we extend
     it to Zeus we have to take into account the fact that "total energy" is
     really gas energy in the Zeus algorithm, and if we do the MHD routines
     we also need to worry about magnetic fields in total and gas energy.
     This is tedious but straightforward, and deferred to the future.
     TODO: UPDATE THIS ROUTINE FOR NON-PPM FLUID SOLVERS. */
  if (HydroMethod != PPM_DirectEuler){
    ENZO_FAIL("Grid::ApplyBoundsToBaryonFields: currently only supports PPM Direct Euler (HydroMethod = 0) at the moment!\n");
  }

  /* Get Units. */
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0,
    VelocityUnits = 1.0, TimeUnits = 1.0;
  double MassUnits = 1.0;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* declarations */
  float orig_ge_from_te, new_ge_from_te, vel_mag, vel_ratio;
  int i, size = 1;

  /* more declarations. TODO: THESE SHOULD BE REPLACED BY INPUT PARAMETERS.
     Density cap: 1e-8 g/cm^3 (10^16 particles per cc)
     Density floor: tiny_number
     Velocity magnitude cap: 3000 km/s
     GE cap corresponds to a temperature of 10^9 K, assuming mean molecular weight is 0.6 (fully ionized)
     GE floor corresponds to a temperature of 1 K, assuming mean molecular weight is 1.2 (fully neutral)
    */
  float rho_floor = tiny_number, rho_ceiling = 1.0e-8/DensityUnits,
    vel_max = 3.0e+8/(VelocityUnits),
    ge_floor = 1.032e8/POW(VelocityUnits,2.0), ge_ceiling = 2.476e17/POW(VelocityUnits,2.0);

  /* Compute the size of the baryon fields given grid dimensions. */
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, internal and total energy, velocity, magnetic fields*/
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
    Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {     
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }

  /* Find Multi-species fields if we need them. */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }
 
  /* Strides used by the non-finite-cell report below. */
  int stride_x = 1;
  int stride_y = GridDimension[0];
  int stride_z = GridDimension[0]*GridDimension[1];

// Loop over the baryon fields, including ghost zones
for (i = 0; i < size; i++) {

  /* In loop, check for NaN/inf and CRASH if you have those
     (we could in principle replace these with a min or max value or the value of a nearby cell,
     but a NaN/inf is usually a sign of a bigger problem).  Note that this check may not work for
     very old versions of the C++ STL; if you run into errors just comment it out and recompile.

     TEMPORARY INSTRUMENTATION (velocity NaN hunt): the velocity components are now checked
     as well.  Previously a NaN velocity passed through here silently -- the clamp below tests
     "vel_mag > vel_max", which is false for NaN -- and was only caught a subcycle later, once
     the hydro solver had laundered it into the total energy.  Checking the velocities catches
     the corruption in the subcycle it is created, and the block below dumps the exact cell.
     REMOVE this instrumentation once the source of the NaN is identified. */

  const char *bad_field = NULL;
  if      (std::isfinite(BaryonField[DensNum][i]) == FALSE) bad_field = "Density";
  else if (std::isfinite(BaryonField[TENum][i])   == FALSE) bad_field = "TotalEnergy";
  else if (DualEnergyFormalism &&
           std::isfinite(BaryonField[GENum][i])   == FALSE) bad_field = "GasEnergy";
  else if (std::isfinite(BaryonField[Vel1Num][i]) == FALSE) bad_field = "Velocity1";
  else if (GridRank > 1 &&
           std::isfinite(BaryonField[Vel2Num][i]) == FALSE) bad_field = "Velocity2";
  else if (GridRank > 2 &&
           std::isfinite(BaryonField[Vel3Num][i]) == FALSE) bad_field = "Velocity3";

  if (bad_field != NULL) {

    int ci =  i % GridDimension[0];
    int cj = (i / GridDimension[0]) % GridDimension[1];
    int ck =  i / stride_z;

    int in_ghost = (ci < GridStartIndex[0] || ci > GridEndIndex[0] ||
                    (GridRank > 1 && (cj < GridStartIndex[1] || cj > GridEndIndex[1])) ||
                    (GridRank > 2 && (ck < GridStartIndex[2] || ck > GridEndIndex[2]))) ? 1 : 0;

    double pos[3] = {0.0, 0.0, 0.0};
    pos[0] = (double)(CellLeftEdge[0][ci] + 0.5*CellWidth[0][ci]);
    if (GridRank > 1) pos[1] = (double)(CellLeftEdge[1][cj] + 0.5*CellWidth[1][cj]);
    if (GridRank > 2) pos[2] = (double)(CellLeftEdge[2][ck] + 0.5*CellWidth[2][ck]);

    /* NB: all integers are cast to long long and printed with %lld -- under
       LARGE_INTS macros_and_parameters.h does "#define int long_int", so a
       plain %d here would be wrong.  GridLevel is deliberately NOT reported:
       it is only assigned along the subgrid-marker/radiation paths, so it is
       not reliable for a general grid.  CellWidth identifies the level. */
    fprintf(stderr, "\n===== NaNReport BEGIN (proc %lld) =====\n", (long long)MyProcessorNumber);
    fprintf(stderr, "NaNReport: first non-finite field = %s\n", bad_field);
    fprintf(stderr, "NaNReport: Time %.17g, dtFixed %.17g\n",
            (double)Time, (double)dtFixed);
    fprintf(stderr, "NaNReport: flat index %lld -> (i,j,k) = (%lld,%lld,%lld), %s zone\n",
            (long long)i, (long long)ci, (long long)cj, (long long)ck,
            in_ghost ? "GHOST" : "active");
    fprintf(stderr, "NaNReport: cell center = (%.17g, %.17g, %.17g)\n",
            pos[0], pos[1], pos[2]);
    fprintf(stderr, "NaNReport: GridDimension = (%lld,%lld,%lld), Start = (%lld,%lld,%lld), End = (%lld,%lld,%lld)\n",
            (long long)GridDimension[0], (long long)GridDimension[1], (long long)GridDimension[2],
            (long long)GridStartIndex[0], (long long)GridStartIndex[1], (long long)GridStartIndex[2],
            (long long)GridEndIndex[0], (long long)GridEndIndex[1], (long long)GridEndIndex[2]);
    fprintf(stderr, "NaNReport: GridLeftEdge = (%.17g, %.17g, %.17g), CellWidth = %.17g\n",
            (double)GridLeftEdge[0], (double)GridLeftEdge[1], (double)GridLeftEdge[2],
            (double)CellWidth[0][ci]);
    fprintf(stderr, "NaNReport: D = %.17g  TE = %.17g  GE = %.17g\n",
            (double)BaryonField[DensNum][i], (double)BaryonField[TENum][i],
            DualEnergyFormalism ? (double)BaryonField[GENum][i] : 0.0);
    fprintf(stderr, "NaNReport: v = (%.17g, %.17g, %.17g), vel_max = %.17g\n",
            (double)BaryonField[Vel1Num][i],
            (GridRank > 1) ? (double)BaryonField[Vel2Num][i] : 0.0,
            (GridRank > 2) ? (double)BaryonField[Vel3Num][i] : 0.0,
            (double)vel_max);

    /* Dump the 5-point stencil along each axis so we can see whether the
       corruption is isolated or a front, and which side it came from. */
    int axis_stride[3] = {stride_x, stride_y, stride_z};
    int axis_index[3]  = {ci, cj, ck};
    const char *axis_name[3] = {"x", "y", "z"};
    for (int dim = 0; dim < GridRank; dim++) {
      for (int off = -2; off <= 2; off++) {
        int idx = axis_index[dim] + off;
        if (idx < 0 || idx >= GridDimension[dim]) continue;
        int n = i + off*axis_stride[dim];
        fprintf(stderr, "NaNReport: %s%+lld  D = %.17g  TE = %.17g  GE = %.17g  v = (%.17g, %.17g, %.17g)\n",
                axis_name[dim], (long long)off,
                (double)BaryonField[DensNum][n], (double)BaryonField[TENum][n],
                DualEnergyFormalism ? (double)BaryonField[GENum][n] : 0.0,
                (double)BaryonField[Vel1Num][n],
                (GridRank > 1) ? (double)BaryonField[Vel2Num][n] : 0.0,
                (GridRank > 2) ? (double)BaryonField[Vel3Num][n] : 0.0);
      }
    }

    /* Nearest particle, and nearest star particle -- to test whether this cell
       is inside the 3x3x3 feedback zone of a star (star_feedback6). */
    if (NumberOfParticles > 0 && ParticlePosition[0] != NULL) {
      double best_any = -1.0, best_star = -1.0;
      int    best_any_n = -1, best_star_n = -1;
      for (int n = 0; n < NumberOfParticles; n++) {
        double d2 = 0.0;
        for (int dim = 0; dim < GridRank; dim++) {
          double dd = (double)ParticlePosition[dim][n] - pos[dim];
          d2 += dd*dd;
        }
        if (best_any < 0.0 || d2 < best_any) { best_any = d2; best_any_n = n; }
        if (ParticleType != NULL && ParticleType[n] == PARTICLE_TYPE_STAR)
          if (best_star < 0.0 || d2 < best_star) { best_star = d2; best_star_n = n; }
      }
      if (best_any_n >= 0)
        fprintf(stderr, "NaNReport: nearest particle: type %lld, dist %.17g cells, mass %.17g\n",
                (long long)((ParticleType != NULL) ? ParticleType[best_any_n] : -1),
                sqrt(best_any)/(double)CellWidth[0][ci],
                (ParticleMass != NULL) ? (double)ParticleMass[best_any_n] : 0.0);
      if (best_star_n >= 0)
        fprintf(stderr, "NaNReport: nearest STAR particle: dist %.17g cells, mass %.17g\n",
                sqrt(best_star)/(double)CellWidth[0][ci],
                (ParticleMass != NULL) ? (double)ParticleMass[best_star_n] : 0.0);
      else
        fprintf(stderr, "NaNReport: no star particles on this grid (%lld particles total)\n",
                (long long)NumberOfParticles);
    } else {
      fprintf(stderr, "NaNReport: no particles on this grid\n");
    }
    fprintf(stderr, "===== NaNReport END (proc %lld) =====\n\n", (long long)MyProcessorNumber);
    fflush(stderr);

    ENZO_VFAIL("ApplyBoundsToBaryonFields: %s is NaN or inf at index %"ISYM" -- see NaNReport above.\n",
               bad_field, i)
  }

  // Check for density floor
  if(RestrictDensity && BaryonField[DensNum][i] < rho_floor){
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Density %"ESYM" g/cm^3 replaced with lower bound %"ESYM" g/cm^3\n", BaryonField[DensNum][i]*DensityUnits, rho_floor*DensityUnits);
    }
    BaryonField[DensNum][i] = rho_floor;

    /* Rescale multispecies fields for the floored density!
       If we are at the density floor we are assuming that the gas is totally ionized.
       This is approximately reasonable given that we're likely either in a cosmic void
       or (more likely if this routine is being called) somewhere around an area with very
       vigorous stellar feedback. Even if it's not reasonable, it will get set correctly
       at the next time step by the chemistry solver.
       Also, note that the 0.76 and 0.24 factors in MS=1 come from the cosmic H and He
       fractions; this is slightly off when metals are added, but only by a couple of percent. */
    if(MultiSpecies>0){ // 6-species non-eq (H, He, H2, e-)
        BaryonField[HINum][i] = 0.0;
        BaryonField[HIINum][i] = 0.76*rho_floor;
        BaryonField[HeINum][i] = 0.0;
        BaryonField[HeIINum][i] = 0.0;
        BaryonField[HeIIINum][i] = 0.24*rho_floor;

        /* Electron density (DeNum) is a scaled density field, so there's not a factor of
           electron mass to proton mass to contend with.  i.e., the equation below is NOT a bug!
           Also, this equation does not need to be changed if the assumptions about species densities
           right above this are changed, as it is generally true. */
        BaryonField[DeNum][i] = BaryonField[HIINum][i] + 0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];

      if(MultiSpecies>1){ // 9-species non-eq (adds H2, H2+, H- to MS=1)
        BaryonField[H2INum][i] = BaryonField[H2IINum][i] = BaryonField[HMNum][i] = 0.0; // assuming molecular fraction, H- fraction are zero.

        /* Update electron density for MS=2. Negative sign for HMnum is because H- removes an electron.
           This equation does not need to be changed if the H2I, H2II and HMN assumptions above are
           modified. */
        BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] - BaryonField[HMNum][i];

        if(MultiSpecies>2){ // 12-species non-eq (adds D, DII, HD to MS=2)
          BaryonField[DIINum][i] = 2.0*3.4e-5*BaryonField[HIINum][i]; // D-to-H ratio from Burles & Tytler 1998; factor of 2 is because it's a mass density.
          BaryonField[DINum][i] = BaryonField[HDINum][i] = 0.0;
          BaryonField[DeNum][i] += 0.5*BaryonField[DIINum][i]; // 0.5 takes into account factor of 2 in deuterium mass
        } // if MS>2
      } // if MS>1
    } // if MS>0

  } // if(RestrictDensity && BaryonField[DensNum][i] < rho_floor)

  // Check for density ceiling
  if(RestrictDensity && BaryonField[DensNum][i] > rho_ceiling){
    if(debug){
      printf("Grid::ApplyBoundsToBaryonFields: Density %"ESYM" g/cm^3 replaced with upper bound %"ESYM" g/cm^3\n", BaryonField[DensNum][i]*DensityUnits, rho_ceiling*DensityUnits);
    }
    BaryonField[DensNum][i] = rho_ceiling;

    /* If we are at the density ceiling we are assuming that the gas is close to being
       totally molecular (or totally atomic in the case of MultiSpecies = 1). This is
       heavily dependent on a whole host of assumptions about primordial gas chemistry at
       very high density, but generally is pretty reasonable and consistent with the
       MultiSpecies solvers in Enzo and Grackle. To that end, we assume that the gas is
       99% molecular and 1% atomic, with no ionized gas at all.  As with the assumptions
       around MultiSpecies for the density floor, this is approximately reasonable but any
       inconsistency will ultimately get corrected by the rate solver in the next time step. */
    if(MultiSpecies>0){ // 6-species non-eq (H, He, H2, e-)
        BaryonField[HINum][i] = 0.76*rho_ceiling;  // Assume fully atomic if MS=1; this will be reset soon if MS > 1
        BaryonField[HIINum][i] = 0.0;
        BaryonField[HeINum][i] = 0.24*rho_ceiling;
        BaryonField[HeIINum][i] = 0.0;
        BaryonField[HeIIINum][i] = 0.0;

        /* Electron density (DeNum) is a scaled density field, so there's not a factor of
           electron mass to proton mass to contend with.  i.e., the equation below is NOT a bug!
           Also, this equation does not need to be changed if the assumptions about species densities
           right above this are changed, as it is generally true. */
        BaryonField[DeNum][i] = BaryonField[HIINum][i] + 0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];

      if(MultiSpecies>1){ // 9-species non-eq (adds H2, H2+, H-)
        BaryonField[HINum][i] = 0.01*0.76*rho_ceiling; // reset HI density to 1% of mass density if MS=2
        BaryonField[H2INum][i] = 0.99*0.76*rho_ceiling;  // H2I density set to 99% of mass density
        BaryonField[H2IINum][i] = BaryonField[HMNum][i] = 0.0; // assume H2II and H- are negligible.

        /* Update electron density for MS=2. Negative sign for HMnum is because H- removes an electron.
           This equation does not need to be changed if the H2I, H2II and HMN assumptions above are
           modified, but if you assume in the MS=1 case that the gas is not fully neutral then you
           need to revisit how the electron density is calculated to make sure it's all consistent. */
        BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] - BaryonField[HMNum][i]; // negative sign for HMnum is because H- removes an electron

        if(MultiSpecies>2){ // 12-species non-eq (adds D, DII, HD)
          BaryonField[DINum][i] = 2*3.4e-5*BaryonField[HINum][i]; // D-to-H ratio from Burles & Tytler 1998; factor of 2 is because it's a mass density.
          BaryonField[DIINum][i] = 0.0;
          BaryonField[HDINum][i] = 1.5*3.4e-5*BaryonField[H2INum][i]; // as above, D-to-H from Burles+Tytler; factor of 1.5 comes from mass of HD compared to H2

          BaryonField[DeNum][i] += 0.5*BaryonField[DIINum][i]; // 0.5 takes into account factor of 2 in deuterium mass
        } // if MS>2
      } // if MS>1
    } // if MS>0

  } // if(RestrictDensity && BaryonField[DensNum][i] > rho_ceiling)

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
