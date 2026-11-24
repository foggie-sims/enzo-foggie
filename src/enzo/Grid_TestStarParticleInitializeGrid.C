/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A STAR PARTICLE TEST)
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass, 
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
					 FLOAT TestStarParticleStarPosition[],
					 float TestStarParticleMetallicity)
{
  /* declarations */

  float CentralMass = 1.0;
  int i, dim;
  float TestInitialdt = *Initialdt;
  
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  
  /* Get Units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1;
  double MassUnits = 1;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Set Central Mass in simulation units */

  CentralMass = TestStarParticleStarMass*SolarMass* POW(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;

  printf("Central Mass: %f \n",CentralMass);

  /* Set number of particles for this grid and allocate space. */

  NumberOfParticles = 1;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    ParticleType[i] = PARTICLE_TYPE_STAR;
  }

  /* Set central particle. */ 
  for (dim = 0; dim < GridRank; dim++) {
    ParticlePosition[dim][0] = TestStarParticleStarPosition[dim];
    ParticleVelocity[dim][0] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
  }
  ParticleMass[0] = CentralMass;
  if (StarMakerStoreInitialMass)
    ParticleInitialMass[0] = CentralMass;
  ParticleAttribute[0][0] = Time+1e-7;
  ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;

  ParticleAttribute[1][0] = StarMakerMinimumDynamicalTime * yr_s/TimeUnits; // default
  if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;
  if (STARFEED_METHOD(MOM_STAR) || STARFEED_METHOD(MECH_STAR))
    if(StarMakerExplosionDelayTime >= 0.0)
      ParticleAttribute[1][0] = 1.0;
  
  ParticleAttribute[2][0] = TestStarParticleMetallicity;  // Metal fraction
  ParticleAttribute[3][0] = 0.0;  // metalfSNIa

  return SUCCESS;
}


/* Optionally split the medium at x_SplitDim = 0.5 (code units) into two
   halves of different density: SplitDensity[0] for cells whose center is
   below 0.5, SplitDensity[1] above. All density-type fields (species,
   metals) are rescaled so their mass fractions are unchanged. If
   PressureEquilibrium is on, the specific internal energy set by
   TestStarParticleEnergy/Temperature is kept on the denser half and the
   other half is heated so both halves have the same pressure; otherwise
   both halves keep the same specific internal energy. Velocities are
   unchanged. */

int grid::TestStarParticleSplitMedium(int SplitDim, float SplitDensity[2],
                                      int PressureEquilibrium)
{
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  if (SplitDim < 0)
    return SUCCESS;
  if (SplitDim >= GridRank)
    ENZO_FAIL("TestStarParticleSplitPlane must be less than TopGridRank.\n");
  if (SplitDensity[0] <= 0 || SplitDensity[1] <= 0)
    ENZO_FAIL("TestStarParticleSplitDensity values must be positive.\n");

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL)
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");

  float dmax = max(SplitDensity[0], SplitDensity[1]);
  int i, j, k, n, field, side, idx[3];
  float ratio, ke, eint, vel[3];
  FLOAT xc;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++) {
        n = i + GridDimension[0]*(j + GridDimension[1]*k);
        idx[0] = i; idx[1] = j; idx[2] = k;
        xc = CellLeftEdge[SplitDim][idx[SplitDim]] +
             0.5*CellWidth[SplitDim][idx[SplitDim]];
        side = (xc < 0.5) ? 0 : 1;

        /* specific kinetic and internal energy before the change */
        vel[0] = BaryonField[Vel1Num][n];
        vel[1] = (GridRank > 1) ? BaryonField[Vel2Num][n] : 0.0;
        vel[2] = (GridRank > 2) ? BaryonField[Vel3Num][n] : 0.0;
        ke = 0.5*(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
        eint = (DualEnergyFormalism) ? BaryonField[GENum][n]
                                     : BaryonField[TENum][n] - ke;

        /* rescale every density-type field (density, species, metals) */
        ratio = SplitDensity[side] / BaryonField[DensNum][n];
        for (field = 0; field < NumberOfBaryonFields; field++)
          if (FieldTypeIsDensity(FieldType[field]))
            BaryonField[field][n] *= ratio;

        /* equal pressure: dense side keeps eint, other side scaled up */
        if (PressureEquilibrium)
          eint *= dmax / SplitDensity[side];

        if (DualEnergyFormalism)
          BaryonField[GENum][n] = eint;
        BaryonField[TENum][n] = eint + ke;
      }

  return SUCCESS;
}