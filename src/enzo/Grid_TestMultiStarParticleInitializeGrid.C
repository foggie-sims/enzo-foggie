/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A MULTI-STAR PARTICLE TEST)
/
/  written by: Cassi Lochhaas
/  date:       2026
/  modified1:
/
/  PURPOSE:
/    Places star particles distributed uniformly throughout a sphere
/    centered at Center[], using concentric Fibonacci shells.
/
/    Shells are placed at r = Spacing, 2*Spacing, 3*Spacing, ..., with
/    each shell independently populated using the Fibonacci spiral for
/    near-uniform surface coverage.  The number of particles on shell k
/    is max(1, round(4*pi*k^2)), matching the shell's surface area, so
/    the inter-particle spacing is Spacing throughout the volume.
/    One particle is placed at the center (k=0).
/
/    NParticles sets the minimum number of particles: shells are added
/    until the total meets or exceeds NParticles, with the last shell
/    always completed to preserve the uniform spacing.  Spacing is
/    derived from the user-supplied SphereRadius:
/      Spacing = SphereRadius / n_shells
/    The actual particle count and derived spacing are printed at runtime.
/
/    If UseSphProfile is set, the gas density and temperature outside
/    SphProfileRadius follow a power-law profile:
/      rho(r) = InnerDensity * (SphProfileRadius / r)^SphProfileAlpha
/    with the internal energy adjusted to maintain constant pressure:
/      epsilon(r) = InnerEnergy * (r / SphProfileRadius)^SphProfileAlpha
/    Inside SphProfileRadius the fields are left at the uniform values
/    set by InitializeUniformGrid.
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

int grid::TestMultiStarParticleInitializeGrid(int NParticles,
                                              FLOAT SphereRadius,
                                              FLOAT Center[],
                                              float StarMass,
                                              FLOAT StarVelocity[],
                                              float StarMetallicity,
                                              int UseSphProfile,
                                              FLOAT SphProfileRadius,
                                              float SphProfileAlpha,
                                              float InnerDensity,
                                              float InnerEnergy,
                                              float *Initialdt)
{
  int i, j, dim;

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

  /* Convert star mass to simulation units (mass per cell volume). */

  float ParticleMassCode = StarMass * SolarMass * POW(LengthUnits * CellWidth[0][0], -3.0) / DensityUnits;

  printf("Multi-star particle mass (code units): %f\n", ParticleMassCode);

  /* Add shells using natural filling (N_k = round(4*pi*k^2)) until the
     running total meets or exceeds NParticles.  The last shell is always
     completed so every shell has uniform surface coverage.
     Spacing is derived from SphereRadius once n_shells is known. */

  int n_shells = 0;
  int actual_N = 1;  /* center particle */
  while (actual_N < NParticles) {
    n_shells++;
    int N_k = (int)(4.0 * M_PI * n_shells * n_shells + 0.5);
    if (N_k < 1) N_k = 1;
    actual_N += N_k;
  }

  FLOAT Spacing = (n_shells > 0) ? SphereRadius / n_shells : SphereRadius;

  if (actual_N != NParticles)
    printf("TestMultiStarParticle: requested %d particles; "
           "placing %d to complete last shell\n", NParticles, actual_N);

  /* Build per-shell counts */
  int *N_per_shell = new int[n_shells];
  for (i = 0; i < n_shells; i++) {
    int k = i + 1;
    N_per_shell[i] = (int)(4.0 * M_PI * k * k + 0.5);
    if (N_per_shell[i] < 1) N_per_shell[i] = 1;
  }

  /* Set number of particles and allocate. */

  NumberOfParticles = actual_N;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles on %d shells\n", NumberOfParticles, n_shells);
  printf("Multi-star sphere radius (code units): %"PSYM"\n", SphereRadius);
  printf("Multi-star inter-particle spacing (code units): %"PSYM"\n", Spacing);

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    ParticleType[i] = PARTICLE_TYPE_STAR;
  }

  /* Place particles on concentric Fibonacci shells.
     Shell k (k >= 1) is at radius r = k * Spacing and receives N_per_shell[k-1]
     particles placed with the Fibonacci spiral, covering the full sphere surface.
     One particle is placed at the center (k=0). */

  const double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;
  int count = 0;

  /* Center particle */
  if (GridRank > 0) ParticlePosition[0][count] = Center[0];
  if (GridRank > 1) ParticlePosition[1][count] = Center[1];
  if (GridRank > 2) ParticlePosition[2][count] = Center[2];
  for (dim = 0; dim < GridRank; dim++)
    ParticleVelocity[dim][count] = StarVelocity[dim] * 1e5 * TimeUnits / LengthUnits;
  ParticleMass[count] = ParticleMassCode;
  if (StarMakerStoreInitialMass)
    ParticleInitialMass[count] = ParticleMassCode;
  ParticleAttribute[0][count] = Time + 1e-7;
  ParticleAttribute[1][count] = StarMakerMinimumDynamicalTime * yr_s / TimeUnits;
  if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][count] = 10.0 * Myr_s / TimeUnits;
  if (STARFEED_METHOD(MOM_STAR) || STARFEED_METHOD(MECH_STAR))
    if (StarMakerExplosionDelayTime >= 0.0)
      ParticleAttribute[1][count] = 1.0;
  ParticleAttribute[2][count] = StarMetallicity;
  ParticleAttribute[3][count] = 0.0;
  count++;

  /* Shells */
  int shell;
  for (shell = 1; shell <= n_shells; shell++) {
    double r = shell * (double)Spacing;
    int N_shell = N_per_shell[shell - 1];

    for (j = 0; j < N_shell; j++) {
      double phi   = acos(1.0 - 2.0 * (j + 0.5) / (double)N_shell);
      double theta = 2.0 * M_PI * (double)j / golden_ratio;

      FLOAT dx = r * sin(phi) * cos(theta);
      FLOAT dy = r * sin(phi) * sin(theta);
      FLOAT dz = r * cos(phi);

      if (GridRank > 0) ParticlePosition[0][count] = Center[0] + dx;
      if (GridRank > 1) ParticlePosition[1][count] = Center[1] + dy;
      if (GridRank > 2) ParticlePosition[2][count] = Center[2] + dz;

      for (dim = 0; dim < GridRank; dim++)
        ParticleVelocity[dim][count] = StarVelocity[dim] * 1e5 * TimeUnits / LengthUnits;

      ParticleMass[count] = ParticleMassCode;
      if (StarMakerStoreInitialMass)
        ParticleInitialMass[count] = ParticleMassCode;

      ParticleAttribute[0][count] = Time + 1e-7;
      ParticleAttribute[1][count] = StarMakerMinimumDynamicalTime * yr_s / TimeUnits;
      if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][count] = 10.0 * Myr_s / TimeUnits;
      if (STARFEED_METHOD(MOM_STAR) || STARFEED_METHOD(MECH_STAR))
        if (StarMakerExplosionDelayTime >= 0.0)
          ParticleAttribute[1][count] = 1.0;

      ParticleAttribute[2][count] = StarMetallicity;
      ParticleAttribute[3][count] = 0.0;
      count++;
    }
  }

  delete[] N_per_shell;

  /* Apply spherical atmosphere profile outside SphProfileRadius. */

  if (UseSphProfile) {

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum);

    int k_idx, j_idx, i_idx, index;
    float delx, dely, delz, r, ratio;

    for (k_idx = GridStartIndex[2]; k_idx <= GridEndIndex[2]; k_idx++) {
      delz = (GridRank > 2) ?
        (float)(CellLeftEdge[2][k_idx] + 0.5*CellWidth[2][k_idx] - Center[2]) : 0.0;
      for (j_idx = GridStartIndex[1]; j_idx <= GridEndIndex[1]; j_idx++) {
        dely = (GridRank > 1) ?
          (float)(CellLeftEdge[1][j_idx] + 0.5*CellWidth[1][j_idx] - Center[1]) : 0.0;
        for (i_idx = GridStartIndex[0]; i_idx <= GridEndIndex[0]; i_idx++) {
          delx = (float)(CellLeftEdge[0][i_idx] + 0.5*CellWidth[0][i_idx] - Center[0]);
          r = sqrt(delx*delx + dely*dely + delz*delz);

          index = (k_idx * GridDimension[1] + j_idx) * GridDimension[0] + i_idx;
          if (r > (float)SphProfileRadius) {
            ratio = r / (float)SphProfileRadius;
            float density   = InnerDensity * pow(ratio, -SphProfileAlpha);
            float internalE = InnerEnergy  * pow(ratio,  SphProfileAlpha);
            float vx = BaryonField[Vel1Num][index];
            float vy = (GridRank > 1) ? BaryonField[Vel2Num][index] : 0.0f;
            float vz = (GridRank > 2) ? BaryonField[Vel3Num][index] : 0.0f;
            float totalE = internalE + 0.5f * (vx*vx + vy*vy + vz*vz);

            BaryonField[DensNum][index] = density;
            BaryonField[TENum][index]   = totalE;
            if (DualEnergyFormalism)
              BaryonField[GENum][index] = internalE;
          }
        }
      }
    }
  }

  return SUCCESS;
}
