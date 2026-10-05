/***********************************************************************
/
/  GRID CLASS (SOLVE FOR OPTICALLY-THIN YOUNG-STAR RADIATION ON A SUBGRID)
/
/  written by: Cameron Trapp
/  date:       Sept, 2026
/  modified1:  Oct, 2026 - all three RT fields; parent subtraction made
/              consistent with the prolonged parent flux
/
/  PURPOSE: LEBRON-like RT.  On entry the flux fields hold the parent-level
/    solution, prolonged in PreparePotentialField.  For the sources S on
/    this grid's mesh we then
/      add      K_L * S                      (this level's view)
/      subtract prolong(K_{L-1} * coarsen(S)) (the parent's view, put through
/                                              the same interpolation as the
/                                              parent flux, so it cancels)
/    where K_l = dx_l^3 / (4 pi max(r^2, (0.38 dx_l)^2)), matching the
/    root-grid Green's function.  Over all levels this telescopes to K_L * S.
/
/  NOTE: Direct sums: O(N_source * N_cell) per grid on the fine mesh, and
/    the same on the (R^3 smaller) coarse mesh.  Source fields are densities
/    (per unit volume), hence the cell volume in K_l.  Assumes a uniform
/    refinement factor RefineBy.
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

extern "C" void FORTRAN_NAME(prolong)(float *source, float *dest, int *ndim,
				      int *sdim1, int *sdim2, int *sdim3,
				      int *ddim1, int *ddim2, int *ddim3,
				      int *start1, int *start2, int *start3,
				      int *refine1, int *refine2,int *refine3);

/* Integer division rounding toward -infinity (mesh offsets can be
   negative where the buffer extends past the domain edge). */

static int FloorDiv(int a, int b)
{
  return (a >= 0) ? a/b : -((-a + b - 1)/b);
}
 
int grid::SolveForLocallyExtinctRadiation(int level)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  float *Source[3] = {kdissH2SourceField, kdetHMSourceField, isrfSourceField};
  float *Flux[3]   = {kdissH2FluxField,   kdetHMFluxField,   isrfFluxField};

  for (int f = 0; f < 3; f++)
    if (Source[f] == NULL || Flux[f] == NULL)
      return SUCCESS;  // nothing deposited / no parent solution yet

  LCAPERF_START("grid_SolveForLocallyExtinctRadiation");
 
  int dim, f, i0, j0, k0, i1, j1, k1, n0, n1;
  int *Dim = GravitatingMassFieldDimension;
  int R = RefineBy;
  int size = Dim[0]*Dim[1]*Dim[2];

  float dx  = GravitatingMassFieldCellSize;
  float dxp = dx * R;  // parent cell width
  float CellVolume  = POW(dx,  GridRank);
  float CellVolumeP = POW(dxp, GridRank);
  float soft2  = (0.38*dx)  * (0.38*dx);
  float softp2 = (0.38*dxp) * (0.38*dxp);
  float r2, kernel;

  /* ------------------------------------------------------------------ */
  /* 1) This level's view: Flux += K_L * S on the fine mesh. */

  for (k0 = 0, n0 = 0; k0 < Dim[2]; k0++)
    for (j0 = 0; j0 < Dim[1]; j0++)
      for (i0 = 0; i0 < Dim[0]; i0++, n0++) {

        if (Source[0][n0] <= 0 && Source[1][n0] <= 0 && Source[2][n0] <= 0)
          continue;  // no young stars in this cell

        for (k1 = 0, n1 = 0; k1 < Dim[2]; k1++)
          for (j1 = 0; j1 < Dim[1]; j1++)
            for (i1 = 0; i1 < Dim[0]; i1++, n1++) {
              r2 = ((k1-k0)*(k1-k0) + (j1-j0)*(j1-j0) + (i1-i0)*(i1-i0)) *
                   dx * dx;
              kernel = CellVolume / (4.0*pi*max(r2, soft2));
              for (f = 0; f < 3; f++)
                Flux[f][n1] += Source[f][n0] * kernel;
            }
      }

  /* ------------------------------------------------------------------ */
  /* 2) The parent's view of the same sources: coarsen S onto a temporary
     parent-resolution mesh laid out exactly as PreparePotentialField lays
     out the parent region it prolongs from, convolve with K_{L-1}, prolong
     back with the same call, and subtract. */

  int Offset[MAX_DIMENSION], CStart[MAX_DIMENSION], CDim[MAX_DIMENSION],
      ProlongStart[MAX_DIMENSION], Refine[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Offset[dim] = CStart[dim] = ProlongStart[dim] = 0;
    CDim[dim] = Refine[dim] = 1;
  }
  for (dim = 0; dim < GridRank; dim++) {
    Offset[dim] = nint((GravitatingMassFieldLeftEdge[dim] -
                        DomainLeftEdge[dim]) / dx);  // global fine index
    CStart[dim] = FloorDiv(Offset[dim], R) - 1;
    CDim[dim]   = FloorDiv(Offset[dim] + Dim[dim] - 1, R) - CStart[dim] + 3;
    ProlongStart[dim] = Offset[dim] - R*CStart[dim];
    Refine[dim] = R;
  }

  int csize = CDim[0]*CDim[1]*CDim[2];
  float *Coarse[3], *CoarseView[3];
  for (f = 0; f < 3; f++) {
    Coarse[f]     = new float[csize]();
    CoarseView[f] = new float[csize]();
  }

  /* Coarsen: parent-cell density = sum of fine densities / R^rank. */

  float inv_R_rank = 1.0 / POW(float(R), GridRank);
  int ic, jc, kc, c0, c1;
  for (k0 = 0, n0 = 0; k0 < Dim[2]; k0++) {
    kc = (GridRank > 2) ? FloorDiv(k0 + Offset[2], R) - CStart[2] : 0;
    for (j0 = 0; j0 < Dim[1]; j0++) {
      jc = (GridRank > 1) ? FloorDiv(j0 + Offset[1], R) - CStart[1] : 0;
      for (i0 = 0; i0 < Dim[0]; i0++, n0++) {
        ic = FloorDiv(i0 + Offset[0], R) - CStart[0];
        c0 = (kc*CDim[1] + jc)*CDim[0] + ic;
        for (f = 0; f < 3; f++)
          Coarse[f][c0] += Source[f][n0] * inv_R_rank;
      }
    }
  }

  /* Convolve with K_{L-1} on the coarse mesh. */

  for (k0 = 0, c0 = 0; k0 < CDim[2]; k0++)
    for (j0 = 0; j0 < CDim[1]; j0++)
      for (i0 = 0; i0 < CDim[0]; i0++, c0++) {

        if (Coarse[0][c0] <= 0 && Coarse[1][c0] <= 0 && Coarse[2][c0] <= 0)
          continue;

        for (k1 = 0, c1 = 0; k1 < CDim[2]; k1++)
          for (j1 = 0; j1 < CDim[1]; j1++)
            for (i1 = 0; i1 < CDim[0]; i1++, c1++) {
              r2 = ((k1-k0)*(k1-k0) + (j1-j0)*(j1-j0) + (i1-i0)*(i1-i0)) *
                   dxp * dxp;
              kernel = CellVolumeP / (4.0*pi*max(r2, softp2));
              for (f = 0; f < 3; f++)
                CoarseView[f][c1] += Coarse[f][c0] * kernel;
            }
      }

  /* Prolong to the fine mesh (same routine and relative offsets as the
     parent flux) and subtract. */

  float *ParentView = new float[size];
  for (f = 0; f < 3; f++) {
    FORTRAN_NAME(prolong)(CoarseView[f], ParentView, &GridRank,
                          CDim, CDim+1, CDim+2,
                          Dim, Dim+1, Dim+2,
                          ProlongStart, ProlongStart+1, ProlongStart+2,
                          Refine, Refine+1, Refine+2);
    for (n1 = 0; n1 < size; n1++)
      Flux[f][n1] -= ParentView[n1];
  }

  delete [] ParentView;
  for (f = 0; f < 3; f++) {
    delete [] Coarse[f];
    delete [] CoarseView[f];
  }
 
  LCAPERF_STOP("grid_SolveForLocallyExtinctRadiation");
  return SUCCESS;
}


/* Copy the LEBRON-like flux fields (GravitatingMassField mesh) into
   BaryonField-shaped arrays for Grackle, converting units.  The flux
   fields hold sum(table [rate*cm^2/Msun] * Msun) / (code length)^2, so
   dividing by LengthUnits^2 gives 1/s (or Habing for the ISRF).  The
   offset into the gravity mesh follows CopyPotentialToBaryonField.
   Output arrays are left untouched where no flux field exists. */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::GetLocallyExtinctRadiationRates(float *kdissH2, float *kdetHM,
                                          float *isrf)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
           &VelocityUnits, Time);

  /* TODO: in comoving runs LengthUnits is comoving; check whether a
     factor of a^2 is needed here. */

  float inv_L2 = 1.0 / (LengthUnits * LengthUnits);
  float *Flux[]   = {kdissH2FluxField, kdetHMFluxField, isrfFluxField};
  float *Out[]    = {kdissH2, kdetHM, isrf};
  float Factor[]  = {TimeUnits*inv_L2, TimeUnits*inv_L2, inv_L2};

  int dim, i, j, k, f, index, n;
  int Off[3] = {0, 0, 0};
  for (dim = 0; dim < GridRank; dim++)
    Off[dim] = (GravitatingMassFieldDimension[dim] - GridDimension[dim])/2;

  for (f = 0; f < 3; f++) {
    if (Flux[f] == NULL || Out[f] == NULL)
      continue;
    n = 0;
    for (k = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++) {
        index = (((k+Off[2])*GravitatingMassFieldDimension[1]) +
                 (j+Off[1]))*GravitatingMassFieldDimension[0] + Off[0];
        for (i = 0; i < GridDimension[0]; i++, index++, n++)
          Out[f][n] = max(Flux[f][index], 0.0) * Factor[f];
      }
  }

  return SUCCESS;
}


/* Allocate (if needed) and zero the LEBRON-like RT source arrays.
   ClearRTSourceParticles mirrors ClearGravitatingMassFieldParticles
   (particle mesh); ClearRTSourceField mirrors ClearGravitatingMassField
   (gravity mesh).  Both require the corresponding mesh to be initialized. */

int grid::ClearRTSourceParticles()
{
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (GravitatingMassFieldParticlesCellSize == FLOAT_UNDEFINED)
    ENZO_FAIL("GravitatingMassFieldParticles uninitialized.\n");

  int dim, i, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldParticlesDimension[dim];

  if (kdissH2SourceParticles == NULL) kdissH2SourceParticles = new float[size];
  if (kdetHMSourceParticles  == NULL) kdetHMSourceParticles  = new float[size];
  if (isrfSourceParticles    == NULL) isrfSourceParticles    = new float[size];

  for (i = 0; i < size; i++) {
    kdissH2SourceParticles[i] = 0.0;
    kdetHMSourceParticles[i]  = 0.0;
    isrfSourceParticles[i]    = 0.0;
  }

  return SUCCESS;
}

int grid::ClearRTSourceField()
{
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (GravitatingMassFieldCellSize == FLOAT_UNDEFINED)
    ENZO_FAIL("GravitatingMassField uninitialized.\n");

  int dim, i, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldDimension[dim];

  if (kdissH2SourceField == NULL) kdissH2SourceField = new float[size];
  if (kdetHMSourceField  == NULL) kdetHMSourceField  = new float[size];
  if (isrfSourceField    == NULL) isrfSourceField    = new float[size];

  for (i = 0; i < size; i++) {
    kdissH2SourceField[i] = 0.0;
    kdetHMSourceField[i]  = 0.0;
    isrfSourceField[i]    = 0.0;
  }

  return SUCCESS;
}
