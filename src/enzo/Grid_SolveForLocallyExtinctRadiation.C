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
/    the same on the (R^3 smaller) coarse mesh, over compact lists of the
/    nonzero source cells with K_l tabulated by integer squared separation.
/    Grids without sources return early (both views vanish).  Source fields
/    are densities (per unit volume), hence the cell volume in K_l.  Assumes
/    a uniform refinement factor RefineBy.
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

double ReturnWallTime(void);
int FastFourierTransform(float *buffer, int Rank, int DimensionReal[],
			 int Dimension[], int direction, int type);

/* Fine-mesh sum method.  The FFT is used when
     ns * N_cell > LERFFTCostRatio * M log2(M),   M = padded FFT size,
   i.e. when the direct sum's pair count exceeds the FFT's cost by this
   factor (~5 transforms at a few flops per point, vs ~1 flop-equivalent
   per pair; tune from the debug1 timing lines).  0 forces the FFT
   whenever there are sources, < 0 forces the direct sum. */

static const float LERFFTCostRatio = 25.0;

/* Integer division rounding toward -infinity (mesh offsets can be
   negative where the buffer extends past the domain edge). */

static int FloorDiv(int a, int b)
{
  return (a >= 0) ? a/b : -((-a + b - 1)/b);
}

/* Tabulate K = CellVolume / (4 pi max(r^2, (0.38 dx)^2)) by the integer
   squared separation m = r^2/dx^2, for every separation on a Dim mesh.
   Only m = 0 is softened, since 0.38^2 < 1. */

static float *BuildKernelTable(int *Dim, float CellVolume, float dx)
{
  int dim, m, mmax = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    mmax += (Dim[dim]-1)*(Dim[dim]-1);

  float *Ktab = new float[mmax+1];
  float norm = CellVolume / (4.0*pi*dx*dx);
  Ktab[0] = norm / (0.38*0.38);
  for (m = 1; m <= mmax; m++)
    Ktab[m] = norm / float(m);

  return Ktab;
}

/* Gather the cells of a Dim mesh where any of the three fields is positive
   into a compact list (indices si,sj,sk and values S[f]).  Returns the
   number of sources; the arrays are allocated here (NULL if none). */

static int GatherSources(int *Dim, float *Field[3], int **si, int **sj,
                         int **sk, float *S[3])
{
  int i, j, k, n, f, ns = 0;
  int size = Dim[0]*Dim[1]*Dim[2];

  for (n = 0; n < size; n++)
    if (Field[0][n] > 0 || Field[1][n] > 0 || Field[2][n] > 0)
      ns++;

  *si = *sj = *sk = NULL;
  for (f = 0; f < 3; f++)
    S[f] = NULL;
  if (ns == 0)
    return 0;

  *si = new int[ns];
  *sj = new int[ns];
  *sk = new int[ns];
  for (f = 0; f < 3; f++)
    S[f] = new float[ns];

  ns = 0;
  for (k = 0, n = 0; k < Dim[2]; k++)
    for (j = 0; j < Dim[1]; j++)
      for (i = 0; i < Dim[0]; i++, n++)
        if (Field[0][n] > 0 || Field[1][n] > 0 || Field[2][n] > 0) {
          (*si)[ns] = i;
          (*sj)[ns] = j;
          (*sk)[ns] = k;
          for (f = 0; f < 3; f++)
            S[f][ns] = Field[f][n];
          ns++;
        }

  return ns;
}

/* Flux[f] += sum_s S[f][s] * Ktab[|target - source|^2] over a Dim mesh.
   Targets are the outer loop so each flux cell is written once; the
   (j,k) part of the separation is hoisted out of each row. */

static void DirectSum1OverR2(int *Dim, int ns, int *si, int *sj, int *sk,
                             float *S[3], float *Ktab, float *Flux[3])
{
  int i, j, k, n, s, di;
  float K, sum0, sum1, sum2;
  int *dyz2 = new int[ns];

  for (k = 0, n = 0; k < Dim[2]; k++)
    for (j = 0; j < Dim[1]; j++) {

      for (s = 0; s < ns; s++)
        dyz2[s] = (j-sj[s])*(j-sj[s]) + (k-sk[s])*(k-sk[s]);

      for (i = 0; i < Dim[0]; i++, n++) {
        sum0 = sum1 = sum2 = 0;
        for (s = 0; s < ns; s++) {
          di = i - si[s];
          K = Ktab[dyz2[s] + di*di];
          sum0 += S[0][s] * K;
          sum1 += S[1][s] * K;
          sum2 += S[2][s] * K;
        }
        Flux[0][n] += sum0;
        Flux[1][n] += sum1;
        Flux[2][n] += sum2;
      }
    }

  delete [] dyz2;
}

static void DeleteSources(int *si, int *sj, int *sk, float *S[3])
{
  delete [] si;
  delete [] sj;
  delete [] sk;
  for (int f = 0; f < 3; f++)
    delete [] S[f];
}

/* Smallest m >= n with no prime factors other than 2, 3 and 5 (sizes
   every Enzo FFT backend accepts, including ffte). */

static int NextFFTSize(int n)
{
  for (int m = max(n, 1); ; m++) {
    int r = m;
    while (r % 2 == 0) r /= 2;
    while (r % 3 == 0) r /= 3;
    while (r % 5 == 0) r /= 5;
    if (r == 1)
      return m;
  }
}

/* |offset| stored at padded index ip (offsets -(D-1)..D-1 wrap to the top
   of the padded axis), or -1 for padding that no source-target pair hits. */

static int WrapOffset(int ip, int D, int P)
{
  if (ip < D)     return ip;
  if (ip > P - D) return P - ip;
  return -1;
}

/* Same sum as DirectSum1OverR2 (to roundoff), by zero-padded FFT
   convolution on a P mesh (P >= 2 Dim - 1, so nothing wraps around).
   The kernel is even in every offset, so its transform is real and
   leaves real and imaginary parts uncoupled: fields 0 and 1 share one
   complex transform (S0 + i S1) and field 2 uses another. */

static int FFTSum1OverR2(int Rank, int *Dim, int *P, float *Source[3],
                         float *Ktab, float *Flux[3])
{
  int i, j, k, n, np, di, dj, dk;
  int M = P[0]*P[1]*P[2];

  float *Kbuf = new float[2*M]();
  float *A    = new float[2*M]();  // S0 + i S1
  float *B    = new float[2*M]();  // S2

  /* Kernel at every wrapped offset. */

  for (k = 0, np = 0; k < P[2]; k++) {
    dk = WrapOffset(k, Dim[2], P[2]);
    for (j = 0; j < P[1]; j++) {
      dj = WrapOffset(j, Dim[1], P[1]);
      for (i = 0; i < P[0]; i++, np++) {
        di = WrapOffset(i, Dim[0], P[0]);
        if (di >= 0 && dj >= 0 && dk >= 0)
          Kbuf[2*np] = Ktab[di*di + dj*dj + dk*dk];
      }
    }
  }

  /* Sources in the low corner of the padded mesh. */

  for (k = 0, n = 0; k < Dim[2]; k++)
    for (j = 0; j < Dim[1]; j++) {
      np = (k*P[1] + j)*P[0];
      for (i = 0; i < Dim[0]; i++, n++, np++) {
        A[2*np  ] = Source[0][n];
        A[2*np+1] = Source[1][n];
        B[2*np  ] = Source[2][n];
      }
    }

  if (FastFourierTransform(Kbuf, Rank, P, P, FFT_FORWARD,
                           COMPLEX_TO_COMPLEX) == FAIL ||
      FastFourierTransform(A, Rank, P, P, FFT_FORWARD,
                           COMPLEX_TO_COMPLEX) == FAIL ||
      FastFourierTransform(B, Rank, P, P, FFT_FORWARD,
                           COMPLEX_TO_COMPLEX) == FAIL) {
    ENZO_FAIL("Error in FastFourierTransform (LER forward).\n");
  }

  /* Multiply by the (real) kernel transform. */

  float Khat;
  for (np = 0; np < M; np++) {
    Khat = Kbuf[2*np];
    A[2*np  ] *= Khat;
    A[2*np+1] *= Khat;
    B[2*np  ] *= Khat;
    B[2*np+1] *= Khat;
  }
  delete [] Kbuf;

  if (FastFourierTransform(A, Rank, P, P, FFT_INVERSE,
                           COMPLEX_TO_COMPLEX) == FAIL ||
      FastFourierTransform(B, Rank, P, P, FFT_INVERSE,
                           COMPLEX_TO_COMPLEX) == FAIL) {
    ENZO_FAIL("Error in FastFourierTransform (LER inverse).\n");
  }

  /* Add the low corner back into the flux fields. */

  for (k = 0, n = 0; k < Dim[2]; k++)
    for (j = 0; j < Dim[1]; j++) {
      np = (k*P[1] + j)*P[0];
      for (i = 0; i < Dim[0]; i++, n++, np++) {
        Flux[0][n] += A[2*np  ];
        Flux[1][n] += A[2*np+1];
        Flux[2][n] += B[2*np  ];
      }
    }

  delete [] A;
  delete [] B;
  return SUCCESS;
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
 
  int dim, f, i0, j0, k0, n0, n1, ns;
  int *Dim = GravitatingMassFieldDimension;
  int R = RefineBy;
  int size = Dim[0]*Dim[1]*Dim[2];

  float dx  = GravitatingMassFieldCellSize;
  float dxp = dx * R;  // parent cell width
  float CellVolume  = POW(dx,  GridRank);
  float CellVolumeP = POW(dxp, GridRank);
  int *si, *sj, *sk;
  float *S[3];

  /* Compact list of source cells.  With none, both views below vanish
     (sources are non-negative), so there is nothing to do. */

  ns = GatherSources(Dim, Source, &si, &sj, &sk, S);
  if (ns == 0) {
    LCAPERF_STOP("grid_SolveForLocallyExtinctRadiation");
    return SUCCESS;
  }

  /* ------------------------------------------------------------------ */
  /* 1) This level's view: Flux += K_L * S on the fine mesh. */

  double WallTime0 = ReturnWallTime();  // for the debug1 timing line below
  int nsFine = ns;

  /* Direct sum or FFT, whichever is estimated cheaper (see
     LERFFTCostRatio). */

  int Pad[MAX_DIMENSION] = {1, 1, 1};
  for (dim = 0; dim < GridRank; dim++)
    Pad[dim] = NextFFTSize(2*Dim[dim] - 1);
  double M = double(Pad[0])*double(Pad[1])*double(Pad[2]);
  int UseFFT = (LERFFTCostRatio >= 0 &&
                double(ns)*double(size) > LERFFTCostRatio*M*log(M)/log(2.0));

  float *Ktab = BuildKernelTable(Dim, CellVolume, dx);
  if (UseFFT) {
    if (FFTSum1OverR2(GridRank, Dim, Pad, Source, Ktab, Flux) == FAIL) {
      ENZO_FAIL("Error in FFTSum1OverR2.\n");
    }
  } else
    DirectSum1OverR2(Dim, ns, si, sj, sk, S, Ktab, Flux);
  delete [] Ktab;
  DeleteSources(si, sj, sk, S);

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
  int ic, jc, kc, c0;
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

  ns = GatherSources(CDim, Coarse, &si, &sj, &sk, S);
  float *KtabP = BuildKernelTable(CDim, CellVolumeP, dxp);
  DirectSum1OverR2(CDim, ns, si, sj, sk, S, KtabP, CoarseView);
  delete [] KtabP;
  DeleteSources(si, sj, sk, S);
  int nsCoarse = ns;

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

  /* Per-grid cost, to tune LERFFTCostRatio: direct-sum time should scale
     as ns * N_cell, FFT time as M log2(M) (Pad = padded FFT dims), plus
     the R^6-cheaper coarse term either way. */

  if (debug1)
    fprintf(stdout, "LER: P%"ISYM" grid %"ISYM" level %"ISYM
            " Dim %"ISYM" %"ISYM" %"ISYM" ns %"ISYM" (coarse %"ISYM")"
            " %s Pad %"ISYM" %"ISYM" %"ISYM" time %"GSYM" s\n",
            MyProcessorNumber, ID, level,
            Dim[0], Dim[1], Dim[2], nsFine, nsCoarse,
            (UseFFT ? "FFT" : "direct"), Pad[0], Pad[1], Pad[2],
            float(ReturnWallTime() - WallTime0));

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
