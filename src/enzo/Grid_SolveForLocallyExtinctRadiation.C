/***********************************************************************
/
/  GRID CLASS (SOLVE FOR OPTICALLY-THIN YOUNG-STAR RADIATION ON A SUBGRID)
/
/  written by: Cameron Trapp
/  date:       Sept, 2026
/  modified1:  Oct, 2026 - all three RT fields; parent subtraction made
/              consistent with the prolonged parent flux
/
/  PURPOSE: Locally-Extinct-Stellar-Radiation. The flux fields initially
/    hold the parent-level solution, prolonged in PreparePotentialField.
/    For the sources S on this grid's mesh we then:
/      add      K_L * S                      (this level's view)
/      subtract prolong(K_{L-1} * coarsen(S)) (the parent's view, put through
/                                              the same interpolation as the
/                                              parent flux, so it cancels)
/    where K_l = dx_l^3 / (4 pi max(r^2, (0.38 dx_l)^2)).
/
/  NOTE: Direct sums: O(N_source * N_cell) per grid on the fine mesh, and
/    the same on the (R^3 smaller) coarse mesh. Source fields are densities,
/    hence the cell volume in K_l. Assumes a uniform refinement factor RefineBy.
/
/    Depending on the number of sources, this will compute this sum in Fourier
/    space, O(M log M), where M~8*N_cell.
/
/    With LocallyExtinctRadiationBoxRadius = rbox > 0, the correction is applied
/    only within rbox fine cells of each source, O(N_source * (2 rbox + 1)^3).
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

static int FloorDiv(int a, int b)
{
  /* Floor division that rounds negative indices down */
  return int(floor(double(a)/double(b)));
}


static float *BuildKernelTable(int *Dim, float CellVolume, float dx)
{
  /* Tabulate r2tab = CellVolume / (4 pi max(r^2, (0.38 dx)^2)) by the integer
   squared separation r2idx = r^2/dx^2, for every separation on a Dim mesh.
   Only r2idx = 0 is softened*/
  int dim, r2idx, max_idx = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    max_idx += (Dim[dim]-1)*(Dim[dim]-1);

  float *r2tab = new float[max_idx+1];
  float norm = CellVolume / (4.0*pi*dx*dx);
  r2tab[0] = norm / (0.38*0.38);
  for (r2idx = 1; r2idx <= max_idx; r2idx++)
    r2tab[r2idx] = norm / float(r2idx);

  return r2tab;
}

static int GatherSources(int *Dim, float *Field[3], int **si, int **sj,
                         int **sk, float *S[3])
{
/* Gather the cells of a Dim mesh where any of the three fields is positive
   into a compact list (indices si,sj,sk and values S[f]).  Returns the
   number of sources; the arrays are allocated here (NULL if none). */

  int i, j, k, n, f, ns = 0;
  int size = Dim[0]*Dim[1]*Dim[2];

  for (n = 0; n < size; n++) //How many sources are there
    if (Field[0][n] > 0 || Field[1][n] > 0 || Field[2][n] > 0)
      ns++;

  *si = *sj = *sk = NULL; //initialize
  for (f = 0; f < 3; f++)
    S[f] = NULL;
  if (ns == 0)
    return 0;

  *si = new int[ns];
  *sj = new int[ns];
  *sk = new int[ns];
  for (f = 0; f < 3; f++)
    S[f] = new float[ns];

  ns = 0; //store source indices and fluxes
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

static void DirectSumOverR2(int *Dim, int ns, int *si, int *sj, int *sk,
                             float *S[3], float *r2tab, float *Flux[3])
{
  /* 
  Directly calculate the 1/r^2 term for all sources and sum
  Flux[field] += Sum_sources S[f][source] * r2tab[(target-source)^2]
  */

  int i, j, k, n, s, di;
  float r2term, sum0, sum1, sum2;
  int *dyz2 = new int[ns];

  for (k = 0, n = 0; k < Dim[2]; k++)
    for (j = 0; j < Dim[1]; j++) { //loop over outer 2 grid dimensions

      for (s = 0; s < ns; s++) //loop over sources
        dyz2[s] = (j-sj[s])*(j-sj[s]) + (k-sk[s])*(k-sk[s]);

      for (i = 0; i < Dim[0]; i++, n++) { //loop over inner grid dimension
        sum0 = sum1 = sum2 = 0;
        for (s = 0; s < ns; s++) {
          di = i - si[s];
          r2term = r2tab[dyz2[s] + di*di];
          sum0 += S[0][s] * r2term; //H2 photodissociation
          sum1 += S[1][s] * r2term; //HM photodetachment
          sum2 += S[2][s] * r2term; //isrf
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
  /* Delete saved source indices and fluxes */
  delete [] si;
  delete [] sj;
  delete [] sk;
  for (int f = 0; f < 3; f++)
    delete [] S[f];
}

static int NextFFTSize(int n)
{
  /* Smallest m >= n with no prime factors other than 2, 3 and 5 (sizes
   every Enzo FFT backend accepts, including ffte). */
  for (int m = max(n, 1); ; m++) {
    int r = m;
    while (r % 2 == 0) r /= 2;
    while (r % 3 == 0) r /= 3;
    while (r % 5 == 0) r /= 5;
    if (r == 1)
      return m;
  }
}

static int WrapOffset(int ip, int D, int P)
{
  /* |offset| stored at padded index ip (offsets -(D-1)..D-1 wrap to the top
   of the padded axis), or -1 for padding that no source-target pair hits. */
  if (ip < D)     return ip;
  if (ip > P - D) return P - ip;
  return -1;
}

static int FFTSumOverR2(int Rank, int *Dim, int *P, float *Source[3],
                         float *r2tab, float *Flux[3])
{
  /* Same sum as DirectSumOverR2 (to roundoff), by zero-padded FFT
   convolution on a P mesh (P >= 2 Dim - 1, so nothing wraps around).
   The kernel is even in every offset, so its transform is real and
   leaves real and imaginary parts uncoupled: fields 0 and 1 share one
   complex transform (S0 + i S1) and field 2 uses another. */
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
          Kbuf[2*np] = r2tab[di*di + dj*dj + dk*dk];
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


static float *BuildBoxCorrection(int rbox, float dx, float CellVolume,
                                 float CellVolumeParent)
{
  /* Table for the box correction (LocallyExtinctRadiationBoxRadius = rbox > 0).
   3D only (checked in ReadParameterFile).

   Each entry is what one fine source cell of unit density adds to a target
   offset = (di, dj, dk) fine cells away, for |di|, |dj|, |dk| <= rbox:

     correction = this level's flux - the parent's flux interpolated to the
                  target
                = normFine / |offset|^2
                  - sum of weight * normParent / |cell|^2 over the 8 parent
                    cells around the target (trilinear, as in prolong.F)

   The parent term depends on where the source sits inside its parent cell
   (subpos = 0..RefineBy-1 per dimension), since that sets how the target
   lines up with the parent cell centers.  So the table has one block of
   (2 rbox + 1)^3 entries for each subpos.  Both terms are softened to
   r = 0.38 cell at zero separation.

   Layout: block for subpos (subpos[0] fastest), then offset (offset[0]
   fastest, each -rbox..rbox). */

  int BoxWidth = 2*rbox + 1;
  int BoxSize  = BoxWidth*BoxWidth*BoxWidth;
  float *BoxTable = new float[RefineBy*RefineBy*RefineBy*BoxSize];

  float dxParent = dx*RefineBy;
  float normFine = CellVolume / (4.0*pi*dx*dx);
  /* A unit-density fine cell is a parent density of 1/RefineBy^3. */
  float normParent = CellVolumeParent / (4.0*pi*dxParent*dxParent) /
                     float(RefineBy*RefineBy*RefineBy);
  float soft2 = 0.38*0.38;  // softening length squared, in cells

  int dim, subpos[3], offset[3], below[3];
  int r2Fine, r2Parent, entry = 0;
  float FineCenter, ParentCoord, wAbove[3], weight, ParentFlux;

  for (subpos[2] = 0; subpos[2] < RefineBy; subpos[2]++) //The 8 positions the fine resolution cell could be in the parent cell
  for (subpos[1] = 0; subpos[1] < RefineBy; subpos[1]++)
  for (subpos[0] = 0; subpos[0] < RefineBy; subpos[0]++)
    for (offset[2] = -rbox; offset[2] <= rbox; offset[2]++) //All cells in the box
    for (offset[1] = -rbox; offset[1] <= rbox; offset[1]++)
    for (offset[0] = -rbox; offset[0] <= rbox; offset[0]++, entry++) {

      /* This level's view. */

      r2Fine = offset[0]*offset[0] + offset[1]*offset[1] + offset[2]*offset[2];
      BoxTable[entry] = normFine / max(float(r2Fine), soft2);

      /* Parent's view.  Locate the target among the parent cell centers,
         one direction at a time, with the source's parent cell as parent
         cell 0. */

      for (dim = 0; dim < 3; dim++) {

        /* Target's center in fine cells, counted from the left edge of
           the source's parent cell (fine cells 0..RefineBy-1 lie in it). */

        FineCenter = subpos[dim] + offset[dim] + 0.5;

        /* The same position in parent cells, shifted by half a cell so
           parent cell p has its center at p. */

        ParentCoord = FineCenter / RefineBy - 0.5;

        /* below: the parent cell centered at or just below the target.
           wAbove: how far the target is toward the next center up (0-1),
           which is the interpolation weight on that cell. */

        below[dim]  = int(floor(ParentCoord));
        wAbove[dim] = ParentCoord - below[dim];
      }

      /* Interpolate the parent's flux to the target (trilinear, as
         prolong.F does). */

      float w[3][2];  // w[dim][0]: lower cell's weight, w[dim][1]: upper cell's
      for (dim = 0; dim < 3; dim++) {
        w[dim][0] = 1 - wAbove[dim]; //Lower cell's weight
        w[dim][1] = wAbove[dim]; //Upper cell's weight
      }

      ParentFlux = 0;
      for (int uz = 0; uz <= 1; uz++)      // 0 = lower, 1 = upper parent cell
      for (int uy = 0; uy <= 1; uy++)
      for (int ux = 0; ux <= 1; ux++) {
        int cx = below[0] + ux;
        int cy = below[1] + uy;
        int cz = below[2] + uz;
        weight   = w[0][ux] * w[1][uy] * w[2][uz];
        r2Parent = cx*cx + cy*cy + cz*cz;
        ParentFlux += weight * normParent / max(float(r2Parent), soft2);
      }

      BoxTable[entry] -= ParentFlux;
    }

  return BoxTable;
}

static void BoxSumOverR2(int *Dim, int *Offset, int rbox,
                          int ns, int *si, int *sj, int *sk, float *S[3],
                          float *BoxTable, float *Flux[3])
{
  /* Add the box correction around each source cell: the table block for
     the source's sub-position, times the source value, over the cells
     within rbox of it (clipped to the mesh). */

  int BoxWidth = 2*rbox + 1;
  int BoxSize  = BoxWidth*BoxWidth*BoxWidth;

  int src, dim, i, j, k, n, ilo, ihi, jlo, jhi, klo, khi;
  int SrcIndex[3], GlobalIndex, ParentIndex, subpos[3], BoxTableStart, entry;
  float S0, S1, S2;

  for (src = 0; src < ns; src++) {

    /* Which block of the table this source uses: its position inside its
       parent cell along each axis (0..RefineBy-1). */

    SrcIndex[0] = si[src];
    SrcIndex[1] = sj[src];
    SrcIndex[2] = sk[src];
    for (dim = 0; dim < 3; dim++) {
      GlobalIndex = SrcIndex[dim] + Offset[dim];       // fine index from the domain's left edge
      ParentIndex = FloorDiv(GlobalIndex, RefineBy);   // parent cell it lies in 
      subpos[dim] = GlobalIndex - RefineBy*ParentIndex; // position inside that parent cell
    }
    BoxTableStart = ((subpos[2]*RefineBy + subpos[1])*RefineBy + subpos[0])*BoxSize;

    /* The box, clipped to the mesh. */

    ilo = max(si[src] - rbox, 0);  ihi = min(si[src] + rbox, Dim[0]-1);
    jlo = max(sj[src] - rbox, 0);  jhi = min(sj[src] + rbox, Dim[1]-1);
    klo = max(sk[src] - rbox, 0);  khi = min(sk[src] + rbox, Dim[2]-1);

    S0 = S[0][src];
    S1 = S[1][src];
    S2 = S[2][src];
    for (k = klo; k <= khi; k++)
      for (j = jlo; j <= jhi; j++)
        for (i = ilo; i <= ihi; i++) {

          /* Target cell in the mesh, and the table entry for its offset
             from the source (each offset -rbox:rbox, shifted to 0:2 rbox). */

          n     = (k*Dim[1] + j)*Dim[0] + i;
          entry = BoxTableStart +
                  ((k - sk[src] + rbox)*BoxWidth + (j - sj[src] + rbox))*BoxWidth +
                  (i - si[src] + rbox);

          Flux[0][n] += S0 * BoxTable[entry];
          Flux[1][n] += S1 * BoxTable[entry];
          Flux[2][n] += S2 * BoxTable[entry];
        }
  }
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
      return SUCCESS;  // nothing deposited / no parent fluxes

  LCAPERF_START("grid_SolveForLocallyExtinctRadiation");
 
  int dim, f, i0, j0, k0, n0, n1, ns;
  int *Dim = GravitatingMassFieldDimension;
  int R = RefineBy;
  int size = Dim[0]*Dim[1]*Dim[2];

  float dx  = GravitatingMassFieldCellSize;
  float dxp = dx * R;  // parent cell width
  float CellVolume  = POW(dx,  GridRank);
  float CellVolumeParent = POW(dxp, GridRank);
  int *si, *sj, *sk;
  float *S[3];

  /* Compact list of source cells. (Luminosities have already been deposited) */

  ns = GatherSources(Dim, Source, &si, &sj, &sk, S);
  if (ns == 0) { //No young stars on this grid
    LCAPERF_STOP("grid_SolveForLocallyExtinctRadiation");
    return SUCCESS;
  }

  double WallTime0 = ReturnWallTime();  // for the debug1 timing lines below
  int nsFine = ns;

  /* Alignment with the parent mesh: the parent region this grid's flux was
     prolonged from, laid out as in PreparePotentialField. */

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

  /* ------------------------------------------------------------------ */
  /* 0) Box-truncated correction instead of steps 1 and 2 (see
     BuildBoxCorrection). */

  if (LocallyExtinctRadiationBoxRadius > 0) {
    int rbox = LocallyExtinctRadiationBoxRadius;
    float *BoxTable = BuildBoxCorrection(rbox, dx, CellVolume,
                                         CellVolumeParent);
    BoxSumOverR2(Dim, Offset, rbox, ns, si, sj, sk, S, BoxTable, Flux);
    delete [] BoxTable;
    DeleteSources(si, sj, sk, S);

    if (debug1)
      fprintf(stdout, "LER: P%"ISYM" grid %"ISYM" level %"ISYM
              " Dim %"ISYM" %"ISYM" %"ISYM" ns %"ISYM" box rbox %"ISYM
              " time %"GSYM" s\n", MyProcessorNumber, ID, level,
              Dim[0], Dim[1], Dim[2], nsFine, rbox,
              float(ReturnWallTime() - WallTime0));

    LCAPERF_STOP("grid_SolveForLocallyExtinctRadiation");
    return SUCCESS;
  }

  /* ------------------------------------------------------------------ */
  /* 1) This level's view: Flux += K_L * S on the fine mesh. */

  /* Direct sum or FFT, whichever is estimated cheaper */

  int Pad[MAX_DIMENSION] = {1, 1, 1};
  for (dim = 0; dim < GridRank; dim++)
    Pad[dim] = NextFFTSize(2*Dim[dim] - 1);
  double M = double(Pad[0])*double(Pad[1])*double(Pad[2]);
  int UseFFT = (LocallyExtinctRadiationFFTRatio >= 0 &&
                double(ns)*double(size) >
                LocallyExtinctRadiationFFTRatio*M*log(M)/log(2.0));

  float *r2tab = BuildKernelTable(Dim, CellVolume, dx);
  if (UseFFT) {
    if (FFTSumOverR2(GridRank, Dim, Pad, Source, r2tab, Flux) == FAIL) {
      ENZO_FAIL("Error in FFTSumOverR2.\n");
    }
  } else
    DirectSumOverR2(Dim, ns, si, sj, sk, S, r2tab, Flux);
  delete [] r2tab;
  DeleteSources(si, sj, sk, S);

  /* ------------------------------------------------------------------ */
  /* 2) The parent's view of the same sources: coarsen S onto a temporary
     parent-resolution mesh laid out exactly as PreparePotentialField lays
     out the parent region it prolongs from, convolve with K_{L-1}, prolong
     back with the same call, and subtract. */

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
  float *r2tabP = BuildKernelTable(CDim, CellVolumeParent, dxp);
  DirectSumOverR2(CDim, ns, si, sj, sk, S, r2tabP, CoarseView);
  delete [] r2tabP;
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

  /* Per-grid cost, to tune LocallyExtinctRadiationFFTRatio: direct-sum time
     should scale
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


/* Copy the LER flux fields (GravitatingMassField mesh) into
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

  float *Flux[]   = {kdissH2FluxField, kdetHMFluxField, isrfFluxField};
  float *Out[]    = {kdissH2, kdetHM, isrf};
  float UnitConversion[]  = {TimeUnits/LengthUnits/LengthUnits, TimeUnits/LengthUnits/LengthUnits, 1.0/LengthUnits/LengthUnits};

  int dim, i, j, k, f, index, n;
  int Off[3] = {0, 0, 0};
  for (dim = 0; dim < GridRank; dim++) //offset between grid and gravitating mass fields
    Off[dim] = (GravitatingMassFieldDimension[dim] - GridDimension[dim])/2;

  for (f = 0; f < 3; f++) { //loop over bands
    if (Flux[f] == NULL || Out[f] == NULL)
      continue;
    n = 0;
    for (k = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++) {
        index = (((k+Off[2])*GravitatingMassFieldDimension[1]) +
                 (j+Off[1]))*GravitatingMassFieldDimension[0] + Off[0];
        for (i = 0; i < GridDimension[0]; i++, index++, n++)
          Out[f][n] = max(Flux[f][index], 0.0) * UnitConversion[f];
      }
  }

  return SUCCESS;
}


/* Allocate (if needed) and zero the Locally-Extinct-Stellar-Radiation source arrays.
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
