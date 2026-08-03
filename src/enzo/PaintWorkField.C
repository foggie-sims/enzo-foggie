/***********************************************************************
/
/  PAINT THE PER-CELL WORK FIELD FOR A WORK-CAPPED BERGER-RIGOUTSOS
/
/  written by: performance audit (T2.1 / T0.3 successor)
/  date:       August 2026
/
/  PURPOSE: Before a level is rebuilt, give each parent grid a per-cell
/           map of what the region cost last time, so that
/           ProtoSubgrid::AcceptableSubgrid can cap a candidate by
/           predicted WORK rather than only by cell count.
/
/           Why this is needed at all: the cell cap works as designed but
/           predicts an individual grid's cost only to within ~2.4x, so
/           the largest grid by work is not the largest by cells, and a
/           few grids end up holding several ranks' worth of a level's
/           work.  A grid cannot be split across ranks, so no assignment
/           policy can repair that - it has to be fixed where grids are
/           made.  Simulated on measured work, capping at one rank's even
/           share takes the achievable imbalance from 3.02x to 1.04x for
/           about 5% more grids.
/
/           Why the reduction is unavoidable: measured work lives on the
/           outgoing subgrids, which are scattered across ranks (only 9
/           of 350 level-10 grids sit on the rank owning the root tile
/           they occupy), while the parent may live anywhere.  Every
/           locally available alternative was regressed against measured
/           work and none beat the cell count - the cooling rate, the
/           physically motivated candidate, was far worse.  Grid
/           structure must also be identical on every rank or the run
/           deadlocks, so a partial local estimate is not an option.
/
/           Grid EXTENTS are already replicated, so only the work values
/           are communicated: each rank reports its own grids and zero
/           elsewhere, then one sum makes the array agree everywhere.
/           This is placed among the collectives the rebuild already
/           performs, where ranks are synchronised anyway, so it should
/           cost message latency rather than a fresh barrier.  That is a
/           claim about placement, so it is timed: set GridWorkMapOutput
/           to report the cost.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"

double ReturnWallTime(void);

/* Running cost of the work-field reduction, reported once per root cycle
   when the work map is enabled.  Nothing depends on these; they exist so
   the overhead is measured rather than assumed. */

double WorkFieldReduceTime = 0.0;
int    WorkFieldReduceCount = 0;

int PaintWorkField(LevelHierarchyEntry *ParentLevel,
		   LevelHierarchyEntry *OldChildLevel, int level)
{

  MaximumSubgridWork = 0.0;

  if (SubgridMaximumWorkFraction <= 0)
    return SUCCESS;
  if (ParentLevel == NULL || OldChildLevel == NULL)
    return SUCCESS;

  /* Count the outgoing children.  Their extents are already known on
     every rank; only their measured cost is not. */

  int nchild = 0;
  LevelHierarchyEntry *Temp;
  for (Temp = OldChildLevel; Temp; Temp = Temp->NextGridThisLevel)
    nchild++;
  if (nchild == 0)
    return SUCCESS;

  float *ChildWork = new float[nchild];
  int i = 0;
  for (Temp = OldChildLevel; Temp; Temp = Temp->NextGridThisLevel, i++)
    ChildWork[i] =
      (Temp->GridData->ReturnProcessorNumber() == MyProcessorNumber)
      ? float(Temp->GridData->ReturnChemWorkTime() +
	      Temp->GridData->ReturnGravWorkTime())
      : 0.0;

  double t0 = ReturnWallTime();
  CommunicationAllSumValues(ChildWork, nchild);
  WorkFieldReduceTime += ReturnWallTime() - t0;
  WorkFieldReduceCount++;

  /* The cap: one rank's even share of this level's work, scaled.  If
     nothing was measured - the first rebuild after a restart, say - the
     cap stays zero and the cell cap alone applies, which is the
     historical behaviour. */

  float total = 0.0;
  for (i = 0; i < nchild; i++)
    total += ChildWork[i];

  if (total <= 0) {
    delete [] ChildWork;
    return SUCCESS;
  }

  MaximumSubgridWork =
    SubgridMaximumWorkFraction * total / float(NumberOfProcessors);

  /* Paint each local parent.  A child's cost is spread evenly over the
     parent cells it covers: the measurement has no structure finer than
     one grid, and the parent's cells are the resolution the cut decision
     is made at, so nothing is gained by pretending otherwise. */

  for (Temp = ParentLevel; Temp; Temp = Temp->NextGridThisLevel) {

    grid *P = Temp->GridData;
    if (P->ReturnProcessorNumber() != MyProcessorNumber)
      continue;

    int Rank, Dim[MAX_DIMENSION];
    FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
    P->ReturnGridInfo(&Rank, Dim, Left, Right);

    int size = 1;
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      size *= Dim[dim];

    float *WF = P->AllocateWorkField(size);

    /* cell width from the active region, which is what GridLeftEdge and
       GridRightEdge bracket */
    FLOAT dx[MAX_DIMENSION];
    int Start[MAX_DIMENSION], End[MAX_DIMENSION];
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      Start[dim] = P->GetGridStartIndex(dim);
      End[dim]   = P->GetGridEndIndex(dim);
      int nactive = End[dim] - Start[dim] + 1;
      dx[dim] = (nactive > 0) ? (Right[dim] - Left[dim]) / FLOAT(nactive) : 1.0;
      if (dx[dim] <= 0)
	dx[dim] = 1.0;
    }

    i = 0;
    for (LevelHierarchyEntry *C = OldChildLevel; C;
	 C = C->NextGridThisLevel, i++) {

      if (ChildWork[i] <= 0)
	continue;

      int cRank, cDim[MAX_DIMENSION];
      FLOAT cLeft[MAX_DIMENSION], cRight[MAX_DIMENSION];
      C->GridData->ReturnGridInfo(&cRank, cDim, cLeft, cRight);

      /* overlap of the child with this parent, in parent cell indices */
      int lo[MAX_DIMENSION], hi[MAX_DIMENSION];
      int covered = 1;
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
	FLOAT a = max(Left[dim], cLeft[dim]);
	FLOAT b = min(Right[dim], cRight[dim]);
	if (b <= a) {
	  covered = 0;
	  break;
	}
	lo[dim] = Start[dim] + int(floor((a - Left[dim]) / dx[dim]));
	hi[dim] = Start[dim] + int(ceil ((b - Left[dim]) / dx[dim])) - 1;
	lo[dim] = max(lo[dim], Start[dim]);
	hi[dim] = min(hi[dim], End[dim]);
	if (hi[dim] < lo[dim]) {
	  covered = 0;
	  break;
	}
	covered *= (hi[dim] - lo[dim] + 1);
      }
      if (!covered)
	continue;

      float per_cell = ChildWork[i] / float(covered);
      for (int k = lo[2]; k <= hi[2]; k++)
	for (int j = lo[1]; j <= hi[1]; j++) {
	  int base = (k * Dim[1] + j) * Dim[0];
	  for (int ii = lo[0]; ii <= hi[0]; ii++)
	    WF[base + ii] += per_cell;
	}
    }
  }

  delete [] ChildWork;

  return SUCCESS;
}
