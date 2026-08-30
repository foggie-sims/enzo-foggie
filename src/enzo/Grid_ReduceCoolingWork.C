/***********************************************************************
/
/  GRID CLASS (REDUCE A COOLING-TIME FIELD TO LOAD-BALANCE WORK PROXIES)
/
/  written by: Claude, 2026-08-30
/
/  PURPOSE:  Turn a per-cell cooling-time array into a few cheap scalars
/     that predict how long the Grackle solve will take on this grid, so
/     the load balancer can weight grids by real work instead of by cell
/     count (CommunicationLoadBalanceGrids.C:84 currently uses
/     ComputeTime[i] = float(NumberOfCells), which is blind to gas state).
/
/  *** THIS IS A DIAGNOSTIC.  IT MUST NOT ENTER ANY PHYSICS PATH. ***
/
/     It reads BaryonField and the supplied cooling_time array and writes
/     ONLY to this grid's CoolingWork* members, which are consumed solely
/     by the load balancer and by the calibration dump.  The quantities
/     here are deliberately approximate and are meaningless as physics:
/     do not feed them to hydro, chemistry, feedback, or the timestep.
/
/  WHY THESE PROXIES:  Grackle sub-cycles a 1D slice until its SLOWEST
/     cell converges -- solve_rate_cool_g.F:817-819 leaves the iteration
/     loop only when ttmin, the least-advanced cell in the slice, reaches
/     dt.  Converged cells are masked (itmask, :813) but the slice is
/     still swept.  So a row's cost has a term set by its single shortest
/     cooling time (its densest cell) and a term set by how many of its
/     cells are actually still integrating:
/
/       cost(row) ~ nx * n_iter_max(row) * w_overhead     <- CoolingWorkRowMax
/                 + sum_cells n_iter(i)  * w_work         <- CoolingWorkCellSum
/
/     A mean over cells would capture neither: it is dominated by
/     volume-filling diffuse gas and under-predicts exactly the dense
/     grids that drive the imbalance.  The two terms are recorded
/     separately so their weights can be fitted against the measured
/     per-grid time rather than guessed.
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

/* Grackle's own sub-cycle ceiling; the estimator saturates here for the
   same reason the solver does, so a single runaway cell cannot make the
   proxy diverge. */
#define COOLING_WORK_ITMAX 10000.0

/* dtit is limited to ~0.1 of the local cooling time (solve_rate_cool_g.F
   :594), so a cell needing to cross dt takes roughly 10*dt/t_cool
   sub-iterations. */
#define COOLING_WORK_ITER_PER_TCOOL 10.0

int grid::ReduceCoolingWork(float *cooling_time)
{

  /* Only the owning rank has the data. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  if (cooling_time == NULL || GridRank < 1)
    return SUCCESS;

  this->CoolingWorkRowMax   = 0.0;
  this->CoolingWorkCellSum  = 0.0;
  this->CoolingWorkDenseCells = 0;

  float dt = this->dtFixed;
  if (dt <= 0.0)
    return SUCCESS;

  int nx = GridDimension[0];
  int i, j, k, index;

  /* Active zone only: ghost cells are not integrated by the solver, so
     including them would inflate both proxies by the boundary volume. */

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {

      float row_iter_max = 0.0;
      index = (k*GridDimension[1] + j)*nx + GridStartIndex[0];

      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

        /* t_cool is negative where the gas is being heated; its
           magnitude still sets the sub-cycle length, so use fabs.  A
           zero entry would divide by zero, so skip it. */

        float tc = fabs(cooling_time[index]);
        if (tc <= tiny_number)
          continue;

        float iter = COOLING_WORK_ITER_PER_TCOOL * dt / tc;
        if (iter > COOLING_WORK_ITMAX)
          iter = COOLING_WORK_ITMAX;
        if (iter < 1.0)
          iter = 1.0;

        this->CoolingWorkCellSum += iter;
        if (iter > row_iter_max)
          row_iter_max = iter;
        if (tc < dt)
          this->CoolingWorkDenseCells++;
      }

      /* The whole row is swept until its laggard finishes. */
      this->CoolingWorkRowMax += row_iter_max * float(GridEndIndex[0] -
                                                     GridStartIndex[0] + 1);
    }

  return SUCCESS;
}
