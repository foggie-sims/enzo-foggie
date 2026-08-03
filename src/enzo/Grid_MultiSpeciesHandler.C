/***********************************************************************
/
/  GRID CLASS (HANDLE CALLING AND SOLVING COOLING/CHEMISTRY)
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: Move logic for chemistry/cooling module selection here
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

double ReturnWallTime(void);

int grid::MultiSpeciesHandler()
{
  if ((!MultiSpecies) && (!RadiativeCooling)) return SUCCESS;
  if (GadgetEquilibriumCooling != 0) return SUCCESS;

  LCAPERF_START("grid_MultiSpeciesHandler");

  /* Measure this grid's chemistry/cooling cost for the T2.1 work map.
     Only the owning rank does real work here (the solvers return
     immediately otherwise), so only it accumulates. */

  double chem_t0 = ReturnWallTime();

#ifdef USE_GRACKLE
  if (grackle_data->use_grackle == TRUE) {
    grackle_data->radiative_transfer_intermediate_step = FALSE;
    if (this->GrackleWrapper() == FAIL) {
      ENZO_FAIL("Error in GrackleWrapper.\n");
    }
    if (ProcessorNumber == MyProcessorNumber)
      this->AddChemWorkTime(ReturnWallTime() - chem_t0);
    return SUCCESS;
  }
#endif

  if (MultiSpecies && RadiativeCooling ) {
    int RTCoupledSolverIntermediateStep = FALSE;
    this->SolveRateAndCoolEquations(RTCoupledSolverIntermediateStep);
  } else {
    if (MultiSpecies)
      this->SolveRateEquations();
    if (RadiativeCooling)
      this->SolveRadiativeCooling();
  }

  if (ProblemType == 62)
    this->CoolingTestResetEnergies();

  if (ProcessorNumber == MyProcessorNumber)
    this->AddChemWorkTime(ReturnWallTime() - chem_t0);

  LCAPERF_STOP("grid_MultiSpeciesHandler");
  return SUCCESS;
}
