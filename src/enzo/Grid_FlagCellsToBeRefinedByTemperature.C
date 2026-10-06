/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED WHERE COOLING TIME < DX/CSOUND)
/
/  written by: Greg Bryan
/  date:       February, 2000
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
         float *TemperatureUnits, float *TimeUnits,
         float *VelocityUnits, FLOAT Time);
 
 
int grid::FlagCellsToBeRefinedByTemperature()
{
  /* declarations */
 
  int i, dim;
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute the temperature field. */
 
  float *temperature = NULL;
  temperature = new float[size];
  if (this->ComputeTemperatureField(temperature) == FAIL)
    ENZO_FAIL("Error in grid->ComputeTemperature.");

  /* Get units. */

  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  /* Loop over grid, checking if temperature > threshold and mass > threshold */

  float cell_mass, cell_vol = pow(CellWidth[0][0] * LengthUnits, 3);
  for (i = 0; i < size; i++) {
    cell_mass = BaryonField[DensNum][i] * DensityUnits * cell_vol / SolarMass;
    if (temperature[i] > MinimumTemperatureForRefinement && cell_mass > TemperatureRefinementStoppingMassMsun)
      FlaggingField[i]++;
  }

  /* clean up */
 
  delete [] temperature;
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)

      NumberOfFlaggedCells++;
 
  return NumberOfFlaggedCells;
 
}
