#include <stdio.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

int TestRefinementSchemeInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData)
{

    char *DensName = "Density";
    char *TEName   = "TotalEnergy";
    char *GEName   = "GasEnergy";
    char *Vel1Name = "x-velocity";
    char *Vel2Name = "y-velocity";
    char *Vel3Name = "z-velocity";

    int i = 0;
    DataLabel[i++] = DensName;
    DataLabel[i++] = TEName;
    DataLabel[i++] = GEName;
    DataLabel[i++] = Vel1Name;
    DataLabel[i++] = Vel2Name;
    DataLabel[i++] = Vel3Name;

    float TestRefinementUniformDensity = 1.0;
    float TestRefinementUniformTotalEnergy = 1.0;
    float TestRefinementUniformInternalEnergy = 1.0;
    float TestRefinementUniformVelocity[MAX_DIMENSION];
    float TestRefinementUniformBField[MAX_DIMENSION];
    float TestRefinementUniformCR = 0;

    /* Set velocity and B-fields to 0 */
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        TestRefinementUniformBField[dim] = 0.0;
        TestRefinementUniformVelocity[dim] = 0.0;   
    }

    if (TopGrid.GridData->InitializeUniformGrid(TestRefinementUniformDensity,
        TestRefinementUniformTotalEnergy, TestRefinementUniformInternalEnergy,
        TestRefinementUniformVelocity, TestRefinementUniformBField, 
        TestRefinementUniformCR) == FAIL) 
    {
        ENZO_FAIL("Error in InitializeUniformGrid.");
    }

    return SUCCESS;
}