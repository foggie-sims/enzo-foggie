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

int grid::TestRefinementSchemeInitializeGrid()
{

    float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];

    /*  Return if this grid is not on this processor. */
    if (ProcessorNumber != MyProcessorNumber)
        return SUCCESS;

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                     Vel3Num, TENum, B1Num, B2Num, B3Num);

                            
    //And then set the values for the fields in the manner that suits your
    //problem.
    int index=0;
    for(int k=0; k<GridDimension[2];k++){
        for(int j=0; j<GridDimension[1];j++){
            for(int i=0; i<GridDimension[0]; i++){
               index = i+GridDimension[0]*(j+GridDimension[1]*k);
               BaryonField[DensNum][index] = 46.2; //define your problem here.
            }
        }
    }

}