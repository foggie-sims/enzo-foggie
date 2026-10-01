/***********************************************************************
/
/  GRID CLASS (Solve For Luminosity from YSOs on each grid)
/
/  written by: Cameron Trapp
/  date:       Sept, 2026
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
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
 

#define TOLERANCE 2.0e-6
#define MAX_ITERATION 20
 
int grid::SolveForLocallyExtinctRadiation(int level)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  //if (GravitatingMassField == NULL)  // if this is not set we have nothing to do.
  //  return SUCCESS;

  LCAPERF_START("grid_SolveForLocallExtinctRadiation");
 
  /* declarations */
 
  int dim, size = 1, i;
  float tol_dim = TOLERANCE * POW(0.1, 3-GridRank);
  //  if (GridRank == 3)
  //    tol_dim = 1.0e-5;
 

 
  /* Compute adot/a at time = t+1/2dt (time-centered). */
  
  /*
  float* r2_kernel = new float[kernel_size];
  int nn=0;
  int kk,jj,ii;
  float r2;


  int kernel_dim = 2*RefineBy+1;
  int kernel_size = 8 * GridDimension[2] * GridDimension[1] * GridDimension[0];
  for (kk = 0; kk < GridDimension[2];kk++){
    for (jj = 0; jj < GridDimension[1]; jj++){
        for (ii = 0; ii < GridDimension[0]; ii++,nn++){
            zpos = float(kk-RefineBy);
            ypos = float(jj-RefineBy);
            xpos = float(ii-RefineBy);
	        r2 = xpos*xpos + ypos*ypos + zpos*zpos;
            r2 = r2*this->CellWidth[0][0]*this->CellWidth[0][0];
	        r2 = max(r2, dilrad2);
            r2_kernel[nn] =  1.0 / (4.0 * PI * r2);
        }
    }
  }
*/


  int k0,j0,i0;
  int k1,j1,i1;
  int k2,j2,i2;
  int n0=0,n1,n2;
  float rsquared;
  float root_rsquared;
  float refine_factor = POW(RefineBy,level);
  float root_cell_width = this->CellWidth[0][0] * refine_factor;
  float root_dilrad2 = 0.1444*root_cell_width*root_cell_width;
  float near_field_term,far_field_term;

  /* Maybe speed up by replacing with r2 kernel? Tough on memory?*/
  for (k0 = 0; k0 < GridDimension[2]; k0++){ 
    for (j0 = 0; j0 < GridDimension[1]; j0++){
      for (i0 = 0; i0 < GridDimension[0]; i0++,n0++){ //Loop over to find young stellar sources
         if (this->kdissH2SourceField[n0]>0){ //YSO Present on grid

           n1=0;
           for (k1=0; k1 < GridDimension[2];k1++){ //loop over 5^3 kernel
            for (j1=0; j1 < GridDimension[1]; j1++){
              for (i1=0; i1 < GridDimension[0]; i1++, n1++){
 
                  rsquared = ((k1-k0)*(k1-k0) + (j1-j0)*(j1-j0) + (i1-i0)*(i1-i0)) * this->CellWidth[0][0] * this->CellWidth[0][0]
                  if rsquared < dilrad2: rsquared = dilrad2
                  near_field_term = this->kdissH2SourceField[n0] / (4*PI*rsquared); //Term from this grid level

                  root_rsquared = floor( ((k1-k0)*(k1-k0) + (j1-j0)*(j1-j0) + (i1-i0)*(i1-i0)) / refine_factor ) * root_cell_width * root_cell_width;
                  if root_rsquared<root_dilrad2: root_rsquared = root_dilrad2; //Avoid double counting from coarser parent
                  
                  far_field_term  = this->kdissH2SourceField[n0] / (4.0 * PI * root_rsquared);

                  this->kdissH2FluxField[n1] += near_field_term - far_field_term;
              }
            }
           }
         }
      }
    }
  }


            /* kernel alternative
           n1=0;
           for (k1 = 0; k1 < GridDimension[2]; k1++){
             for (j1 = 0; j1 < GridDimension[1]; j1++){
               for (i1 = 0; i1 < GridDeminsion[0]; i1++,n1++){ //Loop over convolve
                 k2 = k1 - k0 + GridDimension[2];
                 j2 = j1 - j0 + GridDimension[1];
                 i2 = i1 - i0 + GridDiemnsion[0];
                 n2 = i2 + j2 * 2 * GridDimension[0] + k2 * 4 * GridDimension[1] * GridDimension[0];
                 this->kdissH2FluxField[n1] += this->kdissH2SourceField[n0] * r2_kernel[n2]; //9+ 6* -> ~15



               }
             }
           }
         }
      }
    }
                 */


 
  LCAPERF_STOP("grid_SolveForLocallExtinctRadiation");
  return SUCCESS;
}
 
