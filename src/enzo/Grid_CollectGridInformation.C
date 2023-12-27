/***********************************************************************
/
/  GRID CLASS (COLLECT SUMMARY INFORMATION)
/
/  written by: Greg Bryan
/  date:       September, 1996
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int grid::CollectGridInformation(int &GridMemory, float &GridVolume,
				 int &CellsActive, float &AxialRatio,
				 int &CellsTotal, int &Particles, int level)
{
 
  GridMemory    = NumberOfBaryonFields*sizeof(float);
  GridVolume    = 1;
  CellsTotal    = 1;
  CellsActive   = 1;
  int MaxDim    = 1, MinDim = GridDimension[0];
 
  for (int dim = 0; dim < GridRank; dim++) {
    int DimActive = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    GridMemory    *= GridDimension[dim];
    GridVolume    *= (GridRightEdge[dim]   - GridLeftEdge[dim]  )/
                     (DomainRightEdge[dim] - DomainLeftEdge[dim]);
    CellsActive *= DimActive;
    CellsTotal *= GridDimension[dim];
    MaxDim = max(MaxDim, DimActive);
    MinDim = min(MinDim, DimActive);
  }
 
  if (level > 3) 
    GridMemory = GridMemory + NumberOfParticles*
               (sizeof(float)*(GridRank+2+NumberOfParticleAttributes) + sizeof(float)*GridRank);
 
  AxialRatio = float(MaxDim)/float(MinDim);

  printf("JT Grid_CGI0: The input grid level is: %d\n", level); 
  printf("JT Grid_CGI1: ID=%d GridLevel=%d  #Cells %d CellsTot %d GridMem %d Part %d \n", ID, level, CellsActive, CellsTotal, GridMemory, NumberOfParticles);
  printf("JT Grid_CGI2: ID=%d   LeftEdge: %f %f %f\n", ID, GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
  printf("JT Grid_CGI3: ID=%d  RightEdge: %f %f %f\n", ID, GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);

  return SUCCESS;
 
}
