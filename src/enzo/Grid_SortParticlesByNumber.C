/***********************************************************************
/
/  GRID CLASS (SORT PARTICLES BY PARTICLE NUMBER)
/
/  written by: Greg Bryan
/  date:       Jan, 2001
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
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
void QuickSortAndDrag(PINT List[], int left, int right,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, int   *DragList3[]);
 
 
void grid::SortParticlesByNumber()
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || NumberOfParticles == 0)
    return;
 
  int dim, j;
  int masses = (StarMakerStoreInitialMass) ? 2 : 1;
 
  /* Allocate arrays of pointer, one for float type and one for FLOAT type,
     and file them up with pointers to the particle data. */
 
  float **DragList1 = new float*[GridRank+masses+NumberOfParticleAttributes];
  FLOAT **DragList2 = new FLOAT*[GridRank];
  int   **DragList3 = new int*[1];
  for (dim = 0; dim < GridRank; dim++) {
    DragList2[dim] = ParticlePosition[dim];
    DragList1[dim] = ParticleVelocity[dim];
  }
  DragList1[GridRank] = ParticleMass;
  if (StarMakerStoreInitialMass) DragList1[GridRank+1] = ParticleInitialMass;
  DragList3[0]        = ParticleType;
  for (j = 0; j < NumberOfParticleAttributes; j++)
    DragList1[GridRank+masses+j] = ParticleAttribute[j];
 
  /* Sort by particle index, dragging the data along. */
 
  QuickSortAndDrag(ParticleNumber, 0, NumberOfParticles-1,
		   GridRank+masses+NumberOfParticleAttributes, DragList1,
		   GridRank, DragList2, 1, DragList3);
 
  /* Clean up. */
 
  delete [] DragList1;
  delete [] DragList2;
  delete [] DragList3;
 
  return;
}
