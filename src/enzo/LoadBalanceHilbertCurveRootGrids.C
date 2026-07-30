/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE BY A HILBERT CURVE
/
/  written by: John Wise
/  date:       April, 2010
/  modified1:
/
/  NOTES: For a given level, sort the grids on a 3D Hilbert curve, and
/         then partition the list with equal amounts of work.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "communication.h"
#include "CommunicationUtilities.h"
#include "SortCompareFunctions.h"

hilbert_key HilbertCurve3D(FLOAT *coord);
double HilbertKeyDiff(const hilbert_key &a, const hilbert_key &b);
Eint32 compare_hkey(const void *a, const void *b);

#define FUZZY_BOUNDARY 0.1
#define FUZZY_ITERATIONS 10

int LoadBalanceHilbertCurveRootGrids(FLOAT *GridCenters[], int *CellCount,
				     int NumberOfGrids, int* &RootProcessors)
{

  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;

  /* Initialize */
  
  RootProcessors = new int[NumberOfGrids];
  int *GridWork = new int[NumberOfGrids];
  hilbert_data *HilbertData = new hilbert_data[NumberOfGrids];
  int *BlockDivisions = new int[NumberOfProcessors];
  // 64-bit work accounting: the summed ghost-inclusive cell counts
  // overflow 32 bits above ~2.1e9 cells, silently garbling partitions.
  long long *ProcessorWork = new long long[NumberOfProcessors];

  FLOAT ThisCenter[MAX_DIMENSION];
  long long TotalWork, WorkThisProcessor, WorkPerProcessor, WorkLeft;
  int i, dim, grid_num, Rank, block_num, Dims[MAX_DIMENSION];
  int iter;

  /* Compute the position of each grid on a Hilbert curve */
  // TODO: PARALLELIZE

  for (i = 0; i < NumberOfGrids; i++) {

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ThisCenter[dim] = GridCenters[dim][i];

    HilbertData[i].grid_num = i;
    HilbertData[i].hkey = HilbertCurve3D(ThisCenter);
    
  } // ENDFOR grids

  /* Sort the grids along the curve and partition it into pieces with
     equal amounts of work. */

  //qsort(HilbertData, NumberOfGrids, sizeof(hilbert_data), compare_hkey);
  std::sort(HilbertData, HilbertData+NumberOfGrids, cmp_hkey());
  TotalWork = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    GridWork[i] = CellCount[HilbertData[i].grid_num];
    TotalWork += CellCount[HilbertData[i].grid_num];
  }

  /* Partition into nearly equal workloads */

  grid_num = 0;
  WorkLeft = TotalWork;
  for (i = 0; i < NumberOfProcessors-1; i++) {
    WorkThisProcessor = 0;
    WorkPerProcessor = WorkLeft / (NumberOfProcessors-i);
    while (WorkThisProcessor < WorkPerProcessor &&
	   grid_num < NumberOfGrids) {
      WorkThisProcessor += GridWork[grid_num];
      grid_num++;
    } // ENDWHILE

    // Determine if removing the last grid results in a closer match
    if (grid_num > 1)
      if (ABS(WorkThisProcessor - GridWork[grid_num-1] - WorkPerProcessor) <
	  ABS(WorkThisProcessor - WorkPerProcessor)) {
	WorkThisProcessor -= GridWork[grid_num-1];
	grid_num--;
      }

    // -1 because we advanced grid_num before checking the workload
    BlockDivisions[i] = grid_num-1;  
    ProcessorWork[i] = WorkThisProcessor;
    WorkLeft -= WorkThisProcessor;

  } // ENDFOR processors

  /* Fill in the last entry */

  BlockDivisions[i] = NumberOfGrids-1;
  ProcessorWork[i] = WorkLeft;

  /* Mark the new processor numbers with the above divisions. */

  block_num = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    if (i > BlockDivisions[block_num]) 
      block_num++;
    RootProcessors[HilbertData[i].grid_num] = block_num;
  }

  /* Fuzzy boundaries: There will be some work imbalance if we split
     the grids along the Hilbert curve because the grids where the
     work is split are too large.  So we further refine the balancing
     by moving grids to adjacent processors within some epsilon on the
     Hilbert curve. */

  const double CriticalBalance = 0.1 / NumberOfProcessors;
  hilbert_key div_key;
  double lo_span, hi_span, boundary_off;
  char direction;
  int LoadedBlock, UnloadedBlock;
  long long WorkDifference;
  long long MinWork, MaxWork;
  float WorkImbalance;

  for (iter = 0; iter < FUZZY_ITERATIONS; iter++) {
    MinWork = 0x7FFFFFFFFFFFFFFFLL;
    MaxWork = -1;
    for (i = 0; i < NumberOfProcessors-1; i++) {

      /* Degenerate partitions (a block with zero grids leaves its
	 division at -1, or 0 with nothing before it) have no boundary
	 to fuzz - and indexing HilbertData with them reads out of
	 bounds. */
      if (BlockDivisions[i] <= 0 || BlockDivisions[i+1] < 0 ||
	  (i > 0 && BlockDivisions[i-1] < 0)) {
	MinWork = min(MinWork, ProcessorWork[i]);
	MaxWork = max(MaxWork, ProcessorWork[i]);
	continue;
      }

      /* Curve spans (as double offsets from the division key) of the
	 segments we may move grids across.  Keys are exact integers;
	 the fuzzy interpolation is a heuristic, so double precision on
	 the differences is plenty. */
      div_key = HilbertData[BlockDivisions[i]].hkey;
      if (i == 0)
	lo_span = HilbertKeyDiff(div_key, HilbertData[0].hkey);
      else
	lo_span = HilbertKeyDiff(div_key,
				 HilbertData[BlockDivisions[i-1]].hkey);
      hi_span = HilbertKeyDiff(HilbertData[BlockDivisions[i+1]].hkey,
			       div_key);

      // Which processor has more work?
      if (ProcessorWork[i] > ProcessorWork[i+1]) {
	LoadedBlock = i;
	UnloadedBlock = i+1;
	boundary_off = -FUZZY_BOUNDARY * lo_span;
	direction = -1;
      } else {
	LoadedBlock = i+1;
	UnloadedBlock = i;
	boundary_off = FUZZY_BOUNDARY * hi_span;
	direction = +1;
      }

      /* Move grids from the loaded to unloaded processor until we
	 reach the Hilbert key boundary */

      grid_num = BlockDivisions[i] + direction;
      while (grid_num >= 0 && grid_num < NumberOfGrids &&
	     direction * (boundary_off -
	       HilbertKeyDiff(HilbertData[grid_num].hkey, div_key)) > 0) {
	WorkDifference = 
	  ProcessorWork[LoadedBlock] - ProcessorWork[UnloadedBlock];
	if (2*GridWork[grid_num] < WorkDifference) {
//	  if (debug)
//	    printf("Moving grid %d (work=%d) from P%d -> P%d\n",
//		   grid_num, GridWork[grid_num], LoadedBlock, UnloadedBlock);
	  ProcessorWork[LoadedBlock] -= GridWork[grid_num];
	  ProcessorWork[UnloadedBlock] += GridWork[grid_num];
	  RootProcessors[HilbertData[grid_num].grid_num] = UnloadedBlock;
	}
	grid_num += direction;
      } // ENDWHILE move grids
      MinWork = min(MinWork, ProcessorWork[i]);
      MaxWork = max(MaxWork, ProcessorWork[i]);
    } // ENDFOR processors

    MinWork = min(MinWork, ProcessorWork[NumberOfProcessors-1]);
    MaxWork = max(MaxWork, ProcessorWork[NumberOfProcessors-1]);
    // A block can legitimately hold zero work; never divide by it.
    if (MinWork <= 0)
      continue;
    WorkImbalance = float(MaxWork - MinWork) / float(MinWork);
    if (WorkImbalance < CriticalBalance)
      break;

  } // ENDFOR iterations

  /* Cleanup */
  
  delete [] GridWork;
  delete [] HilbertData;
  delete [] ProcessorWork;
  delete [] BlockDivisions;

  return SUCCESS;

}
