/***********************************************************************
/
/  WRITE THE PER-GRID CHEMISTRY WORK MAP  (audit T2.1 diagnostic)
/
/  written by: performance audit
/  date:       August 2026
/
/  PURPOSE: Records, once per root step, the wall time every locally
/           owned grid spent in chemistry/cooling, together with that
/           grid's level and spatial extent.
/
/           This is measurement only - nothing here feeds the load
/           balancer.  It exists to test the premise that measured-work
/           balancing would rest on: that the chemistry cost of a region
/           of space is predictable from one root step to the next.
/           Grid objects do not survive RebuildHierarchy, so a work
/           estimate can only be carried across a rebuild spatially;
/           this map is what makes that testable before the mapping
/           machinery is built.
/
/           Enabled with GridWorkMapOutput = 1.  Each rank appends to
/           its own file to avoid serialising on a shared one.
/
/  RETURNS: SUCCESS or FAIL
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
#include "Hierarchy.h"
#include "LevelHierarchy.h"

int WriteGridWorkMap(LevelHierarchyEntry *LevelArray[], int cycle,
		     int minlevel)
{

  if (GridWorkMapOutput == 0)
    return SUCCESS;

  /* Records must be written before RebuildHierarchy destroys the grids
     they describe, otherwise the accumulated time dies with the object
     and the deepest levels - which rebuild most often - are sampled
     least.  Every call is collective, so this counter stays in step
     across ranks and orders the intervals globally. */

  static int sequence = 0;
  sequence++;

  char name[MAX_LINE_LENGTH];
  sprintf(name, "gridwork_rank%4.4"ISYM".txt", MyProcessorNumber);

  FILE *fptr = fopen(name, "a");
  if (fptr == NULL) {
    ENZO_VFAIL("Error opening grid work map file %s\n", name)
  }

  /* One header per file, written when it is still empty. */

  if (ftell(fptr) == 0)
    fprintf(fptr, "# seq cycle level rank cells chemwork_sec subcycles "
	    "xl yl zl xr yr zr\n");

  for (int level = minlevel; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (LevelHierarchyEntry *Temp = LevelArray[level]; Temp != NULL;
	 Temp = Temp->NextGridThisLevel) {

      if (Temp->GridData->ReturnProcessorNumber() != MyProcessorNumber)
	continue;

      int Rank, Dims[MAX_DIMENSION];
      FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
      Temp->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);

      fprintf(fptr,
	      "%"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %.10e %"ISYM" "
	      "%.10e %.10e %.10e %.10e %.10e %.10e\n",
	      sequence, cycle, level, MyProcessorNumber,
	      Temp->GridData->GetActiveSize(),
	      Temp->GridData->ReturnChemWorkTime(),
	      Temp->GridData->ReturnChemWorkCalls(),
	      (double) Left[0],  (double) Left[1],  (double) Left[2],
	      (double) Right[0], (double) Right[1], (double) Right[2]);

      /* The accumulator restarts each root step, so consecutive records
	 describe disjoint intervals. */

      Temp->GridData->ClearChemWorkTime();
    }

  fclose(fptr);

  return SUCCESS;
}
