/***
* CHECK FOR LIBYT CALL
* PURPOSE:
*   This routine checks if libyt in situ analysis workflow should be called.
*/

#ifdef USE_LIBYT

#include "preincludes.h"

#include <stdio.h>
#include <stdlib.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#include "CommunicationUtilities.h"

int CallInSitulibyt(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
                    int level, int from_topgrid);

int CheckForLibytCall(LevelHierarchyEntry *LevelArray[], TopGridData &MetaData) {
     /* Check for libyt call: cycle-based. */
  if (MetaData.CycleNumber >= CycleLastLibytCall + CycleSkipLibytCall   &&
      CycleSkipLibytCall > 0) {

    CycleLastLibytCall += CycleSkipLibytCall;

    // Call libyt in situ routine
    CallInSitulibyt(LevelArray, &MetaData, 0, 1);
  }
    return SUCCESS;
}

#endif // #ifdef USE_LIBYT