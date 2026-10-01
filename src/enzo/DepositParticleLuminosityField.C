/***********************************************************************
/
/  DEPOSIT PARTICLES INTO LUMINOSITY FIELDS IN THIS GRID AND ALL SUBGRIDS
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:Cameron Trapp
/  date:       Sept, 2026
/
/  PURPOSE: LEBRON-like implementation of RT
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
#include "TopGridData.h"
#include "LevelHierarchy.h"
 
/* function prototypes */
 
int DepositParticleLuminosityFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT Time);
 
int DepositParticleLuminosityField(HierarchyEntry *Grid, FLOAT TimeMidStep)
{
 
  /* Get the time and dt for this grid.  Compute time+1/2 dt. */
 
  if (TimeMidStep < 0)
    TimeMidStep =     Grid->GridData->ReturnTime() +
                  0.5*Grid->GridData->ReturnTimeStep();
 
  /* Initialize the gravitating mass field only if in send-receive mode
     (i.e. this routine is called only once) or if in the first of the
     three communication modes (post-receive). */


  //I think this part can be piggy backed of the previous DepositParticleMass call
  //if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
  //    CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
 
    /* Initialize the gravitating mass field parameters (if necessary). */
 
  //  if (Grid->GridData->InitializeGravitatingMassFieldParticles(RefineBy)
  //                                                                == FAIL) {
  //    ENZO_FAIL("Error in grid->InitializeGravitatingMassFieldParticles.\n");
  //  }
 
    /* Clear the GravitatingMassFieldParticles. */
 
  //  if (Grid->GridData->ClearGravitatingMassFieldParticles() == FAIL) {
  //    ENZO_FAIL("Error in grid->ClearGravitatingMassFieldParticles.\n");
   // }
 
//  fprintf(stderr, "--DepositParticleMassField (Send) Initialize & Clear\n");
 
  //} // end: if (CommunicationDirection != COMMUNICATION_SEND)
  
  /* Deposit particles to GravitatingMassFieldParticles in this grid. */
 
//  fprintf(stderr, "--DepositParticleMassField Call DepositParticlePositions\n");
 
  if (Grid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
				 KDISSH2_SOURCE_FIELD) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }

  if (Grid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
				 KDETHM_SOURCE_FIELD) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }
 
  if (Grid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
				 ISRF_SOURCE_FIELD) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }

  /* Recursively deposit particles in children (at TimeMidStep). */
 
  if (Grid->NextGridNextLevel != NULL)
    if (DepositParticleLuminosityFieldChildren(Grid, Grid->NextGridNextLevel,
					 TimeMidStep)
	== FAIL) {
      ENZO_FAIL("Error in DepositParticleLuminosityFieldChildren.\n");
    }
 
  return SUCCESS;
}
 
 
 
 
int DepositParticleLuminosityFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT DepositTime)
{
 
  /* Deposit particles in Grid into DepositGrid at the given time. */
 
  if (Grid->GridData->DepositParticlePositions(DepositGrid->GridData,
		     DepositTime, KDISSH2_SOURCE_FIELD) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }
 
  if (Grid->GridData->DepositParticlePositions(DepositGrid->GridData,
		     DepositTime, KDETHM_SOURCE_FIELD) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }

  if (Grid->GridData->DepositParticlePositions(DepositGrid->GridData,
		     DepositTime, ISRF_SOURCE_FIELD) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }

  /* Next grid on this level. */
 
  if (Grid->NextGridThisLevel != NULL)
    if (DepositParticleLuminosityFieldChildren(DepositGrid, Grid->NextGridThisLevel,
					 DepositTime) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassFieldChildren(1).\n");
    }
 
  /* Recursively deposit particles in children. */
 
  if (Grid->NextGridNextLevel != NULL)
    if (DepositParticleLuminosityFieldChildren(DepositGrid, Grid->NextGridNextLevel,
					 DepositTime) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassFieldChildren(2).\n");

    }
 
 
  return SUCCESS;
}
