/***********************************************************************
/
/  ACTIVE PARTICLE HELPER ROUTINE:
/   Return Number of Active Particles with ID = ActiveParticleIDToFind
/
/  written by: Nathan Goldbaum
/  date:       March 2012
/
*************************************************************************/

#ifdef USE_MPI
#endif /* USE_MPI */

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"
#include "SortCompareFunctions.h"


int grid::ReturnNumberOfParticlesOfThisType(int ParticleIDToFind) {

  // Return if this does not concern us
  if (MyProcessorNumber != ProcessorNumber)
    return 0;
  
  int NumberOfParticlesOfThisType = 0;
  for (int j = 0; j<NumberOfParticles; j++) {
    if (this->ParticleType[j] == ParticleIDToFind) {
      NumberOfParticlesOfThisType++;
    }
  }
  return NumberOfParticlesOfThisType;
}
