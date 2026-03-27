/***********************************************************************
/
/  INITIALIZE A MULTI-STAR PARTICLE TEST
/
/  written by: Cassi Lochhaas
/  date:       2026
/  modified1:
/
/  PURPOSE:
/    Initialize a multi-star particle test in a uniform medium.
/    Particles are placed in a roughly spherical formation using the
/    Fibonacci spiral algorithm, with user-controlled count and spacing.
/
/    Parameters:
/      TestMultiStarParticleNParticles    -- number of star particles (default 50)
/      TestMultiStarParticleSpacing       -- approx spacing between adjacent
/                                           particles in code units (default 0.1);
/                                           sets sphere radius as
/                                           R = spacing * cbrt(3*N / (4*pi))
/      TestMultiStarParticleCenter        -- center of the sphere (default 0.5 0.5 0.5)
/      TestMultiStarParticleStarMass      -- mass per particle in solar masses (default 100)
/      TestMultiStarParticleStarVelocity  -- velocity for all particles in km/s (default 0 0 0)
/      TestMultiStarParticleStarMetallicityFraction -- metallicity fraction (default 0)
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <string.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int TestMultiStarParticleInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                                    TopGridData &MetaData, float *Initialdt)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *MetalName = "Metal_Density";
  char *MetalIIName = "MetalSNII_Density";
  char *MetalIaName = "MetalSNIa_Density";
  char *MetalAGBName = "MetalAGB_Density";
  char *MetalNSMName = "MetalNSM_Density";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";

  char *TracerFluidO1Name = "TracerFluid01";
  char *TracerFluidO2Name = "TracerFluid02";
  char *TracerFluidO3Name = "TracerFluid03";
  char *TracerFluidO4Name = "TracerFluid04";
  char *TracerFluidO5Name = "TracerFluid05";
  char *TracerFluidO6Name = "TracerFluid06";
  char *TracerFluidO7Name = "TracerFluid07";
  char *TracerFluidO8Name = "TracerFluid08";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int dim, ret;

  /* set default parameters */

  float TestMultiStarParticleDensity     = 1.0;
  float TestMultiStarParticleEnergy      = -1.0;
  float TestMultiStarParticleTemperature = -1.0;
  float TestMultiStarParticleVelocity[3] = {0.0, 0.0, 0.0};
  float TestMultiStarParticleBField[3]   = {0.0, 0.0, 0.0};

  int   TestMultiStarParticleNParticles   = 50;
  FLOAT TestMultiStarParticleSphereRadius = 0.1;
  FLOAT TestMultiStarParticleCenter[3]    = {0.5, 0.5, 0.5};
  float TestMultiStarParticleStarMass     = 100.0;
  FLOAT TestMultiStarParticleStarVelocity[3] = {0.0, 0.0, 0.0};
  float TestMultiStarParticleStarMetallicityFraction = 0.0;

  int TestProblemUseMetallicityField = 1;
  float TestProblemInitialMetallicityFraction = 2e-3; // 0.1 Zsun

  TestProblemData.MultiSpecies = MultiSpecies;
  TestProblemData.UseMetallicityField = TestProblemUseMetallicityField;
  TestProblemData.MetallicityField_Fraction = TestProblemInitialMetallicityFraction;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TestMultiStarParticleDensity = %"FSYM,
                  &TestMultiStarParticleDensity);
    ret += sscanf(line, "TestMultiStarParticleEnergy = %"FSYM,
                  &TestMultiStarParticleEnergy);
    ret += sscanf(line, "TestMultiStarParticleTemperature = %"FSYM,
                  &TestMultiStarParticleTemperature);
    ret += sscanf(line, "TestMultiStarParticleNParticles = %"ISYM,
                  &TestMultiStarParticleNParticles);
    ret += sscanf(line, "TestMultiStarParticleSphereRadius = %"PSYM,
                  &TestMultiStarParticleSphereRadius);
    ret += sscanf(line, "TestMultiStarParticleCenter = %"PSYM" %"PSYM" %"PSYM,
                  &TestMultiStarParticleCenter[0],
                  &TestMultiStarParticleCenter[1],
                  &TestMultiStarParticleCenter[2]);
    ret += sscanf(line, "TestMultiStarParticleStarMass = %"FSYM,
                  &TestMultiStarParticleStarMass);
    ret += sscanf(line, "TestMultiStarParticleStarVelocity = %"PSYM" %"PSYM" %"PSYM,
                  &TestMultiStarParticleStarVelocity[0],
                  &TestMultiStarParticleStarVelocity[1],
                  &TestMultiStarParticleStarVelocity[2]);
    ret += sscanf(line, "TestMultiStarParticleStarMetallicityFraction = %"FSYM,
                  &TestMultiStarParticleStarMetallicityFraction);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    ret += sscanf(line, "TestProblemInitialHIFraction  = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "TestProblemInitialHIIFraction  = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIFraction  = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIFraction  = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "TestMultiStarParticle")
        && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* use either the internal energy or the temperature parameters */

  if ((TestMultiStarParticleEnergy > 0) && (TestMultiStarParticleTemperature > 0))
    ENZO_FAIL("Error in TestMultiStarParticleInitialize: please specify only one of TestMultiStarParticleEnergy and TestMultiStarParticleTemperature");

  if ((TestMultiStarParticleEnergy < 0) && (TestMultiStarParticleTemperature < 0))
    ENZO_FAIL("Error in TestMultiStarParticleInitialize: please specify either TestMultiStarParticleEnergy or TestMultiStarParticleTemperature");

  if (TestMultiStarParticleTemperature > 0) {
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
                 &VelocityUnits, MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in GetUnits.\n");
      return FAIL;
    }
    TestMultiStarParticleEnergy = TestMultiStarParticleTemperature / TemperatureUnits / ((Gamma - 1.0) * 0.6);
  }

  if (TestMultiStarParticleNParticles < 1)
    ENZO_FAIL("Error in TestMultiStarParticleInitialize: TestMultiStarParticleNParticles must be >= 1");

  if (TestMultiStarParticleSphereRadius <= 0.0)
    ENZO_FAIL("Error in TestMultiStarParticleInitialize: TestMultiStarParticleSphereRadius must be > 0");

  /* add gas velocity to internal energy to find total energy */

  float TotalEnergy = TestMultiStarParticleEnergy;
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    TotalEnergy += 0.5 * POW(TestMultiStarParticleVelocity[dim], 2);

  /* set up uniform grid */

  if (TopGrid.GridData->InitializeUniformGrid(TestMultiStarParticleDensity,
                                              TotalEnergy,
                                              TestMultiStarParticleEnergy,
                                              TestMultiStarParticleVelocity,
                                              TestMultiStarParticleBField) == FAIL)
    ENZO_FAIL("Error in InitializeUniformGrid.");

  if (TopGrid.GridData->
      TestMultiStarParticleInitializeGrid(TestMultiStarParticleNParticles,
                                          TestMultiStarParticleSphereRadius,
                                          TestMultiStarParticleCenter,
                                          TestMultiStarParticleStarMass,
                                          TestMultiStarParticleStarVelocity,
                                          TestMultiStarParticleStarMetallicityFraction,
                                          Initialdt) == FAIL)
    ENZO_FAIL("Error in TestMultiStarParticleInitializeGrid.\n");

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if (TestProblemData.MultiSpecies) {
    DataLabel[count++] = ElectronName;
    DataLabel[count++] = HIName;
    DataLabel[count++] = HIIName;
    DataLabel[count++] = HeIName;
    DataLabel[count++] = HeIIName;
    DataLabel[count++] = HeIIIName;
    if (TestProblemData.MultiSpecies > 1) {
      DataLabel[count++] = HMName;
      DataLabel[count++] = H2IName;
      DataLabel[count++] = H2IIName;
    }
    if (TestProblemData.MultiSpecies > 2) {
      DataLabel[count++] = DIName;
      DataLabel[count++] = DIIName;
      DataLabel[count++] = HDIName;
    }
  }
  if (TestProblemData.UseMetallicityField)
    DataLabel[count++] = MetalName;
    if (StarMakerTypeIaSNe || StarFeedbackTrackMetalSources)
      DataLabel[count++] = MetalIaName;
    if (StarFeedbackTrackMetalSources) {
      DataLabel[count++] = MetalIIName;
      DataLabel[count++] = MetalAGBName;
      DataLabel[count++] = MetalNSMName;
    }

  if (UseTracerFluid) {
    if (NumberOfTracerFluidFields >= 1) DataLabel[count++] = TracerFluidO1Name;
    if (NumberOfTracerFluidFields >= 2) DataLabel[count++] = TracerFluidO2Name;
    if (NumberOfTracerFluidFields >= 3) DataLabel[count++] = TracerFluidO3Name;
    if (NumberOfTracerFluidFields >= 4) DataLabel[count++] = TracerFluidO4Name;
    if (NumberOfTracerFluidFields >= 5) DataLabel[count++] = TracerFluidO5Name;
    if (NumberOfTracerFluidFields >= 6) DataLabel[count++] = TracerFluidO6Name;
    if (NumberOfTracerFluidFields >= 7) DataLabel[count++] = TracerFluidO7Name;
    if (NumberOfTracerFluidFields == 8) DataLabel[count++] = TracerFluidO8Name;
  }

  int j;
  for (j = 0; j < count; j++)
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "TestMultiStarParticleDensity = %"FSYM"\n",
            TestMultiStarParticleDensity);
    fprintf(Outfptr, "TestMultiStarParticleEnergy = %"FSYM"\n",
            TestMultiStarParticleEnergy);
    fprintf(Outfptr, "TestMultiStarParticleNParticles = %"ISYM"\n",
            TestMultiStarParticleNParticles);
    fprintf(Outfptr, "TestMultiStarParticleSphereRadius = %"PSYM"\n",
            TestMultiStarParticleSphereRadius);
    fprintf(Outfptr, "MetallicityField_Fraction = %"FSYM"\n",
            TestProblemData.MetallicityField_Fraction);
  }

  fprintf(stderr, "TestMultiStarParticleDensity = %"FSYM"\n",
          TestMultiStarParticleDensity);
  fprintf(stderr, "TestMultiStarParticleEnergy = %"FSYM"\n",
          TestMultiStarParticleEnergy);
  fprintf(stderr, "TestMultiStarParticleNParticles = %"ISYM"\n",
          TestMultiStarParticleNParticles);
  fprintf(stderr, "TestMultiStarParticleSphereRadius = %"PSYM"\n",
          TestMultiStarParticleSphereRadius);
  fprintf(stderr, "MetallicityField_Fraction = %"FSYM"\n",
          TestProblemData.MetallicityField_Fraction);

  return SUCCESS;
}
