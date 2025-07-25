/***********************************************************************
/
/  INITIALIZE A STAR PARTICLE TEST
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/    Initialize a star particle test in a uniform medium.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

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

int TestStarParticleInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			       TopGridData &MetaData,float *Initialdt)
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

  /* Error check. */


  /* set default parameters */

  float TestStarParticleDensity     = 1.0;
  float TestStarParticleEnergy      = -1.0;
  float TestStarParticleTemperature = -1.0;
  float TestStarParticleVelocity[3] = {0.0, 0.0, 0.0};
  FLOAT TestStarParticleStarVelocity[3] = {0.0, 0.0, 0.0};
  FLOAT TestStarParticleStarPosition[3] = {0.5, 0.5, 0.5};
  float TestStarParticleBField[3]   = {0.0, 0.0, 0.0};
  float TestStarParticleStarMass    = 100.0;
  int TestProblemUseMetallicityField = 1;
  float TestProblemInitialMetallicityFraction = 2e-3; // 0.1 Zsun
  float TestStarParticleStarMetallicityFraction = 0.0;

  TestProblemData.MultiSpecies = MultiSpecies;
  TestProblemData.UseMetallicityField = TestProblemUseMetallicityField;
  TestProblemData.MetallicityField_Fraction = TestProblemInitialMetallicityFraction;

  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TestStarParticleDensity = %"FSYM,
		  &TestStarParticleDensity);
    ret += sscanf(line, "TestStarParticleEnergy = %"FSYM,
		  &TestStarParticleEnergy);
    ret += sscanf(line, "TestStarParticleTemperature = %"FSYM,
		  &TestStarParticleTemperature);
    ret += sscanf(line, "TestStarParticleStarMass = %"FSYM,
		  &TestStarParticleStarMass);
    ret += sscanf(line, "TestStarParticleStarMetallicityFraction = %"FSYM,
      &TestStarParticleStarMetallicityFraction);
    ret += sscanf(line,"TestStarParticleStarVelocity = %"PSYM" %"PSYM" %"PSYM, 
		  &TestStarParticleStarVelocity[0],
		  &TestStarParticleStarVelocity[1],
		  &TestStarParticleStarVelocity[2]);
    ret += sscanf(line,"TestStarParticleStarPosition = %"PSYM" %"PSYM" %"PSYM, 
		  &TestStarParticleStarPosition[0],
		  &TestStarParticleStarPosition[1],
		  &TestStarParticleStarPosition[2]);
    

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction); 

    ret += sscanf(line, "TestProblemInitialHIFraction  = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "TestProblemInitialHIIFraction  = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIFraction  = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIFraction  = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);


    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "TestStarParticle") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* use either the internal energy or the temperature parameters */

  if ((TestStarParticleEnergy > 0) && (TestStarParticleTemperature > 0))
    ENZO_FAIL("Error in TestStarParticleInitialize: please specify only one of TestStarParticleEnergy and TestStarParticleTemperature");

  if ((TestStarParticleEnergy < 0) && (TestStarParticleTemperature < 0))
    ENZO_FAIL("Error in TestStarParticleInitialize: please specify either TestStarParticleEnergy or TestStarParticleTemperature");

  if (TestStarParticleTemperature > 0){
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
    if (GetUnits(&DensityUnits, &LengthUnits,&TemperatureUnits, &TimeUnits,
              &VelocityUnits, MetaData.Time) == FAIL){
      fprintf(stderr, "Error in GetUnits.\n");
      return FAIL;
    }

    TestStarParticleEnergy = TestStarParticleTemperature/TemperatureUnits/((Gamma-1.0)*0.6);
  }

  /* add gas velocity to internal energy to find total energy */

  float TotalEnergy = TestStarParticleEnergy;
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    TotalEnergy += 0.5*POW(TestStarParticleVelocity[dim], 2);
  

  /* set up uniform grid as of before explosion */

  if (TopGrid.GridData->InitializeUniformGrid(TestStarParticleDensity, 
					      TotalEnergy,
					      TestStarParticleEnergy,
					      TestStarParticleVelocity,
					      TestStarParticleBField) == FAIL)
    ENZO_FAIL("Error in InitializeUniformGrid.");

  if (TopGrid.GridData->
      TestStarParticleInitializeGrid(TestStarParticleStarMass,
				     Initialdt, 
				     TestStarParticleStarVelocity,
				     TestStarParticleStarPosition,
             TestStarParticleStarMetallicityFraction) == FAIL)
  ENZO_FAIL("Error in TestStarParticleInitializeGrid.\n");

  /* set up field names and units */
  
  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if (TestProblemData.MultiSpecies){
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

  if(UseTracerFluid){
    if(NumberOfTracerFluidFields >= 1) DataLabel[count++] = TracerFluidO1Name;
    if(NumberOfTracerFluidFields >= 2) DataLabel[count++] = TracerFluidO2Name;
    if(NumberOfTracerFluidFields >= 3) DataLabel[count++] = TracerFluidO3Name;
    if(NumberOfTracerFluidFields >= 4) DataLabel[count++] = TracerFluidO4Name;
    if(NumberOfTracerFluidFields >= 5) DataLabel[count++] = TracerFluidO5Name;
    if(NumberOfTracerFluidFields >= 6) DataLabel[count++] = TracerFluidO6Name;
    if(NumberOfTracerFluidFields >= 7) DataLabel[count++] = TracerFluidO7Name;
    if(NumberOfTracerFluidFields == 8) DataLabel[count++] = TracerFluidO8Name;
  }


  int j;
  for(j=0; j < count; j++)
    DataUnits[j] = NULL;

   
  /* Write parameters to parameter output file */
  
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "TestStarParticleDensity = %"FSYM"\n",
	    TestStarParticleDensity);
    fprintf(Outfptr, "TestStarParticleEnergy = %"FSYM"\n",
	    TestStarParticleEnergy);
    fprintf(Outfptr, "MetallicityField_Fraction = %"FSYM"\n",
            TestProblemData.MetallicityField_Fraction);
  }

  fprintf(stderr, "TestStarParticleDensity = %"FSYM"\n",
	  TestStarParticleDensity);
  fprintf(stderr, "TestStarParticleEnergy = %"FSYM"\n",
	  TestStarParticleEnergy);
  fprintf(stderr, "MetallicityField_Fraction = %"FSYM"\n",
	  TestProblemData.MetallicityField_Fraction);



  return SUCCESS;

}

