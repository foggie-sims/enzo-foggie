#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"

int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData, 
		 ExternalBoundary *Exterior, FLOAT WriteTime);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int HydroShockTubesInitialize(FILE *fptr, FILE *Outfptr, 
			      HierarchyEntry &TopGrid, TopGridData &MetaData) 
{
  char *DensName = "Density";
  char *PresName = "Pressure";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ColourName = "colour";
  char *MachName   = "Mach";
  char *PSTempName = "PreShock_Temperature";
  char *PSDenName  = "PreShock_Density";

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
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = FALSE;
  float  InitialDiscontinuity = 0.5, SecondDiscontinuity = FLOAT_UNDEFINED,
    LeftDensity = 1.0, RightDensity = 1.0, CenterDensity = 1.0, 
    LeftVelocityX = 0.0, RightVelocityX = 0.0, CenterVelocityX = 0.0,
    LeftVelocityY = 0.0, RightVelocityY = 0.0, CenterVelocityY = 0.0,
    LeftVelocityZ = 0.0, RightVelocityZ = 0.0, CenterVelocityZ = 0.0,
    LeftPressure = 1.0, RightPressure = 1.0, CenterPressure = 1.0;
  
  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "HydroShockTubesRefineAtStart = %"ISYM, 
		  &RefineAtStart);
    ret += sscanf(line, "HydroShockTubesInitialDiscontinuity = %"FSYM, 
		  &InitialDiscontinuity);
    ret += sscanf(line, "HydroShockTubesSecondDiscontinuity = %"FSYM, 
		  &SecondDiscontinuity);
    ret += sscanf(line, "HydroShockTubesLeftVelocityX = %"FSYM,
		  &LeftVelocityX);
    ret += sscanf(line, "HydroShockTubesLeftVelocityY = %"FSYM,
		  &LeftVelocityY);
    ret += sscanf(line, "HydroShockTubesLeftVelocityZ = %"FSYM,
		  &LeftVelocityZ);
    ret += sscanf(line, "HydroShockTubesLeftPressure = %"FSYM, 
		  &LeftPressure);
    ret += sscanf(line, "HydroShockTubesLeftDensity = %"FSYM, 
		  &LeftDensity);
    ret += sscanf(line, "HydroShockTubesRightVelocityX = %"FSYM, 
		  &RightVelocityX);
    ret += sscanf(line, "HydroShockTubesRightVelocityY = %"FSYM, 
		  &RightVelocityY);
    ret += sscanf(line, "HydroShockTubesRightVelocityZ = %"FSYM, 
		  &RightVelocityZ);
    ret += sscanf(line, "HydroShockTubesRightPressure = %"FSYM, 
		  &RightPressure);
    ret += sscanf(line, "HydroShockTubesRightDensity = %"FSYM,
                  &RightDensity);
    ret += sscanf(line, "HydroShockTubesCenterVelocityX = %"FSYM, 
		  &CenterVelocityX);
    ret += sscanf(line, "HydroShockTubesCenterVelocityY = %"FSYM, 
		  &CenterVelocityY);
    ret += sscanf(line, "HydroShockTubesCenterVelocityZ = %"FSYM, 
		  &CenterVelocityZ);
    ret += sscanf(line, "HydroShockTubesCenterPressure = %"FSYM, 
		  &CenterPressure);
    ret += sscanf(line, "HydroShockTubesCenterDensity = %"FSYM,
                  &CenterDensity);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "HydroShockTubes") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file
  
  float DensityUnits = 1, LengthUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  /* check if SecondDiscontinuity is defined, if not, then make it
     sit on top of InitialDiscontinuity (to mimic a 2-region shock) */

  if (SecondDiscontinuity == FLOAT_UNDEFINED) {
      SecondDiscontinuity = InitialDiscontinuity;
  }

  /* set up grid */

  TopGrid.GridData->
    HydroShockTubesInitializeGrid(InitialDiscontinuity, 
				  SecondDiscontinuity,
				  LeftDensity, RightDensity, CenterDensity,
				  LeftVelocityX,  RightVelocityX, CenterVelocityX,
				  LeftVelocityY,  RightVelocityY, CenterVelocityY,
				  LeftVelocityZ,  RightVelocityZ, CenterVelocityZ,
				  LeftPressure,   RightPressure,  CenterPressure);

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */

  if (RefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      printf("In level %"ISYM"\n", level);
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	Temp->GridData->
	  HydroShockTubesInitializeGrid
	  (InitialDiscontinuity, SecondDiscontinuity,
	   LeftDensity, RightDensity, CenterDensity,
	   LeftVelocityX,  RightVelocityX, CenterVelocityX,
	   LeftVelocityY,  RightVelocityY, CenterVelocityY,
	   LeftVelocityZ,  RightVelocityZ, CenterVelocityZ,
	   LeftPressure,   RightPressure,  CenterPressure);
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels


    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
				   *LevelArray[level-1]->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

    //WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
    //       &TopGrid, MetaData, Exterior, -1);

  } // end: if (RefineAtStart)


  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = Vel1Name;
  if (TopGrid.GridData->GetGridRank() > 1 || HydroMethod > 2) {
    DataLabel[count++] = Vel2Name;
    if (TopGrid.GridData->GetGridRank() > 2 || HydroMethod > 2) {
      DataLabel[count++] = Vel3Name;
    }
  }
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = GEName;
  }

  if (ShockMethod) {
    DataLabel[count++] = MachName;
    if(StorePreShockFields){
      DataLabel[count++] = PSTempName;
      DataLabel[count++] = PSDenName;
    }
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

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "HydroShockTubesRefineAtStart        = %"ISYM"\n",
	    RefineAtStart);
    fprintf(Outfptr, "HydroShockTubesInitialDiscontinuity = %"FSYM"\n",
	    InitialDiscontinuity);
    fprintf(Outfptr, "HydroShockTubesLeftDensity          = %"FSYM"\n",
	    LeftDensity);
    fprintf(Outfptr, "HydroShockTubesRightDensity         = %"FSYM"\n",
	    RightDensity);
    fprintf(Outfptr, "HydroShockTubesLeftVelocityX        = %"FSYM"\n",
	    LeftVelocityX);
    fprintf(Outfptr, "HydroShockTubesRightVelocityX       = %"FSYM"\n",
            RightVelocityX);
    fprintf(Outfptr, "HydroShockTubesLeftVelocityY        = %"FSYM"\n",
	    LeftVelocityY);
    fprintf(Outfptr, "HydroShockTubesRightVelocityY       = %"FSYM"\n",
            RightVelocityY);
    fprintf(Outfptr, "HydroShockTubesLeftVelocityZ        = %"FSYM"\n",
	    LeftVelocityZ);
    fprintf(Outfptr, "HydroShockTubesRightVelocityZ       = %"FSYM"\n",
            RightVelocityZ);
    fprintf(Outfptr, "HydroShockTubesLeftPressure         = %"FSYM"\n",
            LeftPressure);
    fprintf(Outfptr, "HydroShockTubesRightPressure        = %"FSYM"\n",
            RightPressure);

    fprintf(Outfptr, "HydroShockTubesSecondDiscontinuity = %"FSYM"\n",
	    SecondDiscontinuity);
    fprintf(Outfptr, "HydroShockTubesCenterDensity       = %"FSYM"\n",
	    CenterDensity);
    fprintf(Outfptr, "HydroShockTubesCenterVelocityX     = %"FSYM"\n",
	    CenterVelocityX);
    fprintf(Outfptr, "HydroShockTubesCenterVelocityY     = %"FSYM"\n",
	    CenterVelocityY);
    fprintf(Outfptr, "HydroShockTubesCenterVelocityZ     = %"FSYM"\n",
	    CenterVelocityZ);
    fprintf(Outfptr, "HydroShockTubesCenterPressure      = %"FSYM"\n",
	    CenterPressure);

  }
  //return FAIL;
  return SUCCESS;

}

