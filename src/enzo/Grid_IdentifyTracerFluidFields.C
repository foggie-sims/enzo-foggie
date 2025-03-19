/***********************************************************************
/
/  GRID CLASS (IDENTIFY TRACER FLUID FIELDS)
/
/  written by: Brian O'Shea
/  date:       February 2025
/  modified1:
/
/  PURPOSE:  Identifies field numbers for tracer fluid fields
/
/  NOTE:  If you add additional tracer fluid fields, you have to modify this
/         AND you have to update MAX_NUMBER_OF_TRACER_FIELDS in macros_and_parameters.h
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);

int grid::IdentifyTracerFluidFields(int &TracerFluidField01Num, int &TracerFluidField02Num,
			            int &TracerFluidField03Num, int &TracerFluidField04Num,
			            int &TracerFluidField05Num, int &TracerFluidField06Num,
			            int &TracerFluidField07Num, int &TracerFluidField08Num)
{

  TracerFluidField01Num = TracerFluidField02Num = TracerFluidField03Num = TracerFluidField04Num =
    TracerFluidField05Num = TracerFluidField06Num = TracerFluidField07Num=TracerFluidField08Num = 0;
			       
 
  TracerFluidField01Num = FindField(TracerFluidField01Density, FieldType, NumberOfBaryonFields);
  TracerFluidField02Num = FindField(TracerFluidField02Density, FieldType, NumberOfBaryonFields);
  TracerFluidField03Num = FindField(TracerFluidField03Density, FieldType, NumberOfBaryonFields);
  TracerFluidField04Num = FindField(TracerFluidField04Density, FieldType, NumberOfBaryonFields);
  TracerFluidField05Num = FindField(TracerFluidField05Density, FieldType, NumberOfBaryonFields);
  TracerFluidField06Num = FindField(TracerFluidField06Density, FieldType, NumberOfBaryonFields);
  TracerFluidField07Num = FindField(TracerFluidField07Density, FieldType, NumberOfBaryonFields);
  TracerFluidField08Num = FindField(TracerFluidField08Density, FieldType, NumberOfBaryonFields);

  return SUCCESS;
}

