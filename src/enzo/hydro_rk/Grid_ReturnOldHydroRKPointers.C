/***********************************************************************
/
/  GRID CLASS (RETURNS AN ARRAY OF POINTERS THAT ARE COMPATIBLE WITH 
/              THE HYDRO_RK SOLVERS)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
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
#include "TopGridData.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

int grid::ReturnOldHydroRKPointers(float **Prim, bool ReturnMassFractions)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int i, n, dim, size, nfield, n0;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum, CRNum;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  /* Add the physical quantities */

  if (HydroMethod == HD_RK) {
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum);
    nfield = n0 = NEQ_HYDRO;
  }

  else if (HydroMethod == MHD_RK) {
    if (CRModel)
      this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum, CRNum);
    else
      this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num,
				       PhiNum);
    nfield = n0 = NEQ_MHD;
  }
  
  Prim[iden] = OldBaryonField[DensNum];
  Prim[ivx] = OldBaryonField[Vel1Num];
  if (GridRank > 1)
    Prim[ivy] = OldBaryonField[Vel2Num];
  if (GridRank > 2)
    Prim[ivz] = OldBaryonField[Vel3Num];
  Prim[ietot] = OldBaryonField[TENum];
  if (DualEnergyFormalism)
    Prim[ieint] = OldBaryonField[GENum];

  if (HydroMethod == MHD_RK) {
    Prim[iBx] = OldBaryonField[B1Num];
    Prim[iBy] = OldBaryonField[B2Num];
    Prim[iBz] = OldBaryonField[B3Num];
    Prim[iPhi] = OldBaryonField[PhiNum];
  }
 
  if (CRModel){
    Prim[iCR] = OldBaryonField[CRNum];
  }
  /* Add the species */

  if (MultiSpecies) {

    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
				HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

    //Prim[nfield++] = OldBaryonField[DeNum];
    Prim[nfield++] = OldBaryonField[HINum];
    Prim[nfield++] = OldBaryonField[HIINum];
    Prim[nfield++] = OldBaryonField[HeINum];
    Prim[nfield++] = OldBaryonField[HeIINum];
    Prim[nfield++] = OldBaryonField[HeIIINum];

    if (MultiSpecies > 1) {
      Prim[nfield++] = OldBaryonField[HMNum];
      Prim[nfield++] = OldBaryonField[H2INum];
      Prim[nfield++] = OldBaryonField[H2IINum];
    }

    if (MultiSpecies > 2) {
      Prim[nfield++] = OldBaryonField[DINum];
      Prim[nfield++] = OldBaryonField[DIINum];
      Prim[nfield++] = OldBaryonField[HDINum];
    }

  } // ENDIF MultiSpecies

  /* Add the colours (treat them as species) */

  int SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum, MetalAGBNum, MetalNSMNum; 

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, 
         MetalAGBNum, MetalNSMNum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    return FAIL;
  }
  
  if (MetalNum != -1) {
    Prim[nfield++] = OldBaryonField[MetalNum];
    if (StarMakerTypeIaSNe)
      Prim[nfield++] = OldBaryonField[MetalIaNum];
    else if (StarFeedbackTrackMetalSources) {// mutually exclusive with StarMakerTypeIaSNe
      Prim[nfield++] = OldBaryonField[MetalIaNum];
      Prim[nfield++] = OldBaryonField[MetalIINum];
      Prim[nfield++] = OldBaryonField[MetalAGBNum];
      Prim[nfield++] = OldBaryonField[MetalNSMNum];
    }
    if (StarMakerTypeIISNeMetalField)
      Prim[nfield++] = OldBaryonField[MetalIINum];
    if (MultiMetals || TestProblemData.MultiMetals) {
      Prim[nfield++] = OldBaryonField[MetalNum+1]; // ExtraType0
      Prim[nfield++] = OldBaryonField[MetalNum+2]; // ExtraType1
    }
  }

  if (SNColourNum      != -1) Prim[nfield++] = OldBaryonField[SNColourNum];  
  /*   //##### These fields are currently not being used and only causing interpolation problems
  if (MBHColourNum     != -1) Prim[nfield++] = OldBaryonField[MBHColourNum];
  if (Galaxy1ColourNum != -1) Prim[nfield++] = OldBaryonField[Galaxy1ColourNum];
  if (Galaxy2ColourNum != -1) Prim[nfield++] = OldBaryonField[Galaxy2ColourNum];
  */

  /* Tracer fluid fields */
  int TF01Num, TF02Num, TF03Num, TF04Num, TF05Num, TF06Num, TF07Num, TF08Num;

  if (this->IdentifyTracerFluidFields(TF01Num, TF02Num, TF03Num, TF04Num, TF05Num, TF06Num, TF07Num, TF08Num) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyTracerFluidFields.\n");

  if(TF01Num != -1) Prim[nfield++] = OldBaryonField[TF01Num];
  if(TF02Num != -1) Prim[nfield++] = OldBaryonField[TF02Num];
  if(TF03Num != -1) Prim[nfield++] = OldBaryonField[TF03Num];
  if(TF04Num != -1) Prim[nfield++] = OldBaryonField[TF04Num];
  if(TF05Num != -1) Prim[nfield++] = OldBaryonField[TF05Num];
  if(TF06Num != -1) Prim[nfield++] = OldBaryonField[TF06Num];
  if(TF07Num != -1) Prim[nfield++] = OldBaryonField[TF07Num];
  if(TF08Num != -1) Prim[nfield++] = OldBaryonField[TF08Num];


  /* Convert the species and color fields into mass fractions */

  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (ReturnMassFractions)  
    for (n = n0; n < nfield; n++)
      for (i = 0; i < size; i++)
	Prim[n][i] /= Prim[iden][i];  

  return SUCCESS;

}
