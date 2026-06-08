/***********************************************************************
/
/  GRID CLASS (WRAP THE GRACKLE CHEMISTRY SOLVER)
/
/  written by: Britton Smith
/  date:       April, 2013
/  modified1:
/
/  PURPOSE: Solve chemistry and cooling with grackle.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
int search_lower_bound(float *arr, float value, int low, int high, 
		       int total);
           
int grid::GrackleWrapper()
{

#ifdef USE_GRACKLE
  if (grackle_data->use_grackle == FALSE)
    return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  LCAPERF_START("grid_GrackleWrapper");

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  double dt_cool = dtFixed;
#ifdef TRANSFER
  dt_cool = (grackle_data->radiative_transfer_intermediate_step == TRUE) ? dtPhoton : dtFixed;
#endif

  
  /* Compute the size of the fields. */
 
  int i;
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
  g_grid_dimension = new Eint32[GridRank];
  g_grid_start = new Eint32[GridRank];
  g_grid_end = new Eint32[GridRank];
  for (i = 0; i < GridRank; i++) {
    g_grid_dimension[i] = (Eint32) GridDimension[i];
    g_grid_start[i] = (Eint32) GridStartIndex[i];
    g_grid_end[i] = (Eint32) GridEndIndex[i];
  }
 
  /* Find fields: density, total energy, velocity1-3. */
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* Find Multi-species fields. */

  DeNum = HINum = HIINum = HeINum = HeIINum = HeIIINum = HMNum = H2INum = 
    H2IINum = DINum = DIINum = HDINum = 0;
 
  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }
 
  /* Get easy to handle pointers for each variable. */
 
  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  float *volumetric_heating_rate = NULL;
  float *specific_heating_rate   = NULL;

  /* Compute the cooling time. */

  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dt_cool, &a, &dadt);
 
    aUnits = 1.0/(1.0 + InitialRedshift);
  } else if (RadiationFieldRedshift > -1){
    a        = 1.0 / (1.0 + RadiationFieldRedshift);
    aUnits   = 1.0;
  }
  float afloat = float(a);

  /* Update units. */

  code_units grackle_units;
  grackle_units.comoving_coordinates = (Eint32) ComovingCoordinates;
  grackle_units.density_units        = (double) DensityUnits;
  grackle_units.length_units         = (double) LengthUnits;
  grackle_units.time_units           = (double) TimeUnits;
  grackle_units.velocity_units       = (double) VelocityUnits;
  grackle_units.a_units              = (double) aUnits;
  grackle_units.a_value              = (double) a;

  /* Metal cooling codes. */
 
  int MetalNum = 0, SNColourNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  // Double check if there's a metal field when we have metal cooling
  if (MetalCooling && MetalFieldPresent == FALSE) {
    ENZO_FAIL("Metal cooling is on, but no metal field present.");
  }

  /* If both metal fields (Pop I/II and III) exist, create a field
     that contains their sum */

  float *MetalPointer = NULL;
  float *TotalMetals = NULL;

  if (MetalNum != -1 && SNColourNum != -1) {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
      TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalPointer = TotalMetals;
  } // ENDIF both metal types
  else {
    if (MetalNum != -1)
      MetalPointer = BaryonField[MetalNum];
    else if (SNColourNum != -1)
      MetalPointer = BaryonField[SNColourNum];
  } // ENDELSE both metal types

  /*
    Cap the metal density at 90% of the gas density.
    This was historically in Enzo, then Grackle and was removed in
    Grackle PR #121 (https://github.com/grackle-project/grackle/pull/121)
  */
  if (MetalPointer != NULL) {
    for (i = 0; i < size; i++) {
      MetalPointer[i] = min(MetalPointer[i], 0.9 * BaryonField[DensNum][i]);
    }
  }
 
  int temp_thermal = FALSE;
  float *thermal_energy;
  if ( UseMHD ){
    iBx = FindField(Bfield1, FieldType, NumberOfBaryonFields);
    iBy = FindField(Bfield2, FieldType, NumberOfBaryonFields);
    iBz = FindField(Bfield3, FieldType, NumberOfBaryonFields);  
  }

  if (HydroMethod==Zeus_Hydro) {
    thermal_energy = BaryonField[TENum];
  }
  else if (DualEnergyFormalism) {
    thermal_energy = BaryonField[GENum];
  }
  else {
    temp_thermal = TRUE;
    thermal_energy = new float[size];
    for (i = 0; i < size; i++) {
      thermal_energy[i] = BaryonField[TENum][i] - 
        0.5 * POW(BaryonField[Vel1Num][i], 2.0);
      if(GridRank > 1)
        thermal_energy[i] -= 0.5 * POW(BaryonField[Vel2Num][i], 2.0);
      if(GridRank > 2)
        thermal_energy[i] -= 0.5 * POW(BaryonField[Vel3Num][i], 2.0);

      if( UseMHD ) {
        thermal_energy[i] -= 0.5 * (POW(BaryonField[iBx][i], 2.0) + 
                                    POW(BaryonField[iBy][i], 2.0) + 
                                    POW(BaryonField[iBz][i], 2.0)) / 
          BaryonField[DensNum][i];
      }
    } // for (int i = 0; i < size; i++)
  }

  //
  // Put code here to assign fields to volumetric or specific
  // heating rate pointers
  //

  /* set up grackle fields object */
  grackle_field_data my_fields;

  my_fields.grid_rank = (Eint32) GridRank;
  my_fields.grid_dimension = g_grid_dimension;
  my_fields.grid_start     = g_grid_start;
  my_fields.grid_end       = g_grid_end;
  my_fields.grid_dx        = this->CellWidth[0][0];

  /* now add in the baryon fields */
  my_fields.density         = density;
  my_fields.internal_energy = thermal_energy;
  my_fields.x_velocity      = velocity1;
  my_fields.y_velocity      = velocity2;
  my_fields.z_velocity      = velocity3;
  my_fields.HI_density      = BaryonField[HINum];
  my_fields.HII_density     = BaryonField[HIINum];
  my_fields.HeI_density     = BaryonField[HeINum];
  my_fields.HeII_density    = BaryonField[HeIINum];
  my_fields.HeIII_density   = BaryonField[HeIIINum];
  my_fields.e_density       = BaryonField[DeNum];

  my_fields.HM_density      = BaryonField[HMNum];
  my_fields.H2I_density     = BaryonField[H2INum];
  my_fields.H2II_density    = BaryonField[H2IINum];

  my_fields.DI_density      = BaryonField[DINum];
  my_fields.DII_density     = BaryonField[DIINum];
  my_fields.HDI_density     = BaryonField[HDINum];

  my_fields.metal_density   = MetalPointer;

  my_fields.volumetric_heating_rate = volumetric_heating_rate;
  my_fields.specific_heating_rate   = specific_heating_rate;

#ifdef TRANSFER
  /* Find RT fields */
  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum,
        gammaNum, kphHMNum, kdissH2IINum;

  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum,
                                  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  /* unit conversion from Enzo RT units to CGS */
  float rtunits = erg_eV / TimeUnits;

  if( RadiativeTransfer ){
    my_fields.RT_HI_ionization_rate   = BaryonField[kphHINum];

    if (RadiativeTransferHydrogenOnly == FALSE){
      my_fields.RT_HeI_ionization_rate  = BaryonField[kphHeINum];
      my_fields.RT_HeII_ionization_rate = BaryonField[kphHeIINum];
    }

    if (MultiSpecies > 1)
      my_fields.RT_H2_dissociation_rate = BaryonField[kdissH2INum];

    /* need to convert to CGS units */
    for( i = 0; i < size; i++) BaryonField[gammaNum][i] *= rtunits;

    my_fields.RT_heating_rate = BaryonField[gammaNum];
    

  }
#endif // TRANSFER



  /* Estimate local radiation field from new Stars */
  //Simplest Alternative:
  //Sum all the mass of all young star particles on the grid
  //Get an estimate of the LW photon production from fitting Starbust99 results (w/o loading in a table)
  //Convert to RT_H2_dissociation_rate and add to grackle fields

  //Sum all mass of young star particles on grid
  /*
  static const float metallicity_bins[6] = {0.0, 0.0004, 0.002, 0.006, 0.014, 0.02}; //Metallicity bins for SB99 tables
  static const float age_bins[11] = {0,5e6,1e7,1.5e7,2e7,2.5e7,3e7,3.5e7,4e7,4.5e7,5e7}; //Age bins for SB99 tables. To do: check units, currently years

  float years_to_seconds = 3.15576e7f; //Seconds in a year
  float MassUnits = DensityUnits * POW(LengthUnits,3);

  //Cheap implementation of lookup tables for SB99. Should do full interpolation in tabular_feedback, but this will give a quick proof of concept.
  //I've generated the full tables and have the framework for reading them in following the tabular_feedback scheme, but will leave as this for now.
  //H2 photodissociation rates from Lyman-Werner photons in units of cm^2/s per Msun of stars. To convert to RT units, need to multiply by young_star_mass and divide by LengthUnits^2
  //Each row corresponds to a different metallicity bin, each column corresponds to a different age bin.
  static const float kdiss_H2_sb99[5][10] = {{8.616753475461946e+28, 3.0609842295148783e+28, 1.7793259544275774e+28, 1.5579481975570319e+28, 7.336999024747869e+27, 6.16771377021248e+27, 4.934612144598924e+27, 4.672521186828272e+27, 4.4146802267792734e+27, 4.083479510828315e+27},
      {9.721982223039793e+28, 4.304822507060031e+28, 1.8644395112296062e+28, 1.0619152696364634e+28, 7.848321220290816e+27, 5.642701209907671e+27, 4.4481784324738383e+27, 3.5499680262354453e+27, 2.8571716506680415e+27, 2.3976026941788822e+27},
      {1.0573363485793936e+29, 3.774233520528279e+28, 1.4727367324254517e+28, 8.13310248001581e+27, 5.741325336346668e+27, 4.0256979767777485e+27, 3.126377653994177e+27, 2.3896288079559745e+27, 1.8599668891094716e+27, 1.612425486113449e+27},
      {1.0891678985339411e+29, 3.2456721710048436e+28, 1.1576683096773793e+28, 6.311547803801093e+27, 4.382626900327114e+27, 2.9285407981848826e+27, 2.2622327737699204e+27, 1.6952595798785602e+27, 1.2641920083191978e+27, 1.0863416172338696e+27},
      {1.0801135857923557e+29, 2.8455249259730975e+28, 9.102459204349159e+27, 4.668747480690447e+27, 2.889872536333432e+27, 2.0071649538074096e+27, 1.3899321088430588e+27, 9.873410563441971e+26, 7.167294859311428e+26, 5.888193311217845e+26}};

  //H- photodetachment rates in units of cm^2/s per Msun of stars. To convert to RT units, need to multiply by young_star_mass and divide by LengthUnits^2
  static const float kdet_HM_sb99[5][10] = {{3.8241226738957066e+30, 1.2691075221881906e+30, 6.320187615663541e+29, 4.589734608406248e+29, 2.0886091048841194e+29, 1.7686502309462676e+29, 1.449086479665172e+29, 1.3742786972477078e+29, 1.3156335472850541e+29, 1.248175564829444e+29},
      {5.355080973335116e+30, 1.506101136754258e+30, 7.367051989107666e+29, 4.94826722391895e+29, 4.3672645585508546e+29, 3.646219930630708e+29, 3.2993246094530356e+29, 3.349240499947722e+29, 2.81035243418751e+29, 2.3917918242342782e+29},
      {4.975427798636371e+30, 2.11430381110695e+30, 1.0134575716416472e+30, 7.015201614220793e+29, 6.536613088761112e+29, 5.45449846075009e+29, 4.6504715263353755e+29, 4.8418816203992866e+29, 4.325479957866995e+29, 3.2547242824127256e+29},
      {4.5062467799254436e+30, 2.7441414260306114e+30, 1.5135046492296408e+30, 9.762287625447108e+29, 7.833203390009536e+29, 6.784162146506757e+29, 5.1530282236800166e+29, 5.319482714958312e+29, 5.128378151096622e+29, 3.811304647568074e+29},
      {4.367128430849667e+30, 2.7899493363546306e+30, 1.438774324771867e+30, 8.63880338912953e+29, 6.441637068851158e+29, 5.4212026702962396e+29, 5.100041216387224e+29, 4.6203283682706974e+29, 4.27653791399001e+29, 3.507518956008909e+29}};

  float k_diss_H2 = 0; //Photodissociation rate for H2 from LW band
  float k_det_HM = 0; //Photodetachment rate for H- from photons above 0.755 eV
  for (i = 0; i < this->NumberOfParticles; i++) {
    if (this->ParticleType[i] == PARTICLE_TYPE_STAR) {
      float age = (this->Time - this->ParticleAttribute[0][i]) * TimeUnits / years_to_seconds; //Convert to yr
      if (age < 5e7) { //To do: Double check units are being handled correctly throughout
          //fprintf(stdout, "CWT: Interpolating for star particle with age %"FSYM" yr?\n",age);

          //To do: This stuff shouldn't be hardcoded like this, but I'm assuming we will replace this all with actually reading the tables
          int aa = search_lower_bound((float*)age_bins, age, 0, 11, 11); 
          float t_age=0.5f;
          if (aa>=10){
            aa=9;
            t_age = 1;
          }
          else if (aa<0){
            aa=0;
            t_age=0;
          }
          else{
              t_age = (age - age_bins[aa]) / (age_bins[aa+1] - age_bins[aa]);
          }

          float metallicity = this->ParticleAttribute[2][i];
          int zz = search_lower_bound((float*)metallicity_bins, metallicity, 0, 6, 6);

          float t_z=0.5f;
          if (zz>=5){
            zz=4;
            t_z = 1;
          }
          else if (zz<0){
            zz=0;
            t_z=0;
          }
          else{
              t_z = (metallicity - metallicity_bins[zz]) / (metallicity_bins[zz+1] - metallicity_bins[zz]);
          }

          float k_diss_H2_sb99_interp = (1-t_age) * (1-t_z) * kdiss_H2_sb99[zz][aa] + t_age * (1-t_z) * kdiss_H2_sb99[zz][aa+1] + (1-t_age) * t_z * kdiss_H2_sb99[zz+1][aa] + t_age * t_z * kdiss_H2_sb99[zz+1][aa+1];
          float k_det_HM_sb99_interp = (1-t_age) * (1-t_z) * kdet_HM_sb99[zz][aa] + t_age * (1-t_z) * kdet_HM_sb99[zz][aa+1] + (1-t_age) * t_z * kdet_HM_sb99[zz+1][aa] + t_age * t_z * kdet_HM_sb99[zz+1][aa+1];

          float dx = this->CellWidth[0][0];
          float ParticleMass_Msun = this->ParticleMass[i] * dx * dx * dx * MassUnits / (1.989e33); //Convert from code mass to Msun
          //fprintf(stdout, "CWT: Interpolating for star particle with mass %"FSYM" Msun?\n",ParticleMass_Msun);

          k_diss_H2 += k_diss_H2_sb99_interp * ParticleMass_Msun; //Convert from code mass to Msun
          k_det_HM += k_det_HM_sb99_interp * ParticleMass_Msun; //Convert from code mass to Msun

//ParticleMassCode = StarMass * SolarMass * POW(LengthUnits * CellWidth[0][0], -3.0) / DensityUnits;


      }
    }
  }
  k_diss_H2 = k_diss_H2 * TimeUnits / (LengthUnits * LengthUnits); //Convert from cm^2/s to code units //Correct?
  k_det_HM  = k_det_HM  * TimeUnits / (LengthUnits * LengthUnits); //Convert from cm^2/s to code units
  float grid_dx = this->GridRightEdge[0]-this->GridLeftEdge[0];
  float grid_dy = this->GridRightEdge[1]-this->GridLeftEdge[1];
  float grid_dz = this->GridRightEdge[2]-this->GridLeftEdge[2];
  //This is the most tunable part of this code, as calculating the r^2 for each cell will get expensive
  //Currently estimating as half the average extent of the grid which isn't great
  //To do: Account for any local extinction from unresolved sources around stars?
  float dilutionRadius = 0.5 * (grid_dx + grid_dy + grid_dz)/3.0; //Get Half the average extent of the grid
  //float dilutionRadius = 4.848e-6 * pc_cm / (double) LengthUnits;  // 1 AU //Try an extreme case
  float dilRad2 = dilutionRadius * dilutionRadius;
  k_diss_H2 = k_diss_H2  / (4.0 * 3.14159 * dilRad2);
  k_det_HM = k_det_HM    / (4.0 * 3.14159 * dilRad2);
  
  */
  float *k_diss_H2_grid  = new float[size];
  float *k_det_HM_grid  = new float[size];
  for (int i = 0; i < size; i++){
      k_diss_H2_grid[i] = k_diss_H2I_grid_sum; //From Grid Property written in StarParticleHandler
      k_det_HM_grid[i] = k_det_HM_grid_sum; 
  }
  fprintf(stdout, "CWT: Setting My Fields...\n");
  fprintf(stdout, "CWT: k_diss_H2 = %"FSYM" Hz?\n", k_diss_H2_grid[0] / TimeUnits);
  fprintf(stdout, "CWT: k_det_HM  = %"FSYM" Hz?\n", k_det_HM_grid[0] / TimeUnits);

  my_fields.RT_H2_dissociation_rate =  k_diss_H2_grid;//Already in units of seconds (from table)
  my_fields.RT_HM_detachment_rate =  k_det_HM_grid;  //Feeds in Britton's Grackle Branch (foggie-sf) only

  // Need to set the other fields to the same 0 array for now
  float *EmptyRtArray0  = new float[size];
  float *EmptyRtArray1  = new float[size];
  float *EmptyRtArray2  = new float[size];
  float *EmptyRtArray3  = new float[size];

  for (int i = 0; i < size; i++){
      EmptyRtArray0[i] = 0;
      EmptyRtArray1[i] = 0;
      EmptyRtArray2[i] = 0;
      EmptyRtArray3[i] = 0;
  }

  float *EmptyRTArray  = new float[size];

  for (int i = 0; i < size; i++){
      EmptyRTArray[i] = 0;
  }

  my_fields.RT_HI_ionization_rate   = EmptyRtArray0;
  my_fields.RT_HeI_ionization_rate  = EmptyRtArray1;
  my_fields.RT_HeII_ionization_rate = EmptyRtArray2;
  my_fields.RT_heating_rate = EmptyRtArray3;
  
  /*                                              */

  /* Call the chemistry solver. */
  if (solve_chemistry(&grackle_units, &my_fields, (double) dt_cool) == FAIL){
    fprintf(stderr, "Error in Grackle solve_chemistry.\n");
    return FAIL;
  }

  if (HydroMethod != Zeus_Hydro) {
    for (i = 0; i < size; i++) {
      BaryonField[TENum][i] = thermal_energy[i] +
        0.5 * POW(BaryonField[Vel1Num][i], 2.0);
      if(GridRank > 1)
        BaryonField[TENum][i] += 0.5 * POW(BaryonField[Vel2Num][i], 2.0);
      if(GridRank > 2)
        BaryonField[TENum][i] += 0.5 * POW(BaryonField[Vel3Num][i], 2.0);

      if( UseMHD ) {
        BaryonField[TENum][i] += 0.5 * (POW(BaryonField[iBx][i], 2.0) + 
                                        POW(BaryonField[iBy][i], 2.0) + 
                                        POW(BaryonField[iBz][i], 2.0)) / 
          BaryonField[DensNum][i];
      }

    } // for (int i = 0; i < size; i++)
  } // if (HydroMethod != Zeus_Hydro)

  if (temp_thermal == TRUE) {
    delete [] thermal_energy;
  }

#ifdef TRANSFER
  if (RadiativeTransfer){
    /* convert the RT units back to Enzo */
    for(i = 0; i < size; i ++) BaryonField[gammaNum][i] /= rtunits;

  }
#endif TRANSFER


  delete [] TotalMetals;
  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

  delete[] k_diss_H2_grid;
  delete[] k_det_HM_grid;
  delete[] EmptyRtArray0;
  delete[] EmptyRtArray1;
  delete[] EmptyRtArray2;
  delete[] EmptyRtArray3;


  LCAPERF_STOP("grid_GrackleWrapper");

#endif
  return SUCCESS;
}
