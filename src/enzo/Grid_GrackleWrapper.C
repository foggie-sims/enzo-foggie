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

  float *k_diss_H2_grid = NULL;
  float *k_det_HM_grid = NULL;
  float *k_diss_CO_grid = NULL;
  float *k_ion_CI_grid = NULL;
  float *k_ion_OI_grid = NULL;
  float *EmptyRtArray = NULL;
  if (UseLocalStellarRadiation){
    /* Estimate local radiation field from new Stars */
    //Sum all the mass of all young star particles on the grid
    //Get an estimate of the LW photon production from fitting results 
    //Convert to RT_H2_dissociation_rate and add to grackle fields
    float years_to_seconds = 3.15576e7f; //Seconds in a year
    float MassUnits = DensityUnits * POW(LengthUnits,3);
  
    //n_metal_bins = pSNFBTable.n_met
    //n_age_bins = pSNFBTable.n_age
    //metallicity bins = pSNFBTable.ini_met?
    //age bins = pSNFBTable.pop_age
    //kdiss_H2_sb99 = pSNFBTable.kdiss_H2
    //kdet_HM = pSNFBTable.kdet_HM

    k_diss_H2I_grid_sum = 0; //Grid Attribute
    k_det_HM_grid_sum = 0; //Grid Attribute - Calculated here and used in CoolingRate
    k_diss_COI_grid_sum = 0;
    k_ion_CI_grid_sum = 0;
    k_ion_OI_grid_sum = 0;
    for (i = 0; i < this->NumberOfParticles; i++) {
      if (this->ParticleType[i] == PARTICLE_TYPE_STAR) {
        float age = (this->Time - this->ParticleAttribute[0][i]) * TimeUnits / years_to_seconds; //Convert to yr
        if (age < 5e7) { 
          //int aa = search_lower_bound((float*)pSNFBTable.pop_age, age, 0, pSNFBTable.n_age+1, pSNFBTable.n_age+1);  //+1?
          float dt_table = pSNFBTable.pop_age[1] - pSNFBTable.pop_age[0];
          float t_age = (age - pSNFBTable.pop_age[0]) / dt_table;
          int aa = (int)floor(t_age);
          if (aa>=pSNFBTable.n_age){
            aa=pSNFBTable.n_age-1;
            t_age = 1;
          }
          else if (aa<0){
            aa=0;
            t_age=0;
          }
          else{
              t_age = (age - pSNFBTable.pop_age[aa]) / (pSNFBTable.pop_age[aa+1] - pSNFBTable.pop_age[aa]);
          }

          float metallicity = this->ParticleAttribute[2][i];
         // int zz = search_lower_bound((float*)metallicity_bins, metallicity, 0, 6, 6);
          int zz = search_lower_bound((float*)pSNFBTable.ini_met, metallicity, 0, pSNFBTable.n_met+1, pSNFBTable.n_met+1);

          //fprintf(stdout, "Age %"ESYM", ", age);
          //fprintf(stdout, "aa %"ISYM",", aa);
          //fprintf(stdout, "Z %"ESYM", ", metallicity);
          //fprintf(stdout, "zz %"ISYM"\n", zz);

          float t_z=0.5f;
          if (zz>=pSNFBTable.n_met){
            zz=pSNFBTable.n_met-1;
            t_z = 1;
          }
          else if (zz<0){
            zz=0;
            t_z=0;
          }
          else{
              t_z = (metallicity - pSNFBTable.ini_met[zz]) / (pSNFBTable.ini_met[zz+1] - pSNFBTable.ini_met[zz]);
          }

          int ii0 = zz * pSNFBTable.n_age + aa;
          int ii1 = zz * pSNFBTable.n_age + (aa+1);
          int ii2 = (zz+1) * pSNFBTable.n_age + aa;
          int ii3 = (zz+1) * pSNFBTable.n_age + (aa+1);

          /* In Units Hz/cm^2 per Solar Mass*/
          float k_diss_H2_sb99_interp = (1-t_age) * (1-t_z) * pSNFBTable.kdiss_H2[ii0] + t_age * (1-t_z) * pSNFBTable.kdiss_H2[ii1] + (1-t_age) * t_z * pSNFBTable.kdiss_H2[ii2] + t_age * t_z * pSNFBTable.kdiss_H2[ii3];
          float k_det_HM_sb99_interp = (1-t_age) * (1-t_z) * pSNFBTable.kdet_HM[ii0] + t_age * (1-t_z) * pSNFBTable.kdet_HM[ii1] + (1-t_age) * t_z * pSNFBTable.kdet_HM[ii2] + t_age * t_z * pSNFBTable.kdet_HM[ii3];
          float k_diss_CO_sb99_interp = (1-t_age) * (1-t_z) * pSNFBTable.kdiss_CO[ii0] + t_age * (1-t_z) * pSNFBTable.kdiss_CO[ii1] + (1-t_age) * t_z * pSNFBTable.kdiss_CO[ii2] + t_age * t_z * pSNFBTable.kdiss_CO[ii3];
          float k_ion_CI_sb99_interp = (1-t_age) * (1-t_z) * pSNFBTable.kion_CI[ii0] + t_age * (1-t_z) * pSNFBTable.kion_CI[ii1] + (1-t_age) * t_z * pSNFBTable.kion_CI[ii2] + t_age * t_z * pSNFBTable.kion_CI[ii3];
          float k_ion_OI_sb99_interp = (1-t_age) * (1-t_z) * pSNFBTable.kion_OI[ii0] + t_age * (1-t_z) * pSNFBTable.kion_OI[ii1] + (1-t_age) * t_z * pSNFBTable.kion_OI[ii2] + t_age * t_z * pSNFBTable.kion_OI[ii3];


          float dx = this->CellWidth[0][0];
          float ParticleMass_Msun = this->ParticleMass[i] * dx * dx * dx *  MassUnits / SolarMass; //Convert from code mass to Msun

          k_diss_H2I_grid_sum += k_diss_H2_sb99_interp * ParticleMass_Msun;
          k_det_HM_grid_sum += k_det_HM_sb99_interp * ParticleMass_Msun;
          k_diss_COI_grid_sum += k_diss_CO_sb99_interp * ParticleMass_Msun;
          k_ion_CI_grid_sum += k_ion_CI_sb99_interp * ParticleMass_Msun;
          k_ion_OI_grid_sum += k_ion_OI_sb99_interp * ParticleMass_Msun;


        }
      }
    }

    k_diss_H2I_grid_sum = k_diss_H2I_grid_sum * TimeUnits / (LengthUnits * LengthUnits); //Convert from cm^2/s to code units //Correct?
    k_det_HM_grid_sum  = k_det_HM_grid_sum  * TimeUnits / (LengthUnits * LengthUnits); //Convert from cm^2/s to code units
    k_diss_COI_grid_sum = k_diss_COI_grid_sum * TimeUnits / (LengthUnits * LengthUnits); //Convert from cm^2/s to code units //Correct?
    k_ion_CI_grid_sum = k_ion_CI_grid_sum * TimeUnits / (LengthUnits * LengthUnits); //Convert from cm^2/s to code units //Correct?
    k_ion_OI_grid_sum = k_ion_OI_grid_sum * TimeUnits / (LengthUnits * LengthUnits); //Convert from cm^2/s to code units //Correct?

    float grid_dx = this->GridRightEdge[0]-this->GridLeftEdge[0];
    float grid_dy = this->GridRightEdge[1]-this->GridLeftEdge[1];
    float grid_dz = this->GridRightEdge[2]-this->GridLeftEdge[2];
    //This is the most tunable part of this code, as calculating the r^2 for each cell will get expensive
    //Currently estimating as half the average extent of the grid which isn't great
    //To do: Account for any local extinction from unresolved sources around stars?
    float dilutionRadius = 0.5 * (grid_dx + grid_dy + grid_dz)/3.0; //Get Half the average extent of the grid
    //float dilutionRadius = 4.848e-6 * pc_cm / (double) LengthUnits;  // 1 AU //Try an extreme case
    float dilRad2 = dilutionRadius * dilutionRadius;
    k_diss_H2I_grid_sum = k_diss_H2I_grid_sum  / (4.0 * 3.14159 * dilRad2);
    k_det_HM_grid_sum   = k_det_HM_grid_sum    / (4.0 * 3.14159 * dilRad2);
    k_diss_COI_grid_sum = k_diss_COI_grid_sum  / (4.0 * 3.14159 * dilRad2);
    k_ion_CI_grid_sum = k_ion_CI_grid_sum  / (4.0 * 3.14159 * dilRad2);
    k_ion_OI_grid_sum = k_ion_OI_grid_sum  / (4.0 * 3.14159 * dilRad2);

    /* Define Grid */
    k_diss_H2_grid  = new float[size];
    k_det_HM_grid  = new float[size];
    k_diss_CO_grid  = new float[size];
    k_ion_CI_grid  = new float[size];
    k_ion_OI_grid  = new float[size];

    for (int i = 0; i < size; i++){
      k_diss_H2_grid[i] = k_diss_H2I_grid_sum;
      k_det_HM_grid[i] = k_det_HM_grid_sum; 
      k_diss_CO_grid[i] = k_diss_COI_grid_sum;
      k_ion_CI_grid[i] = k_ion_CI_grid_sum;
      k_ion_OI_grid[i] = k_ion_OI_grid_sum;

    }

    if (k_diss_H2I_grid_sum>0){
      fprintf(stdout, "Grid %"ISYM",", this->ID);
      fprintf(stdout, "Time %"ESYM",", this->Time);
      fprintf(stdout, " k_diss_H2 = %"ESYM" 1/CodeTime,", k_diss_H2_grid[0]);
      fprintf(stdout, " k_det_HM  = %"ESYM" 1/CodeTime", k_det_HM_grid[0]);
      fprintf(stdout, " k_diss_CO  = %"ESYM" 1/CodeTime", k_diss_CO_grid[0]);
      fprintf(stdout, " k_ion_CI  = %"ESYM" 1/CodeTime", k_ion_CI_grid[0]);
      fprintf(stdout, " k_ion_OI  = %"ESYM" 1/CodeTime\n", k_ion_OI_grid[0]);

    }

    //To Do -> Add CO, CI, OI to grackle rates?
    my_fields.RT_H2_dissociation_rate =  k_diss_H2_grid;
#ifdef HM_GRACKLE
    my_fields.RT_HM_detachment_rate   =  k_det_HM_grid; //Feeds in Britton's Grackle Branch (foggie-sf) only
    //The following rates are commented out for now, but will feed into the newchem-cpp branch of grackle when ready
    //my_fields.RT_CO_dissociation_rate =  k_diss_CO_grid; 
    //my_fields.RT_CI_ionization_rate   =  k_ion_CI_grid; 
    //my_fields.RT_OI_ionization_rate   =  k_ion_OI_grid; 
#endif

    // Need to set the other fields to the same 0 array for now
    EmptyRtArray  = new float[size];
    for (int i = 0; i < size; i++){
      EmptyRtArray[i] = 0;
    }

    my_fields.RT_HI_ionization_rate   = EmptyRtArray;
    my_fields.RT_HeI_ionization_rate  = EmptyRtArray;
    my_fields.RT_HeII_ionization_rate = EmptyRtArray;
    my_fields.RT_heating_rate = EmptyRtArray;
  } // UseLocalStellarRadiation
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
#endif //TRANSFER


  delete [] TotalMetals;
  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

    if (UseLocalStellarRadiation){
      delete[] k_diss_H2_grid;
      delete[] k_det_HM_grid;
      delete[] k_diss_CO_grid;
      delete[] k_ion_CI_grid;
      delete[] k_ion_OI_grid;
      delete[] EmptyRtArray;
    }


  LCAPERF_STOP("grid_GrackleWrapper");

#endif
  return SUCCESS;
}
