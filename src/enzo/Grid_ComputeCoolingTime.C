/***********************************************************************
/
/  GRID CLASS (COMPUTE THE COOLING TIME FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Elizabeth Harper-Clark, August 2009
/              added in CoolingModel parameter
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the cooling time

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "phys_constants.h"
 
/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */
 
extern int RadiationFieldRecomputeMetalRates;
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);
int GadgetCoolingTime(float *d, float *e, float *ge, 
		      float *u, float *v, float *w,
		      float *cooltime,
		      int *in, int *jn, int *kn, int *iexpand, 
		      hydro_method *imethod, int *idual, int *idim,
		      int *is, int *js, int *ks, int *ie, int *je, 
		      int *ke, float *dt, float *aye,
		      float *fh, float *utem, float *uxyz, 
		      float *uaye, float *urho, float *utim,
		      float *gamma);
int search_lower_bound(float *arr, float value, int low, int high, 
		       int total);

extern "C" void FORTRAN_NAME(cool_multi_time)(
	float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
	float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
	float *cooltime, int *coolonly,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
	hydro_method *imethod,
        int *idual, int *ispecies, int *imetal, int *imcool, int *idust, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co,
	int *ipiht, int *igammah,
	float *dt, float *aye, float *redshift, float *temstart, float *temend,
	float *utem, float *uxyz, float *uaye, float *urho, float *utim,
	float *eta1, float *eta2, float *gamma, float *z_solar,
	float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, float *ciHeIa,
	float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a,
	float *reHeII2a, float *reHeIIIa, float *brema, float *compa, float *gammaha,
	float *comp_xraya, float *comp_temp, float *piHI, float *piHeI, float *piHeII,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI, float *metal,
	float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
	float *gpldl, float *gphdl, float *HDltea, float *HDlowa,
	float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
	float *gasgra, float *metala, int *n_xe, float *xe_start, float *xe_end,
	float *inutot, int *iradfield, int *nfreq, int *imetalregen,
	int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p,
	int *iradtrans, float *photogamma,
	int *ih2optical, int *iciecool, float *ciecoa,
 	int *icmbTfloor, int *iClHeat,
 	float *clEleFra, int *clGridRank, int *clGridDim,
 	float *clPar1, float *clPar2, float *clPar3, float *clPar4, float *clPar5,
 	int *clDataSize, float *clCooling, float *clHeating);

extern "C" void FORTRAN_NAME(cool_time)(
	float *d, float *e, float *ge, float *u, float *v, float *w,
           float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
	hydro_method *imethod, int *idual, int *idim, int *igammah,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
	float *dt, float *aye, float *temstart, float *temend,
	float *fh, float *utem, float *urho, 
	float *eta1, float *eta2, float *gamma, float *coola, float *gammaha, float *mu);
 
int grid::ComputeCoolingTime(float *cooling_time, int CoolingTimeOnly)
{
 
  /* Return if this doesn't concern us. */

  if (RadiativeCooling == 0) return SUCCESS;
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
 
  /* Compute the size of the fields. */
 
  int i;
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
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
 
  /* Find photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

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
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
 
    aUnits = 1.0/(1.0 + InitialRedshift);
  } else if (RadiationFieldRedshift > -1){
    a       = 1.0 / (1.0 + RadiationFieldRedshift);
    aUnits  = 1.0;
  }
  float afloat = float(a);
 
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
    if (debug)
      fprintf(stderr, "Warning: No metal field found.  Turning OFF MetalCooling.\n");
    MetalCooling = FALSE;
    MetalNum = 0;
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
 
#ifdef USE_GRACKLE
  if (grackle_data->use_grackle == TRUE) {

    Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
    g_grid_dimension = new Eint32[GridRank];
    g_grid_start = new Eint32[GridRank];
    g_grid_end = new Eint32[GridRank];
    for (i = 0; i < GridRank; i++) {
      g_grid_dimension[i] = (Eint32) GridDimension[i];
      g_grid_start[i] = 0;
      g_grid_end[i] = (Eint32) GridDimension[i]-1;
    }

    /* Update units. */

    code_units grackle_units;
    grackle_units.comoving_coordinates = (Eint32) ComovingCoordinates;
    grackle_units.density_units        = (double) DensityUnits;
    grackle_units.length_units         = (double) LengthUnits;
    grackle_units.time_units           = (double) TimeUnits;
    grackle_units.velocity_units       = (double) VelocityUnits;
    grackle_units.a_units              = (double) aUnits;
    grackle_units.a_value              = (double) a;

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

    /* set up the my_fields */
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

    my_fields.volumetric_heating_rate  = volumetric_heating_rate;
    my_fields.specific_heating_rate    = specific_heating_rate;

#ifdef TRANSFER

    /* unit conversion from Enzo RT units to CGS */
    float rtunits = erg_eV / TimeUnits;

    if ( RadiativeTransfer ){
      my_fields.RT_HI_ionization_rate = BaryonField[kphHINum];

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

  //CWT 05072026
  /* Estimate local radiation field from new Stars */
  //Currently copied from what I put in Grid_GrackleWrapper.C to see if this works, but not very efficient!
  //Simplest Alternative:
  //Sum all the mass of all young star particles on the grid
  //Get an estimate of the LW photon production from fitting Starbust99 results (w/o loading in a table)
  //Convert to RT_H2_dissociation_rate and add to grackle fields

  //Sum all mass of young star particles on grid
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
         // fprintf(stdout, "CWT: Interpolating for star particle with mass %"FSYM" Msun?\n",ParticleMass_Msun);

          k_diss_H2 += k_diss_H2_sb99_interp * ParticleMass_Msun; //Convert from code mass to Msun
          k_det_HM += k_det_HM_sb99_interp * ParticleMass_Msun; //Convert from code mass to Msun


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
  //float dilutionRadius = 4.848e-6 * pc_cm / (double) LengthUnits;  // 1 AU //Try an extreme case first
  float dilRad2 = dilutionRadius * dilutionRadius;
  k_diss_H2 = k_diss_H2  / (4.0 * 3.14159 * dilRad2);
  k_det_HM = k_det_HM    / (4.0 * 3.14159 * dilRad2);

  //float *k_diss_H2_grid = new float[size];
  //float *k_det_HM_grid = new float[size];

  //for (int i = 0; i < size; i++){
  //    k_diss_H2_grid[i] = k_diss_H2;
  //    k_det_HM_grid[i] = k_det_HM;
  //}

  float *k_diss_H2_grid  = new float[size];
  float *k_det_HM_grid  = new float[size];
  for (int i = 0; i < size; i++){
      k_diss_H2_grid[i] = k_diss_H2;
      k_det_HM_grid[i] = k_det_HM;
  }
  fprintf(stdout, "CWT: Setting My Fields in Grid_ComputeCoolingTime...\n");

  my_fields.RT_H2_dissociation_rate =  k_diss_H2_grid;//Already in units of seconds (from table)
  my_fields.RT_HM_detachment_rate =  k_det_HM_grid;  //Feeds in Britton's Grackle Branch (foggie-sf) only


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


    if (calculate_cooling_time(&grackle_units, &my_fields, cooling_time) == FAIL) {
      ENZO_FAIL("Error in Grackle calculate_cooling_time.\n");
    }

    for (i = 0; i < size; i++) {
      cooling_time[i] = fabs(cooling_time[i]);
    }

    if (temp_thermal == TRUE) {
      delete [] thermal_energy;
    }

#ifdef TRANSFER
  if (RadiativeTransfer){
    /* convert the RT units back to Enzo */
    for(i = 0; i < size; i ++) BaryonField[gammaNum][i] /= rtunits;

  }
#endif // TRANSFER

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

    return SUCCESS;
  }
#endif // USE_GRACKLE

  /* Calculate the rates due to the radiation field. */
 
  if (RadiationFieldCalculateRates(Time+0.5*dtFixed) == FAIL) {
    ENZO_FAIL("Error in RadiationFieldCalculateRates.\n");
  }
 
  /* Set up information for rates which depend on the radiation field. 
     Precompute factors for self shielding (this is the cross section * dx). */
 
  float HIShieldFactor = RadiationData.HIAveragePhotoHeatingCrossSection *
                         double(LengthUnits) * CellWidth[0][0];
  float HeIShieldFactor = RadiationData.HeIAveragePhotoHeatingCrossSection *
                          double(LengthUnits) * CellWidth[0][0];
  float HeIIShieldFactor = RadiationData.HeIIAveragePhotoHeatingCrossSection *
                           double(LengthUnits) * CellWidth[0][0];

  /* Calculate start and end indices for cooling time calculations and feed 
     them into the appropriate routine(s).  Unlike many other quantities, we
     want to actually calculate the cooling time for the ghost zones.  This is 
     important for grid refinement when refining by local cooling time, because 
     Enzo flags a buffer region of grid cells around those cells that meet the
     various refinement criteria.  Ghost zones need to be examined for refinement
     criteria so this is done correctly, and thus cooling time needs to be calculated
     in the ghost zones.  Note that we calculate these numbers as unique variables 
     and feed them into the various routines for clarity!  */
    
  int GridStartIndexX,GridStartIndexY,GridStartIndexZ,
    GridEndIndexX,GridEndIndexY,GridEndIndexZ;

  GridStartIndexX = 0;
  GridStartIndexY = 0;
  GridStartIndexZ = 0;

  GridEndIndexX=GridDimension[0]-1;  
  GridEndIndexY=GridDimension[1]-1;
  GridEndIndexZ=GridDimension[2]-1;

  /* Call the appropriate FORTRAN routine to do the work. */

  if (MultiSpecies) {
    // printf("Grid_ComputeCoolingTime.C, integer arguments to cool_multi_time:\n");
    // printf("  in =%"ISYM"\n",GridDimension[0]);
    // printf("  jn =%"ISYM"\n",GridDimension[1]);
    // printf("  kn =%"ISYM"\n",GridDimension[2]);
    // printf("  nratec =%"ISYM"\n",CoolData.NumberOfTemperatureBins);
    // printf("  iexpand =%"ISYM"\n",ComovingCoordinates);
    // printf("  imethod =%"ISYM"\n",HydroMethod);
    // printf("  idual =%"ISYM"\n",DualEnergyFormalism);
    // printf("  ispecies =%"ISYM"\n",MultiSpecies);
    // printf("  imetal =%"ISYM"\n",MetalFieldPresent);
    // printf("  imcool =%"ISYM"\n",MetalCooling);
    // printf("  idust =%"ISYM"\n",H2FormationOnDust);
    // printf("  idim =%"ISYM"\n",GridRank);
    // printf("  is =%"ISYM"\n",GridStartIndex[0]);
    // printf("  js =%"ISYM"\n",GridStartIndex[1]);
    // printf("  ks =%"ISYM"\n",GridStartIndex[2]);
    // printf("  ie =%"ISYM"\n",GridEndIndex[0]);
    // printf("  je =%"ISYM"\n",GridEndIndex[1]);
    // printf("  ke =%"ISYM"\n",GridEndIndex[2]);
    // printf("  ih2co =%"ISYM"\n",CoolData.ih2co);
    // printf("  ipiht =%"ISYM"\n",CoolData.ipiht);
    // printf("  igammah =%"ISYM"\n",PhotoelectricHeating);
    // printf("  n_xe =%"ISYM"\n",CoolData.NumberOfElectronFracBins);
    // printf("  iradfield =%"ISYM"\n",RadiationFieldType);
    // printf("  nfreq =%"ISYM"\n",RadiationData.NumberOfFrequencyBins);
    // printf("  imetalregen =%"ISYM"\n",RadiationFieldRecomputeMetalRates);
    // printf("  iradshield =%"ISYM"\n",RadiationData.RadiationShield);
    // printf("  iradtrans =%"ISYM"\n",RadiativeTransfer);
    // printf("  ih2optical =%"ISYM"\n",H2OpticalDepthApproximation);
    // printf("  iciecool =%"ISYM"\n",CIECooling);
    // printf("  icmbTfloor =%"ISYM"\n",CloudyCoolingData.CMBTemperatureFloor);
    // printf("  iClHeat =%"ISYM"\n",CloudyCoolingData.IncludeCloudyHeating);
    // printf("  clGridRank =%"ISYM"\n",CloudyCoolingData.CloudyCoolingGridRank);
    // printf("  clGridDim =%"ISYM", %"ISYM", %"ISYM", %"ISYM", %"ISYM"\n",
    // 	   CloudyCoolingData.CloudyCoolingGridDimension[0],
    // 	   CloudyCoolingData.CloudyCoolingGridDimension[1],
    // 	   CloudyCoolingData.CloudyCoolingGridDimension[2],
    // 	   CloudyCoolingData.CloudyCoolingGridDimension[3],
    // 	   CloudyCoolingData.CloudyCoolingGridDimension[4]);
    // printf("  clDataSize =%"ISYM"\n\n",CloudyCoolingData.CloudyDataSize);
    

    FORTRAN_NAME(cool_multi_time)(
       density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
       BaryonField[DeNum], BaryonField[HINum], BaryonField[HIINum],
       BaryonField[HeINum], BaryonField[HeIINum], BaryonField[HeIIINum],
       cooling_time, &CoolingTimeOnly,
       GridDimension, GridDimension+1, GridDimension+2,
       &CoolData.NumberOfTemperatureBins, &ComovingCoordinates,
       &HydroMethod,
       &DualEnergyFormalism, &MultiSpecies, &MetalFieldPresent, &MetalCooling, 
       &H2FormationOnDust,
       &GridRank, &GridStartIndexX, &GridStartIndexY, &GridStartIndexZ,
       &GridEndIndexX, &GridEndIndexY, &GridEndIndexZ,
       &CoolData.ih2co, &CoolData.ipiht, &PhotoelectricHeating,
       &dtFixed, &afloat, &RadiationFieldRedshift,
       &CoolData.TemperatureStart, &CoolData.TemperatureEnd,
       &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
       &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
       &CoolData.SolarMetalFractionByMass,
       CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI,
       CoolData.ciHeI,
       CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII,
       CoolData.reHeII1,
       CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &CoolData.comp, &CoolData.gammah,
       &CoolData.comp_xray, &CoolData.temp_xray,
       &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
       BaryonField[HMNum], BaryonField[H2INum], BaryonField[H2IINum],
       BaryonField[DINum], BaryonField[DIINum], BaryonField[HDINum],
       MetalPointer,
       CoolData.hyd01k, CoolData.h2k01, CoolData.vibh,
       CoolData.roth, CoolData.rotl,
       CoolData.GP99LowDensityLimit, CoolData.GP99HighDensityLimit,
       CoolData.HDlte, CoolData.HDlow,
       CoolData.GAHI, CoolData.GAH2, CoolData.GAHe, CoolData.GAHp,
       CoolData.GAel, CoolData.gas_grain,
       CoolData.metals, &CoolData.NumberOfElectronFracBins, 
       &CoolData.ElectronFracStart, &CoolData.ElectronFracEnd,
       RadiationData.Spectrum[0], &RadiationFieldType,
       &RadiationData.NumberOfFrequencyBins,
       &RadiationFieldRecomputeMetalRates,
       &RadiationData.RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor,
       &RadiativeTransfer, BaryonField[gammaNum], 
       &H2OpticalDepthApproximation, &CIECooling, CoolData.cieco,
       &CloudyCoolingData.CMBTemperatureFloor,
       &CloudyCoolingData.IncludeCloudyHeating,
       &CloudyCoolingData.CloudyElectronFractionFactor,
       &CloudyCoolingData.CloudyCoolingGridRank,
       CloudyCoolingData.CloudyCoolingGridDimension,
       CloudyCoolingData.CloudyCoolingGridParameters[0],
       CloudyCoolingData.CloudyCoolingGridParameters[1],
       CloudyCoolingData.CloudyCoolingGridParameters[2],
       CloudyCoolingData.CloudyCoolingGridParameters[3],
       CloudyCoolingData.CloudyCoolingGridParameters[4],
       &CloudyCoolingData.CloudyDataSize,
       CloudyCoolingData.CloudyCooling, CloudyCoolingData.CloudyHeating);
  } else if (GadgetEquilibriumCooling==1) {  
    int result = GadgetCoolingTime
      (
       density,totalenergy,gasenergy,velocity1,
       velocity2,velocity3,
       cooling_time,
       GridDimension,GridDimension+1,
       GridDimension+2, &ComovingCoordinates, &HydroMethod,
       &DualEnergyFormalism, &GridRank,
       &GridStartIndexX, &GridStartIndexY, &GridStartIndexZ,
       &GridEndIndexX, &GridEndIndexY, &GridEndIndexZ,&dtFixed,
       &afloat,&CoolData.HydrogenFractionByMass,
       &TemperatureUnits,&LengthUnits,
       &aUnits,&DensityUnits,&TimeUnits,&Gamma);
    if (result == FAIL )  {

      ENZO_FAIL("Error in GadgetCoolingTime.  Exiting.");
    }
  } else { // if not multispecies or Gadget cooling, must be generic cooling.
    FORTRAN_NAME(cool_time)(
       BaryonField[DensNum], BaryonField[TENum], BaryonField[GENum],
          BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
          cooling_time,
       GridDimension, GridDimension+1, GridDimension+2,
          &CoolData.NumberOfTemperatureBins, &ComovingCoordinates,
          &HydroMethod,
       &DualEnergyFormalism, &GridRank, &PhotoelectricHeating,
       &GridStartIndexX, &GridStartIndexY, &GridStartIndexZ,
       &GridEndIndexX, &GridEndIndexY, &GridEndIndexZ,
       &dtFixed, &afloat, &CoolData.TemperatureStart,
          &CoolData.TemperatureEnd, &CoolData.HydrogenFractionByMass,
       &TemperatureUnits, &DensityUnits,
       &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
       CoolData.EquilibriumRate, &CoolData.gammah, &Mu);
  }

  delete [] TotalMetals;
 
  return SUCCESS;
}
