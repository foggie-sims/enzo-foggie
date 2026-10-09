/***********************************************************************
/
/  GRID CLASS (COMPUTE LOCAL STELLAR RADIATION FIELD)
/
/  written by: Cameron Trapp
/  date:       June, 2026
/  modified1:
/
/  PURPOSE: Estimate the grid-averaged photodissociation/ionization rates
/           and ISRF from young star particles on this grid, using the
/           pre-SN feedback (SB99) tables. Results are stored in the grid
/           attributes k_diss_H2I_grid_sum, k_det_HM_grid_sum,
/           k_diss_COI_grid_sum, k_ion_CI_grid_sum, k_ion_OI_grid_sum
/           and isrf_grid_sum (rates in 1/CodeTime, isrf in Habing units).
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

/* function prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int search_lower_bound(float *arr, float value, int low, int high, 
		       int total);

int grid::ComputeLocalStellarRadiation()
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  /* Estimate local radiation field from new Stars - CWT 06/07/26 */
  //Sum all the mass of all young star particles on the grid
  //Get an estimate of the LW photon production from fitting results 
  //Convert to RT_H2_dissociation_rate and add to grackle fields
  float MassUnits = DensityUnits * POW(LengthUnits,3);

  //n_metal_bins = pSNFBTable.n_met
  //n_age_bins = pSNFBTable.n_age
  //metallicity bins = pSNFBTable.ini_met?
  //age bins = pSNFBTable.pop_age
  //kdiss_H2_sb99 = pSNFBTable.kdiss_H2
  //kdet_HM = pSNFBTable.kdet_HM

  k_diss_H2I_grid_sum = 0; //Grid Attributes 
  k_det_HM_grid_sum = 0;
  k_diss_COI_grid_sum = 0;
  k_ion_CI_grid_sum = 0;
  k_ion_OI_grid_sum = 0;
  isrf_grid_sum = 0;
  float dt_table = pSNFBTable.pop_age[1] - pSNFBTable.pop_age[0];
  for (i = 0; i < this->NumberOfParticles; i++) {
    if (this->ParticleType[i] == PARTICLE_TYPE_STAR) {
      float age = (this->Time - this->ParticleAttribute[0][i]) * TimeUnits / yr_s; //Convert to yr
      if (age < 5e7) { 
        //Interpolate in age and metallicity to get the rates
        int aa = (int)floor((age - pSNFBTable.pop_age[0]) / dt_table); //Index of the lower bound for age

        float t_age; //How far between the two indices?
        if (aa>=pSNFBTable.n_age-1){
          aa=pSNFBTable.n_age-2;
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
        int zz = search_lower_bound((float*)pSNFBTable.ini_met, metallicity, 0, pSNFBTable.n_met, pSNFBTable.n_met); //Index of the lower bound for metallicity

        float t_z; //How far between the two indices?
        if (zz>=pSNFBTable.n_met-1){
          zz=pSNFBTable.n_met-2;
          t_z = 1;
        }
        else if (zz<0){
          zz=0;
          t_z=0;
        }
        else{
            t_z = (metallicity - pSNFBTable.ini_met[zz]) / (pSNFBTable.ini_met[zz+1] - pSNFBTable.ini_met[zz]);
        }

        //bilinear interpolation in age and metallicity to get the rates
        int ii0 = zz * pSNFBTable.n_age + aa;
        int ii1 = zz * pSNFBTable.n_age + (aa+1);
        int ii2 = (zz+1) * pSNFBTable.n_age + aa;
        int ii3 = (zz+1) * pSNFBTable.n_age + (aa+1);

        /* In Units Hz cm^2 per Solar Mass*/
        float k_diss_H2_sb99_interp = (1-t_age) * (1-t_z) * pSNFBTable.kdiss_H2[ii0] + t_age * (1-t_z) * pSNFBTable.kdiss_H2[ii1] + (1-t_age) * t_z * pSNFBTable.kdiss_H2[ii2] + t_age * t_z * pSNFBTable.kdiss_H2[ii3];
        float k_det_HM_sb99_interp  = (1-t_age) * (1-t_z) * pSNFBTable.kdet_HM[ii0] + t_age * (1-t_z) * pSNFBTable.kdet_HM[ii1] + (1-t_age) * t_z * pSNFBTable.kdet_HM[ii2] + t_age * t_z * pSNFBTable.kdet_HM[ii3];
        float k_diss_CO_sb99_interp = (1-t_age) * (1-t_z) * pSNFBTable.kdiss_CO[ii0] + t_age * (1-t_z) * pSNFBTable.kdiss_CO[ii1] + (1-t_age) * t_z * pSNFBTable.kdiss_CO[ii2] + t_age * t_z * pSNFBTable.kdiss_CO[ii3];
        float k_ion_CI_sb99_interp  = (1-t_age) * (1-t_z) * pSNFBTable.kion_CI[ii0] + t_age * (1-t_z) * pSNFBTable.kion_CI[ii1] + (1-t_age) * t_z * pSNFBTable.kion_CI[ii2] + t_age * t_z * pSNFBTable.kion_CI[ii3];
        float k_ion_OI_sb99_interp  = (1-t_age) * (1-t_z) * pSNFBTable.kion_OI[ii0] + t_age * (1-t_z) * pSNFBTable.kion_OI[ii1] + (1-t_age) * t_z * pSNFBTable.kion_OI[ii2] + t_age * t_z * pSNFBTable.kion_OI[ii3];
        float isrf_sb99_interp      = (1-t_age) * (1-t_z) * pSNFBTable.isrf[ii0] + t_age * (1-t_z) * pSNFBTable.isrf[ii1] + (1-t_age) * t_z * pSNFBTable.isrf[ii2] + t_age * t_z * pSNFBTable.isrf[ii3];


        k_diss_H2I_grid_sum += k_diss_H2_sb99_interp * this->ParticleInitialMass[i];
        k_det_HM_grid_sum += k_det_HM_sb99_interp * this->ParticleInitialMass[i];
        k_diss_COI_grid_sum += k_diss_CO_sb99_interp * this->ParticleInitialMass[i];
        k_ion_CI_grid_sum += k_ion_CI_sb99_interp * this->ParticleInitialMass[i];
        k_ion_OI_grid_sum += k_ion_OI_sb99_interp * this->ParticleInitialMass[i];
        isrf_grid_sum += isrf_sb99_interp * this->ParticleInitialMass[i];



      }
    }
  }

  float dx = this->CellWidth[0][0];
  float conversion_factor = dx * dx * dx * MassUnits / SolarMass; //Convert particle density to units of Msun
  conversion_factor = conversion_factor * TimeUnits/(LengthUnits*LengthUnits); //Convert from cm^2/s to code units

  k_diss_H2I_grid_sum = k_diss_H2I_grid_sum * conversion_factor;
  k_det_HM_grid_sum  = k_det_HM_grid_sum  * conversion_factor;
  k_diss_COI_grid_sum = k_diss_COI_grid_sum  * conversion_factor;
  k_ion_CI_grid_sum = k_ion_CI_grid_sum * conversion_factor;
  k_ion_OI_grid_sum = k_ion_OI_grid_sum * conversion_factor;
  isrf_grid_sum = isrf_grid_sum * conversion_factor;


  float grid_dx = this->GridRightEdge[0]-this->GridLeftEdge[0];
  float grid_dy = this->GridRightEdge[1]-this->GridLeftEdge[1];
  float grid_dz = this->GridRightEdge[2]-this->GridLeftEdge[2];
  //This is the most tunable part of this code, as calculating the r^2 for each cell will get expensive
  float dilutionRadius = sqrt((grid_dx*grid_dx + grid_dy*grid_dy + grid_dz*grid_dz) / 6.0); //~point by point separation
  float dilRad2 = dilutionRadius * dilutionRadius;
  k_diss_H2I_grid_sum = k_diss_H2I_grid_sum  / (4.0 * 3.14159 * dilRad2);
  k_det_HM_grid_sum   = k_det_HM_grid_sum    / (4.0 * 3.14159 * dilRad2);
  k_diss_COI_grid_sum = k_diss_COI_grid_sum  / (4.0 * 3.14159 * dilRad2);
  k_ion_CI_grid_sum = k_ion_CI_grid_sum  / (4.0 * 3.14159 * dilRad2);
  k_ion_OI_grid_sum = k_ion_OI_grid_sum  / (4.0 * 3.14159 * dilRad2);
  isrf_grid_sum = isrf_grid_sum  / (4.0 * 3.14159 * dilRad2); //Convert to G0

  if (debug){
    if (k_diss_H2I_grid_sum>0){
      fprintf(stdout, "Grid %"ISYM",", this->ID);
      fprintf(stdout, "Time %"ESYM",", this->Time);
      fprintf(stdout, " k_diss_H2 = %"ESYM" 1/CodeTime,", k_diss_H2I_grid_sum);
      fprintf(stdout, " k_det_HM  = %"ESYM" 1/CodeTime,", k_det_HM_grid_sum);
      fprintf(stdout, " k_diss_CO  = %"ESYM" 1/CodeTime,", k_diss_COI_grid_sum);
      fprintf(stdout, " k_ion_CI  = %"ESYM" 1/CodeTime,", k_ion_CI_grid_sum);
      fprintf(stdout, " k_ion_OI  = %"ESYM" 1/CodeTime,", k_ion_OI_grid_sum);
      fprintf(stdout, " isrf  = %"ESYM" (Habing Units)\n", isrf_grid_sum);
    }
  }

  return SUCCESS;
}
