/***********************************************************************
/
/  WRITES GRACKLE PARAMETERS FROM INPUT FILE
/
/  written by: Andrew Emerick
/  date:       April, 2020
/  modified1:
/
/  PURPOSE:
/
/  NOTE: Handles writing of Grackle-specific parameters to file
/
************************************************************************/

#include "preincludes.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* function prototypes */

int GrackleWriteParameters(FILE *fptr)
{

#ifndef USE_GRACKLE

  if (use_grackle == TRUE){
    ENZO_FAIL("Error: Enzo must be compiled with 'make grackle-yes' to run with use_grackle = 1");
  }

#else

  fprintf(fptr, "with_radiative_cooling      = %d\n", grackle_data->with_radiative_cooling);
  fprintf(fptr, "use_volumetric_heating_rate = %d\n", grackle_data->use_volumetric_heating_rate);
  fprintf(fptr, "use_specific_heating_rate   = %d\n", grackle_data->use_specific_heating_rate);
  fprintf(fptr, "self_shielding_method       = %d\n", grackle_data->self_shielding_method);
  fprintf(fptr, "H2_self_shielding           = %d\n", grackle_data->H2_self_shielding);
  fprintf(fptr, "grackle_data_file           = %s\n", grackle_data->grackle_data_file);
  fprintf(fptr, "UVbackground                = %d\n", grackle_data->UVbackground);
  fprintf(fptr, "Compton_xray_heating        = %d\n", grackle_data->Compton_xray_heating);
  fprintf(fptr, "LWbackground_intensity      = %lf\n", grackle_data->LWbackground_intensity);
  fprintf(fptr, "LWbackground_sawtooth_suppression = %d\n", grackle_data->LWbackground_sawtooth_suppression);
  fprintf(fptr, "dust_chemistry              = %d\n",  grackle_data->dust_chemistry);
  fprintf(fptr, "local_dust_to_gas_ratio     = %lf\n", grackle_data->local_dust_to_gas_ratio);
  fprintf(fptr, "interstellar_radiation_field = %lf\n", grackle_data->interstellar_radiation_field);
  fprintf(fptr, "dust_recombination_cooling  = %d\n",  grackle_data->dust_recombination_cooling);
  fprintf(fptr, "use_isrf_field              = %d\n",  grackle_data->use_isrf_field);
  fprintf(fptr, "use_dust_density_field      = %d\n",  grackle_data->use_dust_density_field);

  /* New dust physics parameters (newchemcpp Grackle) */
  fprintf(fptr, "dust_model                  = %d\n",  grackle_data->dust_model);
  fprintf(fptr, "solver_method               = %d\n",  grackle_data->solver_method);
  fprintf(fptr, "use_sne_field               = %d\n",  grackle_data->use_sne_field);
  fprintf(fptr, "use_tau_dest_field          = %d\n",  grackle_data->use_tau_dest_field);
  fprintf(fptr, "dust_destruction_eff        = %lf\n", grackle_data->dust_destruction_eff);
  fprintf(fptr, "sne_coeff                   = %lf\n", grackle_data->sne_coeff);
  fprintf(fptr, "sne_shockspeed              = %lf\n", grackle_data->sne_shockspeed);
  fprintf(fptr, "dust_grainsize              = %lf\n", grackle_data->dust_grainsize);
  fprintf(fptr, "dust_growth_densref         = %lf\n", grackle_data->dust_growth_densref);
  fprintf(fptr, "dust_growth_tauref          = %lf\n", grackle_data->dust_growth_tauref);
  fprintf(fptr, "dust_condensation_eff       = %lf\n", grackle_data->dust_condensation_eff);
  fprintf(fptr, "sne_metal_yield             = %lf\n", grackle_data->sne_metal_yield);

  /* Species-resolved dust tracking (Trayford+2026 MNRAS 545, staf2040).
     dust_species_track is mapped from the Enzo flag UseDustSpeciesTrack,
     written by WriteParameterFile. */
  fprintf(fptr, "dust_growth_sticking_coeff  = %lf\n", grackle_data->dust_growth_sticking_coeff);
  fprintf(fptr, "dust_growth_tauref_silicate = %lf\n", grackle_data->dust_growth_tauref_silicate);
  fprintf(fptr, "dust_growth_tauref_carbon   = %lf\n", grackle_data->dust_growth_tauref_carbon);
  fprintf(fptr, "dust_growth_clumping_factor_max = %lf\n", grackle_data->dust_growth_clumping_factor_max);
  fprintf(fptr, "dust_growth_clumping_nH_min = %lf\n", grackle_data->dust_growth_clumping_nH_min);
  fprintf(fptr, "dust_growth_clumping_nH_max = %lf\n", grackle_data->dust_growth_clumping_nH_max);
  fprintf(fptr, "dust_sputter_tauref         = %lf\n", grackle_data->dust_sputter_tauref);
  fprintf(fptr, "dust_silicate_mg_fraction   = %lf\n", grackle_data->dust_silicate_mg_fraction);

#endif

  return SUCCESS;
}
