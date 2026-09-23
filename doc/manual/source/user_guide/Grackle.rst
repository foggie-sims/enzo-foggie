.. _Grackle:

Running Enzo with Grackle
=========================

The Grackle is an external chemistry and cooling library originally derived from 
Enzo's MultiSpecies chemistry and Cloudy cooling modules.  The non-equilibrium 
primordial chemistry and cooling functionality is essentially identical to the 
MultiSpecies network.  However, significant updates have been made to the treatment 
of metals and UV backgrounds that may make using the Grackle a more attactive option, 
such as:

- UV backgrounds are treated via interpolating from data tables loaded from disk rather 
  than piece-wise polynomial functions, making it somewhat easier to add support for 
  new background models.

- UV background and cooling data are contained within the same input file and are 
  more consistent than the currently available Cloudy cooling data and Enzo UV 
  background models.  Currently, the Grackle distribution comes with UV background 
  and cooling data for two different models:

    1. `Faucher-Giguere et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...703.1416F>`_.

    2. `Haardt & Madau (2012) <http://adsabs.harvard.edu/abs/2012ApJ...746..125H>`_.

- Unlike the original Cloudy cooling which required separate input files for cooling 
  before the UV background turns on and after, all data is contained in a single 
  table.  This means one no longer has to run with one input file to the redshift 
  where the UV background starts, stop the simulation, and restart with another input 
  file.

- When the UV background comes on is controlled by Enzo, not by the redshift range of
  the table.  Enzo passes its ``RadiationRedshiftOn`` (default 7.0) and
  ``RadiationRedshiftFullOn`` (default 6.0) parameters to Grackle, which ramps the
  background from zero to full strength between those two redshifts even if the table
  starts at a higher redshift (the Haardt & Madau 2012 table starts at z = 15.13).
  Set both parameters to the table's maximum redshift to use the whole table.  The
  values in effect are printed at start-up; see the Grackle section of the parameter
  list for details.

- Also unlike the original Cloudy cooling module, Grackle supports the option to also 
  solve the primordial cooling via interpolation from a table.  Thus, one is no longer 
  required to run with the MultiSpecies functionality in order to calculate the 
  primordial component.  This simplified method is somewhat faster and requires fewer 
  baryon fields to be stored, lowering the ram and disk footprint.

For more information on the Grackle library, see the 
`Grackle documentation <https://grackle.readthedocs.org/>`_.

Obtaining and Building the Grackle
----------------------------------

See the `Grackle documentation <https://grackle.readthedocs.org/>`_ for complete 
instruction on how to obtain, compile, and install the Grackle library.

Compiling Enzo with Grackle
---------------------------

In order to compile Enzo with support for Grackle, the following lines need to be added 
to your machine make file in the appropriate places:

::

   LOCAL_GRACKLE_INSTALL = PATH/TO/GRACKLE
   LOCAL_INCLUDES_GRACKLE = -I$(LOCAL_GRACKLE_INSTALL)/include
   MACH_INCLUDES_GRACKLE = $(LOCAL_INCLUDES_GRACKLE)
   LOCAL_LIBS_GRACKLE = -L$(LOCAL_GRACKLE_INSTALL)/lib -lgrackle
   MACH_LIBS_GRACKLE = $(LOCAL_LIBS_GRACKLE)

See the example make file, **Make.mach.unkown**, in the Enzo source for an example.

To configuration Enzo to build with Grackle support, do the following before typing 
"make":

::

   make grackle-yes

Running with the Grackle
------------------------

Grackle parameters should be given in the same parameter file as the rest of the Enzo 
parameters.  Since the Grackle is based on Enzo's MultiSpecies, many of the parameter 
names are the same.  For a full list of Grackle parameters, see 
:ref:`here <grackle_pars>`.
