.. _local_stellar_radiation:

Local Stellar Radiation
=========================

This is a cheap alternative to existing radiation transfer methods
to get local photodissociation/photodetachment/photoionization rates
from young stars. These rates are then fed into grackle. Each star
only affects gas cells on their grid and assumes an average particle
cell separation. Fluxes from stars are loaded from
`input/preSN_feedback_SB99_RT.hdf5` table set by
`StarFeedbackPreSNFilename` and calculated with
`input/combine_pySB_hdf5_with_RT_Fields.py` (or interactively with 
`input/combine_pySB_hdf5_with_RT_Fields.ipynb`). Rates are
pre-mutliplied by the reaction cross section, and thus have units
of cm^2/s. This option sums the relative reaction rates for each star
particle on a grid, and calculates the typical reaction rate as a
grid averaged value. Grids without any young star particles will not
see any radiation. Currently supports H2 Photodissociation,
HM photodetachment, and  the isrf radiation field of grackle. Has
commented out support for CO photodissociation, as well as CI and
OI photoionization for when these are added to grackle. Calculations
are done in `Grid_ComputeLocalStellarRadiation.C`, which is called
in both `Grid_GraclkeWrapper.C` and `Grid_ComputeCoolingTime.C`.


Compiling and Running
=======================
Currently requires the `foggie-sf (https://github.com/foggie-sims/grackle)`
branch of grackle, which includes a field for HM-photodetachment. To compile,
add the `HM_GRACKLE` flag to indicate you are on this branch of grackle.
To run, set `UseLocalStellarRadiation=1` and `StarFeedbackPreSNFilename`
in the parameter file.
