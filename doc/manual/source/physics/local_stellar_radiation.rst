.. _local_stellar_radiation:

Local Stellar Radiation
=========================

This is a cheap alternative to existing radiation transfer methods
to get local photodissociation/photodetachment/photoionization rates
from young stars. These rates are then fed into grackle. Each star
only affects gas cells on their grid and assumes an average particle
cell separation. Luminosities from stars are loaded from
``input/preSN_feedback_SB99_RT.hdf5`` table set by
``StarFeedbackPreSNFilename`` and calculated with
``input/combine_pySB_hdf5_with_RT_Fields.py``. Rates are
pre-mutliplied by the reaction cross section, and thus have units
of cm^2/s. This option sums the relative reaction rates for each star
particle on a grid, and calculates the typical reaction rate as a
grid averaged value. Grids without any young star particles will not
see any radiation. Currently supports H2 Photodissociation,
HM photodetachment, and the isrf radiation field of grackle. Has
commented out support for CO photodissociation, as well as CI and
OI photoionization for when these are added to grackle. Calculations
are done in ``Grid_ComputeLocalStellarRadiation.C``, which is called
in both ``Grid_GraclkeWrapper.C`` and ``Grid_ComputeCoolingTime.C``.

Locally Extinct Radiation
=========================

This is a slightly more extensive, but more robust alternative to
Local Stellar Radiation. This method calculates the H2 photodissociation,
HM photodetachment, and ISRF for each young star particle in the same
manner as Local Stellar Radiation. Luminosities from stars are loaded from
``input/preSN_feedback_SB99_RT.hdf5`` table set by
``StarFeedbackPreSNFilename`` and calculated with
``input/combine_pySB_hdf5_with_RT_Fields.py``. Rates are
pre-mutliplied by the reaction cross section, and thus have units
of cm^2/s. These luminosities are then deposited onto the grid during the
gravity solve. Luminosities are converted to fluxes across the entire
simulation domain by convolving with a 1/4*pi*r^2 kernel, re-using
the existing machinery for the gravity solver as much as possible.
The initial 1/r^2 calculation on the root grid sets the coarsest
radiation field. Deeper layers perform the 1/r^2 calculation using
only star particles on their gravity mesh, and update their radiation fields.
Currently supports H2 Photodissociation, HM photodetachment, and
the isrf radiation field of grackle.

Compiling and Running
=======================
Currently requires the `foggie-sf <https://github.com/foggie-sims/grackle>`_
branch of grackle, which includes a field for HM-photodetachment. To compile,
add the ``HM_GRACKLE`` flag to indicate you are on this branch of grackle
(e.g., add ``-DHM_GRACKLE`` to the ``MACH_DEFINES`` line of your machine file).
To run, set ``UseLocalStellarRadiation=1`` or ``UseLocallyExtinctRadiation=1``.
Requires setting ``StarFeedbackPreSNFilename`` in the parameter file.
