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
pre-multiplied by the reaction cross section, and thus have units
of cm^2/s. This option sums the relative reaction rates for each star
particle on a grid, and calculates the typical reaction rate as a
grid averaged value. Grids without any young star particles will not
see any radiation. Currently supports H2 Photodissociation,
HM photodetachment, and the isrf radiation field of grackle. Has
commented out support for CO photodissociation, as well as CI and
OI photoionization for when these are added to grackle. Calculations
are done in ``Grid_ComputeLocalStellarRadiation.C``, which is called
in both ``Grid_GrackleWrapper.C`` and ``Grid_ComputeCoolingTime.C``.

Locally Extinct Radiation
=========================

This is a slightly more expensive, but more robust alternative to
Local Stellar Radiation. This method calculates the H2 photodissociation,
HM photodetachment, and ISRF for each young star particle in the same
manner as Local Stellar Radiation. Luminosities from stars are loaded from
``input/preSN_feedback_SB99_RT.hdf5`` table set by
``StarFeedbackPreSNFilename`` and calculated with
``input/combine_pySB_hdf5_with_RT_Fields.py``. Rates are calculated from
`Starburst99 <https://github.com/CalumHawcroft/Starburst>`_ models for various
ages and metallicities, pre-multiplied by their reaction cross section, and
integrated across each band (units of cm^2/s). Integrated dust opacities are
also calculated for these rates using tables from
`Draine 2003 <https://www.astro.princeton.edu/~draine/dust/extcurvs/kext_albedo_WD_MW_3.1_60_D03.all>`_.
These rates are then deposited onto the grid during the gravity solve. Rates
are attenuated based on the table's kappa values and the local Sobolev column
density approximation of the gas (done in the same way as in Grackle).
Luminosities are converted to fluxes across the entire
simulation domain by convolving with a 1/(4 pi r^2) kernel, re-using
the existing machinery for the gravity solver as much as possible.
The initial 1/r^2 calculation on the root grid sets the coarsest
radiation field. Deeper layers perform the 1/r^2 calculation using
star particles on their gravity mesh and subtract their parent contribution.
Currently supports H2 Photodissociation, HM photodetachment, and
the isrf radiation field of grackle.

Only star particles younger than 50 Myr are considered, and luminosity
is scaled by their initial mass.

Local extinction is calculated using the Sobolev approximation, and dust
opacities are scaled linearly with metallicity as follows:
``tau = kappa (Z/Zsun) rho L``, where kappa is the tabular dust opacity, Z is the metallicity (capped at 0.9), and L
is the Sobolev length (capped at 1 code length).

Self shielding at absorption is left to Grackle.

Outputs
--------

Both Local Stellar Radiation and Locally Extinct Radiation primarily work by feeding
values to grackle. Local Stellar Radiation does not have any additional fields, but
grid averaged values for radiation can be read in the debug outputs.

Locally Extinct Radiation writes out additional fields for each band:
``LER_kdissH2`` and ``LER_kdetHM`` are in units of 1/s. ``LER_ISRF`` is in Habing units.
These are the rates given to grackle prior to any self-shielding.

Cost and Debugging
-------------------

The Locally Extinct Radiation module essentially adds three more gravity solves. For forced refinement runs,
wherein there are many cells, the time is dominated by the particle deposition step,
and essentially quadruples the time of the gravity solve. In runs with fewer cells,
the time is dominated by the 1/r^2 solves at high refinement levels. This can be alleviated
almost entirely by setting the parameter ``LocallyExtinctRadiationBoxRadius=6``, which
limits the 1/r^2 calculation on finer grids to a 6 cell radius around each star particle.

The timers ``LER_Deposit``, ``LER_RootFFT``, and ``LER_Subgrid`` are available to profile the
deposition, root-grid, and subgrid steps, respectively.

Compiling and Running
----------------------

Currently requires the `foggie-sf <https://github.com/foggie-sims/grackle>`_
branch of grackle, which includes a field for HM-photodetachment. To compile,
add the ``HM_GRACKLE`` flag to indicate you are on this branch of grackle
(e.g., add ``-DHM_GRACKLE`` to the ``MACH_DEFINES`` line of your machine file).

To run, set ``UseLocalStellarRadiation=1`` or ``UseLocallyExtinctRadiation=1``.
Users may additionally set ``LocallyExtinctRadiationBoxRadius=6`` for a
slightly more efficient calculation of the 1/r^2 kernel.
Requires setting ``StarFeedbackPreSNFilename`` in the parameter file.

``UseLocallyExtinctRadiation=1`` also requires the following (each is
checked when the parameter file is read):

* ``UseLocalStellarRadiation=0`` (mutually exclusive).
* ``SelfGravity=1``, the radiation reuses the gravity solver's
  meshes and FFTs.
* ``MaximumGravityRefinementLevel`` >= ``MaximumRefinementLevel``, since no
  radiation is computed on levels above the gravity solve.
* ``StarMakerStoreInitialMass=1``, luminosities are scaled by each
  star particle's initial mass.
* ``TopGridRank=3``.
* No active particles (``AppendActiveParticleType``).
* A ``StarFeedbackPreSNFilename`` table that includes the dust opacity
  (``*_kappa``) datasets, as written by
  ``input/combine_pySB_hdf5_with_RT_Fields.py``.

To have Grackle use the ISRF (e.g. for dust heating), set
``use_isrf_field=1`` in the parameter file (read by
``GrackleReadParameters.C``).

Source Files
-------------

Shared by both methods:

* ``input/combine_pySB_hdf5_with_RT_Fields.py``: builds the
  ``input/preSN_feedback_SB99_RT.hdf5`` table (luminosities and, for
  Locally Extinct Radiation, dust opacities) from Starburst99 outputs.
* ``src/enzo/ReadPreSNFeedbackTable.C``: reads the table and
  interpolates it in age and metallicity (``InterpolatePreSNFeedbackTable``).
* ``src/enzo/Grid_GrackleWrapper.C`` and ``src/enzo/Grid_ComputeCoolingTime.C``:
  hand the rates to Grackle.
* ``src/enzo/GrackleReadParameters.C``: turns on Grackle's radiative
  transfer rates when either method is used.

Local Stellar Radiation (``UseLocalStellarRadiation``):

* ``src/enzo/Grid_ComputeLocalStellarRadiation.C``: grid-averaged rates
  from the young stars on each grid.

Locally Extinct Radiation (``UseLocallyExtinctRadiation``):

* ``src/enzo/DepositParticleLuminosityField.C`` and
  ``src/enzo/Grid_DepositParticlePositions.C``: deposit young-star
  luminosities, with the source-side dust attenuation.
* ``src/enzo/ComputePotentialFieldLevelZero.C`` and
  ``src/enzo/PrepareIsolatedGreensFunction.C``: root-grid 1/r^2
  convolution by FFT.
* ``src/enzo/Grid_SolveForLocallyExtinctRadiation.C``: subgrid correction
  on each finer level, and conversion of the flux to rates for Grackle.
* ``src/enzo/EvolveLevel.C`` and ``src/enzo/PrepareDensityField.C``: call
  the deposit and the subgrid solve alongside gravity.
* Gravity-solver plumbing that carries the source and flux fields between
  grids and processors: ``PrepareGravitatingMassField.C``,
  ``Grid_AddOverlappingParticleMassField.C``, ``Grid_PreparePotentialField.C``,
  ``Grid_CopyPotentialField.C``, ``Grid_CommunicationSendRegion.C``,
  ``Grid_PrepareFFT.C`` and ``Grid_FinishFFT.C``.
* ``src/enzo/New_Grid_WriteGrid.C``: writes the ``LER_kdissH2``,
  ``LER_kdetHM`` and ``LER_ISRF`` output fields.
