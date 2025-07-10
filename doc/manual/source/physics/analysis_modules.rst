.. _analysis_modules:

Analysis Modules
================


Tracer Particles
----------------
Enzo has extensive analysis methods to track tracer particles. Details about this will be updated soon.


Tracer Fluids
----------------
"Tracer fluids"  are a field-based alternative to tracer particles.  This method creates
a user-specified number of "color field" fluids upon the initialization of a simulation
that allow the user to set a "dye" that advects along with the density field and other
color fields.  A working example of this can be found in the routines for the Rotating
Cylinder problem, and cosmology simulations can also be initialized with tracer fluids
(though the fields themselves are currently set to a constant tiny value when they are
instantiated and will need to be modified by users as needed). This can be readily
extended to other types of simulations as well. The tracer fluids can also be manipulated
during simulation runtime, as will be described shortly.

Note that the tracer fluids have no inherent physical meaning, so they can be modified on a
per-simulation and per-fluid basis to represent whatever you want.  For example, they can trace
outputs from specific star or AGN particles over time, they can follow a specific nucleosynthesis
channel from stellar feedback, they can be an instantaneous "dye" to see where a specific volume
of fluid goes, or you can simultaneously do any or all of these things with separate tracer
fluids within a single simulation. Each tracer fluid can be modified separately and they are totally
independent of each other (i.e., unlike the MultiSpecies fields there is no normalization of
density values within the tracer fluid fields), so this is a relatively flexible capability.
Note that tracer fluids are restricted to baryon fields, however - particles (such as star
particles) do NOT have "tracer field" attributes and it's not desirable for them to do so
for both performance and backward compatibility reasons, so if you want to record tracer fluid
information relating to a star particle that will have to be written into a log file of some
sort.

The Python code in the Enzo git repository (enzo-dev/src/field-injector) can be used to
modify the tracer fluid fields in an Enzo restart data output.  There are separate codes to
either (1) create entirely new fields and then manipulate them or (2) modify pre-existing
tracer fluids that were created at the beginning of the simulation or in a previous phase of
editing.  The tracer fluids only need to be added to the code a single time (i.e., at the first
restart output file where you wish to use them), but you can in principle modify the tracer
fluids using the modification code as many times as you want.  Care must be taken to not
over-write tracer fluid information that was previously injected, unless that is a desirable
behavior!

As stated above, the tracer fluids can also be manipulated during simulation runtime, and since
they are a non-normalized "color field" that gets advected along with the fluid density this is
extremely flexible.  There is a heavily annotated star particle-related example in the file
src/enzo/star_maker2.F .  This routine now takes in pointers to the tracer fluids, plus some
new parameters, and modifies the tracer fluids in a user-defined way.  This can easily be replicated
elsewhere. A comment block at the top of the file contains a lot of the necessary details, with
the rest available in comments in the star_maker and star_feedback routines in the appropriate
places to modify the fields.


Inline Halo Finder
------------------

*Source:  FOF.C*

Enzo has in-built Friends-of-Friends halo finder to identify dark matter halos, originally written by Volker Springel. Dark matter halos are identified by linking "friend" dark matter particles that lie within a specified linking length. All the output files are written in the directory FOF/.  See :ref:`inline_analysis` for more details on parameters. 



Embedded Python
---------------

*Source:  InitializePythonInterface.C*

Python can now be embedded inside Enzo, for inline analysis as well as interaction. See :ref:`embedded-python` for the details about enabling it and using it (out of date).



