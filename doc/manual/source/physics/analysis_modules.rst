.. _analysis_modules:

Analysis Modules
================


Tracer Particles
----------------
Enzo has extensive analysis methods to track tracer particles. Details about this will be updated soon.


Tracer Fluids
----------------
"Tracer fluids"  are a field-based alternative to tracer particles.  This method creates a user-specified number of "color field" fluids upon the initialization of a simulation that allow the user to set a "dye" that advects along with the density field and other color fields.  A working example of this can be found in the routines for the Rotating Cylinder problem, and cosmology simulations can also be initialized with tracer fluids (though the fields themselves are currently set to a constant tiny value, and will need to be modified by users as needed). This can be readily extended to other types of calculations as well. The Python code in the Enzo git repository (enzo-dev/src/field-injector) can be used to modify the tracer fields during simulation runtime.


Inline Halo Finder
------------------

*Source:  FOF.C*

Enzo has in-built Friends-of-Friends halo finder to identify dark matter halos, originally written by Volker Springel. Dark matter halos are identified by linking "friend" dark matter particles that lie within a specified linking length. All the output files are written in the directory FOF/.  See :ref:`inline_analysis` for more details on parameters. 



Embedded Python
---------------

*Source:  InitializePythonInterface.C*

Python can now be embedded inside Enzo, for inline analysis as well as interaction. See :ref:`embedded-python` for the details about enabling it and using it (out of date).



