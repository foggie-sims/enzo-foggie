# Enzo field injector

Written by Brian O'Shea, Feb. 2025

## Purpose 

The goal of this code is to add tracer fluids to a three-dimensional Enzo dataset that does not currently have said tracer fluid fields in it. It does so by modifying basically all of the files in a standard Enzo data output (parameter file, hierarchy file, both boundary conditions files, and the actual .cpuNNNN files where the simulation data is kept).

**WARNING: Make sure to create a backup copy of your simulation dataset before you modify it - this type of modification will irrevocably change parts of the data output being modified!**


## The code 

**Running the code:** Type `python field_injector.py` and then step back and watch the fireworks.

The files included in this code are:

`field_injector.py` - the driver code, which calls routines that edit the various Enzo data output files.

`user_inputs.py` - This includes both a dictionary of user inputs that are used throughout the code as well as a routine called `modify_grid_files` that will inevitably need to be edited by users so that it sets the tracer fluid fields to the values they want.

`mod_routines.py` - This includes the routines that modify the parameter file, hierarchy file, and both the ASCII and binary boundary conditions files.  This should not need to be edited by the user.


## Notes/Caveats/known problems

* When adding the tracer field this code goes through grids in numerical order, which means that the code is accessing the binary (.cpuNNNN) files in effectively random order.  This could cause a performance issue; if so, the `modify_grid_files` routine can be modified so that it goes through the grids in a different order so that a single .cpuNNNN file is being edited at a time.
* This code should be easy to parallelize if need be - the .cpuNNNN files can be edited completely independently of each other.
* The code only runs on 3D, cubic datasets at present. This limitation is due to assumptions baked into the code in various places, and is relatively easy to fix if necessary.