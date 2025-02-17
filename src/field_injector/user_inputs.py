'''
This file contains user inputs, which are stored in a dictionary for the sake of flexibility.
Definitions of specific parameters are right above the dictionary entries.
'''

user_inputs = {

    # this is the directory where the dataset that is going to be modified lives.
    # Do not leave a trailing backslash at the end of that directory, and do not
    # add any file names here.
    "dataset_directory":"/Users/bwoshea/Desktop/tests/DD0000",

    # This is the name of the restart parameter file in the dataset directory
    # The code knows how to figure out the names of other files from that.
    "filename_stem":"DD0000",

    # Number of tracer fluid fields.  Must be at least 1 and at most 8
    # (the "at most 8" comes from the Enzo tracer fluid code).
    "NumberOfTracerFluidFields": 4,

    # This is the number of baryon fields that are in the ORIGINAL dataset.
    # If you don't know offhand, look in the dataset's .hierarchy file - each
    # grid entry has a line that says 'NumberOfBaryonFields', and it should be
    # the same for every grid entry.  Use that number.
    "NumberOfOriginalBaryonFields": 6, 

    # This controls the level of verbosity of the outputs.  If you set it to True
    # you will get a lot of output, but it will also tell you what the code is doing.
    "DEBUG_OUTPUTS": True,  # True of False

    # This sets the default values of the tracer fluid density. "tiny_number" is an
    # Enzo internal value that is typically set to 1e-20.  You probably don't need
    # to modify this.
    "tiny_number": 1.0e-20
}


import yt
import h5py
import numpy as np


def modify_grid_files(user_inputs):
    '''
    This is the routine that actually modifies the grid files.  This is much more annoying than any other part of this
    particular code, because what it does is entirely problem-dependent. The general flow is:

    1. Use yt to get the basic grid information (it could be done with the hierarchy file, but no need to reinvent the wheel).
    2. Loop over all grids in the dataset.  For each grid:
         * Get the grid and cell positions and calculate some useful quantities
         * Open the HDF5 file the grid lives in
         * Open up the density dataset (which always must exist)
         * Loop over the number of tracer fields that the user has specified.  For each tracer field:
            * Create an array that's the same size/shape/precision as the density array
            * Set it to tiny_number first
            * Modify it to whatever values the user wants based on each cell's position. This will automatically be
              saved in the file upon closing.
         * Close the HDF5 file, which will save the datasets.

    NOTES FOR USERS:
      * The tracer fluid fields must be added to ALL grids, not just grids where you want to trace something.  This
        is because Enzo requires that all grids have the same set of baryon fields.  Just set values in grids that 
        you aren't interested in to some small value.
      * This routine currently does everything in Enzo's internal coordinate system (which is 0-1 in all three spatial
        dimensions for cosmology simulations).
    '''

    print("******** Modifying the grid files. ********")

    MODIFY_FILE = True  # if True, this will actually write the tracer fields.
                        # if False, it does everything BUT write the tracer fields (dataset is unmodified)
                        # It seems useful to have this feature because adding the fields is a bit tricky with
                        # the various unit conversions, so you might want to do a dry run first.

    # sphere center 
    sph_cen_x = 0.5
    sph_cen_y = 0.5
    sph_cen_z = 0.5

    # sphere radius - will be multiplied by tracer field number as a test
    sph_dr = 0.03125

    # load up the dataset we're interested in (from user inputs)
    enzo_param_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem']
    ds = yt.load(enzo_param_file)

    # Loop over all of the grids and do things.
    # Note that we have to add the tracer fluid fields to all of the grids, even if you only want to
    # trace fluids in some subvolume of the simulations.  This is because Enzo expects that all grids
    # will have the same baryon fields. Just set values of the tracer field to a very small value in
    # uninteresting regions.
    for i in range(len(ds.index.grids)):

        if user_inputs['DEBUG_OUTPUTS']:
            print("working on grid", i, "in file", ds.index.grids[i].filename)

        # grid numbers (in the grid names) are 1-indexed, not zero-indexed
        # also the name is zero-padded to have 8 digits total.  If this
        # is ever changed in an Enzo dataset (to have more padding, for example)
        # this is immediately going to crash.
        grid_name = 'Grid' + '{:08d}'.format(i+1)

        # print out some useful information about this grid
        if user_inputs['DEBUG_OUTPUTS']:
            print("grid name is     ", grid_name)  # the actual grid name (that we created)
            print("Grid left edge:  ", ds.index.grids[i].LeftEdge)   # yt-provided grid left edge
            print("Grid right edge: ", ds.index.grids[i].RightEdge)  # yt-provided grid right edge
            print("Grid dimensions: ", ds.index.grids[i].ActiveDimensions)  # yt-provided grid active dimensions (no ghost zones)
            print("Grid level:      ", ds.index.grids[i].Level)      # yt-provided grid level

        # figure out grid edge size along each dimension - they should ALWAYS be identical.
        dx_each_dim = (ds.index.grids[i].RightEdge.d - ds.index.grids[i].LeftEdge.d)/ds.index.grids[i].ActiveDimensions

        if user_inputs['DEBUG_OUTPUTS']:
            print("Grid dx (per dim):", dx_each_dim)

        # we are going to use the numpy meshgrid functionality to create a 3D grid of x,y,z cell centers.  First we need
        # the cell centers along each dimension so we can fill in the meshgrid.
        xcenters_1D = ds.index.grids[i].LeftEdge.d[0] + (0.5+np.arange(ds.index.grids[i].ActiveDimensions[0]))*dx_each_dim[0]
        ycenters_1D = ds.index.grids[i].LeftEdge.d[1] + (0.5+np.arange(ds.index.grids[i].ActiveDimensions[1]))*dx_each_dim[1]
        zcenters_1D = ds.index.grids[i].LeftEdge.d[2] + (0.5+np.arange(ds.index.grids[i].ActiveDimensions[2]))*dx_each_dim[2]

        if user_inputs['DEBUG_OUTPUTS']:
            print("x,y,z centers (for mesh grid):")
            print(xcenters_1D)
            print(ycenters_1D)
            print(zcenters_1D)

        # now we actually create the 3D mesh grid, which annoyingly can have two different indexing schemes
        # and also is returned as a list.
        mesh_3D = np.meshgrid(xcenters_1D,ycenters_1D,zcenters_1D, indexing='xy')  # indexing can be 'xy' or 'ij'

        # split this out into three 3D arrays, one for each dimension
        xcenters_3D = mesh_3D[0]
        ycenters_3D = mesh_3D[1]
        zcenters_3D = mesh_3D[2]

        # calculate a grid of radius arrays using the sphere center provided by the user.
        # This will have the same dimensions as the various *centers_3D arrays.
        radius = ((xcenters_3D-sph_cen_x)**2 + (ycenters_3D-sph_cen_y)**2 + (zcenters_3D-sph_cen_z)**2 )**0.5

        if user_inputs['DEBUG_OUTPUTS']:
            print("radius min, max:", radius.min(), radius.max(), "Grid, level:", i, ds.index.grids[i].Level)

        # open up HDF5 file
        # The 'r+' option should let me open the file.
        f = h5py.File(ds.index.grids[i].filename,'r+')

        # read density field (which should always be there) to get dataset dimensions
        # (the tracer fluid fields must be the same size as the other baryon fields, so we're
        # just going to be creatively lazy here)
        dens_name = grid_name + "/Density"

        if user_inputs['DEBUG_OUTPUTS']:
            print("density dataset name:", dens_name)

        # actually read the dataset here
        dens_dset = f[dens_name]

        # Add tracer fluids, up to the number the user has specified
        # TracerFluid fields in Enzo are 1-indexed, which is why the range starts with 1
        # Remember that tracer fluids need to be added to ALL grids or else it will break Enzo!
        for tfnum in range(1,user_inputs['NumberOfTracerFluidFields']+1):

            # This will create a tracer fluid dataset name that is aligned with what Enzo expects
            tracer_fluid_name = 'TracerFluid' + '{:02d}'.format(tfnum)
            tf_dset_name = grid_name + '/' + tracer_fluid_name

            if user_inputs['DEBUG_OUTPUTS']:
                print("tracer fluid number, field name:", tfnum, tracer_fluid_name)
                print("tracer fluid dataset name:", tf_dset_name)

            # first create a tracer field of zeros
            this_tracer_field = np.zeros_like(dens_dset)

            # then set it to tiny_number (not necessary, but it's consistent with how Enzo
            # generates initial uniform grids)
            this_tracer_field[...] = user_inputs['tiny_number']

            # ***** And now we actually modify the tracer fluid in some spatially-aware way! *****

            # radius now depends on the tracer fluid number, so the
            # spatial extent of each tracer fluid field is different
            myrad = sph_dr*tfnum

            # if the tracer fluid is within myrad, give it the same value as
            # the density field (this is arbitrary but convenient, you can do whatever
            # you want)
            this_tracer_field[radius<=myrad] = dens_dset[radius<=myrad]

            # Then actually write the dataset, if the user wants you to!
            if MODIFY_FILE:
                if user_inputs['DEBUG_OUTPUTS']:
                    print("writing dataset", tf_dset_name, "for field", tfnum, "in grid", grid_name)
                f.create_dataset(tf_dset_name,data=this_tracer_field)

            # do a bit of housekeeping in case Python is sloppy with memory management
            # this is not always necessary, but when you have a lot of grids/arrays being created
            # sometimes weird and annoying things happen
            del this_tracer_field

        # memory housekeeping, as described immediately above.
        del xcenters_1D, ycenters_1D, zcenters_1D, mesh_3D, xcenters_3D, ycenters_3D, zcenters_3D, radius

        # close HDF5 file, ensuring everything gets saved.
        f.close()

    return
