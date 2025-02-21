'''
This code modifies pre-existing tracer fluid fields. As such, it is intended to be
used on Enzo datasets that already have the tracer fluid fields in them. 

WARNING: This program will modify a dataset in-place, so you should absolutely make a 
         backup copy of your original Enzo dataset before you start and then verify
         that the tracer fluids were correctly edited before you delete said backup copy.

User arguments are immediately below the imports; the modify_tracer_fluids() routine
will need to be edited by the user to make it do what they want it to.

Author:  Brian O'Shea (oshea@msu.edu), Feb. 2025
'''

import yt
import h5py
import numpy as np
import sys

user_inputs = {

    # this is the directory where the dataset that is going to be modified lives.
    # Do not leave a trailing backslash at the end of that directory, and do not
    # add any file names here.
    # NOTE: Python seems to have problems with directories that have spaces in their names, even if
    #       you put backslashes in them.
    "dataset_directory":"/Users/bwoshea/Desktop/tracer_fluid_tests/cosmotest/RD0009",

    # This is the name of the restart parameter file in the dataset directory
    # The code knows how to figure out the names of other files from that.
    "filename_stem":"RD0009",

    # Number of tracer fluid fields.  Must be at least 1 and at most 8
    # (the "at most 8" comes from the Enzo tracer fluid code).
    "NumberOfTracerFluidFields": 4,

    # This controls the level of verbosity of the outputs.  If you set it to True
    # you will get a lot of output, but it will also tell you what the code is doing.
    "DEBUG_OUTPUTS": True,  # True of False

    # if True, this will actually write the tracer fields.
    # if False, it does everything BUT write the tracer fields (dataset is unmodified)
    # It seems useful to have this feature because adding the fields is a bit tricky with
    # the various unit conversions, so you might want to do a dry run first.
    "MODIFY_FILE": True,  

    # This sets the default values of the tracer fluid density. "tiny_number" is an
    # Enzo internal value that is typically set to 1e-20.  You probably don't need
    # to modify this.
    "tiny_number": 1.0e-20
}

def modify_tracer_fields(user_inputs):
    '''
    This is the routine that actually modifies the tracer fluid fields in Enzo grid (.cpu) files.  This is 
    entirely problem-dependent, so this function will need to be modified. The example given here is meant to 
    show users how to work with grids based on cell positions.

    The general process is:

    1. Use yt to get the basic grid information (it could be done with the hierarchy file, but no need to reinvent the wheel).
    2. Loop over all grids in the dataset.  For each grid:
         * Get the grid and cell positions and calculate some useful quantities
         * Open the HDF5 file the grid lives in
         * Loop over the number of tracer fields that the user has specified.  For each tracer field:
            * Read in the tracer fluid field
            * Modify it to whatever values the user wants based on each cell's position. This will then be saved to the file.
         * Close the HDF5 file.

    VERY IMPORTANT NOTES FOR USERS:
      * This routine currently does everything in Enzo's internal coordinate system (which is 0-1 in all three spatial
        dimensions for cosmology simulations).
      * Enzo uses column-major array ordering in memory (z-dimension goes first: k, j, i) due to its solvers being
        in Fortran. Python (and numpy) use row-major array ordering in memory (x-dimension goes first: i, j, k).  So,
        any Enzo array that is read from a .cpu file into a numpy array needs to be transposed so that it is in the order
        that numpy expects. It then needs to be transposed BACK before being written to disk so that Enzo gets the arrays
        in the ordering that it expects. The code below does all of this.
    '''

    print("******** Modifying the tracer fluid fields. ********")

    # sphere center (user sets this)
    sph_cen_x = 0.5
    sph_cen_y = 0.5
    sph_cen_z = 0.5

    # sphere radius - will be multiplied by tracer field number as a test (user sets this)
    sph_dr = 0.03125

    # load up the Enzo dataset we're interested in (from user inputs)
    enzo_param_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem']
    ds = yt.load(enzo_param_file)

    # do some error checking to make sure that this dataset actually has tracer fluids in it
    if ds.parameters['UseTracerFluid'] not in ds.parameters:
        print("The parameter UseTracerFluid does not exist in the parameter file", enzo_param_file)
        print("This dataset may not have tracer fluids in it.  You need to investigate this.")
        print("Exiting.")
        sys.exit()

    # do some error checking to make sure that this dataset actually has tracer fluids in it
    if ds.parameters['NumberOfTracerFluidFields'] not in ds.parameters:
        print("The parameter NumberOfTracerFluidFields does not exist in the parameter file", enzo_param_file)
        print("This dataset may not have tracer fluids in it.  You need to investigate this.")
        print("Exiting.")
        sys.exit()

    # do some error checking to make sure that the number of tracer fluids the user has listed
    # is NO MORE THAN the number of tracer fluids in the dataset.  This assumes that the tracer
    # fields are being algorithmically modified in some way; if the user is doing something more 
    # manual this check may not actually be needed (i.e., if the user_inputs['NumberOfTracerFluidFields'] 
    # is never used, this is not relevant).
    if user_inputs['NumberOfTracerFluidFields'] > ds.parameters['NumberOfTracerFluidFields'] :
        print("The number of tracer fluid fields you specified in user_inputs is more than the number in the dataset.")
        print("Number you specified:", user_inputs['NumberOfTracerFluidFields'])
        print("Number in the dataset:", ds.parameters['NumberOfTracerFluidFields'])
        print("Exiting.")
        sys.exit()

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

        # we are going to use the numpy meshgrid functionality to create a 3D grid of x,y,z cell centers so that we can later
        # modify cell values based on their spatial positions (which seems more intuitive than array indices).  First we need
        # the cell centers along each dimension so we can fill in the meshgrid.
        xcenters_1D = ds.index.grids[i].LeftEdge.d[0] + (0.5+np.arange(ds.index.grids[i].ActiveDimensions[0]))*dx_each_dim[0]
        ycenters_1D = ds.index.grids[i].LeftEdge.d[1] + (0.5+np.arange(ds.index.grids[i].ActiveDimensions[1]))*dx_each_dim[1]
        zcenters_1D = ds.index.grids[i].LeftEdge.d[2] + (0.5+np.arange(ds.index.grids[i].ActiveDimensions[2]))*dx_each_dim[2]

        if user_inputs['DEBUG_OUTPUTS']:
            print("x,y,z centers (for mesh grid):")
            print(xcenters_1D)
            print(ycenters_1D)
            print(zcenters_1D)

        # Now we actually create the 3D mesh grid, which annoyingly can have two different indexing schemes
        # and also is returned as a list.  The 'ij' indexing scheme does things in the way that is aligned
        # with how Enzo works (after the data arrays are transposed, at least), so we use that.
        mesh_3D = np.meshgrid(xcenters_1D,ycenters_1D,zcenters_1D, indexing='ij')

        # split this out into three 3D arrays, one for each dimension.  Each of these arrays
        # now has the x, y, or z cell center for the indexed cell.
        xcenters_3D = mesh_3D[0]
        ycenters_3D = mesh_3D[1]
        zcenters_3D = mesh_3D[2]

        # calculate a grid of radius arrays using the sphere center provided by the user.
        # This will have the same dimensions as the various *centers_3D arrays.
        # This is not required in general, but is an example of something you could do!
        radius = ((xcenters_3D-sph_cen_x)**2 + (ycenters_3D-sph_cen_y)**2 + (zcenters_3D-sph_cen_z)**2 )**0.5

        if user_inputs['DEBUG_OUTPUTS']:
            print("radius min, max:", radius.min(), radius.max(), "Grid, level:", i, ds.index.grids[i].Level)

        # open up HDF5 file
        # The 'r+' option allows both reading and writing to the file.
        f = h5py.File(ds.index.grids[i].filename,'r+')

        # read density field (which should always be there) because we're going to use this
        # to set values in the tracer fluid fields.  This is NOT required if you have other
        # plans, but density is a convenient field to read because it always exists in a
        # simulation that has baryons.
        dens_name = grid_name + "/Density"

        if user_inputs['DEBUG_OUTPUTS']:
            print("density dataset name:", dens_name)

        # actually read the dataset here
        dens_dset = f[dens_name]

        # Enzo uses column-major ordering in the internal datasets, so we have to transpose datasets
        # to work with them in matplotlib.  This means we have to transpose our tracer fields
        # back to the correct ordering when we write them to the files!
        dens_dset = np.transpose(dens_dset)

        # Add tracer fluids, up to the number the user has specified
        # TracerFluid fields in Enzo are 1-indexed, which is why the range starts with 1
        # Remember that tracer fluids need to be added to ALL grids or else it will break Enzo!
        for tfnum in range(1,user_inputs['NumberOfTracerFluidFields']+1):

            # Get the name of the tracer fluid field we're going to work with
            tracer_fluid_name = 'TracerFluid' + '{:02d}'.format(tfnum)
            tf_dset_name = grid_name + '/' + tracer_fluid_name

            if user_inputs['DEBUG_OUTPUTS']:
                print("tracer fluid number, field name:", tfnum, tracer_fluid_name)
                print("tracer fluid dataset name:", tf_dset_name)

            # read this tracer fluid field, which is now a numpy array.
            this_tracer_field = f[tf_dset_name]

            # Enzo uses column-major ordering in the internal datasets, so we have to transpose
            # to work with them in matplotlib.  This means we have to transpose our tracer fields
            # back to the correct ordering when we write them to the files!
            this_tracer_field = np.transpose(this_tracer_field)

            # then set it to tiny_number (not necessary, but this sets all cells to a uniform
            # value to start up with)
            this_tracer_field[...] = user_inputs['tiny_number']

            # ***** And now we actually modify the tracer fluid in some spatially-aware way! *****

            # radius now depends on the tracer fluid number (variable tfnum), so the
            # spatial extent of each tracer fluid field is different
            myrad = sph_dr*tfnum

            # if the tracer fluid is within myrad, give it the same value as
            # the density field (this is arbitrary but convenient, you can do whatever
            # you want)
            this_tracer_field[radius<=myrad] = dens_dset[radius<=myrad]

            # We now take the tracer fluid field and transpose it back into the
            # column-major order that Enzo expects so that we can write it to disk.
            this_tracer_field = np.transpose(this_tracer_field)

            # we now copy the tracer fluid values from the numpy array back into the 
            # HDF5 file's buffers (note that the h5py docs imply this isn't necessary,
            # but the modified datasets do not seem to get set correctly otherwise).
            f[tf_dset_name][...] = this_tracer_field

            # this command forces h5py to flush dataset buffers to disk (i.e., it ensures
            # that the modified dataset actually gets written to disk)
            f.flush()

            # do a bit of housekeeping in case Python is sloppy with memory management
            # this is not always necessary, but when you have a lot of grids/arrays being created
            # sometimes weird and annoying things happen
            del this_tracer_field

        # close HDF5 file, ensuring everything gets written to disk.
        f.close()

        # memory housekeeping, as described immediately above.
        del xcenters_1D, ycenters_1D, zcenters_1D, mesh_3D, xcenters_3D, ycenters_3D, zcenters_3D, radius, dens_dset

    return

def main():
    print("Modifying Enzo data output -- .cpu files only!")
    modify_tracer_fields(user_inputs)

if __name__=="__main__":
    main()
