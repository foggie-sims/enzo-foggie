'''
This script modifies the Enzo parameter file, hierarchy file, and both
the text-based and HDF5 boundary conditions files to add the fields needed
for Enzo's shock finding algorithm, if they were not already present.
This is built on the functionality introduced by Brian O'Shea to inject
tracer fields into Enzo simulations that did not already have them.
Authors: Cassi Lochhaas, Brian O'Shea
Date created: July 29, 2026
'''


import os
import sys
import h5py
import numpy as np
import yt
import time

user_inputs = {

    # this is the directory where the dataset that is going to be modified lives.
    # Do not leave a trailing backslash at the end of that directory, and do not
    # add any file names here.
    # NOTE: Python seems to have problems with directories that have spaces in their names, even if
    #       you put backslashes in them.
    "dataset_directory":"/Users/clochhaas/Documents/Research/FOGGIE/FOGGHORN/shock_finding_test/injected_shock_fields/DD0003",

    # This is the name of the restart parameter file in the dataset directory
    # The code knows how to figure out the names of other files from that.
    "filename_stem": "data0003",

    # Shock finder method.  Must be at least 1 and at most 4
    "ShockMethod": 1,

    # Whether or not to store pre-shock fields.  Must be 0 or 1
    "StorePreShockFields": 1,

    # Temperature floor to use for shock finder.
    "ShockTemperatureFloor": 1.0,

    # When to find shocks. 0 - during EvolveLevel and writing out data, 1 - only
    # when writing out data, 2 - only during EvolveLevel
    "FindShocksOnlyOnOutput": 1,

    # This is the number of baryon fields that are in the ORIGINAL dataset.
    # If you don't know offhand, look in the dataset's .hierarchy file - each
    # grid entry has a line that says 'NumberOfBaryonFields', and it should be
    # the same for every grid entry.  Use that number.
    "NumberOfOriginalBaryonFields": 7,

    # This controls the level of verbosity of the outputs.  If you set it to True
    # you will get a lot of output, but it will also tell you what the code is doing.
    "DEBUG_OUTPUTS": True,  # True of False

    # If True, this will actually write the tracer fields.
    # If False, it does everything BUT write the tracer fields (dataset is unmodified)
    # It seems useful to have this feature because adding the fields is a bit tricky with
    # the various unit conversions, so you might want to do a dry run first.
    "MODIFY_FILES": True,

    # This has been added to fix a problem with the file system on the NASA Pleiades
    # supercomputer's file system where, if a single file is closed and then opened too
    # quickly, it will throw an HDF5 error.  If set to True, if two grids in a row are
    # stored in the same file then it will wait PLEIADES_SLEEP_TIME_SECONDS seconds
    # after it closes the file to try opening it again.
    "PLEIADES_SLEEP": False,

    # Number of seconds to wait before trying to open a file after it has been closed,
    # if PLEAIDES_SLEEP is set to True.
    # Note that this is set to what seems like a reasonable value, but it's possible it
    # could be smaller and still be fine, or may need to be larger if the file system is
    # very laggy.
    "PLEIADES_SLEEP_TIME_SECONDS": 1,

    # This sets the default values of the tracer fluid density. "tiny_number" is an
    # Enzo internal value that is typically set to 1e-20.  You probably don't need
    # to modify this.
    "tiny_number": 1.0e-20
}

################################################################################
def edit_param_file(user_inputs):
    '''
    edit_param_file

    Edits the parameter file to add shock field information. This routine may feel
    like it's doing an awful lot of error-checking, but we're trying to do something
    very invasive so we need to be a bit cautious.

    What this routine actually DOES is go through the parameter file and look for the
    line that starts with "ShockMethod" and, if it exists, set it to '1' (i.e., on).
    It also looks for the lines that start with "StorePreShockFields" and 
    "ShockTemperatureFloor" and sets those to the user-specified values. If those 
    lines do NOT exist then they are created and added to the file.  Then, we add 
    DataLabel lines for each of the new shock fluids.

    Note that we do the check for ShockMethod and StorePreShockFields because
    older simulation datasets (pre implementation of the shock finder methods in Enzo) will
    not have those lines at all.
    '''

    print("******** Editing parameter file. ********")

    # create some parameter files
    orig_param_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem']
    new_param_file = orig_param_file + ".new"
    backup_param_file = orig_param_file + ".orig"

    # print out file names
    if(user_inputs['DEBUG_OUTPUTS']):
        print(orig_param_file)
        print(new_param_file)
        print(backup_param_file)


    did_SM_exist = False   # we will check for "ShockMethod" in the parameter file.
    did_SPSF_exist = False  # we will check for "StorePreShockFields" in the parameter file.

    # does original parameter file exist?  It should.
    if(os.path.exists(orig_param_file)==True):
        print("Original parameter file exists, continuing")

    # does new parameter file exist?  It should NOT at this point.
    if(os.path.exists(new_param_file)==True):
        print("*** New parameter file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No new parameter file exists, continuing.")

    # Does backup parameter file exist?  It should NOT at this point.
    if(os.path.exists(backup_param_file)==True):
        print("*** Backup parameter file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No backup parameter file exists, continuing.")

    # Now we copy the original parameter file into a backup
    mycommand = "cp " + orig_param_file + " " + backup_param_file

    if(user_inputs['DEBUG_OUTPUTS']):
        print(mycommand)

    # do the actual copying here
    errcode = os.system(mycommand)

    # error check copying
    if errcode != 0:
        print("*** System did something fishy (parameter file), quitting.")
        sys.exit(1)

    # Now check to make sure the backup file really exists!
    if(os.path.exists(backup_param_file)==True):
        print("Backup parameter file exists at this point and SHOULD. YAY!")
    else:
        print("*** No new backup parameter file exists, and it should. Something funny is happening.  Quitting.")
        sys.exit(1)

    # And now the magic happens.  We open the original parameter file and our new one,
    # read through the original file and modify the ShockMethod-related lines as we
    # encounter them, and then finally add the extra DataLabel entries at the end of
    # the file.

    # open files
    inputfile = open(orig_param_file,"r")
    outputfile = open(new_param_file,"w")

    # loop over every line in the new file
    for thisline in inputfile:

        # split the string.
        split_line = thisline.split()

        # DO SOME ERROR-CHECKING.
        # This code strongly assumes that the simulation is 3D and that
        # the top grid is a cube, so check to make sure that's true.

        if len(split_line) > 0 and split_line[0] == 'TopGridRank':
            if int(split_line[2]) != 3:
                print("*** I expect this simulation to be 3D! Your simulation has dimensionality ", split_line[2])
                print("*** Exiting.")
                sys.exit(1)

        if len(split_line) > 0 and split_line[0] == 'TopGridDimensions':
            if (int(split_line[2]) != int(split_line[3])) or (int(split_line[2]) != int(split_line[4])):
                print("*** I expect this simulation to be a cube! Your simulation has root grid dimensions ", split_line[2], split_line[3], split_line[4])
                print("*** Exiting.")
                sys.exit(1)


        # Look for ShockMethod line.  Note that some lines have
        # length of zero, so check for that too.
        if len(split_line) > 0 and split_line[0] == 'ShockMethod':
            did_SM_exist = True
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)
            if int(split_line[2]) > 0: # do some error checking - make sure that ShockMethod is not turned on!
                print("*** Wait, this parameter file already has shock finding (ShockMethod > 0). Quitting.")
                sys.exit(1)
            else: # assuming we pass error checking, make our modified line.
                split_line[2] = str(user_inputs['ShockMethod'])
                thisline = '  '.join(split_line) + "\n"

                if(user_inputs['DEBUG_OUTPUTS']):
                    print("****new ShockMethod line:", thisline)

        # As immediately above, but for the StorePreShockFields line
        if len(split_line) > 0 and split_line[0] == 'StorePreShockFields':
            did_SPSF_exist = True
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)
            if int(split_line[2]) == 1:  # do error checking - there should not be any preshock fields in the original file.
                print("*** Wait, this parameter file already has preshock fields (StorePreShockFields = 1). Quitting.")
                sys.exit(1)
            else: # assuming we pass error checking, make our modified line
                split_line[2] = str(user_inputs['StorePreShockFields'])
                thisline = '  '.join(split_line) + "\n"
                if(user_inputs['DEBUG_OUTPUTS']):
                    print("****new StorePreShockFields line:", thisline)

        # As immediately above, but for the ShockTemperatureFloor line
        if len(split_line) > 0 and split_line[0] == 'ShockTemperatureFloor':
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)
            if float(split_line[2]) > 1.0:  # do error checking - there should not be any temperature floor in the original file.
                print("*** Wait, this parameter file already has a shock temperature floor (ShockTemperatureFloor > 1.0). Quitting.")
                sys.exit(1)
            else: # assuming we pass error checking, make our modified line
                split_line[2] = str(user_inputs['ShockTemperatureFloor'])
                thisline = '  '.join(split_line) + "\n"
                if(user_inputs['DEBUG_OUTPUTS']):
                    print("****new ShockTemperatureFloor line:", thisline)

        # As immediately above, but for the FindShocksOnlyOnOutput line
        if len(split_line) > 0 and split_line[0] == 'FindShocksOnlyOnOutput':
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)
            if int(split_line[2]) > 0:  # do error checking - there should not be any FindShocksOnlyOnOutput in the original file.
                print("*** Wait, this parameter file already has FindShocksOnlyOnOutput > 0. Quitting.")
                sys.exit(1)
            else: # assuming we pass error checking, make our modified line
                split_line[2] = str(user_inputs['FindShocksOnlyOnOutput'])
                thisline = '  '.join(split_line) + "\n"
                if(user_inputs['DEBUG_OUTPUTS']):
                    print("****new FindShocksOnlyOnOutput line:", thisline)

        # now we write either the original line or the modified line to the output file
        print(thisline, end = "", file=outputfile)


    print("\n", file=outputfile) # add a newline

    # add a 'ShockMethod' line if it didn't already exist (i.e., if the dataset is older
    # than this functionality)
    if(did_SM_exist == False):
        newline = "ShockMethod = " + str(user_inputs['ShockMethod'])
        if(user_inputs['DEBUG_OUTPUTS']):
            print("ShockMethod line did not exist, creating it. NEW LINE:")
            print(newline)
        print(newline, file=outputfile)

    # as above, but for the StorePreShockFields lines
    if(did_SPSF_exist == False):
        newline = "StorePreShockFields = " + str(user_inputs['StorePreShockFields'])
        if(user_inputs['DEBUG_OUTPUTS']):
            print("StorePreShockFields line did not exist, creating it. NEW LINE:")
            print(newline)
        print(newline, file=outputfile)

    # If we didn't have either of the two lines listed above, the ShockTemperatureFloor
    # and FindShocksOnlyOnOutput also definitely doesn't exist, so add them.
    if(did_SM_exist == False and did_SPSF_exist == False):
        newline = "ShockTemperatureFloor = " + str(user_inputs['ShockTemperatureFloor'])
        if(user_inputs['DEBUG_OUTPUTS']):
            print("ShockTemperatureFloor line almost certainly did not exist, creating it. NEW LINE:")
            print(newline)
        print(newline, file=outputfile)
        newline = "FindShocksOnlyOnOutput = " + str(user_inputs['FindShocksOnlyOnOutput'])
        if(user_inputs['DEBUG_OUTPUTS']):
            print("FindShocksOnlyOnOutput line almost certainly did not exist, creating it. NEW LINE:")
            print(newline)
        print(newline, file=outputfile)

    # now we add DataLabel entries for the shock fields.  This is not
    # necessary for Enzo, but yt (and other analysis codes) need it.
    print(user_inputs['NumberOfOriginalBaryonFields'])

    # we're assuming that the new fields are the last fields in that grid entry (which is true for Cosmology and NestedCosmology simulations as long as there are no tracer fields)
    newline = 'DataLabel[{:d}]             = '.format(user_inputs['NumberOfOriginalBaryonFields']) + "Mach\n"
    print(newline)
    print(newline, end = "", file=outputfile)

    # If storing pre-shock fields, add these data labels too
    if (user_inputs['StorePreShockFields'] == 1):
        newline = 'DataLabel[{:d}]             = '.format(1 + user_inputs['NumberOfOriginalBaryonFields']) + "PreShock_Temperature\n"
        print(newline)
        print(newline, end = "", file=outputfile)

        newline = 'DataLabel[{:d}]             = '.format(2 + user_inputs['NumberOfOriginalBaryonFields']) + "PreShock_Density\n"
        print(newline)
        print(newline, end = "", file=outputfile)

    print("\n", file=outputfile)

    # close the output files.
    inputfile.close()
    outputfile.close()

    # Remove the starting parameter file (no extra extension) and move the .new one to the original name.
    # Note that we still have a backup file with the ".orig" extension!

    mycommand = "mv " + new_param_file + " " + orig_param_file

    if(user_inputs['DEBUG_OUTPUTS']):
        print(mycommand)

    # remove our original parameter file
    os.remove(orig_param_file)

    # do the actual moving here
    errcode = os.system(mycommand)

    # error check file move
    if errcode != 0:
        print("*** system did something fishy (parameter file moving), quitting.")
        sys.exit(1)

    return


################################################################################
def edit_hierarchy_file(user_inputs):
    '''

    edit_hierarchy_file

    This file edits the hierarchy file to add shock field information.  As with the
    parameter file editing routine, it does a lot of error-checking, but that's probably
    a good thing here.

    What this actually does is modify each grid entry in two ways.  First, it updates the
    NumberOfBaryonField lines to increment it by the number of shock fields that have
    been added.  Then, it updates the FieldType line to include the typedefs (from enzo's typedefs.h file)
    for the shock fields that have been added.
    '''

    print("******** Editing hierarchy file. ********")

    orig_hierarchy_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem'] + ".hierarchy"

    new_hierarchy_file = orig_hierarchy_file + ".new"

    backup_hierarchy_file = orig_hierarchy_file + ".orig"

    if(user_inputs['DEBUG_OUTPUTS']):
        print(orig_hierarchy_file)
        print(new_hierarchy_file)
        print(backup_hierarchy_file)

    # does original hierarchy file exist?  It should.
    if(os.path.exists(orig_hierarchy_file)==True):
        print("Original hierarchy file exists, continuing")

    # does new hierarchy file exist?  It should NOT at this point.
    if(os.path.exists(new_hierarchy_file)==True):
        print("*** New hierarchy file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No new hierarchy file exists, continuing.")

    # Does backup hierarchy file exist?  It should NOT at this point.
    if(os.path.exists(backup_hierarchy_file)==True):
        print("*** Backup hierarchy file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No backup hierarchy file exists, continuing.")

    # Now we copy the original hierarchy file into a backup
    mycommand = "cp " + orig_hierarchy_file + " " + backup_hierarchy_file

    if(user_inputs['DEBUG_OUTPUTS']):
        print(mycommand)

    # do the actual copying here
    errcode = os.system(mycommand)

    # error check copying
    if errcode != 0:
        print("*** system did something fishy (hierarchy file), quitting.")
        sys.exit(1)

    # Now check to make sure the backup file really exists!
    if(os.path.exists(backup_hierarchy_file)==True):
        print("Backup hierarchy file exists at this point and SHOULD. YAY!")
    else:
        print("*** No new backup hierarchy file exists, and it should.  Something funny is happening.  Quitting.")
        sys.exit(1)

    # And now the magic happens.  We open the original hierarchy file and our new one,
    # read through the original file and modify the shock-related lines as we
    # encounter them.
    #
    # The lines that we have to modify in each grid entry are the "NumberOfBaryonFields" lines,
    # which needs to be incremented by 1 or 3 if StorePreShockFields is on, and then the FieldType line needs
    # to have an additional NumberOfBaryonFields with the typedefs for each of the tracer fluid fields
    # (as defined in Enzo's typedefs.h file).

    # open files
    inputfile = open(orig_hierarchy_file,"r")
    outputfile = open(new_hierarchy_file,"w")

    # loop over every line in the new file
    for thisline in inputfile:

        # split the string.
        split_line = thisline.split()

        # Look for NumberOfBaryonFields line.  Note that some lines have
        # length of zero, so check for that too.
        if len(split_line) > 0 and split_line[0] == 'NumberOfBaryonFields':
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)

            orig_field_num = int(split_line[2])

            if orig_field_num != user_inputs['NumberOfOriginalBaryonFields']:
                print("*** The number of baryon fields you THINK are in the original dataset are not")
                print("*** what the .hierarchy file thinks they are.  Check NumberOfOriginalBaryonFields")
                print("*** in user_inputs.py.  Exiting!")
                sys.exit(1)

            # If storing preshock density and temperature, there are 3 new fields. If not, only the Mach field is added
            if (user_inputs['StorePreShockFields'] == 1):
                new_field_num = orig_field_num + 3
            else:
                new_field_num = orig_field_num + 1
            split_line[2] = str(new_field_num)
            thisline = ' '.join(split_line) + "\n"

            if(user_inputs['DEBUG_OUTPUTS']):
                print("****new NumberOfBaryonFields line:", thisline)


        # Look for FieldType line.  Note that some lines have
        # length of zero, so check for that too.
        if len(split_line) > 0 and split_line[0] == 'FieldType':
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)

            # These numbers are assigned to the fields Mach, PreShockTemperature, and PreShockDensity in typedefs.h
            split_line.append(str(63))
            if (user_inputs['StorePreShockFields'] == 1):
                split_line.append(str(64))
                split_line.append(str(65))

            thisline = ' '.join(split_line) + "\n"

            if(user_inputs['DEBUG_OUTPUTS']):
                print("****new FieldType line:", thisline)

        # now we write either the original line or the modified line to the output file
        print(thisline, end = "", file=outputfile)


    print("\n", file=outputfile) # add a newline

    # close the output files.
    inputfile.close()
    outputfile.close()

    # Remove the starting hierarchy file (no extra extension) and move the .new one to the original name.
    # Note that we still have a backup file with the ".orig" extension!

    mycommand = "mv " + new_hierarchy_file + " " + orig_hierarchy_file

    if(user_inputs['DEBUG_OUTPUTS']):
        print(mycommand)

    # remove our original parameter file
    os.remove(orig_hierarchy_file)

    # do the actual moving here
    errcode = os.system(mycommand)

    # error check file move
    if errcode != 0:
        print("*** system did something fishy (hierarchy file moving), quitting.")
        sys.exit(1)

    return


################################################################################
def edit_boundary_files(user_inputs):
    '''

    edit_boundary_file

    This file edits both of the boundary files to add shock field information.  As with
    the other file editing routines, it does a lot of error-checking, but that's probably
    a good thing here.


    '''

    print("******** Editing boundary conditions files. ********")

    orig_boundary_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem'] + ".boundary"
    orig_HDF_boundary_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem'] + ".boundary.hdf"

    new_boundary_file = orig_boundary_file + ".new"
    new_HDF_boundary_file = orig_HDF_boundary_file + ".new"

    backup_boundary_file = orig_boundary_file + ".orig"
    backup_HDF_boundary_file = orig_HDF_boundary_file + ".orig"

    if(user_inputs['DEBUG_OUTPUTS']):
        print(orig_boundary_file)
        print(new_boundary_file)
        print(backup_boundary_file)
        print(orig_HDF_boundary_file)
        print(new_HDF_boundary_file)

    # does original boundary file exist?  It should.
    if(os.path.exists(orig_boundary_file)==True):
        print("Original boundary file exists, continuing")

    # does original HDF boundary file exist?  It should.
    if(os.path.exists(orig_HDF_boundary_file)==True):
        print("Original HDF boundary file exists, continuing")

    # does new boundary file exist?  It should NOT at this point.
    if(os.path.exists(new_boundary_file)==True):
        print("*** New boundary file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No new boundary file exists, continuing.")

    # does new HDF boundary file exist?  It should NOT at this point.
    if(os.path.exists(new_HDF_boundary_file)==True):
        print("*** New HDF boundary file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No new HDF boundary file exists, continuing.")

    # Does backup boundary file exist?  It should NOT at this point.
    if(os.path.exists(backup_boundary_file)==True):
        print("*** Backup boundary file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No backup boundary file exists, continuing.")

    # Does backup HDF boundary file exist?  It should NOT at this point.
    if(os.path.exists(backup_HDF_boundary_file)==True):
        print("*** Backup HDF boundary file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No backup HDF boundary file exists, continuing.")

    # Now we copy the original boundary file into a backup.
    mycommand = "cp " + orig_boundary_file + " " + backup_boundary_file

    if(user_inputs['DEBUG_OUTPUTS']):
        print(mycommand)

    # do the actual copying here
    errcode = os.system(mycommand)

    # error check copying
    if errcode != 0:
        print("*** system did something fishy (boundary file), quitting.")
        sys.exit(1)

    # Now check to make sure the backup file really exists!
    if(os.path.exists(backup_boundary_file)==True):
        print("Backup boundary file exists at this point and SHOULD. YAY!")
    else:
        print("*** No new backup boundary file exists, and it should.  Something funny is happening.  Quitting.")
        sys.exit(1)

    # Now we copy the original HDF boundary file into a backup.
    mycommand = "cp " + orig_HDF_boundary_file + " " + backup_HDF_boundary_file

    if(user_inputs['DEBUG_OUTPUTS']):
        print(mycommand)

    # do the actual copying here
    errcode = os.system(mycommand)

    # error check copying
    if errcode != 0:
        print("*** system did something fishy (HDF boundary file), quitting.")
        sys.exit(1)

    # Now check to make sure the backup HDF file really exists!
    if(os.path.exists(backup_HDF_boundary_file)==True):
        print("Backup HDF boundary file exists at this point and SHOULD. YAY!")
    else:
        print("*** No new backup HDF boundary file exists, and it should.  Something funny is happening.  Quitting.")
        sys.exit(1)

    # we're going to need the size of the boundary conditions for when we create our new HDF5
    # boundary files.
    bcx = bcy = bcz = -1


    # And now the magic happens.  We open the original boundary file and our new one,
    # read through the original file and modify the shock-related lines as we
    # encounter them.
    #
    # The lines that we have to modify in each grid entry are the "NumberOfBaryonFields" lines,
    # which needs to be incremented by 1 or 3 depending on if StorePreShockFields is 0 or 1, and then the FieldType line needs
    # to have an additional NumberOfBaryonFields with the typedefs for each of the shock fields
    # (as defined in Enzo's typedefs.h file).

    # open files
    inputfile = open(orig_boundary_file,"r")
    outputfile = open(new_boundary_file,"w")

    # loop over every line in the new file
    for thisline in inputfile:

        # split the string.
        split_line = thisline.split()

        # Look for NumberOfBaryonFields line.  Note that some lines have
        # length of zero, so check for that too.
        if len(split_line) > 0 and split_line[0] == 'NumberOfBaryonFields':
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)

            orig_field_num = int(split_line[2])
            if (user_inputs['StorePreShockFields'] == 1):
                new_field_num = orig_field_num + 3
            else:
                new_field_num = orig_field_num + 1
            split_line[2] = str(new_field_num)
            thisline = ' '.join(split_line) + "\n"

            if(user_inputs['DEBUG_OUTPUTS']):
                print("****new NumberOfBaryonFields line:", thisline)


        # Look for FieldType line.  Note that some lines have
        # length of zero, so check for that too.
        if len(split_line) > 0 and split_line[0] == 'BoundaryFieldType':
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)

            split_line.append(str(63))
            if (user_inputs['StorePreShockFields'] == 1):
                split_line.append(str(64))
                split_line.append(str(65))

            thisline = ' '.join(split_line) + "\n"

            if(user_inputs['DEBUG_OUTPUTS']):
                print("****new FieldType line:", thisline)


        # Look for BoundaryDimension line.  Note that some lines have
        # length of zero, so check for that too.
        if len(split_line) > 0 and split_line[0] == 'BoundaryDimension':
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)

            bcx = int(split_line[2])
            bcy = int(split_line[3])
            bcz = int(split_line[4])

            if bcx != bcy or bcx != bcz:
                print("*** this box is not a cube, lots of assumptions break.  Exiting.")
                sys.exit(1)

            if(user_inputs['DEBUG_OUTPUTS']):
                print("boundary dimensions: ", bcx, bcy, bcz)


        # now we write either the original line or the modified line to the output file
        print(thisline, end = "", file=outputfile)


    print("\n", file=outputfile) # add a newline

    # close the output files.
    inputfile.close()
    outputfile.close()

    # Remove the starting boundary file (no extra extension) and move the .new one to the original name.
    # Note that we still have a backup file with the ".orig" extension!

    mycommand = "mv " + new_boundary_file + " " + orig_boundary_file

    if(user_inputs['DEBUG_OUTPUTS']):
        print(mycommand)

    # remove our original boundary file
    os.remove(orig_boundary_file)

    # do the actual moving here
    errcode = os.system(mycommand)

    # error check file move
    if errcode != 0:
        print("*** system did something fishy (boundary file moving), quitting.")
        sys.exit(1)

    # Now on to the HDF5 file! There's no point in copying this in a streaming way like we did the text files.
    # We're just going to go ahead and create a whole new file with h5py.

    # need total number of baryon fields
    if (user_inputs['StorePreShockFields'] == 1):
        total_num_baryon_fields = user_inputs['NumberOfOriginalBaryonFields'] + 3
    else:
        total_num_baryon_fields = user_inputs['NumberOfOriginalBaryonFields'] + 1

    # The arrays are for the 2D faces of the root grid, including ghost zones, and for each dimension there are
    # two faces and total_num_baryon_fields entries.
    array_size = bcx**2 * 2 * total_num_baryon_fields

    if(user_inputs['DEBUG_OUTPUTS']):
        print("boundary array sizes are: ", array_size)
        print("estimated file size: ", array_size * 4 * 6)  # 4 comes from 4 bytes/float; 6 comes from 6 arrays total.

    # the boundary dimension type arrays are ones
    ones_array = np.ones(array_size, dtype='float32')

    # the boundary VALUE type is zero.
    zeros_array = np.zeros(array_size, dtype='float32')

    # open file
    f = h5py.File(new_HDF_boundary_file,'w')

    # open first boundary dimension dataset and fill it with ones (and the specific type of floats that Enzo wants here)
    # Then add a bunch of attributes that Enzo expects to see in this file, based on the values calculated above.
    # All of the '>f4' and '>i4' stuff is making sure that the floats and ints are 32-bit big-endian, which is what Enzo needs.
    dset = f.create_dataset("BoundaryDimensionType.0",data=ones_array.astype('>f4'))
    dset.attrs.create(name='BoundaryDimension', data=np.array((bcx, bcx, 0), dtype='int32'), shape=(3,), dtype='>i4')
    dset.attrs.create(name='BoundaryRank', data=np.array((3), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='Index', data=np.array((2), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='NumberOfBaryonFields', data=np.array((total_num_baryon_fields), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='size', data=np.array((bcx*bcx), dtype='int32'), shape=(1,), dtype='>i4')

    # do the same for the second boundary dimension
    dset = f.create_dataset("BoundaryDimensionType.1",data=ones_array.astype('>f4'))
    dset.attrs.create(name='BoundaryDimension', data=np.array((bcx, bcx, 0), dtype='int32'), shape=(3,), dtype='>i4')
    dset.attrs.create(name='BoundaryRank', data=np.array((3), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='Index', data=np.array((2), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='NumberOfBaryonFields', data=np.array((total_num_baryon_fields), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='size', data=np.array((bcx*bcx), dtype='int32'), shape=(1,), dtype='>i4')

    # do the same for the third boundary dimension
    dset = f.create_dataset("BoundaryDimensionType.2",data=ones_array.astype('>f4'))
    dset.attrs.create(name='BoundaryDimension', data=np.array((bcx, bcx, 0), dtype='int32'), shape=(3,), dtype='>i4')
    dset.attrs.create(name='BoundaryRank', data=np.array((3), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='Index', data=np.array((2), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='NumberOfBaryonFields', data=np.array((total_num_baryon_fields), dtype='int32'), shape=(1,), dtype='>i4')
    dset.attrs.create(name='size', data=np.array((bcx*bcx), dtype='int32'), shape=(1,), dtype='>i4')

    # now set the actual values for the boundary dimension.  This is all zeros, and Enzo doesn't want any attributes.
    f.create_dataset("BoundaryDimensionValue.0",data=zeros_array.astype('>f4'))
    f.create_dataset("BoundaryDimensionValue.1",data=zeros_array.astype('>f4'))
    f.create_dataset("BoundaryDimensionValue.2",data=zeros_array.astype('>f4'))

    f.close()

    # Remove the starting HDF5 boundary file (no extra extension) and move the .new one to the original name.
    # Note that we still have a backup file with the ".orig" extension!

    mycommand = "mv " + new_HDF_boundary_file + " " + orig_HDF_boundary_file

    if(user_inputs['DEBUG_OUTPUTS']):
        print(mycommand)

    # remove our original HDF5 boundary condition file
    os.remove(orig_HDF_boundary_file)

    # do the actual moving here
    errcode = os.system(mycommand)

    # error check file move
    if errcode != 0:
        print("*** system did something fishy (boundary HDF5 file moving), quitting.")
        sys.exit(1)

    return

def modify_grid_files(user_inputs):
    '''
    This is the routine that actually modifies the grid files.

    The general flow is:

    1. Use yt to get the basic grid information (it could be done with the hierarchy file, but no need to reinvent the wheel).
    2. Loop over all grids in the dataset.  For each grid:
         * Get the grid and cell positions and calculate some useful quantities
         * Open the HDF5 file the grid lives in
         * Open up the density dataset (which always must exist)
         * Loop over the shock fields that the user has specified.  For each shock field:
            * Create an array that's the same size/shape/precision as the density array
            * Set it to tiny_number
         * Close the HDF5 file, which will save the datasets.

    VERY IMPORTANT NOTES FOR USERS:
      * The shock fields must be added to ALL grids, not just grids where you want to find shocks.  This
        is because Enzo requires that all grids have the same set of baryon fields.
      * This routine currently does everything in Enzo's internal coordinate system (which is 0-1 in all three spatial
        dimensions for cosmology simulations).
      * Enzo uses column-major array ordering in memory (z-dimension goes first: k, j, i) due to its solvers being
        in Fortran. Python (and numpy) use row-major array ordering in memory (x-dimension goes first: i, j, k).  So,
        any Enzo array that is read from a .cpu file into a numpy array needs to be transposed so that it is in the order
        that numpy expects. It then needs to be transposed BACK before being written to disk so that Enzo gets the arrays
        in the ordering that it expects. The code below does all of this.
    '''

    print("******** Modifying the grid files. ********")

    # sphere center (user sets this)
    sph_cen_x = 0.52587891
    sph_cen_y = 0.51708984
    sph_cen_z = 0.48486328

    # sphere radius - will be multiplied by tracer field number as a test (user sets this)
    sph_dr = 0.015625

    # load up the Enzo dataset we're interested in (from user inputs)
    enzo_param_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem']
    ds = yt.load(enzo_param_file)

    # keeps track of the last file that was opened so that we can hold off
    # on re-opening it if necessary.
    last_file_opened = None

    # Loop over all of the grids and do things.
    # Note that we have to add the shock fields to all of the grids, even if you only want to
    # find shocks in some subvolume of the simulations.  This is because Enzo expects that all grids
    # will have the same baryon fields.
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

        # If the last grid was in the same file as the current grid we're working on (i.e., if the last
        # file we worked on is the same as this file) then some file systems need a bit of time to
        # realized that the file was recently closed so that HDF5/h5py doesn't throw an error.  So,
        # if this feature is turned on and the last file is the same as this file, then wait for
        # a user-specified number of seconds.
        if last_file_opened == ds.index.grids[i].filename and user_inputs['PLEIADES_SLEEP'] == True:
            if user_inputs['DEBUG_OUTPUTS']:
                print("************  Last file is the same as this file.  Waiting for this many seconds:", user_inputs['PLEIADES_SLEEP_TIME_SECONDS'])
            time.sleep(user_inputs['PLEIADES_SLEEP_TIME_SECONDS'])

        # open up HDF5 file
        # The 'r+' option allows both reading and writing to the file.
        f = h5py.File(ds.index.grids[i].filename,'r+')

        # read density field (which should always be there) to get dataset dimensions
        # (the tracer fluid fields must be the same size as the other baryon fields, so we're
        # just going to be creatively lazy here)
        dens_name = grid_name + "/Density"

        if user_inputs['DEBUG_OUTPUTS']:
            print("density dataset name:", dens_name)

        # actually read the dataset here
        dens_dset = f[dens_name]

        # Enzo uses column-major ordering in the internal datasets, so we have to transpose datasets
        # to work with them in matplotlib.  This means we have to transpose our tracer fields
        # back to the correct ordering when we write them to the files!
        #dens_dset = np.transpose(dens_dset)

        # Add shock fields, either 3 if StorePreShockFields is on or 1 if not
        # Remember that new fields need to be added to ALL grids or else it will break Enzo!
        # This will create a tracer fluid dataset name that is aligned with what Enzo expects

        Mach_name = 'Mach'
        Mach_dset_name = grid_name + '/' + Mach_name

        if user_inputs['DEBUG_OUTPUTS']:
            print("field name:", Mach_name)
            print("dataset name:", Mach_dset_name)

        # first create a field of zeros
        Mach_field = np.zeros_like(dens_dset)

        # then set it to tiny_number
        Mach_field[...] = user_inputs['tiny_number']

        # Then actually write the dataset, if the user wants you to!
        if user_inputs['MODIFY_FILES']:
            if user_inputs['DEBUG_OUTPUTS']:
                print("writing dataset", Mach_dset_name, "in grid", grid_name)
            f.create_dataset(Mach_dset_name,data=Mach_field)

        if (user_inputs['StorePreShockFields'] == 1):
            Temperature_name = 'PreShock_Temperature'
            Temperature_dset_name = grid_name + '/' + Temperature_name

            if user_inputs['DEBUG_OUTPUTS']:
                print("field name:", Temperature_name)
                print("dataset name:", Temperature_dset_name)

            # first create a field of zeros
            Temperature_field = np.zeros_like(dens_dset)

            # then set it to tiny_number
            Temperature_field[...] = user_inputs['tiny_number']

            # Then actually write the dataset, if the user wants you to!
            if user_inputs['MODIFY_FILES']:
                if user_inputs['DEBUG_OUTPUTS']:
                    print("writing dataset", Temperature_dset_name, "in grid", grid_name)
                f.create_dataset(Temperature_dset_name,data=Temperature_field)

            # do a bit of housekeeping in case Python is sloppy with memory management
            # this is not always necessary, but when you have a lot of grids/arrays being created
            # sometimes weird and annoying things happen
            del Temperature_field

            Density_name = 'PreShock_Density'
            Density_dset_name = grid_name + '/' + Density_name

            if user_inputs['DEBUG_OUTPUTS']:
                print("field name:", Density_name)
                print("dataset name:", Density_dset_name)

            # first create a field of zeros
            Density_field = np.zeros_like(dens_dset)

            # then set it to tiny_number
            Density_field[...] = user_inputs['tiny_number']

            # Then actually write the dataset, if the user wants you to!
            if user_inputs['MODIFY_FILES']:
                if user_inputs['DEBUG_OUTPUTS']:
                    print("writing dataset", Density_dset_name, "in grid", grid_name)
                f.create_dataset(Density_dset_name,data=Density_field)

            # do a bit of housekeeping in case Python is sloppy with memory management
            # this is not always necessary, but when you have a lot of grids/arrays being created
            # sometimes weird and annoying things happen
            del Density_field

        # memory housekeeping, as described immediately above.
        del dens_dset

        # close HDF5 file, ensuring everything gets written to disk.
        f.close()

        last_file_opened = ds.index.grids[i].filename

    return

if __name__ == '__main__':

    if user_inputs['DEBUG_OUTPUTS']:
        print("*"*40)
        print("input user dictionary:\n")
        print(user_inputs)
        print("*"*40,"\n")

    # Creates the new parameter file (with shock field contents in it)
    if user_inputs['MODIFY_FILES'] == True:
        edit_param_file(user_inputs)
    else:
        print("Skipping edit_param_file because this is a dry run (MODIFY_FILES = False)\n")

    # Creates the new hierarchy file (with shock field contents in it)
    if user_inputs['MODIFY_FILES'] == True:
        edit_hierarchy_file(user_inputs)
    else:
        print("Skipping edit_hierarchy_file because this is a dry run (MODIFY_FILES = False)\n")

    # Creates the two new boundary conditions files (with shock field contents in them)
    if user_inputs['MODIFY_FILES'] == True:
        edit_boundary_files(user_inputs)
    else:
        print("Skipping edit_boundary_files because this is a dry run (MODIFY_FILES = False)\n")

    # Modifies the existing grid files to add shock fields. Note that we don't use the
    # same logic as the prior function calls (re: MODIFY_FILES) because we often will want to
    # actually go through the logic of modifying the grid files without doing so, but the other
    # files are much more straightforward so that's not as big of a concern.
    modify_grid_files(user_inputs)