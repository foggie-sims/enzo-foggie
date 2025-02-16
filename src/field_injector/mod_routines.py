import os
import sys

################################################################################
def edit_param_file(user_inputs):
    '''
    edit_param_file

    Edits the parameter file to add tracer fluid information. This routine may feel
    like it's doing an awful lot of error-checking, but we're trying to do something
    very invasive so we need to be a bit cautious.

    What this routine actually DOES is go through the parameter file and look for the
    line that starts with "UseTracerFluid" and, if it exists, set it to '1' (i.e., on).
    It also looks for the line that starts with "NumberOfTracerFluidFields" and sets
    that to the user-specified values. If those lines to NOT exist then they are created
    and added to the file.  Then, we add DataLabel lines for each of the new tracer fluids.

    Note that we do the check for UseTracerField and NumberOfTracerFluidFields because 
    older simulation datasets (pre implementation of the tracer fluid methods in Enzo) will
    not have those lines at all.
    '''

    # create some parameter files
    orig_param_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem']
    new_param_file = orig_param_file + ".new"
    backup_param_file = orig_param_file + ".orig"

    # print out file names
    if(user_inputs['DEBUG_OUTPUTS']):
        print(orig_param_file)
        print(new_param_file)
        print(backup_param_file)


    did_UTF_exist = False   # we will check for "UseTracerFluid" in the parameter file.
    did_NOTFF_exist = False  # we will check for "NumberOfTracerFluidFields" in the parameter file.
        
    # does original parameter file exist?  It should.
    if(os.path.exists(orig_param_file)==True):
        print("Original parameter file exists, continuing")

    # does new parameter file exist?  It should NOT at this point.
    if(os.path.exists(new_param_file)==True):
        print("New parameter file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No new parameter file exists, continuing.")

    # Does backup parameter file exist?  It should NOT at this point.
    if(os.path.exists(backup_param_file)==True):
        print("Backup parameter file exists at this point and shouldn't. You may need to delete something. Exiting.")
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
        print("system did something fishy (parameter file), quitting.")
        sys.exit(1)

    # Now check to make sure the backup file really exists!
    if(os.path.exists(backup_param_file)==True):
        print("Backup parameter file exists at this point and SHOULD. YAY!")
    else:
        print("No new backup parameter file exists, and it should. Something funny is happening.  Quitting.")
        sys.exit(1)

    # And now the magic happens.  We open the original parameter file and our new one,
    # read through the original file and modify the TracerFluid-related lines as we
    # encounter them, and then finally add the extra DataLabel entries at the end of
    # the file.

    # open files
    inputfile = open(orig_param_file,"r")
    outputfile = open(new_param_file,"w")

    # loop over every line in the new file
    for thisline in inputfile:

        # split the string.
        split_line = thisline.split()

        # Look for UseTracerFluid line.  Note that some lines have
        # length of zero, so check for that too.
        if len(split_line) > 0 and split_line[0] == 'UseTracerFluid':
            did_UTF_exist = True
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)
            if int(split_line[2]) == 1: # do some error checking - make sure that UseTracerFluid is not turned on!
                print("Wait, this parameter file already has tracer fluids (UseTracerFluid = 1). Quitting.")
                sys.exit(1)
            else: # assuming we pass error checking, make our modified line.
                split_line[2] = '1'
                thisline = '  '.join(split_line) + "\n"

                if(user_inputs['DEBUG_OUTPUTS']):
                    print("****new UseTracerFluid line:", thisline)

        # As immediately above, but for the NumberOfTracerFluidFields line
        if len(split_line) > 0 and split_line[0] == 'NumberOfTracerFluidFields':
            did_NOTFF_exist = True
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)
            if int(split_line[2]) > 0:  # do error checking - there should not be any tracer fluids in the original file.
                print("Wait, this parameter file already has tracer fluids (NumberOfTracerFluidFields > 0). Quitting.")
                sys.exit(1)
            else: # assuming we pass error checking, make our modified line
                split_line[2] = str(user_inputs['NumberOfTracerFluidFields'])
                thisline = '  '.join(split_line) + "\n"
                if(user_inputs['DEBUG_OUTPUTS']):
                    print("****new NumberOfTracerFluidFields line:", thisline)

        # now we write either the original line or the modified line to the output file
        print(thisline, end = "", file=outputfile)


    print("\n", file=outputfile) # add a newline

    # add a 'UserTracerFluid' line if it didn't already exist (i.e., if the dataset is older
    # than this functionality)
    if(did_UTF_exist == False):
        newline = "UseTracerFluid = 1"
        if(user_inputs['DEBUG_OUTPUTS']):
            print("UserTracerFluid line did not exist, creating it. NEW LINE:")
            print(newline)
        print(newline, file=outputfile)

    # as above, but for the NumberOfTracerFluidFields lines
    if(did_NOTFF_exist == False):
        newline = "NumberOfTracerFluidFields = " + str(user_inputs['NumberOfTracerFluidFields'])
        if(user_inputs['DEBUG_OUTPUTS']):
            print("NumberOfTracerFluidFields line did not exist, creating it. NEW LINE:")
            print(newline)
        print(newline, file=outputfile)

    # If we didn't have either of the two lines listed above, the SetTracerFluidFieldsOnStart
    # also definitely doesn't exist, so add it.
    if(did_UTF_exist == False and did_NOTFF_exist == False):
        newline = "SetTracerFluidFieldsOnStart = 0"
        if(user_inputs['DEBUG_OUTPUTS']):
            print("SetTracerFluidFieldsOnStart line almost certainly did not exist, creating it. NEW LINE:")
            print(newline)
        print(newline, file=outputfile)

    # now we add DataLabel entries for the tracer fluids.  This is not
    # necessary for Enzo, but yt (and other analysis codes) need it.
    for i in range(user_inputs['NumberOfTracerFluidFields']):
        print(i, i+user_inputs['NumberOfOriginalBaryonFields'])

        # we're assuming that the new fields are the last fields in that grid entry (which is currently true for tracer fields)
        newline = 'DataLabel[{:d}]             = '.format(i+user_inputs['NumberOfOriginalBaryonFields']) + 'TracerFluid' + '{:02d}'.format(i+1) + "\n"
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
        print("system did something fishy (parameter file moving), quitting.")
        sys.exit(1)
   
    return
    


################################################################################
def edit_hierarchy_file(user_inputs):
    '''

    edit_hierarchy_file

    This file edits the hierarchy file to add tracer fluid information.  As with the
    parameter file editing routine, it does a lot of error-checking, but that's probably
    a good thing here.

    What this 

    '''
    
    orig_hierarchy_file = user_inputs['dataset_directory'] + "/" + user_inputs['filename_stem'] + ".hierarchy"

    new_hierarchy_file = orig_hierarchy_file + ".new"

    backup_hierarchy_file = orig_hierarchy_file + ".orig"

    if(user_inputs['DEBUG_OUTPUTS']):
        print(orig_hierarchy_file)
        print(new_hierarchy_file)
        print(backup_hierarchy_file)


    ################################################################################

    # does original hierarchy file exist?  It should.
    if(os.path.exists(orig_hierarchy_file)==True):
        print("Original hierarchy file exists, continuing")

    # does new hierarchy file exist?  It should NOT at this point.
    if(os.path.exists(new_hierarchy_file)==True):
        print("New hierarchy file exists at this point and shouldn't. You may need to delete something. Exiting.")
        sys.exit(1)
    else:
        print("No new hierarchy file exists, continuing.")

    # Does backup hierarchy file exist?  It should NOT at this point.
    if(os.path.exists(backup_hierarchy_file)==True):
        print("Backup hierarchy file exists at this point and shouldn't. You may need to delete something. Exiting.")
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
        print("system did something fishy (hierarchy file), quitting.")
        sys.exit(1)

    # Now check to make sure the backup file really exists!
    if(os.path.exists(backup_hierarchy_file)==True):
        print("Backup hierarchy file exists at this point and SHOULD. YAY!")
    else:
        print("No new backup hierarchy file exists, and it should.  Something funny is happening.  Quitting.")
        sys.exit(1)

    # And now the magic happens.  We open the original hierarchy file and our new one,
    # read through the original file and modify the TracerFluid-related lines as we
    # encounter them.
    # 
    # The lines that we have to modify in each grid entry are the "NumberOfBaryonFields" lines,
    # which needs to be incremented by NumberOfTracerFluidFields, and then the FieldType line needs
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
            new_field_num = orig_field_num + user_inputs['NumberOfTracerFluidFields']
            split_line[2] = str(new_field_num)
            thisline = ' '.join(split_line) + "\n"

            if(user_inputs['DEBUG_OUTPUTS']):
                print("****new NumberOfBaryonFields line:", thisline)
            

        # Look for FieldType line.  Note that some lines have
        # length of zero, so check for that too.
        if len(split_line) > 0 and split_line[0] == 'FieldType':
            if(user_inputs['DEBUG_OUTPUTS']):
                print(split_line)

            for i in range(user_inputs['NumberOfTracerFluidFields']):
                split_line.append( str(106 + i)  )

            thisline = ' '.join(split_line) + "\n"

            if(user_inputs['DEBUG_OUTPUTS']):
                print("****new FieldType line:", thisline)

        # now we write either the original line or the modified line to the output file
        print(thisline, end = "", file=outputfile)


    print("\n", file=outputfile) # add a newline
        
    # close the output files.
    inputfile.close()
    outputfile.close()

    # TODO -- need to remove the starting hierarchy file (no extra extension) and move the .new one to the original name


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
        print("system did something fishy (hierarchy file moving), quitting.")
        sys.exit(1)
   

    
    ################################################################################

    return

