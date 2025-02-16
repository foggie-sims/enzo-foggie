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
    "DEBUG_OUTPUTS": True  # True of False
}
