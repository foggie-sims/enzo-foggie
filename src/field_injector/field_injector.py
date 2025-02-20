'''
This is the driver routine to inject tracer fluid fields into an existing Enzo dataset.

WARNING: This program will modify a dataset in-place, so you should absolutely make a 
         backup copy of your original Enzo dataset before you start and then verify
         that the tracer fluids were correctly edited before you delete said backup copy.

All of the user inputs, including the routine that actually edits the .cpu files, can be
found in the file user_inputs.py.  You should not have to modify any other files.

Author:  Brian O'Shea (oshea@msu.edu), Feb. 2025
'''

from user_inputs import *
from mod_routines import *

if user_inputs['DEBUG_OUTPUTS']:
    print("*"*40)
    print("input user dictionary:\n")
    print(user_inputs)
    print("*"*40,"\n")

# creates the new parameter file (with tracer fluid contents in it)
edit_param_file(user_inputs)

# creates the new hierarchy file (with tracer fluid contents in it)
edit_hierarchy_file(user_inputs)

# creates the two new boundary conditions files (with tracer fluid contents in them)
edit_boundary_files(user_inputs)

# modifies the existing grid files to add tracer fluid fields.
modify_grid_files(user_inputs)
