#import yt
#import numpy as np
#import h5py
#import os

from user_inputs import *
from mod_routines import *

if user_inputs['DEBUG_OUTPUTS']:
    print("*"*40)
    print("input user dictionary:\n")
    print(user_inputs)
    print("*"*40,"\n")

#edit_param_file(user_inputs)

edit_hierarchy_file(user_inputs)

