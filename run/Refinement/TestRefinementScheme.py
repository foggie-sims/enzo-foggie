import os
import yt
import numpy as np

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

sim_dir = os.path.basename(_dir_name)

def test_multi_refinement_regions():
    ds = yt.load(os.path.join(_dir_name, 'DD0008/DD0008'))
    data  = ds.all_data()
    current_time = ds.current_time

    # Read in track box coordinates for this snapshot
    t,x1,y1,z1,x2,y2,z2,llim = np.loadtxt('RefinementTestTrackFile.txt',skiprows=2,usecols=(1,2,3,4,5,6,7,8),unpack=True)

    # Find expected positions of MRR 1 & MRR 2 at this time
    tfrac = (current_time-t[0])/(t[1]-t[0])
    MRR1_corners = np.array([x1[0]+(x1[1]-x1[0])*tfrac,y1[0]+(y1[1]-y1[0])*tfrac,z1[0]+(z1[1]-z1[0])*tfrac,\
                                  x2[0]+(x2[1]-x2[0])*tfrac,y2[0]+(y2[1]-y2[0])*tfrac,z2[0]+(z2[1]-z2[0])*tfrac])
    MRR2_corners = np.array([x1[2]+(x1[3]-x1[2])*tfrac,y1[2]+(y1[3]-y1[2])*tfrac,z1[2]+(z1[3]-z1[2])*tfrac,\
                                  x2[2]+(x2[3]-x2[2])*tfrac,y2[2]+(y2[3]-y2[2])*tfrac,z2[2]+(z2[3]-z2[2])*tfrac])
    
    # Create boxes that correspond to where our MRRs should be and check that all cells within them meet the minimum
    # refinement criteria specified in the track file
    MRR1 = ds.box(MRR1_corners[:3],MRR1_corners[3:])
    MRR2 = ds.box(MRR2_corners[:3],MRR2_corners[3:])
    cellvolumes = np.sort(np.unique(data['index','cell_volume']))
    maxlevel = ds.parameters['MaximumRefinementLevel']
    assert (np.all(MRR1['index','cell_volume']<=cellvolumes[int(maxlevel-llim[0])]))
    assert (np.all(MRR2['index','cell_volume']<=cellvolumes[int(maxlevel-llim[2])]))

    # We are ONLY refining by MultiRefineRegion, so check that refinement isn't happening outside of these regions
    FullRegion = ds.box([0,0,0],[1,1,1])
    NoMRR = FullRegion-MRR1-MRR2
    NoMRR_refined = NoMRR['index','cell_volume']<cellvolumes[-1]
    assert(len(NoMRR['index','cell_volume'][NoMRR_refined])/len(NoMRR['index','cell_volume'])<0.95)
    
