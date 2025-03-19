'''
Plot to make tracer fluid plots and print out some quantities (from a Sod shock tube with tracer fluids)

The tracer fluid fields should perfectly track the density.  The Sod shock tube setup makes each tracer
fluid field equal to the density divided by the number of the tracer fluid (so tracer fluid #1 is the
value of the density fluid, fluid #2 is 1/2 of the density value, fluid #3 is 1/3 of the density value, 
etc.).  This script plots the density and the *scaled* fractional value of the tracer fluid (so, it is
renormalized by the number of the tracer fluid so that the values of TracerFluid/Density should be 1 
everywhere).

'''

import yt
import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib.pyplot import cm

ds = yt.load("DD0001/data0001")

# Needed for 1D plot (will probably break if you use an AMR simulation instead of a single-grid simulation)
my_ortho_ray = ds.ortho_ray('x',[0.5,0.5])

# Error-check: are there actually tracer fluids?
if(ds.parameters['UseTracerFluid']==0):
    print("There are no tracer fluids in this simulation.  Exiting.")
    sys.exit(1)

# Set up the color map
color = cm.viridis(np.linspace(0,1,ds.parameters['NumberOfTracerFluidFields']))


###### FIGURE #1 -- DENSITY VALUES
# This is mostly just to give people something nice to look at.
plt.figure(clear=True)

# Start by plotting density
dens = plt.plot(my_ortho_ray['x'],my_ortho_ray['Density'],'k', linewidth=4, label='Density')

tf = []

# Loop over tracer fields and plot those.
for i in range(1,ds.parameters['NumberOfTracerFluidFields']+1):
    fieldname = 'TracerFluid' + '{:02d}'.format(i)
    tf.append(plt.plot(my_ortho_ray['x'],my_ortho_ray[fieldname],c=color[i-1],linestyle='--',label=fieldname))

plt.xlabel('x')
plt.ylabel('Values')
plt.title('Sod shock tube density and tracer fluid(s)')
plt.legend()
plt.savefig('TracerFluidDensities.png',bbox_inches='tight')

###### FIGURE #2 -- FRACTIONAL VALUES
# This should be an extremely boring plot - all of the lines should be
# right on top of each other.
plt.figure(clear=True)

tf_frac = []

# Small fraction to return minimum, maxmimum, and absolute value of the fractional difference 
# of the normalized tracer field (i.e., tracer fluid field / density, times the number of the field)
def minmaxfracdiff(field):
    return field.min(), field.max(), np.abs(field.max()-field.min())/(0.5*(field.min()+field.max()) )

# Loop over fields 
for i in range(1,ds.parameters['NumberOfTracerFluidFields']+1):
    fieldname = 'TracerFluid' + '{:02d}'.format(i)

    # normalized fractional value -- note that the float(i) at the end
    # does the normalization.
    fracval = my_ortho_ray[fieldname].d/my_ortho_ray['Density'].d*float(i)

    # get the min, max, abs. fractional diffs
    fracvals = minmaxfracdiff(fracval)

    # add fractional difference to tracer fluid name to make labels for plot
    mylabel = fieldname + ' ' + '({:.4e})'.format(fracvals[2])

    # make the plots
    tf_frac.append(plt.plot(my_ortho_ray['x'],fracval,c=color[i-1],linestyle='--',label=mylabel))

plt.legend()
plt.xlabel('x')
plt.ylabel('Fractional difference')
plt.title('Scaled fractional difference between density and tracer field\n(Value in parentheses should be near floating-point limit)')
plt.savefig('TracerFluidFractionalDifferences.png')

##### NOT A FIGURE, JUST PRINTING THE FRACTIONAL VALUES
for i in range(1,ds.parameters['NumberOfTracerFluidFields']+1):
    fieldname = 'TracerFluid' + '{:02d}'.format(i)
    fracval = my_ortho_ray[fieldname].d/my_ortho_ray['Density'].d*float(i)

    fracvals = minmaxfracdiff(fracval)

    print(fieldname,"min, max, fractional difference:", fracvals)

print("NOTE: Fractional differences should be a very small number, around 10^-16 or 10^-15 if Enzo was run in double precision.")
