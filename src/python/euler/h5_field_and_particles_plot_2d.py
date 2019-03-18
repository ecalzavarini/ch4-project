# This script visualize the temperature field, the streamlines and the particles.

#!/usr/local/bin/python
import h5py
import numpy as np
import sys
import math
import matplotlib
import matplotlib.pyplot as plt


#f = h5py.File('../Data/RUN/field_2000.h5', 'r')
f = h5py.File('../Data/RUN_highRa1/field_5000.h5','r')
t=np.array(f['euler']['temperature'])
vx=np.array(f['euler']['velocity_x'])
vy=np.array(f['euler']['velocity_y'])
x=np.array(f['euler']['position_x'])
y=np.array(f['euler']['position_y'])

NZ = t.shape[0]
NY = t.shape[1]
NX = t.shape[2]
print("NZ= "+str(NZ))
print("NY= "+str(NY))
print("NX= "+str(NX))
tslice=t[0,:,:]
strivex=2
strivey=2
vxslice=vx[0,0:NY:strivey,0:NX:strivex]
vyslice=vy[0,0:NY:strivey,0:NX:strivex]
xslice=x[0,0:NY:strivey,0:NX:strivex]
yslice=y[0,0:NY:strivey,0:NX:strivex]
speed = np.sqrt(vxslice*vxslice + vyslice*vyslice)
lw = 1.5*speed / speed.max()


fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
ax.set_xlim(0, NX)
ax.set_ylim(0, NY)
ax.axis((0,NX,0,NY))
ax.imshow(tslice, cmap='seismic') #for color map https://matplotlib.org/tutorials/colors/colormaps.html
ax.streamplot(xslice,yslice,vxslice,vyslice, density=2,linewidth=lw,color="gray",arrowsize=0.5)


p = h5py.File('../Data/RUN/particle_2000.h5', 'r')
#p = h5py.File('../Data/RUN_highRa1/particle_5000.h5','r')

xpart=np.array(p['lagrange']['x'])
ypart=np.array(p['lagrange']['y'])
pxpart=np.array(p['lagrange']['px'])
pypart=np.array(p['lagrange']['py'])
ar=np.array(p['lagrange']['aspect_ratio'])

npart = ar.shape[0]
ar_types=np.asarray(list(set(ar)))
ar_types=np.sort(ar_types)
n_ar_types=len(ar_types)
print("ar_types = "+str(ar_types))
print("n_ar_types = "+str(n_ar_types))


# Now I add the ellipses
from matplotlib.patches import Ellipse

ells=[]
for i in range(npart):
    if ar[i]==ar_types[0]:
        ells.append(Ellipse(xy=(xpart[i], ypart[i]),
                            width=1.0, height=ar[i]/20,
                            angle=np.arcsin(pxpart[i])/np.pi*180., zorder=10))


for e in ells:
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_facecolor("black")

plt.savefig('visual.eps',dpi=256)    
plt.show()
