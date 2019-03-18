#!/usr/local/bin/python
import h5py
import numpy as np
#import sys
#import math
#import matplotlib
import matplotlib.pyplot as plt


#f = h5py.File('../Data/RUN/field_500.h5', 'r')
f = h5py.File('field_5000.h5','r')
t=np.array(f['euler']['temperature'])
vx=np.array(f['euler']['velocity_x'])
vy=np.array(f['euler']['velocity_y'])

NZ = t.shape[0]
NY = t.shape[1]
NX = t.shape[2]
print("NZ= "+str(NZ))
print("NY= "+str(NY))
print("NX= "+str(NX))
tslice=t[0,:,:]
vxslice=vx[0,:,:]
vyslice=vy[0,:,:]

#for ii in range(0,NZ):
#    for jj in range(0,NY):
#        for kk in range(0,NZ):
#           t[k,j,i]

#print(tslice)
#plt.imshow(tslice, cmap='hot', interpolation='nearest')
plt.imshow(tslice)
plt.show()



