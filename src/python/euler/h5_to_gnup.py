#!/usr/local/bin/python
import h5py
import numpy as np
import sys

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

fname = sys.argv[1]
#f = h5py.File('RUN/field_2000.h5', 'r')
f = h5py.File(fname, 'r')
#f['euler']['position_x']
#f['euler']['position_y']
#f['euler']['position_z']
#field_dtype = np.dtype({ 'names':['x','y','z','field'], 'formats':[float,float,float,float] })
x = np.array( f['euler']['position_x'])
y = np.array( f['euler']['position_y'])
z = np.array( f['euler']['position_z'])
t = np.array( f['euler']['temperature'])
vx = np.array( f['euler']['velocity_x'])
vy = np.array( f['euler']['velocity_y'])
vz = np.array( f['euler']['velocity_z'])

NZ = t.shape[0]
NY = t.shape[1]
NX = t.shape[2]

for i in xrange(0,NX):
    for j in xrange(0,NY):
        for k in xrange(0,NZ):
            print x[k,j,i],y[k,j,i],z[k,j,i],t[k,j,i], vx[k,j,i], vy[k,j,i], vz[k,j,i]
    print 

