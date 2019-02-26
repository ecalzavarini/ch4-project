#!/usr/local/bin/python
import h5py
import numpy as np
import sys
import math

increment_x = int(sys.argv[1])
increment_y = int(sys.argv[2])
increment_z = int(sys.argv[3])

f1 = h5py.File('pop.h5', 'r')
f2 = h5py.File('liquid_frac.h5', 'r')
t = np.array(f1['population']['temperature'])
v = np.array(f1['population']['velocity'])
l = np.array(f2['euler']['liquid_frac'])

NZ = t.shape[0]
NY = t.shape[1]
NX = t.shape[2]

NNZ = NZ + increment_z
NNY = NY + increment_y
NNX = NX + increment_x

new_t = np.zeros((NNZ,NNY,NNX),t.dtype)
new_v = np.zeros((NNZ,NNY,NNX),v.dtype)
new_l = np.zeros((NNZ,NNY,NNX),l.dtype)

for i in xrange(0,NNX):
    for j in xrange(0,NNY):
        for k in xrange(0,NNZ):
        	if i < NX and j < NY and k < NZ:
        		new_t[k,j,i] = t[k,j,i]
        		new_v[k,j,i] = v[k,j,i]
        		new_l[k,j,i] = l[k,j,i]

new_f = h5py.File('new_pop.h5', 'w')
grp = new_f.create_group('population')
dset_t = grp.create_dataset('temperature',(NNZ,NNY,NNX), new_t.dtype)
dset_v = grp.create_dataset('velocity',(NNZ,NNY,NNX), new_v.dtype)
dset_t[:] = new_t
dset_v[:] = new_v

new_f = h5py.File('new_liquid_frac.h5', 'w')
grp = new_f.create_group('euler')
dset_l = grp.create_dataset('liquid_frac',(NNZ,NNY,NNX), new_l.dtype)
dset_l[:] = new_l


