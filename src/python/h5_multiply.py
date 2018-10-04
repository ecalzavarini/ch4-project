#!/usr/local/bin/python
import h5py
import numpy as np
import sys
import math

def wrap(x,s):
    return int( (x) - math.floor((x) / s) * s )

scale_x = float(sys.argv[1])
scale_y = float(sys.argv[2])
scale_z = float(sys.argv[3])

f = h5py.File('pop.h5', 'r')
#t=np.array(f['population']['temperature'])
v=np.array(f['population']['velocity'])

NZ = v.shape[0]
NY = v.shape[1]
NX = v.shape[2]

NNZ = int(scale_z*NZ)
NNY = int(scale_y*NY)
NNX = int(scale_x*NX)

#new_t = np.zeros((NNZ,NNY,NNX),t.dtype)
new_v = np.zeros((NNZ,NNY,NNX),v.dtype)

for ii in xrange(0,NNX):
    i=wrap(ii,NX)
    for jj in xrange(0,NNY):
        j=wrap(jj,NY)
        for kk in xrange(0,NNZ):
            k=wrap(kk,NZ)
#           new_t[kk,jj,ii] = t[k,j,i]
            new_v[kk,jj,ii] = v[k,j,i]

new_f = h5py.File('new_pop.h5', 'w')
grp = new_f.create_group('population')
#dset_t = grp.create_dataset('temperature',(NNZ,NNY,NNX), new_t.dtype)
dset_v = grp.create_dataset('velocity',(NNZ,NNY,NNX), new_v.dtype)
#dset_t[:] = new_t
dset_v[:] = new_v

