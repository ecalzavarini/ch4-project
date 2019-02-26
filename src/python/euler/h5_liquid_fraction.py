#!/usr/local/bin/python
import h5py
import numpy as np
import sys
import math
import os

Dir = 'Data'
if not os.path.exists(Dir):
    os.makedirs(Dir)

for file in os.listdir('RUN'):
    if file.endswith('.h5'):
        print(file)
        try:
            f = h5py.File('RUN/' + file, 'r')
            lf = np.array(f['euler']['liquid_fraction'])
            NZ = lf.shape[0]
            NY = lf.shape[1]
            NX = lf.shape[2]
            new_f = h5py.File(Dir + '/' + file, 'w')
            grp = new_f.create_group('euler')
            dset_v = grp.create_dataset('liquid_fraction',(NZ,NY,NX), lf.dtype)
            dset_v[:] = lf
        except:
            print 'Unexpected error'