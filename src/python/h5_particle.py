#!/usr/local/bin/python
from __future__ import print_function
import h5py
import numpy as np
import sys
import os

fname = sys.argv[1]
#f = h5py.File('/Users/enrico/here/HIT/RUN/particle_5000.h5','r')
f = h5py.File( fname,'r') 

name=np.array( f['lagrange']['name'])
gyro=np.array(f['lagrange']['gyrotaxis_velocity'])
x=np.array( f['lagrange']['x'])
y=np.array( f['lagrange']['y'])
z=np.array( f['lagrange']['z'])

npart = gyro.shape[0]
which_part = float(sys.argv[2])


f = open('myfile','w')

for i in xrange(0,npart):
     if name[i]%13 == which_part :
          print(name[i],gyro[i],x[i],y[i],z[i],file=f)

f.close()


os.system("cat myfile") 





