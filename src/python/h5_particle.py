#!/usr/local/bin/python                                                                                        
from __future__ import print_function
import h5py
import numpy as np
import sys
import os

import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


dirname = sys.argv[1]
filenames = os.listdir(dirname)
filenames = sorted(filenames, key=numericalSort)

for file in filenames:
    if 'particle' in file:
        if 'h5' in file:
            fpath = dirname+"/"+file
            #print(fpath)
            fin = h5py.File( fpath,'r')

            name=np.array( fin['lagrange']['name'])
            x=np.array( fin['lagrange']['x'])
            y=np.array( fin['lagrange']['y'])
            z=np.array( fin['lagrange']['z'])
            vx=np.array( fin['lagrange']['vx'])
            vy=np.array( fin['lagrange']['vy'])
            vz=np.array( fin['lagrange']['vz'])
            ax=np.array( fin['lagrange']['ax'])
            ay=np.array( fin['lagrange']['ay'])
            az=np.array( fin['lagrange']['az'])

            npart = name.shape[0]
            which_part = int(sys.argv[2])

            foutname = 'trajectory_'+str(which_part)+'.dat'
            fout = open(foutname,'a')

            for i in xrange(0,npart):
               if name[i] == which_part :
                    print(name[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],ax[i],ay[i],az[i],file=fout)

            fout.close()


#os.system("cat myfile") 





