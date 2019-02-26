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
            group = fin['lagrange']
            labels = group.keys()
            #nl = len(labels)
            
            name=np.array( group['name'])
            ax=np.array( group['ax'])
            ay=np.array( group['ay'])
            az=np.array( group['az'])
            dx_ux=np.array( group['dx_ux'])
            dx_uy=np.array( group['dx_uy'])
            dx_uz=np.array( group['dx_uz'])
            dy_ux=np.array( group['dy_ux'])
            dy_uy=np.array( group['dy_uy'])
            dy_uz=np.array( group['dy_uz'])
            dz_ux=np.array( group['dz_ux'])
            dz_uy=np.array( group['dz_uy'])
            dz_uz=np.array( group['dz_uz'])

            npart = name.shape[0]
            which_part = int(sys.argv[2])
            ntypes = int(sys.argv[3])

            foutname = 'gradients_'+str(which_part)+'.dat'
            fout = open(foutname,'a')

            for i in xrange(0,npart):
               if name[i]%ntypes == which_part :
                    print(name[i],ax[i],ay[i],az[i],dx_ux[i],dy_uy[i],dz_uz[i],dy_ux[i],dz_ux[i],dx_uy[i],dz_uy[i],dx_uz[i],dy_uz[i],file=fout)
            fout.close()
            fin.close()

#os.system("cat myfile") 





