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
            x=np.array( group['x'])
            y=np.array( group['y'])
            z=np.array( group['z'])
            vx=np.array( group['vx'])
            vy=np.array( group['vy'])
            vz=np.array( group['vz'])
            ax=np.array( group['ax'])
            ay=np.array( group['ay'])
            az=np.array( group['az'])

            if  'px' in labels:
                px=np.array( group['px'])
                py=np.array( group['py'])
                pz=np.array( group['pz'])

            if 'temperature' in labels:
                tt=np.array( group['temperature'])
                
            npart = name.shape[0]
            which_part = int(sys.argv[2])

            foutname = 'trajectory_'+str(which_part)+'.dat'
            fout = open(foutname,'a')

            for i in xrange(0,npart):
               if name[i] == which_part :
                   if  'px' in labels:
                       if 'temperature' in labels:
                           print(name[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],ax[i],ay[i],az[i],px[i],py[i],pz[i],tt[i],file=fout)
                       else:                        
                           print(name[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],ax[i],ay[i],az[i],px[i],py[i],pz[i],file=fout)
                   else:
                       if 'temperature' in labels:
                           print(name[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],ax[i],ay[i],az[i],tt[i],file=fout)
                       else:
                           print(name[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],ax[i],ay[i],az[i],file=fout)                           
            fout.close()
            fin.close()

#os.system("cat myfile") 





