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
dirnumstart = sys.argv[2]
dirnumend = sys.argv[3]
#dirnames = []
#filenames = []
filenames= os.listdir(dirname+str(dirnumstart))
filenames = sorted(filenames, key=numericalSort)
#print(filenames)

nbins=200

# create array
acceleration_x = np.array([])
acceleration_y = np.array([])
acceleration_z = np.array([])
temperature = np.array([])

ntype = int(sys.argv[5])
totpart =  int(sys.argv[6])
print("totpart"+str(totpart))
n=0.0
for j in range(0,totpart,ntype):
    # loop on all files
    for i in range(int(dirnumstart),int(dirnumend)+1):
        for file in filenames:
            if 'particle' in file:
                if 'h5' and 'sort'  in file:            
                    fpath = dirname+str(i)+"/"+file
                    #print(fpath)
                    fin = h5py.File( fpath,'r')

# read the group
                    group = fin['lagrange']
                    labels = group.keys()
#read the particle data
                    name=np.array( group['name'])
                    ax=np.array( group['vx'])
                    ay=np.array( group['vy'])
                    az=np.array( group['vz'])
                    if 'temperature' in labels: 
                        tt=np.array( group['temperature'])

#            npart = name.shape[0]
#            print(ax[0])
                    which_part = int(sys.argv[4])+j           
#                    print(which_part)
                    acceleration_x = np.append(acceleration_x,ax[which_part])
                    acceleration_y = np.append(acceleration_y,ay[which_part])
                    acceleration_z = np.append(acceleration_z,az[which_part])
                    temperature = np.append(temperature,tt[which_part])
                    fin.close()

 
#print(acceleration)
    tmp_acorr_x,bin_edges = np.histogram(acceleration_x,bins=nbins,range=(-.2,.2),normed=True)
    tmp_acorr_y,bin_edges = np.histogram(acceleration_y,bins=nbins,range=(-.2,.2),normed=True)
    tmp_acorr_z,bin_edges = np.histogram(acceleration_z,bins=nbins,range=(-.2,.2),normed=True)
    tmp_acorr_t,bin_edges = np.histogram(temperature,bins=nbins,range=(-.2,.2),normed=True)
    length = len(acceleration_x)
    print("len ",len(acceleration_x))
    print("len ",len(acceleration_y))
    print("len ",len(acceleration_z))            
    acceleration_x = np.delete(acceleration_x,range(len(acceleration_x)))
    acceleration_y = np.delete(acceleration_y,range(len(acceleration_y)))
    acceleration_z = np.delete(acceleration_z,range(len(acceleration_z)))
    temperature = np.delete(temperature,range(len(temperature)))
    if j==0:
        acorr_x = tmp_acorr_x
        acorr_y = tmp_acorr_y
        acorr_z = tmp_acorr_z
        acorr_t = tmp_acorr_t
    else:
        for k in range(0,len(acorr_x)):
            acorr_x[k] += tmp_acorr_x[k]
            acorr_y[k] += tmp_acorr_y[k]
            acorr_z[k] += tmp_acorr_z[k]
            acorr_t[k] += tmp_acorr_t[k]

    n=n+1.0
    print("n =",n)

    foutname = 'histogram_particle_'+sys.argv[4]+'.txt'
    fout = open(foutname,'w')
    for k in range(0,len(acorr_x)):
        print((bin_edges[k]+bin_edges[k+1])/2.,acorr_x[k]/n,acorr_y[k]/n,acorr_z[k]/n,acorr_t[k]/n,file=fout)
    fout.close()



