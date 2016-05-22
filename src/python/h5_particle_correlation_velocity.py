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

#autocorrelation function
def acf(series):
    n = len(series)
    data = np.asarray(series)
    mean = 0.0 #np.mean(data)
    c0 = np.sum((data - mean) ** 2) / float(n)

    def r(h):
        acf_lag = ((data[:n - h] - mean) * (data[h:] - mean)).sum() / float(n-h)
        return acf_lag
    x = np.arange(n) # Avoiding lag 0 calculation
    acf_coeffs = map(r, x)
    return acf_coeffs

#create arrays
acceleration_x = np.array([])
acceleration_y = np.array([])
acceleration_z = np.array([])

#number of particle types 
ntype = int(sys.argv[5])
# total number of particles
totpart = int(sys.argv[6])
n=0.0
for j in xrange(0,totpart,ntype):
    # loop on all files
    for i in xrange(int(dirnumstart),int(dirnumend)+1):
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

#append the data of a given particle
                    which_part = int(sys.argv[4])+j           
#                    print(which_part)
                    acceleration_x = np.append(acceleration_x,ax[which_part])
                    acceleration_y = np.append(acceleration_y,ay[which_part])
                    acceleration_z = np.append(acceleration_z,az[which_part])
                    fin.close()

#compute spectra 
    tmp_acorr_x = np.asarray(acf(acceleration_x))
    tmp_acorr_y = np.asarray(acf(acceleration_y))
    tmp_acorr_z = np.asarray(acf(acceleration_z))
    print("length ",len(acceleration_x))
    print("length ",len(acceleration_y))
    print("length ",len(acceleration_z))

#compute mean
    tmp_amean_x = np.mean(acceleration_x)
    tmp_amean_y = np.mean(acceleration_y)
    tmp_amean_z = np.mean(acceleration_z)

#free the array
    acceleration_x = np.delete(acceleration_x,range(len(acceleration_x)))
    acceleration_y = np.delete(acceleration_y,range(len(acceleration_y)))
    acceleration_z = np.delete(acceleration_z,range(len(acceleration_z)))

    if j==0:
        acorr_x = tmp_acorr_x
        acorr_y = tmp_acorr_y
        acorr_z = tmp_acorr_z

        amean_x = tmp_amean_x
        amean_y = tmp_amean_y
        amean_z = tmp_amean_z
    else:
        for k in xrange(0,len(acorr_x)):
            acorr_x[k] += tmp_acorr_x[k]
            acorr_y[k] += tmp_acorr_y[k]
            acorr_z[k] += tmp_acorr_z[k]

        amean_x += tmp_amean_x
        amean_y += tmp_amean_y
        amean_z += tmp_amean_z

    n=n+1.0
    print("analyzed trajectories -> n = ",n)

    foutname = 'vcorr_'+sys.argv[4]+'.txt'
    fout = open(foutname,'w')
    for k in xrange(0,len(acorr_x)):
#        print(k*10,acorr_x[k]/acorr_x[0],acorr_y[k]/acorr_y[0],acorr_z[k]/acorr_z[0],amean_x/n,amean_y/n,amean_z/n,file=fout)
        print(k*10,(acorr_x[k]-amean_x**2.)/(acorr_x[0]-amean_x**2.),(acorr_y[k]-amean_y**2.)/(acorr_y[0]-amean_y**2.),(acorr_z[k]-amean_z**2.)/(acorr_z[0]-amean_z**2.),amean_x/n,amean_y/n,amean_z/n,file=fout)
    fout.close()

