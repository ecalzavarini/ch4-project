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
#        return round(acf_lag, 3)
        return acf_lag 
    x = np.arange(n) # Avoiding lag 0 calculation
    acf_coeffs = map(r, x)
    return acf_coeffs

# create array
temperature = np.array([])


ntype = int(sys.argv[5])
n=0.0
for j in xrange(0,22000,ntype):
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
                    temp=np.array( group['temperature'])


#            npart = name.shape[0]
#            print(ax[0])
                    which_part = int(sys.argv[4])+j           
#                    print(which_part)
                    temperature = np.append(temperature,temp[which_part])
                    fin.close()


    print("len = ",len(temperature))
    tmp_tcorr = np.asarray(acf(temperature))
    temperature = np.delete(temperature,range(len(temperature)))

    if j==0:
        tcorr = tmp_tcorr

    else:
        for k in xrange(0,len(tcorr)):
            tcorr[k] += tmp_tcorr[k]

    n=n+1.0
    print("n =",n)

    foutname = 'tcorr_'+sys.argv[4]+'.txt'
    fout = open(foutname,'w')
    for k in xrange(0,len(tcorr)):
        print(k*10,tcorr[k]/tcorr[0],file=fout)
    fout.close()

