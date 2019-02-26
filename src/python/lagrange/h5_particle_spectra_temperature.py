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

#spectra
def spt(series):
    n = len(series)
    data = np.asarray(series,dtype=np.float64)
    ft_data = np.array(n,dtype=np.complex128)
    spectrum = np.array(n,dtype=np.float64)
    ft_data = np.fft.rfft(data)
    spectrum  = (np.abs(ft_data))**2.
    return spectrum
                            

#autocorrelation function
def acf(series):
    n = len(series)
    data = np.asarray(series)
    mean = np.mean(data)
    c0 = np.sum((data - mean) ** 2) / float(n)

    def r(h):
        acf_lag = ((data[:n - h] - mean) * (data[h:] - mean)).sum() / float(n) / c0
        return round(acf_lag, 3)
    x = np.arange(n) # Avoiding lag 0 calculation
    acf_coeffs = map(r, x)
    return acf_coeffs

# create array
temperature = np.array([])


ntype = int(sys.argv[5])
totpart =  int(sys.argv[6])
print("totpart "+str(totpart))
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
                    temp=np.array( group['temperature'])


#            npart = name.shape[0]
#            print(ax[0])
                    which_part = int(sys.argv[4])+j           
#                    print(which_part)
                    temperature = np.append(temperature,temp[which_part])
                    fin.close()


                    
    print("len = ",len(temperature))
    tmp_tcorr = np.asarray(spt(temperature))
    length = len(temperature)
    print("len ",len(temperature))
    temperature = np.delete(temperature,range(len(temperature)))
    
    if j==0:
        tcorr = tmp_tcorr

    else:
        for k in xrange(0,len(tcorr)):
            tcorr[k] += tmp_tcorr[k]

    n=n+1.0
    print("n =",n)

    foutname = 'tspectra_'+sys.argv[4]+'.txt'
    fout = open(foutname,'w')
    for k in xrange(0,len(tcorr)):
        print(k/(10.*length),tcorr[k]/n,file=fout)
    fout.close()

