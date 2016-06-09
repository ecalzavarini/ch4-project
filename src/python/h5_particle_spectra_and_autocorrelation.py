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


# create array
acceleration_x = np.array([])
acceleration_y = np.array([])
acceleration_z = np.array([])
acceleration_m = np.array([])
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
                    if 'swim_velocity':
                        ax=np.array( group['vx'])-np.array( group['px'])*np.array( group['swim_velocity'])
                        ay=np.array( group['vy'])-np.array( group['py'])*np.array( group['swim_velocity'])
                        az=np.array( group['vz'])-np.array( group['pz'])*np.array( group['swim_velocity'])
                    am=np.sqrt(np.power(ax,2.) + np.power(ay,2.) + np.power(az,2.))
                    if 'temperature' in labels: 
                        tt=np.array( group['temperature'])

# remove man value
                    ax=ax-np.mean(ax)
                    ay=ay-np.mean(ay)
                    az=az-np.mean(az)
                    am=am-np.mean(am)
                    tt=tt-np.mean(tt)
                    
#            npart = name.shape[0]
#            print(ax[0])
                    which_part = int(sys.argv[4])+j           
#                    print(which_part)
                    acceleration_x = np.append(acceleration_x,ax[which_part])
                    acceleration_y = np.append(acceleration_y,ay[which_part])
                    acceleration_z = np.append(acceleration_z,az[which_part])
                    acceleration_m = np.append(acceleration_m,am[which_part])
                    temperature = np.append(temperature,tt[which_part])
                    fin.close()

 
#print(acceleration)
# first we pad with 0 , because our signal is not periodic                    
    acceleration_x = np.append(acceleration_x,np.zeros(len(acceleration_x)))
    acceleration_y = np.append(acceleration_y,np.zeros(len(acceleration_y)))
    acceleration_z = np.append(acceleration_z,np.zeros(len(acceleration_z)))
    acceleration_m = np.append(acceleration_m,np.zeros(len(acceleration_m)))
    temperature = np.append(temperature,np.zeros(len(temperature)))

    tmp_acorr_x = np.asarray(spt(acceleration_x))
    tmp_acorr_y = np.asarray(spt(acceleration_y))
    tmp_acorr_z = np.asarray(spt(acceleration_z))
    tmp_acorr_m = np.asarray(spt(acceleration_m))
    tmp_acorr_t = np.asarray(spt(temperature))
    length = len(acceleration_x)
    print("len ",len(acceleration_x))
    print("len ",len(acceleration_y))
    print("len ",len(acceleration_z))            
    acceleration_x = np.delete(acceleration_x,range(len(acceleration_x)))
    acceleration_y = np.delete(acceleration_y,range(len(acceleration_y)))
    acceleration_z = np.delete(acceleration_z,range(len(acceleration_z)))
    acceleration_m = np.delete(acceleration_m,range(len(acceleration_m)))
    temperature = np.delete(temperature,range(len(temperature)))
    if j==0:
        acorr_x = tmp_acorr_x
        acorr_y = tmp_acorr_y
        acorr_z = tmp_acorr_z
        acorr_m = tmp_acorr_m
        acorr_t = tmp_acorr_t
    else:
        for k in range(0,len(acorr_x)):
            acorr_x[k] += tmp_acorr_x[k]
            acorr_y[k] += tmp_acorr_y[k]
            acorr_z[k] += tmp_acorr_z[k]
            acorr_m[k] += tmp_acorr_m[k]
            acorr_t[k] += tmp_acorr_t[k]

    n=n+1.0
    print("n =",n)

    foutname = 'spectra_particle_'+sys.argv[4]+'.txt'
    fout = open(foutname,'w')
    for k in range(0,len(acorr_x)):
        print(k/(10.*length),acorr_x[k]/n,acorr_y[k]/n,acorr_z[k]/n,acorr_t[k]/n,acorr_m[k]/n,file=fout)
    fout.close()

    foutname = 'acorr_particle_'+sys.argv[4]+'.txt'
    fout = open(foutname,'w')
    acorr_x_nozero = acorr_x
    acorr_y_nozero = acorr_y
    acorr_z_nozero = acorr_z
    acorr_m_nozero = acorr_m
    acorr_t_nozero = acorr_t

#    acorr_x_nozero[0] = 0.0;
#    acorr_y_nozero[0] = 0.0;
#    acorr_z_nozero[0] = 0.0;
#    acorr_m_nozero[0] = 0.0;
#    acorr_t_nozero[0] = 0.0;

    real_acorr_x = np.fft.irfft(acorr_x_nozero)
    real_acorr_y = np.fft.irfft(acorr_y_nozero)
    real_acorr_z = np.fft.irfft(acorr_z_nozero)
    real_acorr_m = np.fft.irfft(acorr_m_nozero)
    real_acorr_t = np.fft.irfft(acorr_t_nozero)

    for k in range(0,int(len(real_acorr_x)/2)):
        print(k*10., real_acorr_x[k]/real_acorr_x[0], real_acorr_y[k]/real_acorr_y[0], real_acorr_z[k]/real_acorr_z[0],real_acorr_t[k]/real_acorr_t[0],real_acorr_m[k]/real_acorr_m[0],file=fout)
    fout.close()


