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
filenames= os.listdir(dirname+str(dirnumstart))
filenames = sorted(filenames, key=numericalSort)

# Compute structure function of order 2 for a given variable
def structure_function_order2(series):
    n = len(series)
    data = np.asarray(series)
    def r(h):
        acf_lag = ((data[:n - h] - data[h:])**2.).sum() / float(n-h)
        return acf_lag
    x = np.arange(n) 
    acf_coeffs = map(r, x)
    return acf_coeffs

# Create array
position_x = np.array([])
position_y = np.array([])
position_z = np.array([])
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
                    posx=np.array( group['x'])
                    posy=np.array( group['y'])
                    posz=np.array( group['z'])
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

                    
#            npart = name.shape[0]
#            print(ax[0])
                    which_part = int(sys.argv[4])+j           
#                    print(which_part)
                    position_x = np.append(position_x,posx[which_part])
                    position_y = np.append(position_y,posy[which_part])
                    position_z = np.append(position_z,posz[which_part])
                    acceleration_x = np.append(acceleration_x,ax[which_part])
                    acceleration_y = np.append(acceleration_y,ay[which_part])
                    acceleration_z = np.append(acceleration_z,az[which_part])
                    acceleration_m = np.append(acceleration_m,am[which_part])
                    temperature = np.append(temperature,tt[which_part])
                    fin.close()

 
#print(acceleration)
# first we pad with 0 , because our signal is not periodic                    
    position_x = np.append(position_x,np.zeros(len(position_x)))
    position_y = np.append(position_y,np.zeros(len(position_y)))
    position_z = np.append(position_z,np.zeros(len(position_z)))
    acceleration_x = np.append(acceleration_x,np.zeros(len(acceleration_x)))
    acceleration_y = np.append(acceleration_y,np.zeros(len(acceleration_y)))
    acceleration_z = np.append(acceleration_z,np.zeros(len(acceleration_z)))
    acceleration_m = np.append(acceleration_m,np.zeros(len(acceleration_m)))
    temperature = np.append(temperature,np.zeros(len(temperature)))

    tmp_acorr_posx = np.asarray(structure_function_order2(position_x))
    tmp_acorr_posy = np.asarray(structure_function_order2(position_y))
    tmp_acorr_posz = np.asarray(structure_function_order2(position_z))
    tmp_acorr_x = np.asarray(structure_function_order2(acceleration_x))
    tmp_acorr_y = np.asarray(structure_function_order2(acceleration_y))
    tmp_acorr_z = np.asarray(structure_function_order2(acceleration_z))
    tmp_acorr_m = np.asarray(structure_function_order2(acceleration_m))
    tmp_acorr_t = np.asarray(structure_function_order2(temperature))
    length = len(acceleration_x)
    print("len ",len(acceleration_x))
    print("len ",len(acceleration_y))
    print("len ",len(acceleration_z))            
    position_x = np.delete(position_x,range(len(position_x)))
    position_y = np.delete(position_y,range(len(position_y)))
    position_z = np.delete(position_z,range(len(position_z)))
    acceleration_x = np.delete(acceleration_x,range(len(acceleration_x)))
    acceleration_y = np.delete(acceleration_y,range(len(acceleration_y)))
    acceleration_z = np.delete(acceleration_z,range(len(acceleration_z)))
    acceleration_m = np.delete(acceleration_m,range(len(acceleration_m)))
    temperature = np.delete(temperature,range(len(temperature)))
    if j==0:
        acorr_posx = tmp_acorr_posx
        acorr_posy = tmp_acorr_posy
        acorr_posz = tmp_acorr_posz        
        acorr_x = tmp_acorr_x
        acorr_y = tmp_acorr_y
        acorr_z = tmp_acorr_z
        acorr_m = tmp_acorr_m
        acorr_t = tmp_acorr_t
    else:
        for k in range(0,len(acorr_x)):
            acorr_posx[k] += tmp_acorr_posx[k]
            acorr_posy[k] += tmp_acorr_posy[k]
            acorr_posz[k] += tmp_acorr_posz[k]
            acorr_x[k] += tmp_acorr_x[k]
            acorr_y[k] += tmp_acorr_y[k]
            acorr_z[k] += tmp_acorr_z[k]
            acorr_m[k] += tmp_acorr_m[k]
            acorr_t[k] += tmp_acorr_t[k]

    n=n+1.0
    print("n =",n)

    foutname = 'displacement_particle_'+sys.argv[4]+'.txt'
    fout = open(foutname,'w')
    for k in range(0,len(acorr_x)):
        print(k*10.,acorr_posx[k]/n,acorr_posy[k]/n,acorr_posz[k]/n, acorr_x[k]/n,acorr_y[k]/n,acorr_z[k]/n,acorr_t[k]/n,acorr_m[k]/n,file=fout)
    fout.close()



