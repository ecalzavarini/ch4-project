#!/usr/local/bin/python
import h5py
import numpy as np
import sys , os
#import math
from math import *

# defining global types and fields
real_vector_dtype=np.dtype({ 'names':['x','y','z'], 'formats':[float,float,float] })
complex_vector_dtype=np.dtype({ 'names':['x','y','z'], 'formats':[complex,complex,complex] })



fname = sys.argv[1]
f = h5py.File( fname,'r') 
#f = h5py.File('../RUN/field_5000.h5', 'r')

vx=np.array(f['euler']['velocity_x'])
vy=np.array(f['euler']['velocity_y'])
vz=np.array(f['euler']['velocity_z'])

NZ = vx.shape[0]
NY = vx.shape[1]
NX = vx.shape[2]

ur = np.zeros((NX,NY,NZ),real_vector_dtype)
uc = np.zeros((NX,NY,NZ/2+1),complex_vector_dtype)

ur['x'] = vx
ur['y'] = vy
ur['z'] = vz 

###################################################
# Fourier  transform
def fftv(ur):
    uc = np.zeros((NX,NY,NZ/2+1),complex_vector_dtype) #this is an internal variable
    norm=1.#/ur['x'].size
    uc['x'] = np.fft.rfftn(ur['x'])*norm
    uc['y'] = np.fft.rfftn(ur['y'])*norm
    uc['z'] = np.fft.rfftn(ur['z'])*norm
    return uc;

###############################################
# Preliminary to spectra 

kx = np.fft.fftfreq(NX, 1./NX)
ky = np.fft.fftfreq(NY, 1./NY)
#kz = np.linspace(0, NZ/2, NZ/2+1, float)
kz = np.fft.rfftfreq(NZ, 1./NZ)


####################################################
# Compute spectra
def compute_spectra(uc):
    ene=0.0
    spect = np.zeros(NZ/2+1, float)
    for i in xrange(0,NX):
        for j in xrange(0,NY):
            for k in xrange(0,NZ/2+1):
                kabs=sqrt(kx[i]**2. + ky[j]**2. + kz[k]**2.)
                ik=int(ceil(kabs))
                if (k==0 or k==NZ/2):
                    fac=0.5
                else:
                    fac=1.0                    
                if ik < NZ/2+1:
                    ene = (np.abs(uc['x'][i,j,k]))**2. + (np.abs(uc['y'][i,j,k]))**2. + (np.abs(uc['y'][i,j,k]))**2.
                    spect[ik] += fac*ene;
    with open('spectra.dat', 'w') as f:
        for i in xrange(0,NZ/2+1):
            print >> f, i , spect[i]
    f.closed                
##################


uc = fftv(ur) 
compute_spectra(uc)
number=filter(str.isdigit, fname)
fname_new = 'spectra_'+number[:-1]+'.dat'
os.rename("spectra.dat", fname_new)
