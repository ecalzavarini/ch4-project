#!/usr/local/bin/python
import h5py
import numpy as np
import sys
#import math
from math import *

# defining global types and fields
#real_vector_dtype=np.dtype({ 'names':['x','y','z'], 'formats':[float,float,float] })
#complex_vector_dtype=np.dtype({ 'names':['x','y','z'], 'formats':[complex,complex,complex] })

fname = sys.argv[1]
f = h5py.File( fname,'r')
#f = h5py.File('../RUN/field_5000.h5', 'r')    
#f = h5py.File('/home/enrico/HIT/RUN/field_30000.h5', 'r')

t=np.array(f['euler']['temperature'])

NZ = t.shape[0]
NY = t.shape[1]
NX = t.shape[2]

tr = np.zeros((NX,NY,NZ),float)
tc = np.zeros((NX,NY,NZ/2+1),complex)

tr = t

###################################################
# Fourier  transform
def ffts(tr):
    tc = np.zeros((NX,NY,NZ/2+1),complex) #this is internal variable
    norm=1.#/ur['x'].size
    tc = np.fft.rfftn(tr)*norm
    return tc;

##########
# Preliminary to spectra 

kx = np.fft.fftfreq(NX, 1./NX)
ky = np.fft.fftfreq(NY, 1./NY)
kz = np.linspace(0, NZ/2, NZ/2+1, float)

####################################################
# Compute spectra
def compute_spectra_s(tc):
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
                ene = abs(tc[i,j,k])**2. 
                if ik < NZ/2+1:
                    spect[ik] += fac*ene;
    with open('spectra_temperature.dat', 'w') as f:
        for i in xrange(0,NZ/2+1):
            print >> f, i , spect[i]
    f.closed                
##################


tc = ffts(tr) 
compute_spectra_s(tc)

