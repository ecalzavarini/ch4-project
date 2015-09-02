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


# loop on all files
for file in filenames:
    if 'particle' in file:
        if 'h5' in file:
            fpath = dirname+"/"+file
            fin = h5py.File( fpath,'r')
            fpath = dirname+"/"+file[0:-3]+"_sort.h5"
            fout = h5py.File( fpath,'w')


# read the group
            group = fin['lagrange']
            group2 = fout.create_group('lagrange')

# sort name of the particles (and take the argument list)            
            idx = np.argsort( group['name'])

# we sort all the fields
            for key in group:
                buffer=np.array(group[key])
                buffer=buffer[idx]
                dset = group2.create_dataset(key, data=buffer)

            fout.close()
            fin.close()






