#!/usr/local/bin/python                                                                                        
from __future__ import print_function
import h5py
import numpy as np
import sys
import os
import glob
import re

from numpy import linalg as LA

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = list(map(int, parts[1::2]))
    return parts

order = int(sys.argv[5])
dirname = sys.argv[1]
dirnumstart = sys.argv[2]
dirnumend = sys.argv[3]
filenames = glob.glob(dirname+"["+str(dirnumstart)+"-"+str(dirnumend)+"]/particle*sort*h5")
#filenames= os.listdir(dirname+str(dirnumstart))
filenames = sorted(filenames, key=numericalSort)
howmanyfiles = len(filenames)

for file in filenames:
    print(file)

print(howmanyfiles)
    
ntype = int(sys.argv[4])
#totpart =  int(sys.argv[5])
print("ntype = "+str(ntype))
#print("totpart = "+str(totpart))


moment_p_times_gradt = np.zeros(ntype,dtype=np.float64)
moment_gradt = np.zeros(ntype,dtype=np.float64)
moment_p = np.zeros(ntype,dtype=np.float64)
moment_dt_p = np.zeros(ntype,dtype=np.float64)
moment_p_times_gradt = np.zeros(ntype,dtype=np.float64)
moment_tt = np.zeros(ntype,dtype=np.float64)
moment_t = np.zeros(ntype,dtype=np.float64)
moment_dt_t = np.zeros(ntype,dtype=np.float64)
counter = np.zeros(ntype,dtype=np.float64) 
C = np.zeros(ntype,dtype=np.float64)
moment_ax = np.zeros(ntype,dtype=np.float64)
moment_ay = np.zeros(ntype,dtype=np.float64)
moment_az = np.zeros(ntype,dtype=np.float64)
moment_eps = np.zeros(ntype,dtype=np.float64)
moment_enst = np.zeros(ntype,dtype=np.float64)
moment_u2 = np.zeros(ntype,dtype=np.float64)
moment_hyperb = np.zeros(ntype,dtype=np.float64)
moment_hyperb_plus = np.zeros(ntype,dtype=np.float64)
moment_hyperb_minus = np.zeros(ntype,dtype=np.float64)
moment_p_times_vorticity = np.zeros(ntype,dtype=np.float64)
Cvorticity = np.zeros(ntype,dtype=np.float64)
moment_p_times_e1 = np.zeros(ntype,dtype=np.float64)
moment_p_times_e2 = np.zeros(ntype,dtype=np.float64)
moment_p_times_e3 = np.zeros(ntype,dtype=np.float64)
Cpstrain1 = np.zeros(ntype,dtype=np.float64)
Cpstrain2 = np.zeros(ntype,dtype=np.float64)
Cpstrain3 = np.zeros(ntype,dtype=np.float64)
Cvorticitystrain1 = np.zeros(ntype,dtype=np.float64)
Cvorticitystrain2 = np.zeros(ntype,dtype=np.float64)
Cvorticitystrain3 = np.zeros(ntype,dtype=np.float64)
Cgradtstrain1 = np.zeros(ntype,dtype=np.float64)
Cgradtstrain2 = np.zeros(ntype,dtype=np.float64)
Cgradtstrain3 = np.zeros(ntype,dtype=np.float64)
#ntype = int(sys.argv[5])
#totpart =  int(sys.argv[6])
#print("totpart"+str(totpart))


# loop on all files across different directories
for i in range(0,howmanyfiles-1):

    file = filenames[i]
    fpath = file
    fin = h5py.File( fpath,'r')
# Read the group
    group = fin['lagrange']
    labels = group.keys()
# Read the particle data
    name=np.array( group['name'])
    ax=np.array( group['ax'])
    ay=np.array( group['ay'])
    az=np.array( group['az'])
    dx_ux=np.array( group['dx_ux'])
    dx_uy=np.array( group['dx_uy'])
    dx_uz=np.array( group['dx_uz'])
    dy_ux=np.array( group['dy_ux'])
    dy_uy=np.array( group['dy_uy'])
    dy_uz=np.array( group['dy_uz'])
    dz_ux=np.array( group['dz_ux'])
    dz_uy=np.array( group['dz_uy'])
    dz_uz=np.array( group['dz_uz'])                
    if 'swim_velocity' in labels:
        ux=np.array( group['vx'])-np.array( group['px'])*np.array( group['swim_velocity'])
        uy=np.array( group['vy'])-np.array( group['py'])*np.array( group['swim_velocity'])
        uz=np.array( group['vz'])-np.array( group['pz'])*np.array( group['swim_velocity'])                           
        px=np.array( group['px'])
        py=np.array( group['py'])
        pz=np.array( group['pz'])   
        dt_px=np.array( group['dt_px'])
        dt_py=np.array( group['dt_py'])
        dt_pz=np.array( group['dt_pz'])                        
    if 'temperature' in labels: 
        tt=np.array( group['temperature'])
        dx_t=np.array( group['dx_t'])
        dy_t=np.array( group['dy_t'])
        dz_t=np.array( group['dz_t'])                        
    fin.close()

    file = filenames[i+1]
    fpath = file
    fin = h5py.File( fpath,'r')
# Read the group
    group = fin['lagrange']
    labels = group.keys()
# Read the particle data
    name_later=np.array( group['name'])  
    if 'temperature' in labels: 
        tt_later=np.array( group['temperature'])
    fin.close()

#the data have been read the analysis begins

#eigenvalues
    hyperb=np.zeros(len(name))  
    hyperb_plus=np.zeros(len(name))
    hyperb_minus=np.zeros(len(name))    
    for m in range(0,len(name)):
        D=np.array([[dx_ux[m],dx_uy[m],dx_uz[m]],[dy_ux[m],dy_uy[m],dy_uz[m]],[dz_ux[m],dz_uy[m],dz_uz[m]]])
        vals=LA.eigvals(D)
        iswhat=np.iscomplex(vals)
        if(iswhat[0]==True or iswhat[1]==True or iswhat[2]==True):
            hyperb[m] = 1.0
            for mm in range(0,3):
                if(iswhat[mm]!=True and vals[mm]>=0):
                    hyperb_plus[m] = 1.0
                if(iswhat[mm]!=True and vals[mm]<0):
                    hyperb_minus[m] = 1.0

#strain rate
    vecs1x=np.zeros(len(name))
    vecs1y=np.zeros(len(name))
    vecs1z=np.zeros(len(name))
    vecs2x=np.zeros(len(name))
    vecs2y=np.zeros(len(name))
    vecs2z=np.zeros(len(name))
    vecs3x=np.zeros(len(name))
    vecs3y=np.zeros(len(name))
    vecs3z=np.zeros(len(name))            
    for m in range(0,len(name)):
        s12=0.5*(dx_uy[m]+dy_ux[m])
        s13=0.5*(dx_uz[m]+dz_ux[m])
        s23=0.5*(dy_uz[m]+dz_uy[m])        
        S=np.array([[dx_ux[m],s12,s13],[s12,dy_uy[m],s23],[s13,s23,dz_uz[m]]])
        vals_unsort,vecs_unsort=LA.eigh(S)
        vals = np.sort(vals_unsort)
        vecs = vecs_unsort[:, vals_unsort.argsort()]
        vecs1x[m]=vecs[0,0]        
        vecs1y[m]=vecs[1,0]
        vecs1z[m]=vecs[2,0]        
        vecs2x[m]=vecs[0,1]
        vecs2y[m]=vecs[1,1]
        vecs2z[m]=vecs[2,1]
        vecs3x[m]=vecs[0,2]
        vecs3y[m]=vecs[1,2]
        vecs3z[m]=vecs[2,2]

        
# out of this loop
    p_times_vecs1 = px*vecs1x + py*vecs1y + pz*vecs1z
    p_times_vecs2 = px*vecs2x + py*vecs2y + pz*vecs2z
    p_times_vecs3 = px*vecs3x + py*vecs3y + pz*vecs3z        

    vecs1_norm = np.sqrt(np.power(vecs1x,2.0) + np.power(vecs1y,2.0)+ np.power(vecs1z,2.0))
    vecs2_norm = np.sqrt(np.power(vecs2x,2.0) + np.power(vecs2y,2.0)+ np.power(vecs2z,2.0))
    vecs3_norm = np.sqrt(np.power(vecs3x,2.0) + np.power(vecs3y,2.0)+ np.power(vecs3z,2.0))    
    
# gradients
    eps = ( (dx_ux + dx_ux)*(dx_ux + dx_ux) + (dx_uy + dy_ux)*(dx_uy + dy_ux) + (dx_uz + dz_ux)*(dx_uz + dz_ux) + (dy_ux + dx_uy)*(dy_ux + dy_uy) + (dy_uy + dy_uy)*(dy_uy + dy_uy) + (dy_uz + dz_uy)*(dy_uz + dz_uy) + (dz_ux + dx_uz)*(dz_ux + dx_uz) + (dz_uy + dz_uy)*(dz_uy + dz_uy) + (dz_uz + dz_uz)*(dz_uz + dz_uz) )/4.
   
    power_eps = eps

    wx = dy_uz - dz_uy;
    wy = dz_ux - dx_uz;
    wz = dx_uy - dy_ux;
    
    power_enst = (wx*wx+wy*wy+wz*wz)/4.    

# p cdot vorticity
    p_times_vorticity = px*wx + py*wy + pz*wz
    vorticity_norm = np.sqrt(np.power(wx,2.0) + np.power(wy,2.0)+ np.power(wz,2.0))

    vorticity_times_vecs1 = wx*vecs1x + wy*vecs1y + wz*vecs1z
    vorticity_times_vecs2 = wx*vecs2x + wy*vecs2y + wz*vecs2z
    vorticity_times_vecs3 = wx*vecs3x + wy*vecs3y + wz*vecs3z
            
#acceleration
    power_ax  = np.power(ax,order)
    power_ay  = np.power(ay,order)
    power_az  = np.power(az,order)

# p cdot gradt 
    p_times_gradt = px*dx_t + py*dy_t + pz*dz_t 
    power_p_times_gradt  = np.power(p_times_gradt,order)

    gradt_times_vecs1 = dx_t*vecs1x + dy_t*vecs1y + dz_t*vecs1z
    gradt_times_vecs2 = dx_t*vecs2x + dy_t*vecs2y + dz_t*vecs2z
    gradt_times_vecs3 = dx_t*vecs3x + dy_t*vecs3y + dz_t*vecs3z
            
# temperature gradient    
    gradt_norm = np.sqrt(np.power(dx_t,2.0) + np.power(dy_t,2.0)+ np.power(dz_t,2.0))
    power_gradt_norm  = np.power(gradt_norm,order)

# u2
    u2 = (np.power(ux,2.0) + np.power(uy,2.0)+ np.power(uz,2.0))
    power_u2  = np.power(u2,order)
    
# p norm 
    p_norm = np.sqrt(np.power(px,2.0) + np.power(py,2.0)+ np.power(pz,2.0))
    power_p_norm  = np.power(p_norm,order)

# dt_p norm
    dt_p = np.sqrt(np.power(dt_px,2.0) + np.power(dt_py,2.0)+ np.power(dt_pz,2.0))
    power_dt_p  = np.power(dt_p,order)    
    
#   temperature lagrangian derivative (unfiltered for the moment)
    dt_t = (tt_later - tt)/10.
    power_tt  = np.power(tt,order)
    power_dt_t = np.power(dt_t,order)

    totpart=len(name)
    print("totpart = "+str(totpart))
# sum up
    for ipart in range(0,totpart):
        itype = int(name[ipart])%ntype

        moment_hyperb[itype] += hyperb[ipart]
        moment_hyperb_plus[itype] += hyperb_plus[ipart]
        moment_hyperb_minus[itype] += hyperb_minus[ipart]
        
        moment_eps[itype] += power_eps[ipart]
        moment_enst[itype] += power_enst[ipart]

        moment_u2[itype] += power_u2[ipart]
        
        moment_ax[itype] += power_ax[ipart]
        moment_ay[itype] += power_ay[ipart]
        moment_az[itype] += power_az[ipart]
        
        moment_p_times_gradt[itype] += power_p_times_gradt[ipart] 

        moment_gradt[itype] += power_gradt_norm[ipart]

        moment_p[itype] += power_p_norm[ipart] 

        moment_dt_p[itype] += power_dt_p[ipart]

#   instantaneous correlation p \cdot \nabla T / ( p^2 \nabla T)  
        C[itype] += np.fabs(p_times_gradt[ipart])/(p_norm[ipart]*gradt_norm[ipart])

#   instantaneous correlation p \cdot vort / ( p^2 vort)
        Cvorticity[itype] += np.fabs(p_times_vorticity[ipart])/(p_norm[ipart]*vorticity_norm[ipart])

        Cpstrain1[itype] += np.fabs(p_times_vecs1[ipart])/(p_norm[ipart]*vecs1_norm[ipart])
        Cpstrain2[itype] += np.fabs(p_times_vecs2[ipart])/(p_norm[ipart]*vecs2_norm[ipart])
        Cpstrain3[itype] += np.fabs(p_times_vecs3[ipart])/(p_norm[ipart]*vecs3_norm[ipart])
        
        Cvorticitystrain1[itype] += np.fabs(vorticity_times_vecs1[ipart])/(vorticity_norm[ipart]*vecs1_norm[ipart])
        Cvorticitystrain2[itype] += np.fabs(vorticity_times_vecs2[ipart])/(vorticity_norm[ipart]*vecs2_norm[ipart])
        Cvorticitystrain3[itype] += np.fabs(vorticity_times_vecs3[ipart])/(vorticity_norm[ipart]*vecs3_norm[ipart])

        Cgradtstrain1[itype] += np.fabs(gradt_times_vecs1[ipart])/(gradt_norm[ipart]*vecs1_norm[ipart])
        Cgradtstrain2[itype] += np.fabs(gradt_times_vecs2[ipart])/(gradt_norm[ipart]*vecs2_norm[ipart])
        Cgradtstrain3[itype] += np.fabs(gradt_times_vecs3[ipart])/(gradt_norm[ipart]*vecs3_norm[ipart])
        
        moment_tt[itype] += power_tt[ipart]
        moment_dt_t[itype] += power_dt_t[ipart]  
        moment_t[itype] += tt[ipart]

# counter
        counter[itype] += 1 

    foutname = 'variance_particle_'+str(ntype)+'_'+str(order)+'.txt'
    fout = open(foutname,'w')
    for n in range(0,ntype):
        print(n,i,moment_tt[n]/counter[n], moment_dt_t[n]/counter[n], moment_gradt[n]/counter[n], moment_p_times_gradt[n]/counter[n], moment_p[n]/counter[n], moment_dt_p[n]/counter[n],C[n]/counter[n],moment_ax[n]/counter[n],moment_ay[n]/counter[n],moment_az[n]/counter[n],moment_eps[n]/counter[n],moment_enst[n]/counter[n],moment_u2[n]/counter[n],moment_hyperb[n]/counter[n], moment_hyperb_plus[n]/counter[n], moment_hyperb_minus[n]/counter[n],moment_t[n]/counter[n],Cvorticity[n]/counter[n],Cpstrain1[n]/counter[n],Cpstrain2[n]/counter[n],Cpstrain3[n]/counter[n],Cvorticitystrain1[n]/counter[n],Cvorticitystrain2[n]/counter[n],Cvorticitystrain3[n]/counter[n],Cgradtstrain1[n]/counter[n],Cgradtstrain2[n]/counter[n],Cgradtstrain3[n]/counter[n],file=fout)
    fout.close()



