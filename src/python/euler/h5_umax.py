# multiplying the velocity field read in pop.h5 by umax_new/umax_old
# call: python h5_umax.py 'pop.h5' 0.4 0.02
#!/usr/local/bin/python
import h5py
import numpy as np
import sys
import math
import matplotlib.pyplot as plt

spop = sys.argv[1]
umax_old = float(sys.argv[2])
umax_new = float(sys.argv[3])
vscl = umax_new / umax_old

f = h5py.File(spop, 'r')
t=np.array(f['population']['temperature'])
v=np.array(f['population']['velocity'])

NZ = v.shape[0]
NY = v.shape[1]
NX = v.shape[2]

new_v = np.zeros((NZ,NY,NX),v.dtype)
rho = np.zeros((NZ,NY,NX))
unorm = np.zeros((NZ,NY,NX))
mom = np.zeros((NZ,NY,NX,3))
c = np.zeros((3,9))
w = np.zeros(9)
u = np.zeros(3)

# must be consistent with the convention used in lb.c, GRID_POP_D2Q9
c[2,0]= 0.0; c[1,0]= 0.0; c[0,0]= 0.0;
c[2,1]=-1.0; c[1,1]= 1.0; c[0,1]= 0.0;
c[2,2]=-1.0; c[1,2]= 0.0; c[0,2]= 0.0;
c[2,3]=-1.0; c[1,3]=-1.0; c[0,3]= 0.0;
c[2,4]= 0.0; c[1,4]=-1.0; c[0,4]= 0.0;
c[2,5]= 1.0; c[1,5]=-1.0; c[0,5]= 0.0;
c[2,6]= 1.0; c[1,6]= 0.0; c[0,6]= 0.0;
c[2,7]= 1.0; c[1,7]= 1.0; c[0,7]= 0.0;
c[2,8]= 0.0; c[1,8]= 1.0; c[0,8]= 0.0;
w[0] = 4. / 9.;
w[1] = 1. / 36.;
w[2] = 1. / 9.;
w[3] = 1. / 36.;
w[4] = 1. / 9.;
w[5] = 1. / 36.;
w[6] = 1. / 9.;
w[7] = 1. / 36.;
w[8] = 1. / 9.;
css=1./3.

for i in range(0,NX):
    for j in range(0,NY):
        for k in range(0,NZ):
            for l in range(0,9):
                rho[k,j,i] = rho[k,j,i] + v[k,j,i]['p'+str(l)]
                mom[k,j,i,:] = mom[k,j,i,:] + c[:,l]*v[k,j,i]['p'+str(l)]
                
for i in range(0,NX):
    for j in range(0,NY):
        for k in range(0,NZ):            
            u[:] = vscl * mom[k,j,i,:] / rho[k,j,i]
            # unorm[k,j,i] = math.sqrt(mom[k,j,i,1]**2 + mom[k,j,i,2]**2) / rho[k,j,i]
            for l in range(0,9):                
                new_v[k,j,i]['p'+str(l)] = w[l]*rho[k,j,i]*(1. + np.dot(u,c[:,l])/css + (np.dot(u,c[:,l]))**2/(2.*css**2) - np.dot(u,u)/css)

new_f = h5py.File('newumax_pop.h5', 'w')
grp = new_f.create_group('population')
dset_t = grp.create_dataset('temperature',(NZ,NY,NX), t.dtype)
dset_v = grp.create_dataset('velocity',(NZ,NY,NX), new_v.dtype)
dset_t[:] = t
dset_v[:] = new_v

#THE TWO LOOPS ABOVE DO THE JOB, the two ones below are just to look at the new velocity field directly here in the python script

#rhon = np.zeros((NZ,NY,NX))
#momn = np.zeros((NZ,NY,NX,3))
#unormn = np.zeros((NZ,NY,NX))
#for i in range(0,NX):
#    for j in range(0,NY):
#        for k in range(0,NZ):
#            for l in range(0,9):
#                rhon[k,j,i] = rhon[k,j,i] + new_v[k,j,i]['p'+str(l)]
#                momn[k,j,i,:] = momn[k,j,i,:] + c[:,l]*new_v[k,j,i]['p'+str(l)]
                
#for i in range(0,NX):
#    for j in range(0,NY):
#        for k in range(0,NZ):            
#            unormn[k,j,i] = math.sqrt(momn[k,j,i,1]**2 + momn[k,j,i,2]**2) / rhon[k,j,i]

#VISUALIZING RESULTS

#fig_size = plt.rcParams["figure.figsize"]  
#fig_size[0] = 10  
#fig_size[1] = 8  
#plt.rcParams["figure.figsize"] = fig_size

#plt.imshow(unorm[0,:,:], cmap='coolwarm', interpolation='nearest', origin='lower')
#plt.imshow(rho[0,:,:], cmap='coolwarm', interpolation='nearest', origin='lower')
#plt.imshow(mom[0,:,:,2]/rho[0,:,:], cmap='coolwarm', interpolation='nearest', origin='lower')
#plt.colorbar()
#plt.savefig('v_x.png')

#plt.imshow(momn[0,:,:,2]/rhon[0,:,:], cmap='coolwarm', interpolation='nearest', origin='lower')
#plt.savefig('new_v_x.png')

#NX=50; NY=50;
#X = np.arange(0, NX, 1)
#Y = np.arange(0, NY, 1)
#fig, ax = plt.subplots()
#q = ax.quiver(X, Y, mom[0,:NY,:NX,2], mom[0,:NY,:NX,1])
#plt.savefig('v.png')
#q = ax.quiver(X, Y, momn[0,:NY,:NX,2], momn[0,:NY,:NX,1])
#plt.savefig('new_v.png')

#plt.show()

