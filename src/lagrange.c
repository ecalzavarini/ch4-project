#include "common_object.h"

#ifdef LAGRANGE


/* allocate particle containers */
void allocate_particles(){        

  int i;

  if((int)property.particle_number%nprocs ==0 ){

    /* ALL processor will take the same number of particles  */
    npart = property.particle_number/nprocs;


  }else{

   if(ROOT) fprintf(stderr,"Warning : total particle number is different from the process number!\n");

  for (i=0;i<nprocs;i++){

    /* ROOT processor will take just a little bit more particles */
    if(ROOT) npart = (int)property.particle_number - (nprocs-1)*(int)floor((int)property.particle_number/nprocs); else npart = (int)floor((int)property.particle_number /nprocs);  
					       
  }

  }/* end if on PARTICLE_NUMBER */

  fprintf(stderr,"me %d : I have %d particles\n",me, npart);

tracer  = (point_particle*) malloc(sizeof(point_particle)*npart);
if(tracer == NULL){ fprintf(stderr,"Not enough memory to allocate tracer\n"); exit(-1);}

tracer_here  = (point_particle*) malloc(sizeof(point_particle)*npart);
if(tracer_here == NULL){ fprintf(stderr,"Not enough memory to allocate tracer_here\n"); exit(-1);}

tracer_there  = (point_particle*) malloc(sizeof(point_particle)*npart);
if(tracer_there == NULL){ fprintf(stderr,"Not enough memory to allocate tracer_there\n"); exit(-1);}

all_tracer_there  = (point_particle*) malloc(sizeof(point_particle)*npart*nprocs);
if(all_tracer_there == NULL){ fprintf(stderr,"Not enough memory to allocate all_tracer_there\n"); exit(-1);}

}

/* initial conditions for particles */
void initial_conditions_particles(int restart){  

  int i;
  int *rcounts;
  int name_offset = 0;

    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);

    fprintf(stderr,"me : %d , name_offset %d\n", me , name_offset);


/* restart from file */
    if(restart){

      read_point_particle_h5();

    }else{
      /* restart from memory */

for (i=0;i<npart;i++) {

/* name */
(tracer+i)->name = i+name_offset;

/* position: randomly distributed particles */
(tracer+i)->x = LNX_START + drand48()*LNX;
(tracer+i)->y = LNY_START + drand48()*LNY;
(tracer+i)->z = LNZ_START + drand48()*LNZ;

/* velocity: null speed */
(tracer+i)->vx = 0.0;
(tracer+i)->vy = 0.0;
(tracer+i)->vz = 0.0;

 }

    }/* end of if /else on restart */


}/* end of function */



/* this function implement periodic or wall boundary conditions for the hydrodynamics fields 
   it is needed to correctly interpolate these fields at the particle positions */
void boundary_conditions_hydro(){

  int i,j,k;
  my_double fac, T_wall, S_wall;

  /* bc for the velocity field */
#ifdef LB_FLUID
sendrecv_borders_vector(u);

#ifdef LB_FLUID_BC_X
  for (j = 0; j < LNY + TWO_BRD; j++)                     
    for (k = 0; k < LNZ + TWO_BRD; k++){
        if(LNX_START == 0){
            i = BRD; 
            u[IDX(i-1,j,k)].x =  -u[IDX(i,j,k)].x;
            u[IDX(i-1,j,k)].y =  -u[IDX(i,j,k)].y;
            u[IDX(i-1,j,k)].z =  -u[IDX(i,j,k)].z;
        }
	if(LNX_END == NX){
            i = LNX+BRD-1;
            u[IDX(i+1,j,k)].x =  -u[IDX(i,j,k)].x;
            u[IDX(i+1,j,k)].y =  -u[IDX(i,j,k)].y;
            u[IDX(i+1,j,k)].z =  -u[IDX(i,j,k)].z;
        }
}
#endif

#ifdef LB_FLUID_BC_Y
  for (i = 0; i < LNX + TWO_BRD; i++)                     
    for (k = 0; k < LNZ + TWO_BRD; k++){
        if(LNY_START == 0){
            j = BRD; 
            u[IDX(i,j-1,k)].x =  -u[IDX(i,j,k)].x;
            u[IDX(i,j-1,k)].y =  -u[IDX(i,j,k)].y;
            u[IDX(i,j-1,k)].z =  -u[IDX(i,j,k)].z;
        }
	if(LNY_END == NY){
            j = LNY+BRD-1;
            u[IDX(i,j+1,k)].x =  -u[IDX(i,j,k)].x;
            u[IDX(i,j+1,k)].y =  -u[IDX(i,j,k)].y;
            u[IDX(i,j+1,k)].z =  -u[IDX(i,j,k)].z;
        }
}
#endif
#endif /* endif define LB_FLUID */

#ifdef LB_TEMPERATURE
sendrecv_borders_scalar(t);

#ifdef LB_TEMPERATURE_BC_Y
//my_double T_wall;

  for (i = 0; i < LNX + TWO_BRD; i++) 			
    for (k = 0; k < LNZ + TWO_BRD; k++){


if(LNY_START == 0){

	  j = BRD; 

	  T_wall = property.T_bot;

#ifdef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = 0.0;
#endif

	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  t[IDX(i,j-1,k)] =  fac*t[IDX(i,j,k)];
}

if(LNY_END == NY){

 	  j = LNY+BRD-1; 

	  T_wall = property.T_top;

#ifdef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = 0.0;
#endif
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  t[IDX(i,j+1,k)] =  fac*t[IDX(i,j,k)];
 }



      
}
#endif
#endif /* LB_TEMPERATURE */


#ifdef LB_SCALAR
sendrecv_borders_scalar(s);

#ifdef LB_SCALAR_BC_Y
// my_double S_wall;

  for (i = 0; i < LNX + TWO_BRD; i++) 			
    for (k = 0; k < LNZ + TWO_BRD; k++){


if(LNY_START == 0){

	  j = BRD; 

	  S_wall = property.S_bot;

#ifdef LB_SCALAR_FLUCTUATION 
	  S_wall = 0.0;
#endif

	  fac = 2.0*(S_wall-property.S_ref)/s[IDX(i,j,k)] - 1.0;
	  s[IDX(i,j-1,k)] =  fac*s[IDX(i,j,k)];
}

if(LNY_END == NY){

 	  j = LNY+BRD-1; 

	  S_wall = property.S_top;

#ifdef LB_SCALAR_FLUCTUATION 
	  S_wall = 0.0;
#endif
	  fac = 2.0*(S_wall-property.S_ref)/s[IDX(i,j,k)] - 1.0;
	  s[IDX(i,j+1,k)] =  fac*s[IDX(i,j,k)];
 }

      
    }
#endif
#endif /* end of SCALAR */


}


/* Interpolation for the moment only with regular grid */
void interpolate_vector_at_particles(vector *f,char which_vector){

 point_particle part;
 vector v;

 int ipart,im,jm,km,ip,jp,kp;
 double dxm,dxp,dym,dyp,dzm,dzp;
 double vol_ip_jp_kp,vol_im_jp_kp , vol_ip_jm_kp , vol_ip_jp_km , vol_im_jm_kp , vol_ip_jm_km , vol_im_jp_km , vol_im_jm_km; 

int i,j,k;


 for (ipart=0;ipart<npart;ipart++) {

   // fprintf(stderr,"vel %g %g\n", f[IDX(BRD-1,BRD,BRD)].x , f[IDX(BRD,BRD,BRD)].x );

  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

  //fprintf(stderr,"\n part %g %g %g\n",part.x,part.y,part.z);

for (i=0; i<LNX+TWO_BRD-1; i++) if(center_V[IDX(i, BRD, BRD)].x <= part.x && part.x < center_V[IDX(i+1,BRD, BRD)].x) im = i; 
ip =  im + 1;
for (j=0; j<LNY+TWO_BRD-1; j++) if(center_V[IDX(BRD, j, BRD)].y <= part.y && part.y < center_V[IDX(BRD, j+1, BRD)].y) jm = j;
jp =  jm + 1;
for (k=0; k<LNZ+TWO_BRD-1; k++) if(center_V[IDX(BRD, BRD, k)].z <= part.z && part.z < center_V[IDX(BRD, BRD, k+1)].z) km = k;
kp =  km + 1;

//fprintf(stderr,"index %d %d %d %d %d %d\n", im,ip,jm,jp,km,kp);

//for (j=0;j<10;j++) fprintf(stderr,"%d center_V %e\n",j,center_V[IDX(im, j, km)].y);

  dxm = part.x - center_V[IDX(im, BRD, BRD)].x;
  dxp = center_V[IDX(ip, BRD, BRD)].x - part.x;
  dym = part.y - center_V[IDX(BRD, jm, BRD)].y;
  dyp = center_V[IDX(BRD, jp, BRD)].y - part.y;
  dzm = part.z - center_V[IDX(BRD, BRD, km)].z;
  dzp = center_V[IDX(BRD, BRD, kp)].z - part.z;

  //fprintf(stderr,"distance %g %g %g %g %g %g\n", dxm,dxp,dym,dyp,dzm,dzp);

  vol_ip_jp_kp = dxp*dyp*dzp;
  vol_im_jp_kp = dxm*dyp*dzp;
  vol_ip_jm_kp = dxp*dym*dzp;
  vol_ip_jp_km = dxp*dyp*dzm;
  vol_im_jm_kp = dxm*dym*dzp;
  vol_ip_jm_km = dxp*dym*dzm;
  vol_im_jp_km = dxm*dyp*dzm;
  vol_im_jm_km = dxm*dym*dzm;


  v.x = f[IDX(im, jm, km)].x * vol_ip_jp_kp + 
        f[IDX(ip, jm, km)].x * vol_im_jp_kp +
        f[IDX(im, jp, km)].x * vol_ip_jm_kp +
        f[IDX(im, jm, kp)].x * vol_ip_jp_km +
        f[IDX(ip, jp, km)].x * vol_im_jm_kp +
        f[IDX(im, jp, kp)].x * vol_ip_jm_km +
        f[IDX(ip, jm, kp)].x * vol_im_jp_km +
        f[IDX(ip, jp, kp)].x * vol_im_jm_km ;

  v.y = f[IDX(im, jm, km)].y * vol_ip_jp_kp + 
        f[IDX(ip, jm, km)].y * vol_im_jp_kp +
        f[IDX(im, jp, km)].y * vol_ip_jm_kp +
        f[IDX(im, jm, kp)].y * vol_ip_jp_km +
        f[IDX(ip, jp, km)].y * vol_im_jm_kp +
        f[IDX(im, jp, kp)].y * vol_ip_jm_km +
        f[IDX(ip, jm, kp)].y * vol_im_jp_km +
        f[IDX(ip, jp, kp)].y * vol_im_jm_km ;

  v.z = f[IDX(im, jm, km)].z * vol_ip_jp_kp + 
        f[IDX(ip, jm, km)].z * vol_im_jp_kp +
        f[IDX(im, jp, km)].z * vol_ip_jm_kp +
        f[IDX(im, jm, kp)].z * vol_ip_jp_km +
        f[IDX(ip, jp, km)].z * vol_im_jm_kp +
        f[IDX(im, jp, kp)].z * vol_ip_jm_km +
        f[IDX(ip, jm, kp)].z * vol_im_jp_km +
        f[IDX(ip, jp, kp)].z * vol_im_jm_km ;

  /* if it is the velocity */
  if(which_vector == 'u'){
    (tracer+ipart)->ux = v.x;
    (tracer+ipart)->uy = v.y;
    (tracer+ipart)->uz = v.z; 
  }

 }/* end of for on ipart */

}


/* Interpolation for the moment only with regular grid */
void interpolate_scalar_at_particles(my_double *f,char which_scalar){

 point_particle part;
 my_double s;

 int ipart,im,jm,km,ip,jp,kp;
 double dxm,dxp,dym,dyp,dzm,dzp;
 double vol_ip_jp_kp,vol_im_jp_kp , vol_ip_jm_kp , vol_ip_jp_km , vol_im_jm_kp , vol_ip_jm_km , vol_im_jp_km , vol_im_jm_km; 

int i,j,k;

 for (ipart=0;ipart<npart;ipart++) {

   // fprintf(stderr,"vel %g %g\n", f[IDX(BRD-1,BRD,BRD)].x , f[IDX(BRD,BRD,BRD)].x );

  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

  //fprintf(stderr,"\n part %g %g %g\n",part.x,part.y,part.z);

for (i=0; i<LNX+TWO_BRD-1; i++) if(center_V[IDX(i, BRD, BRD)].x <= part.x && part.x < center_V[IDX(i+1,BRD, BRD)].x) im = i; 
ip =  im + 1;
for (j=0; j<LNY+TWO_BRD-1; j++) if(center_V[IDX(BRD, j, BRD)].y <= part.y && part.y < center_V[IDX(BRD, j+1, BRD)].y) jm = j;
jp =  jm + 1;
for (k=0; k<LNZ+TWO_BRD-1; k++) if(center_V[IDX(BRD, BRD, k)].z <= part.z && part.z < center_V[IDX(BRD, BRD, k+1)].z) km = k;
kp =  km + 1;

//fprintf(stderr,"index %d %d %d %d %d %d\n", im,ip,jm,jp,km,kp);

//for (j=0;j<10;j++) fprintf(stderr,"%d center_V %e\n",j,center_V[IDX(im, j, km)].y);

  dxm = part.x - center_V[IDX(im, BRD, BRD)].x;
  dxp = center_V[IDX(ip, BRD, BRD)].x - part.x;
  dym = part.y - center_V[IDX(BRD, jm, BRD)].y;
  dyp = center_V[IDX(BRD, jp, BRD)].y - part.y;
  dzm = part.z - center_V[IDX(BRD, BRD, km)].z;
  dzp = center_V[IDX(BRD, BRD, kp)].z - part.z;

  //fprintf(stderr,"distance %g %g %g %g %g %g\n", dxm,dxp,dym,dyp,dzm,dzp);

  vol_ip_jp_kp = dxp*dyp*dzp;
  vol_im_jp_kp = dxm*dyp*dzp;
  vol_ip_jm_kp = dxp*dym*dzp;
  vol_ip_jp_km = dxp*dyp*dzm;
  vol_im_jm_kp = dxm*dym*dzp;
  vol_ip_jm_km = dxp*dym*dzm;
  vol_im_jp_km = dxm*dyp*dzm;
  vol_im_jm_km = dxm*dym*dzm;


   s  = f[IDX(im, jm, km)] * vol_ip_jp_kp + 
        f[IDX(ip, jm, km)] * vol_im_jp_kp +
        f[IDX(im, jp, km)] * vol_ip_jm_kp +
        f[IDX(im, jm, kp)] * vol_ip_jp_km +
        f[IDX(ip, jp, km)] * vol_im_jm_kp +
        f[IDX(im, jp, kp)] * vol_ip_jm_km +
        f[IDX(ip, jm, kp)] * vol_im_jp_km +
        f[IDX(ip, jp, kp)] * vol_im_jm_km ;

#ifdef LB_TEMPERATURE
   /* if it is temperature */
if(which_scalar == 't')  (tracer+ipart)->t = s;
#endif

#ifdef LB_SCALAR
  /* if it is a scalar */
if(which_scalar == 's')  (tracer+ipart)->s = s;
#endif

  }

}/* end of interp scalar*/


#define H5FILE_NAME_PARTICLE "part.h5"

/* general output function for particles */
void output_particles(){
  int i,j;
  int np = (int)property.particle_number;
  FILE *fout;

#ifdef LAGRANGE_OUTPUT_DEBUG
    if(ROOT){
     for (i=0;i<npart;i++) {
      fprintf(stdout,"%g %e %e %e %e %e %e\n",time_now, (tracer+i)->x,(tracer+i)->y,(tracer+i)->z,(tracer+i)->vx,(tracer+i)->vy,(tracer+i)->vz);
     }
    }
#endif

#ifdef OUTPUT_H5

    int *rcounts;
    int name_offset = 0;

    /* First check how many particles in each processor and compute offset */
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);



    hid_t       file_id, dataset_id, dataspace_id , group ;  /* identifiers */
    hid_t	plist_id;                 /* property list identifier */
    hid_t hdf5_type;
    hid_t xfer_plist, ret, property_id;
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dims[1], offset[1], count[1];
    herr_t      hdf5_status;
    herr_t status;
    int size;    
    int RANK = 1;

    my_double *aux;

    char NEW_H5FILE_NAME[128];
    char XMF_FILE_NAME[128];

    /* then alloc space */
    aux  = (my_double*) malloc(sizeof(my_double)*npart); 
    if(aux == NULL){ fprintf(stderr,"Not enough memory to allocate aux field t\n"); exit(-1);}

    hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);

                /* Create a new file using default properties */
		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);   
		
		file_id = H5Fcreate(H5FILE_NAME_PARTICLE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		group   = H5Gcreate (file_id, "/lagrange", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

		H5Pclose(plist_id);

		property_id  = H5Pcreate(H5P_DATASET_CREATE);


                /* Create the data space for the dataset. */
                dims[0] = (int)property.particle_number;

		filespace = H5Screate_simple(RANK, dims, NULL);   
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
		count[0] = npart;
		offset[0] = name_offset;

	        memspace = H5Screate_simple(RANK, count, NULL); 	

    /*
     * Select hyperslab in the file.
     */
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

		xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);                       

		/* WRITE PARTICLE NAME */
		dataset_id = H5Dcreate(group, "name", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->name;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
                
		/* WRITE PARTICLE POSITIONS */
		dataset_id = H5Dcreate(group, "x", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->x;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "y", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->y;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "z", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->z;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		/* WRITE PARTICLE VELOCITIES */
		dataset_id = H5Dcreate(group, "vx", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->vx;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "vy", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->vy;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "vz", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->vz;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

#ifdef LB_TEMPERATURE
		/* WRITE PARTICLE TEMPERATURE */		
		dataset_id = H5Dcreate(group, "temperature", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->t;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
		
#endif

#ifdef LB_SCALAR
		/* WRITE PARTICLE SCALAR */
		dataset_id = H5Dcreate(group, "scalar", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->s;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
#endif

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

  /* free scalar auxiliary field */
  free(aux);

 /* create the file names */
  sprintf(NEW_H5FILE_NAME,"%s/particle_%d.h5",OutDir,itime);

  /* we rename the file */
  if(ROOT) rename(H5FILE_NAME_PARTICLE, NEW_H5FILE_NAME);


  /* Xml file */
  if(ROOT){
  sprintf(XMF_FILE_NAME,"%s/particle_%d.xmf" ,OutDir,itime);
  sprintf(NEW_H5FILE_NAME,"particle_%d.h5",itime);
  size=sizeof(my_double);
  

                fout = fopen(XMF_FILE_NAME,"w");
		              
                fprintf(fout,"<?xml version=\"1.0\" ?>\n");
                fprintf(fout,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                fprintf(fout,"<Xdmf Version=\"2.0\">\n");
                fprintf(fout,"<Domain>\n");
                fprintf(fout,"<Grid Name=\"Points\" GridType=\"Uniform\">\n");
                fprintf(fout,"<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n",np);
                fprintf(fout,"<Geometry GeometryType=\"X_Y_Z\">\n");
                
                fprintf(fout,"<DataItem Name=\"x\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/x\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                
                fprintf(fout,"<DataItem Name=\"y\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/y\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");

                fprintf(fout,"<DataItem Name=\"z\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/z\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");

                fprintf(fout,"</Geometry>\n");                

		/* velocity as a vector */
                fprintf(fout,"<Attribute Name=\"velocity\" AttributeType=\"Vector\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d 3\" \n   Function=\"JOIN($0 , $1, $2)\">\n",np);
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/vx\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/vy\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/vz\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  

#ifdef LB_TEMPERATURE
                fprintf(fout,"<Attribute Name=\"temperature\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/temperature\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");          
#endif

#ifdef LB_SCALAR
                fprintf(fout,"<Attribute Name=\"scalar\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/scalar\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");          
#endif
                
                fprintf(fout,"</Grid>\n");
                fprintf(fout,"</Domain>\n");
                fprintf(fout,"</Xdmf>\n");
		                
                fclose(fout);
  }/* end of if root */

#endif


}



/* advance in time particles and assign them to the right processors */
void move_particles(){

  int ipart,i;
  int npart_here,npart_there,all_npart_there,all_npart;  
  point_particle part;

  int *displs,*rcounts; 

  //fprintf(stderr,"me %d I am here, npart %d time %g\n",me, npart,time_now);

/* Begin loop on particles */
 for (ipart=0;ipart<npart;ipart++) {


   /* copy fluid velocity into particle velocity NOTE that this is true only for tracers */
   (tracer+ipart)->vx = (tracer+ipart)->ux;
   (tracer+ipart)->vy = (tracer+ipart)->uy;
   (tracer+ipart)->vz = (tracer+ipart)->uz;

   if(itime==0){
   /* Explicit Euler 1st order */
   (tracer+ipart)->x += property.time_dt*(tracer+ipart)->vx;
   (tracer+ipart)->y += property.time_dt*(tracer+ipart)->vy;
   (tracer+ipart)->z += property.time_dt*(tracer+ipart)->vz;
   }else{
   /* Adams-Bashforth 2nd order */
   (tracer+ipart)->x += property.time_dt*0.5*(3.0*(tracer+ipart)->vx - (tracer+ipart)->vx_old);
   (tracer+ipart)->y += property.time_dt*0.5*(3.0*(tracer+ipart)->vy - (tracer+ipart)->vy_old);
   (tracer+ipart)->z += property.time_dt*0.5*(3.0*(tracer+ipart)->vz - (tracer+ipart)->vz_old);
   }
   /* copy particle velocity in old */
   (tracer+ipart)->vx_old = (tracer+ipart)->vx; 
   (tracer+ipart)->vy_old = (tracer+ipart)->vy; 
   (tracer+ipart)->vz_old = (tracer+ipart)->vz;  

}/* end of loop on particles */


/* Here we perform the particle re-arrangement between processors */

 npart_here = 0;
 npart_there= 0;
 all_npart_there = 0;
 all_npart = 0;

for (ipart=0;ipart<npart;ipart++) {

  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

#ifdef LAGRANGE_WRAP /* this really makes particles to stay in the box */
  (tracer+ipart)->x = wrap( (tracer+ipart)->x ,  property.SX);
  (tracer+ipart)->y = wrap( (tracer+ipart)->y ,  property.SY);
  (tracer+ipart)->z = wrap( (tracer+ipart)->z ,  property.SZ); 
#endif


  //fprintf(stderr,"\n part %g %g %g\n",part.x,part.y,part.z);

  /* check how many particles are still in the local domain (here) and how many have to go (there) */

  /*
 if(  part.x >= center_V[IDX(0, BRD, BRD)].x && part.x < center_V[IDX(LNX+TWO_BRD-1,BRD, BRD)].x &&
      part.y >= center_V[IDX(BRD, 0, BRD)].y && part.y < center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y &&
      part.z >= center_V[IDX(BRD, BRD, 0)].z && part.z < center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z ){
  */
if(  part.x >= mesh[IDXG(BRD, BRD, BRD)].x && part.x < mesh[IDXG(LNXG+BRD-1,BRD, BRD)].x &&
     part.y >= mesh[IDXG(BRD, BRD, BRD)].y && part.y < mesh[IDXG(BRD,LNYG+BRD-1, BRD)].y &&
     part.z >= mesh[IDXG(BRD, BRD, BRD)].z && part.z < mesh[IDXG(BRD, BRD,LNZG+BRD-1)].z ){

      npart_here += 1;

      tracer_here  = (point_particle*) realloc( tracer_here, sizeof(point_particle)*npart_here);
      tracer_here[npart_here-1] = tracer[ipart];

      //fprintf(stderr,"Ehi! ipart %d\n",ipart);
  }else{

      npart_there += 1;   
      tracer_there  = (point_particle*) realloc(tracer_there,sizeof(point_particle)*npart_there);
      tracer_there[npart_there-1] = tracer[ipart];

      //fprintf(stderr,"AAA me %d , tracer_there[%d] %g,  tracer[%d] %g \n",me, npart_there-1, tracer_there[npart_there-1].x , ipart, tracer[ipart].x);

  }/* end of if else */

 }/* for loop on ipart */

//fprintf(stderr,"me %d : npart_here %d , npart_there %d\n",me,npart_here, npart_there);

/* first we communicate to other procs how many particles we have to give away and we sum up them among all the processors */

 MPI_Allreduce(&npart_there, &all_npart_there, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

#ifdef LAGRANGE_DEBUG
 if(ROOT)fprintf(stderr,"me %d : all_npart_there %d\n",me, all_npart_there);
#endif

 if(all_npart_there != 0){

    displs  = (int *)malloc(nprocs*sizeof(int)); 
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart_there, 1 , MPI_INT, rcounts, 1 , MPI_INT,MPI_COMM_WORLD);


    //displs[0]=0;
    //for (i=1;i<nprocs;i++) displs[i] += rcounts[i];
    //for (i=0;i<nprocs;i++) displs[i] = rcounts[i];
    displs[0] = 0;
    for (i=1; i<nprocs; i++) displs[i] = displs[i-1] + rcounts[i-1];


    //for (i=0;i<nprocs;i++) fprintf(stderr,"me %d : rcounts[%d] = %d\n",me,i, rcounts[i]);    
 
  /* space is allocate for coming particles */
  all_tracer_there = (point_particle*) realloc(all_tracer_there,sizeof(point_particle)*all_npart_there);
 
 /* Allgather to get all the migrant particles */
 MPI_Allgatherv(tracer_there, npart_there, MPI_point_particle_type, all_tracer_there, rcounts, displs, MPI_point_particle_type,MPI_COMM_WORLD);
 

 /* Begin loop on particles which just arrived */
 for (ipart=0;ipart<all_npart_there;ipart++) {

   //fprintf(stderr,"BBB me %d , part %g\n",me,(all_tracer_there+ipart)->x);

  part.x = wrap( (all_tracer_there+ipart)->x ,  property.SX);
  part.y = wrap( (all_tracer_there+ipart)->y ,  property.SY);
  part.z = wrap( (all_tracer_there+ipart)->z ,  property.SZ); 

  /* check how many particles are still in the local domain (here) and how many have to go (there) */
  /*
  if( part.x >= center_V[IDX(0, BRD, BRD)].x && part.x < center_V[IDX(LNX+TWO_BRD-1, BRD, BRD)].x &&
      part.y >= center_V[IDX(BRD, 0, BRD)].y && part.y < center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y &&
      part.z >= center_V[IDX(BRD, BRD, 0)].z && part.z < center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z ){
  */
  if(  part.x >= mesh[IDXG(BRD, BRD, BRD)].x && part.x < mesh[IDXG(LNXG+BRD-1,BRD, BRD)].x &&
       part.y >= mesh[IDXG(BRD, BRD, BRD)].y && part.y < mesh[IDXG(BRD,LNYG+BRD-1, BRD)].y &&
       part.z >= mesh[IDXG(BRD, BRD, BRD)].z && part.z < mesh[IDXG(BRD, BRD,LNZG+BRD-1)].z ){

      npart_here += 1;

      tracer_here  = (point_particle*) realloc(tracer_here,sizeof(point_particle)*npart_here);
      tracer_here[npart_here-1] = all_tracer_there[ipart];
  }/* end if */

 }/* for on ipart till all_npart_there */


 /* here I realloc tracer and copy the new data into it */

      tracer = (point_particle*) realloc(tracer,sizeof(point_particle)*npart_here);   

      for (ipart=0;ipart<npart_here;ipart++){  tracer[ipart] = tracer_here[ipart];}

       free(rcounts);
       free(displs);
 }/* end on if on all_npart_there */

      npart = npart_here;

      //fprintf(stderr,"me %d : new npart is %d\n",me, npart);


      /* for debug */
      /*
      fprintf(stderr," %g %g\n %g %g\n %g %g\n", center_V[IDX(0, BRD, BRD)].x , center_V[IDX(LNX+TWO_BRD-1,BRD, BRD)].x, 
	                                         center_V[IDX(BRD, 0 , BRD)].y , center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y,
	                                         center_V[IDX(BRD, BRD, 0)].z , center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z );
      */
      /*
     fprintf(stderr," %g %g\n %g %g\n %g %g\n", mesh[IDXG(BRD, BRD, BRD)].x , mesh[IDXG(LNXG+BRD-1,BRD, BRD)].x, 
	                                         mesh[IDXG(BRD, BRD , BRD)].y , mesh[IDXG(BRD,LNYG+BRD-1, BRD)].y,
	                                         mesh[IDXG(BRD, BRD, BRD)].z , mesh[IDXG(BRD, BRD,LNZG+BRD-1)].z );
      */
      /* final check */
       MPI_Allreduce(&npart, &all_npart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
       if(itime%10==0 && ROOT)fprintf(stderr,"------------------ Check npart = %d\n",all_npart);

       if(all_npart != (int)property.particle_number){
         if(ROOT) fprintf(stderr,"Total number of bubbles has changed during run!!!!\n Was %d now is %d\n Exit.\n",(int)property.particle_number, all_npart);
         MPI_Finalize();
         exit(0);
       }

 

}/* end of move particles */



/* here we have the function to write the full point_particle structure and to read it */

/* general output function for particles */

#ifdef OUTPUT_H5

void write_point_particle_h5(){
  int i,j;
  int np = (int)property.particle_number;
  FILE *fout;

    int *rcounts;
    int name_offset = 0;

    /* First check how many particles in each processor and compute offset */
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);


    hid_t       file_id, dataset_id, dataspace_id , group ;  /* identifiers */
    hid_t	plist_id;                 /* property list identifier */
    hid_t hdf5_type;
    hid_t xfer_plist, ret, property_id;
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dims[1], offset[1], count[1];
    herr_t      hdf5_status;
    herr_t status;
    int size;    
    int RANK = 1;

    my_double *aux;

    char NEW_H5FILE_NAME[128];
    char XMF_FILE_NAME[128];
 
    char label[128]; 

    /* create point particle compound */
    hdf5_type = H5Tcreate (H5T_COMPOUND, sizeof(point_particle));

    /* define its offsets */
    sprintf(label,"x");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, x), H5T_NATIVE_DOUBLE);
    sprintf(label,"y");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, y), H5T_NATIVE_DOUBLE);
    sprintf(label,"z");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, z), H5T_NATIVE_DOUBLE);

    sprintf(label,"vx");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vx), H5T_NATIVE_DOUBLE);
    sprintf(label,"vy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vy), H5T_NATIVE_DOUBLE);
    sprintf(label,"vz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vz), H5T_NATIVE_DOUBLE);

    sprintf(label,"vx_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vx_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"vy_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vy_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"vz_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vz_old), H5T_NATIVE_DOUBLE);

#ifdef LB_FLUID
    sprintf(label,"ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uz), H5T_NATIVE_DOUBLE);
#endif

#ifdef LB_TEMPERATURE
    sprintf(label,"t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, t), H5T_NATIVE_DOUBLE);
#endif

#ifdef LB_SCALAR
    sprintf(label,"s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, s), H5T_NATIVE_DOUBLE);
#endif

    /*************************************************************/
  
                /* Create a new file using default properties */
		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);   
		
		file_id = H5Fcreate(H5FILE_NAME_PARTICLE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		group   = H5Gcreate (file_id, "/particles", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

		H5Pclose(plist_id);

		property_id  = H5Pcreate(H5P_DATASET_CREATE);


                /* Create the data space for the dataset. */
                dims[0] = (int)property.particle_number;

		filespace = H5Screate_simple(RANK, dims, NULL);   
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
		count[0] = npart;
		offset[0] = name_offset;

	        memspace = H5Screate_simple(RANK, count, NULL); 	

    /*
     * Select hyperslab in the file.
     */
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

		xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);                       

		/* WRITE POINT_PARTICLE STRUCTURE */
		dataset_id = H5Dcreate(group, "point_particle", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, tracer);
                status = H5Dclose(dataset_id);                

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

 /* create the file names */
  sprintf(NEW_H5FILE_NAME,"part.h5",itime);

  /* we rename the file */
  if(ROOT) rename(H5FILE_NAME_PARTICLE, NEW_H5FILE_NAME);

}/* end of write point particle */


/* Now the function to read point particles */
void read_point_particle_h5(){
  int i,j;
  int np = (int)property.particle_number;
  FILE *fout;

    int *rcounts;
    int name_offset = 0;

    /* First check how many particles in each processor and compute offset */
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);


    hid_t       file_id, dataset_id, dataspace_id , group ;  /* identifiers */
    hid_t	plist_id;                 /* property list identifier */
    hid_t hdf5_type;
    hid_t xfer_plist, ret, property_id;
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dims[1], offset[1], count[1];
    herr_t      hdf5_status;
    herr_t status;
    int size;    
    int RANK = 1;

    my_double *aux;

    char NEW_H5FILE_NAME[128];
    char XMF_FILE_NAME[128];
 
    char label[128]; 

    /* create point particle compound */
    hdf5_type = H5Tcreate (H5T_COMPOUND, sizeof(point_particle));

    /* define its offsets */
    sprintf(label,"x");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, x), H5T_NATIVE_DOUBLE);
    sprintf(label,"y");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, y), H5T_NATIVE_DOUBLE);
    sprintf(label,"z");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, z), H5T_NATIVE_DOUBLE);

    sprintf(label,"vx");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vx), H5T_NATIVE_DOUBLE);
    sprintf(label,"vy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vy), H5T_NATIVE_DOUBLE);
    sprintf(label,"vz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vz), H5T_NATIVE_DOUBLE);

    sprintf(label,"vx_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vx_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"vy_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vy_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"vz_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vz_old), H5T_NATIVE_DOUBLE);

#ifdef LB_FLUID
    sprintf(label,"ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uz), H5T_NATIVE_DOUBLE);
#endif

#ifdef LB_TEMPERATURE
    sprintf(label,"t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, t), H5T_NATIVE_DOUBLE);
#endif

#ifdef LB_SCALAR
    sprintf(label,"s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, s), H5T_NATIVE_DOUBLE);
#endif

    /*************************************************************/
  
                /* Create a new file using default properties */
		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);   
		
		file_id = H5Fopen(H5FILE_NAME_PARTICLE, H5F_ACC_RDONLY, H5P_DEFAULT); 
		group   = H5Gopen(file_id, "/particles", H5P_DEFAULT);

		H5Pclose(plist_id);

		property_id  = H5Pcreate(H5P_DATASET_CREATE);


                /* Create the data space for the dataset. */
                dims[0] = (int)property.particle_number;

		filespace = H5Screate_simple(RANK, dims, NULL);   
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
		count[0] = npart;
		offset[0] = name_offset;

	        memspace = H5Screate_simple(RANK, count, NULL); 	

    /*
     * Select hyperslab in the file.
     */
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

		xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);                       

		/* WRITE POINT_PARTICLE STRUCTURE */
		dataset_id = H5Dopen(group, "point_particle", H5P_DEFAULT);
                ret = H5Dread(dataset_id, hdf5_type, memspace, filespace,  H5P_DEFAULT, tracer);
                status = H5Dclose(dataset_id);                

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

}/* end of write point particle */


#endif


#endif
