#include "common_object.h"

#ifdef LAGRANGE

/*
NPART = total number of particles in the simulation 

NPART_PROC = NPART/nprocs
*/

#define PARTICLE_NUMBER 100



/* allocate particle containers */
void allocate_particles(){        

  int i;

  if(PARTICLE_NUMBER%nprocs ==0 ){

    /* ALL processor will take the same number of particles  */
    npart = PARTICLE_NUMBER/nprocs;


  }else{

   if(ROOT) fprintf(stderr,"Warning : total particle number is different from the process number!\n");

  for (i=0;i<nprocs;i++){

    /* ROOT processor will take just a little bit more particles */
    if(ROOT) npart = PARTICLE_NUMBER - (nprocs-1)*(int)floor(PARTICLE_NUMBER/nprocs); else npart = (int)floor(PARTICLE_NUMBER/nprocs);  
					       
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
void initial_conditions_particles(){  

  int i;
  int *rcounts;
  int name_offset = 0;

    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);

    fprintf(stderr,"me : %d , name_offset %d\n", me , name_offset);

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

}



/* this function implement periodic or wall boundary conditions for the hydrodynamics fields 
   it is needed to correctly interpolate these fields at the particle positions */
void boundary_conditions_hydro(){

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
            f[IDX(i,j-1,k)].x =  -u[IDX(i,j,k)].x;
            f[IDX(i,j-1,k)].y =  -u[IDX(i,j,k)].y;
            f[IDX(i,j-1,k)].z =  -u[IDX(i,j,k)].z;
        }
	if(LNY_END == NY){
            j = LNY+BRD-1;
            f[IDX(i,j+1,k)].x =  -u[IDX(i,j,k)].x;
            f[IDX(i,j+1,k)].y =  -u[IDX(i,j,k)].y;
            f[IDX(i,j+1,k)].z =  -u[IDX(i,j,k)].z;
        }
}
#endif
#endif /* endif define LB_FLUID */

#ifdef LB_TEMPERATURE
sendrecv_borders_vector(t);

#endif


#ifdef LB_SCALAR
sendrecv_borders_vector(t);

#endif


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
    (tracer+ipart)->vx = v.x;
    (tracer+ipart)->vy = v.y;
    (tracer+ipart)->vz = v.z; 
  }

 }/* end of for on ipart */

}


#ifdef LB_TEMPERATURE
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

   /* if it is temperature */
if(which_scalar == 't')  (tracer+ipart)->t = s;

  /* if it is a scalar */
if(which_scalar == 's')  (tracer+ipart)->s = s;


  }

}
#endif

/* general output function for particles */
void output_particles(){
  int i,j;
  FILE *fout;

    if(ROOT){
     for (i=0;i<npart;i++) {
      fprintf(stdout,"%g %e %e %e %e %e %e\n",time_now, (tracer+i)->x,(tracer+i)->y,(tracer+i)->z,(tracer+i)->vx,(tracer+i)->vy,(tracer+i)->vz);
     }
    }


  int *rcounts;
  int name_offset = 0;

    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);


#ifdef OUTPUT_H5_AAA
    hid_t       file_id, dataset_id, dataspace_id;  /* identifiers */
    hsize_t     dims[1], offset, count;
    herr_t      status;
    hid_t hdf5_type;
    my_double *aux;

    FILE *fout;
    double * dset_data = (double *) malloc (npart*sizeof(double));
    char NEW_H5FILE_NAME[128];
    char XMF_FILE_NAME[128];


    aux  = (my_double*) malloc(sizeof(my_double)*npart); 
    if(aux == NULL){ fprintf(stderr,"Not enough memory to allocate aux field t\n"); exit(-1);}

    hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);

    /* create the file names */
                sprintf(NEW_H5FILE_NAME,"%s/particles_%d.h5",OutDir,itime);
                sprintf(XMF_FILE_NAME,"%s/particles_%d.xmf",OutDir,itime);

                /* Create a new file using default properties. */
		plist_id = H5Pcreate(H5P_FILE_ACCESS);

		hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

		file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		group   = H5Gcreate (file_id, "/lagrange", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

		H5Pclose(plist_id);
                
                if( file_id < 0 )
                {
                    fprintf(stderr,"Error file_id while creating the file at timestep%d\n",itime);
                    break;
                }


		xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);

                /* Create the data space for the dataset. */
                dims[0] = PARTICLE_NUMBER;

		//filespace = H5Screate_simple(RANK, dims, NULL);
    

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
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

                                          
                dataspace_id = H5Screate_simple(1, dims, NULL);
                                
                dataset_id = H5Dcreate2(file_id, "name", hdf5_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		for(i=0;i<npart;i++) aux[i]=(particle + i)->name;

                //status = H5Dwrite(dataset_id, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, aux);               
                status = H5Dwrite(dataset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, aux); 
               /* End access to the dataset and release resources used by it. */
                status = H5Dclose(dataset_id);
                
                /* Terminate access to the data space. */
                status = H5Sclose(dataspace_id);
                
                /* Close the file. */
                status = H5Fclose(file_id);
                
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

   /* Adams-Bashforth 2nd order */
   ///*
   (tracer+ipart)->x += property.time_dt*0.5*(3.0*(tracer+ipart)->vx - (tracer+ipart)->vx_old);
   (tracer+ipart)->y += property.time_dt*0.5*(3.0*(tracer+ipart)->vy - (tracer+ipart)->vy_old);
   (tracer+ipart)->z += property.time_dt*0.5*(3.0*(tracer+ipart)->vz - (tracer+ipart)->vz_old);
   
   (tracer+ipart)->vx_old = (tracer+ipart)->vx; 
   (tracer+ipart)->vy_old = (tracer+ipart)->vy; 
   (tracer+ipart)->vz_old = (tracer+ipart)->vz; 
   //*/

}/* end of loop on particles */


/* Here we perform the particle rearrangement between processors */

 npart_here = 0;
 npart_there= 0;
 all_npart_there = 0;
 all_npart = 0;

for (ipart=0;ipart<npart;ipart++) {

  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

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

 if(ROOT)fprintf(stderr,"me %d : all_npart_there %d\n",me, all_npart_there);

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
       if(ROOT)fprintf(stderr,"------------------ Check npart = %d\n",all_npart);


 

}/* end of move particles */


#endif
