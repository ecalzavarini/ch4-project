#include "common_object.h"

#ifdef LAGRANGE

/*
NPART = total number of particles in the simulation 

NPART_PROC = NPART/nprocs
*/

#define PARTICLE_NUMBER 100



/* allocate particle containers */
void allocate_particles(){        

npart = PARTICLE_NUMBER/nprocs; 

 fprintf(stderr,"We have %d particles\n",npart);

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

for (i=0;i<npart;i++) {

/* name */
(tracer+i)->name = i+me*npart;

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



/* Interpolation for the moment only with regular grid */
void interpolate_vector_at_particles(vector *f){

 point_particle part;
 vector v;

 int ipart,im,jm,km,ip,jp,kp;
 double dxm,dxp,dym,dyp,dzm,dzp;
 double vol_ip_jp_kp,vol_im_jp_kp , vol_ip_jm_kp , vol_ip_jp_km , vol_im_jm_kp , vol_ip_jm_km , vol_im_jp_km , vol_im_jm_km; 

int i,j,k;
sendrecv_borders_vector(f);

#ifdef LB_FLUID_BC_X
  for (j = 0; j < LNY + TWO_BRD; j++)                     
    for (k = 0; k < LNZ + TWO_BRD; k++){
        if(LNX_START == 0){
            i = BRD; 
            f[IDX(i-1,j,k)].x =  -f[IDX(i,j,k)].x;
            f[IDX(i-1,j,k)].y =  -f[IDX(i,j,k)].y;
            f[IDX(i-1,j,k)].z =  -f[IDX(i,j,k)].z;
        }
	if(LNX_END == NX){
            i = LNX+BRD-1;
            f[IDX(i+1,j,k)].x =  -f[IDX(i,j,k)].x;
            f[IDX(i+1,j,k)].y =  -f[IDX(i,j,k)].y;
            f[IDX(i+1,j,k)].z =  -f[IDX(i,j,k)].z;
        }
}
#endif

#ifdef LB_FLUID_BC_Y
  for (i = 0; i < LNX + TWO_BRD; i++)                     
    for (k = 0; k < LNZ + TWO_BRD; k++){
        if(LNY_START == 0){
            j = BRD; 
            f[IDX(i,j-1,k)].x =  -f[IDX(i,j,k)].x;
            f[IDX(i,j-1,k)].y =  -f[IDX(i,j,k)].y;
            f[IDX(i,j-1,k)].z =  -f[IDX(i,j,k)].z;
        }
	if(LNY_END == NY){
            j = LNY+BRD-1;
            f[IDX(i,j+1,k)].x =  -f[IDX(i,j,k)].x;
            f[IDX(i,j+1,k)].y =  -f[IDX(i,j,k)].y;
            f[IDX(i,j+1,k)].z =  -f[IDX(i,j,k)].z;
        }
}
#endif

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


(tracer+ipart)->vx = v.x;
(tracer+ipart)->vy = v.y;
(tracer+ipart)->vz = v.z; 

  }

}


#ifdef LB_TEMPERATURE
/* Interpolation for the moment only with regular grid */
void interpolate_temperature_at_particles(my_double *f){

 point_particle part;
 my_double s;

 int ipart,im,jm,km,ip,jp,kp;
 double dxm,dxp,dym,dyp,dzm,dzp;
 double vol_ip_jp_kp,vol_im_jp_kp , vol_ip_jm_kp , vol_ip_jp_km , vol_im_jm_kp , vol_ip_jm_km , vol_im_jp_km , vol_im_jm_km; 

int i,j,k;
sendrecv_borders_vector(f);

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

(tracer+ipart)->t = s;

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

}



/* advance in time particles and assign them to the right processors */
void move_particles(){

  int ipart,i;
  int npart_here,npart_there,all_npart_there;  
  point_particle part;

  int *displs,*rcounts; 

  fprintf(stderr,"I am here, npart %d\n",npart);

/* Begin loop on particles */
 for (ipart=0;ipart<npart;ipart++) {

   /* Adams-Bashforth 2nd order */
   
   (tracer+ipart)->x += property.time_dt*0.5*(3.0*(tracer+ipart)->vx - (tracer+ipart)->vx_old);
   (tracer+ipart)->y += property.time_dt*0.5*(3.0*(tracer+ipart)->vy - (tracer+ipart)->vy_old);
   (tracer+ipart)->z += property.time_dt*0.5*(3.0*(tracer+ipart)->vz - (tracer+ipart)->vz_old);
   
   (tracer+ipart)->vx_old = (tracer+ipart)->vx; 
   (tracer+ipart)->vy_old = (tracer+ipart)->vy; 
   (tracer+ipart)->vz_old = (tracer+ipart)->vz; 
   

}/* end of loop on particles */


/* Here we perform the particle rearrangement between processors */

 npart_here = 0;
 npart_there= 0;
 all_npart_there = 0;

for (ipart=0;ipart<npart;ipart++) {

  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

  //fprintf(stderr,"\n part %g %g %g\n",part.x,part.y,part.z);

  /* check how many particles are still in the local domain (here) and how many have to go (there) */

 if(  part.x >= center_V[IDX(0, BRD, BRD)].x && part.x < center_V[IDX(LNX+TWO_BRD-1,BRD, BRD)].x &&
      part.y >= center_V[IDX(BRD, 0, BRD)].y && part.y < center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y &&
      part.z >= center_V[IDX(BRD, BRD, 0)].z && part.z < center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z ){

      npart_here += 1;

      tracer_here  = (point_particle*) realloc( tracer_here, sizeof(point_particle)*npart_here);
      tracer_here[npart_here-1] = tracer[ipart];

      //fprintf(stderr,"Ehi! ipart %d\n",ipart);
  }else{

      npart_there += 1;   
      tracer_there  = (point_particle*) realloc(tracer_there,sizeof(point_particle)*npart_there);
      tracer_there[npart_there-1] = tracer[ipart];

  }/* end of if else */

 }/* for loop on ipart */

 fprintf(stderr,"me %d : npart_here %d , npart_there %d\n",me,npart_here, npart_there);

/* first we communicate to other procs how many particles we have to give away and we sum up them among all the processors */

 MPI_Allreduce(&npart_there, &all_npart_there, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

 fprintf(stderr,"me %d : all_npart_there %d\n",me, all_npart_there);

 if(all_npart_there != 0){
    displs = (int *)malloc(nprocs*sizeof(int)); 
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart_there, 1 , MPI_INT, rcounts, 1 , MPI_INT,MPI_COMM_WORLD);

    for (i=0;i<nprocs;i++) displs[i]=0.0;

    for (i=0;i<nprocs;i++) fprintf(stderr,"me %d : rcounts[%d] = %d\n",me,i, rcounts[i]);    
 }
 /* space is allocate for coming particles */
 all_tracer_there = (point_particle*) realloc(all_tracer_there,sizeof(point_particle)*all_npart_there);
 
 /* Allgather to get all the migrant particles */
 MPI_Allgatherv(tracer_there, npart_there, MPI_point_particle_type, all_tracer_there, rcounts,displs, MPI_point_particle_type,MPI_COMM_WORLD);


/* Begin loop on particles which just arrived */
 for (ipart=0;ipart<all_npart_there;ipart++) {

  part.x = wrap( (all_tracer_there+ipart)->x ,  property.SX);
  part.y = wrap( (all_tracer_there+ipart)->y ,  property.SY);
  part.z = wrap( (all_tracer_there+ipart)->z ,  property.SZ); 

  /* check how many particles are still in the local domain (here) and how many have to go (there) */
  if( part.x >= center_V[IDX(0, BRD, BRD)].x && part.x < center_V[IDX(LNX+TWO_BRD-1, BRD, BRD)].x &&
      part.y >= center_V[IDX(BRD, 0, BRD)].y && part.y < center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y &&
      part.z >= center_V[IDX(BRD, BRD, 0)].z && part.z < center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z ){

      npart_here += 1;

      tracer_here  = (point_particle*) realloc(tracer_here,sizeof(point_particle)*npart_here);
      tracer_here[npart_here-1] = all_tracer_there[ipart];
  }/* end if */

 }/* for on ipart till all_npart_there */


 /* here I realloc tracer and copy the new data into it */

      tracer = (point_particle*) realloc(tracer,sizeof(point_particle)*npart_here);   

      for (ipart=0;ipart<npart_here;ipart++){  tracer[ipart] = tracer_here[ipart];}

      npart = npart_here;

      fprintf(stderr,"me %d : new npart is %d\n",me, npart);

      /* for debug */
      /*
      fprintf(stderr," %g %g\n %g %g\n %g %g\n", center_V[IDX(0, BRD, BRD)].x , center_V[IDX(LNX+TWO_BRD-1,BRD, BRD)].x, 
	                                         center_V[IDX(BRD, 0 , BRD)].y , center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y,
	                                         center_V[IDX(BRD, BRD, 0)].z , center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z );
      */

}/* end of move particles */


#endif
