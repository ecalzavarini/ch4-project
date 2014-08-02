#include "common_object.h"

#ifdef LAGRANGE

/*
NPART = total number of particles in the simulation 

NPART_PROC = NPART/nprocs
*/

#define NPART_PROC 1000



/* allocate particle containers */
void allocate_particles(){        

//NPART_PROC = NPART/nprocs; 

npart = NPART_PROC;
/*
npart_new = 0;
npart_xp  = 0;
npart_xm  = 0;
npart_yp  = 0;
npart_ym  = 0;
npart_zp  = 0;
npart_zm  = 0;
*/

tracer  = (point_particle*) malloc(sizeof(point_particle)*npart);
if(tracer == NULL){ fprintf(stderr,"Not enough memory to allocate tracer\n"); exit(-1);}

/*
tracer_new  = (point_particle*) malloc(sizeof(point_particle)*npart_new);
tracer_xp  = (point_particle*) malloc(sizeof(point_particle)*npart_xp);
tracer_xm  = (point_particle*) malloc(sizeof(point_particle)*npart_xm);
tracer_yp  = (point_particle*) malloc(sizeof(point_particle)*npart_yp);
tracer_ym  = (point_particle*) malloc(sizeof(point_particle)*npart_ym);
tracer_zp  = (point_particle*) malloc(sizeof(point_particle)*npart_zp);
tracer_zm  = (point_particle*) malloc(sizeof(point_particle)*npart_zm);
*/

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


/* general output function for particles */
void output_particles(){
  int i,j;

  for (j=0;j<nprocs;j++){
    //if(0==me){
    for (i=0;i<npart;i++) {
      if( (tracer+i)->name == 0 )
	fprintf(stdout,"%g %e %e %e %e %e %e\n",time_now, (tracer+i)->x,(tracer+i)->y,(tracer+i)->z,(tracer+i)->vx,(tracer+i)->vy,(tracer+i)->vz);

     }
    //}
  }

}



/* advance in time particles and assign them to the right processors */
void move_particles(){

  int ipart;

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




/* here we perform the MPI send recv*/
}













#endif
