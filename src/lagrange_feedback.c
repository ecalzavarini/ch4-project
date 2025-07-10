#include "common_object.h"

#ifdef LAGRANGE_TWOWAY
/* Extrapolation only with regular grid of unit spacing */
/* The form of the feedback is 
F_p = \Sigma_i^Npart * [ (D_t u - g) + (rho_p/rho_f) (g - a_p) ]

Also written as
F_p = \Sigma_i^Npart * [ (D_t u - rho_p/rho_f a_p) + (rho_p/rho_f - 1) g  ]
*/
void add_particle_feedbacks(){
 my_double fac, fac2, sp;
 vector fp , part , gravity_vector;
 my_double density_ratio;
 
 int ipart,im,jm,km,ip,jp,kp;
 double dxm,dxp,dym,dyp,dzm,dzp;
 double vol_ip_jp_kp,vol_im_jp_kp , vol_ip_jm_kp , vol_ip_jp_km , vol_im_jm_kp , vol_ip_jm_km , vol_im_jp_km , vol_im_jm_km; 
 int i,j,k;
 double f_im_jm_km , f_ip_jm_km , f_im_jp_km , f_im_jm_kp , f_ip_jp_km , f_im_jp_kp , f_ip_jm_kp , f_ip_jp_kp; 
 double dx_jm_km , dx_jp_kp , dx_jp_km , dx_jm_kp; 
 double dy_im_km , dy_ip_kp , dy_ip_km , dy_im_kp; 
 double dz_im_jm , dz_ip_jp , dz_ip_jm , dz_im_jp; 

#ifdef LAGRANGE_TWOWAY_MOMENTUM
 set_to_zero_vector( force_twoway, (LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD) );
#endif
#ifdef LAGRANGE_TWOWAY_TEMPERATURE
 set_to_zero_scalar( t_source_twoway, (LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD) );
#endif

 for (ipart=0;ipart<npart;ipart++) {

  /* get coordinates in the domain */
  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

  /* get indexes of neighboring nodes */
  for (i=0; i<LNX+TWO_BRD-1; i++) if(center_V[IDX(i, BRD, BRD)].x <= part.x && part.x < center_V[IDX(i+1,BRD, BRD)].x){ im = i;} 
  ip =  im + 1;
  for (j=0; j<LNY+TWO_BRD-1; j++) if(center_V[IDX(BRD, j, BRD)].y <= part.y && part.y < center_V[IDX(BRD, j+1, BRD)].y){ jm = j;}
  jp =  jm + 1;
  for (k=0; k<LNZ+TWO_BRD-1; k++) if(center_V[IDX(BRD, BRD, k)].z <= part.z && part.z < center_V[IDX(BRD, BRD, k+1)].z){ km = k;}
  kp =  km + 1;

  /* compute difference segments */
  dxm = part.x - center_V[IDX(im, BRD, BRD)].x;
  dxp = center_V[IDX(ip, BRD, BRD)].x - part.x;
  dym = part.y - center_V[IDX(BRD, jm, BRD)].y;
  dyp = center_V[IDX(BRD, jp, BRD)].y - part.y;
  dzm = part.z - center_V[IDX(BRD, BRD, km)].z;
  dzp = center_V[IDX(BRD, BRD, kp)].z - part.z;

  /* Compute volume coefficients */
  vol_ip_jp_kp = dxp*dyp*dzp;
  vol_im_jp_kp = dxm*dyp*dzp;
  vol_ip_jm_kp = dxp*dym*dzp;
  vol_ip_jp_km = dxp*dyp*dzm;
  vol_im_jm_kp = dxm*dym*dzp;
  vol_ip_jm_km = dxp*dym*dzm;
  vol_im_jp_km = dxm*dyp*dzm;
  vol_im_jm_km = dxm*dym*dzm;

  #ifdef LAGRANGE_TWOWAY_MOMENTUM

  /* compute rho_p /rho_f = (3-beta)/(2*beta) */
  density_ratio = (3.0 - (tracer+ipart)->beta_coeff)/(2.0*(tracer+ipart)->beta_coeff);

  /* this is the particle volume */
  fac =  (4./3.)*one_pi*pow( 3.0*property.nu*(tracer+ipart)->beta_coeff *(tracer+ipart)->tau_drag , 3./2. ) ;   

  /* this is the convention in our code 
   property.gravity_{x,y,z} just indicates the intensity 
   but the orientation along y is downward */
  gravity_vector.x =  property.gravity_x;
  gravity_vector.y = -property.gravity_y;
  gravity_vector.z =  property.gravity_z;

   /* build feedback */
  fp.x = (tracer+ipart)->Dt_ux - gravity_vector.x + density_ratio*(gravity_vector.x - (tracer+ipart)->ax);
  fp.y = (tracer+ipart)->Dt_uy - gravity_vector.y + density_ratio*(gravity_vector.y - (tracer+ipart)->ay);
  fp.z = (tracer+ipart)->Dt_uz - gravity_vector.z + density_ratio*(gravity_vector.z - (tracer+ipart)->az);

//fprintf(stderr, "fac %e\n",fac);
    /* test  */
  //fp.x =  0.0;
  //fp.y = density_ratio*(gravity_vector.y);
  //fp.z = 0.0;
 
  /* feedback in x */
  force_twoway[IDX(im, jm, km)].x +=  fac * fp.x * vol_ip_jp_kp;
  force_twoway[IDX(ip, jm, km)].x +=  fac * fp.x * vol_im_jp_kp;
  force_twoway[IDX(im, jp, km)].x +=  fac * fp.x * vol_ip_jm_kp;
  force_twoway[IDX(im, jm, kp)].x +=  fac * fp.x * vol_ip_jp_km; 
  force_twoway[IDX(ip, jp, km)].x +=  fac * fp.x * vol_im_jm_kp; 
  force_twoway[IDX(im, jp, kp)].x +=  fac * fp.x * vol_ip_jm_km;
  force_twoway[IDX(ip, jm, kp)].x +=  fac * fp.x * vol_im_jp_km; 
  force_twoway[IDX(ip, jp, kp)].x +=  fac * fp.x * vol_im_jm_km; 

  /* feedback in y */
  force_twoway[IDX(im, jm, km)].y +=  fac * fp.y * vol_ip_jp_kp;
  force_twoway[IDX(ip, jm, km)].y +=  fac * fp.y * vol_im_jp_kp; 
  force_twoway[IDX(im, jp, km)].y +=  fac * fp.y * vol_ip_jm_kp; 
  force_twoway[IDX(im, jm, kp)].y +=  fac * fp.y * vol_ip_jp_km;
  force_twoway[IDX(ip, jp, km)].y +=  fac * fp.y * vol_im_jm_kp;
  force_twoway[IDX(im, jp, kp)].y +=  fac * fp.y * vol_ip_jm_km;
  force_twoway[IDX(ip, jm, kp)].y +=  fac * fp.y * vol_im_jp_km;
  force_twoway[IDX(ip, jp, kp)].y +=  fac * fp.y * vol_im_jm_km;

  /* feedback in z */
  force_twoway[IDX(im, jm, km)].z +=  fac * fp.z * vol_ip_jp_kp; 
  force_twoway[IDX(ip, jm, km)].z +=  fac * fp.z * vol_im_jp_kp; 
  force_twoway[IDX(im, jp, km)].z +=  fac * fp.z * vol_ip_jm_kp;
  force_twoway[IDX(im, jm, kp)].z +=  fac * fp.z * vol_ip_jp_km;
  force_twoway[IDX(ip, jp, km)].z +=  fac * fp.z * vol_im_jm_kp;
  force_twoway[IDX(im, jp, kp)].z +=  fac * fp.z * vol_ip_jm_km;
  force_twoway[IDX(ip, jm, kp)].z +=  fac * fp.z * vol_im_jp_km;
  force_twoway[IDX(ip, jp, kp)].z +=  fac * fp.z * vol_im_jm_km;

  #endif /* end of LAGRANGE_TWO_WAY_MOMENTUM */
  #ifdef LAGRANGE_TWOWAY_TEMPERATURE

 /* this is the particle volume */
  fac =  (4./3.)*one_pi*pow( 3.0*property.nu*(tracer+ipart)->beta_coeff *(tracer+ipart)->tau_drag , 3./2. );    

  /* this term is equivalent to  3 * kappa / radius^2  but here expressed in terms of tau_drag and beta */
  fac2 = property.kappa / ( property.nu * (tracer+ipart)->beta_coeff * (tracer+ipart)->tau_drag ); 

  sp = ( (tracer+ipart)->t_p - (tracer+ipart)->t ) * fac2;

  /* feedback in x */
  t_source_twoway[IDX(im, jm, km)] +=  fac * sp * vol_ip_jp_kp;
  t_source_twoway[IDX(ip, jm, km)] +=  fac * sp * vol_im_jp_kp;
  t_source_twoway[IDX(im, jp, km)] +=  fac * sp * vol_ip_jm_kp;
  t_source_twoway[IDX(im, jm, kp)] +=  fac * sp * vol_ip_jp_km; 
  t_source_twoway[IDX(ip, jp, km)] +=  fac * sp * vol_im_jm_kp; 
  t_source_twoway[IDX(im, jp, kp)] +=  fac * sp * vol_ip_jm_km;
  t_source_twoway[IDX(ip, jm, kp)] +=  fac * sp * vol_im_jp_km; 
  t_source_twoway[IDX(ip, jp, kp)] +=  fac * sp * vol_im_jm_km; 

  #endif

 }/* end of loop on ipart */


/* now we communicate and add the infos on the borders */
  #ifdef LAGRANGE_TWOWAY_MOMENTUM
    add_sendrecv_borders_vector(force_twoway);
  #endif
  #ifdef LAGRANGE_TWOWAY_TEMPERATURE
    add_sendrecv_borders_scalar(t_source_twoway);
  #endif

/* now the field are added to the fluid force and source_t */
  for(i=0;i<(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD);i++){
  #ifdef LAGRANGE_TWOWAY_MOMENTUM  
    force[i].x +=  force_twoway[i].x;
    force[i].y +=  force_twoway[i].y;
    force[i].z +=  force_twoway[i].z;
  #endif
  #ifdef LAGRANGE_TWOWAY_TEMPERATURE
    t_source[i] +=  t_source_twoway[i];
  #endif
  }

}/* end of add_particle_feedbaks function */

#endif


/****************************************************************************************************/ 
void add_sendrecv_borders_scalar(my_double *f){
  int i,j,k,brd_size;
  MPI_Status status1;


  /*     BRD|LNX|BRD     */
  /* Copy borders along x */
  brd_size = BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_scalar[IDX_XBRD(i,j,k)] = f[IDX(i+LNX+BRD,j,k)];
      }

  MPI_Sendrecv( xp_scalar, brd_size, MPI_MY_DOUBLE, me_xp, 10,
                xm_scalar, brd_size, MPI_MY_DOUBLE, me_xm, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i+BRD,j,k)] += xm_scalar[IDX_XBRD(i,j,k)];
        xm_scalar[IDX_XBRD(i,j,k)] = f[IDX(i,j,k)];
      }

 MPI_Sendrecv( xm_scalar, brd_size, MPI_MY_DOUBLE, me_xm, 11,
               xp_scalar, brd_size, MPI_MY_DOUBLE, me_xp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX,j,k)] += xp_scalar[IDX_XBRD(i,j,k)];
      }



  /* Copy borders along y */
  brd_size = BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_scalar[IDX_YBRD(i,j,k)] = f[IDX(i,j+LNY+BRD,k)];
      }

  MPI_Sendrecv( yp_scalar, brd_size, MPI_MY_DOUBLE, me_yp, 10,
                ym_scalar, brd_size, MPI_MY_DOUBLE, me_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j+BRD,k)] += ym_scalar[IDX_YBRD(i,j,k)];
	    ym_scalar[IDX_YBRD(i,j,k)] = f[IDX(i,j,k)];
      }
  
 MPI_Sendrecv( ym_scalar, brd_size, MPI_MY_DOUBLE, me_ym, 13,
               yp_scalar, brd_size, MPI_MY_DOUBLE, me_yp, 13, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
       	f[IDX(i,j+LNY,k)] += yp_scalar[IDX_YBRD(i,j,k)];
      }
  

  /* Copy borders along z */
  brd_size = BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        zp_scalar[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+LNZ+BRD)];
      }

  MPI_Sendrecv( zp_scalar, brd_size, MPI_MY_DOUBLE, me_zp, 14,
                zm_scalar, brd_size, MPI_MY_DOUBLE, me_zm, 14, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k+BRD)] += zm_scalar[IDX_ZBRD(i,j,k)];
        zm_scalar[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k)];
      }
 MPI_Sendrecv( zm_scalar, brd_size, MPI_MY_DOUBLE, me_zm, 15,
               zp_scalar, brd_size, MPI_MY_DOUBLE, me_zp, 15, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j,k+LNZ)] += zp_scalar[IDX_ZBRD(i,j,k)];
      }


  /* First we communicate the 8 corner cubes (they are either 1x1x1 or 2x2x2 depending on BRD) */ 
  
  brd_size = BRD*BRD*BRD;

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_zp_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+LNZ)];
        xp_yp_zm_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+BRD)];
        xp_ym_zp_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+LNZ)];
        xm_yp_zp_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+LNZ)];
      }

  MPI_Sendrecv( xp_yp_zp_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xp_yp_zp, 10,
                xm_ym_zm_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xm_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_yp_zm_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xp_yp_zm, 10,
                xm_ym_zp_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xm_ym_zp, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_ym_zp_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xp_ym_zp, 10,
                xm_yp_zm_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xm_yp_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_zp_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xm_yp_zp, 10,
                xp_ym_zm_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xp_ym_zm, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] += xm_ym_zm_corner_scalar[IDX_CORNER(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] += xm_ym_zp_corner_scalar[IDX_CORNER(i,j,k)];
        f[IDX(i,j+LNY+BRD,k)] += xm_yp_zm_corner_scalar[IDX_CORNER(i,j,k)];
        f[IDX(i+LNX+BRD,j,k)] += xp_ym_zm_corner_scalar[IDX_CORNER(i,j,k)];
        xm_ym_zm_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+BRD)];
        xm_ym_zp_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+LNZ)];
        xm_yp_zm_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+BRD)];
        xp_ym_zm_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+BRD)];
      }
 MPI_Sendrecv( xm_ym_zm_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xm_ym_zm, 11,
               xp_yp_zp_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xp_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_ym_zp_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xm_ym_zp, 11,
               xp_yp_zm_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xp_yp_zm, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_yp_zm_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xm_yp_zm, 11,
               xp_ym_zp_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xp_ym_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_zm_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xp_ym_zm, 11,
               xm_yp_zp_corner_scalar, brd_size, MPI_MY_DOUBLE, me_xm_yp_zp, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k+LNZ+BRD)] += xp_yp_zp_corner_scalar[IDX_CORNER(i,j,k)];
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)] += xp_yp_zm_corner_scalar[IDX_CORNER(i,j,k)];
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)] += xp_ym_zp_corner_scalar[IDX_CORNER(i,j,k)];
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)] += xm_yp_zp_corner_scalar[IDX_CORNER(i,j,k)];
      }
 

 /* Then we communicate the 12 edges  */
 
 /* along x */
 
 brd_size = BRD*BRD*(LNX+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_zp_edge_scalar[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+LNZ)];
        yp_zm_edge_scalar[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+BRD)];
      }

  MPI_Sendrecv( yp_zp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_yp_zp, 10,
                ym_zm_edge_scalar, brd_size, MPI_MY_DOUBLE, me_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( yp_zm_edge_scalar, brd_size, MPI_MY_DOUBLE, me_yp_zm, 10,
                ym_zp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_ym_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] += ym_zm_edge_scalar[IDX_EDGE_X(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] += ym_zp_edge_scalar[IDX_EDGE_X(i,j,k)];
        ym_zm_edge_scalar[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+BRD)];
        ym_zp_edge_scalar[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+LNZ)];
      }
 MPI_Sendrecv( ym_zm_edge_scalar, brd_size, MPI_MY_DOUBLE, me_ym_zm, 11,
               yp_zp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( ym_zp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_ym_zp, 11,
               yp_zm_edge_scalar, brd_size, MPI_MY_DOUBLE, me_yp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)] += yp_zp_edge_scalar[IDX_EDGE_X(i,j,k)];
	f[IDX(i,j+LNY+BRD,k)] += yp_zm_edge_scalar[IDX_EDGE_X(i,j,k)];
      }
 
 /* along y */
 
 brd_size = BRD*BRD*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_zp_edge_scalar[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+LNZ)];
        xp_zm_edge_scalar[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+BRD)];
      }

  MPI_Sendrecv( xp_zp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xp_zp, 10,
                xm_zm_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xm_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_zm_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xp_zm, 10,
                xm_zp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xm_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] += xm_zm_edge_scalar[IDX_EDGE_Y(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] += xm_zp_edge_scalar[IDX_EDGE_Y(i,j,k)];
        xm_zm_edge_scalar[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+BRD)];
        xm_zp_edge_scalar[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+LNZ)];
      }
 MPI_Sendrecv( xm_zm_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xm_zm, 11,
               xp_zp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_zp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xm_zp, 11,
               xp_zm_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)] += xp_zp_edge_scalar[IDX_EDGE_Y(i,j,k)];
	f[IDX(i+LNX+BRD,j,k)] += xp_zm_edge_scalar[IDX_EDGE_Y(i,j,k)];
      }
 
 
 /* along z */
 
 brd_size = BRD*BRD*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_edge_scalar[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+LNY,k)];
        xm_yp_edge_scalar[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+LNY,k)];
      }

  MPI_Sendrecv( xp_yp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xp_yp, 10,
                xm_ym_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xm_ym, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xm_yp, 10,
                xp_ym_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xp_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] += xm_ym_edge_scalar[IDX_EDGE_Z(i,j,k)];
        f[IDX(i+LNX+BRD,j,k)] += xp_ym_edge_scalar[IDX_EDGE_Z(i,j,k)];
        xm_ym_edge_scalar[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+BRD,k)];
        xp_ym_edge_scalar[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+BRD,k)];
      }
 MPI_Sendrecv( xm_ym_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xm_ym, 11,
               xp_yp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xp_yp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xp_ym, 11,
               xm_yp_edge_scalar, brd_size, MPI_MY_DOUBLE, me_xm_yp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)] += xp_yp_edge_scalar[IDX_EDGE_Z(i,j,k)];
	f[IDX(i,j+LNY+BRD,k)] += xm_yp_edge_scalar[IDX_EDGE_Z(i,j,k)];
      }
 

}/* end scalar send rcv function */




/****************************************************************************************************/ 
void add_sendrecv_borders_vector(vector *f){
  int i,j,k,brd_size;
  MPI_Status status1;

  /*     BRD|LNX|BRD     */
  /* Copy borders along x */
  brd_size = BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_vector[IDX_XBRD(i,j,k)] = f[IDX(i+LNX,j,k)];
      }

  MPI_Sendrecv( xp_vector, brd_size, MPI_vector_type, me_xp, 10,
                xm_vector, brd_size, MPI_vector_type, me_xm, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)].x += xm_vector[IDX_XBRD(i,j,k)].x;
        f[IDX(i,j,k)].y += xm_vector[IDX_XBRD(i,j,k)].y;
        f[IDX(i,j,k)].z += xm_vector[IDX_XBRD(i,j,k)].z;
        xm_vector[IDX_XBRD(i,j,k)] = f[IDX(i+BRD,j,k)];
      }
 MPI_Sendrecv( xm_vector, brd_size, MPI_vector_type, me_xm, 11,
               xp_vector, brd_size, MPI_vector_type, me_xp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k)].x += xp_vector[IDX_XBRD(i,j,k)].x;
  f[IDX(i+LNX+BRD,j,k)].y += xp_vector[IDX_XBRD(i,j,k)].y;
  f[IDX(i+LNX+BRD,j,k)].z += xp_vector[IDX_XBRD(i,j,k)].z;
      }



  /* Copy borders along y */
  brd_size = BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_vector[IDX_YBRD(i,j,k)] = f[IDX(i,j+LNY,k)];
      }

  MPI_Sendrecv( yp_vector, brd_size, MPI_vector_type, me_yp, 10,
                ym_vector, brd_size, MPI_vector_type, me_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)].x += ym_vector[IDX_YBRD(i,j,k)].x;
        f[IDX(i,j,k)].y += ym_vector[IDX_YBRD(i,j,k)].y;
        f[IDX(i,j,k)].z += ym_vector[IDX_YBRD(i,j,k)].z;
	ym_vector[IDX_YBRD(i,j,k)] = f[IDX(i,j+BRD,k)];
      }
  
 MPI_Sendrecv( ym_vector, brd_size, MPI_vector_type, me_ym, 13,
               yp_vector, brd_size, MPI_vector_type, me_yp, 13, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
       	f[IDX(i,j+LNY+BRD,k)].x += yp_vector[IDX_YBRD(i,j,k)].x;
        f[IDX(i,j+LNY+BRD,k)].y += yp_vector[IDX_YBRD(i,j,k)].y;
        f[IDX(i,j+LNY+BRD,k)].z += yp_vector[IDX_YBRD(i,j,k)].z;
      }
  

  /* Copy borders along z */
  brd_size = BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        zp_vector[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+LNZ)];
      }

  MPI_Sendrecv( zp_vector, brd_size, MPI_vector_type, me_zp, 14,
                zm_vector, brd_size, MPI_vector_type, me_zm, 14, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)].x += zm_vector[IDX_ZBRD(i,j,k)].x;
        f[IDX(i,j,k)].y += zm_vector[IDX_ZBRD(i,j,k)].y;
        f[IDX(i,j,k)].z += zm_vector[IDX_ZBRD(i,j,k)].z;
        zm_vector[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+BRD)];
      }
 MPI_Sendrecv( zm_vector, brd_size, MPI_vector_type, me_zm, 15,
               zp_vector, brd_size, MPI_vector_type, me_zp, 15, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j,k+LNZ+BRD)].x += zp_vector[IDX_ZBRD(i,j,k)].x;
  f[IDX(i,j,k+LNZ+BRD)].y += zp_vector[IDX_ZBRD(i,j,k)].y;
  f[IDX(i,j,k+LNZ+BRD)].z += zp_vector[IDX_ZBRD(i,j,k)].z;
      }


#ifdef METHOD_EDGES_AND_CORNERS

  /* First we communicate the 8 corner cubes (they are either 1x1x1 or 2x2x2 depending on BRD) */ 
  
  brd_size = BRD*BRD*BRD;

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_zp_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+LNZ)];
        xp_yp_zm_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+BRD)];
        xp_ym_zp_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+LNZ)];
        xm_yp_zp_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+LNZ)];
      }

  MPI_Sendrecv( xp_yp_zp_corner_vector, brd_size, MPI_vector_type, me_xp_yp_zp, 10,
                xm_ym_zm_corner_vector, brd_size, MPI_vector_type, me_xm_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_yp_zm_corner_vector, brd_size, MPI_vector_type, me_xp_yp_zm, 10,
                xm_ym_zp_corner_vector, brd_size, MPI_vector_type, me_xm_ym_zp, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_ym_zp_corner_vector, brd_size, MPI_vector_type, me_xp_ym_zp, 10,
                xm_yp_zm_corner_vector, brd_size, MPI_vector_type, me_xm_yp_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_zp_corner_vector, brd_size, MPI_vector_type, me_xm_yp_zp, 10,
                xp_ym_zm_corner_vector, brd_size, MPI_vector_type, me_xp_ym_zm, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)].x += xm_ym_zm_corner_vector[IDX_CORNER(i,j,k)].x;
        f[IDX(i,j,k)].y += xm_ym_zm_corner_vector[IDX_CORNER(i,j,k)].y;
        f[IDX(i,j,k)].z += xm_ym_zm_corner_vector[IDX_CORNER(i,j,k)].z;
        f[IDX(i,j,k+LNZ+BRD)].x += xm_ym_zp_corner_vector[IDX_CORNER(i,j,k)].x;
        f[IDX(i,j,k+LNZ+BRD)].y += xm_ym_zp_corner_vector[IDX_CORNER(i,j,k)].y;
        f[IDX(i,j,k+LNZ+BRD)].z += xm_ym_zp_corner_vector[IDX_CORNER(i,j,k)].z;
        f[IDX(i,j+LNY+BRD,k)].x += xm_yp_zm_corner_vector[IDX_CORNER(i,j,k)].x;
        f[IDX(i,j+LNY+BRD,k)].y += xm_yp_zm_corner_vector[IDX_CORNER(i,j,k)].y;
        f[IDX(i,j+LNY+BRD,k)].z += xm_yp_zm_corner_vector[IDX_CORNER(i,j,k)].z;
        f[IDX(i+LNX+BRD,j,k)].x += xp_ym_zm_corner_vector[IDX_CORNER(i,j,k)].x;
        f[IDX(i+LNX+BRD,j,k)].y += xp_ym_zm_corner_vector[IDX_CORNER(i,j,k)].y;
        f[IDX(i+LNX+BRD,j,k)].z += xp_ym_zm_corner_vector[IDX_CORNER(i,j,k)].z;
        xm_ym_zm_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+BRD)];
        xm_ym_zp_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+LNZ)];
        xm_yp_zm_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+BRD)];
        xp_ym_zm_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+BRD)];
      }
 MPI_Sendrecv( xm_ym_zm_corner_vector, brd_size, MPI_vector_type, me_xm_ym_zm, 11,
               xp_yp_zp_corner_vector, brd_size, MPI_vector_type, me_xp_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_ym_zp_corner_vector, brd_size, MPI_vector_type, me_xm_ym_zp, 11,
               xp_yp_zm_corner_vector, brd_size, MPI_vector_type, me_xp_yp_zm, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_yp_zm_corner_vector, brd_size, MPI_vector_type, me_xm_yp_zm, 11,
               xp_ym_zp_corner_vector, brd_size, MPI_vector_type, me_xp_ym_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_zm_corner_vector, brd_size, MPI_vector_type, me_xp_ym_zm, 11,
               xm_yp_zp_corner_vector, brd_size, MPI_vector_type, me_xm_yp_zp, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k+LNZ+BRD)].x += xp_yp_zp_corner_vector[IDX_CORNER(i,j,k)].x;
  f[IDX(i+LNX+BRD,j+LNY+BRD,k+LNZ+BRD)].y += xp_yp_zp_corner_vector[IDX_CORNER(i,j,k)].y;
  f[IDX(i+LNX+BRD,j+LNY+BRD,k+LNZ+BRD)].z += xp_yp_zp_corner_vector[IDX_CORNER(i,j,k)].z;
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)].x += xp_yp_zm_corner_vector[IDX_CORNER(i,j,k)].x;
  f[IDX(i+LNX+BRD,j+LNY+BRD,k)].y += xp_yp_zm_corner_vector[IDX_CORNER(i,j,k)].y;
  f[IDX(i+LNX+BRD,j+LNY+BRD,k)].z += xp_yp_zm_corner_vector[IDX_CORNER(i,j,k)].z;
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)].x += xp_ym_zp_corner_vector[IDX_CORNER(i,j,k)].x;
  f[IDX(i+LNX+BRD,j,k+LNZ+BRD)].y += xp_ym_zp_corner_vector[IDX_CORNER(i,j,k)].y;
  f[IDX(i+LNX+BRD,j,k+LNZ+BRD)].z += xp_ym_zp_corner_vector[IDX_CORNER(i,j,k)].z;
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)].x += xm_yp_zp_corner_vector[IDX_CORNER(i,j,k)].x;
  f[IDX(i,j+LNY+BRD,k+LNZ+BRD)].y += xm_yp_zp_corner_vector[IDX_CORNER(i,j,k)].y;
  f[IDX(i,j+LNY+BRD,k+LNZ+BRD)].z += xm_yp_zp_corner_vector[IDX_CORNER(i,j,k)].z;
      }
 

 /* Then we communicate the 12 edges  */
 
 /* along x */
 
 brd_size = BRD*BRD*(LNX+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_zp_edge_vector[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+LNZ)];
        yp_zm_edge_vector[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+BRD)];
      }

  MPI_Sendrecv( yp_zp_edge_vector, brd_size, MPI_vector_type, me_yp_zp, 10,
                ym_zm_edge_vector, brd_size, MPI_vector_type, me_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( yp_zm_edge_vector, brd_size, MPI_vector_type, me_yp_zm, 10,
                ym_zp_edge_vector, brd_size, MPI_vector_type, me_ym_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)].x += ym_zm_edge_vector[IDX_EDGE_X(i,j,k)].x;
        f[IDX(i,j,k)].y += ym_zm_edge_vector[IDX_EDGE_X(i,j,k)].y;
        f[IDX(i,j,k)].z += ym_zm_edge_vector[IDX_EDGE_X(i,j,k)].z;
        f[IDX(i,j,k+LNZ+BRD)].x += ym_zp_edge_vector[IDX_EDGE_X(i,j,k)].x;
        f[IDX(i,j,k+LNZ+BRD)].y += ym_zp_edge_vector[IDX_EDGE_X(i,j,k)].y;
        f[IDX(i,j,k+LNZ+BRD)].z += ym_zp_edge_vector[IDX_EDGE_X(i,j,k)].z;
        ym_zm_edge_vector[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+BRD)];
        ym_zp_edge_vector[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+LNZ)];
      }
 MPI_Sendrecv( ym_zm_edge_vector, brd_size, MPI_vector_type, me_ym_zm, 11,
               yp_zp_edge_vector, brd_size, MPI_vector_type, me_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( ym_zp_edge_vector, brd_size, MPI_vector_type, me_ym_zp, 11,
               yp_zm_edge_vector, brd_size, MPI_vector_type, me_yp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)].x += yp_zp_edge_vector[IDX_EDGE_X(i,j,k)].x;
  f[IDX(i,j+LNY+BRD,k+LNZ+BRD)].y += yp_zp_edge_vector[IDX_EDGE_X(i,j,k)].y;
  f[IDX(i,j+LNY+BRD,k+LNZ+BRD)].z += yp_zp_edge_vector[IDX_EDGE_X(i,j,k)].z;
	f[IDX(i,j+LNY+BRD,k)].x += yp_zm_edge_vector[IDX_EDGE_X(i,j,k)].x;
  f[IDX(i,j+LNY+BRD,k)].y += yp_zm_edge_vector[IDX_EDGE_X(i,j,k)].y;
  f[IDX(i,j+LNY+BRD,k)].z += yp_zm_edge_vector[IDX_EDGE_X(i,j,k)].z;
      }
 
 /* along y */
 
 brd_size = BRD*BRD*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_zp_edge_vector[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+LNZ)];
        xp_zm_edge_vector[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+BRD)];
      }

  MPI_Sendrecv( xp_zp_edge_vector, brd_size, MPI_vector_type, me_xp_zp, 10,
                xm_zm_edge_vector, brd_size, MPI_vector_type, me_xm_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_zm_edge_vector, brd_size, MPI_vector_type, me_xp_zm, 10,
                xm_zp_edge_vector, brd_size, MPI_vector_type, me_xm_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)].x += xm_zm_edge_vector[IDX_EDGE_Y(i,j,k)].x;
        f[IDX(i,j,k)].y += xm_zm_edge_vector[IDX_EDGE_Y(i,j,k)].y;
        f[IDX(i,j,k)].z += xm_zm_edge_vector[IDX_EDGE_Y(i,j,k)].z;
        f[IDX(i,j,k+LNZ+BRD)].x += xm_zp_edge_vector[IDX_EDGE_Y(i,j,k)].x;
        f[IDX(i,j,k+LNZ+BRD)].y += xm_zp_edge_vector[IDX_EDGE_Y(i,j,k)].y;
        f[IDX(i,j,k+LNZ+BRD)].z += xm_zp_edge_vector[IDX_EDGE_Y(i,j,k)].z;
        xm_zm_edge_vector[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+BRD)];
        xm_zp_edge_vector[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+LNZ)];
      }
 MPI_Sendrecv( xm_zm_edge_vector, brd_size, MPI_vector_type, me_xm_zm, 11,
               xp_zp_edge_vector, brd_size, MPI_vector_type, me_xp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_zp_edge_vector, brd_size, MPI_vector_type, me_xm_zp, 11,
               xp_zm_edge_vector, brd_size, MPI_vector_type, me_xp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)].x += xp_zp_edge_vector[IDX_EDGE_Y(i,j,k)].x;
  f[IDX(i+LNX+BRD,j,k+LNZ+BRD)].y += xp_zp_edge_vector[IDX_EDGE_Y(i,j,k)].y;
  f[IDX(i+LNX+BRD,j,k+LNZ+BRD)].z += xp_zp_edge_vector[IDX_EDGE_Y(i,j,k)].z;
	f[IDX(i+LNX+BRD,j,k)].x += xp_zm_edge_vector[IDX_EDGE_Y(i,j,k)].x;
  f[IDX(i+LNX+BRD,j,k)].y += xp_zm_edge_vector[IDX_EDGE_Y(i,j,k)].y;
  f[IDX(i+LNX+BRD,j,k)].z += xp_zm_edge_vector[IDX_EDGE_Y(i,j,k)].z;
      }
 
 
 /* along z */
 
 brd_size = BRD*BRD*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_edge_vector[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+LNY,k)];
        xm_yp_edge_vector[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+LNY,k)];
      }

  MPI_Sendrecv( xp_yp_edge_vector, brd_size, MPI_vector_type, me_xp_yp, 10,
                xm_ym_edge_vector, brd_size, MPI_vector_type, me_xm_ym, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_edge_vector, brd_size, MPI_vector_type, me_xm_yp, 10,
                xp_ym_edge_vector, brd_size, MPI_vector_type, me_xp_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)].x += xm_ym_edge_vector[IDX_EDGE_Z(i,j,k)].x;
        f[IDX(i,j,k)].y += xm_ym_edge_vector[IDX_EDGE_Z(i,j,k)].y;
        f[IDX(i,j,k)].z += xm_ym_edge_vector[IDX_EDGE_Z(i,j,k)].z;
        f[IDX(i+LNX+BRD,j,k)].x += xp_ym_edge_vector[IDX_EDGE_Z(i,j,k)].x;
        f[IDX(i+LNX+BRD,j,k)].y += xp_ym_edge_vector[IDX_EDGE_Z(i,j,k)].y;
        f[IDX(i+LNX+BRD,j,k)].z += xp_ym_edge_vector[IDX_EDGE_Z(i,j,k)].z;
        xm_ym_edge_vector[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+BRD,k)];
        xp_ym_edge_vector[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+BRD,k)];
      }
 MPI_Sendrecv( xm_ym_edge_vector, brd_size, MPI_vector_type, me_xm_ym, 11,
               xp_yp_edge_vector, brd_size, MPI_vector_type, me_xp_yp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_edge_vector, brd_size, MPI_vector_type, me_xp_ym, 11,
               xm_yp_edge_vector, brd_size, MPI_vector_type, me_xm_yp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)].x += xp_yp_edge_vector[IDX_EDGE_Z(i,j,k)].x;
  f[IDX(i+LNX+BRD,j+LNY+BRD,k)].y += xp_yp_edge_vector[IDX_EDGE_Z(i,j,k)].y;
  f[IDX(i+LNX+BRD,j+LNY+BRD,k)].z += xp_yp_edge_vector[IDX_EDGE_Z(i,j,k)].z;
	f[IDX(i,j+LNY+BRD,k)].x += xm_yp_edge_vector[IDX_EDGE_Z(i,j,k)].x;
  f[IDX(i,j+LNY+BRD,k)].y += xm_yp_edge_vector[IDX_EDGE_Z(i,j,k)].y;
  f[IDX(i,j+LNY+BRD,k)].z += xm_yp_edge_vector[IDX_EDGE_Z(i,j,k)].z;
      }
 
#endif /* METHOD_EDGES_AND_CORNERS */

}/* end vector send rcv function */