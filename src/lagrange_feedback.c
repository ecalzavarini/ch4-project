
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

 for (ipart=0;ipart<npart;ipart++) {

  fac = 1.0; //(property.SX*property.SY*property.SZ)/property.particle_number;   

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
  density_ratio = (3.0 - (tracer+ipart)->beta)/(2.0*(tracer+ipart)->beta_coeff);

  /* this is the convention in our code 
   property.gravity_{x,y,z} just indicates the intensity 
   but the orientation along y is downward */
  gravity_vector.x =  property.gravity_x
  gravity_vector.y = -property.gravity_y
  gravity_vector.z =  property.gravity_z

   /* build feedback */
  fp.x = (tracer+ipart)->Dt_ux - gravity_vector.x + density_ratio*(gravity_vector.x - (tracer+ipart)->ax);
  fp.y = (tracer+ipart)->Dt_uy - gravity_vector.y + density_ratio*(gravity_vector.y - (tracer+ipart)->ay);
  fp.z = (tracer+ipart)->Dt_uz - gravity_vector.z + density_ratio*(gravity_vector.z - (tracer+ipart)->az);
 
  /* feedback in x */
  force[IDX(im, jm, km)].x +=  fac * fp.x * vol_ip_jp_kp;
  force[IDX(ip, jm, km)].x +=  fac * fp.x * vol_im_jp_kp;
  force[IDX(im, jp, km)].x +=  fac * fp.x * vol_ip_jm_kp;
  force[IDX(im, jm, kp)].x +=  fac * fp.x * vol_ip_jp_km; 
  force[IDX(ip, jp, km)].x +=  fac * fp.x * vol_im_jm_kp; 
  force[IDX(im, jp, kp)].x +=  fac * fp.x * vol_ip_jm_km;
  force[IDX(ip, jm, kp)].x +=  fac * fp.x * vol_im_jp_km; 
  force[IDX(ip, jp, kp)].x +=  fac * fp.x * vol_im_jm_km; 

  /* feedback in y */
  force[IDX(im, jm, km)].y +=  fac * fp.y * vol_ip_jp_kp;
  force[IDX(ip, jm, km)].y +=  fac * fp.y * vol_im_jp_kp; 
  force[IDX(im, jp, km)].y +=  fac * fp.y * vol_ip_jm_kp; 
  force[IDX(im, jm, kp)].y +=  fac * fp.y * vol_ip_jp_km;
  force[IDX(ip, jp, km)].y +=  fac * fp.y * vol_im_jm_kp;
  force[IDX(im, jp, kp)].y +=  fac * fp.y * vol_ip_jm_km;
  force[IDX(ip, jm, kp)].y +=  fac * fp.y * vol_im_jp_km;
  force[IDX(ip, jp, kp)].y +=  fac * fp.y * vol_im_jm_km;

  /* feedback in z */
  force[IDX(im, jm, km)].z +=  fac * fp.z * vol_ip_jp_kp; 
  force[IDX(ip, jm, km)].z +=  fac * fp.z * vol_im_jp_kp; 
  force[IDX(im, jp, km)].z +=  fac * fp.z * vol_ip_jm_kp;
  force[IDX(im, jm, kp)].z +=  fac * fp.z * vol_ip_jp_km;
  force[IDX(ip, jp, km)].z +=  fac * fp.z * vol_im_jm_kp;
  force[IDX(im, jp, kp)].z +=  fac * fp.z * vol_ip_jm_km;
  force[IDX(ip, jm, kp)].z +=  fac * fp.z * vol_im_jp_km;
  force[IDX(ip, jp, kp)].z +=  fac * fp.z * vol_im_jm_km;

  #endif /* end of LAGRANGE_TWO_WAY_MOMENTUM */
  #ifdef LAGRANGE_TWOWAY_TEMPERATURE

  /* this term is equivalent to  3 * kappa / radius^2  but here expressed in terms of tau_drag and beta */
  fac2 = property.kappa / ( property.nu * (tracer+ipart)->beta_coeff * (tracer+ipart)->tau_drag ); 

  sp = ( (tracer+ipart)->t_p - (tracer+ipart)->t ) * fac2;

  /* feedback in x */
  source_t[IDX(im, jm, km)] +=  fac * sp * vol_ip_jp_kp;
  source_t[IDX(ip, jm, km)] +=  fac * sp * vol_im_jp_kp;
  source_t[IDX(im, jp, km)] +=  fac * sp * vol_ip_jm_kp;
  source_t[IDX(im, jm, kp)] +=  fac * sp * vol_ip_jp_km; 
  source_t[IDX(ip, jp, km)] +=  fac * sp * vol_im_jm_kp; 
  source_t[IDX(im, jp, kp)] +=  fac * sp * vol_ip_jm_km;
  source_t[IDX(ip, jm, kp)] +=  fac * sp * vol_im_jp_km; 
  source_t[IDX(ip, jp, kp)] +=  fac * sp * vol_im_jm_km; 

  #endif

 }/* end of loop on ipart */

}/* end of add_particle_feedbaks function */

#endif