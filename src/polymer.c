#include "common_object.h"

#ifdef LAGRANGE_POLYMER
 #ifdef LAGRANGE_POLYMER_FEEDBACK
/* Extrapolation, for the moment only with regular grid of unit spacing */
void add_lagrangian_polymer_feedback_on_the_flow(){
 point_particle part;
 my_double fac;
 my_double cxx,cyy,czz,cxy,cyz,cxz;

 int ipart,im,jm,km,ip,jp,kp;
 double dxm,dxp,dym,dyp,dzm,dzp;
 double vol_ip_jp_kp,vol_im_jp_kp , vol_ip_jm_kp , vol_ip_jp_km , vol_im_jm_kp , vol_ip_jm_km , vol_im_jp_km , vol_im_jm_km; 
 int i,j,k;
 double f_im_jm_km , f_ip_jm_km , f_im_jp_km , f_im_jm_kp , f_ip_jp_km , f_im_jp_kp , f_ip_jm_kp , f_ip_jp_kp; 

 double dx_jm_km , dx_jp_kp , dx_jp_km , dx_jm_kp; 
 double dy_im_km , dy_ip_kp , dy_ip_km , dy_im_kp; 
 double dz_im_jm , dz_ip_jp , dz_ip_jm , dz_im_jp; 

 for (ipart=0;ipart<npart;ipart++) {

   /* get tensor value */
  cxx = (tracer+ipart)->cxx;
  cyy = (tracer+ipart)->cyy;
  czz = (tracer+ipart)->czz;
  cxy = (tracer+ipart)->cxy;
  cyz = (tracer+ipart)->cyz;
  cxz = (tracer+ipart)->cxz;

  /* C_ij - \delta_ij  */
  cxx = cxx -1.0;
  cyy = cyy -1.0;
  czz = czz -1.0;

  fac = property.nu_polymer / property.tau_polymer;

   /* get coordinates in the domain */
  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

  /* get indexes of neighboring nodes */

for (i=0; i<LNX+TWO_BRD-1; i++) if(center_V[IDX(i, BRD, BRD)].x <= part.x && part.x < center_V[IDX(i+1,BRD, BRD)].x) im = i; 
ip =  im + 1;
for (j=0; j<LNY+TWO_BRD-1; j++) if(center_V[IDX(BRD, j, BRD)].y <= part.y && part.y < center_V[IDX(BRD, j+1, BRD)].y) jm = j;
jp =  jm + 1;
for (k=0; k<LNZ+TWO_BRD-1; k++) if(center_V[IDX(BRD, BRD, k)].z <= part.z && part.z < center_V[IDX(BRD, BRD, k+1)].z) km = k;
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

  /* extrapolate value to field (just change name) */

  f_im_jm_km = vol_ip_jp_kp; 
  f_ip_jm_km = vol_im_jp_kp;
  f_im_jp_km = vol_ip_jm_kp;
  f_im_jm_kp = vol_ip_jp_km;
  f_ip_jp_km = vol_im_jm_kp;
  f_im_jp_kp = vol_ip_jm_km;
  f_ip_jm_kp = vol_im_jp_km;
  f_ip_jp_kp = vol_im_jm_km;


  /* field x gradients  */
  dx_jm_km = f_ip_jm_km - f_im_jm_km;
  dx_jp_kp = f_ip_jp_kp - f_im_jp_kp;
  dx_jp_km = f_ip_jp_km - f_im_jp_km;
  dx_jm_kp = f_ip_jm_kp - f_im_jm_kp;

  /* field y gradients  */
  dy_im_km = f_im_jp_km - f_im_jm_km;
  dy_im_km = f_im_jp_km - f_im_jm_km;
  dy_im_km = f_im_jp_km - f_im_jm_km;
  dy_im_km = f_im_jp_km - f_im_jm_km;

  /* field z gradients  */
  dz_im_jm = f_im_jm_kp - f_im_jm_km;
  dz_im_jm = f_im_jm_kp - f_im_jm_km;
  dz_im_jm = f_im_jm_kp - f_im_jm_km;
  dz_im_jm = f_im_jm_kp - f_im_jm_km;

  /* we want to compute : 
     (1) dx cxx + dy cyx + dz czx
     (2) dx cxy + dy cyy + dz czy
     (3) dx cxz + dy cyz + dz czz
     due to th symmetry :
     (1) dx cxx + dy cxy + dz cxz
     (2) dx cxy + dy cyy + dz cyz
     (3) dx cxz + dy cyz + dz czz
  */

  // fprintf(stderr,"%d %d %d %d %d %d\n",im,ip,jm,jp,km,kp);
 
  /* divergence x */
  force[IDX(im, jm, km)].x +=  fac*( cxx * dx_jm_km + cxy * dy_im_km + cxz * dz_im_jm); 
  force[IDX(ip, jm, km)].x +=  fac*( cxx * dx_jm_km + cxy * dy_ip_km + cxz * dz_ip_jm); 
  force[IDX(im, jp, km)].x +=  fac*( cxx * dx_jp_km + cxy * dy_im_km + cxz * dz_im_jp); 
  force[IDX(im, jm, kp)].x +=  fac*( cxx * dx_jm_kp + cxy * dy_im_kp + cxz * dz_im_jm); 
  force[IDX(ip, jp, km)].x +=  fac*( cxx * dx_jp_km + cxy * dy_ip_km + cxz * dz_ip_jp); 
  force[IDX(im, jp, kp)].x +=  fac*( cxx * dx_jp_kp + cxy * dy_im_kp + cxz * dz_im_jp); 
  force[IDX(ip, jm, kp)].x +=  fac*( cxx * dx_jm_kp + cxy * dy_ip_kp + cxz * dz_ip_jm); 
  force[IDX(ip, jp, kp)].x +=  fac*( cxx * dx_jp_kp + cxy * dy_ip_kp + cxz * dz_ip_jp); 

  /* divergence y */
  force[IDX(im, jm, km)].y +=  fac*( cxy * dx_jm_km + cyy * dy_im_km + cyz * dz_im_jm); 
  force[IDX(ip, jm, km)].y +=  fac*( cxy * dx_jm_km + cyy * dy_ip_km + cyz * dz_ip_jm); 
  force[IDX(im, jp, km)].y +=  fac*( cxy * dx_jp_km + cyy * dy_im_km + cyz * dz_im_jp); 
  force[IDX(im, jm, kp)].y +=  fac*( cxy * dx_jm_kp + cyy * dy_im_kp + cyz * dz_im_jm); 
  force[IDX(ip, jp, km)].y +=  fac*( cxy * dx_jp_km + cyy * dy_ip_km + cyz * dz_ip_jp); 
  force[IDX(im, jp, kp)].y +=  fac*( cxy * dx_jp_kp + cyy * dy_im_kp + cyz * dz_im_jp); 
  force[IDX(ip, jm, kp)].y +=  fac*( cxy * dx_jm_kp + cyy * dy_ip_kp + cyz * dz_ip_jm); 
  force[IDX(ip, jp, kp)].y +=  fac*( cxy * dx_jp_kp + cyy * dy_ip_kp + cyz * dz_ip_jp); 

  /* divergence z */
  force[IDX(im, jm, km)].z +=  fac*( cxz * dx_jm_km + cyz * dy_im_km + czz * dz_im_jm); 
  force[IDX(ip, jm, km)].z +=  fac*( cxz * dx_jm_km + cyz * dy_ip_km + czz * dz_ip_jm); 
  force[IDX(im, jp, km)].z +=  fac*( cxz * dx_jp_km + cyz * dy_im_km + czz * dz_im_jp); 
  force[IDX(im, jm, kp)].z +=  fac*( cxz * dx_jm_kp + cyz * dy_im_kp + czz * dz_im_jm); 
  force[IDX(ip, jp, km)].z +=  fac*( cxz * dx_jp_km + cyz * dy_ip_km + czz * dz_ip_jp); 
  force[IDX(im, jp, kp)].z +=  fac*( cxz * dx_jp_kp + cyz * dy_im_kp + czz * dz_im_jp); 
  force[IDX(ip, jm, kp)].z +=  fac*( cxz * dx_jm_kp + cyz * dy_ip_kp + czz * dz_ip_jm); 
  force[IDX(ip, jp, kp)].z +=  fac*( cxz * dx_jp_kp + cyz * dy_ip_kp + czz * dz_ip_jp); 

 }/* end of loop on ipart */

}/* end of add_... function */
 #endif /* LAGRANGE_POLYMER_FEEDBACK */

void evolve_lagrangian_polymer_conformation_tensor(int ipart){

  int i,j,k;
  my_double invtau;
  my_double matI[3][3];
  my_double matC[3][3];
  my_double matF[3][3];
  my_double matFold[3][3];
  my_double matA[3][3];

  int dim = 3;  /* grid dimension */
#ifdef GRID_POP_D2Q9
  dim = 2;
#endif
              /* polymer relaxation time */
              invtau = 1./property.tau_polymer;

              /* Define the identity matrix */
              matI[0][0] = matI[1][1] = matI[2][2] = 1.0;
              matI[0][1] = matI[0][2] = matI[1][2] = 0.0;
              matI[1][0] = matI[2][0] = matI[2][1] = 0.0;

              /* velocity gradient matrix */
              matA[0][0]=(tracer+ipart)->dx_ux ; matA[0][1]=(tracer+ipart)->dy_ux; matA[0][2]=(tracer+ipart)->dz_ux;
              matA[1][0]=(tracer+ipart)->dx_uy ; matA[1][1]=(tracer+ipart)->dy_uy; matA[1][2]=(tracer+ipart)->dz_uy;
              matA[2][0]=(tracer+ipart)->dx_uz ; matA[2][1]=(tracer+ipart)->dy_uz; matA[2][2]=(tracer+ipart)->dz_uz;   

#ifdef DEBUG
	      if( (tracer+ipart)->name ==0)fprintf(stdout,"before %e %e %e %e %e %e\n", (tracer+ipart)->cxx, (tracer+ipart)->cyy, (tracer+ipart)->czz ,
						                                 (tracer+ipart)->cxy, (tracer+ipart)->cyz, (tracer+ipart)->cxz );
#endif

              /* conformation */
              matC[0][0]=(tracer+ipart)->cxx ; matC[0][1]=(tracer+ipart)->cxy; matC[0][2]=(tracer+ipart)->cxz;
              matC[1][0]=(tracer+ipart)->cxy ; matC[1][1]=(tracer+ipart)->cyy; matC[1][2]=(tracer+ipart)->cyz;
              matC[2][0]=(tracer+ipart)->cxz ; matC[2][1]=(tracer+ipart)->cyz; matC[2][2]=(tracer+ipart)->czz;   

              /* d_t conformation */
              matFold[0][0]=(tracer+ipart)->dt_cxx ; matFold[0][1]=(tracer+ipart)->dt_cxy; matFold[0][2]=(tracer+ipart)->dt_cxz;
              matFold[1][0]=(tracer+ipart)->dt_cxy ; matFold[1][1]=(tracer+ipart)->dt_cyy; matFold[1][2]=(tracer+ipart)->dt_cyz;
              matFold[2][0]=(tracer+ipart)->dt_cxz ; matFold[2][1]=(tracer+ipart)->dt_cyz; matFold[2][2]=(tracer+ipart)->dt_czz;   

	      /* Can be optimized by taking into account symmetry , but be careful! */
	      /* compute ( grad_u * c + c * grad_u^T) + (1-c)/\tau_polym */
              for (i=0; i<dim; i++)
                for (j=0; j<dim; j++){
		  matF[i][j] = 0.0;
		  for (k=0; k<dim; k++){
		  matF[i][j] +=  matA[k][i]*matC[k][j] + matC[i][k]*matA[k][j]; 
		  }
		  matF[i][j] += (matI[i][j] - matC[i][j])*invtau;
		}


	      /* NOTE THAT the loop here below are only on i<=j indexes due to the symmetry */
	      /* if restart Euler 1st order C = C0 + (DT)*F  */
	      if(itime==0 && resume==0){
              for (j=0; j<dim; j++)
		for (i=0; i<=j; i++){              
                    matC[i][j] =  matC[i][j] + property.time_dt*matF[i][j];        
                }
	      }else{
		/* AB 2nd order C = C0 + (DT/2)*(3*F - Fold)  */
              for (j=0; j<dim; j++)
		for (i=0; i<=j; i++){
                  matC[i][j] =  matC[i][j] + 0.5*property.time_dt*( 3.*matF[i][j] -  matFold[i][j] );
                }
	      }

              /* assign C tensor */
	      (tracer+ipart)->cxx = matC[0][0];
	      (tracer+ipart)->cyy = matC[1][1];
	      (tracer+ipart)->czz = matC[2][2]; 
	      (tracer+ipart)->cxy = matC[0][1];
	      (tracer+ipart)->cyz = matC[1][2];
	      (tracer+ipart)->cxz = matC[0][2]; 

              /* assign the just computed dC/dt symmetric tensor */ 
	      (tracer+ipart)->dt_cxx = matF[0][0];
	      (tracer+ipart)->dt_cyy = matF[1][1];
	      (tracer+ipart)->dt_czz = matF[2][2]; 
	      (tracer+ipart)->dt_cxy = matF[0][1];
	      (tracer+ipart)->dt_cyz = matF[1][2];
	      (tracer+ipart)->dt_cxz = matF[0][2]; 

#ifdef DEBUG
	      if( (tracer+ipart)->name ==0)fprintf(stdout,"after %e %e %e %e %e %e\n", (tracer+ipart)->cxx, (tracer+ipart)->cyy, (tracer+ipart)->czz ,
						                                 (tracer+ipart)->cxy, (tracer+ipart)->cyz, (tracer+ipart)->cxz );
#endif

}/* end of evolve_lagrangian_conformation_tensor */

#endif /* LAGRANGE_POLYMER */
