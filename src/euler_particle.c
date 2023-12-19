#include "common_object.h"


/* function to copy a scalar field */
void copy_scalar(my_double *ff, my_double *ff_copy){
  int i,j,k;
  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){

	ff_copy[IDX(i,j,k)] = ff[IDX(i,j,k)]; 	

	}
}


/* function to copy a vector field */
void copy_vector(vector *ff, vector *ff_copy){
  int i,j,k;
  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){

	ff_copy[IDX(i,j,k)] = ff[IDX(i,j,k)]; 	

	}
}

#ifdef EULER_PARTICLE
    #ifdef EULER_PARTICLE_CONCENTRATION
/* Minmod limiter function */
my_double minmod(my_double a, my_double b) {
     //my_double sgn_a = (a > 0) - (a < 0);
     //my_double sgn_b = (b > 0) - (b < 0);
     my_double sgn_a = (a>0)? 1.0 : -1.0;
     my_double sgn_b = (b>0)? 1.0 : -1.0;
     return 0.5 * (sgn_a + sgn_b) * fmin(fabs(a), fabs(b));
}
#ifdef AAA
/* Minmod Van Leer limiter function */
my_double minmod_theta(my_double a, my_double b, my_double c){
     if (a > 0 && b > 0 && c > 0) {
          return fmin(fmin(a, b), c);
     } else if (a < 0 && b < 0 && c < 0) {
          return fmax(fmax(a, b), c);
          } else {
          return 0.0;
          }
}
#endif
    #endif
#endif

/* IMPORTANT: We always assume here that the grid and time spacing are = 1 */
#ifdef EULER_PARTICLE
    #ifdef EULER_PARTICLE_CONCENTRATION
/* which method? */
#define TEST
//#define UPWIND
//#define QUICK
#define KT
//#define WENO

/* compute the right hand side of the concentration equation */
/*   \partial_t conc = - \nabla \cdot ( u * conc)   */
void compute_rhs_conc(){
 int i,j,k;
 my_double dx, dy, dz, dt;
 dx = dy = dz = dt = 1;
 int ip,jp,kp,im,km,jm,ipp,jpp,kpp,imm,jmm,kmm, ippp,jppp,kppp,immm,jmmm,kmmm;
 /* upwind method coefficients */
 #ifdef UPWIND
 my_double uface;
 my_double uface_xp, uface_yp, uface_zp;
 my_double uface_xm, uface_ym, uface_zm;
 #endif
 /* quick method coefficients */
 #ifdef QUICK
 my_double uface;
 my_double coeff_D, coeff_U, coeff_UU;
 coeff_D = 3.0/8.0;
 coeff_U = 6.0/8.0;
 coeff_UU = -1.0/8.0;
 #endif
 /* KT method coefficients */
 #ifdef KT
 my_double uL, uR, vL, vR, wL, wR, uL1, uR1, vL1, vR1, wL1, wR1, cL, cR, dL, dR, eL, eR, cL1, cR1, dL1, dR1, eL1, eR1;
 my_double ux1, ux2, uy1, uy2, uz1, uz2, ux11, uy11, uz11, ux22, uy22, uz22, cx1, cx2, cy1, cy2, cz1, cz2, cx11, cy11, cz11, cx22, cy22, cz22;
 my_double a, b, c, d, e, f, H1, H2, H3, H4, H5, H6;
 my_double fL1, fR1, fL2, fR2, fL3, fR3, fL4, fR4, fL5, fR5, fL6, fR6;
 #endif
 /* WENO method coefficients  */
 #ifdef WENO
   vector a;
   vector uL,uLm,uLmm,uLp,uLpp;
   vector uR,uRm,uRmm,uRp,uRpp;
   vector p0n,p1n,p2n, p0p,p1p,p2p;
   vector B0n,B1n,B2n, B0p,B1p,B2p;
   vector alpha0n,alpha1n,alpha2n,alphasumn,alpha0p,alpha1p,alpha2p,alphasump;
   vector w0n,w1n,w2n, w0p,w1p,w2p;
   vector hn, hp,hnm,hpm;
   my_double d0n,d1n,d2n,d0p,d1p,d2p,epsilon;
 #endif

/* we store the previous rhs */
copy_scalar(old_rhs_conc, old_old_rhs_conc);
copy_scalar(rhs_conc, old_rhs_conc);

/* apply the boundary condition for the fluid hydrodynamics fields */
boundary_conditions_hydro(); 

#ifdef TEST
	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {
          u[IDX(i, j, k)].x = 0.001;
          u[IDX(i, j, k)].y = 0.0;
          u[IDX(i, j, k)].z = 0.0;
               }
     boundary_conditions_hydro();
#endif


	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {
               
               /* precompute the indexes */
               ip  = i+1; jp  = j+1; kp  = k+1;
               im  = i-1; jm  = j-1; km  = k-1;
               ipp = i+2; jpp = j+2; kpp = k+2;
               imm = i-2; jmm = j-2; kmm = k-2;

               rhs_conc[IDX(i,j,k)] = 0.0;

#ifdef UPWIND 
            /*  we implement the linear upwind scheme for the Finite Volume method*/     
            /* NOTE : this scheme is TOO diffusive! */

            /* x + 1/2 face */
            uface =  0.5*( u[IDX(i,j,k)].x + u[IDX(ip,j,k)].x);
            if( uface  >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(ip,j,k)];

            /* x - 1/2 face */
            uface = 0.5 * (u[IDX(i,j,k)].x + u[IDX(im,j,k)].x);
	        if (uface >  0)
		        rhs_conc[IDX(i,j,k)] += uface * conc[IDX(im,j,k)];
            else
		        rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];
            
            /* y + 1/2 face */
	       uface =  0.5*( u[IDX(i,j,k)].y + u[IDX(i,jp,k)].y);
            if( uface >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,jp,k)];
            
            /* y - 1/2 face */
	        uface = 0.5 * (u[IDX(i,j,k)].y + u[IDX(i,jm,k)].y);
	        if( uface >0)
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,jm,k)];
            else
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];	

            /* z + 1/2 face */
	        uface = 0.5 * (u[IDX(i,j,k)].z + u[IDX(i,j,kp)].z);
            if( uface >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,kp)];
            
            /* z - 1/2 face */
    	        uface = 0.5 * (u[IDX(i,j,k)].z + u[IDX(i,j,km)].z);
	        if (uface > 0)
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,km)];
	        else
                rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];		    
#endif /* UPWIND */

#ifdef QUICK
        /*  we implement the QUICK scheme for the Finite Volume method*/

            /* x + 1/2 face */
            uface =  0.5*( u[IDX(i,j,k)].x + u[IDX(ip,j,k)].x);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(im, j, k)] + coeff_D * conc[IDX(ip, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(ip, j, k)] + coeff_UU * conc[IDX(ipp, j, k)] + coeff_D * conc[IDX(i, j, k)]);

            /* x - 1/2 face */
            uface = 0.5 * (u[IDX(i,j,k)].x + u[IDX(im,j,k)].x);
	        if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(im, j, k)] + coeff_UU * conc[IDX(imm, j, k)] + coeff_D * conc[IDX(i, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(ip, j, k)] + coeff_D * conc[IDX(im, j, k)]);
            
            /* y + 1/2 face */
            uface = 0.5 * (u[IDX(i, j, k)].y + u[IDX(1, jp, k)].y);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, jm, k)] + coeff_D * conc[IDX(i, jp, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, jp, k)] + coeff_UU * conc[IDX(i, jpp, k)] + coeff_D * conc[IDX(i, j, k)]);
            /* y - 1/2 face*/ 
            uface = 0.5 * (u[IDX(i, j, k)].y + u[IDX(i, jm, k)].y);
             if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, jm, k)] + coeff_UU * conc[IDX(i, jmm, k)] + coeff_D * conc[IDX(i, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, jp, k)] + coeff_D * conc[IDX(i, jm, k)]);
            
            /* z + 1/2 face */
            uface = 0.5 * (u[IDX(i, j, k)].z + u[IDX(i, j, kp)].z);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, j, km)] + coeff_D * conc[IDX(i, j, kp)]);
            else
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, kp)] + coeff_UU * conc[IDX(i, j, kpp)] + coeff_D * conc[IDX(i, j, k)]);
            /* z - 1/2 face */ 
            uface = 0.5 * (u[IDX(i, j, k)].z + u[IDX(i, j, km)].z);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, km)] + coeff_UU * conc[IDX(i, j, kmm)] + coeff_D * conc[IDX(i, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, j, kp)] + coeff_D * conc[IDX(i, j, km)]);		    
#endif /* QUICK */

#ifdef KT
          /*at face H_i+1/2, H_j+1/2, H_k+1/2*/

          ux1 = minmod((u[IDX(i,  j,  k )].x - u[IDX(im, j,  k )].x) / dx, (u[IDX(ip,  j,   k  )].x - u[IDX(i,  j,  k )].x)/ dx);
          ux2 = minmod((u[IDX(ip, j,  k )].x - u[IDX(i,  j,  k )].x) / dx, (u[IDX(ipp, j,   k  )].x - u[IDX(ip, j,  k )].x)/ dx);
          uy1 = minmod((u[IDX(i,  j,  k )].y - u[IDX(i,  jm, k )].y) / dy, (u[IDX(i,   jp,  k  )].y - u[IDX(i,  j,  k )].y)/ dy);
          uy2 = minmod((u[IDX(i,  jp, k )].y - u[IDX(i,  j,  k )].y) / dy, (u[IDX(i,   jpp, k  )].y - u[IDX(i,  jp, k )].y)/ dy);
          uz1 = minmod((u[IDX(i,  j,  k )].z - u[IDX(i,  j,  km)].z) / dz, (u[IDX(i,   j,   kp )].z - u[IDX(i,  j,  k )].z)/ dz);
          uz2 = minmod((u[IDX(i,  j,  kp)].z - u[IDX(i,  j,  k )].z) / dz, (u[IDX(i,   j,   kpp)].z - u[IDX(i,  j,  kp)].z)/ dz);
          
          /*Velocity at x+1/2 slope limiting*/

          uL = u[IDX(i,  j, k)].x + 0.5 * dx * ux1;
          uR = u[IDX(ip, j, k)].x - 0.5 * dx * ux2;

          vL = u[IDX(i, j,  k)].y + 0.5 * dy * uy1;
          vR = u[IDX(i, jp, k)].y - 0.5 * dy * uy2;

          wL = u[IDX(i, j, k )].z + 0.5 * dz * uz1;
          wR = u[IDX(i, j, kp)].z - 0.5 * dz * uz2;

          cx1 = minmod((conc[IDX(i,  j,  k )] - conc[IDX(im, j,  k )]) / dx, (conc[IDX(ip,  j,   k  )] - conc[IDX(i,  j,  k )])/ dx);
          cx2 = minmod((conc[IDX(ip, j,  k )] - conc[IDX(i,  j,  k )]) / dx, (conc[IDX(ipp, j,   k  )] - conc[IDX(ip, j,  k )])/ dx);
          cy1 = minmod((conc[IDX(i,  j,  k )] - conc[IDX(i,  jm, k )]) / dy, (conc[IDX(i,   jp,  k  )] - conc[IDX(i,  j,  k )])/ dy);
          cy2 = minmod((conc[IDX(i,  jp, k )] - conc[IDX(i,  j,  k )]) / dy, (conc[IDX(i,   jpp, k  )] - conc[IDX(i,  jp, k )])/ dy);
          cz1 = minmod((conc[IDX(i,  j,  k )] - conc[IDX(i,  j,  km)]) / dz, (conc[IDX(i,   j,   kp )] - conc[IDX(i,  j,  k )])/ dz);
          cz2 = minmod((conc[IDX(i,  j,  kp)] - conc[IDX(i,  j,  k )]) / dz, (conc[IDX(i,   j,   kpp)] - conc[IDX(i,  j,  kp)])/ dz);

          cL = conc[IDX(i,  j, k)] + 0.5 * dx * cx1;
          cR = conc[IDX(ip, j, k)] - 0.5 * dx * cx2;

          dL = conc[IDX(i, j, k )] + 0.5 * dy * cy1;
          dR = conc[IDX(i, jp, k)] - 0.5 * dy * cy2;

          eL = conc[IDX(i, j, k )] + 0.5 * dz * cz1;
          eR = conc[IDX(i, j, kp)] - 0.5 * dz * cz2;

          fL1 = cL * uL;
          fR1 = cR * uR;

          fL2 = dL * vL;
          fR2 = dR * vR;

          fL3 = eL * wL;
          fR3 = eR * wR;

          a = fmax(fabs(uL), fabs(uR));
          b = fmax(fabs(vL), fabs(vR));
          c = fmax(fabs(wL), fabs(wR));

          /*Numerical fluxes at x+1/2*/
          H1 = 0.5 * (fL1 + fR1 - a * (cR - cL));
          H3 = 0.5 * (fL2 + fR2 - b * (dR - dL));
          H5 = 0.5 * (fL3 + fR3 - c * (eR - eL));

         /*at face H_i-1/2, H_j-1/2, H_k-1/2*/

          ux1 = minmod((u[IDX(im, j,  k )].x - u[IDX(imm, j,   k  )].x) / dx, (u[IDX(i,   j,  k )].x - u[IDX(im, j,   k )].x)/ dx);
          ux2 = minmod((u[IDX(i,  j,  k )].x - u[IDX(im,  j,   k  )].x) / dx, (u[IDX(ip,  j,  k )].x - u[IDX(i,  j,   k )].x)/ dx);
          uy1 = minmod((u[IDX(i,  jm, k )].y - u[IDX(i,   jmm, k  )].y) / dy, (u[IDX(i,   j,  k )].y - u[IDX(i,  jm,  k )].y)/ dy);
          uy2 = minmod((u[IDX(i,  j,  k )].y - u[IDX(i,   jm,  k  )].y) / dy, (u[IDX(i,   jp, k )].y - u[IDX(i,  j,   k )].y)/ dy);
          uz1 = minmod((u[IDX(i,  j,  km)].z - u[IDX(i,   j,   kmm)].z) / dz, (u[IDX(i,   j,  k )].z - u[IDX(i,  j,   km)].z)/ dz);
          uz2 = minmod((u[IDX(i,  j,  k )].z - u[IDX(i,   j,   km )].z) / dz, (u[IDX(i,   j,  kp)].z - u[IDX(i,  j,   k )].z)/ dz);
          
          /*Velocity at x-1/2 slope limiting*/

          uL = u[IDX(im, j, k )].x + 0.5 * dx * ux1;
          uR = u[IDX(i,  j, k )].x - 0.5 * dx * ux2;

          vL = u[IDX(i, jm, k )].y + 0.5 * dy * uy1;
          vR = u[IDX(i, j,  k )].y - 0.5 * dy * uy2;

          wL = u[IDX(i, j,  km)].z + 0.5 * dz * uz1;
          wR = u[IDX(i, j,  k )].z - 0.5 * dz * uz2;

          cx1 = minmod((conc[IDX(im, j,  k )] - conc[IDX(imm, j,   k  )]) / dx, (conc[IDX(i,   j,  k )] - conc[IDX(im, j,  k )])/ dx);
          cx2 = minmod((conc[IDX(i,  j,  k )] - conc[IDX(im,  j,   k  )]) / dx, (conc[IDX(ip,  j,  k )] - conc[IDX(i,  j,  k )])/ dx);
          cy1 = minmod((conc[IDX(i,  jm, k )] - conc[IDX(i,   jmm, k  )]) / dy, (conc[IDX(i,   j,  k )] - conc[IDX(i,  jm, k )])/ dy);
          cy2 = minmod((conc[IDX(i,  j,  k )] - conc[IDX(i,   jm,  k  )]) / dy, (conc[IDX(i,   jp, k )] - conc[IDX(i,  j,  k )])/ dy);
          cz1 = minmod((conc[IDX(i,  j,  km)] - conc[IDX(i,   j,   kmm)]) / dz, (conc[IDX(i,   j,  k )] - conc[IDX(i,  j,  km)])/ dz);
          cz2 = minmod((conc[IDX(i,  j,  k )] - conc[IDX(i,   j,   km )]) / dz, (conc[IDX(i,   j,  kp)] - conc[IDX(i,  j,  k )])/ dz);

          cL = conc[IDX(im, j, k)] + 0.5 * dx * cx1;
          cR = conc[IDX(i,  j, k)] - 0.5 * dx * cx2;

          dL = conc[IDX(i, jm, k)] + 0.5 * dy * cy1;
          dR = conc[IDX(i, j,  k)] - 0.5 * dy * cy2;

          eL = conc[IDX(i, j, km)] + 0.5 * dz * cz1;
          eR = conc[IDX(i, j, k )] - 0.5 * dz * cz2;

          fL1 = cL * uL;
          fR1 = cR * uR;

          fL2 = dL * vL;
          fR2 = dR * vR;

          fL3 = eL * wL;
          fR3 = eR * wR;

          a = fmax(fabs(uL), fabs(uR));
          b = fmax(fabs(vL), fabs(vR));
          c = fmax(fabs(wL), fabs(wR));

          /*Numerical fluxes at x-1/2*/
          H2 = 0.5 * (fL1 + fR1 - a * (cR - cL));
          H4 = 0.5 * (fL2 + fR2 - b * (dR - dL));
          H6 = 0.5 * (fL3 + fR3 - c * (eR - eL));

         /* Final formula of KT scheme*/
          rhs_conc[IDX(i,j,k)] =  -(H1-H2)/dx -(H3-H4)/dy -(H5-H6)/dz;
#endif

#ifdef WENO

    // Lax-Friedrichs Flux Splitting
    //a.x  = fmax(fabs(u[IDX(i,j,k)].x));
    //a.y  = fmax(fabs(u[IDX(i,j,k)].y));
    //a.z  = fmax(fabs(u[IDX(i,j,k)].z));
    a.x  = fabs(u[IDX(i,j,k)].x);
    a.y  = fabs(u[IDX(i,j,k)].y);
    a.z  = fabs(u[IDX(i,j,k)].z);

    uL.x = 0.5 * ( conc[IDX(i,  j, k )]*u[IDX(i, j, k )].x + a.x * conc[IDX(i,  j, k  )]);
    uR.x = 0.5 * ( conc[IDX(ip, j, k )]*u[IDX(ip,j, k )].x - a.x * conc[IDX(ip, j, k  )]);

    uL.y = 0.5 * ( conc[IDX(i, j,  k )]*u[IDX(i, j, k )].y + a.y * conc[IDX(i,  j,  k )]);
    uR.y = 0.5 * ( conc[IDX(i, jp, k )]*u[IDX(i, jp,k )].y - a.y * conc[IDX(i,  jp, k )]);

    uL.z = 0.5 * ( conc[IDX(i, j,  k )]*u[IDX(i, j, k )].z + a.z * conc[IDX(i,  j,  k )]);
    uR.z = 0.5 * ( conc[IDX(i, j,  kp)]*u[IDX(i, j, kp)].z - a.z * conc[IDX(i,  j,  kp)]);

    // Right Flux
    uLmm.x = 0.5 * ( conc[IDX(imm, j, k)]*u[IDX(imm, j, k)].x + a.x * conc[IDX(imm, j, k)]);
    uLm.x  = 0.5 * ( conc[IDX(im,  j, k)]*u[IDX(im,  j, k)].x + a.x * conc[IDX(im,  j, k)]);
    uLp.x  = 0.5 * ( conc[IDX(ip,  j, k)]*u[IDX(ip,  j, k)].x + a.x * conc[IDX(ip,  j, k)]);
    uLpp.x = 0.5 * ( conc[IDX(ipp, j, k)]*u[IDX(ipp, j, k)].x + a.x * conc[IDX(ipp, j, k)]);

    uLmm.y = 0.5 * ( conc[IDX(i, jmm, k)]*u[IDX(i, jmm, k)].y + a.y * conc[IDX(i, jmm,  k)]);
    uLm.y  = 0.5 * ( conc[IDX(i, jm,  k)]*u[IDX(i, jm,  k)].y + a.y * conc[IDX(i, jm,   k)]);
    uLp.y  = 0.5 * ( conc[IDX(i, jp,  k)]*u[IDX(i, jp,  k)].y + a.y * conc[IDX(i, jp,   k)]);
    uLpp.y = 0.5 * ( conc[IDX(i, jpp, k)]*u[IDX(i, jpp, k)].y + a.y * conc[IDX(i, jpp,  k)]);

    uLmm.z = 0.5 * ( conc[IDX(i, j, kmm)]*u[IDX(i, j, kmm)].z + a.z * conc[IDX(i, j, kmm)]);
    uLm.z  = 0.5 * ( conc[IDX(i, j, km )]*u[IDX(i, j, km )].z + a.z * conc[IDX(i, j, km)]);
    uLp.z  = 0.5 * ( conc[IDX(i, j, kp )]*u[IDX(i, j, kp )].z + a.z * conc[IDX(i, j, kp)]);
    uLpp.z = 0.5 * ( conc[IDX(i, j, kpp)]*u[IDX(i, j, kpp)].z + a.z * conc[IDX(i, j, kpp)]);

    // flux calculations
    p0n.x = (2 * uLmm.x - 7 * uLm.x + 11 * uL.x) / 6;
    p1n.x = (-uLm.x + 5 * uL.x + 2 * uLp.x) / 6;
    p2n.x = (2 * uL.x + 5 * uLp.x - uLpp.x) / 6;

    p0n.y = (2 * uLmm.y - 7 * uLm.y + 11 * uL.y) / 6;
    p1n.y = (-uLm.y + 5 * uL.y + 2 * uLp.y) / 6;
    p2n.y = (2 * uL.y + 5 * uLp.y - uLpp.y) / 6;

    p0n.z = (2 * uLmm.z - 7 * uLm.z + 11 * uL.z) / 6;
    p1n.z = (-uLm.z + 5 * uL.z + 2 * uLp.z) / 6;
    p2n.z = (2 * uL.z + 5 * uLp.z - uLpp.z) / 6;

    // smoothness indicators 
    B0n.x = (13 / 12) * pow((uLmm.x - 2 * uLm.x + uL.x),2) + (1 / 4) * pow((uLmm.x - 4 * uLm.x + 3 * uL.x),2);
    B1n.x = (13 / 12) * pow((uLm.x - 2 * uL.x + uLp.x),2) + (1 / 4) * pow((uLm.x - uLp.x),2);
    B2n.x = (13 / 12) * pow((uL.x - 2 * uLp.x + uLpp.x),2) + (1 / 4) * pow((3 * uL.x - 4 * uLp.x + uLpp.x),2);

    B0n.y = (13 / 12) * pow((uLmm.y - 2 * uLm.y + uL.y),2) + (1 / 4) * pow((uLmm.y - 4 * uLm.y + 3 * uL.y),2);
    B1n.y = (13 / 12) * pow((uLm.y - 2 * uL.y + uLp.y),2) + (1 / 4) * pow((uLm.y - uLp.y),2);
    B2n.y = (13 / 12) * pow((uL.y - 2 * uLp.y + uLpp.y),2) + (1 / 4) * pow((3 * uL.y - 4 * uLp.y + uLpp.y),2);

    B0n.z = (13 / 12) * pow((uLmm.z - 2 * uLm.z + uL.z),2) + (1 / 4) * pow((uLmm.z - 4 * uLm.z + 3 * uL.z),2);
    B1n.z = (13 / 12) * pow((uLm.z - 2 * uL.z + uLp.z),2) + (1 / 4) * pow((uLm.z - uLp.z),2);
    B2n.z = (13 / 12) * pow((uL.z - 2 * uLp.z + uLpp.z),2) + (1 / 4) * pow((3 * uL.z - 4 * uLp.z + uLpp.z),2);

    d0n = 1. / 10.;
    d1n = 6. / 10.;
    d2n = 3. / 10.;
    epsilon = 1e-6;

    alpha0n.x = d0n / pow((epsilon + B0n.x),2);
    alpha1n.x = d1n / pow((epsilon + B1n.x),2);
    alpha2n.x = d2n / pow((epsilon + B2n.x),2);
    alphasumn.x = alpha0n.x + alpha1n.x + alpha2n.x;

    alpha0n.y = d0n / pow((epsilon + B0n.y),2);
    alpha1n.y = d1n / pow((epsilon + B1n.y),2);
    alpha2n.y = d2n / pow((epsilon + B2n.y),2);
    alphasumn.y = alpha0n.y + alpha1n.y + alpha2n.y;

    alpha0n.z = d0n / pow((epsilon + B0n.z),2);
    alpha1n.z = d1n / pow((epsilon + B1n.z),2);
    alpha2n.z = d2n / pow((epsilon + B2n.z),2);
    alphasumn.z = alpha0n.z + alpha1n.z + alpha2n.z;

    w0n.x = alpha0n.x / alphasumn.x;
    w1n.x = alpha1n.x / alphasumn.x;
    w2n.x = alpha2n.x / alphasumn.x;

    w0n.y = alpha0n.y / alphasumn.y;
    w1n.y = alpha1n.y / alphasumn.y;
    w2n.y = alpha2n.y / alphasumn.y;

    w0n.z = alpha0n.z / alphasumn.z;
    w1n.z = alpha1n.z / alphasumn.z;
    w2n.z = alpha2n.z / alphasumn.z;

    hn.x = w0n.x * p0n.x + w1n.x * p1n.x + w2n.x * p2n.x;
    hn.y = w0n.y * p0n.y + w1n.y * p1n.y + w2n.y * p2n.y;
    hn.z = w0n.z * p0n.z + w1n.z * p1n.z + w2n.z * p2n.z;

    // Left Flux
    uRmm.x = 0.5 * ( conc[IDX(im,   j, k)]*u[IDX(im,   j, k)].x - a.x * conc[IDX(im,   j, k)]);
    uRm.x  = 0.5 * ( conc[IDX(i,    j, k)]*u[IDX(i,    j, k)].x - a.x * conc[IDX(i,    j, k)]);
    uRp.x  = 0.5 * ( conc[IDX(ipp,  j, k)]*u[IDX(ipp,  j, k)].x - a.x * conc[IDX(ipp,  j, k)]);
    uRpp.x = 0.5 * ( conc[IDX(ippp, j, k)]*u[IDX(ippp, j, k)].x - a.x * conc[IDX(ippp, j, k)]);

    uRmm.y = 0.5 * ( conc[IDX(i, jm,   k)]*u[IDX(i, jm,   k)].y - a.y * conc[IDX(i, jm,   k)]);
    uRm.y  = 0.5 * ( conc[IDX(i, j,    k)]*u[IDX(i, j,    k)].y - a.y * conc[IDX(i, j,    k)]);
    uRp.y  = 0.5 * ( conc[IDX(i, jpp,  k)]*u[IDX(i, jpp,  k)].y - a.y * conc[IDX(i, jpp,  k)]);
    uRpp.y = 0.5 * ( conc[IDX(i, jppp, k)]*u[IDX(i, jppp, k)].y - a.y * conc[IDX(i, jppp, k)]);

    uRmm.z = 0.5 * ( conc[IDX(i, j, km  )]*u[IDX(i, j, km  )].z - a.z * conc[IDX(i, j, km  )]);
    uRm.z  = 0.5 * ( conc[IDX(i, j, k   )]*u[IDX(i, j, k   )].z - a.z * conc[IDX(i, j, k   )]);
    uRp.z  = 0.5 * ( conc[IDX(i, j, kpp )]*u[IDX(i, j, kpp )].z - a.z * conc[IDX(i, j, kpp )]);
    uRpp.z = 0.5 * ( conc[IDX(i, j, kppp)]*u[IDX(i, j, kppp)].z - a.z * conc[IDX(i, j, kppp)]);
 
    p0p.x = (-uRmm.x + 5 * uRm.x + 2 * uR.x) / 6;
    p1p.x = (2 * uRm.x + 5 * uR.x - uRp.x) / 6;
    p2p.x = (11 * uR.x - 7 * uRp.x + 2 * uRpp.x) / 6;

    p0p.y = (-uRmm.y + 5 * uRm.y + 2 * uR.y) / 6;
    p1p.y = (2 * uRm.y + 5 * uR.y - uRp.y) / 6;
    p2p.y = (11 * uR.y - 7 * uRp.y + 2 * uRpp.y) / 6;

    p0p.z = (-uRmm.z + 5 * uRm.z + 2 * uR.z) / 6;
    p1p.z = (2 * uRm.z + 5 * uR.z - uRp.z) / 6;
    p2p.z = (11 * uR.z - 7 * uRp.z + 2 * uRpp.z) / 6;

    B0p.x = (13 / 12) * pow((uRmm.x - 2 * uRm.x + uR.x),2) + (1 / 4) * pow((uRmm.x - 4 * uRm.x + 3 * uR.x),2);
    B1p.x = (13 / 12) * pow((uRm.x - 2 * uR.x + uRp.x), 2) + (1 / 4) * pow((uRm.x - uRp.x), 2);
    B2p.x = (13 / 12) * pow((uR.x - 2 * uRp.x + uRpp.x),2) + (1 / 4) * pow((3 * uR.x - 4 * uRp.x + uRpp.x), 2);

    B0p.y = (13 / 12) * pow((uRmm.y - 2 * uRm.y + uR.y),2) + (1 / 4) * pow((uRmm.y - 4 * uRm.y + 3 * uR.y),2);
    B1p.y = (13 / 12) * pow((uRm.y - 2 * uR.y + uRp.y), 2) + (1 / 4) * pow((uRm.y - uRp.y), 2);
    B2p.y = (13 / 12) * pow((uR.y - 2 * uRp.y + uRpp.y),2) + (1 / 4) * pow((3 * uR.y - 4 * uRp.y + uRpp.y), 2);

    B0p.z = (13 / 12) * pow((uRmm.z - 2 * uRm.z + uR.z),2) + (1 / 4) * pow((uRmm.z - 4 * uRm.z + 3 * uR.z),2);
    B1p.z = (13 / 12) * pow((uRm.z - 2 * uR.z + uRp.z), 2) + (1 / 4) * pow((uRm.z - uRp.z), 2);
    B2p.z = (13 / 12) * pow((uR.z - 2 * uRp.z + uRpp.z),2) + (1 / 4) * pow((3 * uR.z - 4 * uRp.z + uRpp.z), 2);

    d0p = 3 / 10;
    d1p = 6 / 10;
    d2p = 1 / 10;
    epsilon = 1e-6;

    alpha0p.x = d0p / pow((epsilon + B0p.x),2);
    alpha1p.x = d1p / pow((epsilon + B1p.x),2);
    alpha2p.x = d2p / pow((epsilon + B2p.x),2);
    alphasump.x = alpha0p.x + alpha1p.x + alpha2p.x;

    alpha0p.y = d0p / pow((epsilon + B0p.y),2);
    alpha1p.y = d1p / pow((epsilon + B1p.y),2);
    alpha2p.y = d2p / pow((epsilon + B2p.y),2);
    alphasump.y = alpha0p.y + alpha1p.y + alpha2p.y;

    alpha0p.z = d0p / pow((epsilon + B0p.z),2);
    alpha1p.z = d1p / pow((epsilon + B1p.z),2);
    alpha2p.z = d2p / pow((epsilon + B2p.z),2);
    alphasump.z = alpha0p.z + alpha1p.z + alpha2p.z;

    w0p.x = alpha0p.x / alphasump.x;
    w1p.x = alpha1p.x / alphasump.x;
    w2p.x = alpha2p.x / alphasump.x;

    w0p.y = alpha0p.y / alphasump.y;
    w1p.y = alpha1p.y / alphasump.y;
    w2p.y = alpha2p.y / alphasump.y;

    w0p.z = alpha0p.z / alphasump.z;
    w1p.z = alpha1p.z / alphasump.z;
    w2p.z = alpha2p.z / alphasump.z;

    hp.x = w0p.x * p0p.x + w1p.x * p1p.x + w2p.x * p2p.x;
    hp.y = w0p.y * p0p.y + w1p.y * p1p.y + w2p.y * p2p.y;
    hp.z = w0p.z * p0p.z + w1p.z * p1p.z + w2p.z * p2p.z;

    // From here all is the same but with -1 in the index for i or j or k in respectively x, y and z components
    // Lax-Friedrichs Flux Splitting
    
    a.x  = fabs(u[IDX(im,j,k)].x);
    a.y  = fabs(u[IDX(i,jm,k)].y);
    a.z  = fabs(u[IDX(i,j,km)].z);

    uL.x = 0.5 * ( conc[IDX(im,  j, k )]*u[IDX(im, j, k )].x + a.x * conc[IDX(im,  j, k  )]);
    uR.x = 0.5 * ( conc[IDX(i, j, k )]*u[IDX(i,j, k )].x - a.x * conc[IDX(i, j, k  )]);

    uL.y = 0.5 * ( conc[IDX(i, jm,  k )]*u[IDX(i, jm, k )].y + a.y * conc[IDX(i,  jm,  k )]);
    uR.y = 0.5 * ( conc[IDX(i, j, k )]*u[IDX(i, j,k )].y - a.y * conc[IDX(i,  j, k )]);

    uL.z = 0.5 * ( conc[IDX(i, j,  km )]*u[IDX(i, j, km )].z + a.z * conc[IDX(i,  j,  km )]);
    uR.z = 0.5 * ( conc[IDX(i, j,  k)]*u[IDX(i, j, k)].z - a.z * conc[IDX(i,  j,  k)]);

    // Right Flux
    uLmm.x = 0.5 * ( conc[IDX(immm, j, k)]*u[IDX(immm, j, k)].x + a.x * conc[IDX(immm, j, k)]);
    uLm.x  = 0.5 * ( conc[IDX(imm,  j, k)]*u[IDX(imm, j,k)].x + a.x * conc[IDX(imm,  j, k)]);
    uLp.x  = 0.5 * ( conc[IDX(i,    j, k)]*u[IDX(i,   j,k)].x + a.x * conc[IDX(i,  j, k)]);
    uLpp.x = 0.5 * ( conc[IDX(ip,   j, k)]*u[IDX(ip,  j,k)].x + a.x * conc[IDX(ip, j, k)]);

    uLmm.y = 0.5 * ( conc[IDX(i, jmmm, k)]*u[IDX(i,jmmm,k)].y + a.y * conc[IDX(i,jmmm,  k)]);
    uLm.y  = 0.5 * ( conc[IDX(i, jmm,  k)]*u[IDX(i,jmm, k)].y + a.y * conc[IDX(i,jmm,   k)]);
    uLp.y  = 0.5 * ( conc[IDX(i, j,    k)]*u[IDX(i,j, k)].y + a.y * conc[IDX(i,j,   k)]);
    uLpp.y = 0.5 * ( conc[IDX(i, jp,   k)]*u[IDX(i,jp,k)].y + a.y * conc[IDX(i,jp,  k)]);

    uLmm.z = 0.5 * ( conc[IDX(i,j,kmmm)]*u[IDX(i,j,kmmm)].z + a.z * conc[IDX(i,j,kmmm)]);
    uLm.z  = 0.5 * ( conc[IDX(i,j,kmm)]*u[IDX(i,j,kmm)].z + a.z * conc[IDX(i,j,kmm)]);
    uLp.z  = 0.5 * ( conc[IDX(i,j,k)]*u[IDX(i,j,k)].z + a.z * conc[IDX(i,j,k)]);
    uLpp.z = 0.5 * ( conc[IDX(i,j,kp)]*u[IDX(i,j,kp)].z + a.z * conc[IDX(i,j,kp)]);

    // flux calculations
    p0n.x = (2 * uLmm.x - 7 * uLm.x + 11 * uL.x) / 6;
    p1n.x = (-uLm.x + 5 * uL.x + 2 * uLp.x) / 6;
    p2n.x = (2 * uL.x + 5 * uLp.x - uLpp.x) / 6;

    p0n.y = (2 * uLmm.y - 7 * uLm.y + 11 * uL.y) / 6;
    p1n.y = (-uLm.y + 5 * uL.y + 2 * uLp.y) / 6;
    p2n.y = (2 * uL.y + 5 * uLp.y - uLpp.y) / 6;

    p0n.z = (2 * uLmm.z - 7 * uLm.z + 11 * uL.z) / 6;
    p1n.z = (-uLm.z + 5 * uL.z + 2 * uLp.z) / 6;
    p2n.z = (2 * uL.z + 5 * uLp.z - uLpp.z) / 6;

    // smoothness indicators 
    B0n.x = (13 / 12) * pow((uLmm.x - 2 * uLm.x + uL.x),2) + (1 / 4) * pow((uLmm.x - 4 * uLm.x + 3 * uL.x),2);
    B1n.x = (13 / 12) * pow((uLm.x - 2 * uL.x + uLp.x),2) + (1 / 4) * pow((uLm.x - uLp.x),2);
    B2n.x = (13 / 12) * pow((uL.x - 2 * uLp.x + uLpp.x),2) + (1 / 4) * pow((3 * uL.x - 4 * uLp.x + uLpp.x),2);

    B0n.y = (13 / 12) * pow((uLmm.y - 2 * uLm.y + uL.y),2) + (1 / 4) * pow((uLmm.y - 4 * uLm.y + 3 * uL.y),2);
    B1n.y = (13 / 12) * pow((uLm.y - 2 * uL.y + uLp.y),2) + (1 / 4) * pow((uLm.y - uLp.y),2);
    B2n.y = (13 / 12) * pow((uL.y - 2 * uLp.y + uLpp.y),2) + (1 / 4) * pow((3 * uL.y - 4 * uLp.y + uLpp.y),2);

    B0n.z = (13 / 12) * pow((uLmm.z - 2 * uLm.z + uL.z),2) + (1 / 4) * pow((uLmm.z - 4 * uLm.z + 3 * uL.z),2);
    B1n.z = (13 / 12) * pow((uLm.z - 2 * uL.z + uLp.z),2) + (1 / 4) * pow((uLm.z - uLp.z),2);
    B2n.z = (13 / 12) * pow((uL.z - 2 * uLp.z + uLpp.z),2) + (1 / 4) * pow((3 * uL.z - 4 * uLp.z + uLpp.z),2);

    d0n = 1 / 10;
    d1n = 6 / 10;
    d2n = 3 / 10;
    epsilon = 1e-6;

    alpha0n.x = d0n / pow((epsilon + B0n.x),2);
    alpha1n.x = d1n / pow((epsilon + B1n.x),2);
    alpha2n.x = d2n / pow((epsilon + B2n.x),2);
    alphasumn.x = alpha0n.x + alpha1n.x + alpha2n.x;

    alpha0n.y = d0n / pow((epsilon + B0n.y),2);
    alpha1n.y = d1n / pow((epsilon + B1n.y),2);
    alpha2n.y = d2n / pow((epsilon + B2n.y),2);
    alphasumn.y = alpha0n.y + alpha1n.y + alpha2n.y;

    alpha0n.z = d0n / pow((epsilon + B0n.z),2);
    alpha1n.z = d1n / pow((epsilon + B1n.z),2);
    alpha2n.z = d2n / pow((epsilon + B2n.z),2);
    alphasumn.z = alpha0n.z + alpha1n.z + alpha2n.z;

    w0n.x = alpha0n.x / alphasumn.x;
    w1n.x = alpha1n.x / alphasumn.x;
    w2n.x = alpha2n.x / alphasumn.x;

    w0n.y = alpha0n.y / alphasumn.y;
    w1n.y = alpha1n.y / alphasumn.y;
    w2n.y = alpha2n.y / alphasumn.y;

    w0n.z = alpha0n.z / alphasumn.z;
    w1n.z = alpha1n.z / alphasumn.z;
    w2n.z = alpha2n.z / alphasumn.z;

    hnm.x = w0n.x * p0n.x + w1n.x * p1n.x + w2n.x * p2n.x;
    hnm.y = w0n.y * p0n.y + w1n.y * p1n.y + w2n.y * p2n.y;
    hnm.z = w0n.z * p0n.z + w1n.z * p1n.z + w2n.z * p2n.z;

    // Left Flux
    uRmm.x = 0.5 * ( conc[IDX(imm, j, k)]*u[IDX(imm, j, k)].x - a.x * conc[IDX(imm, j, k)]);
    uRm.x  = 0.5 * ( conc[IDX(im,  j, k)]*u[IDX(im,  j, k)].x - a.x * conc[IDX(im,  j, k)]);
    uRp.x  = 0.5 * ( conc[IDX(ip,  j, k)]*u[IDX(ip,  j, k)].x - a.x * conc[IDX(ip,  j, k)]);
    uRpp.x = 0.5 * ( conc[IDX(ipp, j, k)]*u[IDX(ipp, j, k)].x - a.x * conc[IDX(ipp, j, k)]);

    uRmm.y = 0.5 * ( conc[IDX(i, jmm, k)]*u[IDX(i, jmm, k)].y - a.y * conc[IDX(i, jmm, k)]);
    uRm.y  = 0.5 * ( conc[IDX(i, jm,  k)]*u[IDX(i, jm,  k)].y - a.y * conc[IDX(i, jm,  k)]);
    uRp.y  = 0.5 * ( conc[IDX(i, jp,  k)]*u[IDX(i, jp,  k)].y - a.y * conc[IDX(i, jp,  k)]);
    uRpp.y = 0.5 * ( conc[IDX(i, jpp, k)]*u[IDX(i, jpp, k)].y - a.y * conc[IDX(i, jpp, k)]);

    uRmm.z = 0.5 * ( conc[IDX(i, j, kmm)]*u[IDX(i, j, kmm)].z - a.z * conc[IDX(i, j, kmm)]);
    uRm.z  = 0.5 * ( conc[IDX(i, j, km )]*u[IDX(i, j, km )].z - a.z * conc[IDX(i, j, km )]);
    uRp.z  = 0.5 * ( conc[IDX(i, j, kp )]*u[IDX(i, j, kp )].z - a.z * conc[IDX(i, j, kp )]);
    uRpp.z = 0.5 * ( conc[IDX(i, j, kpp)]*u[IDX(i, j, kpp)].z - a.z * conc[IDX(i, j, kpp)]);
 
    p0p.x = (-uRmm.x + 5 * uRm.x + 2 * uR.x) / 6;
    p1p.x = (2 * uRm.x + 5 * uR.x - uRp.x) / 6;
    p2p.x = (11 * uR.x - 7 * uRp.x + 2 * uRpp.x) / 6;

    p0p.y = (-uRmm.y + 5 * uRm.y + 2 * uR.y) / 6;
    p1p.y = (2 * uRm.y + 5 * uR.y - uRp.y) / 6;
    p2p.y = (11 * uR.y - 7 * uRp.y + 2 * uRpp.y) / 6;

    p0p.z = (-uRmm.z + 5 * uRm.z + 2 * uR.z) / 6;
    p1p.z = (2 * uRm.z + 5 * uR.z - uRp.z) / 6;
    p2p.z = (11 * uR.z - 7 * uRp.z + 2 * uRpp.z) / 6;

    B0p.x = (13 / 12) * pow((uRmm.x - 2 * uRm.x + uR.x),2) + (1 / 4) * pow((uRmm.x - 4 * uRm.x + 3 * uR.x),2);
    B1p.x = (13 / 12) * pow((uRm.x - 2 * uR.x + uRp.x), 2) + (1 / 4) * pow((uRm.x - uRp.x), 2);
    B2p.x = (13 / 12) * pow((uR.x - 2 * uRp.x + uRpp.x),2) + (1 / 4) * pow((3 * uR.x - 4 * uRp.x + uRpp.x), 2);

    B0p.y = (13 / 12) * pow((uRmm.y - 2 * uRm.y + uR.y),2) + (1 / 4) * pow((uRmm.y - 4 * uRm.y + 3 * uR.y),2);
    B1p.y = (13 / 12) * pow((uRm.y - 2 * uR.y + uRp.y), 2) + (1 / 4) * pow((uRm.y - uRp.y), 2);
    B2p.y = (13 / 12) * pow((uR.y - 2 * uRp.y + uRpp.y),2) + (1 / 4) * pow((3 * uR.y - 4 * uRp.y + uRpp.y), 2);

    B0p.z = (13 / 12) * pow((uRmm.z - 2 * uRm.z + uR.z),2) + (1 / 4) * pow((uRmm.z - 4 * uRm.z + 3 * uR.z),2);
    B1p.z = (13 / 12) * pow((uRm.z - 2 * uR.z + uRp.z), 2) + (1 / 4) * pow((uRm.z - uRp.z), 2);
    B2p.z = (13 / 12) * pow((uR.z - 2 * uRp.z + uRpp.z),2) + (1 / 4) * pow((3 * uR.z - 4 * uRp.z + uRpp.z), 2);

    d0p = 3 / 10;
    d1p = 6 / 10;
    d2p = 1 / 10;
    epsilon = 1e-6;

    alpha0p.x = d0p / pow((epsilon + B0p.x),2);
    alpha1p.x = d1p / pow((epsilon + B1p.x),2);
    alpha2p.x = d2p / pow((epsilon + B2p.x),2);
    alphasump.x = alpha0p.x + alpha1p.x + alpha2p.x;

    alpha0p.y = d0p / pow((epsilon + B0p.y),2);
    alpha1p.y = d1p / pow((epsilon + B1p.y),2);
    alpha2p.y = d2p / pow((epsilon + B2p.y),2);
    alphasump.y = alpha0p.y + alpha1p.y + alpha2p.y;

    alpha0p.z = d0p / pow((epsilon + B0p.z),2);
    alpha1p.z = d1p / pow((epsilon + B1p.z),2);
    alpha2p.z = d2p / pow((epsilon + B2p.z),2);
    alphasump.z = alpha0p.z + alpha1p.z + alpha2p.z;

    w0p.x = alpha0p.x / alphasump.x;
    w1p.x = alpha1p.x / alphasump.x;
    w2p.x = alpha2p.x / alphasump.x;

    w0p.y = alpha0p.y / alphasump.y;
    w1p.y = alpha1p.y / alphasump.y;
    w2p.y = alpha2p.y / alphasump.y;

    w0p.z = alpha0p.z / alphasump.z;
    w1p.z = alpha1p.z / alphasump.z;
    w2p.z = alpha2p.z / alphasump.z;

    hpm.x = w0p.x * p0p.x + w1p.x * p1p.x + w2p.x * p2p.x;
    hpm.y = w0p.y * p0p.y + w1p.y * p1p.y + w2p.y * p2p.y;
    hpm.z = w0p.z * p0p.z + w1p.z * p1p.z + w2p.z * p2p.z;

    // Compute finite volume residual term, df/dx.
    rhs_conc[IDX(i,j,k)] = (hp.x + hn.x - hpm.x - hnm.x) / dx + 
                           (hp.y + hn.y - hpm.y - hnm.y) / dy + 
                           (hp.z + hn.z - hpm.z - hnm.z) / dz;

#endif

          } /* i, j, k */
            
          sendrecv_borders_scalar(rhs_conc);

}/* rhs_conc */

    #endif
#endif


#ifdef EULER_PARTICLE
void time_stepping_scalar(my_double *ff, my_double *rhs_ff, my_double *old_rhs_ff){
    int i,j,k;
    int dt=1.0;
	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {

            /* 1st order Euler */
            //ff[IDX(i,j,k)] += rhs_ff[IDX(i,j,k)];
            
            /* this is 2nd order Adams-Bashforth */
            ff[IDX(i,j,k)] += 1.5*rhs_ff[IDX(i,j,k)] - 0.5*old_rhs_ff[IDX(i,j,k)];

            /* this is 3rd order Adams-Bashforth */
            //ff[IDX(i,j,k)] += (23.*rhs_ff[IDX(i,j,k)] - 16.*old_rhs_ff[IDX(i,j,k)] +  5.*old_old_rhs_conc[IDX(i,j,k)])/12. ;
            
     /* this is 3rd order Runge-Kutta */
    /*
    copy_scalar(conc, old_conc);
    //u = uo - dt * dF
    ff[IDX(i,j,k)] = old_conc[IDX(i,j,k)] - dt*rhs_ff[IDX(i,j,k)];
    compute_rhs_conc();
    //u = 0.75 * uo + 0.25 * (u - dt * dF)
    ff[IDX(i,j,k)] = 0.75 * old_conc[IDX(i,j,k)] + 0.25 * (ff[IDX(i,j,k)] - dt*rhs_ff[IDX(i,j,k)]);
    compute_rhs_conc();
    //u = (uo + 2 * (u - dt * dF)) / 3
    ff[IDX(i,j,k)] = (old_conc[IDX(i,j,k)] + 2.0 *(ff[IDX(i,j,k)] - dt*rhs_ff[IDX(i,j,k)]))/3.;
     */
    
            }
            /* copy borders */
            sendrecv_borders_scalar(ff);
}
#endif




