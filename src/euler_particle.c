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
     my_double sgn_a = (a > 0) - (a < 0);
     my_double sgn_b = (b > 0) - (b < 0);

     return 0.5 * (sgn_a + sgn_b) * fmin(fabs(a), fabs(b));
     }
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

/* compute the right hand side of the concentration equation */
/*   \partial_t conc = - \nabla \cdot ( u * conc)   */
void compute_rhs_conc(){
 int i,j,k;
 my_double uface;
 my_double uface_xp, uface_yp, uface_zp;
 my_double uface_xm, uface_ym, uface_zm;
 /* quick method coefficients */
 my_double coeff_D, coeff_U, coeff_UU;
 coeff_D = 3.0/8.0;
 coeff_U = 6.0/8.0;
 coeff_UU = -1.0/8.0;
 /* KT method coefficients */
 my_double dx, dy, dz, dt;
 dx = dy = dz = dt = 1;
 my_double uL, uR, vL, vR, wL, wR, uL1, uR1, vL1, vR1, wL1, wR1, cL, cR, dL, dR, eL, eR, cL1, cR1, dL1, dR1, eL1, eR1;
 my_double ux1, ux2, uy1, uy2, uz1, uz2, ux11, uy11, uz11, ux22, uy22, uz22, cx1, cx2, cy1, cy2, cz1, cz2, cx11, cy11, cz11, cx22, cy22, cz22;
 my_double a, b, c, d, e, f, H1, H2, H3, H4, H5, H6;
 my_double fL1, fR1, fL2, fR2, fL3, fR3, fL4, fR4, fL5, fR5, fL6, fR6;
 int ip,jp,kp,im,km,jm,ipp,jpp,kpp,imm,jmm,kmm;

/* we store the previous rhs */
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

               ip  = i+1; jp  = j+1; kp  = k+1;
               im  = i-1; jm  = j-1; km  = k-1;
               ipp = i+2; jpp = j+2; kpp = k+2;
               imm = i-2; jmm = j-2; kmm = k-2;

               rhs_conc[IDX(i,j,k)] = 0.0;

#ifdef UPWIND 
            /*  we implement the linear upwind scheme for the Finite Volume method*/     
            /* NOTE : this scheme is TOO diffusive! */

            /* x + 1/2 face */
            uface =  0.5*( u[IDX(i,j,k)].x + u[IDX(i+1,j,k)].x);
            if( uface  >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i+1,j,k)];

            /* x - 1/2 face */
            uface = 0.5 * (u[IDX(i,j,k)].x + u[IDX(i-1,j,k)].x);
	        if (uface >  0)
		        rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i-1,j,k)];
            else
		        rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];
            
            /* y + 1/2 face */
	       uface =  0.5*( u[IDX(i,j,k)].y + u[IDX(i,j+1,k)].y);
            if( uface >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j+1,k)];
            
            /* y - 1/2 face */
	        uface = 0.5 * (u[IDX(i,j,k)].y + u[IDX(i,j-1,k)].y);
	        if( uface >0)
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j-1,k)];
            else
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];	

            /* z + 1/2 face */
	        uface = 0.5 * (u[IDX(i,j,k)].z + u[IDX(i,j,k+1)].z);
            if( uface >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k+1)];
            
            /* z - 1/2 face */
    	        uface = 0.5 * (u[IDX(i,j,k)].z + u[IDX(i,j,k-1)].z);
	        if (uface > 0)
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k-1)];
	        else
                rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];		    
#endif /* UPWIND */

#ifdef QUICK
        /*  we implement the QUICK scheme for the Finite Volume method*/

            /* x + 1/2 face */
            uface =  0.5*( u[IDX(i,j,k)].x + u[IDX(i+1,j,k)].x);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i-1, j, k)] + coeff_D * conc[IDX(i+1, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i+1, j, k)] + coeff_UU * conc[IDX(i+2, j, k)] + coeff_D * conc[IDX(i, j, k)]);

            /* x - 1/2 face */
            uface = 0.5 * (u[IDX(i,j,k)].x + u[IDX(i-1,j,k)].x);
	        if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i-1, j, k)] + coeff_UU * conc[IDX(i-2, j, k)] + coeff_D * conc[IDX(i, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i+1, j, k)] + coeff_D * conc[IDX(i-1, j, k)]);
            
            /* y + 1/2 face */
            uface = 0.5 * (u[IDX(i, j, k)].y + u[IDX(1, j+1, k)].y);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, j-1, k)] + coeff_D * conc[IDX(i, j+1, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j+1, k)] + coeff_UU * conc[IDX(i, j+2, k)] + coeff_D * conc[IDX(i, j, k)]);
            /* y - 1/2 face*/ 
            uface = 0.5 * (u[IDX(i, j, k)].y + u[IDX(i, j-1, k)].y);
             if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j-1, k)] + coeff_UU * conc[IDX(i, j-2, k)] + coeff_D * conc[IDX(i, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, j+1, k)] + coeff_D * conc[IDX(i, j-1, k)]);
            
            /* z + 1/2 face */
            uface = 0.5 * (u[IDX(i, j, k)].z + u[IDX(i, j, k+1)].z);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, j, k-1)] + coeff_D * conc[IDX(i, j, k+1)]);
            else
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k+1)] + coeff_UU * conc[IDX(i, j, k+2)] + coeff_D * conc[IDX(i, j, k)]);
            /* z - 1/2 face */ 
            uface = 0.5 * (u[IDX(i, j, k)].z + u[IDX(i, j, k-1)].z);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k-1)] + coeff_UU * conc[IDX(i, j, k-2)] + coeff_D * conc[IDX(i, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, j, k+1)] + coeff_D * conc[IDX(i, j, k-1)]);		    
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
          uR = u[IDX(i, j,  k )].x - 0.5 * dx * ux2;

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

          cL = conc[IDX(im,  j, k)] + 0.5 * dx * cx1;
          cR = conc[IDX(i , j, k)] - 0.5 * dx * cx2;

          dL = conc[IDX(i, jm,  k)] + 0.5 * dy * cy1;
          dR = conc[IDX(i, j , k)] - 0.5 * dy * cy2;

          eL = conc[IDX(i, j, km )] + 0.5 * dz * cz1;
          eR = conc[IDX(i, j, k  )] - 0.5 * dz * cz2;

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

         /* final formula of KT scheme*/
          rhs_conc[IDX(i,j,k)] =  ( -(H1-H2)-(H3-H4)-(H5-H6) )/dx;
#endif
          } /* i, j, k */
            
          sendrecv_borders_scalar(rhs_conc);

}/* rhs_conc */

    #endif
#endif


#ifdef EULER_PARTICLE
void time_stepping_scalar(my_double *ff, my_double *rhs_ff, my_double *old_rhs_ff){
    int i,j,k;
	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {
            /* this is 2nd order Adams-Bashforth */
            ff[IDX(i,j,k)] += 1.5*rhs_ff[IDX(i,j,k)] - 0.5*old_rhs_ff[IDX(i,j,k)];
            /* 1st order Euler */
            //ff[IDX(i,j,k)] += rhs_ff[IDX(i,j,k)];
    
            }
            /* copy borders */
            sendrecv_borders_scalar(ff);
}
#endif




