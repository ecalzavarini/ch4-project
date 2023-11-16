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



/* IMPORTANT: We always assume here that the grid and time spacing are = 1 */
#ifdef EULER_PARTICLE
    #ifdef EULER_PARTICLE_CONCENTRATION

/* compute the right hand side of the concentration equation */
/*   \partial_t conc = - \nabla \cdot ( u * conc)   */
void compute_rhs_conc(){
 int i,j,k;
 my_double uface;
 /* quick method coefficients */
 my_double coeff_D, coeff_U, coeff_UU;
 coeff_D = 3.0/8.0;
 coeff_U = 6.0/8.0;
 coeff_UU = -1.0/8.0;

/* we store the previous rhs */
copy_scalar(rhs_conc, old_rhs_conc);
/* apply the boundary condition for the fluid hydrodynamics fields */
boundary_conditions_hydro(); 

	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {

            rhs_conc[IDX(i,j,k)] = 0.0;

           #ifdef UPWIND 
            /*  we implement the linear upwind scheme for the Finite Volume method*/     
            /* NOTE : this scheme is TOO diffusive! */

            /* x + 1/2 face */
            uface = 0.01;
            //uface =  0.5*( u[IDX(i,j,k)].x + u[IDX(i+1,j,k)].x);
            if( uface  >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i+1,j,k)];

            /* x - 1/2 face */
            //uface = 0.5 * (u[IDX(i,j,k)].x + u[IDX(i-1,j,k)].x);
	        if (uface >  0)
		        rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i-1,j,k)];
            else
		        rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];
            
            /* y + 1/2 face */
            uface =  0.01;
	        //uface =  0.5*( u[IDX(i,j,k)].y + u[IDX(i,j+1,k)].y);
            if( uface >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j+1,k)];
            
            /* y - 1/2 face */
	        //uface = 0.5 * (u[IDX(i,j,k)].y + u[IDX(i,j-1,k)].y);
	        if( uface >0)
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j-1,k)];
            else
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];	

            /* z + 1/2 face */
            uface =  0.0;
	        //uface = 0.5 * (u[IDX(i,j,k)].z + u[IDX(i,j,k+1)].z);
            if( uface >0)
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k)];
            else
                rhs_conc[IDX(i,j,k)] +=  -uface * conc[IDX(i,j,k+1)];
            
            /* z - 1/2 face */
    	    //uface = 0.5 * (u[IDX(i,j,k)].z + u[IDX(i,j,k-1)].z);
	        if (uface > 0)
	            rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k-1)];
	        else
                rhs_conc[IDX(i,j,k)] += uface * conc[IDX(i,j,k)];		    
            #endif /* UPWIND */

        /*  we implement the QUICK scheme for the Finite Volume method*/

            //uface = 0.01;
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
            
            //uface = 0.01;
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
            
            uface = 0.0;
            /* z + 1/2 face */
            //uface = 0.5 * (u[IDX(i, j, k)].z + u[IDX(i, j, k+1)].z);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, j, k-1)] + coeff_D * conc[IDX(i, j, k+1)]);
            else
                 rhs_conc[IDX(i, j, k)] += -uface * (coeff_U * conc[IDX(i, j, k+1)] + coeff_UU * conc[IDX(i, j, k+2)] + coeff_D * conc[IDX(i, j, k)]);
            /* z - 1/2 face */ 
            //uface = 0.5 * (u[IDX(i, j, k)].z + u[IDX(i, j, k-1)].z);
            if (uface > 0)
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k-1)] + coeff_UU * conc[IDX(i, j, k-2)] + coeff_D * conc[IDX(i, j, k)]);
            else
                 rhs_conc[IDX(i, j, k)] += uface * (coeff_U * conc[IDX(i, j, k)] + coeff_UU * conc[IDX(i, j, k+1)] + coeff_D * conc[IDX(i, j, k-1)]);		    
         

            }/* for i,j,k */

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
            //f[IDX(i,j,k)] += rhs_f[IDX(i,j,k)];
    
            }
            /* copy borders */
            sendrecv_borders_scalar(ff);
}
#endif




