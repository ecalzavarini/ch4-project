#include "common_object.h"


/* function to copy a scalar field */
void copy_scalar(my_double *f, my_double *f_copy){
  int i,j,k;
  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){

	f_copy[IDX(i,j,k)] = f[IDX(i,j,k)]; 	

	}
}


/* function to copy a vector field */
void copy_vector(vector *f, vector *f_copy){
  int i,j,k;
  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){

	f_copy[IDX(i,j,k)] = f[IDX(i,j,k)]; 	

	}
}



/* IMPORTANT: We always assume here that the grid and time spacing are = 1 */
#ifdef EULER_PARTICLE
    #ifdef EULER_PARTICLE_CONCENTRATION

/* compute the right hand side of the concentration equation */
/*   \partial_t c = - u \nabla c   */
void compute_rhs_c()
{

	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {

             /*  we implement the linear upwind scheme for FV */
            uface =  0.5*( u[IDX(i,j,k)].x + u[IDX(i+1,j,k)].x )
            if( uface  >0){
                rhs_c[IDX(i,j,k)] =  -uface * c[IDX(i,j,k)];
            }else{
                rhs_c[IDX(i,j,k)] =  -uface * c[IDX(i+1,j,k)];
            }
            uface =  0.5*( u[IDX(i,j,k)].x + u[IDX(i+1,j,k)].x )
            if( u[IDX(i,j,k)].y >0){
                rhs_c[IDX(i,j,k)] +=  -u[IDX(i,j,k)].y*( c[IDX(i,j,k)] - c[IDX(i,j-1,k)] );
            }else{
                rhs_c[IDX(i,j,k)] +=  -u[IDX(i,j,k)].y*( c[IDX(i,j+1,k)] - c[IDX(i,j,k)] );
            }
            if( u[IDX(i,j,k)].z >0){
                rhs_c[IDX(i,j,k)] =  -u[IDX(i,j,k)].z*( c[IDX(i,j,k)] - c[IDX(i,j,k-1)] );
            }else{
                rhs_c[IDX(i,j,k)] =  -u[IDX(i,j,k)].z*( c[IDX(i,j,k+1)] - c[IDX(i,j,k)] );
            }

            /*  we implement the linear upwind scheme for FD */
            /*
            if( u[IDX(i,j,k)].x >0){
                rhs_c[IDX(i,j,k)] =  -u[IDX(i,j,k)].x*( c[IDX(i,j,k)] - c[IDX(i-1,j,k)] );
            }else{
                rhs_c[IDX(i,j,k)] =  -u[IDX(i,j,k)].x*( c[IDX(i+1,j,k)] - c[IDX(i,j,k)] );
            }
            if( u[IDX(i,j,k)].y >0){
                rhs_c[IDX(i,j,k)] +=  -u[IDX(i,j,k)].y*( c[IDX(i,j,k)] - c[IDX(i,j-1,k)] );
            }else{
                rhs_c[IDX(i,j,k)] +=  -u[IDX(i,j,k)].y*( c[IDX(i,j+1,k)] - c[IDX(i,j,k)] );
            }
            if( u[IDX(i,j,k)].z >0){
                rhs_c[IDX(i,j,k)] =  -u[IDX(i,j,k)].z*( c[IDX(i,j,k)] - c[IDX(i,j,k-1)] );
            }else{
                rhs_c[IDX(i,j,k)] =  -u[IDX(i,j,k)].z*( c[IDX(i,j,k+1)] - c[IDX(i,j,k)] );
            }
            */
            }/* for i,j,k */
            sendrecv_borders_scalar(rhs_c);

}/* rhs_c */

    #endif
#endif



void time_stepping_scalar(my_double *f, my_double *rhs_f, my_double *old_rhs_f)
{
	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {
            /* this is 2nd order Adams-Bashforth */
            f[IDX(i,j,k)] += 1.5*rhs_f[IDX(i,j,k)] - 0.5*old_rhs_f[IDX(i,j,k)];
    
            }
}



copy_scalar(rhs_c, old_rhs_c);
compute_rhs_c();
time_stepping_scalar(c,rhs_c,old_rhs_c);