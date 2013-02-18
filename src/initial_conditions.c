#include "common_object.h"

void initial() 
{
  int i,j,k, pp;
  my_double nu,gradP;
  my_double fn,kn;

#ifdef LB_FLUID
  nu = (my_double)property.nu;
#ifdef LB_FLUID_FORCING_POISEUILLE
  gradP = (my_double)property.gradP;
#endif
#endif

  /* initialize random seed: */
  //  srand( time(NULL) );

  /* this is the loop to run over the bulk of all the vertices */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

#ifdef LB_FLUID
	/* constant density */
  for (pp = 0; pp < NPOP; pp++)  p[IDX(y,x)].p[pp] = wgt[pp];
  
#ifdef LB_FLUID_INITIAL_KOLMOGOROV 
    fn=0.00000000001;
    kn=10.0;
    j_global= j-BRD+LNY_START;

	/* along x */  
        p[IDX(i,j,k)].p[1]  +=  wgt[1]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[2]  += -wgt[2]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[7]  +=  wgt[7]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[8]  +=  wgt[8]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[9]  += -wgt[9]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[10] += -wgt[10]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[11] +=  wgt[11]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[12] += -wgt[12]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[13] +=  wgt[13]*fn*sin(kn*2.*3.14*j_global/NY); 
	p[IDX(i,j,k)].p[14] += -wgt[14]*fn*sin(kn*2.*3.14*j_global/NY); 
#endif  

 }/* for i,j,k */

}
