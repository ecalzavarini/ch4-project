#include "common_object.h"

void initial_conditions() 
{
  int i,j,k, pp;
  my_double nu,gradP;
  my_double fn,kn;
  my_double y;
  my_double L;

#ifdef LB_FLUID
  nu = (my_double)property.nu;
#ifdef LB_FLUID_INITIAL_POISEUILLE
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
	for (pp = 0; pp < NPOP; pp++)  p[IDX(i,j,k)].p[pp] = wgt[pp];
  
#ifdef LB_FLUID_INITIAL_KOLMOGOROV 
    fn=0.01;
    kn=1.0;
    L=(my_double)(NY);

    fn*=3.0; /* to get the correct amplitude */  
    y = (my_double)center_V[IDX(i,j,k)].y;

	/* along x */  
        p[IDX(i,j,k)].p[1]  +=  wgt[1]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[2]  += -wgt[2]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[7]  +=  wgt[7]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[8]  +=  wgt[8]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[9]  += -wgt[9]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[10] += -wgt[10]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[11] +=  wgt[11]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[12] += -wgt[12]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[13] +=  wgt[13]*fn*sin(kn*two_pi*y/L); 
	p[IDX(i,j,k)].p[14] += -wgt[14]*fn*sin(kn*two_pi*y/L); 
#endif  

#ifdef LB_FLUID_INITIAL_POISEUILLE 
    L=(my_double)(NY);
    y = (my_double)center_V[IDX(i,j,k)].y;
    fn=gradP/(2.0*property.nu);

	/* along x */  
        p[IDX(i,j,k)].p[1]  +=  wgt[1]*fn*y*(L-y);
	p[IDX(i,j,k)].p[2]  += -wgt[2]*fn*y*(L-y);
	p[IDX(i,j,k)].p[7]  +=  wgt[7]*fn*y*(L-y);
	p[IDX(i,j,k)].p[8]  +=  wgt[8]*fn*y*(L-y);
	p[IDX(i,j,k)].p[9]  += -wgt[9]*fn*y*(L-y);
	p[IDX(i,j,k)].p[10] += -wgt[10]*fn*y*(L-y);
	p[IDX(i,j,k)].p[11] +=  wgt[11]*fn*y*(L-y);
	p[IDX(i,j,k)].p[12] += -wgt[12]*fn*y*(L-y);
	p[IDX(i,j,k)].p[13] +=  wgt[13]*fn*y*(L-y);
	p[IDX(i,j,k)].p[14] += -wgt[14]*fn*y*(L-y);
#endif  

#endif
 }/* for i,j,k */


   /* communicate borders for populations */
   sendrecv_borders_pop(p);
}
