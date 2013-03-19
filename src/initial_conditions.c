#include "common_object.h"

void initial_conditions() 
{
  int i,j,k, pp;
  my_double nu,Amp_x;
  my_double fn,kn;
  my_double y;
  my_double L;

#ifdef LB_FLUID
  nu = (my_double)property.nu;
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
    y = (my_double)center_V[IDX(i,j,k)].y;

        /* along x */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] +=  6.0*wgt[pp]*c[pp].x*fn*sin(kn*two_pi*y/L);

#endif  

#ifdef LB_FLUID_INITIAL_POISEUILLE 
    L=(my_double)(NY);
    y = (my_double)center_V[IDX(i,j,k)].y;
    Amp_x = (my_double)property.Amp_x;
    fn=-Amp_x*(4.0*nu)*pow(L,-2.0);

        /* along x */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] +=  6.0*wgt[pp]*c[pp].x*fn*y*(y-L);

#endif  

#endif
 }/* for i,j,k */


   /* communicate borders for populations */
   sendrecv_borders_pop(p);
}
