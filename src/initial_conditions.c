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

#ifdef LB_FLUID

  /* this is the loop to run over the bulk of all the vertices */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

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

}/* for i,j,k */

   /* communicate borders for populations */
   sendrecv_borders_pop(p);
#endif


#ifdef LB_TEMPERATURE
  /* this is the loop to run over the bulk of all the vertices */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 


#ifdef LB_TEMPERATURE_INITIAL_LINEAR
	/* linear temperature gradient */
	y = (my_double)center_V[IDX(i,j,k)].y;
	t[IDX(i,j,k)] = ( (property.T_bot-property.T_ref) - (property.deltaT/NY)*y );
#endif 

#ifdef LB_TEMPERATURE_INITIAL_ADD_PERTURBATION	 
	if(x<NX/2){ t[IDX(i,j,k)] += 1.e-4; }else{ t[IDX(i,j,k)] -= 1.e-4; }
#endif

#ifdef LB_TEMPERATURE_INITIAL_CONSTANT
        /* constant temperature */
        t[IDX(i,j,k)] = property.T_top;
#endif
	/* on the populations */
	for (pp = 0; pp < NPOP; pp++) 
	  g[IDX(i,j,k)].p[pp] = wgt[pp]*t[IDX(i,j,k)];

 }/* for i,j,k */

   /* communicate borders for populations */
   sendrecv_borders_pop(g);
#endif

}
