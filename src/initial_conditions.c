#include "common_object.h"

void initial_conditions(int restart) 
{
  int i,j,k, pp;
  my_double nu,Amp_x;
  my_double fn,kn;
  my_double y,x;
  my_double L,LX,LY;


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


#if (defined LB_TEMPERATURE_BUOYANCY && defined LB_INITIAL_BAROMETRIC)	
   L=(my_double)property.SY; //NY;
   y = (my_double)center_V[IDX(i,j,k)].y;
  /* hydrostatic density profile -  barometric formula dP/P = -\rho g dy , P=\rho cs^2 , \rho = \beta g T_{Lin}*/
   for (pp = 0; pp < NPOP; pp++) 
     // p[IDX(i,j,k)].p[pp] = wgt[pp]* (exp(property.beta_t*property.gravity_y*y*( (property.T_bot-property.T_ref) - 0.5*(property.deltaT/L)*y )/cs2 ));
     /* the good one */
     p[IDX(i,j,k)].p[pp] = wgt[pp]* (exp(property.beta_t*property.gravity_y*y*( 0.5*property.deltaT - 0.5*(property.deltaT/L)*y )/cs2 ));
     //p[IDX(i,j,k)].p[pp] = wgt[pp];
     //p[IDX(i,j,k)].p[pp] = wgt[pp]*(exp(property.beta_t*y/cs2));
#endif

  
#ifdef LB_FLUID_INITIAL_VORTICES 
    fn=0.0001;
    kn=10.0;
    LY=(my_double)property.SY; //NY;
    y = (my_double)center_V[IDX(i,j,k)].y;
    LX=(my_double)property.SX; //NX;
    x = (my_double)center_V[IDX(i,j,k)].x;
        
       	for (pp = 0; pp < NPOP; pp++)  
	  //p[IDX(i,j,k)].p[pp] +=  3.0*fn*wgt[pp]*( c[pp].y*sin(0.5*kn*two_pi*x/LX) + c[pp].y*sin(kn*two_pi*x/LX) );
	  p[IDX(i,j,k)].p[pp] +=  3.0*fn*wgt[pp]*( c[pp].y*sin(0.5*kn*two_pi*x/LX) - c[pp].x*cos(kn*two_pi*y/LY) );	
#endif  

#ifdef LB_FLUID_INITIAL_KOLMOGOROV 
    fn=0.01;
    kn=1.0;
    L=(my_double)property.SY; //NY;
    y = (my_double)center_V[IDX(i,j,k)].y;

        /* along x */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] +=  3.0*wgt[pp]*c[pp].x*fn*sin(kn*two_pi*y/L);
#endif  


#ifdef LB_FLUID_INITIAL_POISEUILLE 
    L=(my_double)property.SY; //NY;
    y = (my_double)center_V[IDX(i,j,k)].y;
    Amp_x = (my_double)property.Amp_x;
    fn=-3.0*Amp_x*(4.0*nu)*pow(L,-2.0);

        /* along x */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] +=  3.0*wgt[pp]*c[pp].x*fn*y*(y-L);

#endif  


#ifdef LB_FLUID_INITIAL_PERTURBATION 
    fn=1.e-2;
    kn=16;
    //Lx=(my_double)property.SY; //NX;
    //x = (my_double)center_V[IDX(i,j,k)].x;
    L=(my_double)property.SX; //NX;
    x = (my_double)center_V[IDX(i,j,k)].x;

        /* along x */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] += fn*wgt[pp]*( (2.0*drand48()-1.0)*c[pp].x + (2.0*drand48()-1.0)*c[pp].y + (2.0*drand48()-1.0)*c[pp].z );  
	    //3.0*wgt[pp]*c[pp].y*fn*sin(kn*two_pi*x/L);

#endif 

}/* for i,j,k */


   /* We set the global density to be =1 , this is an optional step */
   my_double norm, norm_all;
   norm = norm_all = 0.0;

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++)
	for (pp = 0; pp < NPOP; pp++){ 
	  norm += p[IDX(i,j,k)].p[pp];
      }

  MPI_Allreduce(&norm, &norm_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

  norm_all /= property.NX* property.NY* property.NZ;
  /*
  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++)
	for (pp = 0; pp < NPOP; pp++){ 
	  p[IDX(i,j,k)].p[pp] /= norm_all;
      }
  */
  /* end of density normalization*/

#ifdef METHOD_LOG
  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++)
	for (pp = 0; pp < NPOP; pp++) 
	  p[IDX(i,j,k)].p[pp]=property.tau_u*log(p[IDX(i,j,k)].p[pp]);
#endif

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
	L=(my_double)property.SY; //LY;
	y = (my_double)center_V[IDX(i,j,k)].y;
	t[IDX(i,j,k)] = ( (property.T_bot-property.T_ref) - (property.deltaT/L)*y );
#endif 

#ifdef LB_TEMPERATURE_INITIAL_CONSTANT
        /* constant temperature */
        t[IDX(i,j,k)] = 0.5*(property.T_top + property.T_bot);
#endif

#ifdef LB_TEMPERATURE_INITIAL_BL
	/* BOUNDARY LAYER INITIALIZATION */
	L=(my_double)property.SY; //LY;
	y = (my_double)center_V[IDX(i,j,k)].y;
	fn=20.; 
	if(y>L/fn && y<L*(1.-1./fn) ) t[IDX(i,j,k)] = 0.5*(property.T_top + property.T_bot);
	if(y<=L/fn)  t[IDX(i,j,k)] =  (property.T_bot) - (0.5*fn*property.deltaT/L)*y ;
	if(y>=L*(1.-1./fn))  t[IDX(i,j,k)] =  (property.deltaT*0.5*(1.-1./fn)*fn) - (0.5*fn*property.deltaT/L)*y ;
	  
#endif

#ifdef LB_TEMPERATURE_INITIAL_SPOT
  /* impose a mean temperature profile , note that bc for temp shall be set to 0 */
      my_double spot;
      spot = pow(center_V[IDX(i,j,k)].x-property.SX/2.0, 2.0)+pow(center_V[IDX(i,j,k)].y-property.SY/2.0, 2.0);
      if( spot < 10.0 ) t[IDX(i,j,k)] = property.T_bot; else  t[IDX(i,j,k)] = property.T_top;
#endif

#ifdef LB_TEMPERATURE_INITIAL_ADD_PERTURBATION	 
      if(NZ==1){
        if(center_V[IDX(i, j, k)].x<property.SX/2){ t[IDX(i,j,k)] += 1.e-2; }else{ t[IDX(i,j,k)] -= 1.e-2; }
      }else{
	if(center_V[IDX(i, j, k)].x<property.SX/2 && center_V[IDX(i, j, k)].z<property.SZ/2){ t[IDX(i,j,k)] += 1.e-1; }else{ t[IDX(i,j,k)] -= 0.25e-1; }
	//t[IDX(i,j,k)] += 1.e-1*sin(10.0*two_pi*center_V[IDX(i, j, k)].x/property.SX)*sin(10.0*two_pi*center_V[IDX(i, j, k)].z/property.SZ);
      }
#endif

	/* on the populations */
	for (pp = 0; pp < NPOP; pp++) 
	  g[IDX(i,j,k)].p[pp] = wgt[pp]*t[IDX(i,j,k)];

 }/* for i,j,k */

   /* communicate borders for populations */
   sendrecv_borders_pop(g);

   /* Melting fields */
#ifdef LB_TEMPERATURE_MELTING
  /* 1 is fluid , 0 is solid */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
#ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID
      liquid_frac[IDX(i, j, k)]=liquid_frac_old[IDX(i, j, k)]=1.0;
#else
      liquid_frac[IDX(i, j, k)]=liquid_frac_old[IDX(i, j, k)]=0.0;
#endif
    }  
#endif

#endif


#ifdef OUTPUT_H5
   if(restart) read_pop_h5();

   //#define PERTURBATE
#ifdef PERTURBATE
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	for (pp = 0; pp < NPOP; pp++)  p[IDX(i,j,k)].p[pp] += wgt[pp]*(2.0*drand48()-1.0);
	for (pp = 0; pp < NPOP; pp++)  g[IDX(i,j,k)].p[pp] += wgt[pp]*(2.0*drand48()-1.0);
      }
#endif

#ifdef LB_FLUID
  sendrecv_borders_pop(p);
#endif
#ifdef LB_TEMPERATURE
  sendrecv_borders_pop(g);
#endif
#ifdef LB_SCALAR
  sendrecv_borders_pop(h);
#endif
#endif

}
