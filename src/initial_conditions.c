#include "common_object.h"

void initial_conditions(int restart) 
{
  int i,j,k, pp;
  my_double nu,Amp_x;
  my_double fn,kn;
  my_double y,x,z;
  my_double L,LX,LY,LZ;
  my_double LY_half;
  my_double val;
  int iy,iyp,iym;
  my_double norm, norm_all;
  char fnamein[128];
  FILE * fin;
  


#ifdef LB_FLUID
  nu = (my_double)property.nu;
#endif

  /* initialize random seed: */
  //  srand( time(NULL) );

#ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_LAYER
  my_double lambda = compute_lambda_from_stefan(property.deltaT*property.specific_heat/property.latent_heat);
  fprintf(stderr, "Lambda = %e\n", lambda);
#endif

#ifdef LB_FLUID

  /* this is the loop to run over the bulk of all the vertices */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

	/* constant density */
	for (pp = 0; pp < NPOP; pp++)  p[IDX(i,j,k)].p[pp] = wgt[pp];


//#if (defined LB_TEMPERATURE_BUOYANCY && defined LB_INITIAL_BAROMETRIC)	
#ifdef LB_TEMPERATURE_BUOYANCY
 #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_LAYER
   L=(my_double)property.SY; 
   y = (my_double)center_V[IDX(i,j,k)].y;
     /* the good one for linear temperature profile*/
   if ( y <= property.layer_height )
   {
    my_double parameter = lambda*y/property.layer_height;
    my_double int_y = y*(erf(parameter) + exp(-parameter*parameter))/(parameter*sqrt(one_pi)*erf(lambda));
    for (pp = 0; pp < NPOP; pp++) 
     //p[IDX(i,j,k)].p[pp] = wgt[pp]* (1-property.beta_t*property.deltaT*(1-erf(parameter)/erf(lambda)));
     p[IDX(i,j,k)].p[pp] = wgt[pp] * exp( (property.beta_t*property.gravity_y/cs2) * ( property.deltaT*y - property.deltaT*int_y ) );
   }
   else
   {
    my_double parameter = lambda;
    my_double int_y = property.layer_height*(erf(parameter) + exp(-parameter*parameter))/(parameter*sqrt(one_pi)*erf(lambda));
    for (pp = 0; pp < NPOP; pp++) 
     //p[IDX(i,j,k)].p[pp] = wgt[pp]* (1-property.beta_t*property.deltaT*(1-erf(parameter)/erf(lambda)));
     p[IDX(i,j,k)].p[pp] = wgt[pp] * exp( (property.beta_t*property.gravity_y/cs2) * ( property.deltaT*property.layer_height - property.deltaT*int_y ) );

   }
 #endif
 #ifdef LB_INITIAL_BAROMETRIC
  /* hydrostatic density profile  (kind of barometric formula) :  
     if velocity is u = 0  ==>  0 = -\grad P + \rho \beta g T \hat{y}
     Since P=\rho cs^2  and  T is a function of height T=T(y)
     we get \rho = \rho_0 exp( \beta g cs^{-2} \Int_0^y T(y) dy )
     rho_0 is here assumed to be = 1 for T=T_mean
  */
   L=(my_double)property.SY; 
   y = (my_double)center_V[IDX(i,j,k)].y;
     /* the good one for linear temperature profile*/
   for (pp = 0; pp < NPOP; pp++) 
     p[IDX(i,j,k)].p[pp] = wgt[pp]* (exp(property.beta_t*property.gravity_y*y*( 0.5*property.deltaT - 0.5*(property.deltaT/L)*y )/cs2 ));
 #endif
 #ifdef LB_INITIAL_BULK
   /*for quadratic profile */
	L=(my_double)property.SY; 
	y = (my_double)center_V[IDX(i,j,k)].y;
        val = (my_double)property.Amp_t*0.5/property.kappa;
   for (pp = 0; pp < NPOP; pp++) 
     p[IDX(i,j,k)].p[pp] = wgt[pp]*(exp( property.beta_t*property.gravity_y*((property.T_bot-property.T_ref)*y + (val*L -  property.deltaT/L)*y*y/2. - val*y*y*y/3.)/cs2 ));
 #endif
 #ifdef LB_INITIAL_CONSTANT_T_TOP
     p[IDX(i,j,k)].p[pp] = wgt[pp]* (exp(property.beta_t*property.gravity_y*y*( -property.T_top )/cs2 ));
 #endif
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
    //Amp_x = (my_double)property.Amp_x;
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
    //fn=-3.0*Amp_x*(4.0*nu)*pow(L,-2.0);
    fn=-Amp_x*4.0*pow(L,-2.0);
  
        /* along x */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] +=  3.0*wgt[pp]*c[pp].x*fn*y*(y-L);
#endif  

#ifdef LB_FLUID_INITIAL_POISEUILLE_HALF 
    L=2.0*(my_double)property.SY; //NY;
    y = (my_double)center_V[IDX(i,j,k)].y;
    Amp_x = (my_double)property.Amp_x;
    //fn=-3.0*Amp_x*(4.0*nu)*pow(L,-2.0);
    fn=-Amp_x*4.0*pow(L,-2.0);
  
        /* along x */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] +=  3.0*wgt[pp]*c[pp].x*fn*y*(y-L);
#endif  


#ifdef LB_FLUID_INITIAL_CHANNEL 
    /* first we load a turbulent channel profile*/
    turbulent_channel_profile();

    LY=(my_double)property.SY;
    LY_half = LY/2.0;
    y = (my_double)center_V[IDX(i,j,k)].y/LY_half;
    if(y>1.0) y=(2.0-y);
    fn = (my_double)property.Amp_x;
    
    //fprintf(stderr,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA fn = %e\n",fn);

    for (iy = 1; iy < 65; iy++){ 
      if( y >= channel_y[iy-1]  &&  y < channel_y[iy]){
       iym = iy-1;
       iyp = iy;
      }
      //fprintf(stderr,"%d %d \n",iym,iyp);
    }  

    val = fn*( channel_u[iym] + (channel_u[iyp] - channel_u[iym])*(y - channel_y[iym])/(channel_y[iyp] - channel_y[iym]) );

    //fprintf(stdout,"%e %e \n",center_V[IDX(i,j,k)].y, val);
        /* Profile along x */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] +=  3.0*wgt[pp]*c[pp].x*val;  //*sqrt(fn/LY)

	 free(channel_y); 
	 free(channel_u); 
#endif  

#ifdef LB_FLUID_INITIAL_PERTURBATION 
    fn=0.2*property.Amp_x;
        /* random component of amplitude 0.1 U_max in each direction */
       	for (pp = 0; pp < NPOP; pp++)  
	  p[IDX(i,j,k)].p[pp] += 3.0*fn*wgt[pp]*( (2.0*myrand()-1.0)*c[pp].x + (2.0*myrand()-1.0)*c[pp].y + (2.0*myrand()-1.0)*c[pp].z );  
#endif 

#ifdef LB_FLUID_INITIAL_ADD_NOISE
	/* add random noise of amplitude +noise_u or -noise_u only on the x velocity component */ 
    val = 2.0*myrand()-1.0;
    #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_LAYER
        y = center_V[IDX(i, j, k)].y;
        if ( y <= property.layer_height )
    #endif
    for (pp = 0; pp < NPOP; pp++)
        p[IDX(i,j,k)].p[pp] += 6.0*wgt[pp]*property.noise_u*c[pp].x*((val > 0) - (val < 0));
#endif

   
#ifdef LB_FLUID_INITIAL_LANDSCAPE
	/* assigning density value 1 in solid regions (landscape = 1) */  
	if(landscape[IDX(i, j, k)] == 1.0){  	
	  for (pp = 0; pp < NPOP; pp++)  p[IDX(i,j,k)].p[pp] = wgt[pp];
	}
#endif


}/* for i,j,k */


 #ifdef LB_FLUID_INITIAL_UNIT_DENSITY
   /* We set the global density to be =1 , this is an optional step */
   /* ATTENTION: in the present formulation it just works for NX=SX, etc. */
   norm = norm_all = 0.0;

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++)
	for (pp = 0; pp < NPOP; pp++){ 
	  norm += p[IDX(i,j,k)].p[pp];
      }

  MPI_Allreduce(&norm, &norm_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  norm_all /= property.NX* property.NY* property.NZ;

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++)
	for (pp = 0; pp < NPOP; pp++){ 
	  p[IDX(i,j,k)].p[pp] /= norm_all;
      }

  /* end of density normalization*/
 #endif

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

	/* first the temperature is put to zero */
	t[IDX(i,j,k)] = 0.0;

 #ifdef LB_TEMPERATURE_INITIAL_LINEAR
	/* linear temperature gradient */
	L=(my_double)property.SY; //LY;
	y = (my_double)center_V[IDX(i,j,k)].y;
	t[IDX(i,j,k)] = ( (property.T_bot-property.T_ref) - (property.deltaT/L)*y );
 #endif 

 #ifdef LB_TEMPERATURE_INITIAL_CONSTANT
  #ifdef LB_TEMPERATURE_INITIAL_CONSTANT_MEAN
        /* constant mean temperature */
	  t[IDX(i,j,k)] = 0.5*(property.T_top + property.T_bot);
  #endif
  #ifdef LB_TEMPERATURE_INITIAL_CONSTANT_BOT
        /* constant bottom temperature */
          t[IDX(i,j,k)] = property.T_bot;
  #endif
  #ifdef LB_TEMPERATURE_INITIAL_CONSTANT_TOP
        /* constant top temperature */
         t[IDX(i,j,k)] = property.T_top;
  #endif
 #endif

 #ifdef LB_TEMPERATURE_FORCING
  #ifdef LB_TEMPERATURE_INITIAL_BULK
	 //if( abs(center_V[IDX(i,j,k)].y - property.SY/2.0)<=1   ) t[IDX(i,j,k)] = property.T_bot; else  t[IDX(i,j,k)] = property.T_top;
	/* Quadratic temperature gradient due to temperature and bulk forcing */
	L=(my_double)property.SY; 
	y = (my_double)center_V[IDX(i,j,k)].y;
        val = (my_double)property.Amp_t*0.5/property.kappa;
	t[IDX(i,j,k)] = ( (property.T_bot-property.T_ref) + (val*L -  property.deltaT/L)*y - val*y*y);
  #endif
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
      val = 1.e-3; 
      if(NZ==1){ 
        if(center_V[IDX(i, j, k)].x<property.SX/2){ t[IDX(i,j,k)] += val; }else{ t[IDX(i,j,k)] -= val; }
      }else{
        //if(center_V[IDX(i, j, k)].x<property.SX/2){ t[IDX(i,j,k)] += val; }else{ t[IDX(i,j,k)] -= val; }
	if(center_V[IDX(i, j, k)].x<property.SX/2 && center_V[IDX(i, j, k)].z<property.SZ/2){ t[IDX(i,j,k)] += val; }else{ t[IDX(i,j,k)] -= val/4.0; }
	//t[IDX(i,j,k)] += 1.e-1*sin(10.0*two_pi*center_V[IDX(i, j, k)].x/property.SX)*sin(10.0*two_pi*center_V[IDX(i, j, k)].z/property.SZ);
      }
 #endif
 #ifdef LB_TEMPERATURE_INITIAL_ADD_NOISE
      val = property.noise_t;
      t[IDX(i,j,k)] += val*(2.0*myrand()-1.0);
 #endif

#ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_LAYER
  /* non-linear temperature gradient */
  y = center_V[IDX(i, j, k)].y; 
  L=(my_double)property.SY; //LY;
  if ( y <= property.layer_height )
    t[IDX(i,j,k)] = ( (property.T_bot-property.T_ref) - (property.deltaT*erf(lambda*y/property.layer_height)/erf(lambda)) );
#endif  
 
#ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_CAVITY_TEMPERATURE
  x = abs(center_V[IDX(i, j, k)].x - property.SX/2.0);
  y = center_V[IDX(i, j, k)].y; 
  z = abs(center_V[IDX(i, j, k)].z - property.SZ/2.0);
 #ifdef GRID_POP_D2Q9  /* 2D simulations */
  if( x < property.cavity_SX/2 && y < property.cavity_SY )
  {
    my_double DT,A,mpx;
    DT = property.T_bot - property.T_top;
    for (int m=1;m<100;m++)
    {
    	mpx = m*one_pi/property.cavity_SX;
    	A = 2*DT*(1-cos(m*one_pi))/(m*one_pi*sinh(mpx*property.cavity_SY));
    	t[IDX(i,j,k)] += A*sin(mpx*(x+property.cavity_SX/2))*sinh(mpx*(property.cavity_SY-y));
    }
  }
 #endif
 #ifdef GRID_POP_D3Q19 /* 3D simulations */
  if( x < property.cavity_SX/2 && y < property.cavity_SY && z < property.cavity_SZ/2)
  {
    my_double DT,mpx,npz,K,A;
    DT = property.T_bot - property.T_top;
    for (int m=1;m<10;m++)
    	for (int n=1;n<10;n++){
    		mpx = m*one_pi/property.cavity_SX;
    		npz = n*one_pi/property.cavity_SZ;
    		K = sqrt(mpx*mpx + npz*npz);
    		A = 4*DT*(1-cos(n*one_pi))*(1-cos(m*one_pi))/(m*n*one_pi*one_pi*sinh(K*property.cavity_SY));
    		t[IDX(i,j,k)] += A*sin(mpx*(x+property.cavity_SX/2))*sin(npz*(z+property.cavity_SZ/2))*sinh(K*(property.cavity_SY-y));
    	}
  }
 #endif
#endif //LB_TEMPERATURE_MELTING_INITIAL_LIQUID_CAVITY

#ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_SEMISPHERE_TEMPERATURE   
  x = center_V[IDX(i, j, k)].x - property.SX/2.0;
  y = center_V[IDX(i, j, k)].y; 
  z = center_V[IDX(i, j, k)].z - property.SZ/2.0;
  val = sqrt(x*x + y*y + z*z);
  if( val < property.cavity_radius ){
    my_double DT,radius;
    DT = property.T_bot - property.T_top;
    radius = val / property.cavity_radius;
    t[IDX(i,j,k)] = property.T_bot;
    for (int m=1;m<100;m++)
    	t[IDX(i,j,k)] += (4*DT/one_pi)*( pow(radius,2*m-1) * sin( (2*m-1) * atan2(-y,sqrt(x*x+z*z)) ) ) / (2*m-1);
  }
#endif //LB_TEMPERATURE_MELTING_INITIAL_LIQUID_SEMISPHERE



	/* on the populations */
	for (pp = 0; pp < NPOP; pp++) 
	  g[IDX(i,j,k)].p[pp] = wgt[pp]*t[IDX(i,j,k)];

 }/* for i,j,k */

   /* communicate borders for populations */
   sendrecv_borders_pop(g);


 #ifdef LB_TEMPERATURE_FORCING
   /* Melting fields */
  #ifdef LB_TEMPERATURE_MELTING
  /* 1 is fluid , 0 is solid */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
   #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID
      liquid_frac[IDX(i, j, k)]=liquid_frac_old[IDX(i, j, k)]=1.0;
    #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_HALF
      if(center_V[IDX(i, j, k)].y > property.SY/2.0) liquid_frac[IDX(i, j, k)]=liquid_frac_old[IDX(i, j, k)]=0.0;
    #endif  
    #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_SEMISPHERE      
      x = center_V[IDX(i, j, k)].x - property.SX/2.0;
      y = center_V[IDX(i, j, k)].y; 
      z = center_V[IDX(i, j, k)].z - property.SZ/2.0;
      val = sqrt(x*x + y*y + z*z);
      if( val > property.cavity_radius ) liquid_frac[IDX(i, j, k)]=liquid_frac_old[IDX(i, j, k)]=0.0;
    #endif //LB_TEMPERATURE_MELTING_INITIAL_LIQUID_SEMISPHERE
    #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_CAVITY      
      x = abs(center_V[IDX(i, j, k)].x - property.SX/2.0);
      y = center_V[IDX(i, j, k)].y;
      z = abs(center_V[IDX(i, j, k)].z - property.SZ/2.0);
      #ifdef GRID_POP_D2Q9
      if( x > property.cavity_SX/2 || y > property.cavity_SY )
        liquid_frac[IDX(i, j, k)] = liquid_frac_old[IDX(i, j, k)] = 0.0;
      #endif
      #ifdef GRID_POP_D3Q19
      if( x > property.cavity_SX/2 || y > property.cavity_SY ||  z > property.cavity_SZ/2)
        liquid_frac[IDX(i, j, k)] = liquid_frac_old[IDX(i, j, k)] = 0.0;
      #endif
    #endif //LB_TEMPERATURE_MELTING_INITIAL_LIQUID_CAVITY
    #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_LAYER
      y = center_V[IDX(i, j, k)].y; 
      if( y > property.layer_height ) liquid_frac[IDX(i, j, k)]=liquid_frac_old[IDX(i, j, k)]=0.0;
    #endif  
   #else
      liquid_frac[IDX(i, j, k)]=liquid_frac_old[IDX(i, j, k)]=0.0;
   #endif
      }  
  #endif
 #endif

#endif /* end of LB_TEMPERATURE */


#ifdef LB_SCALAR
  /* this is the loop to run over the bulk of all the vertices */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

#ifdef LB_SCALAR_INITIAL_LINEAR
	/* linear temperature gradient */
	L=(my_double)property.SY; //LY;
	y = (my_double)center_V[IDX(i,j,k)].y;
	s[IDX(i,j,k)] = ( (property.S_bot-property.S_ref) - (property.deltaS/L)*y );			
#endif 

#ifdef LB_SCALAR_INITIAL_ADD_PERTURBATION	 
      if(NZ==1){
        if(center_V[IDX(i, j, k)].x<property.SX/2){ s[IDX(i,j,k)] += 1.e-2; }else{ s[IDX(i,j,k)] -= 1.e-2; }
      }else{
	if(center_V[IDX(i, j, k)].x<property.SX/2 && center_V[IDX(i, j, k)].z<property.SZ/2){ s[IDX(i,j,k)] += 1.e-1; }else{ s[IDX(i,j,k)] -= 0.25e-1; }
      }
#endif

#ifdef LB_SCALAR_INITIAL_CONSTANT
#ifdef LB_SCALAR_INITIAL_CONSTANT_MEAN
        /* constant mean scalar */
	  s[IDX(i,j,k)] = 0.5*(property.S_top + property.S_bot);
#endif
#ifdef LB_SCALAR_INITIAL_CONSTANT_BOT
        /* constant bottom scalar */
          s[IDX(i,j,k)] = property.S_bot;
#endif
#ifdef LB_SCALAR_INITIAL_CONSTANT_TOP
        /* constant top scalar */
          s[IDX(i,j,k)] = property.S_top;
#endif
#endif


#ifdef LB_SCALAR_INITIAL_BULK
	  if( abs(center_V[IDX(i,j,k)].x - property.SX/2.0)<10   ) s[IDX(i,j,k)] = property.S_bot; else  s[IDX(i,j,k)] = property.S_top;
#endif

	/* on the populations */
	for (pp = 0; pp < NPOP; pp++) 
	  h[IDX(i,j,k)].p[pp] = wgt[pp]*s[IDX(i,j,k)];
	
      }/* for i,j,k */

   /* communicate borders for populations */
   sendrecv_borders_pop(h);
#endif


#ifdef OUTPUT_H5
   /*  if requested in param.in , read the populations from file pop.h5*/
  sprintf(fnamein,"pop.h5");  
  fin = fopen("pop.h5","r"); 
  if(restart && fin != NULL){
    read_pop_h5();
      if(ROOT) fprintf(stderr,"The %s file is present!\n Populations are initialized from file.\n",fnamein);
      fclose(fin);
  }else{   
    if(restart) if(ROOT)  fprintf(stderr,"Warning message -> %s file is missing!\n Populations are initialized from memory.\n",fnamein);
    if(!restart) if(ROOT) fprintf(stderr,"Warning message -> %s file not requested!\n Populations are initialized from memory.\n",fnamein);
  }

 #ifdef LB_TEMPERATURE_MELTING
   /*  if requested melting is active, and resume is asked for, read it from file liquid_frac.h5 */
  sprintf(fnamein,"liquid_frac.h5");  
  fin = fopen("liquid_frac.h5","r"); 
  if(restart && resume_t && fin != NULL){
    read_scalar_h5(liquid_frac,"liquid_frac");
      if(ROOT) fprintf(stderr,"The %s file is present!\n melting interface is initialized from file.\n",fnamein);
      fclose(fin);
  }else{   
      if(ROOT) fprintf(stderr,"Warning message -> %s file is missing!\n melting interface is initialized from memory.\n",fnamein);
  }
 #endif

#ifdef LB_FLUID
 #ifdef LB_FLUID_AFTER_INIT_PERTURBATION 
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	for (pp = 0; pp < NPOP; pp++)  p[IDX(i,j,k)].p[pp] += wgt[pp]*(2.0*myrand()-1.0);
      } 
 #endif
#endif

#ifdef LB_TEMPERATURE
 #ifdef LB_TEMPERATURE_AFTER_INIT_PERTURBATION 
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
       	for (pp = 0; pp < NPOP; pp++)  g[IDX(i,j,k)].p[pp] += wgt[pp]*(2.0*myrand()-1.0);
      } 
 #endif
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



#ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_LAYER
my_double compute_lambda_from_stefan(my_double St)
{
  my_double lambda_new,lambda_old=1;
  my_double F,dF,error=1;
  int count=0;
  while (error>1e-6 && count<100)
  {
    F=sqrt(one_pi)*lambda_old*exp(lambda_old*lambda_old)*erf(lambda_old)-St;
    dF=sqrt(one_pi)*exp(lambda_old*lambda_old)*(2*lambda_old*lambda_old+1)*erf(lambda_old)+2*lambda_old;
    lambda_new=lambda_old - F/dF;
    error=fabs(lambda_new - lambda_old);
    count++;
    lambda_old=lambda_new;
  }
  return lambda_new;
}
#endif


void turbulent_channel_profile(){

int iy;

channel_y  = (my_double*) malloc(sizeof(my_double)*65);
channel_u  = (my_double*) malloc(sizeof(my_double)*65); 
  /*
# Authors: Moser, Kim & Mansour
# Reference: DNS of Turbulent Channel Flow up to Re_tau=590, 1999,
#            Physics of Fluids, vol 11, 943-945.
# Numerical Method: Kim, Moin & Moser, 1987, J. Fluid Mech. vol 177, 133-166 
# Re_tau = 178.12
# Normalization: U_max, h 
# Description: Mean velocities 
# Filename: chan180/profiles/chan180.means 
# Date File Created: Mar 01, 2000 
#---------------- 
# ny = 129,  Re = 178.12  
  */
/*  Number of points  N = 65 */   
/*  y coordinate [0,1]    |   U(y) profile [0,1]  */
channel_y[0]=0.000000e+00; channel_u[0]=0.000000e+00;
channel_y[1]=3.011800e-04; channel_u[1]=2.930933e-03;
channel_y[2]=1.204500e-03; channel_u[2]=1.171685e-02;
channel_y[3]=2.709500e-03; channel_u[3]=2.633572e-02;
channel_y[4]=4.815300e-03; channel_u[4]=4.674881e-02;
channel_y[5]=7.520500e-03; channel_u[5]=7.288673e-02;
channel_y[6]=1.082300e-02; channel_u[6]=1.046282e-01;
channel_y[7]=1.472200e-02; channel_u[7]=1.417354e-01;
channel_y[8]=1.921500e-02; channel_u[8]=1.837714e-01;
channel_y[9]=2.429800e-02; channel_u[9]=2.300148e-01;
channel_y[10]=2.996900e-02; channel_u[10]=2.794000e-01;
channel_y[11]=3.622400e-02; channel_u[11]=3.305448e-01;
channel_y[12]=4.306000e-02; channel_u[12]=3.819026e-01;
channel_y[13]=5.047200e-02; channel_u[13]=4.319545e-01;
channel_y[14]=5.845600e-02; channel_u[14]=4.794328e-01;
channel_y[15]=6.700700e-02; channel_u[15]=5.234140e-01;
channel_y[16]=7.612000e-02; channel_u[16]=5.634118e-01;
channel_y[17]=8.579000e-02; channel_u[17]=5.992569e-01;
channel_y[18]=9.601100e-02; channel_u[18]=6.311131e-01;
channel_y[19]=1.067800e-01; channel_u[19]=6.593082e-01;
channel_y[20]=1.180800e-01; channel_u[20]=6.841156e-01;
channel_y[21]=1.299100e-01; channel_u[21]=7.060270e-01;
channel_y[22]=1.422700e-01; channel_u[22]=7.254248e-01;
channel_y[23]=1.551500e-01; channel_u[23]=7.425824e-01;
channel_y[24]=1.685300e-01; channel_u[24]=7.578821e-01;
channel_y[25]=1.824200e-01; channel_u[25]=7.715972e-01;
channel_y[26]=1.967900e-01; channel_u[26]=7.840555e-01;
channel_y[27]=2.116500e-01; channel_u[27]=7.954210e-01;
channel_y[28]=2.269900e-01; channel_u[28]=8.059669e-01;
channel_y[29]=2.427900e-01; channel_u[29]=8.158571e-01;
channel_y[30]=2.590500e-01; channel_u[30]=8.251462e-01;
channel_y[31]=2.757500e-01; channel_u[31]=8.340528e-01;
channel_y[32]=2.928900e-01; channel_u[32]=8.425223e-01;
channel_y[33]=3.104600e-01; channel_u[33]=8.507185e-01;
channel_y[34]=3.284400e-01; channel_u[34]=8.586416e-01;
channel_y[35]=3.468300e-01; channel_u[35]=8.663461e-01;
channel_y[36]=3.656100e-01; channel_u[36]=8.738867e-01;
channel_y[37]=3.847700e-01; channel_u[37]=8.812633e-01;
channel_y[38]=4.043000e-01; channel_u[38]=8.884760e-01;
channel_y[39]=4.241900e-01; channel_u[39]=8.955248e-01;
channel_y[40]=4.444300e-01; channel_u[40]=9.024097e-01;
channel_y[41]=4.650000e-01; channel_u[41]=9.090760e-01;
channel_y[42]=4.859000e-01; channel_u[42]=9.155784e-01;
channel_y[43]=5.071000e-01; channel_u[43]=9.219168e-01;
channel_y[44]=5.286000e-01; channel_u[44]=9.280914e-01;
channel_y[45]=5.503900e-01; channel_u[45]=9.340473e-01;
channel_y[46]=5.724400e-01; channel_u[46]=9.398394e-01;
channel_y[47]=5.947600e-01; channel_u[47]=9.454128e-01;
channel_y[48]=6.173200e-01; channel_u[48]=9.507677e-01;
channel_y[49]=6.401000e-01; channel_u[49]=9.559040e-01;
channel_y[50]=6.631100e-01; channel_u[50]=9.608765e-01;
channel_y[51]=6.863200e-01; channel_u[51]=9.656303e-01;
channel_y[52]=7.097200e-01; channel_u[52]=9.702202e-01;
channel_y[53]=7.332900e-01; channel_u[53]=9.745369e-01;
channel_y[54]=7.570200e-01; channel_u[54]=9.786897e-01;
channel_y[55]=7.809000e-01; channel_u[55]=9.825146e-01;
channel_y[56]=8.049100e-01; channel_u[56]=9.860117e-01;
channel_y[57]=8.290400e-01; channel_u[57]=9.891809e-01;
channel_y[58]=8.532700e-01; channel_u[58]=9.919677e-01;
channel_y[59]=8.775900e-01; channel_u[59]=9.943719e-01;
channel_y[60]=9.019800e-01; channel_u[60]=9.963936e-01;
channel_y[61]=9.264400e-01; channel_u[61]=9.979783e-01;
channel_y[62]=9.509300e-01; channel_u[62]=9.991257e-01;
channel_y[63]=9.754600e-01; channel_u[63]=9.997814e-01;
channel_y[64]=1.000000e+00; channel_u[64]=1.000000e+00;

#ifdef DEBUG
for (iy = 0; iy < 65; iy++) 
  fprintf(stderr,"%d %e %e\n",iy, channel_y[iy],channel_y[iy]);
#endif

}



void initialization_forcing(){
  int ii;
  my_double LX,LY,LZ, U;

  LX = (my_double)(property.SX);
  LY = (my_double)(property.SY);
  LZ = (my_double)(property.SZ);
  U = 0.1;

#ifdef LB_FLUID_FORCING_HIT  /* initialize HOMOGENEOUS ISOTROPIC TURBULENCE */
  if(ROOT) fprintf(stderr,"initialization of LB_FLUID_FORCING_HIT\n");

  nk = 6;
  if(phi_u == NULL) phi_u  = (vector*) malloc(sizeof(vector)*nk);
  vk  = (vector*) malloc(sizeof(vector)*nk);
  vk2  = (my_double*) malloc(sizeof(my_double)*nk);
  randomization_itime = (int)floor( ( sqrt( pow(LX,2.0) +  pow(LY,2.0) +  pow(LZ,2.0) ) / U ) / property.time_dt );  
  /*  fprintf(stderr,"random %d\n",randomization_itime); */

    /*  k vectors such as k^2 = kx^2 + ky^2 +kz^2 <= 2 */
    vk[0].x=1; vk[0].y=0; vk[0].z=0; 
    vk[1].x=0; vk[1].y=1; vk[1].z=0; 
    vk[2].x=0; vk[2].y=0; vk[2].z=1; 
    vk[3].x=1; vk[3].y=1; vk[3].z=0; 
    vk[4].x=0; vk[4].y=1; vk[4].z=1; 
    vk[5].x=1; vk[5].y=0; vk[5].z=1;
    /* compute the square */
    for (ii=0; ii<nk; ii++) vk2[ii] = vk[ii].x*vk[ii].x + vk[ii].y*vk[ii].y + vk[ii].z*vk[ii].z;

    if(ROOT){
      fprintf(stderr,"Forced modes\n"); 
      for (ii=0; ii<nk; ii++) fprintf(stderr,"%d : %e %e %e squared norm %e\n",ii+1, vk[ii].x,vk[ii].y,vk[ii].z,vk2[ii]); 
    }

    if(!resume_u){
  /* initialize phases */
   for (ii=0; ii<nk; ii++) phi_u[ii].x = phi_u[ii].y = phi_u[ii].z = 0.0;

      if(ROOT){ 
        for (ii=0; ii<nk; ii++){
          phi_u[ii].x = myrand();
          phi_u[ii].y = myrand();
          phi_u[ii].z = myrand();
        }
      }
      MPI_Bcast(phi_u, nk, MPI_vector_type, 0, MPI_COMM_WORLD);
    }
#endif

#ifdef LB_TEMPERATURE_FORCING_HIT
  if(ROOT) fprintf(stderr,"initialization of LB_TEMPERATURE_FORCING_HIT\n");
  nk_t = 6;
  if(phi_t == NULL) phi_t  = (vector*) malloc(sizeof(vector)*nk_t);
  vk_t  = (vector*) malloc(sizeof(vector)*nk_t);
  vk2_t  = (my_double*) malloc(sizeof(my_double)*nk_t);
  randomization_itime_t = (int)floor( ( sqrt( pow(LX,2.0) +  pow(LY,2.0) +  pow(LZ,2.0) ) / U ) / property.time_dt );  
  /*  fprintf(stderr,"random %d\n",randomization_itime_t); */

    /*  k vectors such as k^2 = kx^2 + ky^2 +kz^2 <= 2 */
    vk_t[0].x=1; vk_t[0].y=0; vk_t[0].z=0; 
    vk_t[1].x=0; vk_t[1].y=1; vk_t[1].z=0; 
    vk_t[2].x=0; vk_t[2].y=0; vk_t[2].z=1; 
    vk_t[3].x=1; vk_t[3].y=1; vk_t[3].z=0; 
    vk_t[4].x=0; vk_t[4].y=1; vk_t[4].z=1; 
    vk_t[5].x=1; vk_t[5].y=0; vk_t[5].z=1;
    /* compute the square */
    for (ii=0; ii<nk_t; ii++) vk2_t[ii] = vk_t[ii].x*vk_t[ii].x + vk_t[ii].y*vk_t[ii].y + vk_t[ii].z*vk_t[ii].z;

    if(!resume_t){
  /* initialize phases */
   for (ii=0; ii<nk_t; ii++) phi_t[ii].x = phi_t[ii].y = phi_t[ii].z = 0.0;

      if(ROOT){ 
        for (ii=0; ii<nk_t; ii++){
          phi_t[ii].x = myrand();
          phi_t[ii].y = myrand();
          phi_t[ii].z = myrand();
        }
     }
    MPI_Bcast(phi_t, nk_t, MPI_vector_type, 0, MPI_COMM_WORLD);
    }
#endif

#ifdef LB_SCALAR_FORCING_HIT
  if(ROOT) fprintf(stderr,"initialization of LB_SCALAR_FORCING_HIT\n");
  nk_s = 6;
  if(phi_s == NULL) phi_s  = (vector*) malloc(sizeof(vector)*nk_s);
  vk_s  = (vector*) malloc(sizeof(vector)*nk_s);
  vk2_s  = (my_double*) malloc(sizeof(my_double)*nk_s);
  randomization_itime_s = (int)floor( ( sqrt( pow(LX,2.0) +  pow(LY,2.0) +  pow(LZ,2.0) ) / U ) / property.time_dt );  
  /*  fprintf(stderr,"random %d\n",randomization_itime); */

    /*  k vectors such as k^2 = kx^2 + ky^2 +kz^2 <= 2 */
    vk_s[0].x=1; vk_s[0].y=0; vk_s[0].z=0; 
    vk_s[1].x=0; vk_s[1].y=1; vk_s[1].z=0; 
    vk_s[2].x=0; vk_s[2].y=0; vk_s[2].z=1; 
    vk_s[3].x=1; vk_s[3].y=1; vk_s[3].z=0; 
    vk_s[4].x=0; vk_s[4].y=1; vk_s[4].z=1; 
    vk_s[5].x=1; vk_s[5].y=0; vk_s[5].z=1;
    /* compute the square */
    for (ii=0; ii<nk_s; ii++) vk2_s[ii] = vk_s[ii].x*vk_s[ii].x + vk_s[ii].y*vk_s[ii].y + vk_s[ii].z*vk_s[ii].z;

    if(!resume_s){
  /* initialize phases */
   for (ii=0; ii<nk_s; ii++) phi_s[ii].x = phi_s[ii].y = phi_s[ii].z = 0.0;

   randomization_itime_s = (int)floor( ( sqrt( pow(LX,2.0) +  pow(LY,2.0) +  pow(LZ,2.0) ) / U ) / property.time_dt );  

      if(ROOT){ 
        for (ii=0; ii<nk_s; ii++){
          phi_s[ii].x = myrand();
          phi_s[ii].y = myrand();
          phi_s[ii].z = myrand();
        }
      }
    MPI_Bcast(phi_s, nk_s, MPI_vector_type, 0, MPI_COMM_WORLD);
    }
#endif
} 
