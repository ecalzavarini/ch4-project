#include "common_object.h"

#ifdef NOSLIP
void noslip_yp(pop*f){

  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){
      
if(LNY_END == NY){
	j = LNY+BRD-1; 

       	for(pp=0;pp<NPOP;pp++) p[IDX(i,j+1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
       	for(pp=0;pp<NPOP;pp++) p[IDX(i,j+2,k)].p[pp] = p[IDX(i,j-1,k)].p[inv[pp]];
#endif
}

}
#endif



#ifdef LB_FLUID_BC
void boundary_conditions(){

  int i,j,k,pp;
  vector vel;
  my_double rho;
  pop p_eq,p_neq;


	/* X direction */	
#ifdef LB_FLUID_BC_X

  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){
      for(pp=0;pp<NPOP;pp++){

if(LNX_END == NX){
	i = LNX+BRD-1;

#ifdef LB_FLUID_BC_XP_OUTLET
	if(pp==0){
	  vel.x = 0.0;
	  vel.y = 0.0;
	  vel.z = 0.0;
	    rho = 1.0;
	  p[IDX(i+1,j,k)] = equilibrium_given_velocity(vel,rho);
	}
#else 
#ifdef LB_FLUID_BC_XP_SLIP

	p[IDX(i+1,j,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].z != 0.0 ) p[IDX(i+1,j,k)].p[pp] = p[IDX(i,j,k)].p[pp];
#else
	/* NOSLIP */
	p[IDX(i+1,j,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#endif
#endif

 }/* if */

if(LNX_START == 0){
  i = BRD; 

#ifdef LB_FLUID_BC_XM_INLET
	if(pp==0){
	  vel.x = 0.01;
	  vel.y = 0.0;
	  vel.z = 0.0;
	    rho = 1.0;
	  p[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho);
	}
#else
#ifdef LB_FLUID_BC_XM_SLIP
		
	p[IDX(i-1,j,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].z != 0.0 ) f[IDX(i-1,j,k)].p[pp] = p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
	p[IDX(i-1,j,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#endif
#endif

 }/* if */

      }/* for pp */
    }/* for j,k */
#endif


  /************************************/

	/* Y direction */	
#ifdef LB_FLUID_BC_Y

  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){
      for(pp=0;pp<NPOP;pp++){

if(LNY_END == NY){
	j = LNY+BRD-1; 
#ifdef LB_FLUID_BC_YP_SLIP
	/*
       	if(norm_yp_pop[IDX_Y(i, k)].p[pp] > 0.0){ 
	  f[IDX(i,j+1,k)].p[inv[pp]] = norm_yp_pop[IDX_Y(i, k)].p[pp]*f[IDX(i,j,k)].p[pp];
	  }else{
	  f[IDX(i,j+1,k)].p[pp] = 0.0;
	  }
	*/
	//f[IDX(i,j+1,k)].p[pp] = wgt[pp];
	p[IDX(i,j+1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) p[IDX(i,j+1,k)].p[pp] = p[IDX(i,j,k)].p[pp];
#ifdef METHOD_MYQUICK
       	p[IDX(i,j+2,k)].p[pp] = p[IDX(i,j-1,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) p[IDX(i,j+2,k)].p[pp] = p[IDX(i,j-1,k)].p[pp];
#endif
#else
	/* NOSLIP */
	  vel.x = u[IDX(i,j,k)].x;
	  vel.y = u[IDX(i,j,k)].y;
	  vel.z = u[IDX(i,j,k)].z;
	    rho = 1.0;
	  p_neq = equilibrium_given_velocity(vel,rho);
	  p_neq.p[pp] -= p[IDX(i,j,k)].p[pp];
	  vel.x *= -1.0;
	  vel.y *= -1.0;
	  vel.z *= -1.0;
	  p_eq = equilibrium_given_velocity(vel,rho);
       	  //f[IDX(i,j+1,k)].p[pp] = 0.5*(f_neq.p[inv[pp]]+f_eq.p[pp]);

       	p[IDX(i,j+1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
       	p[IDX(i,j+2,k)].p[pp] = p[IDX(i,j-1,k)].p[inv[pp]];
#endif

#endif
}

if(LNY_START == 0){
  j = BRD; 
#ifdef LB_FLUID_BC_YM_SLIP	
	/*
        if(norm_ym_pop[IDX_Y(i, k)].p[pp] > 0.0){ 
	  f[IDX(i,j-1,k)].p[inv[pp]] = norm_ym_pop[IDX_Y(i,k)].p[pp]*f[IDX(i,j,k)].p[pp];
	  }else{
	  f[IDX(i,j-1,k)].p[pp] = 0.0;
	  }
	*/	
	p[IDX(i,j-1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) p[IDX(i,j-1,k)].p[pp] = p[IDX(i,j,k)].p[pp];
#ifdef METHOD_MYQUICK
       	  p[IDX(i,j-2,k)].p[pp] = p[IDX(i,j+1,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) p[IDX(i,j-2,k)].p[pp] = p[IDX(i,j+1,k)].p[pp];
#endif
#else
	/* NO SLIP */

	  vel.x = u[IDX(i,j,k)].x;
	  vel.y = u[IDX(i,j,k)].y;
	  vel.z = u[IDX(i,j,k)].z;
	    rho = 1.0;
	  p_neq = equilibrium_given_velocity(vel,rho);
	  p_neq.p[pp] -= p[IDX(i,j,k)].p[pp];
	  vel.x *= -1.0;
	  vel.y *= -1.0;
	  vel.z *= -1.0;
	  p_eq = equilibrium_given_velocity(vel,rho);
	  //f[IDX(i,j-1,k)].p[pp] = 0.5*(f_neq.p[inv[pp]]+f_eq.p[pp]);
	  
	  p[IDX(i,j-1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
       	  p[IDX(i,j-2,k)].p[pp] = p[IDX(i,j+1,k)].p[inv[pp]];
#endif
  
#endif


 }

      }/* for pp */
    }/* for i,k */
#endif


	/* Y direction */	
#ifdef LB_TEMPERATURE_BC_Y

  pop g_eq, g_eq_w;
  my_double effDT, rho2;
  my_double T_wall;
  my_double fac;

#ifdef KALYAN_BC
  my_double a,b,c,d,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3;
#endif


  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNY_END == NY){

 	  j = LNY+BRD-1; 

	  rho = t[IDX(i,j,k)];
#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#else
	  T_wall = 0.0;
#endif

	  //effDT = ( (T_wall-property.T_ref) - rho )*2.0 +  rho;
	  //for(pp=0;pp<NPOP;pp++) g[IDX(i,j+1,k)].p[pp] =  (effDT/rho)*g[IDX(i,j,k)].p[pp];
	  fac = 2.0*(T_wall-property.T_ref)/rho - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j+1,k)].p[pp] =  fac*g[IDX(i,j,k)].p[pp];

	  /* imposed heat flux */
	  //for(pp=0;pp<NPOP;pp++) g[IDX(i,j+1,k)].p[pp] =  g[IDX(i,j,k)].p[inv[pp]]+(-property.kappa*property.deltaT/property.SY)*wgt[pp];

#ifdef KALYAN_BC  /* Works only single processor */
	  a = -1.0/interp4_yp[IDX(i,j,k)];
	  b = (property.T_bot-property.T_ref) - ((property.deltaT/LNY)*(my_double)grid_ruler_y[LNY-1]);	
	  c = interp3_yp[IDX(i,j,k)];
	  d = 1.0 - interp3_yp[IDX(i,j,k)] + interp4_yp[IDX(i,j,k)];

	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j+1,k)].p[pp] = (a * (b- (c * t[IDX(i,j-1,k)]) - (d * t[IDX(i,j,k)]))) * wgt[pp];	
#endif 


#ifdef METHOD_MYQUICK
	  
	  rho = t[IDX(i,j-1,k)];
#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#else
	  T_wall = 0.0;
#endif 	
	  //effDT = ( (T_wall-property.T_ref) - rho )*2.0 +  rho;
	  //for(pp=0;pp<NPOP;pp++) g[IDX(i,j+2,k)].p[pp] =   (effDT/rho)*g[IDX(i,j-1,k)].p[pp];
	  fac = 2.0*(T_wall-property.T_ref)/rho - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j+2,k)].p[pp] =  fac*g[IDX(i,j-1,k)].p[pp];

	  /* imposed heat flux */
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j+2,k)].p[pp] =  g[IDX(i,j-1,k)].p[inv[pp]]+(-property.kappa*property.deltaT/property.SY)*wgt[pp];

#ifdef KALYAN_BC
	  a1 = -1.0/interp4_yp[IDX(i,j+1,k)];
	  b1 = property.T_top-property.T_ref;	 
	  c1 = interp3_yp[IDX(i,j+1,k)];
	  d1 = 1.0 - interp3_yp[IDX(i,j+1,k)] + interp4_yp[IDX(i,j+1,k)];

	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j+2,k)].p[pp] = (a1* (b1 - (c1 * t[IDX(i,j,k)]) - (d1*  (a * (b - (c * t[IDX(i,j-1,k)]) - (d * t[IDX(i,j,k)])))))) * wgt[pp];
#endif

#endif
 }

if(LNY_START == 0){

	  j = BRD; 

	  rho  =  t[IDX(i,j,k)]; 
#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif
	  // effDT = ( (T_wall-property.T_ref) - rho )*2.0 +  rho;
	  //effDT = ( (T_wall-property.T_ref) - rho ) +  (T_wall-property.T_ref);
	  //effDT = 2.0*(T_wall-property.T_ref) - rho;
 	  //for(pp=0;pp<NPOP;pp++) g[IDX(i,j-1,k)].p[pp] =  (effDT/rho)*g[IDX(i,j,k)].p[pp];	  
	  fac = 2.0*(T_wall-property.T_ref)/rho - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j-1,k)].p[pp] =  fac*g[IDX(i,j,k)].p[pp];
	  //fprintf(stderr, "1 fac %lf\n", fac );
	  /* imposed heat flux */
	  //for(pp=0;pp<NPOP;pp++) g[IDX(i,j-1,k)].p[pp] =  g[IDX(i,j,k)].p[inv[pp]]+(-property.kappa*property.deltaT/property.SY)*wgt[pp];

     

#ifdef KALYAN_BC
	  a2 = -1.0/interp2_yp[IDX(i,j,k)];
	  b2 = (property.T_bot-property.T_ref) - ((property.deltaT/LNY)*(my_double)grid_ruler_y[1]);	 
	  c2 = interp_yp[IDX(i,j,k)];
	  d2 = 1.0 - interp_yp[IDX(i,j,k)] + interp2_yp[IDX(i,j,k)];

	   for(pp=0;pp<NPOP;pp++) g[IDX(i,j-1,k)].p[pp] = (a2 * ( b2 - (c2 * t[IDX(i,j+1,k)]) - (d2 * t[IDX(i,j,k)]))) * wgt[pp];
#endif


#ifdef METHOD_MYQUICK 
	  
	  rho =  t[IDX(i,j+1,k)];
#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif   
	  //effDT = ( (T_wall-property.T_ref) - rho )*2.0 +  rho;	
	  //effDT = ( (T_wall-property.T_ref) - rho ) +  (T_wall-property.T_ref);
	  //effDT = 2.0*(T_wall-property.T_ref) - rho;
	  //for(pp=0;pp<NPOP;pp++) g[IDX(i,j-2,k)].p[pp] =  (effDT/rho)*g[IDX(i,j+1,k)].p[pp];
	  fac = 2.0*(T_wall-property.T_ref)/rho - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j-2,k)].p[pp] =  fac*g[IDX(i,j+1,k)].p[pp];
	  //fprintf(stderr, "2 fac %lf\n", fac );

	  /* imposed heat flux */
	  //for(pp=0;pp<NPOP;pp++) g[IDX(i,j-2,k)].p[pp] =  g[IDX(i,j+1,k)].p[inv[pp]]+(-property.kappa*property.deltaT/property.SY)*wgt[pp];

#ifdef KALYAN_BC
	  a3 = -1.0/interp2_yp[IDX(i,j-1,k)];
	  b3 = property.T_bot-property.T_ref;	 
	  c3 = interp_yp[IDX(i,j-1,k)];
	  d3 = 1.0 - interp_yp[IDX(i,j-1,k)] + interp2_yp[IDX(i,j-1,k)];

       for(pp=0;pp<NPOP;pp++) g[IDX(i,j-2,k)].p[pp] = (a3 * (b3 - (c3 * t[IDX(i,j,k)]) - (d3 * (a2 * ( b2 - (c2 * t[IDX(i,j+1,k)]) - (d2 * t[IDX(i,j,k)])))))) * wgt[pp];
#endif

#endif
}

      
    }
#endif



  /*****************************************************************************************/
  /* Z direction */	

#ifdef LB_FLUID_BC_Z

  for (j = BRD; j < LNY + BRD; j++) 			
    for (i = BRD; i < LNX + BRD; i++){
      for(pp=0;pp<NPOP;pp++){

	if(LNZ_END == NZ){
	k = LNZ+BRD-1;

#ifdef LB_FLUID_BC_ZP_SLIP

	p[IDX(i,j,k+1)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].x != 0.0 ) p[IDX(i,j,k+1)].p[pp] = p[IDX(i,j,k)].p[pp];
#else
	/* NOSLIP */
	p[IDX(i,j,k+1)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#endif

 }/* if */

	if(LNZ_START == 0){
	  k = BRD; 

#ifdef LB_FLUID_BC_ZM_SLIP
		
	p[IDX(i,j,k-1)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].x != 0.0 ) p[IDX(i,j,k-1)].p[pp] = p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
	p[IDX(i,j,k-1)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#endif

 }/* if */

      }/* for pp */
    }/* for j,i */
#endif


}
#endif
