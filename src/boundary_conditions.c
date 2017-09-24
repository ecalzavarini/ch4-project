#include "common_object.h"


#ifdef LB_FLUID_BC
void boundary_conditions(){

  int i,j,k,pp;
  vector vel;
  my_double rho,fac;
  pop p_eq;
  pop p_neq;

	/* X direction */	
#ifdef LB_FLUID_BC_X

  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){
      for(pp=0;pp<NPOP;pp++){

if(LNX_END == NX){
	i = LNX+BRD-1;
	
#ifdef LB_FLUID_BC_X_P_OUTLET
	
	if(pp==0){
	  vel.x = u[IDX(i, j, k)].x;
	  vel.y = u[IDX(i, j, k)].y;
	  vel.z = u[IDX(i, j, k)].z;
	  rho = dens[IDX(i, j, k)];
	  p[IDX(i+1,j,k)] = equilibrium_given_velocity(vel,rho);
	  //p[IDX(i+1,j,k)] = p[IDX(i,j,k)];
#ifdef METHOD_MYQUICK
	  p[IDX(i+2,j,k)] = p[IDX(i+1,j,k)];
	  //p[IDX(i+2,j,k)] = p[IDX(i-1,j,k)];
#endif
	/*
          rho = 1.0;
	  fac = 2.0*(rho)/dens[IDX(i,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) p[IDX(i+1,j,k)].p[pp] =  fac*p[IDX(i,j,k)].p[pp];
#ifdef METHOD_MYQUICK
	  fac = 2.0*(rho)/dens[IDX(i-1,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) p[IDX(i+2,j,k)].p[pp] =  fac*p[IDX(i-1,j,k)].p[pp];

#endif
	*/
	  
	}	  
#else 
#ifdef LB_FLUID_BC_X_P_SLIP

	p[IDX(i+1,j,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].z != 0.0 ) p[IDX(i+1,j,k)].p[pp] = p[IDX(i,j,k)].p[pp];
#else
	/* NOSLIP */
	p[IDX(i+1,j,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
        p[IDX(i+2,j,k)].p[pp] = p[IDX(i-1,j,k)].p[inv[pp]];
#endif
#endif
#endif

 }/* if */

if(LNX_START == 0){
  i = BRD; 

#ifdef LB_FLUID_BC_X_M_INLET
  /* the following if is to compute the equilibrium only one time */
	if(pp==0){
	  vel.x = property.Amp_x;
	  vel.y = property.Amp_y;
	  vel.z = property.Amp_z;
	    rho = 1.0;
	  p[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho);
#ifdef METHOD_MYQUICK
	  p[IDX(i-2,j,k)] = p[IDX(i-1,j,k)];
#endif
	}
	
#else
#ifdef LB_FLUID_BC_X_M_SLIP
		
	p[IDX(i-1,j,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].z != 0.0 ) p[IDX(i-1,j,k)].p[pp] = p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
	p[IDX(i-1,j,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
        p[IDX(i-2,j,k)].p[pp] = p[IDX(i+1,j,k)].p[inv[pp]];
#endif
#endif
#endif

 }/* if */

      }/* for pp */
    }/* for j,k */
#endif


  /************************************/

	/* Y direction */	
#ifdef LB_FLUID_BC_Y
 for (i = 0; i < LNX + TWO_BRD; i++) 			
    for (k = 0; k < LNZ + TWO_BRD; k++){
      // this is the good one
      //  for (i = BRD; i < LNX + BRD; i++) 			
      //    for (k = BRD; k < LNZ + BRD; k++){
      for(pp=0;pp<NPOP;pp++){

if(LNY_END == NY){
	j = LNY+BRD-1; 

#ifdef LB_FLUID_BC_Y_P_OUTLET
	if(pp==0){
	  vel.x =  u[IDX(i, j, k)].x;
	  vel.y =  u[IDX(i, j, k)].y;
	  vel.z =  u[IDX(i, j, k)].z;
	  rho = dens[IDX(i, j, k)];
	  p[IDX(i,j+1,k)] = equilibrium_given_velocity(vel,rho);
#ifdef METHOD_MYQUICK
          p[IDX(i,j+2,k)] = p[IDX(i,j+1,k)];
#endif
	}
#else

#ifdef LB_FLUID_BC_Y_P_SLIP
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
	/*
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
	*/
       	  //f[IDX(i,j+1,k)].p[pp] = 0.5*(f_neq.p[inv[pp]]+f_eq.p[pp]);
	p[IDX(i,j+1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
       //p[IDX(i,j+1,k)].p[pp] = wgt[pp]*dens[IDX(i,j,k)];
       //if(c[pp].y < 0) p[IDX(i,j+1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
       	p[IDX(i,j+2,k)].p[pp] = p[IDX(i,j-1,k)].p[inv[pp]];
        //p[IDX(i,j+2,k)].p[pp] = wgt[pp]*dens[IDX(i,j,k)];
	//if(c[pp].y < 0) p[IDX(i,j+2,k)].p[pp] = p[IDX(i,j-1,k)].p[inv[pp]];
#endif

#endif
#endif
 }

if(LNY_START == 0){
  j = BRD; 
#ifdef LB_FLUID_BC_Y_M_SLIP	
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
	/*
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
	*/
	  //f[IDX(i,j-1,k)].p[pp] = 0.5*(f_neq.p[inv[pp]]+f_eq.p[pp]);
	  
	  p[IDX(i,j-1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	//if(c[pp].y < 0 ) p[IDX(i,j-1,k)].p[pp] = p[IDX(i,j,k)].p[pp];
	//if(c[pp].y > 0 ) p[IDX(i,j-1,k)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
       	  p[IDX(i,j-2,k)].p[pp] = p[IDX(i,j+1,k)].p[inv[pp]];	
	//if(c[pp].y < 0.0 ) p[IDX(i,j-2,k)].p[pp] =  p[IDX(i,j,k)].p[pp];
	//if(c[pp].y > 0.0 ) p[IDX(i,j-2,k)].p[pp] =  p[IDX(i,j+1,k)].p[inv[pp]];	 
#endif
  
#endif


 }

      }/* for pp */
    }/* for i,k */
#endif



#ifdef LB_TEMPERATURE_BC_KEEP_WITHIN

  my_double T_max, T_min;

  /* Rather DIRTY TRICK , in order to keep the temperature inside the boundaries, in the whole bulk */
  T_max = property.T_bot;
  T_min = property.T_top;
  if( property.T_top > property.T_bot){ T_max = property.T_top;  T_min = property.T_bot;}

  for (i = BRD; i < LNX + BRD; i++) 
    for (j = BRD; j < LNY + BRD; j++)
      for (k = BRD; k < LNZ + BRD; k++){

	if(m(g[IDX(i,j,k)]) < T_min) for(pp=0;pp<NPOP;pp++) g[IDX(i,j,k)].p[pp] *= T_min/m(g[IDX(i,j,k)]);
	if(m(g[IDX(i,j,k)]) > T_max) for(pp=0;pp<NPOP;pp++) g[IDX(i,j,k)].p[pp] *= T_max/m(g[IDX(i,j,k)]); 
      }
#endif


	/* Y direction */	
#ifdef LB_TEMPERATURE
 #ifdef LB_TEMPERATURE_BC_Y

  pop g_eq, g_eq_w;
  my_double effDT, rho2;
  my_double T_wall;
  // my_double fac;

  /************************/

  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNY_END == NY){

 	  j = LNY+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
	  
 #ifdef LB_TEMPERATURE_BC_Y_P_VARIABLE
	  /* half fixed temperature, half  fixed-flux */
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_top; else T_wall = 0.5*property.grad_T_top + t[IDX(i,j,k)] + property.T_ref;
 #endif

 #ifdef LB_TEMPERATURE_BC_Y_P_FLUX
	  /* we fix the temperature gradient at the wall */
	  T_wall = 0.5*property.grad_T_top + t[IDX(i,j,k)] + property.T_ref;
 #endif
	  
#else
	  T_wall = 0.0;
#endif

	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j+1,k)].p[pp] =  fac*g[IDX(i,j,k)].p[pp];	  

#ifdef METHOD_MYQUICK
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j-1,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j+2,k)].p[pp] =  fac*g[IDX(i,j-1,k)].p[pp];
#endif
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
 #ifdef LB_TEMPERATURE_BC_Y_M_VARIABLE
	  //if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/2.0), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 2.0)  T_wall = property.T_bot; else T_wall = property.T_top;
  //  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = 0.0;  
	  /* half fixed temperature, half  fixed-flux */
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = -0.5*property.grad_T_bot + t[IDX(i,j,k)] + property.T_ref;
 #endif

 #ifdef LB_TEMPERATURE_BC_Y_M_FLUX
	  /* we fix the temperature gradient at the wall */
	  T_wall = -0.5*property.grad_T_bot + t[IDX(i,j,k)] + property.T_ref;
 #endif

#else
	  T_wall = 0.0;
#endif

	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j-1,k)].p[pp] =  fac*g[IDX(i,j,k)].p[pp];
     
#ifdef METHOD_MYQUICK 
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j+1,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j-2,k)].p[pp] =  fac*g[IDX(i,j+1,k)].p[pp];
#endif
}

      
    }
 #endif
 
#ifdef LB_TEMPERATURE_BC_X
#ifdef LB_TEMPERATURE_BC_X_NOFLUX
 /* now along X , the default is insultaing BC*/
  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNX_END == NX){

 	  i = LNX+BRD-1; 

	  for(pp=0;pp<NPOP;pp++) g[IDX(i+1,j,k)].p[pp] =  g[IDX(i,j,k)].p[pp];
#ifdef METHOD_MYQUICK	  
	  for(pp=0;pp<NPOP;pp++) g[IDX(i+2,j,k)].p[pp] =  g[IDX(i-1,j,k)].p[pp];
#endif
 }

if(LNX_START == 0){

	  i = BRD; 
	
	  for(pp=0;pp<NPOP;pp++) g[IDX(i-1,j,k)].p[pp] =  g[IDX(i,j,k)].p[pp];     
#ifdef METHOD_MYQUICK 	  
	  for(pp=0;pp<NPOP;pp++) g[IDX(i-2,j,k)].p[pp] =  g[IDX(i+1,j,k)].p[pp];
#endif
 }      
    }
#endif
#endif

#endif


  /*****************************************************************************************/
  /* Z direction */	

#ifdef LB_FLUID_BC_Z

  for (j = BRD; j < LNY + BRD; j++) 			
    for (i = BRD; i < LNX + BRD; i++){
      for(pp=0;pp<NPOP;pp++){

	if(LNZ_END == NZ){
	k = LNZ+BRD-1;

#ifdef LB_FLUID_BC_Z_P_SLIP

	p[IDX(i,j,k+1)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].x != 0.0 ) p[IDX(i,j,k+1)].p[pp] = p[IDX(i,j,k)].p[pp];
#else
	/* NOSLIP */
	p[IDX(i,j,k+1)].p[pp] = p[IDX(i,j,k)].p[inv[pp]];
#endif

 }/* if */

	if(LNZ_START == 0){
	  k = BRD; 

#ifdef LB_FLUID_BC_Z_M_SLIP
		
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


} /* end of void boundary conditions */
#endif

/**********************************************************************************************************/
/* When STREAMING is active bc shall be implemented in a different way*/
//#ifdef LB_FLUID_BC
#ifdef METHOD_STREAMING
void boundary_and_pbc_conditions_for_streaming(){

  int i,j,k,pp;
  vector vel;
  my_double rho;
  pop p_eq,p_neq;
  int ii,jj,kk;	
  my_double fac;

  /* communications to be done in any case  (especially for pbc)*/

#ifdef LB_FLUID
sendrecv_borders_pop(rhs_p);
#endif

#ifdef LB_TEMPERATURE
sendrecv_borders_pop(rhs_g);
//sendrecv_borders_scalar(t);
#endif

#ifdef LB_SCALAR
sendrecv_borders_pop(rhs_h);
#endif

/***************** Bounce-back BC for the complex boundaries  ********************/
  int iii, jjj, kkk;
#ifdef LB_FLUID_FORCING_LANDSCAPE
 for(k=BRD;k<LNZ+BRD;k++){
   for(j=BRD;j<LNY+BRD;j++){
      for(i=BRD;i<LNX+BRD;i++){ 

     for(pp=0;pp<NPOP;pp++){
          iii = i+(int)c[pp].x;
	  jjj = j+(int)c[pp].y;
	  kkk = k+(int)c[pp].z;
	  /* landscape = 1 is a solid region ,  landscape = 0 is a fluid  region */
	  if(landscape[IDX(i, j, k)] == 0.0 && landscape[IDX(iii, jjj, kkk)] == 1.0 ) rhs_p[IDX(iii,jjj,kkk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
     }/* pp */   
      }/* i */
   }/* j */
 }/* k */
#endif

 pop dummy_here, dummy_there;

#ifdef LB_TEMPERATURE_MELTING_BOUNCEBACK
 /* Results given by this bc are not yet satisfactory (total mass is not conserved). Penalization method works better */
 for(k=BRD;k<LNZ+BRD;k++){
   for(j=BRD;j<LNY+BRD;j++){
      for(i=BRD;i<LNX+BRD;i++){ 

     for(pp=0;pp<NPOP;pp++){
          iii = i+(int)c[pp].x;
	  jjj = j+(int)c[pp].y;
	  kkk = k+(int)c[pp].z;

	  /* in place bounce-back for solid boundary regions */
	  if(liquid_frac[IDX(i, j, k)] < 0.5 && liquid_frac[IDX(iii, jjj, kkk)] >= 0.5 ){	 
	    dummy_here.p[inv[pp]] = rhs_p[IDX(i,j,k)].p[inv[pp]];
	    rhs_p[IDX(i,j,k)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
	    rhs_p[IDX(i,j,k)].p[pp] = dummy_here.p[inv[pp]];
	  }  	 
	  
	  /* liquid_frac = 0 is a solid region ,  liquid_frac = 1 is a fluid  region */
	  // good
	  //if(liquid_frac[IDX(i, j, k)] == 0.0 && liquid_frac[IDX(iii, jjj, kkk)] > 0.0 ) rhs_p[IDX(iii,jjj,kkk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
	  //if(liquid_frac[IDX(i, j, k)] < 0.5 && liquid_frac[IDX(iii, jjj, kkk)] >= 0.5 ) rhs_p[IDX(iii,jjj,kkk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
	  // bad
	  //if(liquid_frac[IDX(i, j, k)] > 0.5 && liquid_frac[IDX(iii, jjj, kkk)] <= 0.5 ) rhs_p[IDX(iii,jjj,kkk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
	  //if(liquid_frac[IDX(i, j, k)] == 1.0 && liquid_frac[IDX(iii, jjj, kkk)] < 1.0 ) rhs_p[IDX(iii,jjj,kkk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
     }/* pp */   
      }/* i */
   }/* j */
 }/* k */
#endif

 //#define LB_TEMPERATURE_MELTING_FIXNEGATIVE
#ifdef LB_TEMPERATURE_MELTING_FIXNEGATIVE
 /*  A rather ad-hoc way to fix the negative temperature generation problem */
 for(k=BRD;k<LNZ+BRD;k++){
   for(j=BRD;j<LNY+BRD;j++){
      for(i=BRD;i<LNX+BRD;i++){ 
	/* if the temperature has been decreased below the inital value then raise it up to that minimal value */
	  if(t[IDX(i,j,k)] < property.T_solid)
	    for(pp=0;pp<NPOP;pp++) rhs_g[IDX(i,j,k)].p[pp] *= property.T_solid/m(rhs_g[IDX(i,j,k)]);
	  //}/* pp */   
      }/* i */
   }/* j */
 }/* k */
#endif

/************************************/
	/* X direction */ 
#ifdef LB_FLUID_BC_X

   for (j = 0; j < LNY + TWO_BRD; j++) 			
   for (k = 0; k < LNZ + TWO_BRD; k++){

      for(pp=0;pp<NPOP;pp++){

	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;	

if(LNX_END == NX){

  	i = LNX+BRD-1;

#ifdef LB_FLUID_BC_X_P_OUTLET
	if(pp==0){
	 /*borrowing the density value from the nearest fluid node*/
	  if (j >= BRD && j< LNY+BRD && k >= BRD && k < LNZ+BRD) {rho = dens[IDX(i, j, k)];} else {rho = 1.0;}

	 /*definition of all velocity components by Dirichlet boundary condition method */
#ifdef GRID_POP_D2Q9	  
	  vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[8]+2*(rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[5]))/rho;
	  vel.y = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[6]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[7]))/rho;
#endif

#ifdef GRID_POP_D3Q19
	 vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[18]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[13]))/rho;
	 vel.y = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[7]))/rho;	 
	 vel.z = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[10]+2*(rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[11]))/rho;
#endif
	 /*non-reflecting pressure boundary condition Ref: Finck et al*/
       	 rho = (rho+rho*sqrt(3.0)*(vel.x-u[IDX(i,j,k)].x)+1.0*1.0*1.0)/(1.0+1.0*1.0);

	 rhs_p[IDX(i+1,j,k)] = equilibrium_given_velocity(vel,rho); 
	} 
#else

#ifdef LB_FLUID_BC_X_P_SLIP
	  /* SLIP */
          rhs_p[IDX(i+1,j,k)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#else	 
	  /* NO SLIP */
	  if (jj>= 0 && jj< LNY+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(i+1,jj,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
  
#endif
#endif
}
// }

if(LNX_START == 0){

       i = BRD;

#ifdef LB_FLUID_BC_X_M_OUTLET
	if(pp==0){

	    if (j >= BRD && j< LNY+BRD && k >= BRD && k < LNZ+BRD) {rho = dens[IDX(i, j, k)];} else {rho = 1.0;}

#ifdef GRID_POP_D2Q9
	  vel.x = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[8]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]))/rho;
	  vel.y = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[6]+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]))/rho;
#endif
	  
#ifdef GRID_POP_D3Q19 
	  vel.x = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[18]+2*(rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]))/rho; 
	  vel.y = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[8]))/rho;
	  vel.z = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[10]+2*(rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[13]))/rho;
#endif

      	  rho = (rho+rho*sqrt(3.0)*(vel.x-u[IDX(i,j,k)].x)+1.0*1.0*1.0)/(1.0+1.0*1.0);

	  rhs_p[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho); 
	}
#else
 
#ifdef LB_FLUID_BC_X_M_INLET
	if(pp==0){

#ifdef LB_FLUID_BC_X_M_INLET_POISEUILLE
	  /*Poiseuille flow velocity profile - inlet*/
	     vel.x = -(4.0*(my_double)property.Amp_x*(pow(((my_double)property.SY),-2.0)))*((my_double)center_V[IDX(i,j,k)].y)*(((my_double)center_V[IDX(i,j,k)].y)-((my_double)property.SY));
	  vel.y = 0.0;
          vel.z = 0.0;
#endif

#ifdef LB_FLUID_BC_X_M_INLET_POISEUILLE_HALF
	  /* Half Poiseuille flow velocity profile - inlet*/
	  vel.x = -(4.0*(my_double)property.Amp_x*(pow((2.0*(my_double)property.SY),-2.0)))*((my_double)center_V[IDX(i,j,k)].y)*(((my_double)center_V[IDX(i,j,k)].y)-(2.0*(my_double)property.SY));
	  vel.y = 0.0;
          vel.z = 0.0;
#endif

#ifdef LB_FLUID_BC_X_M_INLET_CONSTANT
	  /* Constant velocity = Amp_x*/        
	  vel.x = property.Amp_x;
	  vel.y = property.Amp_y;
	  vel.z = property.Amp_z;
#endif

#ifdef LB_FLUID_BC_X_M_INLET_JET
	  /*Hyperbolic cosine - inlet*/
	  vel.x = ((my_double)property.Amp_x)/(pow(cosh(0.02*((my_double)center_V[IDX(i,j,k)].y-((my_double)property.SY/2.0))),2.0));
#endif

	  /* If the inlet velocity profile contains only x-velocity component*/   
#ifdef GRID_POP_D2Q9 
	  rho = (1.0/(1.0-vel.x))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[8])+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]));
#endif

#ifdef GRID_POP_D3Q19
	  rho = (1.0/(1.0-vel.x))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[18])+2*(rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]));
#endif
	  rhs_p[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho); 
	}
		
#else
#ifdef LB_FLUID_BC_X_M_SLIP

	/* SLIP */
          rhs_p[IDX(i-1,j,k)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
         if (jj>= 0 && jj< LNY+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(i-1,jj,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif
#endif
#endif 
 }

      }/* for pp */
    }/* for j,k */
#endif


/************************************/
	/* Y direction */ 
#ifdef LB_FLUID_BC_Y

   for (i = 0; i < LNX + TWO_BRD; i++) 			
   for (k = 0; k < LNZ + TWO_BRD; k++){

         for(pp=0;pp<NPOP;pp++){

	 ii = i+(int)c[pp].x;
	 kk = k+(int)c[pp].z;
	
if(LNY_END == NY){

               j = LNY+BRD-1;	 

#ifdef LB_FLUID_BC_Y_P_OUTLET
	if(pp==0){

      	  if (i >= BRD && i< LNX+BRD && k >= BRD && k < LNZ+BRD) {rho = dens[IDX(i, j, k)];} else {rho=1.0;}

#ifdef GRID_POP_D2Q9
	   vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[8]+2*(rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[5]))/rho;
	   vel.y = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[6]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[7]))/rho;
#endif

#ifdef GRID_POP_D3Q19
	 vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[18]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[13]))/rho;
	 vel.y = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[7]))/rho;	 
	 vel.z = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[11]))/rho;
#endif

          rho = (rho+rho*sqrt(3.0)*(vel.y-u[IDX(i,j,k)].y)+1.0*1.0*1.0)/(1.0+1.0*1.0);

	  rhs_p[IDX(i,j+1,k)] = equilibrium_given_velocity(vel,rho); 
}
	   	   		   		   		   
#else

#ifdef LB_FLUID_BC_Y_P_SLIP
	/* SLIP */
          rhs_p[IDX(i,j+1,k)].p[inv_y[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#elif defined(LB_FLUID_BC_Y_P_VELOCITY)
	  // IMPOSED VELOCITY AT WALL SLIP //
	  fac=1.0;
 #ifdef LB_TEMPERATURE_MELTING
	  /* no moving wall if we are in the solid */
	  if(liquid_frac[IDX(i, j, k)] < 1.0) fac = 0.0;
 #endif
	  if (ii >= 0 && ii < LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j+1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp] - fac*6.0*wgt[pp]*( c[pp].x*property.yp_wall_velocity_x  + c[pp].z*property.yp_wall_velocity_z );
#else
	  /* NO SLIP */
	  if (ii>= 0 && ii< LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j+1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#endif
#endif
}

if(LNY_START == 0){
  j = BRD; 
  
#ifdef LB_FLUID_BC_Y_M_SLIP
           // SLIP //
          rhs_p[IDX(i,j-1,k)].p[inv_y[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#elif defined(LB_FLUID_BC_Y_M_VELOCITY)
	  // IMPOSED VELOCITY AT WALL SLIP //
	  fac=1.0;
 #ifdef LB_TEMPERATURE_MELTING
	  /* no moving wall if we are in the solid */
	  if(liquid_frac[IDX(i, j, k)] < 1.0) fac = 0.0;
 #endif
	  if (ii >= 0 && ii < LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j-1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp] - fac*6.0*wgt[pp]*( c[pp].x*property.ym_wall_velocity_x  + c[pp].z*property.ym_wall_velocity_z );
#else
	  // NO SLIP //
	  if (ii >= 0 && ii < LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j-1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif
  

#ifdef LB_FLUID_BC_Y_M_JET
  // inlet geometry - circle: POSITION OF CENTER COORDINATES OF THE CIRCLE SPECIFIED HERE ITSELF //
	  if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/3.158), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 2.0){    	  
    if(pp==0){         

      /*Poiseuille flow velocity profile - inlet*/
      //   vel.y = -(4.0*(my_double)property.Amp_y*(pow(21.85,-2.0)))*((my_double)center_V[IDX(i,j,k)].x-(((property.SX/2.0)+1.0)-10.925))*(((my_double)center_V[IDX(i,j,k)].x-(((property.SX/2.0)+1.0)-10.925))-21.85);

	  /*constant flow velocity profile - inlet*/
	           vel.y =  0.1;

	           vel.x =  0.0;    //property.Amp_y;
	       	   vel.z =  0.0;    //property.Amp_z;

#ifdef GRID_POP_D2Q9 
	  rho = (1.0/(1.0-vel.y))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[6])+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]));
#endif
#ifdef GRID_POP_D3Q19 
	  rho = (1.0/(1.0-vel.y))*(rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[8])); 
#endif
	  vel.y *=  -1.0; /* we invert velocity and apply a bounce back */

	  //	  rho = (rho + 50.0 * dens[IDX(i,j,k)])/(1.0+50.0);
	  p_eq = equilibrium_given_velocity(vel,rho); 
	}

          rhs_p[IDX(ii,j-1,kk)].p[inv[pp]]   = p_eq.p[pp];

      }
#endif /* end of LB_FLUID_BC_Y_M_JET */

#ifdef LB_FLUID_BC_Y_M_OUTLET
	if(pp==0){

      	  if (i >= BRD && i< LNX+BRD && k >= BRD && k < LNZ+BRD) {rho = dens[IDX(i, j, k)];} else {rho=1.0;}

#ifdef GRID_POP_D2Q9
	  vel.x = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[8]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]))/rho;
	  vel.y = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[6]+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]))/rho;
#endif
	  
#ifdef GRID_POP_D3Q19 
	  vel.x = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[18]+2*(rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]))/rho; 
	  vel.y = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[8]))/rho;
	  vel.z = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[10]+2*(rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[13]))/rho;
#endif

      	  rho = (rho+rho*sqrt(3.0)*(vel.y-u[IDX(i,j,k)].y)+1.0*1.0*1.0)/(1.0+1.0*1.0);

	  rhs_p[IDX(i,j-1,k)] = equilibrium_given_velocity(vel,rho); 
}
#endif


 }
     }/* for pp */
    }/* for i,k */
#endif



/************************************/
	/* Z direction */ 
#ifdef LB_FLUID_BC_Z

   for (i = 0; i < LNX + TWO_BRD; i++) 			
   for (j = 0; j < LNY + TWO_BRD; j++){

      for(pp=0;pp<NPOP;pp++){

	  ii = i+(int)c[pp].x;
	  jj = j+(int)c[pp].y;	

if(LNZ_END == NZ){
	k = LNZ+BRD-1;


#ifdef LB_FLUID_BC_Z_P_OUTLET

	if(pp==0){
      	  if (i >= BRD && i< LNX+BRD && j >= BRD && j < LNY+BRD) {rho = dens[IDX(i, j, k)];} else {rho=1.0;}

#ifdef GRID_POP_D3Q19
	 vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[18]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[13]))/rho;
	 vel.y = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[7]))/rho;	 
	 vel.z = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[11]))/rho;
#endif

      	  rho = (rho+rho*sqrt(3.0)*(vel.z-u[IDX(i,j,k)].z)+1.0*1.0*1.0)/(1.0+1.0*1.0);

	  rhs_p[IDX(i,j,k+1)] = equilibrium_given_velocity(vel,rho);
}

#else
#ifdef LB_FLUID_BC_Z_P_SLIP
	  /* SLIP */
          rhs_p[IDX(i,j,k+1)].p[inv_z[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else	 
	  /* NO SLIP */
	  if (ii>= 0 && ii< LNX+TWO_BRD && jj >= 0 && jj < LNY+TWO_BRD) rhs_p[IDX(ii,jj,k+1)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif
#endif	
}

if(LNZ_START == 0){
  k = BRD; 

#ifdef LB_FLUID_BC_Z_M_OUTLET

	if(pp==0){
      	  if (i >= BRD && i< LNX+BRD && j >= BRD && j < LNY+BRD) {rho = dens[IDX(i, j, k)];} else {rho=1.0;}

#ifdef GRID_POP_D3Q19 
	  vel.x = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[18]+2*(rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]))/rho; 
	  vel.y = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[14]+2*(rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[8]))/rho;
	  vel.z = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[10]+2*(rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[13]))/rho;
#endif

      	  rho = (rho+rho*sqrt(3.0)*(vel.z-u[IDX(i,j,k)].z)+1.0*1.0*1.0)/(1.0+1.0*1.0);

	  rhs_p[IDX(i,j,k-1)] = equilibrium_given_velocity(vel,rho);
}

#else 
#ifdef LB_FLUID_BC_Z_M_SLIP
	  /* SLIP */
          rhs_p[IDX(i,j,k-1)].p[inv_z[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
          /* NO SLIP */
	  if (ii>= 0 && ii< LNX+TWO_BRD && jj >= 0 && jj < LNY+TWO_BRD) rhs_p[IDX(ii,jj,k-1)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif
#endif
 }

      }/* for pp */
    }/* for j,k */
#endif

/************************************/


#ifdef LB_TEMPERATURE	
  pop g_eq, g_eq_w;
  my_double effDT, rho2;
  my_double T_wall;
  //  my_double fac;
  

  /************************/
	/* Y direction */
#ifdef LB_TEMPERATURE_BC_Y


   for (i = BRD; i < LNX + BRD; i++)                     
    for (k = BRD; k < LNZ + BRD; k++){
  //for (i = 0; i < LNX + TWO_BRD; i++)                     
  // for (k = 0; k < LNZ + TWO_BRD; k++){

if(LNY_END == NY){

          j = LNY+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#ifdef LB_TEMPERATURE_BC_Y_P_VARIABLE
	  //	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_top; else T_wall = 0.0;  
	  //T_wall = property.T_top;
	  /* half fixed temperature, half  fixed-flux */
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_top; else T_wall = 0.5*property.grad_T_top + t[IDX(i,j,k)] + property.T_ref;
#endif
#else
	  T_wall = 0.0;
#endif

#ifdef LB_TEMPERATURE_BC_Y_P_OUTLET
	  /* this is the wall temperature in case of outlet */
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i,j-1,k)]);
#endif


#ifdef LB_TEMPERATURE_BC_Y_P_FLUX
	  /* we fix the temperature gradient at the wall */
	  //my_double gradT_top = - property.deltaT/property.SY;
	  T_wall = 0.5*property.grad_T_top + t[IDX(i,j,k)] + property.T_ref;
#endif

	  /* this is a linear interpolation */
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_Y_P_NOFLUX
	  fac = 0.0;
#endif

	      for(pp=0;pp<NPOP;pp++){ 
      	         if(c[pp].y>0){
	      ii = i+(int)c[pp].x;
	      kk = k+(int)c[pp].z;
	              #ifdef LB_TEMPERATURE_BC_Y_P_OUTLET
	     	         if (ii==LNX+TWO_BRD-1 || kk==LNZ+TWO_BRD-1 || ii==0 || kk==0) fac=0.0;
	     	      #endif
	      rhs_g[IDX(ii,j+1,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	    }	   
	  }
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;

#ifdef LB_TEMPERATURE_BC_Y_M_VARIABLE
	  //  if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/2.0), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 2.0)  T_wall = property.T_bot; else T_wall = property.T_top;
	  //if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/3.158), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 2.0) T_wall = property.T_bot; else T_wall = property.T_top; 
	  //if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = 0.0;  
	  /* half fixed temperature, half  fixed-flux */
	  //if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = -0.5*property.grad_T_bot + t[IDX(i,j,k)] + property.T_ref;
#endif

#else
	  T_wall = 0.0;
#endif

#ifdef LB_TEMPERATURE_BC_Y_M_OUTLET
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i,j+1,k)]);
#endif

#ifdef LB_TEMPERATURE_BC_Y_M_FLUX
	  /* we fix the temperature gradient at the wall */
	  //my_double gradT_bot = - property.deltaT/property.SY;
	  T_wall = -0.5*property.grad_T_bot + t[IDX(i,j,k)] + property.T_ref;
#endif
	  fac = 2.0*((T_wall-property.T_ref)-t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_Y_M_NOFLUX
	  fac = 0.0;
#endif

	      for(pp=0;pp<NPOP;pp++){
	        if(c[pp].y<0){
	      ii = i+(int)c[pp].x;
	      kk = k+(int)c[pp].z;
                   #ifdef LB_TEMPERATURE_BC_Y_M_OUTLET
	              if (ii==LNX+TWO_BRD-1 || kk==LNZ+TWO_BRD-1 || ii==0 || kk==0) fac=0.0;
                   #endif
	      rhs_g[IDX(ii,j-1,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;	
	  }
	  }
 }
      
 }
#endif


 /************************/
	/* X direction */
#ifdef LB_TEMPERATURE_BC_X

    for (j = BRD; j < LNY + BRD; j++) 			
     for (k = BRD; k < LNZ + BRD; k++){
       // for (j = 0; j < LNY + TWO_BRD; j++) 			
       //    for (k = 0; k < LNZ + TWO_BRD; k++){


if(LNX_END == NX){

          i = LNX+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#else
	  T_wall = 0.0;
#endif

#ifdef LB_TEMPERATURE_BC_X_P_OUTLET
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i-1,j,k)]);
#endif
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_X_P_NOFLUX
	  fac = 0.0;
#endif	  

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].x>0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;
	          #ifdef LB_TEMPERATURE_BC_X_P_OUTLET
		    if (jj==LNY+TWO_BRD-1 || kk==LNZ+TWO_BRD-1 || jj==0 || kk==0) fac=0.0;
	 	  #endif
	//rhs_g[IDX(i+1,jj,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	  rhs_g[IDX(i+1,jj,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[pp] + wgt[inv[pp]]*fac;
	    }
	  }
 }

if(LNX_START == 0){

	  i = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif

#ifdef LB_TEMPERATURE_BC_X_M_OUTLET
	  /* this is outlet */
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i+1,j,k)]);
#endif
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_X_M_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].x<0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;
	          #ifdef LB_TEMPERATURE_BC_X_M_OUTLET	
	 	    if (jj==LNY+TWO_BRD-1 || kk==LNZ+TWO_BRD-1 || jj==0 || kk==0) fac=0.0;
	          #endif
	// rhs_g[IDX(i-1,jj,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	 rhs_g[IDX(i-1,jj,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[pp] + wgt[inv[pp]]*fac;
	    }
	  }
 }
      
}
#endif

 /************************/
	/* Z direction */
#ifdef LB_TEMPERATURE_BC_Z


    for (i = BRD; i < LNX + BRD; i++) 			
     for (j = BRD; j < LNY + BRD; j++){
       // for (i = 0; i < LNX + TWO_BRD; i++) 			
       //for (j = 0; j < LNY + TWO_BRD; j++){


if(LNZ_END == NZ){

          k = LNZ+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#else
	  T_wall = 0.0;
#endif


#ifdef LB_TEMPERATURE_BC_Z_P_OUTLET
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i,j,k-1)]);
#endif
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_Z_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].z>0){
	  ii = i+(int)c[pp].x;	
	  jj = j+(int)c[pp].y;
               #ifdef LB_TEMPERATURE_BC_Z_P_OUTLET
         	    if (ii==LNX+TWO_BRD-1 || jj==LNY+TWO_BRD-1 || ii==0 || jj==0) fac=0.0;
               #endif
	  rhs_g[IDX(ii,jj,k+1)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	    }
	  }
 }

if(LNZ_START == 0){

	  k = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif

#ifdef LB_TEMPERATURE_BC_Z_M_OUTLET
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i,j,k+1)]);
#endif
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_Z_M_NOFLUX
	  fac = 0.0;
#endif
	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].z<0){
	  ii = i+(int)c[pp].x;
	  jj = j+(int)c[pp].y;
               #ifdef LB_TEMPERATURE_BC_Z_M_OUTLET
	          if (ii==LNX+TWO_BRD-1 || jj==LNY+TWO_BRD-1 || ii==0 || jj==0) fac=0.0;
               #endif
	  rhs_g[IDX(ii,jj,k-1)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	    }
	  }
 }
      
}
#endif

#endif

  /************************************************************************************/
  /* from here the scalar part */

#ifdef LB_SCALAR	
  pop h_eq, h_eq_w;
  my_double effDS;
  my_double S_wall;

  
  /************************/
	/* Y direction */
#ifdef LB_SCALAR_BC_Y

  for (i = BRD; i < LNX + BRD; i++)                     
  for (k = BRD; k < LNZ + BRD; k++){

if(LNY_END == NY){

          j = LNY+BRD-1; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_top;
#else
	  S_wall = 0.0;
#endif

#ifdef LB_SCALAR_BC_Y_P_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i,j-1,k)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_Y_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
      	    if(c[pp].y>0){
	  ii = i+(int)c[pp].x;
	  kk = k+(int)c[pp].z;
               #ifdef LB_SCALAR_BC_Y_P_OUTLET	
         	    if (ii==LNX+TWO_BRD-1 || kk==LNZ+TWO_BRD-1 || ii==0 || jj==0) fac=0.0;
               #endif	  
	  rhs_h[IDX(ii,j+1,kk)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	    }	   
	  }
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_bot;

#ifdef LB_SCALAR_BC_Y_M_VARIABLE
  if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/2.0), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 2.0)  S_wall = property.S_bot; else S_wall = property.S_top;
#endif

#else
	  S_wall = 0.0;
#endif

#ifdef LB_SCALAR_BC_Y_M_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i,j+1,k)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]); 

#ifdef LB_SCALAR_BC_Y_M_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].y<0){
	  ii = i+(int)c[pp].x;
	  kk = k+(int)c[pp].z;
               #ifdef LB_SCALAR_BC_Y_M_OUTLET
	          if (ii==LNX+TWO_BRD-1 || kk==LNZ+TWO_BRD-1 || ii==0 || kk==0) fac=0.0;
               #endif
	  rhs_h[IDX(ii,j-1,kk)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;	
	  }
	  }
 }
      
 }
#endif

 /************************/
	/* X direction */
#ifdef LB_SCALAR_BC_X

  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNX_END == NX){

          i = LNX+BRD-1; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_top;
#else
	  S_wall = 0.0;
#endif


#ifdef LB_SCALAR_BC_X_P_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i-1,j,k)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_X_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].x>0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;
               #ifdef LB_SCALAR_BC_X_P_OUTLET
         	    if (jj==LNY+TWO_BRD-1 || kk==LNZ+TWO_BRD-1 || jj==0 || kk==0) fac=0.0;
               #endif	
	  rhs_h[IDX(i+1,jj,kk)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	    }
	  }
 }

if(LNX_START == 0){

	  i = BRD; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_bot;
#else
	  S_wall = 0.0;
#endif

#ifdef LB_SCALAR_BC_X_M_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i+1,j,k)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_X_M_NOFLUX
	  fac = 0.0;
#endif
	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].x<0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;
               #ifdef LB_SCALAR_BC_X_M_OUTLET
	          if (jj==LNY+TWO_BRD-1 || kk==LNZ+TWO_BRD-1 || jj==0 || kk==0) fac=0.0;
               #endif
	  rhs_h[IDX(i-1,jj,kk)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	    }
	  }
 }
      
}
#endif

 /************************/
	/* Z direction */
#ifdef LB_SCALAR_BC_Z

  for (i = BRD; i < LNX + BRD; i++) 			
    for (j = BRD; j < LNY + BRD; j++){


if(LNZ_END == NZ){

          k = LNZ+BRD-1; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_top;
#else
	  S_wall = 0.0;
#endif

#ifdef LB_SCALAR_BC_Z_P_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i,j,k-1)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_Z_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].z>0){
	  ii = i+(int)c[pp].x;
	  jj = j+(int)c[pp].y;
               #ifdef LB_SCALAR_BC_Z_P_OUTLET
	          if (ii==LNX+TWO_BRD-1 || jj==LNY+TWO_BRD-1 || ii==0 || jj==0) fac=0.0;
               #endif	
	  rhs_h[IDX(ii,jj,k+1)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	    }
	  }
 }

if(LNZ_START == 0){

	  k = BRD; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_bot;
#else
	  S_wall = 0.0;
#endif

#ifdef LB_SCALAR_BC_Z_M_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i,j,k+1)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_Z_M_NOFLUX
	  fac = 0.0;
#endif
	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].z<0){
	  ii = i+(int)c[pp].x;
	  jj = j+(int)c[pp].y;
               #ifdef LB_SCALAR_BC_Z_M_OUTLET
        	    if (ii==LNX+TWO_BRD-1 || jj==LNY+TWO_BRD-1 || ii==0 || jj==0) fac=0.0;
               #endif
	  rhs_h[IDX(ii,jj,k-1)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[inv[pp]]*fac;
	    }
	  }
 }
      
}
#endif

#endif /* end of scalar part */

}/* end of bc for streaming */
#endif
/**********************************************************************************************************/
/* When STREAMING is active bc shall be implemented in a different way*/
//#ifdef LB_FLUID_BC
#ifdef METHOD_STREAMING
void boundary_and_pbc_conditions_for_streaming_old(){

  int i,j,k,pp;
  vector vel;
  my_double rho;
  pop p_eq,p_neq;
  int ii,jj,kk;	
  my_double fac;

  /* communications to be done in any case  (especially for pbc)*/

#ifdef LB_FLUID
sendrecv_borders_pop(rhs_p);
#endif

#ifdef LB_TEMPERATURE
sendrecv_borders_pop(rhs_g);
//sendrecv_borders_scalar(t);
#endif

#ifdef LB_SCALAR
sendrecv_borders_pop(rhs_h);
#endif

/***************** Bounce-back BC for the complex boundaries  ********************/
  int iii, jjj, kkk;
#ifdef LB_FLUID_FORCING_LANDSCAPE
 for(k=BRD;k<LNZ+BRD;k++){
   for(j=BRD;j<LNY+BRD;j++){
      for(i=BRD;i<LNX+BRD;i++){ 

     for(pp=0;pp<NPOP;pp++){
          iii = i+(int)c[pp].x;
	  jjj = j+(int)c[pp].y;
	  kkk = k+(int)c[pp].z;

	  if(landscape[IDX(i, j, k)] == 0.0 && landscape[IDX(iii, jjj, kkk)] == 1.0 ) rhs_p[IDX(iii,jjj,kkk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
     }/* pp */   
      }/* i */
   }/* j */
 }/* k */
#endif

/************************************/
	/* X direction */ 
#ifdef LB_FLUID_BC_X

   for (j = 0; j < LNY + TWO_BRD; j++) 			
   for (k = 0; k < LNZ + TWO_BRD; k++){

      for(pp=0;pp<NPOP;pp++){

	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;	

if(LNX_END == NX){
	i = LNX+BRD-1;

#ifdef LB_FLUID_BC_X_P_OUTLET
	if(pp==0){
	  rho = 1;
#ifdef GRID_POP_D2Q9
	  vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[4]+2*(rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[5]))/rho;
#endif

#ifdef GRID_POP_D3Q19
	 vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16]+2*(rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[1]))/rho; 
#endif

	 vel.y = 0.0; //u[IDX(i, j, k)].y; 
	 vel.z = 0.0; //u[IDX(i, j, k)].z;
	 rhs_p[IDX(i+1,j,k)] = equilibrium_given_velocity(vel,rho); 
	} 
#else

#ifdef LB_FLUID_BC_X_P_SLIP
	  /* SLIP */
          rhs_p[IDX(i+1,j,k)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#else	 
	  /* NO SLIP */
	  if (jj>= 0 && jj< LNY+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(i+1,jj,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#endif
#endif
}

if(LNX_START == 0){
  i = BRD;

#ifdef LB_FLUID_BC_X_M_OUTLET
	if(pp==0){
	  rho = 1;
#ifdef GRID_POP_D2Q9
	  vel.x = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[4]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]))/rho;
#endif
	  
#ifdef GRID_POP_D3Q19 
	  vel.x = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16]+2*(rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[2]))/rho; 
#endif
	  vel.y = 0.0; //u[IDX(i, j, k)].y;
	  vel.z = 0.0; //u[IDX(i, j, k)].z; 
	  rhs_p[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho); 
	}
#else 
#ifdef LB_FLUID_BC_X_M_INLET
	if(pp==0){

#ifdef LB_FLUID_BC_X_M_INLET_POISEUILLE
	  /*Poiseuille flow velocity profile - inlet*/
	     vel.x = -(4.0*(my_double)property.Amp_x*(pow(((my_double)property.SY),-2.0)))*((my_double)center_V[IDX(i,j,k)].y)*(((my_double)center_V[IDX(i,j,k)].y)-((my_double)property.SY));
#endif

#ifdef LB_FLUID_BC_X_M_INLET_POISEUILLE_HALF
	  /* Half Poiseuille flow velocity profile - inlet*/
	  vel.x = -(4.0*(my_double)property.Amp_x*(pow((2.0*(my_double)property.SY),-2.0)))*((my_double)center_V[IDX(i,j,k)].y)*(((my_double)center_V[IDX(i,j,k)].y)-(2.0*(my_double)property.SY));
#endif

#ifdef LB_FLUID_BC_X_M_INLET_CONSTANT
	  /* Constant velocity = Amp_x*/        
	  vel.x = property.Amp_x;
	  vel.y = property.Amp_y;
	  vel.z = property.Amp_z;
#endif

#ifdef GRID_POP_D2Q9 
	  rho = (1.0/(1.0-vel.x))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[4])+2*(rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[3]));
#endif

#ifdef GRID_POP_D3Q19
	  rho = (1.0/(1.0-vel.x))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16])+2*(rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[2]));
#endif

	  rhs_p[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho); 
	}		
#else
#ifdef LB_FLUID_BC_X_M_SLIP
	/* SLIP */
          rhs_p[IDX(i-1,j,k)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
         if (jj>= 0 && jj< LNY+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(i-1,jj,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif
#endif
#endif 
 }

      }/* for pp */
    }/* for j,k */
#endif

/************************************/
	/* Y direction */ 
#ifdef LB_FLUID_BC_Y

   for (i = 0; i < LNX + TWO_BRD; i++) 			
   for (k = 0; k < LNZ + TWO_BRD; k++){

         for(pp=0;pp<NPOP;pp++){

	 ii = i+(int)c[pp].x;
	 kk = k+(int)c[pp].z;
	
if(LNY_END == NY){
          j = LNY+BRD-1;	 

#ifdef LB_FLUID_BC_Y_P_OUTLET
	if(pp==0){

	  rho = 1.0; //dens[IDX(i, j, k)];

#ifdef GRID_POP_D2Q9 
	  vel.y = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[6]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[7]))/rho;
#endif

#ifdef GRID_POP_D3Q19
	  vel.y = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[11]+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[7]))/rho; 
#endif
 
	  vel.x = 0.0; //u[IDX(i, j, k)].x; 
	  vel.z = 0.0; //u[IDX(i, j, k)].z; 
	  rhs_p[IDX(i,j+1,k)] = equilibrium_given_velocity(vel,rho); 
}		   	   		   		   		   
#else

#ifdef LB_FLUID_BC_Y_P_SLIP
	/* SLIP */
          rhs_p[IDX(i,j+1,k)].p[inv_y[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	  /* NO SLIP */
	  if (ii>= 0 && ii< LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j+1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#endif
#endif
}

if(LNY_START == 0){
  j = BRD; 

#ifdef LB_FLUID_BC_Y_M_SLIP
        /* SLIP */
          rhs_p[IDX(i,j-1,k)].p[inv_y[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
	  if (ii >= 0 && ii < LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j-1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif


#ifdef LB_FLUID_BC_Y_M_JET
  /* inlet geometry - circle: POSITION OF CENTER COORDINATES OF THE CIRCLE SPECIFIED HERE ITSELF */
  if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/2.0), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 7.0){    
	  
    if(pp==0){         
          vel.x =  0.0;    //property.Amp_y;
	  vel.y =  0.1;    //property.Amp_x;
	  vel.z =  0.0;    //property.Amp_z;
#ifdef GRID_POP_D2Q9 
	  rho = (1.0/(1.0-vel.y))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[2])+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]));
#endif
#ifdef GRID_POP_D3Q19 
	  rho = (1.0/(1.0-vel.y))*(rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[11]+2*(rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[8])); 
#endif
	  vel.y *=  -1.0; /* we invert velocity and apply a bounce back */
	  p_eq = equilibrium_given_velocity(vel,rho); 
	}
          rhs_p[IDX(ii,j-1,kk)].p[inv[pp]]   = p_eq.p[pp];

 	
      }
#endif /* end of LB_FLUID_BC_Y_M_JET */

 }

     }/* for pp */
    }/* for i,k */
#endif

/************************************/
	/* Z direction */ 
#ifdef LB_FLUID_BC_Z

   for (i = 0; i < LNX + TWO_BRD; i++) 			
   for (j = 0; j < LNY + TWO_BRD; j++){

      for(pp=0;pp<NPOP;pp++){

	  ii = i+(int)c[pp].x;
	  jj = j+(int)c[pp].y;	

if(LNZ_END == NZ){
	k = LNZ+BRD-1;


#ifdef LB_FLUID_BC_Z_P_OUTLET

	if(pp==0){
	  rho = 1.0;

#ifdef GRID_POP_D3Q19
	 vel.z = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[7]+2*(rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[11]))/rho; 
#endif
	 vel.x = 0.0; //u[IDX(i, j, k)].x; //0.0;
	 vel.y = 0.0; //u[IDX(i, j, k)].y; //0.0;
	  rhs_p[IDX(i,j,k+1)] = equilibrium_given_velocity(vel,rho);
}

#else
#ifdef LB_FLUID_BC_Z_P_SLIP
	//          rhs_p[IDX(i+1,jj,kk)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];
 
          rhs_p[IDX(i,j,k+1)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#else	 
	  /* NO SLIP */
	  //rhs_p[IDX(i+1,j,k)].p[pp] = rhs_p[IDX(i,jj,kk)].p[inv[pp]];	
	  if (ii>= 0 && ii< LNX+TWO_BRD && jj >= 0 && jj < LNY+TWO_BRD) rhs_p[IDX(ii,jj,k+1)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#endif
#endif	
}

if(LNZ_START == 0){
  k = BRD; 

#ifdef LB_FLUID_BC_Z_M_OUTLET

	if(pp==0){
	  rho = 1.0;

#ifdef GRID_POP_D3Q19
	 vel.z = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[7]+2*(rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[13]))/rho; 
#endif

	  vel.x = 0.0; //u[IDX(i, j, k)].x; //0.0;
	  vel.y = 0.0; //u[IDX(i, j, k)].y; //0.0;
	  rhs_p[IDX(i,j,k-1)] = equilibrium_given_velocity(vel,rho);
}

#else 
#ifdef LB_FLUID_BC_Z_M_SLIP

          rhs_p[IDX(i,j,k-1)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
         //rhs_p[IDX(i-1,j,k)].p[pp] = rhs_p[IDX(i,jj,kk)].p[inv[pp]];
	  if (ii>= 0 && ii< LNX+TWO_BRD && jj >= 0 && jj < LNY+TWO_BRD) rhs_p[IDX(ii,jj,k-1)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif
#endif
 }

      }/* for pp */
    }/* for j,k */
#endif

/************************************/



#ifdef LB_TEMPERATURE	
  pop g_eq, g_eq_w;
  my_double effDT, rho2;
  my_double T_wall;
  //  my_double fac;
  

  /************************/
	/* Y direction */
#ifdef LB_TEMPERATURE_BC_Y


  for (i = BRD; i < LNX + BRD; i++)                     
  for (k = BRD; k < LNZ + BRD; k++){
  //for (i = 0; i < LNX + TWO_BRD; i++) 			
  //for (k = 0; k < LNZ + TWO_BRD; k++){


if(LNY_END == NY){

          j = LNY+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#ifdef LB_TEMPERATURE_BC_Y_P_VARIABLE
	  //	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_top; else T_wall = 0.0;  
	  T_wall = property.T_top;
#endif
#else
	  T_wall = 0.0;
#endif



#ifdef LB_TEMPERATURE_BC_Y_P_OUTLET
	  /* this is the wall temperature in case of outlet */
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i,j-1,k)]);
#endif
	  /* this is a linear interpolation */
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_Y_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
      	    if(c[pp].y>0){
	      ii = i+(int)c[pp].x;
	      kk = k+(int)c[pp].z;		  
	      rhs_g[IDX(ii,j+1,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }	   
	  }
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#ifdef LB_TEMPERATURE_BC_Y_M_VARIABLE
  if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/2.0), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 2.0)  T_wall = property.T_bot; else T_wall = property.T_top;
  //if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = 0.0;  
#endif
#else
	  T_wall = 0.0;
#endif



#ifdef LB_TEMPERATURE_BC_Y_M_OUTLET
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i,j+1,k)]);
#endif

	  fac = 2.0*((T_wall-property.T_ref)-t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_Y_M_NOFLUX
	  fac = 0.0;
#endif
	  

	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].y<0){
	      ii = i+(int)c[pp].x;
	      kk = k+(int)c[pp].z;
	      rhs_g[IDX(ii,j-1,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;	
	  }
	  }
 }
      
 }
#endif



 /************************/
	/* X direction */
#ifdef LB_TEMPERATURE_BC_X

  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNX_END == NX){

          i = LNX+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#else
	  T_wall = 0.0;
#endif


#ifdef LB_TEMPERATURE_BC_X_P_OUTLET
	  /* this is outlet */
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i-1,j,k)]);
#endif
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_X_P_NOFLUX
	  fac = 0.0;
#endif	  


	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].x>0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;	
	  rhs_g[IDX(i+1,jj,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }
	  }
 }

if(LNX_START == 0){

	  i = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif

#ifdef LB_TEMPERATURE_BC_X_M_OUTLET
	  /* this is outlet */
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i+1,j,k)]);
#endif
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_X_M_NOFLUX
	  fac = 0.0;
#endif

	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].x<0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;	
	  rhs_g[IDX(i-1,jj,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }
	  }
 }
      
}
#endif

 /************************/
	/* Z direction */
#ifdef LB_TEMPERATURE_BC_Z

  for (i = BRD; i < LNX + BRD; i++) 			
    for (j = BRD; j < LNY + BRD; j++){


if(LNZ_END == NZ){

          k = LNZ+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#else
	  T_wall = 0.0;
#endif


#ifdef LB_TEMPERATURE_BC_Z_P_OUTLET
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i,j,k-1)]);
#endif
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_Z_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].z>0){
	  ii = i+(int)c[pp].x;	
	  jj = j+(int)c[pp].y;
	  rhs_g[IDX(ii,jj,k+1)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }
	  }
 }

if(LNZ_START == 0){

	  k = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif


#ifdef LB_TEMPERATURE_BC_Z_M_OUTLET
	  T_wall =  t[IDX(i,j,k)] + 0.5*( t[IDX(i,j,k)] - t[IDX(i,j,k+1)]);
#endif
	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

#ifdef LB_TEMPERATURE_BC_Z_M_NOFLUX
	  fac = 0.0;
#endif
	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].z<0){
	  ii = i+(int)c[pp].x;
	  jj = j+(int)c[pp].y;
	  rhs_g[IDX(ii,jj,k-1)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }
	  }
 }
      
}
#endif

#endif


  /* from here the scalar part */

#ifdef LB_SCALAR	
  pop h_eq, h_eq_w;
  my_double effDS;
  my_double S_wall;

  
  /************************/
	/* Y direction */
#ifdef LB_SCALAR_BC_Y


  for (i = BRD; i < LNX + BRD; i++)                     
  for (k = BRD; k < LNZ + BRD; k++){


if(LNY_END == NY){

          j = LNY+BRD-1; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_top;
#else
	  S_wall = 0.0;
#endif


#ifdef LB_SCALAR_BC_Y_P_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i,j-1,k)]);
#endif

	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_Y_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
      	    if(c[pp].y>0){
	      ii = i+(int)c[pp].x;
	      kk = k+(int)c[pp].z;	
	  	rhs_h[IDX(ii,j+1,kk)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }	   
	  }
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_bot;
#else
	  S_wall = 0.0;
#endif

#ifdef LB_SCALAR_BC_Y_M_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i,j+1,k)]);
#endif

	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]); 

#ifdef LB_SCALAR_BC_Y_M_NOFLUX
	  fac = 0.0;
#endif


	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].y<0){
	      ii = i+(int)c[pp].x;
	      kk = k+(int)c[pp].z;
	      rhs_h[IDX(ii,j-1,kk)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;	
	  }
	  }
 }
      
 }
#endif

 /************************/
	/* X direction */
#ifdef LB_SCALAR_BC_X

  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNX_END == NX){

          i = LNX+BRD-1; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_top;
#else
	  S_wall = 0.0;
#endif


#ifdef LB_SCALAR_BC_X_P_OUTLET
	  /* this is outlet */
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i-1,j,k)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_X_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].x>0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;	
	  rhs_h[IDX(i+1,jj,kk)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }
	  }
 }

if(LNX_START == 0){

	  i = BRD; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_bot;
#else
	  S_wall = 0.0;
#endif


#ifdef LB_SCALAR_BC_X_M_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i+1,j,k)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_X_M_NOFLUX
	  fac = 0.0;
#endif
	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].x<0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;	
	  rhs_h[IDX(i-1,jj,kk)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }
	  }
 }
      
}
#endif

 /************************/
	/* Z direction */
#ifdef LB_SCALAR_BC_Z

  for (i = BRD; i < LNX + BRD; i++) 			
    for (j = BRD; j < LNY + BRD; j++){


if(LNZ_END == NZ){

          k = LNZ+BRD-1; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_top;
#else
	  S_wall = 0.0;
#endif


#ifdef LB_SCALAR_BC_Z_P_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i,j,k-1)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_Z_P_NOFLUX
	  fac = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].z>0){
	  ii = i+(int)c[pp].x;
	  jj = j+(int)c[pp].y;	
	  rhs_h[IDX(ii,jj,k+1)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }
	  }
 }

if(LNZ_START == 0){

	  k = BRD; 

#ifndef LB_SCALAR_FLUCTUATION 
	  S_wall = property.S_bot;
#else
	  S_wall = 0.0;
#endif


#ifdef LB_SCALAR_BC_Z_M_OUTLET
	  S_wall =  s[IDX(i,j,k)] + 0.5*( s[IDX(i,j,k)] - s[IDX(i,j,k+1)]);
#endif
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);

#ifdef LB_SCALAR_BC_Z_M_NOFLUX
	  fac = 0.0;
#endif
	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].z<0){
	  ii = i+(int)c[pp].x;
	  jj = j+(int)c[pp].y;	
	  rhs_h[IDX(ii,jj,k-1)].p[inv[pp]] =  rhs_h[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }
	  }
 }
      
}
#endif

#endif /* end of scalar part */

}/* end of bc for streaming */
#endif
//#endif



/***********************************************/

#ifdef LB_FLUID_BC
void boundary_conditions_for_equilibrium(char which_pop){

  int i,j,k,pp;
  vector vel;
  my_double rho;
  pop f_eq;

  if(which_pop == 'p'){
  /* y direction  velocity */
#ifdef LB_FLUID_BC_Y
 
  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){
      for(pp=0;pp<NPOP;pp++){

      /*  top  */
	if(LNY_END == NY){
	   j = LNY+BRD-1;

#ifdef LB_FLUID_BC_Y_P_SLIP
	   /* free slip */
	p_eq[IDX(i,j+1,k)].p[pp] = p_eq[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) p_eq[IDX(i,j+1,k)].p[pp] = p_eq[IDX(i,j,k)].p[pp];
       	p_eq[IDX(i,j+2,k)].p[pp] = p_eq[IDX(i,j-1,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) p_eq[IDX(i,j+2,k)].p[pp] = p_eq[IDX(i,j-1,k)].p[pp];
#else
	/* no slip is the default */
	   p_eq[IDX(i,j+1,k)].p[pp] = p_eq[IDX(i,j,k)].p[inv[pp]];
	   p_eq[IDX(i,j+2,k)].p[pp] = p_eq[IDX(i,j-1,k)].p[inv[pp]];
#endif

	}

     /* bottom  */
	if(LNY_START == 0){
           j = BRD; 

#ifdef LB_FLUID_BC_Y_M_SLIP
	   /* free slip */
	p_eq[IDX(i,j-1,k)].p[pp] = p_eq[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) p_eq[IDX(i,j-1,k)].p[pp] = p_eq[IDX(i,j,k)].p[pp];
       	p_eq[IDX(i,j-2,k)].p[pp] = p_eq[IDX(i,j+1,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) p_eq[IDX(i,j-2,k)].p[pp] = p_eq[IDX(i,j+1,k)].p[pp];
#else
	/* no slip is the default */
             p_eq[IDX(i,j-1,k)].p[pp] =  p_eq[IDX(i,j,k)].p[inv[pp]];
             p_eq[IDX(i,j-2,k)].p[pp] =  p_eq[IDX(i,j+1,k)].p[inv[pp]]; 
#endif	                 

	}

      }
    }
#endif
  }/* end of if which pop */

  if(which_pop == 'g'){
  /* y direction temperature */	
#ifdef LB_TEMPERATURE
#ifdef LB_TEMPERATURE_BC_Y

  my_double effDT, rho2;
  my_double T_wall,fac;

  /**********************************************************/

  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNY_END == NY){

 	  j = LNY+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#else
	  T_wall = 0.0;
#endif

	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g_eq[IDX(i,j+1,k)].p[pp] =  fac*g_eq[IDX(i,j,k)].p[pp];

#ifdef METHOD_MYQUICK
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j-1,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g_eq[IDX(i,j+2,k)].p[pp] =  fac*g_eq[IDX(i,j-1,k)].p[pp];

#endif
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g_eq[IDX(i,j-1,k)].p[pp] =  fac*g_eq[IDX(i,j,k)].p[pp];

     
#ifdef METHOD_MYQUICK 	  
#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif   
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j+1,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g_eq[IDX(i,j-2,k)].p[pp] =  fac*g_eq[IDX(i,j+1,k)].p[pp];
#endif
 }

      
    }
#endif
#endif
  }/* end of if which pop g */
/***************************************************************************************************/
}
#endif
/* end of function boundary_conditions_for_equilibrium */


/****************************************************************************************************/
/* BC for advection (in the FV method )*/
#ifdef LB_FLUID_BC
void boundary_conditions_for_advection(pop * f, char which_pop){

  int i,j,k,pp;
  vector vel;
  my_double rho,fac;
  pop p_eq;
  pop p_neq;

  if(which_pop == 'p'){

	/* X direction */	
#ifdef LB_FLUID_BC_X

  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){
      for(pp=0;pp<NPOP;pp++){

if(LNX_END == NX){
	i = LNX+BRD-1;
	
#ifdef LB_FLUID_BC_X_P_OUTLET
	
	if(pp==0){
	  vel.x = u[IDX(i, j, k)].x;
	  vel.y = u[IDX(i, j, k)].y;
	  vel.z = u[IDX(i, j, k)].z;
	  rho = dens[IDX(i, j, k)];
	  f[IDX(i+1,j,k)] = equilibrium_given_velocity(vel,rho);
#ifdef METHOD_MYQUICK
	  f[IDX(i+2,j,k)] = f[IDX(i+1,j,k)];
#endif	  
	}	  
#else 
#ifdef LB_FLUID_BC_X_P_SLIP

	f[IDX(i+1,j,k)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].z != 0.0 ) p[IDX(i+1,j,k)].p[pp] = p[IDX(i,j,k)].p[pp];
#else
	/* NOSLIP */
	f[IDX(i+1,j,k)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
        f[IDX(i+2,j,k)].p[pp] = f[IDX(i-1,j,k)].p[inv[pp]];
#endif
#endif
#endif

 }/* if */

if(LNX_START == 0){
  i = BRD; 

#ifdef LB_FLUID_BC_X_M_INLET
  /* the following if is to compute the equilibrium only one time */
	if(pp==0){
	  vel.x = property.Amp_x;
	  vel.y = property.Amp_y;
	  vel.z = property.Amp_z;
	    rho = 1.0;
	  f[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho);
#ifdef METHOD_MYQUICK
	  f[IDX(i-2,j,k)] = f[IDX(i-1,j,k)];
#endif
	}
	
#else
#ifdef LB_FLUID_BC_X_M_SLIP
		
	f[IDX(i-1,j,k)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].z != 0.0 ) f[IDX(i-1,j,k)].p[pp] = f[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
	f[IDX(i-1,j,k)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
        f[IDX(i-2,j,k)].p[pp] = f[IDX(i+1,j,k)].p[inv[pp]];
#endif
#endif
#endif

 }/* if */

      }/* for pp */
    }/* for j,k */
#endif


  /************************************/

	/* Y direction */	
#ifdef LB_FLUID_BC_Y
 for (i = 0; i < LNX + TWO_BRD; i++) 			
    for (k = 0; k < LNZ + TWO_BRD; k++){
      for(pp=0;pp<NPOP;pp++){

if(LNY_END == NY){
	j = LNY+BRD-1; 

#ifdef LB_FLUID_BC_Y_P_OUTLET
	if(pp==0){
	  vel.x =  u[IDX(i, j, k)].x;
	  vel.y =  u[IDX(i, j, k)].y;
	  vel.z =  u[IDX(i, j, k)].z;
	  rho = dens[IDX(i, j, k)];
	  f[IDX(i,j+1,k)] = equilibrium_given_velocity(vel,rho);
#ifdef METHOD_MYQUICK
          f[IDX(i,j+2,k)] = f[IDX(i,j+1,k)];
#endif
	}
#else

#ifdef LB_FLUID_BC_Y_P_SLIP
	f[IDX(i,j+1,k)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) f[IDX(i,j+1,k)].p[pp] = f[IDX(i,j,k)].p[pp];
#ifdef METHOD_MYQUICK
       	f[IDX(i,j+2,k)].p[pp] = f[IDX(i,j-1,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) f[IDX(i,j+2,k)].p[pp] = f[IDX(i,j-1,k)].p[pp];
#endif
#else
	f[IDX(i,j+1,k)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
       	f[IDX(i,j+2,k)].p[pp] = f[IDX(i,j-1,k)].p[inv[pp]];
#endif

#endif
#endif
 }

if(LNY_START == 0){
  j = BRD; 
#ifdef LB_FLUID_BC_Y_M_SLIP	
	f[IDX(i,j-1,k)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) f[IDX(i,j-1,k)].p[pp] = f[IDX(i,j,k)].p[pp];
#ifdef METHOD_MYQUICK
       	  f[IDX(i,j-2,k)].p[pp] = f[IDX(i,j+1,k)].p[inv[pp]];
	if(c[pp].x != 0.0 || c[pp].z != 0.0 ) f[IDX(i,j-2,k)].p[pp] = f[IDX(i,j+1,k)].p[pp];
#endif
#else
	  f[IDX(i,j-1,k)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
#ifdef METHOD_MYQUICK
       	  f[IDX(i,j-2,k)].p[pp] = f[IDX(i,j+1,k)].p[inv[pp]];	
#endif  
#endif
 }

      }/* for pp */
    }/* for i,k */
#endif


  /*****************************************************************************************/
  /* Z direction */	

#ifdef LB_FLUID_BC_Z

  for (j = BRD; j < LNY + BRD; j++) 			
    for (i = BRD; i < LNX + BRD; i++){
      for(pp=0;pp<NPOP;pp++){

	if(LNZ_END == NZ){
	k = LNZ+BRD-1;

#ifdef LB_FLUID_BC_Z_P_SLIP
	f[IDX(i,j,k+1)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].x != 0.0 ) f[IDX(i,j,k+1)].p[pp] = f[IDX(i,j,k)].p[pp];
#else
	/* NOSLIP */
	f[IDX(i,j,k+1)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
#endif
 }/* if */

	if(LNZ_START == 0){
	  k = BRD; 

#ifdef LB_FLUID_BC_Z_M_SLIP		
	f[IDX(i,j,k-1)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
	if(c[pp].y != 0.0 || c[pp].x != 0.0 ) f[IDX(i,j,k-1)].p[pp] = f[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
	f[IDX(i,j,k-1)].p[pp] = f[IDX(i,j,k)].p[inv[pp]];
#endif
 }/* if */

      }/* for pp */
    }/* for j,i */
#endif
  }/* end of which_pop on p */


  if(which_pop == 'g'){
#ifdef LB_TEMPERATURE_BC_KEEP_WITHIN

  my_double T_max, T_min;

  /* Rather DIRTY TRICK , in order to keep the temperature inside the boundaries, in the whole bulk */
  T_max = property.T_bot;
  T_min = property.T_top;
  if( property.T_top > property.T_bot){ T_max = property.T_top;  T_min = property.T_bot;}

  for (i = BRD; i < LNX + BRD; i++) 
    for (j = BRD; j < LNY + BRD; j++)
      for (k = BRD; k < LNZ + BRD; k++){

	if(m(f[IDX(i,j,k)]) < T_min) for(pp=0;pp<NPOP;pp++) f[IDX(i,j,k)].p[pp] *= T_min/m(f[IDX(i,j,k)]);
	if(m(f[IDX(i,j,k)]) > T_max) for(pp=0;pp<NPOP;pp++) f[IDX(i,j,k)].p[pp] *= T_max/m(f[IDX(i,j,k)]); 
      }
#endif


	/* Y direction */	
#ifdef LB_TEMPERATURE
#ifdef LB_TEMPERATURE_BC_Y

  pop g_eq, g_eq_w;
  my_double effDT, rho2;
  my_double T_wall;
  /************************/

  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNY_END == NY){

 	  j = LNY+BRD-1; 


#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
	  
 #ifdef LB_TEMPERATURE_BC_Y_P_VARIABLE
	  /* half fixed temperature, half  fixed-flux */
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_top; else T_wall = 0.5*property.grad_T_top + t[IDX(i,j,k)] + property.T_ref;
 #endif

 #ifdef LB_TEMPERATURE_BC_Y_P_FLUX
	  /* we fix the temperature gradient at the wall */
	  T_wall = 0.5*property.grad_T_top + t[IDX(i,j,k)] + property.T_ref;
 #endif
	  
#else
	  T_wall = 0.0;
#endif


	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) f[IDX(i,j+1,k)].p[pp] =  fac*f[IDX(i,j,k)].p[pp];

#ifdef METHOD_MYQUICK
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j-1,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) f[IDX(i,j+2,k)].p[pp] =  fac*f[IDX(i,j-1,k)].p[pp];

#endif
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
 #ifdef LB_TEMPERATURE_BC_Y_M_VARIABLE
	  //if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/2.0), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 2.0)  T_wall = property.T_bot; else T_wall = property.T_top;
  //  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = 0.0;  
	  /* half fixed temperature, half  fixed-flux */
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = -0.5*property.grad_T_bot + t[IDX(i,j,k)] + property.T_ref;
 #endif

 #ifdef LB_TEMPERATURE_BC_Y_M_FLUX
	  /* we fix the temperature gradient at the wall */
	  T_wall = -0.5*property.grad_T_bot + t[IDX(i,j,k)] + property.T_ref;
 #endif

#else
	  T_wall = 0.0;
#endif


	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) f[IDX(i,j-1,k)].p[pp] =  fac*f[IDX(i,j,k)].p[pp];

     
#ifdef METHOD_MYQUICK 	  
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j+1,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) f[IDX(i,j-2,k)].p[pp] =  fac*f[IDX(i,j+1,k)].p[pp];
#endif
}

      
    }
#endif
 
#ifdef LB_TEMPERATURE_BC_X
#ifdef LB_TEMPERATURE_BC_X_NOFLUX
 /* now along X , the default is insultaing BC*/
  for (j = BRD; j < LNY + BRD; j++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNX_END == NX){

 	  i = LNX+BRD-1; 

	  for(pp=0;pp<NPOP;pp++) f[IDX(i+1,j,k)].p[pp] =  f[IDX(i,j,k)].p[pp];
#ifdef METHOD_MYQUICK	  
	  for(pp=0;pp<NPOP;pp++) f[IDX(i+2,j,k)].p[pp] =  f[IDX(i-1,j,k)].p[pp];
#endif
 }

if(LNX_START == 0){

	  i = BRD; 
	
	  for(pp=0;pp<NPOP;pp++) f[IDX(i-1,j,k)].p[pp] =  f[IDX(i,j,k)].p[pp];     
#ifdef METHOD_MYQUICK 	  
	  for(pp=0;pp<NPOP;pp++) f[IDX(i-2,j,k)].p[pp] =  f[IDX(i+1,j,k)].p[pp];
#endif
 }      
    }
#endif
#endif

#endif
  }/* end of which_pop on g */


} /* end of void boundary conditions */
#endif



/****************************************************************/

/* this function implement periodic or wall boundary conditions for the hydrodynamics fields 
   it is needed to correctly interpolate these fields at the particle positions (i.e. when LAGRANGE is active) */
void boundary_conditions_hydro(){

  int i,j,k;
  my_double fac, T_wall, S_wall;
  vector cp,cm; /* just two useful coefficients */
 
  /* bc for the velocity field */
#ifdef LB_FLUID
sendrecv_borders_vector(u);

 #ifdef LB_FLUID_BC
 
  /* X direction */
  #ifdef LB_FLUID_BC_X

   #ifdef LB_FLUID_BC_X_M_SLIP
    cm.x = -1.0; 
    cm.y = cm.z = 1.0;
   #else
    cm.x = cm.y = cm.z = -1.0; 
   #endif

   for (j = 0; j < LNY + TWO_BRD; j++)                     
    for (k = 0; k < LNZ + TWO_BRD; k++){
     if(LNX_START == 0){
            i = BRD; 
            u[IDX(i-1,j,k)].x =  cm.x*u[IDX(i,j,k)].x;
            u[IDX(i-1,j,k)].y =  cm.y*u[IDX(i,j,k)].y;
            u[IDX(i-1,j,k)].z =  cm.z*u[IDX(i,j,k)].z;
        }
    }

   #ifdef LB_FLUID_BC_X_P_SLIP
    cp.x = -1.0; 
    cp.y = cp.z = 1.0;
   #else
    cp.x = cp.y = cp.z = -1.0; 
   #endif
  
  for (j = 0; j < LNY + TWO_BRD; j++)                     
   for (k = 0; k < LNZ + TWO_BRD; k++){
    if(LNX_END == NX){
            i = LNX+BRD-1;
            u[IDX(i+1,j,k)].x =  cp.x*u[IDX(i,j,k)].x;
            u[IDX(i+1,j,k)].y =  cp.y*u[IDX(i,j,k)].y;
            u[IDX(i+1,j,k)].z =  cp.z*u[IDX(i,j,k)].z;
      }
  }
  #endif /* end ifdef LB_FLUID_BC_X */

  /* Y direction */
  #ifdef LB_FLUID_BC_Y

   #ifdef LB_FLUID_BC_Y_M_SLIP
    cm.y = -1.0; 
    cm.x = cm.z = 1.0;
   #else
    cm.x = cm.y = cm.z = -1.0; 
   #endif

   for (i = 0; i < LNX + TWO_BRD; i++)                     
    for (k = 0; k < LNZ + TWO_BRD; k++){
     if(LNY_START == 0){
            j = BRD; 
            u[IDX(i,j-1,k)].x =  cm.x*u[IDX(i,j,k)].x;
            u[IDX(i,j-1,k)].y =  cm.y*u[IDX(i,j,k)].y;
            u[IDX(i,j-1,k)].z =  cm.z*u[IDX(i,j,k)].z;
        }
    }

   #ifdef LB_FLUID_BC_Y_P_SLIP
    cp.y = -1.0; 
    cp.x = cp.z = 1.0;
   #else
    cp.x = cp.y = cp.z = -1.0; 
   #endif
  
  for (i = 0; i < LNX + TWO_BRD; i++)                     
   for (k = 0; k < LNZ + TWO_BRD; k++){
    if(LNY_END == NY){
            j = LNY+BRD-1;
            u[IDX(i,j+1,k)].x =  cp.x*u[IDX(i,j,k)].x;
            u[IDX(i,j+1,k)].y =  cp.y*u[IDX(i,j,k)].y;
            u[IDX(i,j+1,k)].z =  cp.z*u[IDX(i,j,k)].z;
      }
  }
  #endif /* end ifdef LB_FLUID_BC_Y */

  /* Z direction */
  #ifdef LB_FLUID_BC_Z

   #ifdef LB_FLUID_BC_Z_M_SLIP
    cm.z = -1.0; 
    cm.y = cm.x = 1.0;
   #else
    cm.x = cm.y = cm.z = -1.0; 
   #endif

   for (j = 0; j < LNY + TWO_BRD; j++)                     
    for (i = 0; i < LNX + TWO_BRD; i++){
     if(LNZ_START == 0){
            k = BRD; 
            u[IDX(i,j,k-1)].x =  cm.x*u[IDX(i,j,k)].x;
            u[IDX(i,j,k-1)].y =  cm.y*u[IDX(i,j,k)].y;
            u[IDX(i,j,k-1)].z =  cm.z*u[IDX(i,j,k)].z;
        }
    }

   #ifdef LB_FLUID_BC_Z_P_SLIP
    cp.z = -1.0; 
    cp.y = cp.x = 1.0;
   #else
    cp.x = cp.y = cp.z = -1.0; 
   #endif
  
  for (j = 0; j < LNY + TWO_BRD; j++)                     
   for (i = 0; i < LNX + TWO_BRD; i++){
    if(LNZ_END == NZ){
            k = LNZ+BRD-1;
            u[IDX(i,j,k+1)].x =  cp.x*u[IDX(i,j,k)].x;
            u[IDX(i,j,k+1)].y =  cp.y*u[IDX(i,j,k)].y;
            u[IDX(i,j,k+1)].z =  cp.z*u[IDX(i,j,k)].z;
      }
  }
  #endif /* end ifdef LB_FLUID_BC_Z */
 
 #endif /* endif ifdef LB_FLUID_BC */
#endif /* endif ifdef LB_FLUID */

#ifdef LB_TEMPERATURE
sendrecv_borders_scalar(t);

 #ifdef LB_TEMPERATURE_BC
 /* Y direction */
  #ifdef LB_TEMPERATURE_BC_Y
  // my_double T_wall;

   for (i = 0; i < LNX + TWO_BRD; i++) 			
     for (k = 0; k < LNZ + TWO_BRD; k++){

  if(LNY_START == 0){
	  j = BRD; 
	  T_wall = property.T_bot;
   #ifdef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = 0.0;
   #endif

	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  t[IDX(i,j-1,k)] =  fac*t[IDX(i,j,k)];
  }

  if(LNY_END == NY){
 	  j = LNY+BRD-1; 
	  T_wall = property.T_top;
   #ifdef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = 0.0;
   #endif
 	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  t[IDX(i,j+1,k)] =  fac*t[IDX(i,j,k)];
  }
      
    }
  #endif /* end of ifdef LB_TEMPERATURE_BC_Y */
 #endif /* end of ifdef LB_TEMPERATURE_BC*/

#endif /* LB_TEMPERATURE */


#ifdef LB_SCALAR
sendrecv_borders_scalar(s);

 #ifdef LB_SCALAR_BC
 /* Y direction */
  #ifdef LB_SCALAR_BC_Y
  // my_double S_wall;

   for (i = 0; i < LNX + TWO_BRD; i++) 			
     for (k = 0; k < LNZ + TWO_BRD; k++){

  if(LNY_START == 0){
	  j = BRD; 
	  S_wall = property.S_bot;
   #ifdef LB_SCALAR_FLUCTUATION 
	  S_wall = 0.0;
   #endif

	  fac = 2.0*(S_wall-property.S_ref)/s[IDX(i,j,k)] - 1.0;
	  s[IDX(i,j-1,k)] =  fac*s[IDX(i,j,k)];
  }

  if(LNY_END == NY){
 	  j = LNY+BRD-1; 
	  S_wall = property.S_top;
   #ifdef LB_SCALAR_FLUCTUATION 
	  S_wall = 0.0;
   #endif
 	  fac = 2.0*(S_wall-property.S_ref)/s[IDX(i,j,k)] - 1.0;
	  s[IDX(i,j+1,k)] =  fac*s[IDX(i,j,k)];
  }
      
    }
  #endif /* end of ifdef LB_SCALAR_BC_Y */
 #endif /* end of ifdef LB_SCALAR_BC*/
#endif /* end of ifdef LB_SCALAR */


#ifdef LB_FLUID
#ifdef LB_FLUID_FORCING_LANDSCAPE
  sendrecv_borders_scalar(landscape);
#endif
#endif
}

/* end of function boundary_conditions_hydro() */
