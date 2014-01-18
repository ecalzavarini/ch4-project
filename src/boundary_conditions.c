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
	
#ifdef LB_FLUID_BC_XP_OUTLET
	
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
#ifdef LB_FLUID_BC_XP_SLIP

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

#ifdef LB_FLUID_BC_XM_INLET
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
#ifdef LB_FLUID_BC_XM_SLIP
		
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

  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){
      for(pp=0;pp<NPOP;pp++){

if(LNY_END == NY){
	j = LNY+BRD-1; 

#ifdef LB_FLUID_BC_YP_OUTLET
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
#ifdef METHOD_MYQUICK
       	p[IDX(i,j+2,k)].p[pp] = p[IDX(i,j-1,k)].p[inv[pp]];
#endif

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
#ifdef METHOD_MYQUICK
       	  p[IDX(i,j-2,k)].p[pp] = p[IDX(i,j+1,k)].p[inv[pp]];
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
#else
	  T_wall = 0.0;
#endif

	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,j,k)] - 1.0;
	  for(pp=0;pp<NPOP;pp++) g[IDX(i,j-1,k)].p[pp] =  fac*g[IDX(i,j,k)].p[pp];

     
#ifdef METHOD_MYQUICK 
	  
#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif   
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

  /* communications to be done in any case  (especially for pbc)*/

#ifdef LB_FLUID
sendrecv_borders_pop(rhs_p);
#endif

#ifdef LB_TEMPERATURE
sendrecv_borders_pop(rhs_g);
#endif

#ifdef LB_SCALAR
sendrecv_borders_pop(rhs_h);
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
	  /* NO SLIP */
	  rhs_p[IDX(i+1,j,k)].p[pp] = rhs_p[IDX(i,jj,kk)].p[inv[pp]];	
}

if(LNX_START == 0){
  i = BRD; 
	/* NO SLIP */
         rhs_p[IDX(i-1,j,k)].p[pp] = rhs_p[IDX(i,jj,kk)].p[inv[pp]];
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
	  /* NO SLIP */
	  rhs_p[IDX(i,j+1,k)].p[pp] = rhs_p[IDX(ii,j,kk)].p[inv[pp]];	
}

if(LNY_START == 0){
  j = BRD; 
	/* NO SLIP */
         rhs_p[IDX(i,j-1,k)].p[pp] = rhs_p[IDX(ii,j,kk)].p[inv[pp]];
 }

      }/* for pp */
    }/* for i,k */
#endif





#ifdef LB_TEMPERATURE	
  pop g_eq, g_eq_w;
  my_double effDT, rho2;
  my_double T_wall;
  my_double fac;

  /************************/
	/* Y direction */
#ifdef LB_TEMPERATURE_BC_Y


  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){


if(LNY_END == NY){

          j = LNY+BRD-1; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_top;
#else
	  T_wall = 0.0;
#endif

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].y<0){
	  ii = i+(int)c[pp].x;
	  kk = k+(int)c[pp].z;	
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(ii,j,kk)] - 1.0;
	  rhs_g[IDX(i,j+1,k)].p[pp] =  rhs_g[IDX(ii,j,kk)].p[inv[pp]] + wgt[pp]*fac*T_wall;//t[IDX(ii,j,kk)];
	    }
	  }
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#else
	  T_wall = 0.0;
#endif

	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].y>0){
	  ii = i+(int)c[pp].x;
	  kk = k+(int)c[pp].z;	
	  fac = 2.0*(T_wall-property.T_ref)/t[IDX(ii,j,kk)] - 1.0;
	  rhs_g[IDX(i,j-1,k)].p[pp] =  rhs_g[IDX(ii,j,kk)].p[inv[pp]] + wgt[pp]*fac*T_wall;//t[IDX(ii,j,kk)];	    
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

	  for(pp=0;pp<NPOP;pp++){ 
	    if(c[pp].x<0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;	
	  //fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,jj,kk)] - 1.0;
	  rhs_g[IDX(i+1,j,k)].p[pp] =  rhs_g[IDX(i,jj,kk)].p[pp]; //rhs_g[IDX(i,jj,kk)].p[inv[pp]] + wgt[pp]*fac*T_wall;
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

	  
	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].x>0){
	  jj = j+(int)c[pp].y;
	  kk = k+(int)c[pp].z;	
	  //fac = 2.0*(T_wall-property.T_ref)/t[IDX(i,jj,kk)] - 1.0;
	  rhs_g[IDX(i-1,j,k)].p[pp] =  rhs_g[IDX(i,jj,kk)].p[pp]; //rhs_g[IDX(i,jj,kk)].p[inv[pp]] + wgt[pp]*fac*T_wall;
	    }
	  }
 }
      
}
#endif

#endif



}/* end of bc for streaming */
#endif
//#endif



/***********************************************/

#ifdef LB_FLUID_BC
void boundary_conditions_for_equilibrium(){

  int i,j,k;
  vector vel;
  my_double rho;


  /* y direction */
#ifdef LB_FLUID_BC_Y
 
  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){

      /*  top  */
if(LNY_END == NY){
	j = LNY+BRD-1;

	/* no slip is the default */
	  vel.x = -u[IDX(i, j, k)].x;
	  vel.y = -u[IDX(i, j, k)].y;
	  vel.z = -u[IDX(i, j, k)].z;
	  rho = dens[IDX(i, j, k)];
	  p_eq[IDX(i,j+1,k)] = equilibrium_given_velocity(vel,rho);

	  vel.x = -u[IDX(i, j-1, k)].x;
	  vel.y = -u[IDX(i, j-1, k)].y;
	  vel.z = -u[IDX(i, j-1, k)].z;
	  rho = dens[IDX(i, j-1, k)];
	  p_eq[IDX(i,j+2,k)] = equilibrium_given_velocity(vel,rho);
 }

    /*  bottom  */
if(LNY_START == 0){
       j = BRD; 

	/* no slip is the default */
	  vel.x = -u[IDX(i, j, k)].x;
	  vel.y = -u[IDX(i, j, k)].y;
	  vel.z = -u[IDX(i, j, k)].z;
	  rho = dens[IDX(i, j, k)];
	  p_eq[IDX(i,j-1,k)] = equilibrium_given_velocity(vel,rho);

	  vel.x = -u[IDX(i, j+1, k)].x;
	  vel.y = -u[IDX(i, j+1, k)].y;
	  vel.z = -u[IDX(i, j+1, k)].z;
	  rho = dens[IDX(i, j+1, k)];
	  p_eq[IDX(i,j-2,k)] = equilibrium_given_velocity(vel,rho); 
 }

    }
#endif

}

#endif
