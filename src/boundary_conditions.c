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
 for (i = 0; i < LNX + TWO_BRD; i++) 			
    for (k = 0; k < LNZ + TWO_BRD; k++){
      // this is the good one
      //  for (i = BRD; i < LNX + BRD; i++) 			
      //    for (k = BRD; k < LNZ + BRD; k++){
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
#ifdef LB_TEMPERATURE_BC_Y_VARIABLE
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_top; else T_wall = 0.0;  
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
#ifdef LB_TEMPERATURE_BC_Y_VARIABLE
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = 0.0;  
#endif
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

#ifdef LB_FLUID_BC_XP_OUTLET
	if (j>= BRD && j< LNY+BRD && k >= BRD && k < LNZ+BRD){ 
 
	if(pp==0){
	  rho = 1;
#ifdef GRID_POP_D2Q9
	  vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[4]+2*(rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[5]))/rho;
#endif

#ifdef GRID_POP_D3Q19
	 vel.x = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16]+2*(rhs_p[IDX(i,j,k)].p[7]+rhs_p[IDX(i,j,k)].p[11]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[13]+rhs_p[IDX(i,j,k)].p[1]))/rho; 
#endif
	 vel.y = u[IDX(i, j, k)].y; //0.0;
	 vel.z = u[IDX(i, j, k)].z; //0.0;
	  rhs_p[IDX(i+1,j,k)] = equilibrium_given_velocity(vel,rho);
}
	}  
#else

#ifdef LB_FLUID_BC_XP_SLIP
	//          rhs_p[IDX(i+1,jj,kk)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];

          rhs_p[IDX(i+1,j,k)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#else	 
	  /* NO SLIP */
	  //rhs_p[IDX(i+1,j,k)].p[pp] = rhs_p[IDX(i,jj,kk)].p[inv[pp]];	
	  if (jj>= 0 && jj< LNY+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(i+1,jj,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#endif
#endif
}

if(LNX_START == 0){
  i = BRD;

#ifdef LB_FLUID_BC_XM_OUTLET
	if (j>= BRD && j< LNY+BRD && k >= BRD && k < LNZ+BRD){ 
 
	if(pp==0){
	  rho = 1;
#ifdef GRID_POP_D2Q9
	  vel.x = 1.0- (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[4]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[3]))/rho;
#endif
	  //for GRID_POP_D3Q19 - PENDING
	  vel.y = u[IDX(i, j, k)].y; //0.0;
	  vel.z = u[IDX(i, j, k)].z; //0.0;
	  rhs_p[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho);
}
	}

#else 
#ifdef LB_FLUID_BC_XM_INLET
	if (j>= BRD && j< LNY+BRD && k >= BRD && k < LNZ+BRD){ 
	if(pp==0){

	  /*Poiseuille flow velocity profile - inlet*/
	  // vel.x = -(4.0*(my_double)property.Amp_x*(pow(((my_double)property.SY),-2.0)))*((my_double)center_V[IDX(i,j,k)].y)*(((my_double)center_V[IDX(i,j,k)].y)-((my_double)property.SY));


	  /* Half Poiseuille flow velocity profile - inlet*/
	     vel.x = -(4.0*(my_double)property.Amp_x*(pow((2.0*(my_double)property.SY),-2.0)))*((my_double)center_V[IDX(i,j,k)].y)*(((my_double)center_V[IDX(i,j,k)].y)-(2.0*(my_double)property.SY));

	  /* Constant velocity = Amp_x*/        
	  // vel.x = property.Amp_x;

	  vel.y = property.Amp_y;
	  vel.z = property.Amp_z;

#ifdef GRID_POP_D2Q9 
	  rho = (1.0/(1.0-vel.x))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[4])+2*(rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[3]));
#endif

#ifdef GRID_POP_D3Q19
	  rho = (1.0/(1.0-vel.x))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[15]+rhs_p[IDX(i,j,k)].p[5]+rhs_p[IDX(i,j,k)].p[17]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[18]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[16])+2*(rhs_p[IDX(i,j,k)].p[9]+rhs_p[IDX(i,j,k)].p[12]+rhs_p[IDX(i,j,k)].p[10]+rhs_p[IDX(i,j,k)].p[14]+rhs_p[IDX(i,j,k)].p[2]));
#endif

	  rhs_p[IDX(i-1,j,k)] = equilibrium_given_velocity(vel,rho);
                 }		
	}

#else
#ifdef LB_FLUID_BC_XM_SLIP
  //          rhs_p[IDX(i-1,jj,kk)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];

          rhs_p[IDX(i-1,j,k)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
         //rhs_p[IDX(i-1,j,k)].p[pp] = rhs_p[IDX(i,jj,kk)].p[inv[pp]];
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

#ifdef LB_FLUID_BC_YP_OUTLET
	if (i>= BRD && i< LNX+BRD && k >= BRD && k < LNZ+BRD){ 
	   		   		   		   	   
	if(pp==0){
	  rho = 1.0; //dens[IDX(i, j, k)];//1.0;

#ifdef GRID_POP_D2Q9 
	  	  vel.y = -1.0+ (rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[2]+rhs_p[IDX(i,j,k)].p[6]+2*(rhs_p[IDX(i,j,k)].p[1]+rhs_p[IDX(i,j,k)].p[8]+rhs_p[IDX(i,j,k)].p[7]))/rho;
#endif

		  //for GRID_POP_D3Q19 - PENDING
 
	  vel.x = u[IDX(i, j, k)].x; //0.0;
	  vel.z = u[IDX(i, j, k)].z; //0.0;
	  rhs_p[IDX(i,j+1,k)] = equilibrium_given_velocity(vel,rho);
}		   	   		   		   		   
	}
#else

#ifdef LB_FLUID_BC_YP_SLIP
	  //	  rhs_p[IDX(ii,j+1,kk)].p[inv_y[pp]] = rhs_p[IDX(i,j,k)].p[pp];

          rhs_p[IDX(i,j+1,k)].p[inv_y[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#else
	  /* NO SLIP */
	  if (ii>= 0 && ii< LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j+1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
	  //      rhs_p[IDX(ii,j+1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#endif
#endif
}

if(LNY_START == 0){
  j = BRD; 

  #ifdef LB_FLUID_BC_YM_JET

  if (sqrt(pow(center_V[IDX(i,j,k)].x-(property.SX/4.0), 2.0)+pow(center_V[IDX(i,j,k)].z-(property.SZ/2.0), 2.0)) <= 5.0){     //inlet geometry - circle: POSITION OF CENTER COORDINATES OF THE CIRCLE SPECIFIED HERE ITSELF

    if(pp==0){         
          vel.x =  0.0;    //property.Amp_y;
	  vel.y =  0.1;    //property.Amp_x;
	  vel.z =  0.0;    //property.Amp_z;
	  rho = (1.0/(1.0-vel.y))*((rhs_p[IDX(i,j,k)].p[0]+rhs_p[IDX(i,j,k)].p[6]+rhs_p[IDX(i,j,k)].p[2])+2*(rhs_p[IDX(i,j,k)].p[3]+rhs_p[IDX(i,j,k)].p[4]+rhs_p[IDX(i,j,k)].p[5]));
	  rhs_p[IDX(i,j-1,k)] = equilibrium_given_velocity(vel,rho);
                 } 	
	}
   else
   {
#ifdef LB_FLUID_BC_YM_SLIP

          rhs_p[IDX(i,j-1,k)].p[inv_y[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	// NO SLIP
          if (ii >= 0 && ii < LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j-1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
	  //        rhs_p[IDX(ii,j-1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif
  	   }
#else

#ifdef LB_FLUID_BC_YM_SLIP

          rhs_p[IDX(i,j-1,k)].p[inv_y[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */

          if (ii >= 0 && ii < LNX+TWO_BRD && kk >= 0 && kk < LNZ+TWO_BRD) rhs_p[IDX(ii,j-1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
	  //        rhs_p[IDX(ii,j-1,kk)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#endif

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

#ifdef LB_FLUID_BC_ZP_SLIP
	//          rhs_p[IDX(i+1,jj,kk)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];
 
          rhs_p[IDX(i,j,k+1)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#else	 
	  /* NO SLIP */
	  //rhs_p[IDX(i+1,j,k)].p[pp] = rhs_p[IDX(i,jj,kk)].p[inv[pp]];	
	  if (ii>= 0 && ii< LNX+TWO_BRD && jj >= 0 && jj < LNY+TWO_BRD) rhs_p[IDX(ii,jj,k+1)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];

#endif	
}

if(LNZ_START == 0){
  k = BRD; 

#ifdef LB_FLUID_BC_ZM_SLIP

          rhs_p[IDX(i,j,k-1)].p[inv_x[pp]] = rhs_p[IDX(i,j,k)].p[pp];
#else
	/* NO SLIP */
         //rhs_p[IDX(i-1,j,k)].p[pp] = rhs_p[IDX(i,jj,kk)].p[inv[pp]];
	  if (ii>= 0 && ii< LNX+TWO_BRD && jj >= 0 && jj < LNY+TWO_BRD) rhs_p[IDX(ii,jj,k-1)].p[inv[pp]] = rhs_p[IDX(i,j,k)].p[pp];
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
#ifdef LB_TEMPERATURE_BC_Y_VARIABLE
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_top; else T_wall = 0.0;  
#endif
#else
	  T_wall = 0.0;
#endif

	  fac = 2.0*((T_wall-property.T_ref)- t[IDX(i,j,k)]);

	  for(pp=0;pp<NPOP;pp++){ 
      	    if(c[pp].y>0){
	      ii = i+(int)c[pp].x;
	      kk = k+(int)c[pp].z;	
	  
	      //rhs_g[IDX(ii,j+1,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[pp] + wgt[pp]*fac;
	      //rhs_g[IDX(ii,j+1,kk)].p[inv[pp]] =  wgt[pp]*(fac+t[IDX(i,j,k)]);	  
	      rhs_g[IDX(ii,j+1,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]] + wgt[pp]*fac;
	    }	   
	  }
 }

if(LNY_START == 0){

	  j = BRD; 

#ifndef LB_TEMPERATURE_FLUCTUATION 
	  T_wall = property.T_bot;
#ifdef LB_TEMPERATURE_BC_Y_VARIABLE
	  if(center_V[IDX(i,j,k)].x < property.SX/2.0)  T_wall = property.T_bot; else T_wall = 0.0;  
#endif
#else
	  T_wall = 0.0;
#endif

	  fac = 2.0*((T_wall-property.T_ref)-t[IDX(i,j,k)]);
	  
	  //effDT=0.0;
	  //for(pp=0;pp<NPOP;pp++) if(c[pp].y<0)  effDT += wgt[inv[pp]];

	  for(pp=0;pp<NPOP;pp++){
	    if(c[pp].y<0){
	      ii = i+(int)c[pp].x;
	      kk = k+(int)c[pp].z;
            //rhs_g[IDX(ii,j-1,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[pp] + wgt[pp]*fac;
	    //  rhs_g[IDX(ii,j-1,kk)].p[inv[pp]] =  wgt[pp]*(fac+t[IDX(i,j,k)]);	
	    //  rhs_g[IDX(ii,j-1,kk)].p[inv[pp]] =  rhs_g[IDX(i,j,k)].p[inv[pp]]*(fac/t[IDX(i,j,k)] + 1.0);
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


#ifdef LB_SCALAR_BC_YP_OUTLET
	  /* this is the default fixed-at-wall bc */
	  fac = 0.0;
#else
	  /* this is the default fixed-at-wall bc */
	  fac = 2.0*((S_wall-property.S_ref)- s[IDX(i,j,k)]);
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
	  fac = 2.0*((S_wall-property.S_ref)-s[IDX(i,j,k)]);
	  

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
#endif /* end of scalar part */

}/* end of bc for streaming */
#endif
//#endif



/***********************************************/

#ifdef LB_FLUID_BC
void boundary_conditions_for_equilibrium(){

  int i,j,k,pp;
  vector vel;
  my_double rho;
  pop f_eq;

  /* y direction  velocity */
#ifdef LB_FLUID_BC_Y
 
  for (i = BRD; i < LNX + BRD; i++) 			
    for (k = BRD; k < LNZ + BRD; k++){
      for(pp=0;pp<NPOP;pp++){

      /*  top  */
	if(LNY_END == NY){
	   j = LNY+BRD-1;

#ifdef LB_FLUID_BC_YP_SLIP
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

#ifdef LB_FLUID_BC_YM_SLIP
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
/***************************************************************************************************/


}
#endif

