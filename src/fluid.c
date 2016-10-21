#include "common_object.h"


void compute_advection(pop *f, pop *rhs_f, my_double tau, pop *f_eq, char which_pop){

  int i,j,k,pp;
#ifndef METHOD_STREAMING /* if streaming is not defined skip all this */
  my_double adv,aux;
#ifdef DEBUG
	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;
#endif

   my_double fac, dt_over_tau;
   pop ff_eq;
#ifdef METHOD_REDEFINED_POP
  my_double fac2;
  pop f0;
  pop fxp1, fxm1 ,fxm2 ,fxp2;
  pop fyp1, fym1, fym2, fyp2; 
  pop fzp1, fzm1 , fzm2, fzp2;
#ifdef METHOD_REDEFINED_POP_GUO
  vector d;
  my_double rho, ux,uy,uz,cu;
#endif
#endif

#ifdef DEBUG
	/* Each processor prints its mesh */
	sprintf(fnamein, "pop_p.%d.out", me);
	fout = fopen(fnamein, "w");
	
	for (k = 0; k < LNZ + TWO_BRD; k++)
	     for (j = 0; j < LNY + TWO_BRD; j++)
	       for (i = 0; i < LNX + TWO_BRD; i++){
		    //for(pp=0;pp<NPOP;pp++)
		    pp=7;
	fprintf(fout, "%d %d %d %d %e\n", i, j, k, pp, f[IDX(i,j,k)].p[pp] );
	       }
	fclose(fout);
#endif

#ifdef METHOD_REDEFINED_POP
	/* set to zero f_aux , not necessary for the moment*/
	/*
 for(k=0;k<LNZ+TWO_BRD;k++)
   for(j=0;j<LNY+TWO_BRD;j++)
     for(i=0;i<LNX+TWO_BRD;i++){  
       for(pp=0;pp<NPOP;pp++)        
	 f_aux[IDX(i,j,k)].p[pp]  = 0.0;
    }
	*/
      /* We store the equilibrium distribution in all points */
 for(k=BRD;k<LNZ+BRD;k++)
   for(j=BRD;j<LNY+BRD;j++)
    for(i=BRD;i<LNX+BRD;i++){ 
      	f_eq[IDX(i,j,k)]=equilibrium(f,i,j,k);
      }
  /* send the borders, needed to compute the advection */ 
 // OUT     sendrecv_borders_pop(f_eq);

 //#ifdef LB_FLUID_BC
 // OUT    boundary_conditions_for_equilibrium(which_pop);
 //#endif

 /* The population to be advected is different (we call it here f_aux)*/
 /* it is:  f_aux = f + (dt/(2*tau))*(f_eq-f) */
 /* we prepare such a population here */
 fac = 0.5*(property.time_dt/tau); 

 //OUT for(k=0;k<LNZ+TWO_BRD;k++)
 //OUT   for(j=0;j<LNY+TWO_BRD;j++)
 //OUT    for(i=0;i<LNX+TWO_BRD;i++){  
 for(k=BRD;k<LNZ+BRD;k++)
   for(j=BRD;j<LNY+BRD;j++)
    for(i=BRD;i<LNX+BRD;i++){ 
       for(pp=0;pp<NPOP;pp++)        
	 f_aux[IDX(i,j,k)].p[pp]  = f[IDX(i,j,k)].p[pp] + fac*( f_eq[IDX(i,j,k)].p[pp] - f[IDX(i,j,k)].p[pp] );
    }
#ifdef METHOD_REDEFINED_POP_GUO
#ifdef LB_FLUID
#ifdef LB_FLUID_FORCING
if(which_pop == 'p'){
  /* this method is going to be very expensive */

fac = 0.5*property.time_dt;

  /* we update density and velocity (in the bulk) */
 for(k=BRD;k<LNZ+BRD;k++)
   for(j=BRD;j<LNY+BRD;j++)
    for(i=BRD;i<LNX+BRD;i++){ 
   dens[IDX(i, j, k)] = m(p[IDX(i, j, k)]);
   /* TO BE USED if the force is a force */
   //u[IDX(i, j, k)].x= ( mvx(p[IDX(i, j, k)]) + 0.5*property.time_dt*force[IDX(i, j, k)].x )/dens[IDX(i, j, k)];
   //u[IDX(i, j, k)].y= ( mvy(p[IDX(i, j, k)]) + 0.5*property.time_dt*force[IDX(i, j, k)].y )/dens[IDX(i, j, k)];
   //u[IDX(i, j, k)].z= ( mvz(p[IDX(i, j, k)]) + 0.5*property.time_dt*force[IDX(i, j, k)].z )/dens[IDX(i, j, k)];
   /* TO BE USED if the force is a force per unit mass i.e. acceleration */
   u[IDX(i, j, k)].x= mvx(p[IDX(i, j, k)])/dens[IDX(i, j, k)] + fac*force[IDX(i, j, k)].x ;
   u[IDX(i, j, k)].y= mvy(p[IDX(i, j, k)])/dens[IDX(i, j, k)] + fac*force[IDX(i, j, k)].y ;
   u[IDX(i, j, k)].z= mvz(p[IDX(i, j, k)])/dens[IDX(i, j, k)] + fac*force[IDX(i, j, k)].z ;	
   }
  
  fac = 0.5*property.time_dt;

 for(k=0;k<LNZ+TWO_BRD;k++)
   for(j=0;j<LNY+TWO_BRD;j++)
     for(i=0;i<LNX+TWO_BRD;i++){  
	      rho = dens[IDX(i,j,k)];
	      ux=u[IDX(i,j,k)].x;
	      uy=u[IDX(i,j,k)].y;
	      uz=u[IDX(i,j,k)].z;
	        for(pp=0;pp<NPOP;pp++){    
                cu = (c[pp].x*ux + c[pp].y*uy + c[pp].z*uz);
                d.x = (c[pp].x-ux)*invcs2 + c[pp].x*cu*invcs4;
                d.y = (c[pp].y-uy)*invcs2 + c[pp].y*cu*invcs4;
                d.z = (c[pp].z-uz)*invcs2 + c[pp].z*cu*invcs4;
       f_aux[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].x*d.x;
       f_aux[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].y*d.y;
       f_aux[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].z*d.z;  
		}
       }
 }
#endif /* end LB_FLUID_FORCING */
#endif /* end LB_FLUID */
#endif /* end METHOD_REDEFINED_POP_GUO */
/* first PBC */
 sendrecv_borders_pop(f_aux);
 /* then real BC */
#ifdef LB_FLUID_BC
 boundary_conditions_for_advection(f_aux, which_pop);
#endif
#endif /* end METHOD_REDEFINED_POP */

#ifdef METHOD_COLLISION_IMPLICIT  /* this is for Euler implicit time stepping */
      /* We change f in  ( f + (dt/tau)f_eq ) /(1+dt/tau)  */
	dt_over_tau = property.time_dt/tau;  
	fac = 1./(1. + dt_over_tau);

/* We store the equilibrium distribution in all points */
 for(k=BRD;k<LNZ+BRD;k++)
   for(j=BRD;j<LNY+BRD;j++)
    for(i=BRD;i<LNX+BRD;i++){ 
      	f_eq[IDX(i,j,k)]=equilibrium(f,i,j,k);
      }

#ifdef LB_FLUID_BC
	/* bc for equilibrium function */
        boundary_conditions_for_equilibrium(which_pop);
#endif
  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	for(pp=1;pp<NPOP;pp++) f[IDX(i,j,k)].p[pp] = fac*(f[IDX(i,j,k)].p[pp] + dt_over_tau* f_eq[IDX(i,j,k)].p[pp] );
      }
  /* send the borders, needed to compute the advection */ 
        sendrecv_borders_pop(f);
#endif

  /* check this index */
  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

	/* population 0 does not advect anyhow*/
	rhs_f[IDX(i,j,k)].p[0] = 0.0;
	/* now taking care of pop 1 to 18 */
 	for(pp=1;pp<NPOP;pp++){

#ifdef METHOD_CENTERED
#ifdef METHOD_REDEFINED_POP
 /* The auxiliary population is advected f_aux = f + (dt/(2*tau))*(f_eq-f) */
	  adv=0.0;
	  adv = compute_flux_with_central_difference(f_aux,i,j,k,pp);
#else		  
	  /* central difference scheme */
	  adv=0.0;
	  adv = compute_flux_with_central_difference(f,i,j,k,pp);
#endif
#endif

#ifdef METHOD_UPWIND
#ifdef METHOD_REDEFINED_POP
 /* The auxiliary population is advected f_aux = f + (dt/(2*tau))*(f_eq-f) */
	  adv=0.0;
	  adv = compute_flux_with_upwind_first(f_aux,i,j,k,pp);
#else	
 /* first order upwind scheme */
	  adv=0.0;
	  adv = compute_flux_with_upwind_first(f,i,j,k,pp);	  
#endif
#endif

 /* QUICK */
#ifdef METHOD_MYQUICK
#ifdef METHOD_REDEFINED_POP
 /* The auxiliary population is advected f_aux = f + (dt/(2*tau))*(f_eq-f) */
	  adv=0.0;
	  adv = compute_flux_with_quick(f_aux,i,j,k,pp);
#else	
/* quick scheme */
	   adv=0.0;
	   adv = compute_flux_with_quick(f,i,j,k,pp);
          #ifdef METHOD_UPWIND_LINEAR
          adv = compute_flux_with_upwind_linear(f,i,j,k,pp);
          #endif	   

#ifdef METHOD_MYQUICK_LIMITER
	   adv=0.0;
	   adv = compute_flux_with_limiters(f,i,j,k,pp);
#endif

#endif /* end of redefined pop method */
#endif /* end of quick */


#ifdef METHOD_MIXED
adv=0.0;

if((LNY_END == NY && j==LNY+BRD-1) || (LNY_START == 0 && j==BRD)){
//if((LNY_END == NY && j>=LNY+BRD-2) || (LNY_START == 0 && j<=BRD+1)){
//if((LNY_END == NY && j==LNY+BRD-1  && c[pp].y < 0) || (LNY_START == 0 && j==BRD  && c[pp].y > 0)){
//if(c[pp].x*c[pp].x + c[pp].y*c[pp].y + c[pp].z*c[pp].z ==1.0){
#ifndef METHOD_REDEFINED_POP /* redefined pop is off */
#ifdef METHOD_UPWIND
       adv = compute_flux_with_upwind_first(f,i,j,k,pp);
#endif
#ifdef METHOD_UPWIND_LINEAR
       adv = compute_flux_with_upwind_linear(f,i,j,k,pp);
#endif
	 }else{
#ifdef METHOD_CENTERED
        adv = compute_flux_with_central_difference(f,i,j,k,pp);
#endif
#ifdef METHOD_MYQUICK
        adv = compute_flux_with_quick(f,i,j,k,pp);
#endif
#else  /* redefined pop is on */
#ifdef METHOD_UPWIND
       adv = compute_flux_with_upwind_first(f_aux,i,j,k,pp);
#endif
#ifdef METHOD_UPWIND_LINEAR
       adv = compute_flux_with_upwind_linear(f_aux,i,j,k,pp);
#endif
	 }else{
#ifdef METHOD_CENTERED
        adv = compute_flux_with_central_difference(f_aux,i,j,k,pp);
#endif
#ifdef METHOD_MYQUICK
        adv = compute_flux_with_quick(f_aux,i,j,k,pp);
#endif
#endif
 }
#endif

/* with minus sign because we add it to the right hand side */
 rhs_f[IDX(i,j,k)].p[pp] = -adv;

 #ifdef DEBUG_HARD
 if(ROOT && itime ==1000) fprintf(stderr, " %d %d %d %d adv %e\n",i,j,k,pp,adv);
 #endif 

	}/* for pp */
      }/* for i, j , k */

#endif /* end of  ifndef METHOD_STREAMING */
#ifdef METHOD_STREAMING
#ifndef METHOD_STREAMING_INTERPOLATE
 for(k=BRD;k<LNZ+BRD;k++){
   for(j=BRD;j<LNY+BRD;j++){
      for(i=BRD;i<LNX+BRD;i++){ 

	/* in case of normal streaming we just make a copy of the field */
       if(property.time_dt == 1.0){ 
	 /* when the time step is = 1 */
	   rhs_f[IDX(i,j,k)]  = f[IDX(i,j,k)];
       }else{	 
	 /* when the time step is not 1 */  
	   for(pp=0;pp<NPOP;pp++) rhs_f[IDX(i,j,k)].p[pp]  = f[IDX(i,j,k)].p[pp]/property.time_dt;
       }

      }/* i */
   }/* j */
 }/* k */
#else
 /* here METHOD_STREAMING_INTERPOLATE is active */
 /* Interpolation method when the grid is cartesian rectangular, the quick algorithm can be written in a more compact form */ 
 if(coeff_xp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_xp[IDX(i,j,k)].p[pp] > 0.0){
   aux= 1.0 + interp2_xp[IDX(i,j,k)] - interp3_xm[IDX(i,j,k)];
   adv += coeff_xp[IDX(i,j,k)].p[pp]*( 
				      interp_xp[IDX(i,j,k)]*f[IDX(i+1,j,k)].p[pp] 
				      + (aux - interp_xp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] 
				      - (aux + interp4_xm[IDX(i,j,k)])*f[IDX(i-1,j,k)].p[pp] 
				      + interp4_xm[IDX(i,j,k)]*f[IDX(i-2,j,k)].p[pp] 
				       );
   }else{
   aux = 1.0 + interp2_xm[IDX(i,j,k)] - interp3_xp[IDX(i,j,k)];
   adv += coeff_xp[IDX(i,j,k)].p[pp]*( 
				      (aux + interp4_xp[IDX(i,j,k)])*f[IDX(i+1,j,k)].p[pp] 
				      - interp4_xp[IDX(i,j,k)]*f[IDX(i+2,j,k)].p[pp] 
				      - interp_xm[IDX(i,j,k)]*f[IDX(i-1,j,k)].p[pp] 
				      - (aux - interp_xm[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] 
				       );
   }
 }


 if(coeff_yp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_yp[IDX(i,j,k)].p[pp] > 0.0){
   aux = 1.0  + interp2_yp[IDX(i,j,k)] - interp3_ym[IDX(i,j,k)];
   adv += coeff_yp[IDX(i,j,k)].p[pp]*( 
				      interp_yp[IDX(i,j,k)]*f[IDX(i,j+1,k)].p[pp] 
				      + (aux - interp_yp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] 
				      - (aux + interp4_ym[IDX(i,j,k)])*f[IDX(i,j-1,k)].p[pp] 
				      + interp4_ym[IDX(i,j,k)]*f[IDX(i,j-2,k)].p[pp] 
				       );
   }else{
   aux = 1.0 + interp2_ym[IDX(i,j,k)] - interp3_yp[IDX(i,j,k)];
   adv += coeff_yp[IDX(i,j,k)].p[pp]*( 
				      (aux + interp4_yp[IDX(i,j,k)])*f[IDX(i,j+1,k)].p[pp] 
				      - interp4_yp[IDX(i,j,k)]*f[IDX(i,j+2,k)].p[pp] 
				      - interp_ym[IDX(i,j,k)]*f[IDX(i,j-1,k)].p[pp] 
				      - (aux - interp_ym[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] 
				       );
   }
 }


 if(coeff_zp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_zp[IDX(i,j,k)].p[pp] > 0.0){
   aux = 1.0 + interp2_zp[IDX(i,j,k)] - interp3_zm[IDX(i,j,k)];
   adv += coeff_zp[IDX(i,j,k)].p[pp]*( 
				      interp_zp[IDX(i,j,k)]*f[IDX(i,j,k+1)].p[pp] 
				      + (aux - interp_zp[IDX(i,j,k)] )*f[IDX(i,j,k)].p[pp] 
				      - (aux + interp4_zm[IDX(i,j,k)])*f[IDX(i,j,k-1)].p[pp] 
				      + interp4_zm[IDX(i,j,k)]*f[IDX(i,j,k-2)].p[pp]
				       );
   }else{
   aux = 1.0 + interp2_zm[IDX(i,j,k)] - interp3_zp[IDX(i,j,k)];
   adv += coeff_zp[IDX(i,j,k)].p[pp]*( 
				      (aux + interp4_zp[IDX(i,j,k)])*f[IDX(i,j,k+1)].p[pp] 
				      - interp4_zp[IDX(i,j,k)]*f[IDX(i,j,k+2)].p[pp] 
				      -interp_zm[IDX(i,j,k)]*f[IDX(i,j,k-1)].p[pp] 
				      - (aux - interp_zm[IDX(i,j,k)] )*f[IDX(i,j,k)].p[pp] 
				       );
   }
 }
#endif
#endif
}/* end of compute_advection */

/* Advection methods */

#ifdef METHOD_CENTERED
my_double compute_flux_with_central_difference(pop * f, int i, int j, int k, int pp){

  my_double adv=0.0;

 /* centered difference scheme */
  
  /* for equispaced grid */
  /*	  
	  adv += coeff_xp[IDX(i,j,k)].p[pp]*0.5*(f[IDX(i+1,j,k)].p[pp] + f[IDX(i,j,k)].p[pp])+ 
	         coeff_xm[IDX(i,j,k)].p[pp]*0.5*(f[IDX(i-1,j,k)].p[pp] + f[IDX(i,j,k)].p[pp]);
	  
	  adv += coeff_yp[IDX(i,j,k)].p[pp]*0.5*(f[IDX(i,j+1,k)].p[pp] + f[IDX(i,j,k)].p[pp])+
	         coeff_ym[IDX(i,j,k)].p[pp]*0.5*(f[IDX(i,j-1,k)].p[pp] + f[IDX(i,j,k)].p[pp]);

          adv +=  coeff_zp[IDX(i,j,k)].p[pp]*0.5*(f[IDX(i,j,k+1)].p[pp] + f[IDX(i,j,k)].p[pp])+
                  coeff_zm[IDX(i,j,k)].p[pp]*0.5*(f[IDX(i,j,k-1)].p[pp] + f[IDX(i,j,k)].p[pp]);
	  */	  

  /* for stretched cartesian grid */	  	  
	  adv += coeff_xp[IDX(i,j,k)].p[pp]*( (1.0 - interp_xp[IDX(i,j,k)])* f[IDX(i+1,j,k)].p[pp] + interp_xp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] )+ 
	         coeff_xm[IDX(i,j,k)].p[pp]*( (1.0 - interp_xm[IDX(i,j,k)])* f[IDX(i-1,j,k)].p[pp] + interp_xm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] );
	  
	  adv += coeff_yp[IDX(i,j,k)].p[pp]*( (1.0 - interp_yp[IDX(i,j,k)])* f[IDX(i,j+1,k)].p[pp] + interp_yp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] )+
	         coeff_ym[IDX(i,j,k)].p[pp]*( (1.0 - interp_ym[IDX(i,j,k)])* f[IDX(i,j-1,k)].p[pp] + interp_ym[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] );

          adv += coeff_zp[IDX(i,j,k)].p[pp]*( (1.0 - interp_zp[IDX(i,j,k)])* f[IDX(i,j,k+1)].p[pp] + interp_zp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] )+
                 coeff_zm[IDX(i,j,k)].p[pp]*( (1.0 - interp_zm[IDX(i,j,k)])* f[IDX(i,j,k-1)].p[pp] + interp_zm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] );

	  return adv;
}	 
#endif


#ifdef METHOD_UPWIND
my_double compute_flux_with_upwind_first(pop * f, int i, int j, int k, int pp){

  my_double adv=0.0;
  my_double dist,fac;
  pop c2;

#ifndef METHOD_UPWIND_SKEW

 /* first order upwind scheme */   
 if(coeff_xp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_xp[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_xp[IDX(i,j,k)].p[pp]*f[IDX(i+1,j,k)].p[pp];

 if(coeff_xm[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_xm[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_xm[IDX(i,j,k)].p[pp]*f[IDX(i-1,j,k)].p[pp];

 if(coeff_yp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_yp[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_yp[IDX(i,j,k)].p[pp]*f[IDX(i,j+1,k)].p[pp];

 if(coeff_ym[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_ym[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_ym[IDX(i,j,k)].p[pp]*f[IDX(i,j-1,k)].p[pp];

 if(coeff_zp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_zp[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_zp[IDX(i,j,k)].p[pp]*f[IDX(i,j,k+1)].p[pp];

 if(coeff_zm[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_zm[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_zm[IDX(i,j,k)].p[pp]*f[IDX(i,j,k-1)].p[pp];

#else

  /* skew upwind */
  c2.p[pp]=c[pp].x*c[pp].x + c[pp].y*c[pp].y + c[pp].z*c[pp].z;

  if(c2.p[pp]==1){

 /* first order upwind scheme */   
 if(coeff_xp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_xp[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_xp[IDX(i,j,k)].p[pp]*f[IDX(i+1,j,k)].p[pp];

 if(coeff_xm[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_xm[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_xm[IDX(i,j,k)].p[pp]*f[IDX(i-1,j,k)].p[pp];

 if(coeff_yp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_yp[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_yp[IDX(i,j,k)].p[pp]*f[IDX(i,j+1,k)].p[pp];

 if(coeff_ym[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_ym[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_ym[IDX(i,j,k)].p[pp]*f[IDX(i,j-1,k)].p[pp];

 if(coeff_zp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_zp[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_zp[IDX(i,j,k)].p[pp]*f[IDX(i,j,k+1)].p[pp];

 if(coeff_zm[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_zm[IDX(i,j,k)].p[pp]*f[IDX(i,j,k)].p[pp];
 else
   adv += coeff_zm[IDX(i,j,k)].p[pp]*f[IDX(i,j,k-1)].p[pp];
  }
  
  if(c2.p[pp]==2.0){

    fac=1.0; //0.5*(sqrt(5.)-1.); // why?

    if(coeff_xp[IDX(i,j,k)].p[pp] >= 0.0){
      dist = interp_xp[IDX(i,j,k)];      
      if( c[pp].y != 0.0 && c[pp].y < 0.0)  dist /= interp3_yp[IDX(i,j,k)]; 
      if( c[pp].y != 0.0 && c[pp].y > 0.0)  dist /= interp3_ym[IDX(i,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z < 0.0)  dist /= interp3_zp[IDX(i,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z > 0.0)  dist /= interp3_zm[IDX(i,j,k)];                       
   adv += fac*coeff_xp[IDX(i,j,k)].p[pp]*((1.0 -  dist)*f[IDX(i,j,k)].p[pp] + dist*f[IDX(i,j+c[pp].y,k+c[pp].z)].p[pp] );
    }else{
      dist = interp_xm[IDX(i+1,j,k)];
      if( c[pp].y != 0.0 && c[pp].y < 0.0)  dist /= interp3_yp[IDX(i+1,j,k)]; 
      if( c[pp].y != 0.0 && c[pp].y > 0.0)  dist /= interp3_ym[IDX(i+1,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z < 0.0)  dist /= interp3_zp[IDX(i+1,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z > 0.0)  dist /= interp3_zm[IDX(i+1,j,k)]; 
   adv += fac*coeff_xp[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i+1,j,k)].p[pp] + dist*f[IDX(i+1,j+c[pp].y,k+c[pp].z)].p[pp] );
    }

    if(coeff_xm[IDX(i,j,k)].p[pp] >= 0.0){
      dist = interp_xm[IDX(i,j,k)];
      if( c[pp].y != 0.0 && c[pp].y < 0.0)  dist /= interp3_yp[IDX(i,j,k)]; 
      if( c[pp].y != 0.0 && c[pp].y > 0.0)  dist /= interp3_ym[IDX(i,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z < 0.0)  dist /= interp3_zp[IDX(i,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z > 0.0)  dist /= interp3_zm[IDX(i,j,k)];    
   adv += fac*coeff_xm[IDX(i,j,k)].p[pp]*((1.0 -  dist)*f[IDX(i,j,k)].p[pp] + dist*f[IDX(i,j+c[pp].y,k+c[pp].z)].p[pp] );
    }else{
      dist =  interp_xp[IDX(i-1,j,k)];
     if( c[pp].y != 0.0 && c[pp].y < 0.0)  dist /= interp3_yp[IDX(i-1,j,k)]; 
     if( c[pp].y != 0.0 && c[pp].y > 0.0)  dist /= interp3_ym[IDX(i-1,j,k)]; 
     if( c[pp].z != 0.0 && c[pp].z < 0.0)  dist /= interp3_zp[IDX(i-1,j,k)]; 
     if( c[pp].z != 0.0 && c[pp].z > 0.0)  dist /= interp3_zm[IDX(i-1,j,k)]; 
   adv += fac*coeff_xm[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i-1,j,k)].p[pp] + dist*f[IDX(i-1,j+c[pp].y,k+c[pp].z)].p[pp] );
    }

    if(coeff_yp[IDX(i,j,k)].p[pp] >= 0.0){
      dist = interp_yp[IDX(i,j,k)];
      if( c[pp].x != 0.0 && c[pp].x < 0.0)  dist /= interp3_xp[IDX(i,j,k)]; 
      if( c[pp].x != 0.0 && c[pp].x > 0.0)  dist /= interp3_xm[IDX(i,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z < 0.0)  dist /= interp3_zp[IDX(i,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z > 0.0)  dist /= interp3_zm[IDX(i,j,k)]; 
      adv += fac*coeff_yp[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i,j,k)].p[pp] + dist*f[IDX(i+c[pp].x,j,k+c[pp].z)].p[pp] );
    }else{
      dist = interp_ym[IDX(i,j+1,k)];
      if( c[pp].x != 0.0 && c[pp].x < 0.0)  dist /= interp3_xp[IDX(i,j+1,k)]; 
      if( c[pp].x != 0.0 && c[pp].x > 0.0)  dist /= interp3_xm[IDX(i,j+1,k)]; 
      if( c[pp].z != 0.0 && c[pp].z < 0.0)  dist /= interp3_zp[IDX(i,j+1,k)]; 
      if( c[pp].z != 0.0 && c[pp].z > 0.0)  dist /= interp3_zm[IDX(i,j+1,k)];   
  adv += fac*coeff_yp[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i,j+1,k)].p[pp] + dist*f[IDX(i+c[pp].x,j+1,k+c[pp].z)].p[pp] );
    }

    if(coeff_ym[IDX(i,j,k)].p[pp] >= 0.0){
      dist = interp_ym[IDX(i,j,k)];
      if( c[pp].x != 0.0 && c[pp].x < 0.0)  dist /= interp3_xp[IDX(i,j,k)]; 
      if( c[pp].x != 0.0 && c[pp].x > 0.0)  dist /= interp3_xm[IDX(i,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z < 0.0)  dist /= interp3_zp[IDX(i,j,k)]; 
      if( c[pp].z != 0.0 && c[pp].z > 0.0)  dist /= interp3_zm[IDX(i,j,k)];   
  adv += fac*coeff_ym[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i,j,k)].p[pp] + dist*f[IDX(i+c[pp].x,j,k+c[pp].z)].p[pp] );
    }else{
      dist = interp_yp[IDX(i,j-1,k)];
      if( c[pp].x != 0.0 && c[pp].x < 0.0)  dist /= interp3_xp[IDX(i,j-1,k)]; 
      if( c[pp].x != 0.0 && c[pp].x > 0.0)  dist /= interp3_xm[IDX(i,j-1,k)]; 
      if( c[pp].z != 0.0 && c[pp].z < 0.0)  dist /= interp3_zp[IDX(i,j-1,k)]; 
      if( c[pp].z != 0.0 && c[pp].z > 0.0)  dist /= interp3_zm[IDX(i,j-1,k)];  
  adv += fac*coeff_ym[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i,j-1,k)].p[pp] + dist*f[IDX(i+c[pp].x,j-1,k+c[pp].z)].p[pp] );
    }

    if(coeff_zp[IDX(i,j,k)].p[pp] >= 0.0){
      dist = interp_zp[IDX(i,j,k)];
      if( c[pp].y != 0.0 && c[pp].y < 0.0)  dist /= interp3_yp[IDX(i,j,k)]; 
      if( c[pp].y != 0.0 && c[pp].y > 0.0)  dist /= interp3_ym[IDX(i,j,k)]; 
      if( c[pp].x != 0.0 && c[pp].x < 0.0)  dist /= interp3_xp[IDX(i,j,k)]; 
      if( c[pp].x != 0.0 && c[pp].x > 0.0)  dist /= interp3_xm[IDX(i,j,k)]; 
  adv += fac*coeff_zp[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i,j,k)].p[pp] + dist*f[IDX(i+c[pp].x,j+c[pp].y,k)].p[pp] );
    }else{
     dist = interp_zm[IDX(i,j,k+1)];
      if( c[pp].y != 0.0 && c[pp].y < 0.0)  dist /= interp3_yp[IDX(i,j,k+1)]; 
      if( c[pp].y != 0.0 && c[pp].y > 0.0)  dist /= interp3_ym[IDX(i,j,k+1)]; 
      if( c[pp].x != 0.0 && c[pp].x < 0.0)  dist /= interp3_xp[IDX(i,j,k+1)]; 
      if( c[pp].x != 0.0 && c[pp].x > 0.0)  dist /= interp3_xm[IDX(i,j,k+1)]; 
  adv += fac*coeff_zp[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i,j,k+1)].p[pp] + dist*f[IDX(i+c[pp].x,j+c[pp].y,k+1)].p[pp] );
    }

    if(coeff_zm[IDX(i,j,k)].p[pp] >= 0.0){
      dist = interp_zm[IDX(i,j,k)];
      if( c[pp].y != 0.0 && c[pp].y < 0.0)  dist /= interp3_yp[IDX(i,j,k)]; 
      if( c[pp].y != 0.0 && c[pp].y > 0.0)  dist /= interp3_ym[IDX(i,j,k)]; 
      if( c[pp].x != 0.0 && c[pp].x < 0.0)  dist /= interp3_xp[IDX(i,j,k)]; 
      if( c[pp].x != 0.0 && c[pp].x > 0.0)  dist /= interp3_xm[IDX(i,j,k)]; 
  adv += fac*coeff_zm[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i,j,k)].p[pp] + dist*f[IDX(i+c[pp].x,j+c[pp].y,k)].p[pp] );
    }else{
     dist = interp_zp[IDX(i,j,k-1)];
      if( c[pp].y != 0.0 && c[pp].y < 0.0)  dist /= interp3_yp[IDX(i,j,k-1)]; 
      if( c[pp].y != 0.0 && c[pp].y > 0.0)  dist /= interp3_ym[IDX(i,j,k-1)]; 
      if( c[pp].x != 0.0 && c[pp].x < 0.0)  dist /= interp3_xp[IDX(i,j,k-1)]; 
      if( c[pp].x != 0.0 && c[pp].x > 0.0)  dist /= interp3_xm[IDX(i,j,k-1)]; 
  adv += fac*coeff_zm[IDX(i,j,k)].p[pp]*((1.0 - dist)*f[IDX(i,j,k-1)].p[pp] + dist*f[IDX(i+c[pp].x,j+c[pp].y,k-1)].p[pp] );
    }

  }/* end if c2 == 2*/
#endif

 return adv;
}	  
#endif

#ifdef METHOD_UPWIND_LINEAR
my_double compute_flux_with_upwind_linear(pop * f, int i, int j, int k, int pp){

  my_double adv=0.0;
  //fprintf(stderr,"interp %e %e\n",interp_zp[IDX(i,j,k)], interp2_zp[IDX(i,j,k)]);

#ifndef METHOD_UPWIND_LINEAR_IMPROVED
 /* linear upwind scheme */   
 if(coeff_xp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_xp[IDX(i,j,k)].p[pp]*(interp5_xp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_xp[IDX(i,j,k)])*f[IDX(i-1,j,k)].p[pp] );
 else
   adv += coeff_xp[IDX(i,j,k)].p[pp]*(interp6_xp[IDX(i,j,k)]*f[IDX(i+1,j,k)].p[pp] + (1.0 - interp6_xp[IDX(i,j,k)])*f[IDX(i+2,j,k)].p[pp] ); 

 if(coeff_xm[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_xm[IDX(i,j,k)].p[pp]*(interp5_xm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_xm[IDX(i,j,k)])*f[IDX(i+1,j,k)].p[pp] );
 else
   adv += coeff_xm[IDX(i,j,k)].p[pp]*(interp6_xm[IDX(i,j,k)]*f[IDX(i-1,j,k)].p[pp] + (1.0 - interp6_xm[IDX(i,j,k)])*f[IDX(i-2,j,k)].p[pp] ); 

 if(coeff_yp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_yp[IDX(i,j,k)].p[pp]*(interp5_yp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_yp[IDX(i,j,k)])*f[IDX(i,j-1,k)].p[pp] ); 
 else
   adv += coeff_yp[IDX(i,j,k)].p[pp]*(interp6_yp[IDX(i,j,k)]*f[IDX(i,j+1,k)].p[pp] + (1.0 - interp6_yp[IDX(i,j,k)])*f[IDX(i,j+2,k)].p[pp] );

 if(coeff_ym[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_ym[IDX(i,j,k)].p[pp]*(interp5_ym[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_ym[IDX(i,j,k)])*f[IDX(i,j+1,k)].p[pp] );
 else
   adv += coeff_ym[IDX(i,j,k)].p[pp]*(interp6_ym[IDX(i,j,k)]*f[IDX(i,j-1,k)].p[pp] + (1.0 - interp6_ym[IDX(i,j,k)])*f[IDX(i,j-2,k)].p[pp] );

 if(coeff_zp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_zp[IDX(i,j,k)].p[pp]*(interp5_zp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_zp[IDX(i,j,k)])*f[IDX(i,j,k-1)].p[pp] );
 else
   adv += coeff_zp[IDX(i,j,k)].p[pp]*(interp6_zp[IDX(i,j,k)]*f[IDX(i,j,k+1)].p[pp] + (1.0 - interp6_zp[IDX(i,j,k)])*f[IDX(i,j,k+2)].p[pp] );

 if(coeff_zm[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_zm[IDX(i,j,k)].p[pp]*(interp5_zm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_zm[IDX(i,j,k)])*f[IDX(i,j,k+1)].p[pp] );
 else
   adv += coeff_zm[IDX(i,j,k)].p[pp]*(interp6_zm[IDX(i,j,k)]*f[IDX(i,j,k-1)].p[pp] + (1.0 - interp6_zm[IDX(i,j,k)])*f[IDX(i,j,k-2)].p[pp] ); 

 //fprintf(stderr,"adv %e, %d %d %d\n",adv,i,j,k);
#else 
 /* IMPROVED LINEAR METHOD is defined */

if(coeff_xp[IDX(i,j,k)].p[pp] >= 0.0)
 adv += coeff_xp[IDX(i,j,k)].p[pp]*(f[IDX(i,j,k)].p[pp]+interp5_xp[IDX(i,j,k)]*(f[IDX(i+1,j,k)].p[pp]-f[IDX(i-1,j,k)].p[pp]));
 else
 adv += coeff_xp[IDX(i,j,k)].p[pp]*(f[IDX(i+1,j,k)].p[pp]+interp6_xp[IDX(i,j,k)]*(f[IDX(i,j,k)].p[pp]-f[IDX(i+2,j,k)].p[pp]));

 if(coeff_xm[IDX(i,j,k)].p[pp] >= 0.0)
 adv += coeff_xm[IDX(i,j,k)].p[pp]*(f[IDX(i,j,k)].p[pp]+interp5_xm[IDX(i,j,k)]*(f[IDX(i-1,j,k)].p[pp]-f[IDX(i+1,j,k)].p[pp]));
 else
 adv += coeff_xm[IDX(i,j,k)].p[pp]*(f[IDX(i-1,j,k)].p[pp]+interp6_xm[IDX(i,j,k)]*(f[IDX(i,j,k)].p[pp]-f[IDX(i-2,j,k)].p[pp]));

 if(coeff_yp[IDX(i,j,k)].p[pp] >= 0.0)
 adv += coeff_yp[IDX(i,j,k)].p[pp]*(f[IDX(i,j,k)].p[pp]+interp5_yp[IDX(i,j,k)]*(f[IDX(i,j+1,k)].p[pp]-f[IDX(i,j-1,k)].p[pp]));
 else
 adv += coeff_yp[IDX(i,j,k)].p[pp]*(f[IDX(i,j+1,k)].p[pp]+interp6_yp[IDX(i,j,k)]*(f[IDX(i,j,k)].p[pp]-f[IDX(i,j+2,k)].p[pp]));

 if(coeff_ym[IDX(i,j,k)].p[pp] >= 0.0)
 adv += coeff_ym[IDX(i,j,k)].p[pp]*(f[IDX(i,j,k)].p[pp]+interp5_ym[IDX(i,j,k)]*(f[IDX(i,j-1,k)].p[pp]-f[IDX(i,j+1,k)].p[pp]));
 else
 adv += coeff_ym[IDX(i,j,k)].p[pp]*(f[IDX(i,j-1,k)].p[pp]+interp6_ym[IDX(i,j,k)]*(f[IDX(i,j,k)].p[pp]-f[IDX(i,j-2,k)].p[pp]));

 if(coeff_zp[IDX(i,j,k)].p[pp] >= 0.0)
 adv += coeff_zp[IDX(i,j,k)].p[pp]*(f[IDX(i,j,k)].p[pp]+interp5_zp[IDX(i,j,k)]*(f[IDX(i,j,k+1)].p[pp]-f[IDX(i,j,k-1)].p[pp]));
 else
 adv += coeff_zp[IDX(i,j,k)].p[pp]*(f[IDX(i,j,k+1)].p[pp]+interp6_zp[IDX(i,j,k)]*(f[IDX(i,j,k)].p[pp]-f[IDX(i,j,k+2)].p[pp]));

 if(coeff_zm[IDX(i,j,k)].p[pp] >= 0.0)
 adv += coeff_zm[IDX(i,j,k)].p[pp]*(f[IDX(i,j,k)].p[pp]+interp5_zm[IDX(i,j,k)]*(f[IDX(i,j,k-1)].p[pp]-f[IDX(i,j,k+1)].p[pp]));
 else
 adv += coeff_zm[IDX(i,j,k)].p[pp]*(f[IDX(i,j,k-1)].p[pp]+interp6_zm[IDX(i,j,k)]*(f[IDX(i,j,k)].p[pp]-f[IDX(i,j,k-2)].p[pp]));
#endif

return adv;
}	  
#endif


#ifdef METHOD_MYQUICK
my_double compute_flux_with_quick(pop * f, int i, int j, int k, int pp){
 my_double aux;
 my_double adv=0.0;


/* QUICK for a regular grid */
/*
 my_double coeff_d, coeff_u, coeff_uu;
 coeff_d = 3.0/8.0;
 coeff_u = 6.0/8.0;
 coeff_uu = -1.0/8.0;
	  
 if(coeff_xp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_xp[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i+1,j,k)].p[pp] + coeff_u*f[IDX(i,j,k)].p[pp] + coeff_uu*f[IDX(i-1,j,k)].p[pp] );
 else
   adv += coeff_xp[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j,k)].p[pp] + coeff_u*f[IDX(i+1,j,k)].p[pp] + coeff_uu*f[IDX(i+2,j,k)].p[pp] );

 if(coeff_xm[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_xm[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i-1,j,k)].p[pp] + coeff_u*f[IDX(i,j,k)].p[pp] + coeff_uu*f[IDX(i+1,j,k)].p[pp] );
 else
   adv += coeff_xm[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j,k)].p[pp] + coeff_u*f[IDX(i-1,j,k)].p[pp] + coeff_uu*f[IDX(i-2,j,k)].p[pp] );

 if(coeff_yp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_yp[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j+1,k)].p[pp] + coeff_u*f[IDX(i,j,k)].p[pp] + coeff_uu*f[IDX(i,j-1,k)].p[pp] );
 else
   adv += coeff_yp[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j,k)].p[pp] + coeff_u*f[IDX(i,j+1,k)].p[pp] + coeff_uu*f[IDX(i,j+2,k)].p[pp] );

 if(coeff_ym[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_ym[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j-1,k)].p[pp] + coeff_u*f[IDX(i,j,k)].p[pp] + coeff_uu*f[IDX(i,j+1,k)].p[pp] );
 else
   adv += coeff_ym[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j,k)].p[pp] + coeff_u*f[IDX(i,j-1,k)].p[pp] + coeff_uu*f[IDX(i,j-2,k)].p[pp] );

 if(coeff_zp[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_zp[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j,k+1)].p[pp] + coeff_u*f[IDX(i,j,k)].p[pp] + coeff_uu*f[IDX(i,j,k-1)].p[pp] );
 else
   adv += coeff_zp[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j,k)].p[pp] + coeff_u*f[IDX(i,j,k+1)].p[pp] + coeff_uu*f[IDX(i,j,k+2)].p[pp] );

 if(coeff_zm[IDX(i,j,k)].p[pp] >= 0.0)
   adv += coeff_zm[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j,k-1)].p[pp] + coeff_u*f[IDX(i,j,k)].p[pp] + coeff_uu*f[IDX(i,j,k+1)].p[pp] );
 else
   adv += coeff_zm[IDX(i,j,k)].p[pp]*( coeff_d*f[IDX(i,j,k)].p[pp] + coeff_u*f[IDX(i,j,k-1)].p[pp] + coeff_uu*f[IDX(i,j,k-2)].p[pp] );	  
 */

#ifdef METHOD_MYQUICK_CARTESIAN
 /* when the grid is cartesian rectangular, the quick algorithm can be written in a more compact form */ 
 //if(pp>0){
 if(coeff_xp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_xp[IDX(i,j,k)].p[pp] > 0.0){
   aux= 1.0 + interp2_xp[IDX(i,j,k)] - interp3_xm[IDX(i,j,k)];
   adv += coeff_xp[IDX(i,j,k)].p[pp]*( 
				      interp_xp[IDX(i,j,k)]*f[IDX(i+1,j,k)].p[pp] 
				      + (aux - interp_xp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] 
				      - (aux + interp4_xm[IDX(i,j,k)])*f[IDX(i-1,j,k)].p[pp] 
				      + interp4_xm[IDX(i,j,k)]*f[IDX(i-2,j,k)].p[pp] 
				       );
   }else{
   aux = 1.0 + interp2_xm[IDX(i,j,k)] - interp3_xp[IDX(i,j,k)];
   adv += coeff_xp[IDX(i,j,k)].p[pp]*( 
				      (aux + interp4_xp[IDX(i,j,k)])*f[IDX(i+1,j,k)].p[pp] 
				      - interp4_xp[IDX(i,j,k)]*f[IDX(i+2,j,k)].p[pp] 
				      - interp_xm[IDX(i,j,k)]*f[IDX(i-1,j,k)].p[pp] 
				      - (aux - interp_xm[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] 
				       );
   }
 }


 if(coeff_yp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_yp[IDX(i,j,k)].p[pp] > 0.0){
   aux = 1.0  + interp2_yp[IDX(i,j,k)] - interp3_ym[IDX(i,j,k)];
   adv += coeff_yp[IDX(i,j,k)].p[pp]*( 
				      interp_yp[IDX(i,j,k)]*f[IDX(i,j+1,k)].p[pp] 
				      + (aux - interp_yp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] 
				      - (aux + interp4_ym[IDX(i,j,k)])*f[IDX(i,j-1,k)].p[pp] 
				      + interp4_ym[IDX(i,j,k)]*f[IDX(i,j-2,k)].p[pp] 
				       );
   }else{
   aux = 1.0 + interp2_ym[IDX(i,j,k)] - interp3_yp[IDX(i,j,k)];
   adv += coeff_yp[IDX(i,j,k)].p[pp]*( 
				      (aux + interp4_yp[IDX(i,j,k)])*f[IDX(i,j+1,k)].p[pp] 
				      - interp4_yp[IDX(i,j,k)]*f[IDX(i,j+2,k)].p[pp] 
				      - interp_ym[IDX(i,j,k)]*f[IDX(i,j-1,k)].p[pp] 
				      - (aux - interp_ym[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] 
				       );
   }
 }


 if(coeff_zp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_zp[IDX(i,j,k)].p[pp] > 0.0){
   aux = 1.0 + interp2_zp[IDX(i,j,k)] - interp3_zm[IDX(i,j,k)];
   adv += coeff_zp[IDX(i,j,k)].p[pp]*( 
				      interp_zp[IDX(i,j,k)]*f[IDX(i,j,k+1)].p[pp] 
				      + (aux - interp_zp[IDX(i,j,k)] )*f[IDX(i,j,k)].p[pp] 
				      - (aux + interp4_zm[IDX(i,j,k)])*f[IDX(i,j,k-1)].p[pp] 
				      + interp4_zm[IDX(i,j,k)]*f[IDX(i,j,k-2)].p[pp]
				       );
   }else{
   aux = 1.0 + interp2_zm[IDX(i,j,k)] - interp3_zp[IDX(i,j,k)];
   adv += coeff_zp[IDX(i,j,k)].p[pp]*( 
				      (aux + interp4_zp[IDX(i,j,k)])*f[IDX(i,j,k+1)].p[pp] 
				      - interp4_zp[IDX(i,j,k)]*f[IDX(i,j,k+2)].p[pp] 
				      -interp_zm[IDX(i,j,k)]*f[IDX(i,j,k-1)].p[pp] 
				      - (aux - interp_zm[IDX(i,j,k)] )*f[IDX(i,j,k)].p[pp] 
				       );
   }
 }
 //}/*end of if pp > 0 */
#else

 /* The good old one quick, less compact but well tested */
	  //if(pp>0){
 if(coeff_xp[IDX(i,j,k)].p[pp] != 0.0){
 if(coeff_xp[IDX(i,j,k)].p[pp] > 0.0)
   adv += coeff_xp[IDX(i,j,k)].p[pp]*( interp_xp[IDX(i,j,k)]*f[IDX(i+1,j,k)].p[pp] + (1.0 - interp_xp[IDX(i,j,k)] + interp2_xp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_xp[IDX(i,j,k)]*f[IDX(i-1,j,k)].p[pp] );
 else
   adv += coeff_xp[IDX(i,j,k)].p[pp]*( interp3_xp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_xp[IDX(i,j,k)] + interp4_xp[IDX(i,j,k)])*f[IDX(i+1,j,k)].p[pp] - interp4_xp[IDX(i,j,k)]*f[IDX(i+2,j,k)].p[pp] );
 }

 if(coeff_xm[IDX(i,j,k)].p[pp] != 0.0){
 if(coeff_xm[IDX(i,j,k)].p[pp] > 0.0)
   adv += coeff_xm[IDX(i,j,k)].p[pp]*( interp_xm[IDX(i,j,k)]*f[IDX(i-1,j,k)].p[pp] + (1.0 - interp_xm[IDX(i,j,k)] + interp2_xm[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_xm[IDX(i,j,k)]*f[IDX(i+1,j,k)].p[pp] );
 else
   adv += coeff_xm[IDX(i,j,k)].p[pp]*( interp3_xm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_xm[IDX(i,j,k)] + interp4_xm[IDX(i,j,k)])*f[IDX(i-1,j,k)].p[pp] - interp4_xm[IDX(i,j,k)] *f[IDX(i-2,j,k)].p[pp] );
 }

 if(coeff_yp[IDX(i,j,k)].p[pp] != 0.0){
 if(coeff_yp[IDX(i,j,k)].p[pp] > 0.0)
   adv += coeff_yp[IDX(i,j,k)].p[pp]*( interp_yp[IDX(i,j,k)]*f[IDX(i,j+1,k)].p[pp] + (1.0 - interp_yp[IDX(i,j,k)] + interp2_yp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_yp[IDX(i,j,k)]*f[IDX(i,j-1,k)].p[pp] );
 else
   adv += coeff_yp[IDX(i,j,k)].p[pp]*( interp3_yp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_yp[IDX(i,j,k)] + interp4_yp[IDX(i,j,k)])*f[IDX(i,j+1,k)].p[pp] - interp4_yp[IDX(i,j,k)]*f[IDX(i,j+2,k)].p[pp] );
 }

 if(coeff_ym[IDX(i,j,k)].p[pp] != 0.0){
 if(coeff_ym[IDX(i,j,k)].p[pp] > 0.0)
   adv += coeff_ym[IDX(i,j,k)].p[pp]*( interp_ym[IDX(i,j,k)]*f[IDX(i,j-1,k)].p[pp] + (1.0 - interp_ym[IDX(i,j,k)] + interp2_ym[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_ym[IDX(i,j,k)]*f[IDX(i,j+1,k)].p[pp] );
 else
   adv += coeff_ym[IDX(i,j,k)].p[pp]*( interp3_ym[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_ym[IDX(i,j,k)] + interp4_ym[IDX(i,j,k)])*f[IDX(i,j-1,k)].p[pp] - interp4_ym[IDX(i,j,k)]*f[IDX(i,j-2,k)].p[pp] );
 }

 if(coeff_zp[IDX(i,j,k)].p[pp] != 0.0){
 if(coeff_zp[IDX(i,j,k)].p[pp] > 0.0)
   adv += coeff_zp[IDX(i,j,k)].p[pp]*( interp_zp[IDX(i,j,k)]*f[IDX(i,j,k+1)].p[pp] + (1.0 - interp_zp[IDX(i,j,k)] + interp2_zp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_zp[IDX(i,j,k)]*f[IDX(i,j,k-1)].p[pp] );
 else
   adv += coeff_zp[IDX(i,j,k)].p[pp]*( interp3_zp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_zp[IDX(i,j,k)] + interp4_zp[IDX(i,j,k)])*f[IDX(i,j,k+1)].p[pp] - interp4_zp[IDX(i,j,k)]*f[IDX(i,j,k+2)].p[pp] );
 }

 if(coeff_zm[IDX(i,j,k)].p[pp] != 0.0){
 if(coeff_zm[IDX(i,j,k)].p[pp] > 0.0)
   adv += coeff_zm[IDX(i,j,k)].p[pp]*( interp_zm[IDX(i,j,k)]*f[IDX(i,j,k-1)].p[pp] + (1.0 - interp_zm[IDX(i,j,k)] + interp2_zm[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_zm[IDX(i,j,k)]*f[IDX(i,j,k+1)].p[pp] );
 else
   adv += coeff_zm[IDX(i,j,k)].p[pp]*( interp3_zm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_zm[IDX(i,j,k)] + interp4_zm[IDX(i,j,k)])*f[IDX(i,j,k-1)].p[pp] -interp4_zm[IDX(i,j,k)]*f[IDX(i,j,k-2)].p[pp] );	  
   }
 //}/*end of if pp > 0 */
#endif

 return adv;
}
#endif


/* flux limiters */
#ifdef METHOD_MYQUICK
#ifdef METHOD_MYQUICK_LIMITER
/* see http://en.wikipedia.org/wiki/Flux_limiter */
my_double limiter(my_double r){

  my_double phi,b,r3,r1,r2;
  /*
  if(r2==0.0 || r1*r2 < 0.0){  
    r=0.0;
  }else{
    r=r1/r2;
  }  
  */
  //if(r<=0)phi=0;else phi=1.0;
  //phi=0.0;

  /* mine TVD */ 
  b=1.035; // poiseuille value 
  //b=1.14; // kolmogorov value
  /*
  if(r<0.5){ phi=2.*r; }else{
    if(r<1.0){ phi = 1.0; }else{
      if(r<b) phi=r; else phi=b;
    }
  }
  */
  /* osher */
  /*
  b=3.0;
  r3=MIN(r,b);
  phi=MAX(0,r3);
  */

  /* sweeby symmetric*/
  
  b=1.13; //poiseulle
  //b = 1.22; //kolmogorov
  /*
  r1=MIN(b*r,1.0);
  r2=MIN(r,b);
  r3=MAX(r1,r2); 
  phi=MAX(0,r3);
  */

  /* charm */
  //if(r>0) phi = r*(3.*r+1.)/pow((r+1.0),2.0); else phi=0.0;

  /* HQUICK */
  //phi=2.0*(r+fabs(r))/(r+3.0);

  /* smart */
  /*
  r1=2.*r;
  r2=(0.25+0.75*r);
  r3=MIN(r1,r2);
  r1=MIN(r3,4.0);
  phi=MAX(0.0,r1); 
  */

  /* superbee */
  /*     
  r1 = MIN(2.0*r, 1.0);
  r2 = MIN(r, 2.0);
  r3 = MAX(r1,r2);
  phi = MAX(0.,r3);
  */

  /* minmod */
  ///*    
  r1 = MIN(1,r);
  phi = MAX(0,r1);
  //*/

  /* van Leer */  
  //phi = (r+fabs(r))/(1.0+fabs(r));
  //phi = 2.*r/(r*r + 1.0);
  //phi = (r*r + r)/(r*r + 1.0);
  // phi=0.1;
  //fprintf(stderr,"r=%e\n",r);
  //phi=1.2;
  return phi;  
}


my_double gradient_ratio(pop * f, int i, int j, int k, int pp,int dir){
  my_double r,r1,r2;

  if(dir==0){
    r1 =  (f[IDX(i , j,k)].p[pp] - f[IDX(i-1,j,k)].p[pp])/(center_V[IDX(i  ,j,k)].x - center_V[IDX(i-1,j,k)].x);
    r2 =  (f[IDX(i+1,j,k)].p[pp] - f[IDX(i  ,j,k)].p[pp])/(center_V[IDX(i+1,j,k)].x - center_V[IDX(i  ,j,k)].x);
  }
  if(dir==1){
    r1 =  (f[IDX(i , j,k)].p[pp] - f[IDX(i,j-1,k)].p[pp])/(center_V[IDX(i  ,j,k)].y - center_V[IDX(i,j-1,k)].y);
    r2 =  (f[IDX(i,j+1,k)].p[pp] - f[IDX(i  ,j,k)].p[pp])/(center_V[IDX(i,j+1,k)].y - center_V[IDX(i  ,j,k)].y);
  }
  if(dir==2){
    r1 =  (f[IDX(i , j,k)].p[pp] - f[IDX(i,j,k-1)].p[pp])/(center_V[IDX(i  ,j,k)].z - center_V[IDX(i,j,k-1)].z); 
    r2 =  (f[IDX(i,j,k+1)].p[pp] - f[IDX(i  ,j,k)].p[pp])/(center_V[IDX(i,j,k+1)].z - center_V[IDX(i  ,j,k)].z);
  }
     if(r2!=0) r=r1/r2; else r=0.0;
     return r;
}

/************ Advection with flux limiter ****************/

my_double compute_flux_with_limiters(pop * f, int i, int j, int k, int pp){

  my_double advL,advH,adv,r,r1,r2;

  adv=0.0;

/* x */

 if(coeff_xp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_xp[IDX(i,j,k)].p[pp] > 0.0){
   r=gradient_ratio(f,i,j,k,pp,0);
   advL = coeff_xp[IDX(i,j,k)].p[pp]*(interp5_xp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_xp[IDX(i,j,k)])*f[IDX(i-1,j,k)].p[pp] );
   advH = coeff_xp[IDX(i,j,k)].p[pp]*( interp_xp[IDX(i,j,k)]*f[IDX(i+1,j,k)].p[pp] + (1.0 - interp_xp[IDX(i,j,k)] + interp2_xp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_xp[IDX(i,j,k)]*f[IDX(i-1,j,k)].p[pp] );
   }else{
    r=gradient_ratio(f,i,j,k,pp,0); 
   advL = coeff_xp[IDX(i,j,k)].p[pp]*(interp6_xp[IDX(i,j,k)]*f[IDX(i+1,j,k)].p[pp] + (1.0 - interp6_xp[IDX(i,j,k)])*f[IDX(i+2,j,k)].p[pp] ); 
   advH = coeff_xp[IDX(i,j,k)].p[pp]*( interp3_xp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_xp[IDX(i,j,k)] + interp4_xp[IDX(i,j,k)])*f[IDX(i+1,j,k)].p[pp] - interp4_xp[IDX(i,j,k)]*f[IDX(i+2,j,k)].p[pp] );
   }
 
   adv += advL - limiter(r)*(advL - advH);
}



 if(coeff_xm[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_xm[IDX(i,j,k)].p[pp] > 0.0){
   r=gradient_ratio(f,i-1,j,k,pp,0);
   advL = coeff_xm[IDX(i,j,k)].p[pp]*(interp5_xm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_xm[IDX(i,j,k)])*f[IDX(i+1,j,k)].p[pp] );
   advH = coeff_xm[IDX(i,j,k)].p[pp]*( interp_xm[IDX(i,j,k)]*f[IDX(i-1,j,k)].p[pp] + (1.0 - interp_xm[IDX(i,j,k)] + interp2_xm[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_xm[IDX(i,j,k)]*f[IDX(i+1,j,k)].p[pp] );
   }else{
   r=gradient_ratio(f,i-1,j,k,pp,0);
   advL = coeff_xm[IDX(i,j,k)].p[pp]*(interp6_xm[IDX(i,j,k)]*f[IDX(i-1,j,k)].p[pp] + (1.0 - interp6_xm[IDX(i,j,k)])*f[IDX(i-2,j,k)].p[pp] ); 
   advH = coeff_xm[IDX(i,j,k)].p[pp]*( interp3_xm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_xm[IDX(i,j,k)] + interp4_xm[IDX(i,j,k)])*f[IDX(i-1,j,k)].p[pp] - interp4_xm[IDX(i,j,k)] *f[IDX(i-2,j,k)].p[pp] );
   }
 
   adv += advL - limiter(r)*(advL - advH);
 }


/* y */

 if(coeff_yp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_yp[IDX(i,j,k)].p[pp] > 0.0){
   r=gradient_ratio(f,i,j,k,pp,1);
   advL = coeff_yp[IDX(i,j,k)].p[pp]*(interp5_yp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_yp[IDX(i,j,k)])*f[IDX(i,j-1,k)].p[pp] ); 
   advH = coeff_yp[IDX(i,j,k)].p[pp]*( interp_yp[IDX(i,j,k)]*f[IDX(i,j+1,k)].p[pp] + (1.0 - interp_yp[IDX(i,j,k)] + interp2_yp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_yp[IDX(i,j,k)]*f[IDX(i,j-1,k)].p[pp] );
   }else{
   r=gradient_ratio(f,i,j,k,pp,1);
  advL = coeff_yp[IDX(i,j,k)].p[pp]*(interp6_yp[IDX(i,j,k)]*f[IDX(i,j+1,k)].p[pp] + (1.0 - interp6_yp[IDX(i,j,k)])*f[IDX(i,j+2,k)].p[pp] );
   advH = coeff_yp[IDX(i,j,k)].p[pp]*( interp3_yp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_yp[IDX(i,j,k)] + interp4_yp[IDX(i,j,k)])*f[IDX(i,j+1,k)].p[pp] - interp4_yp[IDX(i,j,k)]*f[IDX(i,j+2,k)].p[pp] );
   }
 
   adv += advL - limiter(r)*(advL - advH);
 }


 if(coeff_ym[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_ym[IDX(i,j,k)].p[pp] > 0.0){
   r=gradient_ratio(f,i,j-1,k,pp,1);
   advL = coeff_ym[IDX(i,j,k)].p[pp]*(interp5_ym[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_ym[IDX(i,j,k)])*f[IDX(i,j+1,k)].p[pp] );
   advH = coeff_ym[IDX(i,j,k)].p[pp]*( interp_ym[IDX(i,j,k)]*f[IDX(i,j-1,k)].p[pp] + (1.0 - interp_ym[IDX(i,j,k)] + interp2_ym[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_ym[IDX(i,j,k)]*f[IDX(i,j+1,k)].p[pp] );
   }else{
   r=gradient_ratio(f,i,j-1,k,pp,1);
   advL = coeff_ym[IDX(i,j,k)].p[pp]*(interp6_ym[IDX(i,j,k)]*f[IDX(i,j-1,k)].p[pp] + (1.0 - interp6_ym[IDX(i,j,k)])*f[IDX(i,j-2,k)].p[pp] );
   advH = coeff_ym[IDX(i,j,k)].p[pp]*( interp3_ym[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_ym[IDX(i,j,k)] + interp4_ym[IDX(i,j,k)])*f[IDX(i,j-1,k)].p[pp] - interp4_ym[IDX(i,j,k)]*f[IDX(i,j-2,k)].p[pp] );
   }
 
   adv += advL - limiter(r)*(advL - advH);
 }



/* z */

 if(coeff_zp[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_zp[IDX(i,j,k)].p[pp] > 0.0){
   r=gradient_ratio(f,i,j,k,pp,2);
   advL = coeff_zp[IDX(i,j,k)].p[pp]*(interp5_zp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_zp[IDX(i,j,k)])*f[IDX(i,j,k-1)].p[pp] );
   advH = coeff_zp[IDX(i,j,k)].p[pp]*( interp_zp[IDX(i,j,k)]*f[IDX(i,j,k+1)].p[pp] + (1.0 - interp_zp[IDX(i,j,k)] + interp2_zp[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_zp[IDX(i,j,k)]*f[IDX(i,j,k-1)].p[pp] );
   }else{
   r=gradient_ratio(f,i,j,k,pp,2);
  advL = coeff_zp[IDX(i,j,k)].p[pp]*(interp6_zp[IDX(i,j,k)]*f[IDX(i,j,k+1)].p[pp] + (1.0 - interp6_zp[IDX(i,j,k)])*f[IDX(i,j,k+2)].p[pp] );
  advH = coeff_zp[IDX(i,j,k)].p[pp]*( interp3_zp[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_zp[IDX(i,j,k)] + interp4_zp[IDX(i,j,k)])*f[IDX(i,j,k+1)].p[pp] - interp4_zp[IDX(i,j,k)]*f[IDX(i,j,k+2)].p[pp] );
   }
 
   adv += advL - limiter(r)*(advL - advH);
 }


 if(coeff_zm[IDX(i,j,k)].p[pp] != 0.0){
   if(coeff_zm[IDX(i,j,k)].p[pp] > 0.0){
   r=gradient_ratio(f,i,j,k-1,pp,2);
   advL = coeff_zm[IDX(i,j,k)].p[pp]*(interp5_zm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp5_zm[IDX(i,j,k)])*f[IDX(i,j,k+1)].p[pp] );
   advH = coeff_zm[IDX(i,j,k)].p[pp]*( interp_zm[IDX(i,j,k)]*f[IDX(i,j,k-1)].p[pp] + (1.0 - interp_zm[IDX(i,j,k)] + interp2_zm[IDX(i,j,k)])*f[IDX(i,j,k)].p[pp] - interp2_zm[IDX(i,j,k)]*f[IDX(i,j,k+1)].p[pp] );
   }else{
  r=gradient_ratio(f,i,j,k-1,pp,2);
   advL = coeff_zm[IDX(i,j,k)].p[pp]*(interp6_zm[IDX(i,j,k)]*f[IDX(i,j,k-1)].p[pp] + (1.0 - interp6_zm[IDX(i,j,k)])*f[IDX(i,j,k-2)].p[pp] ); 
   advH = coeff_zm[IDX(i,j,k)].p[pp]*( interp3_zm[IDX(i,j,k)]*f[IDX(i,j,k)].p[pp] + (1.0 - interp3_zm[IDX(i,j,k)] + interp4_zm[IDX(i,j,k)])*f[IDX(i,j,k-1)].p[pp] -interp4_zm[IDX(i,j,k)]*f[IDX(i,j,k-2)].p[pp] );
   }
 
   adv += advL - limiter(r)*(advL - advH);
 }


 return adv;
}
#endif
#endif


/*  Here Collision */

void add_collision(pop *f, pop *rhs_f, my_double tau,pop *f_eq,char which_pop){
  int i, j, k, pp;
  my_double invtau;
  pop ff_eq;
  pop f_eq_xp,f_eq_xm,f_eq_yp,f_eq_ym,f_eq_zp,f_eq_zm;
  pop fcoll, fcoll_xp,fcoll_xm,fcoll_yp,fcoll_ym,fcoll_zp,fcoll_zm;
  my_double fac;
  my_double tau_les;
  vector rho_dt_u , u_dt_rho;
  my_double wgt2;
#ifdef METHOD_TRT
 /* This is Two Relaxation Time (TRT) */
 my_double magic_gamma;
 my_double invtau_minus;
 pop ff_eq_plus, ff_eq_minus, ff_plus, ff_minus;
#endif

 /* Inverse of the relaxation time */ 
 invtau = 1.0/tau;


#ifndef METHOD_COLLISION_IMPLICIT

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

/* Here we begin the LES part , we compute tau for every i,j,k*/
#ifdef LB_FLUID_LES
  if( which_pop == 'p' ){
    tau_les = tau_u_les(i,j,k);
    invtau = 1.0/ tau_les;
    // fprintf(stderr,"%e %e\n",tau,tau_les);
  }
#endif
#ifdef LB_TEMPERATURE_LES
  if( which_pop == 'g' ){
    tau_les = tau_t_les(i,j,k);
    invtau = 1.0/tau_les;
  }
#endif
#ifdef LB_SCALAR_LES
  if( which_pop == 'h' ){
    tau_les = tau_s_les(i,j,k);
    invtau = 1.0/tau_les;
  }
#endif
/* end of LES part */


/* Melting */
#ifdef  LB_TEMPERATURE_MELTING
 #ifdef  LB_TEMPERATURE_MELTING_SOLID_DIFFUSIVITY
  /* Here we implement the variable thermal diffusivity term */
  /* Following "A one-dimensional enthalpy model of sea ice" by Dirk NOTZ, M. Grae WORSTER on Annals of Glaciology 44 (2006) */
  /* see equation (3) : k_effective = lf * k_liquid + (1-lf) * k_solid */
  if( which_pop == 'g' ){
      invtau = 1./(liquid_frac[IDX(i,j,k)]*property.tau_t  + (1.-liquid_frac[IDX(i,j,k)])*property.tau_solid);
  }
 #endif
#endif


#ifdef METHOD_REDEFINED_POP
	ff_eq = f_eq[IDX(i,j,k)];
#else     
	ff_eq=equilibrium(f,i,j,k);
#endif


#ifndef METHOD_EXPONENTIAL
  #ifdef METHOD_TRT
	/* This is two relaxation time TRT */
	magic_gamma = 0.25; //3./16.;
	/* see http://arxiv.org/abs/1508.07982v1 */
	invtau_minus = (4. - 2.*invtau)/(2.+(4.*magic_gamma -1.)*invtau);
	for (pp=0; pp<NPOP; pp++){
	    ff_eq_plus.p[pp]  = 0.5*(ff_eq.p[pp] + ff_eq.p[inv[pp]]);
	    ff_eq_minus.p[pp] = 0.5*(ff_eq.p[pp] - ff_eq.p[inv[pp]]);
	    ff_plus.p[pp]  = 0.5*(f[IDX(i,j,k)].p[pp] + f[IDX(i,j,k)].p[inv[pp]]);
	    ff_minus.p[pp] = 0.5*(f[IDX(i,j,k)].p[pp] - f[IDX(i,j,k)].p[inv[pp]]);		
	fcoll.p[pp] = -invtau * (ff_plus.p[pp] - ff_eq_plus.p[pp]) -invtau_minus * (ff_minus.p[pp] - ff_eq_minus.p[pp]);
      }
  #else
	/* This is single relaxation time SRT */
       	for (pp=0; pp<NPOP; pp++) fcoll.p[pp] = -invtau * (f[IDX(i,j,k)].p[pp] - ff_eq.p[pp]);
  #endif

#ifdef METHOD_LOG
	for (pp=0; pp<NPOP; pp++) fcoll.p[pp] =  ( exp(invtau*(ff_eq.p[pp] - f[IDX(i,j,k)].p[pp]) ) - 1.0);
#endif

#endif
     
#if (defined METHOD_MYQUICK && defined METHOD_TRAPEZOID)
	f_eq_xp = equilibrium(f,i+1,j,k); 
        f_eq_xm = equilibrium(f,i-1,j,k);
	f_eq_yp = equilibrium(f,i,j+1,k); 
        f_eq_ym = equilibrium(f,i,j-1,k);
	f_eq_zp = equilibrium(f,i,j,k+1); 
        f_eq_zm = equilibrium(f,i,j,k-1);

	for (pp=0; pp<NPOP; pp++){
	  fcoll_xp.p[pp] = -invtau * (f[IDX(i+1,j,k)].p[pp] - f_eq_xp.p[pp]);
	  fcoll_xm.p[pp] = -invtau * (f[IDX(i-1,j,k)].p[pp] - f_eq_xm.p[pp]);
	  fcoll_yp.p[pp] = -invtau * (f[IDX(i,j+1,k)].p[pp] - f_eq_yp.p[pp]);
	  fcoll_ym.p[pp] = -invtau * (f[IDX(i,j-1,k)].p[pp] - f_eq_ym.p[pp]);
	  fcoll_zp.p[pp] = -invtau * (f[IDX(i,j,k+1)].p[pp] - f_eq_zp.p[pp]);
	  fcoll_zm.p[pp] = -invtau * (f[IDX(i,j,k-1)].p[pp] - f_eq_zm.p[pp]);
	}
#endif
	
	for (pp=0; pp<NPOP; pp++){

	  //#define ONLY_COLLISION
#ifdef ONLY_COLLISION
	  /* set to zero , just for check to eliminate advection */
	  rhs_f[IDX(i,j,k)].p[pp] = 0.0;
#endif

	/* collision */
#ifdef METHOD_EXPONENTIAL
	  	rhs_f[IDX(i,j,k)].p[pp] +=   invtau * ff_eq.p[pp];
#else

		
#if (defined METHOD_MYQUICK && defined METHOD_TRAPEZOID)
   rhs_f[IDX(i,j,k)].p[pp] +=  0.25*fcoll.p[pp];

   rhs_f[IDX(i,j,k)].p[pp] += 0.125*( interp_xp[IDX(i,j,k)]*fcoll_xp.p[pp] + (1.0 - interp_xp[IDX(i,j,k)] + interp2_xp[IDX(i,j,k)])*fcoll.p[pp] - interp2_xp[IDX(i,j,k)]*fcoll_xm.p[pp] );
   rhs_f[IDX(i,j,k)].p[pp] += 0.125*( interp_xm[IDX(i,j,k)]*fcoll_xm.p[pp] + (1.0 - interp_xm[IDX(i,j,k)] + interp2_xm[IDX(i,j,k)])*fcoll.p[pp] - interp2_xm[IDX(i,j,k)]*fcoll_xp.p[pp] );

   rhs_f[IDX(i,j,k)].p[pp] += 0.125*( interp_yp[IDX(i,j,k)]*fcoll_yp.p[pp] + (1.0 - interp_yp[IDX(i,j,k)] + interp2_yp[IDX(i,j,k)])*fcoll.p[pp] - interp2_yp[IDX(i,j,k)]*fcoll_ym.p[pp] );
   rhs_f[IDX(i,j,k)].p[pp] += 0.125*( interp_ym[IDX(i,j,k)]*fcoll_ym.p[pp] + (1.0 - interp_ym[IDX(i,j,k)] + interp2_ym[IDX(i,j,k)])*fcoll.p[pp] - interp2_ym[IDX(i,j,k)]*fcoll_yp.p[pp] );

   rhs_f[IDX(i,j,k)].p[pp] += 0.125*( interp_zp[IDX(i,j,k)]*fcoll_zp.p[pp] + (1.0 - interp_zp[IDX(i,j,k)] + interp2_zp[IDX(i,j,k)])*fcoll.p[pp] - interp2_zp[IDX(i,j,k)]*fcoll_zm.p[pp] );
   rhs_f[IDX(i,j,k)].p[pp] += 0.125*( interp_zm[IDX(i,j,k)]*fcoll_zm.p[pp] + (1.0 - interp_zm[IDX(i,j,k)] + interp2_zm[IDX(i,j,k)])*fcoll.p[pp] - interp2_zm[IDX(i,j,k)]*fcoll_zp.p[pp] );

#else
   //rhs_f[IDX(i,j,k)].p[pp] +=  -invtau * (f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
			rhs_f[IDX(i,j,k)].p[pp] +=  fcoll.p[pp];
#endif
		
  //rhs_f[IDX(i,j,k)].p[pp] +=  -invtau * (f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
#endif

	}/* pp */
      }/* i,j,k */
#endif

 #ifdef METHOD_CORRECTION_LATT
     // Jonas LATT correction in advection-diffusion equation (see his PhD thesis)
  #ifdef LB_FLUID_PAST      
  if( which_pop == 'g' ||  which_pop == 'h'){
    //fprintf(stderr,"AAA\n");
  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){
      rho_dt_u.x = dens[IDX(i, j, k)]*(u[IDX(i, j, k)].x  - old_u[IDX(i, j, k)].x)*invtau;
      rho_dt_u.y = dens[IDX(i, j, k)]*(u[IDX(i, j, k)].y  - old_u[IDX(i, j, k)].y)*invtau;
      rho_dt_u.z = dens[IDX(i, j, k)]*(u[IDX(i, j, k)].z  - old_u[IDX(i, j, k)].z)*invtau;
      u_dt_rho.x = u[IDX(i, j, k)].x*(dens[IDX(i, j, k)] - old_dens[IDX(i, j, k)])*invtau;
      u_dt_rho.y = u[IDX(i, j, k)].y*(dens[IDX(i, j, k)] - old_dens[IDX(i, j, k)])*invtau;
      u_dt_rho.z = u[IDX(i, j, k)].z*(dens[IDX(i, j, k)] - old_dens[IDX(i, j, k)])*invtau;

  for (pp=0; pp<NPOP; pp++){
      wgt2=c[pp].x*(rho_dt_u.x-u_dt_rho.x) + c[pp].y*(rho_dt_u.y-u_dt_rho.y) + c[pp].z*(rho_dt_u.z-u_dt_rho.z);
       if( which_pop == 'g' ) rhs_f[IDX(i,j,k)].p[pp] += wgt[pp]*(1.0-0.5*property.time_dt*invtau)*invcs2*wgt2;
       if( which_pop == 'h' ) rhs_f[IDX(i,j,k)].p[pp] += wgt[pp]*(1.0-0.5*property.time_dt*invtau)*invcs2*wgt2;
	}/* pp */
      }/* i,j,k */
}/* end of if on g or h */
  #endif
 #endif

}/* end of add_collision */



