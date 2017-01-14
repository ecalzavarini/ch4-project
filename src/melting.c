#include "common_object.h"


#ifdef LB_TEMPERATURE_MELTING
void melting(){
  int i,j,k;
  my_double fl , Tl;
  my_double hs, hl;
  my_double Ts,Cp,Lf,enthalpy;
  my_double fac1, fac2;
  my_double eps = 1.0;
  my_double L,y;

  Ts = (my_double) property.T_solid;
  Cp = (my_double) property.specific_heat;
  Lf = (my_double) property.latent_heat;
  /* solidification enthalpy */
  hs = Cp*Ts;
  hl = hs+Lf;
  Tl = hl/Cp;

  fac1 = (Lf/Cp)/property.time_dt;


   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){

	/* store previous fluid fraction */
	liquid_frac_old[IDX(i,j,k)] = liquid_frac[IDX(i,j,k)];

	/* here we can define a SOLIDUS (or LIQUIDUS) slope */

 #ifdef LB_TEMPERATURE_MELTING_SOLIDUS
   #ifdef LB_TEMPERATURE_MELTING_SOLIDUS_LINEAR
        /* linear gradient */
     L=(my_double)property.SY; 
     y = (my_double)center_V[IDX(i,j,k)].y;
     Ts = ( (property.T_bot-property.T_ref) - (property.deltaT/L)*y );
   #endif
     hs = Cp*Ts;
     hl = hs+Lf;
 #endif

      /* subtract the freezing point depression */
      /*	
      #ifdef LB_SCALAR
      Ts = property.liquidus_slope * s[IDX(i,j,k)];
      hs = Cp*Ts;
      hl = hs+Lf;
      #endif
      */

      /* compute Entalphy */
	enthalpy = Cp*t[IDX(i,j,k)] + Lf*liquid_frac[IDX(i,j,k)];

      /* compute new fluid fraction */
      if(enthalpy < hs) liquid_frac[IDX(i,j,k)]=0.0;
      else if(enthalpy > hl) liquid_frac[IDX(i,j,k)]=1.0;
      else liquid_frac[IDX(i,j,k)]= (enthalpy - hs)/(hl - hs);

      /* add melting term to the temperature field */
      t_source[IDX(i,j,k)] -= fac1*(liquid_frac[IDX(i,j,k)]-liquid_frac_old[IDX(i,j,k)]);
    

#ifndef LB_TEMPERATURE_MELTING_BOUNCEBACK
      /* If melting bounce-back conditions are not defined */
      /* add drag term to the force structure */
      fl = liquid_frac[IDX(i,j,k)]; 

      /* Here we weight the force with the local liquid fraction ( so there will be zero momentum force in the solid ) */
      force[IDX(i,j,k)].x *= fl;
      force[IDX(i,j,k)].y *= fl;
      force[IDX(i,j,k)].z *= fl;

      /* Here different penalization schemes */
      //fac2 = (1.0 - fl*fl)/(eps + fl*fl*fl);
      fac2 = (1.0 - fl*fl);
      //if(fl<0.5) fac2 = 1.0; else fac2 = 0.0;
      force[IDX(i,j,k)].x -= fac2*u[IDX(i,j,k)].x;
      force[IDX(i,j,k)].y -= fac2*u[IDX(i,j,k)].y;  
      force[IDX(i,j,k)].z -= fac2*u[IDX(i,j,k)].z;  
#endif

      }/* end of i,j,k for loop */

}
#endif
