#include "common_object.h"


#ifdef LB_TEMPERATURE_MELTING
void melting(){
  int i,j,k, pp;
  my_double  fl;
  my_double hs, hl , temp , Tl;
  my_double eps = 1.0;
  my_double Ts,Cp,Lf;

  Ts = (my_double) property.T_solid;
  Cp = (my_double) property.specific_heat;
  Lf = (my_double) property.latent_heat;
  /* solidification temperature */
  hs = Cp*Ts;
  hl = hs+Lf;
  Tl = hl/Cp;

  //  if(i==0) fprintf(stderr,"Ts = %g \n Tl = %g\n",(double)Ts, (double)Tl); 

  /* store previous fluid fraction */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){
	liquid_frac_old[IDX(i,j,k)] = liquid_frac[IDX(i,j,k)];
      }

   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){

      //#ifdef SALT
      /* subtract the freezing point depression */
      //Ts = property.liquidus_slope * ss[IDX(y,x)];
      //hs = Cp*Ts;
      //hl = hs+Lf;
      //#endif

      /* compute Entalphy */
      enthalpy[IDX(i,j,k)] = Cp*t[IDX(i,j,k)] + Lf*liquid_frac[IDX(i,j,k)];
      /* compute new fluid fraction */
      if(enthalpy[IDX(i,j,k)] < hs) liquid_frac[IDX(i,j,k)]=0.0;
      else if(enthalpy[IDX(i,j,k)] > hl) liquid_frac[IDX(i,j,k)]=1.0;
      else liquid_frac[IDX(i,j,k)]= (enthalpy[IDX(i,j,k)] - hs)/(hl - hs);
    }

  /* add melting term to the temperature field */
  temp = (Lf/Cp);
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){
	//for (pp=0; pp<NPOP; pp++) g[IDX(i,j,k)].p[pp] -= wgt[pp]*temp*(liquid_frac[IDX(i,j,k)]-liquid_frac_old[IDX(i,j,k)]);
        t_source[IDX(i,j,k)] -= temp*(liquid_frac[IDX(i,j,k)]-liquid_frac_old[IDX(i,j,k)]);
    }


  /* add drag term to the force structure */
   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){
      fl = liquid_frac[IDX(i,j,k)]; 
      temp = (1.0 - fl*fl)/(eps + fl*fl*fl);
      force[IDX(i,j,k)].x -= temp*u[IDX(i,j,k)].x;
      force[IDX(i,j,k)].y -= temp*u[IDX(i,j,k)].y;  
      force[IDX(i,j,k)].z -= temp*u[IDX(i,j,k)].z;  
    }


   //#ifdef THERMIC_SOLID
   //  for (y=1; y<NY+1; y++)
   //    for (x=1; x<NX+1; x++){
   //      fl = liquid_frac[IDX(y,x)];
   //      tauT[IDX(y,x)] = fl*tau2 + (1.0-fl)*tau4; 
   //    }
   //#endif

}
#endif
