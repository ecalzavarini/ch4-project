#include "common_object.h"


#ifdef LB_TEMPERATURE_MELTING 
void melting(){
  int i,j,k;
  my_double fl, Tl;
  my_double hs, hl;
  my_double Ts,Cp,Lf,enthalpy;
  my_double fac1, fac2;
  my_double eps = 1.0;
  my_double L,y;
  #ifdef LB_TEMPERATURE_ZIQI
  my_double Cp_ice; //ZIQI added
  #endif
  my_double liquid_frac_status;
  my_double liquid_frac_status_all;
  
 #ifdef LB_TEMPERATURE_MELTING_UNDEFORMABLE 
  my_double *mean_y , *mean_y_local;
  my_double lx, lz, norm;
 #endif /* end of LB_TEMPERATURE_MELTING_UNDEFORMABLE */

  Ts = (my_double) property.T_solid; 
  Cp = (my_double) property.specific_heat;//heat capacity (of the liquid)
  #ifdef LB_TEMPERATURE_ZIQI
  Cp_ice = (my_double) property.specific_heat_ice;//heat capacity (of the ice), ZIQI added
  #endif
  Lf = (my_double) property.latent_heat;
  /* solidification enthalpy */
  hs = Cp*Ts;
  hl = hs+Lf;
  Tl = hl/Cp;
  
#ifdef LB_TEMPERATURE_NEW_ENTHALPY_ZIQI
  //compute enthalpy based on paper 2019-Energy
  hs = Cp_ice*Ts;
  hl = hs+Lf;
  Tl = hl/Cp;
#endif

  fac1 = (Lf/Cp)/property.time_dt;

 #ifdef LB_TEMPERATURE_MELTING_UNDEFORMABLE 
	/* alloc lf_y_local and lf_y and set them to zero */
	mean_y  = (my_double*) malloc(sizeof(my_double)*NY);
	mean_y_local  = (my_double*) malloc(sizeof(my_double)*NY);
	set_to_zero_scalar(mean_y_local,NY);
	set_to_zero_scalar(mean_y,NY);

	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {
			  /* computing surface element from local from mesh */			  
			  lx = (mesh[IDXG(i+1, j, k)].x - mesh[IDXG(i, j, k)].x);
			  lz = (mesh[IDXG(i, j, k+1)].z - mesh[IDXG(i, j, k)].z);			  			  
			  //mean_y_local[j - BRD + LNY_START] += liquid_frac[IDX(i, j, k)]*lx*lz;
			  mean_y_local[j - BRD + LNY_START] += ( Cp*t[IDX(i,j,k)] + Lf*liquid_frac[IDX(i,j,k)] )*lx*lz;
        }

  	MPI_Allreduce(mean_y_local, mean_y, NY, MPI_MY_DOUBLE, MPI_SUM_my_double, MPI_COMM_WORLD );
	
	norm = 1.0/(my_double)(property.SX*property.SZ);
	for (i = 0; i < NY; i++){
		mean_y[i] *= norm;
  	}
 #endif /* end of LB_TEMPERATURE_MELTING_UNDEFORMABLE */


   for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++)
	  {//begin of i, j, k, for loop
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

 #ifndef LB_TEMPERATURE_MELTING_UNDEFORMABLE 
      /* This is the normal way, we compute the local Enthalpy*/
      /* compute Enthalphy */
	enthalpy = Cp*t[IDX(i,j,k)] + Lf*liquid_frac[IDX(i,j,k)];/* Here we assumes the same heath capacity in both solid and liquid (ice and water) */
 #endif /* end of if not def LB_TEMPERATURE_MELTING_UNDEFORMABLE */

 #ifdef LB_TEMPERATURE_MELTING_UNDEFORMABLE 
       /* This is for the flat interface case, we insert the x,z averaged enthalpy (or fluid fraction, first method) */
        //enthalpy = Cp*t[IDX(i,j,k)] + Lf*mean_y[j - BRD + LNY_START];
	enthalpy = mean_y[j - BRD + LNY_START];	
 #endif /* end of LB_TEMPERATURE_MELTING_UNDEFORMABLE */

 #ifdef LB_TEMPERATURE_NEW_ENTHALPY_ZIQI 
	my_double fi, fw;
	//compute enthalpy based on paper 2019-Energy
	if (t[IDX(i,j,k)]<Ts) {fi = 1.0; fw = 0.0;};
	if (t[IDX(i,j,k)]>Ts) {fi = 0.0; fw = 1.0;};
	if (t[IDX(i,j,k)]==Ts) {fi = 0.5; fw = 0.5;};
	
	enthalpy = Cp_ice*t[IDX(i,j,k)]*fi +
		+ (Cp*(t[IDX(i,j,k)]-Ts)+Cp_ice*Ts)*fw
		+ Lf*liquid_frac[IDX(i,j,k)]; //no solidus or liquidus, here we assume Ts=Tl
 #endif //end of LB_TEMPERATURE_NEW_ENTHALPY_ZIQI 

      /* compute new fluid fraction */ 
      if(enthalpy < hs) liquid_frac[IDX(i,j,k)]=0.0; //in solid (ice)
      else if(enthalpy > hl) liquid_frac[IDX(i,j,k)]=1.0; //in liquid
      else liquid_frac[IDX(i,j,k)]= (enthalpy - hs)/(hl - hs); //at interface

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

 #ifdef LB_TEMPERATURE_MELTING_UNDEFORMABLE
	/* free arrays */
	free(mean_y);
	free(mean_y_local);
 #endif /* end of LB_TEMPERATURE_MELTING_UNDEFORMABLE */


 #ifdef LB_TEMPERATURE_MELTING_CHECK_REACH_YP
	liquid_frac_status = liquid_frac_status_all = 0.0; 
	j = LNY+BRD-1;

	for (i = BRD; i < LNX+BRD; i++)
	  for (k = BRD; k < LNZ+BRD; k++) {
	    if(LNY_END == NY){
	     if( liquid_frac[IDX(i,j,k)] == 1.0 ){ liquid_frac_status = 1.0;
	      //	      fprintf(stderr,"i j k lf %d %d %d %e\n",i,j,k,liquid_frac[IDX(i,j,k)]);	
	     }
	    }
        }

  	MPI_Allreduce(&liquid_frac_status, &liquid_frac_status_all, 1 , MPI_MY_DOUBLE, MPI_SUM_my_double, MPI_COMM_WORLD );

	// if(ROOT) fprintf(stderr,"status %e\n",liquid_frac_status_all);	
	if(liquid_frac_status_all != 0.0){
	  if(ROOT) fprintf(stderr,"Liquid melt has reach the y top wall! Exit.\n");
	  MPI_Finalize();
	  exit(-2);
	}
   
 #endif /* end of LB_TEMPERATURE_MELTING_CHECK_YP */


}
#endif
