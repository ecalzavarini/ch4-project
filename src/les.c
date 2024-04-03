#include "common_object.h"

#ifdef LB_FLUID
 #ifdef LB_FLUID_LES
my_double tau_u_les(int i , int j , int k){

  tensor grad_u, SS, QQ;
  vector grad_t; /*Added by Luz*/
  my_double s_norm, nu_les, tau_les; 
  my_double C_smag = 0.18; /* Leveque , Toschi JFM on shear improved smagorinski model uses instead 0.16*/
  my_double lx,ly,lz,vol;
  my_double fac;
  my_double q_norm;
  /* for WALE */
  tensor WW,JJ;
  my_double iso_const, s_norm2, j_norm2, C_wale, op_wale;


 /* Norm of S  can be computed in two ways */ 
 #ifdef LB_FLUID_LES_LOCAL_STRAIN
/* 1) local computation of strain rate tensor and its norm from the populations */ 
  QQ=stress_tensor(p, i, j, k);
  /* the following factor is used to transform the strss tensor to the rate-of-stain tensor 
  note that our trick , we use the previous value of the relaxation time except at the first time step */
  if (itime == 1) tau_u_les_total[IDX(i,j,k)] = property.tau_u;
  fac = -1.0/( 2*cs2*dens[IDX(i,j,k)]*tau_u_les_total[IDX(i,j,k)]*property.time_dt );
  //fac = -1.0/( 2*cs2*dens[IDX(i,j,k)]*property.tau_u*property.time_dt );
  s_norm = fabs(fac) * sqrt(2.0 * (QQ.xx*QQ.xx + QQ.xy*QQ.xy + QQ.xz*QQ.xz +
                       QQ.yx*QQ.yx + QQ.yy*QQ.yy + QQ.yz*QQ.yz +
                       QQ.zx*QQ.zx + QQ.zy*QQ.zy + QQ.zz*QQ.zz) );                     
 #else
/* 2) calculation from finite difference velocity gradients */
  grad_u = gradient_vector(u,i,j,k);
/* compute |S| = sqrt(2 S_ij S_ij) */
  s_norm = sqrt( 0.5*( (grad_u.xx + grad_u.xx)*(grad_u.xx + grad_u.xx) + 
		        (grad_u.xy + grad_u.yx)*(grad_u.xy + grad_u.yx) +
		        (grad_u.xz + grad_u.zx)*(grad_u.xz + grad_u.zx) +
		        (grad_u.yx + grad_u.xy)*(grad_u.yx + grad_u.xy) + 
		        (grad_u.yy + grad_u.yy)*(grad_u.yy + grad_u.yy) +
		        (grad_u.yz + grad_u.zy)*(grad_u.yz + grad_u.zy) +
		        (grad_u.zx + grad_u.xz)*(grad_u.zx + grad_u.xz) + 
		        (grad_u.zy + grad_u.yz)*(grad_u.zy + grad_u.yz) +
		        (grad_u.zz + grad_u.zz)*(grad_u.zz + grad_u.zz) ) );
 #endif

   #ifdef LB_FLUID_LES_SISM
   my_double L, f_v, c_exp,var_mean;
   my_double s_norm_mean;
   tensor grad_u_mean;

   L=(my_double)property.SY/2.0;
   //   f_v = 0.0043;
   f_v = 0.0033626;  /* For Kalayn: please put a reference or explain the origin of this number */
   c_exp = 3.628 * (f_v/L) * property.time_dt;
   c_exp = property.time_dt/100; /*enrico test because dt/Tmemory */

    /* Initialize */
				if(itime==1) {
				  u_mean[IDX(i, j, k)].x = u[IDX(i, j, k)].x;
				  u_mean[IDX(i, j, k)].y = u[IDX(i, j, k)].y;
				  u_mean[IDX(i, j, k)].z = u[IDX(i, j, k)].z;
				}

    #ifndef LB_FLUID_LES_SISM_KALMAN
				/* EXPONENTIAL FILTER APPROACH */
				u_mean[IDX(i, j, k)].x = (1.0 - c_exp) * u_mean[IDX(i, j, k)].x + c_exp * u[IDX(i, j, k)].x;
				u_mean[IDX(i, j, k)].y = (1.0 - c_exp) * u_mean[IDX(i, j, k)].y + c_exp * u[IDX(i, j, k)].y;
				u_mean[IDX(i, j, k)].z = (1.0 - c_exp) * u_mean[IDX(i, j, k)].z + c_exp * u[IDX(i, j, k)].z;

    #else
				/* KALMAN FILTER APPROACH */				
    var_mean = 3.628 * (f_v/L) * property.time_dt * f_v;

    /* Initialize */
				if(itime==1) {

				  P_kalman[IDX(i, j, k)].x = var_mean * var_mean;
				  P_kalman[IDX(i, j, k)].y = var_mean * var_mean;
				  P_kalman[IDX(i, j, k)].z = var_mean * var_mean;

				  sqr_var[IDX(i, j, k)].x = var_mean * var_mean;
				  sqr_var[IDX(i, j, k)].y = var_mean * var_mean;
				  sqr_var[IDX(i, j, k)].z = var_mean * var_mean;
				}
  /* Predict */
  u_mean_kalman_pre[IDX(i, j, k)].x = u_mean[IDX(i, j, k)].x;
  u_mean_kalman_pre[IDX(i, j, k)].y = u_mean[IDX(i, j, k)].y;
  u_mean_kalman_pre[IDX(i, j, k)].z = u_mean[IDX(i, j, k)].z;

  P_kalman_pre[IDX(i, j, k)].x = P_kalman[IDX(i, j, k)].x + var_mean * var_mean;
  P_kalman_pre[IDX(i, j, k)].y = P_kalman[IDX(i, j, k)].y + var_mean * var_mean;
  P_kalman_pre[IDX(i, j, k)].z = P_kalman[IDX(i, j, k)].z + var_mean * var_mean;

  /* Update */
  K_kalman[IDX(i, j, k)].x = P_kalman_pre[IDX(i, j, k)].x * (1.0/(P_kalman_pre[IDX(i, j, k)].x + sqr_var[IDX(i, j, k)].x));
  K_kalman[IDX(i, j, k)].y = P_kalman_pre[IDX(i, j, k)].y * (1.0/(P_kalman_pre[IDX(i, j, k)].y + sqr_var[IDX(i, j, k)].y));
  K_kalman[IDX(i, j, k)].z = P_kalman_pre[IDX(i, j, k)].z * (1.0/(P_kalman_pre[IDX(i, j, k)].z + sqr_var[IDX(i, j, k)].z));

  u_mean[IDX(i, j, k)].x = u_mean_kalman_pre[IDX(i, j, k)].x +  K_kalman[IDX(i, j, k)].x * (u[IDX(i, j, k)].x - u_mean_kalman_pre[IDX(i, j, k)].x);
  u_mean[IDX(i, j, k)].y = u_mean_kalman_pre[IDX(i, j, k)].y +  K_kalman[IDX(i, j, k)].y * (u[IDX(i, j, k)].y - u_mean_kalman_pre[IDX(i, j, k)].y);
  u_mean[IDX(i, j, k)].z = u_mean_kalman_pre[IDX(i, j, k)].z +  K_kalman[IDX(i, j, k)].z * (u[IDX(i, j, k)].z - u_mean_kalman_pre[IDX(i, j, k)].z);

  P_kalman[IDX(i, j, k)].x = P_kalman_pre[IDX(i, j, k)].x - (K_kalman[IDX(i, j, k)].x * P_kalman_pre[IDX(i, j, k)].x);
  P_kalman[IDX(i, j, k)].y = P_kalman_pre[IDX(i, j, k)].y - (K_kalman[IDX(i, j, k)].y * P_kalman_pre[IDX(i, j, k)].y);
  P_kalman[IDX(i, j, k)].z = P_kalman_pre[IDX(i, j, k)].z - (K_kalman[IDX(i, j, k)].z * P_kalman_pre[IDX(i, j, k)].z);

  sqr_var_kalman[IDX(i, j, k)].x = f_v * fabs(u_mean[IDX(i, j, k)].x - u[IDX(i, j, k)].x);
  sqr_var_kalman[IDX(i, j, k)].y = f_v * fabs(u_mean[IDX(i, j, k)].y - u[IDX(i, j, k)].y);
  sqr_var_kalman[IDX(i, j, k)].z = f_v * fabs(u_mean[IDX(i, j, k)].z - u[IDX(i, j, k)].z);

  sqr_var[IDX(i, j, k)].x = MAX(sqr_var_kalman[IDX(i, j, k)].x, 0.1*f_v*f_v);
  sqr_var[IDX(i, j, k)].y = MAX(sqr_var_kalman[IDX(i, j, k)].y, 0.1*f_v*f_v);
  sqr_var[IDX(i, j, k)].z = MAX(sqr_var_kalman[IDX(i, j, k)].z, 0.1*f_v*f_v);
  #endif

   /* Calculating gradient vector*/
  grad_u_mean = gradient_vector(u_mean,i,j,k);

  /* Norm of S_mean */
  s_norm_mean = sqrt( 0.5*( (grad_u_mean.xx + grad_u_mean.xx)*(grad_u_mean.xx + grad_u_mean.xx) + 
		        (grad_u_mean.xy + grad_u_mean.yx)*(grad_u_mean.xy + grad_u_mean.yx) +
		        (grad_u_mean.xz + grad_u_mean.zx)*(grad_u_mean.xz + grad_u_mean.zx) +
		        (grad_u_mean.yx + grad_u_mean.xy)*(grad_u_mean.yx + grad_u_mean.xy) + 
		        (grad_u_mean.yy + grad_u_mean.yy)*(grad_u_mean.yy + grad_u_mean.yy) +
		        (grad_u_mean.yz + grad_u_mean.zy)*(grad_u_mean.yz + grad_u_mean.zy) +
		        (grad_u_mean.zx + grad_u_mean.xz)*(grad_u_mean.zx + grad_u_mean.xz) + 
		        (grad_u_mean.zy + grad_u_mean.yz)*(grad_u_mean.zy + grad_u_mean.yz) +
			(grad_u_mean.zz + grad_u_mean.zz)*(grad_u_mean.zz + grad_u_mean.zz) ) );
  if(s_norm > s_norm_mean) {
    s_norm =  s_norm - s_norm_mean;
  }
  #endif

  //C_smag = 0.18;

 #ifdef LB_FLUID_LES_VANDRIEST
  if(j < (property.SY/2.0)){C_smag = 0.18 * (1.0 - exp((-(property.Amp_x * center_V[IDX(i,j,k)].y) / property.nu)/26.0));}
  else {C_smag = 0.18 * (1.0 - exp((-(property.Amp_x * ((property.SY - center_V[IDX(i,j,k)].y)) / property.nu)/26.0)));}
  /*
  if (i==64 && k==64)
  {fprintf(stderr," i j k  C_Smag %ld %ld %ld %lf \n",i,j,k,C_smag);}*/
 #endif

#ifdef METHOD_STREAMING
		          vol = 1.0;
#else  		 

			  /* computing volume local from mesh */			  
			  lx = (mesh[IDXG(i+1, j, k)].x - mesh[IDXG(i, j, k)].x);
			  ly = (mesh[IDXG(i, j+1, k)].y - mesh[IDXG(i, j, k)].y);
			  lz = (mesh[IDXG(i, j, k+1)].z - mesh[IDXG(i, j, k)].z);			  
                          vol = lx*ly*lz;
			  vol = pow(vol,1./3.);
#endif

  #ifdef LB_FLUID_LES_SMAGORINSKY_LILLY /* This section was added by Luz */
   my_double N2; /* Brunt Väisälä Frequanecy*/
   my_double Ri; /* Richardson number N²/|S|², |S| = (2SijSij), Sij: Strain rate tensor */
   my_double Cb; /* Buoyancy correction term */
   my_double prandtl = 1.0; /* usual choice is pr=O(1) */ 

   grad_t = gradient_scalar(t,i,j,k); /*Added by Luz*/
 
   if(grad_t.y < 0){N2 = 0;}
   else{N2 = property.beta_t * property.gravity_y * grad_t.y;}
 
   Ri = N2 / pow( s_norm , 2.0 ); 
 
   if(Ri < 1){Cb = sqrt(1 - (Ri/prandtl));}
   else{Cb = 0;}
   
   s_norm = s_norm * Cb;
 
   //fprintf(stderr,"vol C_Smag grad_t.y gravity_y %lf %lf %lf %lf\n",vol, C_smag, grad_t.y, property.gravity_y);
   //fprintf(stderr,"s_norm_lilly %lf \n", s_norm);
  #endif

	/* eddy viscosity */
  //fprintf(stderr,"s_norm_ent %lf \n", s_norm);
	nu_les = s_norm * pow( C_smag * vol, 2.0 ); 		 
		    		    
  #ifdef LB_FLUID_LES_WALE
  /* 
  WALE : Wall-Adapting Local Eddy-Viscosity Model
  Nicoud, F. & Ducros, F. 1999 
  Subgrid-scale stress modelling based on the square of the velocity gradient tensor. 
  Flow, Turbulence and Combustion 62, 183–200.
  */
  /* rate of strain tensor */
  SS.xx = grad_u.xx;
  SS.yy = grad_u.yy;
  SS.zz = grad_u.zz;
  SS.xy = SS.yx = 0.5*(grad_u.xy + grad_u.yx); 
  SS.xz = SS.zx = 0.5*(grad_u.xz + grad_u.zx); 
  SS.zy = SS.yz = 0.5*(grad_u.zy + grad_u.yz); 
  /* rate of rotation tensor */
  WW.xx = WW.yy = WW.zz = 0.0;
  WW.xy = 0.5*(grad_u.xy - grad_u.yx); 
  WW.xz = 0.5*(grad_u.xz - grad_u.zx); 
  WW.yz = 0.5*(grad_u.zy - grad_u.yz); 
  WW.yx = - WW.xy;
  WW.zx = - WW.xz;
  WW.zy = - WW.yz;
/* compute the scalar term SS_ij*SS_ij - WW_ij * WW_ij */
  iso_const = SS.xx*SS.xx + SS.xy*SS.xy + SS.xz*SS.xz 
            + SS.yx*SS.yx + SS.yy*SS.yy + SS.yz*SS.yz 
            + SS.zx*SS.zx + SS.zy*SS.zy + SS.zz*SS.zz 
            - WW.xx*WW.xx - WW.xy*WW.xy - WW.xz*WW.xz 
            - WW.yx*WW.yx - WW.yy*WW.yy - WW.yz*WW.yz 
            - WW.zx*WW.zx - WW.zy*WW.zy - WW.zz*WW.zz;
/* JJ_{ij} = SS_{ik}*SS_{kj} + WW_{ik}*WW_{kj} - delta_{ij}*iso_const/3  */            
  JJ.xx = SS.xx*SS.xx + SS.xy*SS.yx + SS.xz*SS.zx 
        + WW.xx*WW.xx + WW.xy*WW.yx + WW.xz*WW.zx - iso_const/3.0;
  JJ.xy = SS.xx*SS.xy + SS.xy*SS.yy + SS.xz*SS.zy 
        + WW.xx*WW.xy + WW.xy*WW.yy + WW.xz*WW.zy;
  JJ.xz = SS.xx*SS.xz + SS.xy*SS.yz + SS.xz*SS.zz 
        + WW.xx*WW.xz + WW.xy*WW.yz + WW.xz*WW.zz;
  JJ.yx = SS.yx*SS.xx + SS.yy*SS.yx + SS.yz*SS.zx 
        + WW.yx*WW.xx + WW.yy*WW.yx + WW.yz*WW.zx;
  JJ.yy = SS.yx*SS.xy + SS.yy*SS.yy + SS.yz*SS.zy 
        + WW.yx*WW.xy + WW.yy*WW.yy + WW.yz*WW.zy  - iso_const/3.0;
  JJ.yz = SS.yx*SS.xz + SS.yy*SS.yz + SS.yz*SS.zz 
        + WW.yx*WW.xz + WW.yy*WW.yz + WW.yz*WW.zz;
  JJ.zx = SS.zx*SS.xx + SS.zy*SS.yx + SS.zz*SS.zx 
        + WW.zx*WW.xx + WW.zy*WW.yx + WW.zz*WW.zx;
  JJ.zy = SS.zx*SS.xy + SS.zy*SS.yy + SS.zz*SS.zy 
        + WW.zx*WW.xy + WW.zy*WW.yy + WW.zz*WW.zy; 
  JJ.zz = SS.zx*SS.xz + SS.zy*SS.yz + SS.zz*SS.zz 
        + WW.zx*WW.xz + WW.zy*WW.yz + WW.zz*WW.zz - iso_const/3.0; 
   /* compute  S_ij * S_ij */
  /*
  s_norm2 = SS.xx*SS.xx + SS.xy*SS.xy + SS.xz*SS.xz
          + SS.yx*SS.yx + SS.yy*SS.yy + SS.yz*SS.yz
          + SS.zx*SS.zx + SS.zy*SS.zy + SS.zz*SS.zz;
  */
  /* faster */
  s_norm2 = pow(s_norm,2.0)/2.0;

  /* This is J_ij * J_ij */
  j_norm2 = JJ.xx*JJ.xx + JJ.xy*JJ.xy + JJ.xz*JJ.xz
          + JJ.yx*JJ.yx + JJ.yy*JJ.yy + JJ.yz*JJ.yz
          + JJ.zx*JJ.zx + JJ.zy*JJ.zy + JJ.zz*JJ.zz;

  /* op_wale = (J J)^3/2 / ( (S S)^5/2 + (J J)^5/4 )  */
  op_wale = pow(j_norm2,1.5)/( pow(s_norm2,2.5) +  pow(j_norm2,1.25) );
  C_wale = 0.5; 
  nu_les = op_wale * pow( C_wale * vol, 2.0 );
  #endif


#ifdef METHOD_STREAMING
//fprintf(stderr,"nu_les_entrada %lf \n", nu_les);
  tau_les = 3.0*(nu_les + property.nu) + 0.5 ;		 
  //fprintf(stdout,"nu_les %e %d %d %d %d\n", nu_les,i,j,k,itime);	
#else
#ifdef METHOD_REDEFINED_POP
  tau_les = 3.0*(nu_les + property.nu) + 0.5*property.time_dt ;
#else
  tau_les = 3.0*(nu_les + property.nu) ;
#endif
#endif

  //fprintf(stderr,"%e\n",s_norm);

  /* the data are stored in a field , to be reused in other parts of the code */
  tau_u_les_total[IDX(i, j, k)] = tau_les;

return tau_les;
}		 	  
#endif
#endif

#ifdef LB_TEMPERATURE
#if defined(LB_TEMPERATURE_LES) && defined(LB_FLUID_LES)
my_double tau_t_les(int i , int j , int k){

  my_double tau_temp,tau_les;
  my_double prandtl = property.nu/property.kappa;
  /* Usual choice is Pr=O(1) */
  prandtl = 1.0;  //0.85;

    //tau_les = tau_u_les(i , j , k);
    tau_les = tau_u_les_total[IDX(i, j, k)];

#ifdef METHOD_STREAMING
    tau_temp = 3.0*(property.kappa + (tau_les-0.5-3*property.nu)/(3*prandtl)) + 0.5 ; /* Luz modified */
    //tau_temp = (tau_les-0.5)/prandtl + 0.5 ;
#else
#ifdef METHOD_REDEFINED_POP
    tau_temp = (tau_les-0.5*property.time_dt)/prandtl + 0.5*property.time_dt ;
#else
    tau_temp = tau_les/prandtl ;
#endif
#endif

    return tau_temp;
}
#endif
#endif


#ifdef LB_SCALAR
/* to be corrected as the temperature one */
#ifdef LB_SCALAR_LES
my_double tau_s_les(int i , int j , int k){

  my_double tau_scal,tau_les;
  my_double schmidt = property.nu/property.chi;
  schmidt = 1.0; /* same as for the turbulent Prandtl */

  tau_les = tau_u_les(i , j , k);

#ifdef METHOD_STREAMING
    //tau_scal = (tau_les-0.5)/schmidt + 0.5 ;
    tau_scal = 3.0*(property.chi + (tau_les-0.5-3*property.nu)/(3*schmidt)) + 0.5 ; /* same as modified by Luz */
#else
#ifdef METHOD_REDEFINED_POP
    tau_scal = (tau_les-0.5*property.time_dt)/schmidt + 0.5*property.time_dt ;
#else
    tau_scal = tau_les/schmidt ;
#endif
#endif

    return tau_scal;
}
#endif
#endif





