#include "common_object.h"

#ifdef LB_FLUID
 #ifdef LB_FLUID_LES
my_double tau_u_les(int i , int j , int k){

  tensor grad_u;
  my_double s_norm, nu_les, tau_les; 
  my_double C_smag = 0.16; /* See Leveque , Toschi JFM on shear improved smagorinski model */
  my_double lx,ly,lz,vol;

  
  grad_u = gradient_vector(u,i,j,k);

  /* Norm of S */
  s_norm = sqrt( 0.5*( (grad_u.xx + grad_u.xx)*(grad_u.xx + grad_u.xx) + 
		        (grad_u.xy + grad_u.yx)*(grad_u.xy + grad_u.yx) +
		        (grad_u.xz + grad_u.zx)*(grad_u.xz + grad_u.zx) +
		        (grad_u.yx + grad_u.xy)*(grad_u.yx + grad_u.xy) + 
		        (grad_u.yy + grad_u.yy)*(grad_u.yy + grad_u.yy) +
		        (grad_u.yz + grad_u.zy)*(grad_u.yz + grad_u.zy) +
		        (grad_u.zx + grad_u.xz)*(grad_u.zx + grad_u.xz) + 
		        (grad_u.zy + grad_u.yz)*(grad_u.zy + grad_u.yz) +
		        (grad_u.zz + grad_u.zz)*(grad_u.zz + grad_u.zz) ) );


   #ifdef LB_FLUID_LES_SISM
   my_double L, f_v, c_exp,var_mean;
   my_double s_norm_mean;
   tensor grad_u_mean;

   L=(my_double)property.SY/2.0;
   //   f_v = 0.0043;
   f_v = 0.0033626;  /* For Kalayn: please put a reference or explain the origin of this number */
   c_exp = 3.628 * (f_v/L) * property.time_dt;

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

  C_smag = 0.18;

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

		 /* eddy viscosity */
		    nu_les = s_norm * pow( C_smag * vol, 2.0 ); 		 
		    		    

#ifdef METHOD_STREAMING
  tau_les = 3.0*(nu_les + property.nu) + 0.5 ;		 
#else
#ifdef METHOD_REDEFINED_POP
  tau_les = 3.0*(nu_les + property.nu) + 0.5*property.time_dt ;
#else
  tau_les = 3.0*(nu_les + property.nu) ;
#endif
#endif

  //fprintf(stderr,"%e\n",s_norm);
return tau_les;
}		 	  
#endif
#endif

#ifdef LB_TEMPERATURE
#ifdef LB_TEMPERATURE_LES
my_double tau_t_les(int i , int j , int k){

  my_double tau_temp,tau_les;
  my_double prandtl = property.nu/property.kappa;

    tau_les = tau_u_les(i , j , k);

#ifdef METHOD_STREAMING
    tau_temp = (tau_les-0.5)/prandtl + 0.5 ;
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
#ifdef LB_SCALAR_LES
my_double tau_s_les(int i , int j , int k){

  my_double tau_scal,tau_les;
    my_double schmidt = property.nu/property.chi;

    tau_les = tau_u_les(i , j , k);

#ifdef METHOD_STREAMING
    tau_scal = (tau_les-0.5)/schmidt + 0.5 ;
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





