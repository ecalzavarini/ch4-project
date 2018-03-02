#include "common_object.h"


#if (defined LB_FLUID_FORCING || defined LB_TEMPERATURE_FORCING || defined LB_SCALAR_FORCING)
void build_forcing(){
  int i, j, k , ii,jj;
  my_double fnx,fny,fnz,kn,omega_t;
  my_double x,y,z;
  my_double LX,LY,LZ,nu;
  my_double temp, fac, val,fac1,fac2; 
  vector dist,vel;
  vector u0, u0_all;
  vector u2,u2_all;
  vector f0,f0_all;
  vector coeff,vec;
  my_double norm, k_turb, k0_turb , k_ratio;
  my_double t0,t0_all , temp_linear, t_turnover,phi_diffusivity,t0_coeff;
  my_double s0,s0_all,s0_coeff;
  my_double local_depth, radiation_at_bottom, reflection_ceff,lf;
  my_double radiation_at_bottom1,radiation_at_bottom2,radiation_at_bottom3,radiation_at_bottom4;
  FILE *fout;

   LX = (my_double)(property.SX);
   LY = (my_double)(property.SY);
   LZ = (my_double)(property.SZ);

   nu = property.nu;

#ifdef LB_FLUID_FORCING_HIT  /* useful for HOMOGENEOUS ISOTROPIC TURBULENCE */
 #ifdef LB_FLUID_FORCING_HIT_RANDOM
   /* the phases are random */
      if(ROOT){ 
	for (ii=0; ii<nk; ii++){
	  phi_u[ii].x = myrand();
	  phi_u[ii].y = myrand();
	  phi_u[ii].z = myrand();
	}
      }
 #else
   /* the phases make a random walk */
      //fac = sqrt(property.time_dt * 2.0 * (property.SY/0.1)/pow(two_pi,2.0) ); 
      t_turnover = (property.SY/0.1);
      phi_diffusivity = 0.01/t_turnover;
      fac1 = property.time_dt/t_turnover;
      fac2 = sqrt(2.0*property.time_dt*phi_diffusivity);
      if(ROOT){ 
	for (ii=0; ii<nk; ii++){
	  /*
	  val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
	  phi_u[ii].x += val*fac;
	  val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
	  phi_u[ii].y += val*fac;
	  val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
	  phi_u[ii].z += val*fac;
	  */
	  val=random_gauss(0.0,1.0);
	  phi_u[ii].x += -phi_u[ii].x*fac1 + val*fac2;
	  val=random_gauss(0.0,1.0);
	  phi_u[ii].y += -phi_u[ii].y*fac1 + val*fac2;
	  val=random_gauss(0.0,1.0);
	  phi_u[ii].z += -phi_u[ii].z*fac1 + val*fac2;
	}
	/* OUTPUT PHASES */ 
  #ifdef LB_FLUID_FORCING_HIT_DEBUG
	fout = fopen("phases.dat","a");
	fprintf(fout," %e %e %e %e\n", time_now, phi_u[0].x, phi_u[0].y ,phi_u[0].z);
	fclose(fout); 
  #endif	
      }
 #endif
    /* the phases are broadcasted */
   MPI_Bcast(phi_u, nk, MPI_vector_type, 0, MPI_COMM_WORLD);

 #ifdef LB_FLUID_FORCING_HIT_ZEROMODE 
    /* compute the zero mode intensity (the mean velocity) */
    u0.x = u0_all.x = 0.0;
    u0.y = u0_all.y = 0.0;
    u0.z = u0_all.z = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
	for(i=BRD;i<LNX+BRD;i++){ 
	    u0.x += u[IDX(i,j,k)].x;
            u0.y += u[IDX(i,j,k)].y;
            u0.z += u[IDX(i,j,k)].z;
	  }

    MPI_Allreduce(&u0, &u0_all, 1, MPI_vector_type, MPI_SUM_vector, MPI_COMM_WORLD );

    norm = sqrt(u0_all.x*u0_all.x + u0_all.y*u0_all.y + u0_all.z*u0_all.z);
    if(norm !=0.0){
      u0_all.x /= norm;
      u0_all.y /= norm;
      u0_all.z /= norm;
    }
 #endif
 #ifdef LB_FLUID_FORCING_HIT_LINEAR 
    /* compute the total turbulent kinetic energy */

    u0.x = u0_all.x = 0.0;
    u0.y = u0_all.y = 0.0;
    u0.z = u0_all.z = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
	for(i=BRD;i<LNX+BRD;i++){ 
	    u0.x += u[IDX(i,j,k)].x;
            u0.y += u[IDX(i,j,k)].y;
            u0.z += u[IDX(i,j,k)].z;
	  }

    MPI_Allreduce(&u0, &u0_all, 1, MPI_vector_type, MPI_SUM_vector, MPI_COMM_WORLD );

    norm = 1.0/(property.NX*property.NY*property.NZ);

      u0_all.x *= norm;
      u0_all.y *= norm;
      u0_all.z *= norm;

    u2.x = u2_all.x = 0.0;
    u2.y = u2_all.y = 0.0;
    u2.z = u2_all.z = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
	for(i=BRD;i<LNX+BRD;i++){ 

	  u2.x += (u[IDX(i,j,k)].x - u0_all.x)*(u[IDX(i,j,k)].x - u0_all.x);
          u2.y += (u[IDX(i,j,k)].y - u0_all.y)*(u[IDX(i,j,k)].y - u0_all.y);
          u2.z += (u[IDX(i,j,k)].z - u0_all.z)*(u[IDX(i,j,k)].z - u0_all.z);
	  }

   MPI_Allreduce(&u2, &u2_all, 1, MPI_vector_type, MPI_SUM_vector, MPI_COMM_WORLD );

      u2_all.x *= norm;
      u2_all.y *= norm;
      u2_all.z *= norm;

    k_turb =  0.5*( u2_all.x + u2_all.y + u2_all.z );
    k0_turb = 0.5*0.01;  /* assuming a velocity of 0.1 per grid point */
    if(k_turb != 0.0) k_ratio = k0_turb/k_turb; else k_ratio = 1.0;
 #endif
#endif

#ifdef LB_TEMPERATURE_FORCING_HIT
 #ifdef LB_TEMPERATURE_FORCING_HIT_RANDOM
   /* the phases are random */
      if(ROOT){ 
	for (ii=0; ii<nk; ii++){
	  phi_t[ii].x = myrand();
	  phi_t[ii].y = myrand();
	  phi_t[ii].z = myrand();
	}
      }
 #else
      /* the phases do random walk */
      //fac = sqrt(1.0/(my_double)randomization_itime_t);
      //      fac = sqrt(property.time_dt * 2.0 * (property.SY/0.1)/pow(two_pi,2.0) ); 
      //fac = sqrt(property.time_dt * 2.0  * (0.1/property.SY) );
      t_turnover = (property.SY/0.1);
      phi_diffusivity = 0.01/t_turnover;
      fac1 = property.time_dt/t_turnover;
      fac2 = sqrt(2.0*property.time_dt*phi_diffusivity);
      if(ROOT){ 	
	for (ii=0; ii<nk_t; ii++){
	/*
	  val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
          phi_t[ii].x += val*fac;
          val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
          phi_t[ii].y += val*fac;
          val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
          phi_t[ii].z += val*fac;
	*/
	  val=random_gauss(0.0,1.0);
	  phi_t[ii].x += -phi_t[ii].x*fac1 + val*fac2; 
	  val=random_gauss(0.0,1.0);
	  phi_t[ii].y += -phi_t[ii].y*fac1 + val*fac2; 
	  val=random_gauss(0.0,1.0);
	  phi_t[ii].z += -phi_t[ii].z*fac1 + val*fac2; 
	}
      }
#endif
      /* phases are boradcasted */
    MPI_Bcast(phi_t, nk_t, MPI_vector_type, 0, MPI_COMM_WORLD);

 #ifdef LB_TEMPERATURE_FORCING_HIT_ZEROMODE 
    /* compute the zero mode intensity (the mean temperature) */
    t0 = t0_all = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
	for(i=BRD;i<LNX+BRD;i++){ 
	    t0 += t[IDX(i,j,k)];
	  }

    MPI_Allreduce(&t0, &t0_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   norm = 1.0/(my_double)(property.SX*property.SY*property.SZ);
    if(norm !=0.0){
     t0_all *= norm;
    }
     t0_all=(t0_all > 0.0)?1.0:-1.0;
    //    if(ROOT)fprintf(stderr,"t0_all = %e\n",t0_all);
 #endif
#endif

#ifdef LB_SCALAR_FORCING_HIT
  #ifdef LB_SCALAR_FORCING_HIT_RANDOM
   /* the phases are random */
      if(ROOT){ 
	for (ii=0; ii<nk; ii++){
	  phi_s[ii].x = myrand();
	  phi_s[ii].y = myrand();
	  phi_s[ii].z = myrand();
	}
      }
 #else
      /* the phases do random walk */
      //fac = sqrt(1.0/(my_double)randomization_itime_s);
      //fac = sqrt(property.time_dt * 2.0 * (0.1 *property.SY) ); 
      //fac = sqrt(property.time_dt * 2.0  * (0.1/property.SY) );
      t_turnover = (property.SY/0.1);
      phi_diffusivity = 0.01/t_turnover;
      fac1 = property.time_dt/t_turnover;
      fac2 = sqrt(2.0*property.time_dt*phi_diffusivity);
      if(ROOT){ 
	for (ii=0; ii<nk_s; ii++){
	  /*
	  val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
          phi_s[ii].x += val*fac;
          val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
          phi_s[ii].y += val*fac;
          val=(2.0*(myrand()-0.5) > 0.0)?1.0:-1.0;
          phi_s[ii].z += val*fac;
	  */
	  val=random_gauss(0.0,1.0);
	  phi_s[ii].x += -phi_s[ii].x*fac1 + val*fac2; 
	  val=random_gauss(0.0,1.0);
	  phi_s[ii].y += -phi_s[ii].y*fac1 + val*fac2; 
	  val=random_gauss(0.0,1.0);
	  phi_s[ii].z += -phi_s[ii].z*fac1 + val*fac2; 
	}
      }
#endif
      /* phases are boradcasted */
    MPI_Bcast(phi_s, nk_s, MPI_vector_type, 0, MPI_COMM_WORLD);

 #ifdef LB_SCALAR_FORCING_HIT_ZEROMODE 
    /* compute the zero mode intensity (the mean temperature) */
    s0 = s0_all = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
	for(i=BRD;i<LNX+BRD;i++){ 
	    s0 += s[IDX(i,j,k)];
	  }

    MPI_Allreduce(&s0, &s0_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   norm = 1.0/(my_double)(property.SX*property.SY*property.SZ);
    if(norm !=0.0){
     s0_all *= norm;
    }
     s0_all=(s0_all > 0.0)?1.0:-1.0;
    //    if(ROOT)fprintf(stderr,"s0_all = %e\n",s0_all);
 #endif
#endif

     /* Communication for laplacian computation */
#ifdef LB_FLUID_FORCING_LAPLACIAN /* communicate boundaries for laplacian computation */
     sendrecv_borders_vector(u);
#endif
#ifdef LB_TEMPERATURE_FORCING_LAPLACIAN /* communicate boundaries for laplacian computation */
     sendrecv_borders_scalar(t);
#endif

     /* BEGIN LOOP on GRID */ 
 for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

   x = (my_double)center_V[IDX(i,j,k)].x;
   y = (my_double)center_V[IDX(i,j,k)].y;
   z = (my_double)center_V[IDX(i,j,k)].z;

#ifdef LB_FLUID_FORCING
	force[IDX(i,j,k)].x = 0.0;
	force[IDX(i,j,k)].y = 0.0;
	force[IDX(i,j,k)].z = 0.0;

 #ifdef LB_FLUID_FORCING_POISEUILLE 
	/* note that the LX,LY,LZ dependence (i.e. the non-homogeneous direction) is here just an arbitrary choice */
	/* Here Amp indicates the maximal velocity at the center of the channel , note that in poiseuille flow U_max = Force * L^2/(8 \nu)  */
	force[IDX(i,j,k)].x += 2.0*property.Amp_x*((4.*nu)*pow(LY,-2.0));  
	force[IDX(i,j,k)].y += 2.0*property.Amp_y*((4.*nu)*pow(LX,-2.0));  
	force[IDX(i,j,k)].z += 2.0*property.Amp_z*((4.*nu)*pow(LY,-2.0));  
 #endif

 #ifdef LB_FLUID_FORCING_CHANNEL
	/* note that the LX,LY,LZ dependence (i.e. the non-homogeneous direction) is here just an arbitrary choice */
	/* Here Amp indicates the maximal velocity at the center of the turbulent channel , note that in this case we write U_max = sqrt(Force*L/2)  */
	force[IDX(i,j,k)].x += 2.0*pow(property.Amp_x,2.0)/LY;
	force[IDX(i,j,k)].y += 2.0*pow(property.Amp_y,2.0)/LX;
	force[IDX(i,j,k)].z += 2.0*pow(property.Amp_z,2.0)/LY;
 #endif

 #ifdef LB_FLUID_FORCING_CHANNEL_CONSTANT_POWER
	/* Here Amp indicates the power Force*velocity  */
	/* note that out_all.uxt is computed only when we dumpe the averages */	
	if(out_all.ux != 0) force[IDX(i,j,k)].x += property.Amp_x/out_all.ux; 
	if(out_all.uy != 0) force[IDX(i,j,k)].y += property.Amp_y/out_all.uy;
	if(out_all.uz != 0) force[IDX(i,j,k)].z += property.Amp_z/out_all.uz;
 #endif

 #ifdef LB_FLUID_FORCING_KOLMOGOROV 
    kn=1.0;
    fnx=nu*pow(two_pi/LX,2.0);
    fny=nu*pow(two_pi/LY,2.0);
    //y = (my_double)center_V[IDX(i,j,k)].y;
    //x = (my_double)center_V[IDX(i,j,k)].x;

	/* along x */  
        force[IDX(i,j,k)].x += property.Amp_x*fny*sin(kn*two_pi*y/LY); 
	force[IDX(i,j,k)].y += property.Amp_y*fnx*sin(kn*two_pi*x/LX); 
	force[IDX(i,j,k)].z += property.Amp_z*fny*sin(kn*two_pi*y/LY);  

	//   fprintf(stderr,"force[IDX(i,j,k)].x %e force[IDX(i,j,k)].y %e force[IDX(i,j,k)].z %e\n",force[IDX(i,j,k)].x, force[IDX(i,j,k)].y,force[IDX(i,j,k)].z);
	//  fprintf(stderr,"property.Amp_x %e property.Amp_y %e property.Amp_z %e\n",property.Amp_x, property.Amp_y,property.Amp_z);
	//  exit(1);
 #endif  


 #ifdef LB_FLUID_FORCING_CELLULAR
  #ifndef LB_FLUID_FORCING_CELLULAR_UNSTEADY 
    kn=0.5;
	/* along x */  
        force[IDX(i,j,k)].x += property.Amp_x*sin(kn*two_pi*x/LX)*cos(kn*two_pi*y/LY); 
	force[IDX(i,j,k)].y -= property.Amp_x*cos(kn*two_pi*x/LX)*sin(kn*two_pi*y/LY); 
	force[IDX(i,j,k)].z += 0.0; 
 
  #else
    kn=0.5;
    omega_t = two_pi*property.time_dt*(0.1/LX);
	/* along x */  
        force[IDX(i,j,k)].x += property.Amp_x*sin(kn*two_pi*x/LX)*cos(kn*two_pi*y/LY)*sin(omega_t*itime); 
	force[IDX(i,j,k)].y -= property.Amp_x*cos(kn*two_pi*x/LX)*sin(kn*two_pi*y/LY)*sin(omega_t*itime); 
	force[IDX(i,j,k)].z += 0.0; 

  #endif
 #endif  

 #ifdef LB_FLUID_FORCING_LAPLACIAN /* the forcing term has the form nu_add*laplacian( u ), where nu_add is an additional viscosity */
	vec = laplacian_vector(u, i, j, k);
        force[IDX(i,j,k)].x += property.nu_add * vec.x; 
	force[IDX(i,j,k)].y += property.nu_add * vec.y; 
	force[IDX(i,j,k)].z += property.nu_add * vec.z; 
 #endif

 #ifdef LB_FLUID_FORCING_HIT  /* for HOMOGENEOUS ISOTROPIC TURBULENCE */     
  #ifdef LB_FLUID_FORCING_HIT_LINEAR
   /* As in Phares L. Carroll and G. Blanquart PHYSICS OF FLUIDS 25, 105114 (2013) */
	//fac = 0.5*((out_all.ux2 - out_all.ux*out_all.ux) + (out_all.uy2 - out_all.uy*out_all.uy) + (out_all.uz2 - out_all.uz*out_all.uz));	
	//if(fac != 0.0) fac = 1.0/fac; else fac = 1.0;	
	fac = k_ratio;
	force[IDX(i,j,k)].x += fac*property.Amp_x*(u[IDX(i,j,k)].x - u0_all.x);
	force[IDX(i,j,k)].y += fac*property.Amp_y*(u[IDX(i,j,k)].y - u0_all.y);
	force[IDX(i,j,k)].z += fac*property.Amp_z*(u[IDX(i,j,k)].z - u0_all.z);
  #else
  #ifdef LB_FLUID_FORCING_HIT_ZEROMODE 
      /* the zero mode */
      //fac = sqrt(out_all.ux*out_all.ux + out_all.uy*out_all.uy + out_all.uz*out_all.uz);
      //if(fac != 0.0) fac = 1./fac; else fac = 1.0;
      force[IDX(i,j,k)].x += property.Amp_x*(- u0_all.x);
      force[IDX(i,j,k)].y += property.Amp_y*(- u0_all.y);
      force[IDX(i,j,k)].z += property.Amp_z*(- u0_all.z);
  #endif
   #ifdef LB_FLUID_FORCING_HIT_TYPE2
      /* Type 2 forcing , taken from Computers and Mathematics with Applications vol. 58 (2009) pag. 1055-1061
	 "Lattice Boltzmann simulations of homogeneous isotropic turbulence" by Waleed Abdel Kareema, Seiichiro Izawa , Ao-Kui Xiong , Yu Fukunishi  */
    for(ii=0; ii<nk; ii++){
      //fac = pow(vk2[ii],-5./6.) * property.Amp_x * sin(two_pi*(vk[ii].x*x/LX  + vk[ii].y*y/LY + vk[ii].z*z/LZ ) + two_pi*phi_u[ii].x) / vk2[ii]; 
      fac = property.Amp_x * sin(two_pi*(vk[ii].x*x/LX  + vk[ii].y*y/LY + vk[ii].z*z/LZ ) + two_pi*phi_u[ii].x) / vk2[ii];      
      force[IDX(i,j,k)].x += 2.0*fac*(vk[ii].y * vk[ii].z); 
      force[IDX(i,j,k)].y -=     fac*(vk[ii].x * vk[ii].z); 
      force[IDX(i,j,k)].z -=     fac*(vk[ii].x * vk[ii].y);
    }
   #else 
     /* Force at large scale, similar to Federico, Prasad forcing */
    for(ii=0; ii<nk; ii++){
      fac = pow(vk2[ii],-2./6.);
      force[IDX(i,j,k)].x += fac*property.Amp_x* ( sin(two_pi*(vk[ii].y*y/LY + phi_u[ii].y)) + sin(two_pi*(vk[ii].z*z/LZ + phi_u[ii].z)) );
      force[IDX(i,j,k)].y += fac*property.Amp_y* ( sin(two_pi*(vk[ii].x*x/LX + phi_u[ii].x)) + sin(two_pi*(vk[ii].z*z/LZ + phi_u[ii].z)) );
      force[IDX(i,j,k)].z += fac*property.Amp_z* ( sin(two_pi*(vk[ii].y*y/LY + phi_u[ii].y)) + sin(two_pi*(vk[ii].x*x/LX + phi_u[ii].x)) );
      /* diagonal part: be careful it shall satisfy the incompressibility condition */      	
      /*
      force[IDX(i,j,k)].x +=  2.0*fac*property.Amp_x * sin(two_pi*(vk[ii].x*x/LX + phi_u[0].x)); 
      force[IDX(i,j,k)].y +=  -fac*property.Amp_x * sin(two_pi*(vk[ii].y*y/LY + phi_u[0].x)); 
      force[IDX(i,j,k)].z +=  -fac*property.Amp_x * sin(two_pi*(vk[ii].z*z/LZ + phi_u[0].x));       
      */
    }
   #endif

  #endif  
 #endif

 #ifdef LB_FLUID_FORCING_ABSORB  /* attempt to implement an absorbing layer */
    fac = 0.8;
    dist.x = (x/LX - fac)/(1.0-fac);
    dist.y = (y/LY - fac)/(1.0-fac);
    dist.z = (z/LZ - fac)/(1.0-fac);

    if(dist.y>0 ){

    vel.x = 0.0; //out_all.ux;
    vel.y = out_all.uy;
    vel.z = 0.0; //out_all.uz;

    dist.x = 0.5 + pow(dist.x,2.0);
    //    dist.y = 1.0 / ( 0.5 + pow(dist.y,2.0) );
    //dist.y = 1.0;
    dist.z = 0.5 + pow(dist.z,2.0);  

    //fprintf(stderr,"%d %d %d %d %e\n",me, i,j,k, u[IDX(i,j,k)].y);
        
    //force[IDX(i,j,k)].x += -dist.y*(vel.x - u[IDX(i,j,k)].x);  
    //force[IDX(i,j,k)].y += -u[IDX(i,j,k)].y;  // dist.y*(vel.y  -u[IDX(i,j,k)].y);
    //force[IDX(i,j,k)].z += -dist.y*(vel.z  -u[IDX(i,j,k)].z);
    
     }
 #endif

 #ifdef LB_TEMPERATURE
  #ifdef LB_TEMPERATURE_BUOYANCY
	//my_double temp, fac;
   #ifdef LB_TEMPERATURE_BUOYANCY_T0_REF
  temp = (t[IDX(i,j,k)] - property.T_ref);
   #elif defined(LB_TEMPERATURE_BUOYANCY_T0_TOP)
  temp = (t[IDX(i,j,k)] - property.T_top);
   #elif defined(LB_TEMPERATURE_BUOYANCY_T0_BOT)
  temp = (t[IDX(i,j,k)] - property.T_bot);
   #elif defined(LB_TEMPERATURE_BUOYANCY_T0_GRAD)
  /* subtract to the temperature the linear profile */
  temp_linear = 0.0; //-(property.deltaT/property.SY)*center_V[IDX(i,j,k)].y + 0.5*property.deltaT; 
  temp = (t[IDX(i,j,k)] - temp_linear );
   #elif defined(LB_TEMPERATURE_BUOYANCY_T0_REF2)
  temp = (t[IDX(i,j,k)] - property.T_ref2);
   #else
  /* the good one for RB , T0 = T_mean*/  
  temp =  t[IDX(i,j,k)] - 0.5*(property.T_bot + property.T_top);
   #endif

  //temp =  t[IDX(i,j,k)] - (-(property.deltaT/property.SY)*center_V[IDX(i,j,k)].y + property.T_bot) ;
  fac = property.beta_t*temp; 
  if(property.beta2_t != 0.0) fac += property.beta2_t*temp*temp;
  //fac = property.beta_t;
  
     /* This way if we are in a solid the force is not applied */
  /*  Unnecessary , this is done already in melting.c
     #ifdef LB_TEMPERATURE_MELTING
       fac *= liquid_frac[IDX(i,j,k)]; 
     #endif
  */
      force[IDX(i,j,k)].x += fac*property.gravity_x;
      force[IDX(i,j,k)].y += fac*property.gravity_y;
      force[IDX(i,j,k)].z += fac*property.gravity_z;

      //fprintf(stderr, "fy %e\n",property.gravity_y);
  #endif /* LB_TEMPERATURE_BUOYANCY */
 #endif /* LB_TEMPERATURE */

 #ifdef LB_SCALAR_BUOYANCY
      //my_double temp, fac;
  #ifdef LB_SCALAR_BUOYANCY_SREF
      temp = (s[IDX(i,j,k)] - property.S_top);
  #else
      temp =  s[IDX(i,j,k)] - 0.5*(property.S_bot + property.S_top);
  #endif

      fac = property.beta_s*temp;

      force[IDX(i,j,k)].x += fac*property.gravity_x;
      force[IDX(i,j,k)].y += fac*property.gravity_y;
      force[IDX(i,j,k)].z += fac*property.gravity_z;

 #endif


 #ifdef  LB_FLUID_FORCING_PENALIZATION
      my_double mask;
      /* penalization of a cube */
      /*  
      if( fabs(center_V[IDX(i,j,k)].x-property.SX/2.0) < 10 &&
	  fabs(center_V[IDX(i,j,k)].y) < 10 && 
	  fabs(center_V[IDX(i,j,k)].z-property.SZ/2.0) < 10  ) 
	mask=1.0; 
      else 
	mask=0.0;

      if( mask == 1.0 ){
	force[IDX(i,j,k)].x = -u[IDX(i,j,k)].x;  
	force[IDX(i,j,k)].y = -u[IDX(i,j,k)].y;
	force[IDX(i,j,k)].z = -u[IDX(i,j,k)].z;
	  }
      */

  #ifdef LB_FLUID_FORCING_LANDSCAPE
       if(landscape[IDX(i, j, k)]>0.0){
	force[IDX(i,j,k)].x = -u[IDX(i,j,k)].x;  
	force[IDX(i,j,k)].y = -u[IDX(i,j,k)].y;
	force[IDX(i,j,k)].z = -u[IDX(i,j,k)].z;	
      }    
  #else
      /* small central spot penalization */
      mask = pow(center_V[IDX(i,j,k)].x-property.SX/2.0, 2.0)+pow(center_V[IDX(i,j,k)].y-property.SY/2.0, 2.0);
      if( mask < 10.0 ){       
	force[IDX(i,j,k)].x = -u[IDX(i,j,k)].x;  
	force[IDX(i,j,k)].y = -u[IDX(i,j,k)].y;
	force[IDX(i,j,k)].z = -u[IDX(i,j,k)].z;	
      }     
  #endif
 #endif /* endif of LB_FLUID_FORCING_PENALIZATION */


      /* Set forcing to zero if the domain has no thickness (e.g. in 2 dimensions , just 1 grid points in z-direction ) */
      /* TO BE REMOVED, causes problems to 1,50,1 poiseuillle flow */
      /*
      if(NX==1) force[IDX(i,j,k)].x = 0.0;
      if(NY==1) force[IDX(i,j,k)].y = 0.0;
      if(NZ==1) force[IDX(i,j,k)].z = 0.0;
      */

#endif /* end of LB_FLUID_FORCING */

#ifdef LB_TEMPERATURE
/* From here HERE SOURCE TERM ON TEMPERATURE FIELD */
 #ifdef LB_TEMPERATURE_FORCING
      /* set to zero */ 
      t_source[IDX(i,j,k)] = 0.0;

   /* here we can for instance impose a temperature profile , or add a thermal source or make the field reactive*/
 
 #ifdef LB_TEMPERATURE_FORCING_SOURCE
  /* mimic source term  */
      my_double spot;
      /* penalization of a cube */
      /*        
      if( fabs(center_V[IDX(i,j,k)].x-property.SX/2.0) < 10 &&
	  fabs(center_V[IDX(i,j,k)].y) < 10 && 
	  fabs(center_V[IDX(i,j,k)].z-property.SZ/2.0) < 10  ) 
	spot=1.0; 
      else 
	spot=0.0;

      if( spot == 1.0 ) t_source[IDX(i,j,k)] = -(t[IDX(i,j,k)] - property.T_bot);
      */
      /* small central spot penalization */
      spot = pow(center_V[IDX(i,j,k)].x-property.SX/2.0, 2.0)+pow(center_V[IDX(i,j,k)].y-property.SY/2.0, 2.0);
      if( spot < 1.0 ) t_source[IDX(i,j,k)] = -(t[IDX(i,j,k)] - property.T_bot);
 #endif


 #ifdef LB_TEMPERATURE_FORCING_ABSORB  /* attempt to implement an absorbing layer for temperature */
    fac = 0.9;
    dist.x = (x/LX - fac)/(1.0-fac);
    dist.y = (y/LY - fac)/(1.0-fac);
    dist.z = (z/LZ - fac)/(1.0-fac);

    if(dist.y>0){

    //dist.x = 0.5 + pow(dist.x,2.0);
    // dist.y = 1.0 / ( 0.5 + pow(dist.y,2.0) );
    //dist.z = 0.5 + pow(dist.z,2.0);  
        
     t_source[IDX(i,j,k)] = -(t[IDX(i,j,k)] - property.T_top);
       }
 #endif


 #ifdef LB_TEMPERATURE_FORCING_PROFILE
  /* impose a mean linear temperature profile , note that bc for temp shall be set to 0 */
  t_source[IDX(i,j,k)] = (property.deltaT/property.SY)*u[IDX(i,j,k)].y;
 #endif

 #ifdef LB_TEMPERATURE_FORCING_REACTION
  /* make the field reactive */
  #ifdef LB_TEMPERATURE_FORCING_REACTION_FKPP
  /* implement Fisher-KPP */
  /* we use property.Amp_t as reaction rate ,  and property.T_top as the carrying capacity */
   #ifdef LB_TEMPERATURE_FORCING_REACTION_FKPP_FLUCTUATION
    /* the field is inteded as a fluctuation respect to the global mean value "carrying capacity"/2 = property.T_top/2  */
    t_source[IDX(i,j,k)] = property.Amp_t*(property.T_top/4.0 - (t[IDX(i,j,k)]*t[IDX(i,j,k)])/property.T_top );
   #else
    t_source[IDX(i,j,k)] = property.Amp_t*t[IDX(i,j,k)]*(1.0 - t[IDX(i,j,k)]/property.T_top);
   #endif
  #endif 
  #ifdef LB_TEMPERATURE_FORCING_REACTION_ORDER1
    /* implement a first order chemical reaction (or Malthusian growth) */
    /* NOTE : be careful that if other temperature forcings are activated the same parameter property.Amp_t might be in use also there */
    t_source[IDX(i,j,k)] = property.Amp_t*t[IDX(i,j,k)];
  #endif  
 #endif


 #ifdef LB_TEMPERATURE_FORCING_MONOD
  #ifdef LB_SCALAR
  /* make the field react to the scalar concentration */
  t_source[IDX(i,j,k)] = property.Amp_t * s[IDX(i,j,k)]/(0.5 + s[IDX(i,j,k)]) * t[IDX(i,j,k)];
  #endif
 #endif

 #ifdef LB_TEMPERATURE_FORCING_BULK
   t_source[IDX(i,j,k)] = property.Amp_t;
  #ifdef LB_TEMPERATURE_FORCING_BULK_VARIABLE
   if (center_V[IDX(i,j,k)].x > property.SX/2.) t_source[IDX(i,j,k)] = 0.0;
  #endif
  #ifdef LB_TEMPERATURE_MELTING
   /* This way if we are in a solid the bulk source is not applied */
   if (liquid_frac[IDX(i,j,k)]>0.1){
     t_source[IDX(i,j,k)] = liquid_frac[IDX(i,j,k)]*property.Amp_t;
   }else{
     t_source[IDX(i,j,k)] = 0.0;
   }
  #endif
 #endif

 #ifdef LB_TEMPERATURE_FORCING_LAPLACIAN /* the forcing term has the form kappa_add*laplacian( t ), where kappa_add is an additional viscosity */
	val = laplacian_scalar(t, i, j, k);
        t_source[IDX(i,j,k)] += property.kappa_add * val; 
 #endif

 #ifdef LB_TEMPERATURE_FORCING_RADIATION 
  #ifdef LB_TEMPERATURE_FORCING_RADIATION_SOLAR
 /*
   Electromagnetic spectrum (see wikipedia):
   lambda: UV<400nm; VIS: 400nm<lambda<800nm; IR(near+far): lambda>800nm 
   
   Typical values of solar flux (Skyllingstad and Paulson (2007))
   lambda:   (350-700)nm  (700-900)nm   (900-1100)nm   >1100nm
   band index:         m=1          m=2           m=3          m=4
   Pm:                 0.481        0.194         0.123        0.202
   Km:                 0.180        3.250        27.500      300.000

   Where: 
   Pm is the fraction of solar intensity in the given band (The sum of all Pm is therefore =1)
   The actual intensity is I0=Frn*Pm with Frn=160 W/m^2 (order 100 W/m^2)
   Km = alpha , the attenuation coefficient. 
  */
  /*
  t_source[IDX(i,j,k)]  = 0.481 * property.Amp_t* 0.180 *exp(-0.180  * center_V[IDX(i,j,k)].y); 
  t_source[IDX(i,j,k)] += 0.194 * property.Amp_t* 3.250 *exp(-3.250  * center_V[IDX(i,j,k)].y); 
  t_source[IDX(i,j,k)] += 0.123 * property.Amp_t* 27.50 *exp(-27.50  * center_V[IDX(i,j,k)].y); 
  t_source[IDX(i,j,k)] += 0.202 * property.Amp_t* 300.0 *exp(-300.0  * center_V[IDX(i,j,k)].y);  
  */
  t_source[IDX(i,j,k)]  =                 property.Amp_t*property.attenuation*                exp(-property.attenuation * center_V[IDX(i,j,k)].y); 
  t_source[IDX(i,j,k)] += (0.194/0.481) * property.Amp_t*property.attenuation* (3.250/0.180) *exp(-(3.250/0.180)  * center_V[IDX(i,j,k)].y); 
  t_source[IDX(i,j,k)] += (0.123/0.481) * property.Amp_t*property.attenuation* (27.50/0.180) *exp(-(27.50/0.180)  * center_V[IDX(i,j,k)].y); 
  t_source[IDX(i,j,k)] += (0.202/0.481) * property.Amp_t*property.attenuation* (300.0/0.180) *exp(-(300.0/0.180)  * center_V[IDX(i,j,k)].y); 

  #else
  /* just a monocromatic light source */
  t_source[IDX(i,j,k)] = property.Amp_t*property.attenuation*exp(-property.attenuation*center_V[IDX(i,j,k)].y); 
  #endif
  #ifdef LB_TEMPERATURE_FORCING_RADIATION_REFLECTION
    #ifdef LB_TEMPERATURE_MELTING
   /* In general the depth of the fluid layer is < property.SY 
     therefore a fraction of the radiation is reflected and a part is transmitted through the solid medium */
   /* The "local_depth" is computed from the melt fraction */
   lf = 0.0;
   for (jj = BRD; jj < LNY+BRD; jj++) lf += liquid_frac[IDX(i, jj, k)]*(mesh[IDXG(i, jj+1, k)].y - mesh[IDXG(i, jj, k)].y);
   MPI_Allreduce(&lf, &local_depth, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   //if(ROOT) fprintf(stderr,"i %d k %d local_depth %e\n",i,k,local_depth); 

   radiation_at_bottom = property.Amp_t*property.attenuation*exp(-property.attenuation*local_depth); /* the radiation intensity at the bottom of the fluid layer */

   reflection_ceff = 1.0; /* determines the albedo , it shall be given as an external parameter */

   if(center_V[IDX(i,j,k)].y > local_depth){
     /* remove the radiation going through the bottom wall */
     t_source[IDX(i,j,k)] -= property.Amp_t*property.attenuation*exp(-property.attenuation*center_V[IDX(i,j,k)].y); 
     /* add the radiation reflected from bottom */
     t_source[IDX(i,j,k)] += reflection_ceff*radiation_at_bottom*exp(-property.attenuation*(local_depth - center_V[IDX(i,j,k)].y)); 
     /* add the radiation the transmitted through the bottom */    
     t_source[IDX(i,j,k)] += (1.0-reflection_ceff)*radiation_at_bottom*exp(-property.attenuation*center_V[IDX(i,j,k)].y); 
   }

    #else
   /* NO MELTING : in this case the depth of the fluid layer is = property.SY */
   /* the radiation is reflected, the transmitted one just goes out of the system */
     #ifdef LB_TEMPERATURE_FORCING_RADIATION_SOLAR
     local_depth = property.SY;  /* the depth of the fluid layer */

     reflection_ceff = 1.0; /* determines the albedo , it shall be given as an external parameter */

     /* the radiation intensity at the bottom of the fluid layer */
     radiation_at_bottom1 =                 property.Amp_t*property.attenuation*                exp(-property.attenuation * local_depth); 
     radiation_at_bottom2 = (0.194/0.481) * property.Amp_t*property.attenuation* (3.250/0.180) *exp(-(3.250/0.180)  * local_depth); 
     radiation_at_bottom3 = (0.123/0.481) * property.Amp_t*property.attenuation* (27.50/0.180) *exp(-(27.50/0.180)  * local_depth); 
     radiation_at_bottom4 = (0.202/0.481) * property.Amp_t*property.attenuation* (300.0/0.180) *exp(-(300.0/0.180)  * local_depth); 

     t_source[IDX(i,j,k)] += reflection_ceff*radiation_at_bottom1*exp(-property.attenuation*(local_depth - center_V[IDX(i,j,k)].y)); 
     t_source[IDX(i,j,k)] += reflection_ceff*radiation_at_bottom2*exp(-(3.250/0.180)       *(local_depth - center_V[IDX(i,j,k)].y)); 
     t_source[IDX(i,j,k)] += reflection_ceff*radiation_at_bottom3*exp(-(27.50/0.180)       *(local_depth - center_V[IDX(i,j,k)].y)); 
     t_source[IDX(i,j,k)] += reflection_ceff*radiation_at_bottom4*exp(-(300.0/0.180)       *(local_depth - center_V[IDX(i,j,k)].y)); 

     #else
     local_depth = property.SY;  /* the depth of the fluid layer */

     radiation_at_bottom = property.Amp_t*property.attenuation*exp(-property.attenuation*local_depth); /* the radiation intensity at the bottom of the fluid layer */

     reflection_ceff = 1.0; /* determines the albedo , it shall be given as an external parameter */

     t_source[IDX(i,j,k)] += reflection_ceff*radiation_at_bottom*exp(-property.attenuation*(local_depth - center_V[IDX(i,j,k)].y)); 
     #endif
    #endif
  #endif
 #endif 

 #ifdef LB_TEMPERATURE_FORCING_HIT /* for HOMOGENEOUS ISOTROPIC TURBULENCE on TEMPERATURE */     
  #ifdef LB_TEMPERATURE_FORCING_HIT_LINEAR
   /* Inspired by Phares L. Carroll and G. Blanquart PHYSICS OF FLUIDS 25, 105114 (2013) */	
      if(out_all.t2 != 0.0) fac = 1.0/out_all.t2; else fac = 1.0;
      t_source[IDX(i,j,k)] += fac*property.Amp_t*t[IDX(i,j,k)];
  #else
  #ifdef LB_TEMPERATURE_FORCING_HIT_ZEROMODE 
      /* the zero mode */
      t_source[IDX(i,j,k)] += property.Amp_t*(-t0_all);
  #endif
      /*the other modes */
    for(ii=0; ii<nk_t; ii++){
      fac = pow(vk2_t[ii],-2./6.);
      t_source[IDX(i,j,k)] += fac*property.Amp_t*( sin(two_pi*(vk_t[ii].x*x/LX + phi_t[ii].x)) + sin(two_pi*(vk_t[ii].y*y/LY + phi_t[ii].y)) + sin(two_pi*(vk_t[ii].z*z/LZ + phi_t[ii].z)) );
      }
  #endif
 #endif

 #ifdef LB_TEMPERATURE_FORCING_GRAD /* force temperature with constant gradient (along y) */
    fac = property.deltaT/property.SY;
    //t_source[IDX(i,j,k)] += -property.Amp_t*u[IDX(i,j,k)].y;
    t_source[IDX(i,j,k)] += fac*u[IDX(i,j,k)].y;
 #endif

 #endif /* LB_TEMPERATURE_FORCING */
#endif /* LB_TEMPERATURE */

/* From here HERE SOURCE TERM ON SCALAR FIELD */
#ifdef LB_SCALAR_FORCING
      /* set to zero */ 
      s_source[IDX(i,j,k)] = 0.0;

 #ifdef LB_SCALAR_FORCING_REACTION
  /* make the field reactive */
  #ifdef LB_SCALAR_FORCING_REACTION_FKPP
  /* implement Fisher-KPP */
  /* we use property.Amp_s as reaction rate ,  and property.S_top as the carrying capacity */
   #ifdef LB_SCALAR_FORCING_REACTION_FKPP_FLUCTUATION
    /* the field is inteded as a fluctuation respect to the global mean value "carrying capacity"/2 = property.S_top/2  */
    s_source[IDX(i,j,k)] = property.Amp_s*(property.S_top/4.0 - (s[IDX(i,j,k)]*s[IDX(i,j,k)])/property.S_top );
   #else
    s_source[IDX(i,j,k)] = property.Amp_s*s[IDX(i,j,k)]*(1.0 - s[IDX(i,j,k)]/property.S_top);
   #endif
  #endif 
  #ifdef LB_SCALAR_FORCING_REACTION_ORDER1
    /* implement a first order chemical reaction (or Malthusian growth) */
    /* NOTE : be careful that if other temperature forcings are activated the same parameter property.Amp_t might be in use also there */
    s_source[IDX(i,j,k)] = property.Amp_s*s[IDX(i,j,k)];
  #endif  
 #endif


#ifdef LB_SCALAR_FORCING_MONOD
  #ifdef LB_TEMPERATURE
  /* make the field react to negatively to the t (phythoplancton) concentration */
  t_source[IDX(i,j,k)] = -property.Amp_s * s[IDX(i,j,k)]/(0.5 + s[IDX(i,j,k)]) * t[IDX(i,j,k)];
  #endif
#endif


#ifdef LB_SCALAR_FORCING_HIT /* for HOMOGENEOUS ISOTROPIC TURBULENCE on SCALAR */     
 #ifdef LB_SCALAR_FORCING_HIT_LINEAR
   /* Inspired by Phares L. Carroll and G. Blanquart PHYSICS OF FLUIDS 25, 105114 (2013) */	
      if(out_all.s2 != 0.0) fac = 1.0/out_all.s2; else fac = 1.0;
      s_source[IDX(i,j,k)] += fac*property.Amp_s*s[IDX(i,j,k)];
 #else
 #ifdef LB_SCALAR_FORCING_HIT_ZEROMODE 
      /* the zero mode */
      s_source[IDX(i,j,k)] += property.Amp_s*(-s0_all);
 #endif
      /*the other modes */
    for(ii=0; ii<nk_s; ii++){
      fac = pow(vk2_s[ii],-2./6.);
      s_source[IDX(i,j,k)] += fac*property.Amp_s*( sin(two_pi*(vk_s[ii].x*x/LX + phi_s[ii].x)) + sin(two_pi*(vk_s[ii].y*y/LY + phi_s[ii].y)) + sin(two_pi*(vk_s[ii].z*z/LZ + phi_s[ii].z)) );
    }
 #endif
#endif

#ifdef LB_SCALAR_FORCING_GRAD /* force scalar with constant gradient (along y) */
    s_source[IDX(i,j,k)] += -property.Amp_s*u[IDX(i,j,k)].y;
#endif

#endif

      }/* i,j,k */


#ifdef LB_FLUID_FORCING_NOZEROMODE
    /* compute the mean forcing value */
    f0.x = f0_all.x = 0.0;
    f0.y = f0_all.y = 0.0;
    f0.z = f0_all.z = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
        for(i=BRD;i<LNX+BRD;i++){ 
            f0.x += force[IDX(i,j,k)].x;
            f0.y += force[IDX(i,j,k)].y;
            f0.z += force[IDX(i,j,k)].z;
          }

    MPI_Allreduce(&f0, &f0_all, 1, MPI_vector_type, MPI_SUM_vector, MPI_COMM_WORLD );

    norm = 1.0/(property.NX*property.NY*property.NZ);

      f0_all.x *= norm;
      f0_all.y *= norm;
      f0_all.z *= norm;

     for (i = BRD; i < LNX+BRD; i++)
       for (j = BRD; j < LNY+BRD; j++)
	 for (k = BRD; k < LNZ+BRD; k++) {
	   force[IDX(i, j, k)].x -= f0_all.x;
	   force[IDX(i, j, k)].y -= f0_all.y;
	   force[IDX(i, j, k)].z -= f0_all.z;
	 }/* for i,j,k */
#endif

#ifdef LB_FLUID_FORCING_CONSTANT_POWER
     /* compute the total mean power value per unit volume, <F*u>_V / L^3  */
    f0.x = f0_all.x = 0.0;
    f0.y = f0_all.y = 0.0;
    f0.z = f0_all.z = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
        for(i=BRD;i<LNX+BRD;i++){ 
            f0.x += force[IDX(i,j,k)].x*u[IDX(i,j,k)].x;
            f0.y += force[IDX(i,j,k)].y*u[IDX(i,j,k)].y;
            f0.z += force[IDX(i,j,k)].z*u[IDX(i,j,k)].z;
          }

    MPI_Allreduce(&f0, &f0_all, 1, MPI_vector_type, MPI_SUM_vector, MPI_COMM_WORLD );

    norm = 1.0/(property.NX*property.NY*property.NZ);

    /* this is the global mean power per unit volume */
    fac = (f0_all.x + f0_all.y + f0_all.z)*norm;

    if(fac != 0.0){
      coeff.x = property.Amp_x/fac;
      coeff.y = property.Amp_y/fac;
      coeff.z = property.Amp_z/fac;
    }else{
      coeff.x = coeff.y = coeff.z = 1.0;
    }

     for (i = BRD; i < LNX+BRD; i++)
       for (j = BRD; j < LNY+BRD; j++)
	 for (k = BRD; k < LNZ+BRD; k++) {
	   force[IDX(i, j, k)].x *= coeff.x;
	   force[IDX(i, j, k)].y *= coeff.y;
	   force[IDX(i, j, k)].z *= coeff.z;
	 }/* for i,j,k */
#endif

#ifdef LB_TEMPERATURE
 #ifdef LB_FLUID_TEMPERATURE_NOZEROMODE
    /* compute the mean forcing value */
    t0 = t0_all = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
        for(i=BRD;i<LNX+BRD;i++){ 
            t0 += t_source[IDX(i,j,k)];
          }

    MPI_Allreduce(&t0, &t0_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    norm = 1.0/(property.NX*property.NY*property.NZ);

      t0_all *= norm;

     for (i = BRD; i < LNX+BRD; i++)
       for (j = BRD; j < LNY+BRD; j++)
	 for (k = BRD; k < LNZ+BRD; k++) {
	   t_source[IDX(i, j, k)] -= t0_all;
	 }/* for i,j,k */
 #endif
 #ifdef LB_TEMPERATURE_FORCING_CONSTANT_POWER
     /* compute the total mean power value per unit volume, <F*t>_V / L^3  */
    t0 = t0_all = 0.0;

    for(k=BRD;k<LNZ+BRD;k++)
      for(j=BRD;j<LNY+BRD;j++)
        for(i=BRD;i<LNX+BRD;i++){ 
            t0 += t_source[IDX(i,j,k)]*t[IDX(i,j,k)];
          }

    MPI_Allreduce(&t0, &t0_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    norm = 1.0/(property.NX*property.NY*property.NZ);

    /* this is the global mean power per unit volume */
    fac = t0_all *norm;

    if(fac != 0.0){
      t0_coeff = property.Amp_t/fac;
    }else{
      t0_coeff = 1.0;
    }

     for (i = BRD; i < LNX+BRD; i++)
       for (j = BRD; j < LNY+BRD; j++)
	 for (k = BRD; k < LNZ+BRD; k++) {
	   t_source[IDX(i, j, k)] *= t0_coeff;
	 }/* for i,j,k */
 #endif
#endif

}/* end of build_forcing */
#endif

#if (defined LB_FLUID_FORCING || defined LB_TEMPERATURE_FORCING || defined LB_SCALAR_FORCING)
void add_forcing(){
  int i, j, k, pp;
  my_double invtau,invtau_t,invtau_s;
  pop f_eq;
  my_double fac , fac_t, fac_s;
  vector d;
  my_double ux,uy,uz,cu,u2;
  vector vel;
  my_double rho ,temp , wgt2;
  pop p_eq , g_eq , h_eq;
  my_double mask;
  my_double magic_gamma;

#ifdef LB_FLUID
invtau = 1.0/property.tau_u;
 #ifdef METHOD_TRT
        /* This is two relaxation time TRT */
 magic_gamma = 0.25; //3./16.;
        /* see http://arxiv.org/abs/1508.07982v1 */
	/* this is invtau_minus */
        invtau = (4. - 2.*invtau)/(2.+(4.*magic_gamma -1.)*invtau);
 #endif

 fac = 1.0; /* this is for streaming if GUO is not defined */
#ifdef METHOD_FORCING_GUO
  fac = (1.0-0.5*property.time_dt*invtau);
#endif
#ifdef METHOD_REDEFINED_POP
  fac = 1.0;
#ifdef METHOD_REDEFINED_POP_GUO
  fac = (1.0-0.5*property.time_dt*invtau);
#endif
#endif
#endif /* LB_FLUID */

#ifdef LB_TEMPERATURE
invtau_t = 1.0/property.tau_t;

 #ifdef METHOD_TRT
        /* This is two relaxation time TRT */
        magic_gamma = 0.25;
        /* see http://arxiv.org/abs/1508.07982v1 */
	/* this is invtau_minus */
        invtau_t = (4. - 2.*invtau_t)/(2.+(4.*magic_gamma -1.)*invtau_t);
 #endif

#ifdef METHOD_FORCING_MALASPINAS
  fac_t = (1.0-0.5*property.time_dt*invtau_t);
#endif
#ifdef METHOD_REDEFINED_POP
  fac_t = 1.0;
#ifdef METHOD_REDEFINED_POP_GUO
  fac_t = (1.0-0.5*property.time_dt*invtau_t);
#endif
#endif
#endif /* LB_TEMPERATURE */


#ifdef LB_SCALAR
invtau_s = 1.0/property.tau_s;

 #ifdef METHOD_TRT
        /* This is two relaxation time TRT */
        magic_gamma = 0.25;
        /* see http://arxiv.org/abs/1508.07982v1 */
	/* this is invtau_minus */
        invtau_s = (4. - 2.*invtau_s)/(2.+(4.*magic_gamma -1.)*invtau_s);
 #endif

#ifdef METHOD_FORCING_MALASPINAS
  fac_s = (1.0-0.5*property.time_dt*invtau_s);
#endif
#ifdef METHOD_REDEFINED_POP
  fac_s = 1.0;
#ifdef METHOD_REDEFINED_POP_GUO
  fac_s = (1.0-0.5*property.time_dt*invtau_s);
#endif
#endif
#endif /* LB_SCALAR */


  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

	/* Here I prepare equilibrium pop for the direct forcing */
#ifdef  LB_FLUID_FORCING_DIRECT
      /* small central spot with velocity u=0 */
      mask = pow(center_V[IDX(i,j,k)].x-property.SX/2.0, 2.0)+pow(center_V[IDX(i,j,k)].y-property.SY/2.0, 2.0);
      vel.x = 0.0;
      vel.y = 0.0;
      vel.z = 0.0;
      rho = 1.0;
      p_eq = equilibrium_given_velocity(vel,rho);
#endif      

#ifdef  LB_TEMPERATURE_FORCING_DIRECT
      /* small central spot with velocity u=0 */
      mask = pow(center_V[IDX(i,j,k)].x-property.SX/2.0, 2.0)+pow(center_V[IDX(i,j,k)].y-property.SY/2.0, 2.0);
      vel.x = u[IDX(i,j,k)].x;
      vel.y = u[IDX(i,j,k)].y;
      vel.z = u[IDX(i,j,k)].z;
      temp = property.T_bot;
      g_eq = equilibrium_given_velocity(vel,temp);
#endif

      /* start loop on populations */
	for (pp=0; pp<NPOP; pp++){
	/* forcing */

#ifdef LB_FLUID_FORCING

#ifndef METHOD_FORCING_GUO  /* not defined METHOD_FORCING_GUO */
		  	  
	  /* Old version , simple but probably incomplete 
	    rho = dens[IDX(i,j,k)];	  
	    rhs_p[IDX(i,j,k)].p[pp] += invcs2*wgt[pp]*rho*force[IDX(i,j,k)].x*c[pp].x;
            rhs_p[IDX(i,j,k)].p[pp] += invcs2*wgt[pp]*rho*force[IDX(i,j,k)].y*c[pp].y;
            rhs_p[IDX(i,j,k)].p[pp] += invcs2*wgt[pp]*rho*force[IDX(i,j,k)].z*c[pp].z;
	  */

	  /* New version , like in GUO or in PhD thesis EPFL MALASPINAS */	  
	  
	rho = dens[IDX(i,j,k)];
	ux=u[IDX(i,j,k)].x;
	uy=u[IDX(i,j,k)].y;
	uz=u[IDX(i,j,k)].z;
        cu = (c[pp].x*ux + c[pp].y*uy + c[pp].z*uz);
        d.x = (c[pp].x-ux)*invcs2 + c[pp].x*cu*invcs4;
        d.y = (c[pp].y-uy)*invcs2 + c[pp].y*cu*invcs4;
        d.z = (c[pp].z-uz)*invcs2 + c[pp].z*cu*invcs4;

       rhs_p[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].x*d.x;
       rhs_p[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].y*d.y;
       rhs_p[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].z*d.z;   

	  
#ifdef METHOD_LOG
	    fac = 3.0*wgt[pp]*property.tau_u*exp(-p[IDX(i,j,k)].p[pp]*invtau);
	    rhs_p[IDX(i,j,k)].p[pp] += fac*force[IDX(i,j,k)].x*c[pp].x;
            rhs_p[IDX(i,j,k)].p[pp] += fac*force[IDX(i,j,k)].y*c[pp].y;
            rhs_p[IDX(i,j,k)].p[pp] += fac*force[IDX(i,j,k)].z*c[pp].z;
#endif
	    
#else   

	    /* This is METHOD_FORCING_GUO */    
	rho = dens[IDX(i,j,k)];  
	ux=u[IDX(i,j,k)].x;
	uy=u[IDX(i,j,k)].y;
	uz=u[IDX(i,j,k)].z;
        cu = (c[pp].x*ux + c[pp].y*uy + c[pp].z*uz);
        d.x = (c[pp].x-ux)*invcs2 + c[pp].x*cu*invcs4;
        d.y = (c[pp].y-uy)*invcs2 + c[pp].y*cu*invcs4;
        d.z = (c[pp].z-uz)*invcs2 + c[pp].z*cu*invcs4;
		
       rhs_p[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].x*d.x;
       rhs_p[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].y*d.y;
       rhs_p[IDX(i,j,k)].p[pp] += fac*wgt[pp]*rho*force[IDX(i,j,k)].z*d.z;
	    
#endif

#ifdef  LB_FLUID_FORCING_DIRECT

#ifdef LB_FLUID_FORCING_LANDSCAPE
       if(landscape[IDX(i, j, k)]>0.0){
#else
      if( sqrt(mask) < 10.0 ){
#endif
	/* this implementation works only when METHOD_EULER is used */
	rhs_p[IDX(i,j,k)].p[pp] =  (p_eq.p[pp] - p[IDX(i,j,k)].p[pp] )/property.time_dt;
	  }      
#endif
#endif

#ifdef LB_TEMPERATURE
#ifdef LB_TEMPERATURE_FORCING
 #ifndef METHOD_FORCING_MALASPINAS

  #ifdef LB_TEMPERATURE_FORCING_PAST
	/* a new version */
        ux = u[IDX(i, j, k)].x;
        uy = u[IDX(i, j, k)].y;
        uz = u[IDX(i, j, k)].z;                         
        u2 = (ux * ux + uy * uy + uz * uz);
	wgt2 = (1.0 + invcs2 * cu *((property.tau_t-0.5)/property.tau_t)  );
	rhs_g[IDX(i,j,k)].p[pp] += 1.5 * wgt2*wgt[pp]*t_source[IDX(i,j,k)];
        ux = old_u[IDX(i, j, k)].x;
        uy = old_u[IDX(i, j, k)].y;
        uz = old_u[IDX(i, j, k)].z;                         
        cu = (c[pp].x * ux + c[pp].y * uy + c[pp].z * uz);	
	wgt2 = (1.0 + invcs2 * cu *((property.tau_t-0.5)/property.tau_t)  );
	rhs_g[IDX(i,j,k)].p[pp] -= 0.5 * wgt2*wgt[pp]*old_t_source[IDX(i,j,k)];

	old_t_source[IDX(i,j,k)] = t_source[IDX(i,j,k)];
  #else
       /* Not Guo here */
       /* This is the simplest method and the one that should normally be used */
      rhs_g[IDX(i,j,k)].p[pp] += wgt[pp]*t_source[IDX(i,j,k)];      
  #endif

 #else
      /* The forcing is as in Malaspinas PhD pag. 93 (bottom) , with a correction factor 2 (found empirically)*/	    
      //ux = u[IDX(i, j, k)].x;
      //uy = u[IDX(i, j, k)].y;
      //uz = u[IDX(i, j, k)].z;                         
      //u2 = (ux * ux + uy * uy + uz * uz);
      //cu = (c[pp].x * ux + c[pp].y * uy + c[pp].z * uz);
      //wgt2=(1.0 + invcs2 * cu + invtwocs4 * cu * cu - invtwocs2 * u2);
      ////rhs_g[IDX(i,j,k)].p[pp] += fac_t*wgt2*wgt[pp]*t_source[IDX(i,j,k)];
      //rhs_g[IDX(i,j,k)].p[pp] += 2.0*fac_t*wgt2*wgt[pp]*t_source[IDX(i,j,k)];
			         
      /* not as in eq (8.40) , pp. 310 of Timm Krueger at al. book ---> there is an error */
      /* now as in Takeshi Seta, PHYSICAL REVIEW E 87, 063304 (2013) eq (22)  */ 
	rhs_g[IDX(i,j,k)].p[pp] += fac_t*wgt[pp]*t_source[IDX(i,j,k)];
 #endif


 #ifdef  LB_TEMPERATURE_FORCING_DIRECT
      if( sqrt(mask) < 10.0 ){
	/* this implementation works only when METHOD_EULER is used */
	rhs_g[IDX(i,j,k)].p[pp] =  (g_eq.p[pp] - g[IDX(i,j,k)].p[pp] )/property.time_dt;
	  }      
 #endif
#endif
#endif /* LB_TEMPERATURE */

#ifdef LB_SCALAR_FORCING
#ifndef METHOD_FORCING_MALASPINAS
       /* not Guo here  */
	    rhs_h[IDX(i,j,k)].p[pp] += wgt[pp]*s_source[IDX(i,j,k)];
#else
      //ux = u[IDX(i, j, k)].x;
      //uy = u[IDX(i, j, k)].y;
      //uz = u[IDX(i, j, k)].z;                         
      //u2 = (ux * ux + uy * uy + uz * uz);
      //cu = (c[pp].x * ux + c[pp].y * uy + c[pp].z * uz);
      //wgt2=(1.0 + invcs2 * cu + invtwocs4 * cu * cu - invtwocs2 * u2);
      ////rhs_h[IDX(i,j,k)].p[pp] += fac_s*wgt[pp]*wgt2*s_source[IDX(i,j,k)];
      //rhs_h[IDX(i,j,k)].p[pp] += 2.0*fac_s*wgt[pp]*wgt2*s_source[IDX(i,j,k)];

     /* not as in eq (8.40) , pp. 310 of Timm Krueger at al. book ---> there is an error */
     /* now as in Takeshi Seta, PHYSICAL REVIEW E 87, 063304 (2013) eq (22)  */ 
	rhs_h[IDX(i,j,k)].p[pp] += fac_t*wgt[pp]*s_source[IDX(i,j,k)];
#endif
#endif


	}/* pp */
      }/* i,j,k */
}
#endif


