#include "common_object.h"

my_double read_parameter(char * variable){

  char fnamein[256], fnameout[256];
  char name[256] = "NULL";
  double val;
  FILE *fin, *fout;
  int i;
  int cmp=1;

  sprintf(fnamein,"param.in");
  sprintf(fnameout,"param.dat");
  fin = fopen(fnamein,"r");
  if(fin == NULL){
	 fprintf(stderr,"Error message -> %s file is missing! Exit.\n",fnamein);
	 fclose(fin);
	 exit(0);
  }

  while(cmp!=0){  
    val=0.0;
    fscanf(fin,"%s %lf",&name,&val);
    cmp=strcmp(name,variable);
       if(cmp==0){
	 fout = fopen(fnameout,"a");
	 fprintf(fout,"%s %g\n",name, val);
	 fclose(fout);
	 return (my_double)val;
       }
       if(feof(fin)){ 
	 fprintf(stderr,"Error message -> %s not found. End of File. Exit.\n",variable);
	 fclose(fin);
	 exit(0);
       }    
  }
  fclose(fin);
}

void assign_parameters(){
  char name[256];
  int i;
  FILE *fout;

  if(ROOT){
    sprintf(OutDir,"RUN");
    mkdir(OutDir, S_IWUSR|S_IXUSR|S_IRUSR);
    fprintf(stderr,"OutDir is %s\n",OutDir);   

    remove("param.dat");
    remove("numbers.dat");

    /* write the active define flags  */
    /* Note works only if the executable is where the code is */
    system("cat define.h | grep '#def' | grep -v '//#def' > define.dat");

    /* From here on we read parameters from file */

    /* resume flags */
    sprintf(name,"resume");
    resume = (int)read_parameter(name);
#ifdef LB_FLUID
    sprintf(name,"resume_u");
    resume_u = (int)read_parameter(name);
#endif
#ifdef LB_TEMPERATURE
    sprintf(name,"resume_t");
    resume_t = (int)read_parameter(name);
#endif
#ifdef LB_SCALAR
    sprintf(name,"resume_s");
    resume_s = (int)read_parameter(name);
#endif

    /* size of center nodes grid */
    sprintf(name,"NX");
    property.NX = (double)read_parameter(name);
    NX = (int)property.NX;
    sprintf(name,"NY");
    property.NY = (double)read_parameter(name);
    NY = (int)property.NY;
    sprintf(name,"NZ");
    property.NZ = (double)read_parameter(name);
    NZ = (int)property.NZ;
    fprintf(stderr,"System Size:\nNX %d \nNY %d \nNZ %d\n", NX , NY, NZ);


    /* Size of simulation domain */
    sprintf(name,"SX");
    property.SX = (double)read_parameter(name);
    sprintf(name,"SY");
    property.SY = (double)read_parameter(name);
    sprintf(name,"SZ");
    property.SZ = (double)read_parameter(name);
    fprintf(stderr,"System Size:\nSX %d \nSY %d \nSZ %d\n", (int)property.SX , (int)property.SY, (int)property.SZ);
#ifdef METHOD_STREAMING
    /* check if \Delta x = 1 */
    //if(property.SX != property.NX ){ fprintf(stderr," WARNING!  property.SX != property.NX with STREAMING ON, please change it.\n"); exit(-1);}
    //if(property.SY != property.NY ){ fprintf(stderr," WARNING!  property.SY != property.NY with STREAMING ON, please change it.\n"); exit(-1);}
    //if(property.SZ != property.NZ ){ fprintf(stderr," WARNING!  property.SZ != property.NZ with STREAMING ON, please change it.\n"); exit(-1);}
    fprintf(stderr," Delta X = %e \n",property.SX/property.NX);
    fprintf(stderr," Delta Y = %e \n",property.SY/property.NY);
    fprintf(stderr," Delta Z = %e \n",property.SZ/property.NZ);
    //if(property.SX/property.NX != property.SY/property.NY || property.SX/property.NX != property.SZ/property.NZ || property.SY/property.NY != property.SZ/property.NZ )
    //{ fprintf(stderr," WARNING!  DX != DY != DZ with STREAMING ON, please change it.\n"); exit(-1);}
#endif

    /* time stepping parameters */
    sprintf(name,"time_dt");
    property.time_dt = (double)read_parameter(name); 
    fprintf(stderr," Time step: %g\n",property.time_dt);
#ifdef METHOD_STREAMING
    /* check if \Delta t = 1 */
    //if(property.time_dt != 1.0 ){ fprintf(stderr," WARNING! property.time_dt != 1 with STREAMING ON, please change it.\n"); exit(-1);}
   if(property.time_dt != property.SX/property.NX ){ fprintf(stderr," WARNING! property.time_dt != Delta X or Y or Z, with STREAMING ON, please change it.\n"); exit(-1);}

#endif
    sprintf(name,"time_start");
    property.time_start = (double)read_parameter(name); 
    fprintf(stderr,"Time start: %g\n",property.time_start);

    sprintf(name,"time_end");
    property.time_end = (double)read_parameter(name); 
    fprintf(stderr,"Time end: %g\n",property.time_end);
  
    sprintf(name,"time_dump_field");
    property.time_dump_field = (double)read_parameter(name);
    fprintf(stderr,"Time dump fields: %g\n",property.time_dump_field);
    if(property.time_dump_field < property.time_dt){ fprintf(stderr," WARNING! property.time_dump_field < property.time_dt , please change it.\n"); exit(-1);}

    sprintf(name,"time_dump_diagn");
    property.time_dump_diagn = (double)read_parameter(name);
    fprintf(stderr,"Time dump diagnostic: %g\n",property.time_dump_diagn);
    if(property.time_dump_diagn < property.time_dt){ fprintf(stderr," WARNING! property.time_dump_diagn < property.time_dt , please change it.\n"); exit(-1);}

#ifdef LAGRANGE
    sprintf(name,"time_dump_lagr");
    property.time_dump_lagr = (double)read_parameter(name);
    fprintf(stderr,"Time dump lagrange: %g\n",property.time_dump_lagr);
    if(property.time_dump_lagr < property.time_dt){ fprintf(stderr," WARNING! property.time_dump_lagr < property.time_dt , please change it.\n"); exit(-1);}
#endif

#ifdef LB_FLUID
    /* relaxation time and viscosity for fluid */
  fprintf(stderr,"YES <- LB_FLUID\n");
  sprintf(name,"tau_u");
  property.tau_u = read_parameter(name);
  fprintf(stderr,"Properties:\ntau_u %g\n",(double)property.tau_u);
 #ifdef METHOD_STREAMING
  //property.nu = (property.tau_u-0.5)/3.0;
  property.nu = (property.tau_u-0.5*property.time_dt)/3.0;
 #else
  #ifdef METHOD_REDEFINED_POP
  property.nu = (property.tau_u-0.5*property.time_dt)/3.0;
  #else
  property.nu = property.tau_u/3.0;
  #endif
 #endif
  fprintf(stderr,"viscosity %g\n",(double)property.nu);
 #ifdef LB_FLUID_FORCING
  /* forcing Amplitude */
  fprintf(stderr,"YES <- LB_FLUID_FORCING\n");
  sprintf(name,"Amp_x");
  property.Amp_x = read_parameter(name);
  sprintf(name,"Amp_y");
  property.Amp_y = read_parameter(name);
  sprintf(name,"Amp_z");
  property.Amp_z = read_parameter(name);
  fprintf(stderr,"Properties: Amp_x %g\n",(double)property.Amp_x);
  fprintf(stderr,"Properties: Amp_y %g\n",(double)property.Amp_y);
  fprintf(stderr,"Properties: Amp_z %g\n",(double)property.Amp_z);
  #ifdef LB_FLUID_FORCING_POISEUILLE
  fprintf(stderr,"Reynolds (Poiseuille) Number is -> Re= %e\n", property.Amp_x * property.SY /property.nu);
  fout = fopen("numbers.dat","a");
  fprintf(fout,"Reynolds (Poiseuille) %e\n", property.Amp_x * property.SY /property.nu);
  fclose(fout); 
  #endif
  #ifdef LB_FLUID_FORCING_CHANNEL
  fprintf(stderr,"Reynolds (turbulent shear) Number is -> Re= %e\n", property.Amp_x * (property.SY/2.) /property.nu);
  fout = fopen("numbers.dat","a");
  fprintf(fout,"Reynolds (turbulent shear) %e\n", property.Amp_x * (property.SY/2.) /property.nu);
  fclose(fout); 
  #endif
 #endif
#endif

#ifdef LB_TEMPERATURE
    /* relaxation time and viscosity for temperature */
  fprintf(stderr,"YES <- LB_TEMPERATURE\n");
  sprintf(name,"tau_t");
  property.tau_t = read_parameter(name);
  fprintf(stderr,"Properties:\ntau_t %g\n",(double)property.tau_t);
 #ifdef METHOD_STREAMING
  //property.kappa = (property.tau_t-0.5)/3.0;
  property.kappa = (property.tau_t-0.5*property.time_dt)/3.0;
 #else
  #ifdef METHOD_REDEFINED_POP
  property.kappa = (property.tau_t-0.5*property.time_dt)/3.0;
  #else
  property.kappa = property.tau_t/3.0;
  #endif
 #endif
  fprintf(stderr,"thermal diffusivity %g\n",(double)property.kappa);
  /* wall temperature values */
  sprintf(name,"T_bot");
  property.T_bot = read_parameter(name);
  sprintf(name,"T_top");
  property.T_top = read_parameter(name);
  sprintf(name,"T_ref");
  property.T_ref = read_parameter(name);
  property.deltaT = property.T_bot-property.T_top;
  fprintf(stderr,"T_bot %g , T_top %g , deltaT %g\n",(double)property.T_bot, (double)property.T_top, (double)property.deltaT);
  /* wall temperature gradients */
  sprintf(name,"grad_T_bot");
  property.grad_T_bot = read_parameter(name);
  sprintf(name,"grad_T_top");
  property.grad_T_top = read_parameter(name);
  fprintf(stderr,"grad_T_bot %g , grad_T_top %g\n",(double)property.grad_T_bot, (double)property.grad_T_top);
 #ifdef LB_TEMPERATURE_BUOYANCY
  fprintf(stderr,"YES <- LB_TEMPERATURE_BUOYANCY\n");
  sprintf(name,"beta_t");
  property.beta_t = read_parameter(name);
  fprintf(stderr,"linear volume expansion coefficient %g\n",(double)property.beta_t);
  sprintf(name,"beta2_t");
  property.beta2_t = read_parameter(name);
  fprintf(stderr,"quadratic volume expansion coefficient %g\n",(double)property.beta2_t);
  sprintf(name,"gravity_x");
  property.gravity_x = read_parameter(name);
  sprintf(name,"gravity_y");
  property.gravity_y = read_parameter(name);
  sprintf(name,"gravity_z");
  property.gravity_z = read_parameter(name);
  fprintf(stderr,"gravity_x %g, gravity_y %g, gravity_z %g\n",(double)property.gravity_x, (double)property.gravity_y, (double)property.gravity_z);

  fprintf(stderr,"Rayleigh Number is -> Ra = %e\n", property.beta_t*property.gravity_y*property.deltaT*pow(property.SY,3.0)/(property.nu*property.kappa) );
  fprintf(stderr,"Prandtl Number is -> Pr = %e\n", property.nu/property.kappa);
  fout = fopen("numbers.dat","a");
  fprintf(fout,"Rayleigh %e\n",property.beta_t*property.gravity_y*property.deltaT*pow(property.SY,3.0)/(property.nu*property.kappa) );
  fprintf(fout,"Prandtl %e\n", property.nu/property.kappa);   
  fclose(fout);
 #endif 
 #ifdef LB_TEMPERATURE_FORCING
  sprintf(name,"Amp_t");
  property.Amp_t = read_parameter(name);
  #ifdef LB_TEMPERATURE_FORCING_BULK
   #ifdef LB_TEMPERATURE_BUOYANCY
  fprintf(stderr,"Internal Rayleigh Number is -> Ra_{Int} = %e\n", property.beta_t*property.gravity_y*property.Amp_t*pow(property.SY,5.0)/(property.nu*pow(property.kappa,2.0)) ); 
  fout = fopen("numbers.dat","a");
  fprintf(fout,"Internal Rayleigh %e\n",property.beta_t*property.gravity_y*property.Amp_t*pow(property.SY,5.0)/(property.nu*pow(property.kappa,2.0)) );
  fclose(fout);
   #endif
  #endif
  #ifdef LB_TEMPERATURE_FORCING_RADIATION
  sprintf(name,"attenuation");
  property.attenuation = read_parameter(name);
  fprintf(stderr,"Bo Number is -> Bo = %e\n", property.Amp_t/(property.kappa*property.deltaT/property.SY) ); 
  fprintf(stderr,"Er Number is -> Er = %e\n", property.attenuation*property.SY );
  //  fprintf(stderr,"Radiation Rayleigh Number is -> Ra_{rad} = %e\n", property.attenuation*property.beta_t*property.gravity_y*property.Amp_t*pow(property.SY,5.0)/(property.nu*pow(property.kappa,2.0)) );  
  fprintf(stderr,"Radiation Rayleigh Number (new definition) is -> Ra_{rad} = %e\n", property.beta_t*property.gravity_y*property.Amp_t*pow(property.SY,4.0)/(property.nu*pow(property.kappa,2.0)) );
  fout = fopen("numbers.dat","a");
  fprintf(fout,"Bo Number is -> Bo = %e\n", property.Amp_t/(property.kappa*property.deltaT/property.SY) ); 
  fprintf(fout,"Er Number is -> Er = %e\n", property.attenuation*property.SY ); 
  fprintf(fout,"Radiation Rayleigh Number (new definition) is -> Ra_{rad} = %e\n", property.attenuation*property.beta_t*property.gravity_y*property.Amp_t*pow(property.SY,5.0)/(property.nu*pow(property.kappa,2.0)) );  
  fclose(fout);
  #endif
  #ifdef LB_TEMPERATURE_MELTING
  sprintf(name,"T_solid");
  property.T_solid = read_parameter(name);
  sprintf(name,"specific_heat");
  property.specific_heat = read_parameter(name);
  sprintf(name,"latent_heat");
  property.latent_heat = read_parameter(name);
  fprintf(stderr,"Stefan Number is -> Ste = %e\n", property.deltaT*property.specific_heat/property.latent_heat);
  fout = fopen("numbers.dat","a");
  fprintf(fout,"Stefan = %e\n", property.deltaT*property.specific_heat/property.latent_heat);
  fclose(fout);
  #endif
 #endif /* end of ifdef LB_TEMPERATURE_FORCING */
#endif /* end of ifdef LB_TEMPERATURE */

#ifdef LB_SCALAR
    /* relaxation time and viscosity for temperature */
  fprintf(stderr,"YES <- LB_SCALAR\n");
  sprintf(name,"tau_s");
  property.tau_s = read_parameter(name);
  fprintf(stderr,"Properties:\ntau_s %g\n",(double)property.tau_s);
 #ifdef METHOD_STREAMING
  //property.chi = (property.tau_s-0.5)/3.0;
  property.chi = (property.tau_s-0.5*property.time_dt)/3.0;
 #else
  #ifdef METHOD_REDEFINED_POP
  property.chi = (property.tau_s-0.5*property.time_dt)/3.0;
  #else
  property.chi = property.tau_s/3.0;
  #endif
 #endif
  fprintf(stderr,"mass diffusivity %g\n",(double)property.chi);

  sprintf(name,"S_bot");
  property.S_bot = read_parameter(name);
  sprintf(name,"S_top");
  property.S_top = read_parameter(name);
  sprintf(name,"S_ref");
  property.S_ref = read_parameter(name);
  property.deltaS = property.S_bot-property.S_top;
  fprintf(stderr,"S_bot %g , S_top %g , deltaS %g\n",(double)property.S_bot, (double)property.S_top, (double)property.deltaS);
 #ifdef LB_SCALAR_BUOYANCY
  fprintf(stderr,"YES <- LB_SCALAR_BUOYANCY\n");
  sprintf(name,"beta_s");
  property.beta_s = read_parameter(name);
  fprintf(stderr,"linear volume expansion coefficient for scalar %g\n",(double)property.beta_s);
 #endif
 #ifdef LB_SCALAR_FORCING
  sprintf(name,"Amp_s");
  property.Amp_s = read_parameter(name);
  fprintf(stderr,"Amplitude forcing scalar -> Amp_s = %e\n", property.Amp_s ); 
 #endif
#endif /* end of ifded LB_SCALAR */


#ifdef LAGRANGE
    /* total number of particles or tracers */
  fprintf(stderr,"YES <- LAGRANGE\n");
  sprintf(name,"particle_number");
  property.particle_number = read_parameter(name);
  fprintf(stderr,"Properties: particle_number %g\n",(double)property.particle_number);

  /* drag response time tau_drag */  
  sprintf(name,"tau_drag_types");
  property.tau_drag_types = read_parameter(name); 
  sprintf(name,"tau_drag_min");
  property.tau_drag_min = read_parameter(name);
  sprintf(name,"tau_drag_max");
  property.tau_drag_max = read_parameter(name);
  fprintf(stderr,"tau_drag_types %g , tau_drag_max %g , tau_drag_min %g\n",(double)property.tau_drag_types, (double)property.tau_drag_min, (double)property.tau_drag_max);
  if( property.tau_drag_types <1 || property.tau_drag_max < property.tau_drag_min ){ fprintf(stderr,"Error in tau_drag parameters\n Exit.\n"); exit(0);}

  /* total number of particles types up to now */
  property.particle_types = property.tau_drag_types;

 #ifdef LAGRANGE_GRADIENT
  #ifdef LAGRANGE_ADDEDMASS
  /* added mass  */  
  sprintf(name,"beta_coeff_types");
  property.beta_coeff_types = read_parameter(name); 
  sprintf(name,"beta_coeff_min");
  property.beta_coeff_min = read_parameter(name);
  sprintf(name,"beta_coeff_max");
  property.beta_coeff_max = read_parameter(name);
  fprintf(stderr,"beta_coeff_types %g , beta_coeff_max %g , beta_coeff_min %g\n",(double)property.beta_coeff_types, (double)property.beta_coeff_min, (double)property.beta_coeff_max);
  if( property.beta_coeff_types <1 || property.beta_coeff_max < property.beta_coeff_min ){ fprintf(stderr,"Error in beta_coeff parameters\n Exit.\n"); exit(0);}

  /* total number of particles types up to now */
  property.particle_types *= property.beta_coeff_types;
  #endif

  #ifdef LAGRANGE_ORIENTATION
   #ifdef LAGRANGE_ORIENTATION_JEFFREY
   /* aspect ratio for rotation  */ 
   sprintf(name,"aspect_ratio_types");
   property.aspect_ratio_types = read_parameter(name); 
   sprintf(name,"aspect_ratio_min");
   property.aspect_ratio_min = read_parameter(name);
   sprintf(name,"aspect_ratio_max");
   property.aspect_ratio_max = read_parameter(name);
   fprintf(stderr,"aspect_ratio_types %g , aspect_ratio_max %g , aspect_ratio_min %g\n",(double)property.aspect_ratio_types, (double)property.aspect_ratio_min, (double)property.aspect_ratio_max);
   if( property.aspect_ratio_types <1 || property.aspect_ratio_max < property.aspect_ratio_min ){ fprintf(stderr,"Error in aspect_ratio parameters\n Exit.\n"); exit(0);}   
  
   /* total number of particles types up to now */
   property.particle_types *= property.aspect_ratio_types;
 
    #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
    /* gyrotaxis_velocity for rotation  */ 
    sprintf(name,"gyrotaxis_velocity_types");
    property.gyrotaxis_velocity_types = read_parameter(name); 
    sprintf(name,"gyrotaxis_velocity_min");
    property.gyrotaxis_velocity_min = read_parameter(name);
    sprintf(name,"gyrotaxis_velocity_max");
    property.gyrotaxis_velocity_max = read_parameter(name);
    fprintf(stderr,"gyrotaxis_velocity_types %g , gyrotaxis_velocity_max %g , gyrotaxis_velocity_min %g\n",(double)property.gyrotaxis_velocity_types, (double)property.gyrotaxis_velocity_min, (double)property.gyrotaxis_velocity_max);
    if( property.gyrotaxis_velocity_types <1 || property.gyrotaxis_velocity_max < property.gyrotaxis_velocity_min ){ fprintf(stderr,"Error in gyrotaxis_velocity parameters\n Exit.\n"); exit(0);}   
  
    /* total number of particles types up to now */
    property.particle_types *= property.gyrotaxis_velocity_types;
    #endif
   #endif 

   #ifdef LAGRANGE_ORIENTATION_DIFFUSION
   /* rotational_diffusion for rotation  */ 
   sprintf(name,"rotational_diffusion_types");
   property.rotational_diffusion_types = read_parameter(name); 
   sprintf(name,"rotational_diffusion_min");
   property.rotational_diffusion_min = read_parameter(name);
   sprintf(name,"rotational_diffusion_max");
   property.rotational_diffusion_max = read_parameter(name);
   fprintf(stderr,"rotational_diffusion_types %g , rotational_diffusion_max %g , rotational_diffusion_min %g\n",(double)property.rotational_diffusion_types, (double)property.rotational_diffusion_min, (double)property.rotational_diffusion_max);
   if( property.rotational_diffusion_types <1 || property.rotational_diffusion_max < property.rotational_diffusion_min ){ fprintf(stderr,"Error in rotational_diffusion parameters\n Exit.\n"); exit(0);}   
  
  /* total number of particles types up to now */
   property.particle_types *= property.rotational_diffusion_types;
   #endif

   #ifdef LAGRANGE_ORIENTATION_ACTIVE
  /* swim_velocity for rotation  */ 
  sprintf(name,"swim_velocity_types");
  property.swim_velocity_types = read_parameter(name); 
  sprintf(name,"swim_velocity_min");
  property.swim_velocity_min = read_parameter(name);
  sprintf(name,"swim_velocity_max");
  property.swim_velocity_max = read_parameter(name);
  fprintf(stderr,"swim_velocity_types %g , swim_velocity_max %g , swim_velocity_min %g\n",(double)property.swim_velocity_types, (double)property.swim_velocity_min, (double)property.swim_velocity_max);
  if( property.swim_velocity_types <1 || property.swim_velocity_max < property.swim_velocity_min ){ fprintf(stderr,"Error in swim_velocity parameters\n Exit.\n"); exit(0);}   
  
  /* total number of particles types up to now */
  property.particle_types *= property.swim_velocity_types;
    #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP

  /* critical_shear_rate */ 
  sprintf(name,"critical_shear_rate_types");
  property.critical_shear_rate_types = read_parameter(name); 
  sprintf(name,"critical_shear_rate_min");
  property.critical_shear_rate_min = read_parameter(name);
  sprintf(name,"critical_shear_rate_max");
  property.critical_shear_rate_max = read_parameter(name);
  fprintf(stderr,"critical_shear_rate_types %g , critical_shear_rate_max %g , critical_shear_rate_min %g\n",(double)property.critical_shear_rate_types, (double)property.critical_shear_rate_min, (double)property.critical_shear_rate_max);
  if( property.critical_shear_rate_types <1 || property.critical_shear_rate_max < property.critical_shear_rate_min ){ fprintf(stderr,"Error in critical_shear_rate parameters\n Exit.\n"); exit(0);}   
  
  /* total number of particles types up to now */
  property.particle_types *= property.critical_shear_rate_types;

  /* jump_time */ 
  sprintf(name,"jump_time_types");
  property.jump_time_types = read_parameter(name); 
  sprintf(name,"jump_time_min");
  property.jump_time_min = read_parameter(name);
  sprintf(name,"jump_time_max");
  property.jump_time_max = read_parameter(name);
  fprintf(stderr,"jump_time_types %g , jump_time_max %g , jump_time_min %g\n",(double)property.jump_time_types, (double)property.jump_time_min, (double)property.jump_time_max);
  if( property.jump_time_types <1 || property.jump_time_max < property.jump_time_min ){ fprintf(stderr,"Error in jump_time parameters\n Exit.\n"); exit(0);}   
  
  /* total number of particles types up to now */
  property.particle_types *= property.jump_time_types;

    #endif /* LAGRANGE_ORIENTATION_ACTIVE_JUMP  */
   #endif /* LAGRANGE_ORIENTATION_ACTIVE */
  #endif /* LAGRANGE_ORIENTATION */
 #endif /* LAGRANGE_GRADIENT */

  /* total number of particles types */
   fprintf(stderr,"Properties: particle_types %g\n",(double)property.particle_types);   
#endif /* LAGRANGE */

  /* size of types, just for a check */
    fprintf(stderr,"Size of float %d\n",sizeof(float));
    fprintf(stderr,"Size of double %d\n",sizeof(double));
    fprintf(stderr,"Size of long double %d\n",sizeof(long double));
    fprintf(stderr,"Size of my_double %d\n",sizeof(my_double));
    fprintf(stderr,"\n");
  }/* end of if ROOT*/

 /* Now broadcast all properties */
 MPI_Bcast(&resume, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef LB_FLUID
 MPI_Bcast(&resume_u, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef LB_TEMPERATURE
 MPI_Bcast(&resume_t, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef LB_SCALAR
 MPI_Bcast(&resume_s, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

 MPI_Bcast(&NX, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&NY, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&property, 1,MPI_property_type, 0, MPI_COMM_WORLD);
 

#ifdef DEBUG
 for(i=0;i<nprocs;i++){
   if(i==me){ fprintf(stderr,"me %d , System Size: NX %d NY %d NZ %d\n", me , NX , NY, NZ);
             fprintf(stderr,"me %d , time_dt %g\n", me , property.time_dt);
   }
 }
#endif

}/* end of assign_parameters */


void set_to_zero_vector( vector *f,int size){
  int i;
  for(i=0;i<size;i++) f[i].x = f[i].y = f[i].z = 0.0;
}

void set_to_zero_scalar( my_double *f,int size){
  int i;
  for(i=0;i<size;i++) f[i] = 0.0;
}

void set_to_zero_int( int *f,int size){
  int i;
  for(i=0;i<size;i++) f[i] = 0;
}

void set_to_zero_my_double( my_double *f,int size){
  int i;
  for(i=0;i<size;i++) f[i] = 0.0;
}

void set_to_zero_pop(pop  *f,int size){
  int i,pp;
  for(i=0;i<size;i++)
    for (pp=0; pp<NPOP; pp++)
      f[i].p[pp] = 0.0;
}

void set_to_zero_output(output  *f,int size){
  int i,pp;
  for(i=0;i<size;i++){
#ifdef LB_FLUID
     f[i].x = f[i].y = f[i].z = 0.0; 
     f[i].ux = f[i].uy = f[i].uz = 0.0; 
     f[i].ux2 = f[i].uy2 = f[i].uz2 = 0.0; 
     f[i].rho = f[i].ene = f[i].eps = 0.0;
     f[i].rho2 = 0.0;
#endif
#ifdef LB_TEMPERATURE
    f[i].dxt = f[i].dyt = f[i].dzt = 0.0; 
    f[i].uxt = f[i].uyt = f[i].uzt = 0.0; 
    f[i].nux = f[i].nuy = f[i].nuz = 0.0; 
    f[i].t = f[i].t2 = f[i].epst = f[i].lb = 0.0;
#ifdef LB_TEMPERATURE_MELTING
    f[i].lf = f[i].dtlf = f[i].enth = 0.0; 
#endif
#endif
#ifdef LB_SCALAR
    f[i].dxs = f[i].dys = f[i].dzs = 0.0; 
    f[i].uxs = f[i].uys = f[i].uzs = 0.0; 
    f[i].nusx = f[i].nusy = f[i].nusz = 0.0; 
    f[i].s = f[i].s2 = f[i].epss = 0.0;
#endif
  }
}


void allocate_fields(){
 mesh  = (vector*) malloc(sizeof(vector)*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 if(mesh == NULL){ fprintf(stderr,"Not enough memory to allocate mesh\n"); exit(-1);}
 set_to_zero_vector( mesh,(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD));

 /*
 mesh_flag  = (int*) malloc(sizeof(int)*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 if(mesh_flag == NULL){ fprintf(stderr,"Not enough memory to allocate mesh_flag\n"); exit(-1);}
 set_to_zero_int( mesh_flag,(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD));
 */

 center_V = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(center_V == NULL){ fprintf(stderr,"Not enough memory to allocate center_V\n"); exit(-1);}
 set_to_zero_vector( center_V,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));


#ifdef GRID_REFINED
 grid_ruler_x  = (my_double*) malloc(sizeof(my_double)*NXG);
 grid_ruler_y  = (my_double*) malloc(sizeof(my_double)*NYG);
 grid_ruler_z  = (my_double*) malloc(sizeof(my_double)*NZG);
#endif

#ifdef LB
 coeff_xp = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_xp == NULL){ fprintf(stderr,"Not enough memory to allocate coeff_xp\n"); exit(-1);}
 set_to_zero_pop( coeff_xp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 coeff_xm = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_xm == NULL){ fprintf(stderr,"Not enough memory to allocate coeff_xm\n"); exit(-1);}
 set_to_zero_pop( coeff_xm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 coeff_yp = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_yp == NULL){ fprintf(stderr,"Not enough memory to allocate coeff_yp\n"); exit(-1);}
 set_to_zero_pop( coeff_yp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 coeff_ym = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_ym == NULL){ fprintf(stderr,"Not enough memory to allocate coeff_ym\n"); exit(-1);}
 set_to_zero_pop( coeff_ym,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 coeff_zp = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_zp == NULL){ fprintf(stderr,"Not enough memory to allocate coeff_zp\n"); exit(-1);}
 set_to_zero_pop( coeff_zp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 coeff_zm = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_zm == NULL){ fprintf(stderr,"Not enough memory to allocate coeff_zm\n"); exit(-1);}
 set_to_zero_pop( coeff_zm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

 
#if (defined METHOD_CENTERED || defined METHOD_MYQUICK || defined METHOD_STREAMING || defined METHOD_UPWIND)
 /* note that the first block of fields here below is not necessary for upwind method*/
 interp_xp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp_xp == NULL){ fprintf(stderr,"Not enough memory to allocate interp_xp\n"); exit(-1);}
 set_to_zero_my_double( interp_xp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp_xm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp_xm == NULL){ fprintf(stderr,"Not enough memory to allocate interp_xm\n"); exit(-1);}
 set_to_zero_my_double( interp_xm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp_yp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp_yp == NULL){ fprintf(stderr,"Not enough memory to allocate interp_yp\n"); exit(-1);}
 set_to_zero_my_double( interp_yp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp_ym = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp_ym == NULL){ fprintf(stderr,"Not enough memory to allocate interp_ym\n"); exit(-1);}
 set_to_zero_my_double( interp_ym,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp_zp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp_zp == NULL){ fprintf(stderr,"Not enough memory to allocate interp_zp\n"); exit(-1);}
 set_to_zero_my_double( interp_zp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp_zm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp_zm == NULL){ fprintf(stderr,"Not enough memory to allocate interp_zm\n"); exit(-1);}
 set_to_zero_my_double( interp_zm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 /* borders my_double */
 xp_scalar  = (my_double*) malloc(sizeof(my_double)*BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 xm_scalar  = (my_double*) malloc(sizeof(my_double)*BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 if(xp_scalar == NULL || xm_scalar == NULL){ fprintf(stderr,"Not enough memory to allocate x{p,m}_scalar\n"); exit(-1);}
 yp_scalar  = (my_double*) malloc(sizeof(my_double)*BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
 ym_scalar  = (my_double*) malloc(sizeof(my_double)*BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(yp_scalar == NULL || ym_scalar == NULL){ fprintf(stderr,"Not enough memory to allocate y{p,m}_scalar\n"); exit(-1);}
 zp_scalar  = (my_double*) malloc(sizeof(my_double)*BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
 zm_scalar  = (my_double*) malloc(sizeof(my_double)*BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
 if(zp_scalar == NULL || zm_scalar == NULL){ fprintf(stderr,"Not enough memory to allocate z{p,m}_scalar\n"); exit(-1);}

/* borders vector */
 xp_vector  = (vector*) malloc(sizeof(vector)*BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 xm_vector  = (vector*) malloc(sizeof(vector)*BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 if(xp_vector == NULL || xm_vector == NULL){ fprintf(stderr,"Not enough memory to allocate x{p,m}_vector\n"); exit(-1);}
 yp_vector  = (vector*) malloc(sizeof(vector)*BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
 ym_vector  = (vector*) malloc(sizeof(vector)*BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(yp_vector == NULL || ym_vector == NULL){ fprintf(stderr,"Not enough memory to allocate y{p,m}_vector\n"); exit(-1);}
 zp_vector  = (vector*) malloc(sizeof(vector)*BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
 zm_vector  = (vector*) malloc(sizeof(vector)*BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
 if(zp_vector == NULL || zm_vector == NULL){ fprintf(stderr,"Not enough memory to allocate z{p,m}_vector\n"); exit(-1);}
#endif
 
 //#ifdef METHOD_MYQUICK
#if (defined METHOD_MYQUICK || defined METHOD_UPWIND_SKEW)
 interp2_xp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp2_xp == NULL){ fprintf(stderr,"Not enough memory to allocate interp2_xp\n"); exit(-1);}
 set_to_zero_my_double( interp2_xp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp2_xm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp2_xm == NULL){ fprintf(stderr,"Not enough memory to allocate interp2_xm\n"); exit(-1);}
 set_to_zero_my_double( interp2_xm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp2_yp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp2_yp == NULL){ fprintf(stderr,"Not enough memory to allocate interp2_yp\n"); exit(-1);}
 set_to_zero_my_double( interp2_yp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp2_ym = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp2_ym == NULL){ fprintf(stderr,"Not enough memory to allocate interp2_ym\n"); exit(-1);}
 set_to_zero_my_double( interp2_ym,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp2_zp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp2_zp == NULL){ fprintf(stderr,"Not enough memory to allocate interp2_zp\n"); exit(-1);}
 set_to_zero_my_double( interp2_zp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp2_zm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp2_zm == NULL){ fprintf(stderr,"Not enough memory to allocate interp2_zm\n"); exit(-1);}
 set_to_zero_my_double( interp2_zm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 interp3_xp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp3_xp == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_xp\n"); exit(-1);}
 set_to_zero_my_double( interp3_xp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp3_xm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp3_xm == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_xm\n"); exit(-1);}
 set_to_zero_my_double( interp3_xm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp3_yp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp3_yp == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_yp\n"); exit(-1);}
 set_to_zero_my_double( interp3_yp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp3_ym = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp3_ym == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_ym\n"); exit(-1);}
 set_to_zero_my_double( interp3_ym,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp3_zp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp3_zp == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_zp\n"); exit(-1);}
 set_to_zero_my_double( interp3_zp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp3_zm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp3_zm == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_zm\n"); exit(-1);}
 set_to_zero_my_double( interp3_zm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 interp4_xp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp4_xp == NULL){ fprintf(stderr,"Not enough memory to allocate interp4_xp\n"); exit(-1);}
 set_to_zero_my_double( interp4_xp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp4_xm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp4_xm == NULL){ fprintf(stderr,"Not enough memory to allocate interp4_xm\n"); exit(-1);}
 set_to_zero_my_double( interp4_xm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp4_yp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp4_yp == NULL){ fprintf(stderr,"Not enough memory to allocate interp4_yp\n"); exit(-1);}
 set_to_zero_my_double( interp4_yp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp4_ym = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp4_ym == NULL){ fprintf(stderr,"Not enough memory to allocate interp4_ym\n"); exit(-1);}
 set_to_zero_my_double( interp4_ym,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp4_zp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp4_zp == NULL){ fprintf(stderr,"Not enough memory to allocate interp4_zp\n"); exit(-1);}
 set_to_zero_my_double( interp4_zp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp4_zm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp4_zm == NULL){ fprintf(stderr,"Not enough memory to allocate interp4_zm\n"); exit(-1);}
 set_to_zero_my_double( interp4_zm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

#ifdef METHOD_UPWIND_LINEAR
interp5_xp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp5_xp == NULL){ fprintf(stderr,"Not enough memory to allocate interp5_xp\n"); exit(-1);}
 set_to_zero_my_double( interp5_xp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp5_xm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp5_xm == NULL){ fprintf(stderr,"Not enough memory to allocate interp5_xm\n"); exit(-1);}
 set_to_zero_my_double( interp5_xm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp5_yp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp5_yp == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_yp\n"); exit(-1);}
 set_to_zero_my_double( interp5_yp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp5_ym = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp5_ym == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_ym\n"); exit(-1);}
 set_to_zero_my_double( interp5_ym,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp5_zp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp5_zp == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_zp\n"); exit(-1);}
 set_to_zero_my_double( interp5_zp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp5_zm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp5_zm == NULL){ fprintf(stderr,"Not enough memory to allocate interp3_zm\n"); exit(-1);}
 set_to_zero_my_double( interp5_zm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

interp6_xp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp6_xp == NULL){ fprintf(stderr,"Not enough memory to allocate interp6_xp\n"); exit(-1);}
 set_to_zero_my_double( interp6_xp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp6_xm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp6_xm == NULL){ fprintf(stderr,"Not enough memory to allocate interp6_xm\n"); exit(-1);}
 set_to_zero_my_double( interp6_xm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp6_yp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp6_yp == NULL){ fprintf(stderr,"Not enough memory to allocate interp6_yp\n"); exit(-1);}
 set_to_zero_my_double( interp6_yp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp6_ym = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp6_ym == NULL){ fprintf(stderr,"Not enough memory to allocate interp6_ym\n"); exit(-1);}
 set_to_zero_my_double( interp6_ym,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp6_zp = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp6_zp == NULL){ fprintf(stderr,"Not enough memory to allocate interp6_zp\n"); exit(-1);}
 set_to_zero_my_double( interp6_zp,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 interp6_zm = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(interp6_zm == NULL){ fprintf(stderr,"Not enough memory to allocate interp6_zm\n"); exit(-1);}
 set_to_zero_my_double( interp6_zm,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif


#ifdef LB_FLUID
 p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(p == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 set_to_zero_pop( p,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 rhs_p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(rhs_p == NULL){ fprintf(stderr,"Not enough memory to allocate rhs_p\n"); exit(-1);}
 set_to_zero_pop( rhs_p,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

#ifdef METHOD_REDEFINED_POP
 /* auxiliary population field for computing the advection term */
 f_aux  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(f_aux == NULL){ fprintf(stderr,"Not enough memory to allocate f_aux\n"); exit(-1);}
 set_to_zero_pop( f_aux,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

#if (defined METHOD_REDEFINED_POP || defined METHOD_COLLISION_IMPLICIT)
 p_eq  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(p_eq == NULL){ fprintf(stderr,"Not enough memory to allocate p_eq\n"); exit(-1);}
 set_to_zero_pop( p_eq,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

#if (defined METHOD_STEPPING_AB2 || defined METHOD_STEPPING_AB3 || defined METHOD_HEUN)
 old_rhs_p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_rhs_p == NULL){ fprintf(stderr,"Not enough memory to allocate old_rhs_p\n"); exit(-1);}
 set_to_zero_pop( old_rhs_p,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

#ifdef METHOD_STEPPING_AB3
 old_old_rhs_p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_old_rhs_p == NULL){ fprintf(stderr,"Not enough memory to allocate old_old_rhs_p\n"); exit(-1);}
 set_to_zero_pop(old_old_rhs_p,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

 u = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(u == NULL){ fprintf(stderr,"Not enough memory to allocate u\n"); exit(-1);}
 set_to_zero_vector( u,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 dens  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(dens == NULL){ fprintf(stderr,"Not enough memory to allocate dens\n"); exit(-1);}
 set_to_zero_my_double( dens,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 #ifdef LB_FLUID_PAST
 old_u = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_u == NULL){ fprintf(stderr,"Not enough memory to allocate old_u\n"); exit(-1);}
 set_to_zero_vector( old_u,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 old_dens  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_dens == NULL){ fprintf(stderr,"Not enough memory to allocate old_dens\n"); exit(-1);}
 set_to_zero_my_double( old_dens,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 #endif

 #ifdef LB_FLUID_LES
  #ifdef LB_FLUID_LES_SISM 
  u_mean = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
  if(u_mean == NULL){ fprintf(stderr,"Not enough memory to allocate u_mean\n"); exit(-1);}
  set_to_zero_vector( u_mean,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
   #ifdef LB_FLUID_LES_SISM_KALMAN
   u_mean_kalman_pre = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
   if(u_mean_kalman_pre == NULL){ fprintf(stderr,"Not enough memory to allocate  u_mean_kalman_pre\n"); exit(-1);}
   set_to_zero_vector( u_mean_kalman_pre,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

   sqr_var = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
   if(sqr_var == NULL){ fprintf(stderr,"Not enough memory to allocate  sqr_var\n"); exit(-1);}
   set_to_zero_vector( sqr_var,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

   sqr_var_kalman = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
   if(sqr_var_kalman == NULL){ fprintf(stderr,"Not enough memory to allocate sqr_var_kalman\n"); exit(-1);}
   set_to_zero_vector( sqr_var_kalman,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

   K_kalman = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
   if(K_kalman == NULL){ fprintf(stderr,"Not enough memory to allocate  K_kalman\n"); exit(-1);}
   set_to_zero_vector( K_kalman,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

   P_kalman = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
   if(P_kalman == NULL){ fprintf(stderr,"Not enough memory to allocate  P_kalman\n"); exit(-1);}
   set_to_zero_vector(P_kalman,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

   P_kalman_pre = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
   if(P_kalman_pre == NULL){ fprintf(stderr,"Not enough memory to allocate P_kalman_pre\n"); exit(-1);}
   set_to_zero_vector( P_kalman_pre,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
   #endif
  #endif
 #endif

#ifdef LB_FLUID_FORCING
 force  = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(force == NULL){ fprintf(stderr,"Not enough memory to allocate force\n"); exit(-1);}
 set_to_zero_vector( force,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 #ifdef LB_FLUID_FORCING_LANDSCAPE
 landscape  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(landscape == NULL){ fprintf(stderr,"Not enough memory to allocate dens\n"); exit(-1);}
 set_to_zero_my_double( landscape,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 #endif
#endif

 /* borders pop */
 xp_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 xm_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 if(xp_pop == NULL || xm_pop == NULL){ fprintf(stderr,"Not enough memory to allocate x{p,m}_pop\n"); exit(-1);}
 yp_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
 ym_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(yp_pop == NULL || ym_pop == NULL){ fprintf(stderr,"Not enough memory to allocate y{p,m}_pop\n"); exit(-1);}
 zp_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
 zm_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
 if(zp_pop == NULL || zm_pop == NULL){ fprintf(stderr,"Not enough memory to allocate z{p,m}_pop\n"); exit(-1);}


#ifdef METHOD_EDGES_AND_CORNERS
/* 8 corners */
 xp_yp_zp_corner_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*BRD);
 xp_yp_zm_corner_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*BRD);
 xp_ym_zp_corner_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*BRD);
 xp_ym_zm_corner_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*BRD);
 xm_yp_zp_corner_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*BRD);
 xm_yp_zm_corner_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*BRD);
 xm_ym_zp_corner_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*BRD);
 xm_ym_zm_corner_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*BRD);

/* 12 edges */
xp_yp_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNZ+TWO_BRD)); 
xp_ym_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNZ+TWO_BRD));
xm_yp_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNZ+TWO_BRD));
xm_ym_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNZ+TWO_BRD));

xp_zp_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNY+TWO_BRD));
xp_zm_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNY+TWO_BRD));
xm_zp_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNY+TWO_BRD));
xm_zm_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNY+TWO_BRD));

yp_zp_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNX+TWO_BRD));
yp_zm_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNX+TWO_BRD));
ym_zp_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNX+TWO_BRD));
ym_zm_edge_pop = (pop*) malloc(sizeof(pop)*BRD*BRD*(LNX+TWO_BRD));

/* and for vector */
/* 8 corners */
 xp_yp_zp_corner_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*BRD);
 xp_yp_zm_corner_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*BRD);
 xp_ym_zp_corner_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*BRD);
 xp_ym_zm_corner_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*BRD);
 xm_yp_zp_corner_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*BRD);
 xm_yp_zm_corner_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*BRD);
 xm_ym_zp_corner_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*BRD);
 xm_ym_zm_corner_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*BRD);

/* 12 edges */
xp_yp_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNZ+TWO_BRD)); 
xp_ym_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNZ+TWO_BRD));
xm_yp_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNZ+TWO_BRD));
xm_ym_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNZ+TWO_BRD));

xp_zp_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNY+TWO_BRD));
xp_zm_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNY+TWO_BRD));
xm_zp_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNY+TWO_BRD));
xm_zm_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNY+TWO_BRD));

yp_zp_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNX+TWO_BRD));
yp_zm_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNX+TWO_BRD));
ym_zp_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNX+TWO_BRD));
ym_zm_edge_vector = (vector*) malloc(sizeof(vector)*BRD*BRD*(LNX+TWO_BRD));

/* and for scalar */
/* 8 corners */
 xp_yp_zp_corner_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*BRD);
 xp_yp_zm_corner_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*BRD);
 xp_ym_zp_corner_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*BRD);
 xp_ym_zm_corner_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*BRD);
 xm_yp_zp_corner_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*BRD);
 xm_yp_zm_corner_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*BRD);
 xm_ym_zp_corner_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*BRD);
 xm_ym_zm_corner_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*BRD);

/* 12 edges */
xp_yp_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNZ+TWO_BRD)); 
xp_ym_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNZ+TWO_BRD));
xm_yp_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNZ+TWO_BRD));
xm_ym_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNZ+TWO_BRD));

xp_zp_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNY+TWO_BRD));
xp_zm_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNY+TWO_BRD));
xm_zp_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNY+TWO_BRD));
xm_zm_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNY+TWO_BRD));

yp_zp_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNX+TWO_BRD));
yp_zm_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNX+TWO_BRD));
ym_zp_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNX+TWO_BRD));
ym_zm_edge_scalar = (my_double*) malloc(sizeof(my_double)*BRD*BRD*(LNX+TWO_BRD));

#endif


#endif

#ifdef LB_FLUID
 ruler_x_local  = (output*) malloc(sizeof(output)*NX);
 ruler_y_local  = (output*) malloc(sizeof(output)*NY);
 ruler_z_local  = (output*) malloc(sizeof(output)*NZ);
 ruler_x  = (output*) malloc(sizeof(output)*NX);
 ruler_y  = (output*) malloc(sizeof(output)*NY);
 ruler_z  = (output*) malloc(sizeof(output)*NZ);
 ruler_x_running  = (output*) malloc(sizeof(output)*NX);
 ruler_y_running  = (output*) malloc(sizeof(output)*NY);
 ruler_z_running  = (output*) malloc(sizeof(output)*NZ);
 set_to_zero_output(ruler_x_local,NX);
 set_to_zero_output(ruler_y_local,NY);
 set_to_zero_output(ruler_z_local,NZ);
 set_to_zero_output(ruler_x,NX);
 set_to_zero_output(ruler_y,NY);
 set_to_zero_output(ruler_z,NZ);
 set_to_zero_output(ruler_x_running,NX);
 set_to_zero_output(ruler_y_running,NY);
 set_to_zero_output(ruler_z_running,NZ);
#endif

#ifdef LB_TEMPERATURE
 g = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(g == NULL){ fprintf(stderr,"Not enough memory to allocate g\n"); exit(-1);}
 set_to_zero_pop( g,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 rhs_g  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(rhs_g == NULL){ fprintf(stderr,"Not enough memory to allocate rhs_g\n"); exit(-1);}
 set_to_zero_pop( rhs_g,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 //#ifdef METHOD_REDEFINED_POP
#if (defined METHOD_REDEFINED_POP || defined METHOD_COLLISION_IMPLICIT)
 g_eq  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(g_eq == NULL){ fprintf(stderr,"Not enough memory to allocate g_eq\n"); exit(-1);}
 set_to_zero_pop( g_eq,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

#if (defined METHOD_STEPPING_AB2 || defined METHOD_STEPPING_AB3 || defined METHOD_HEUN)
 old_rhs_g  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_rhs_g == NULL){ fprintf(stderr,"Not enough memory to allocate old_rhs_g\n"); exit(-1);}
 set_to_zero_pop( old_rhs_g,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

#ifdef METHOD_STEPPING_AB3
 old_old_rhs_g  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_old_rhs_g == NULL){ fprintf(stderr,"Not enough memory to allocate old_old_rhs_g\n"); exit(-1);}
 set_to_zero_pop( old_old_rhs_g,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

 t  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(t == NULL){ fprintf(stderr,"Not enough memory to allocate t\n"); exit(-1);}
 set_to_zero_my_double( t,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

#ifdef LB_TEMPERATURE_PAST
 old_t  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_t == NULL){ fprintf(stderr,"Not enough memory to allocate old_t\n"); exit(-1);}
 set_to_zero_my_double( old_t,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

 #ifdef LB_TEMPERATURE_FORCING
 t_source  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(t_source == NULL){ fprintf(stderr,"Not enough memory to allocate source_t\n"); exit(-1);}
 set_to_zero_scalar( t_source,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

  #ifdef LB_TEMPERATURE_FORCING_PAST
  old_t_source  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
  if(old_t_source == NULL){ fprintf(stderr,"Not enough memory to allocate old_source_t\n"); exit(-1);}
  set_to_zero_scalar( old_t_source,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
  #endif 
  #ifdef LB_TEMPERATURE_MELTING
  liquid_frac  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
  if(liquid_frac == NULL){ fprintf(stderr,"Not enough memory to allocate liquid_frac\n"); exit(-1);}
  set_to_zero_scalar( liquid_frac,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
  liquid_frac_old  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
  if(liquid_frac_old == NULL){ fprintf(stderr,"Not enough memory to allocate liquid_frac_old\n"); exit(-1);}
  set_to_zero_scalar( liquid_frac_old,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
  //enthalpy  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
  //if(enthalpy == NULL){ fprintf(stderr,"Not enough memory to allocate enthalpy\n"); exit(-1);}
  //set_to_zero_scalar( enthalpy,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
  #endif
 #endif
#endif

#ifdef LB_SCALAR
 h = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(h == NULL){ fprintf(stderr,"Not enough memory to allocate h\n"); exit(-1);}
 set_to_zero_pop( h,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 rhs_h  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(rhs_h == NULL){ fprintf(stderr,"Not enough memory to allocate rhs_h\n"); exit(-1);}
 set_to_zero_pop( rhs_h,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 //#ifdef METHOD_REDEFINED_POP
#if (defined METHOD_REDEFINED_POP || defined METHOD_COLLISION_IMPLICIT)
 h_eq  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(h_eq == NULL){ fprintf(stderr,"Not enough memory to allocate h_eq\n"); exit(-1);}
 set_to_zero_pop( h_eq,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

#if (defined METHOD_STEPPING_AB2 || defined METHOD_STEPPING_AB3 || defined METHOD_HEUN)
 old_rhs_h  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_rhs_h == NULL){ fprintf(stderr,"Not enough memory to allocate old_rhs_h\n"); exit(-1);}
 set_to_zero_pop( old_rhs_h,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif 

#ifdef METHOD_STEPPING_AB3
 old_old_rhs_h  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_old_rhs_h == NULL){ fprintf(stderr,"Not enough memory to allocate old_old_rhs_h\n"); exit(-1);}
 set_to_zero_pop( old_old_rhs_h,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif 

 s  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(s == NULL){ fprintf(stderr,"Not enough memory to allocate s\n"); exit(-1);}
 set_to_zero_my_double( s,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 #ifdef LB_SCALAR_FORCING
 s_source  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(s_source == NULL){ fprintf(stderr,"Not enough memory to allocate s_source\n"); exit(-1);}
 set_to_zero_scalar( s_source,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 #endif
#endif

}



/* Free fields */
void free_fields(){

  free(mesh);
  /*
  free(mesh_flag);
  */
  free(center_V);

#ifdef GRID_REFINED
  free(grid_ruler_x);  
  free(grid_ruler_y);  
  free(grid_ruler_z);  
#endif

#ifdef LB
 free(coeff_xp); 
 free(coeff_xm); 
 free(coeff_yp); 
 free(coeff_ym); 
 free(coeff_zp); 
 free(coeff_zm); 
#endif


#if (defined METHOD_CENTERED || defined METHOD_MYQUICK || defined METHOD_STREAMING || defined METHOD_UPWIND)
 free(interp_xp); 
 free(interp_xm); 
 free(interp_yp); 
 free(interp_ym); 
 free(interp_zp); 
 free(interp_zm); 
 
 /* borders my_double */
 free(xp_scalar);  
 free(xm_scalar);  
 free(yp_scalar);  
 free(ym_scalar);  
 free(zp_scalar);  
 free(zm_scalar);  

/* borders vector */
 free(xp_vector);  
 free(xm_vector);  
 free(yp_vector);  
 free(ym_vector);  
 free(zp_vector);  
 free(zm_vector);  
#endif
 
#ifdef METHOD_MYQUICK
 free(interp2_xp); 
 free(interp2_xm); 
 free(interp2_yp); 
 free(interp2_ym); 
 free(interp2_zp); 
 free(interp2_zm); 

 free(interp3_xp); 
 free(interp3_xm);
 free(interp3_yp);
 free(interp3_ym);
 free(interp3_zp); 
 free(interp3_zm); 

 free(interp4_xp); 
 free(interp4_xm); 
 free(interp4_yp);
 free(interp4_ym); 
 free(interp4_zp); 
 free(interp4_zm); 
#endif



#ifdef LB_FLUID
 free(p); 
 free(rhs_p); 
 free(old_rhs_p); 
 free(old_old_rhs_p); 
 free(u); 
 free(dens); 
 #ifdef LB_FLUID_PAST 
 free(old_u); 
 free(old_dens); 
 #endif
 #ifdef LB_FLUID_LES_SISM 
 free(u_mean);
 #ifdef LB_FLUID_LES_SISM_KALMAN
  free(u_mean_kalman_pre);
  free(sqr_var);
  free(sqr_var_kalman);
  free(K_kalman);
  free(P_kalman);
  free(P_kalman_pre);
  #endif
 #endif
#ifdef LB_FLUID_FORCING
 free(force); 
#ifdef LB_FLUID_FORCING_LANDSCAPE
 free(landscape);
#endif 
#endif

 /* borders pop */
 free(xp_pop); 
 free(xm_pop);
 free(yp_pop); 
 free(ym_pop); 
 free(zp_pop); 
 free(zm_pop); 

#ifdef METHOD_EDGES_AND_CORNERS
/* 8 corners */
 free(xp_yp_zp_corner_pop);
 free(xp_yp_zm_corner_pop);
 free(xp_ym_zp_corner_pop);
 free(xp_ym_zm_corner_pop);
 free(xm_yp_zp_corner_pop);
 free(xm_yp_zm_corner_pop);
 free(xm_ym_zp_corner_pop);
 free(xm_ym_zm_corner_pop);

/* 12 edges */
free(xp_yp_edge_pop); 
free(xp_ym_edge_pop);
free(xm_yp_edge_pop);
free(xm_ym_edge_pop);

free(xp_zp_edge_pop);
free(xp_zm_edge_pop);
free(xm_zp_edge_pop);
free(xm_zm_edge_pop);

free(yp_zp_edge_pop);
free(yp_zm_edge_pop);
free(ym_zp_edge_pop);
free(ym_zm_edge_pop);
#endif


#endif

#ifdef LB_FLUID
 free(ruler_x_local);
 free(ruler_y_local);
 free(ruler_z_local);
 free(ruler_x);
 free(ruler_y);
 free(ruler_z);
 free(ruler_x_running);
 free(ruler_y_running);
 free(ruler_z_running);

#endif

#ifdef LB_TEMPERATURE
 free(g);
 free(rhs_g);
 free(old_rhs_g);
 free(old_old_rhs_g);
 free(t);
 #ifdef LB_TEMPERATURE_PAST
  free(old_t);
 #endif

 #ifdef LB_TEMPERATURE_FORCING
 free(t_source);
  #ifdef LB_TEMPERATURE_FORCING_PAST
   free(old_t_source);
  #endif
  #ifdef LB_TEMPERATURE_MELTING
  free(liquid_frac);  
  free(liquid_frac_old);  
  #endif
 #endif
#endif

#ifdef LB_SCALAR
 free(h); 
 free(rhs_h);  
 free(old_rhs_h);  
 free(old_old_rhs_h);  
 free(s);  

 #ifdef LB_SCALAR_FORCING
 free(s_source);
 #endif
#endif

}



