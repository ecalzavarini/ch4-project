#include "common_object.h"

#ifdef LAGRANGE


/* allocate particle containers */
void allocate_particles(){        

  int i;

  if((int)property.particle_number%nprocs ==0 ){

    /* ALL processor will take the same number of particles  */
    npart = property.particle_number/nprocs;


  }else{

   if(ROOT) fprintf(stderr,"Warning : total particle number is different from the process number!\n");

  for (i=0;i<nprocs;i++){

    /* ROOT processor will take just a little bit more particles */
    if(ROOT) npart = (int)property.particle_number - (nprocs-1)*(int)floor((int)property.particle_number/nprocs); else npart = (int)floor((int)property.particle_number /nprocs);  
					       
  }

  }/* end if on particle_number */

  fprintf(stderr,"me %d : I have %d particles\n",me, npart);

tracer  = (point_particle*) malloc(sizeof(point_particle)*npart);
if(tracer == NULL){ fprintf(stderr,"Not enough memory to allocate tracer\n"); exit(-1);}

tracer_here  = (point_particle*) malloc(sizeof(point_particle)*npart);
if(tracer_here == NULL){ fprintf(stderr,"Not enough memory to allocate tracer_here\n"); exit(-1);}

tracer_there  = (point_particle*) malloc(sizeof(point_particle)*npart);
if(tracer_there == NULL){ fprintf(stderr,"Not enough memory to allocate tracer_there\n"); exit(-1);}

all_tracer_there  = (point_particle*) malloc(sizeof(point_particle)*npart*nprocs);
if(all_tracer_there == NULL){ fprintf(stderr,"Not enough memory to allocate all_tracer_there\n"); exit(-1);}

}

/* initial conditions for particles */
void initial_conditions_particles(int restart){  

  int i,j,k,n;
  int *rcounts;
  int name_offset = 0;
  int ip,im,jp,jm,kp,km;
  my_double solid,theta,phi;
  char fnamein[128];
  FILE *fin;

  int type;
  my_double cycles, repetitions,step, *tau_drag, *beta_coeff, *aspect_ratio, *gyrotaxis_velocity, *rotational_diffusion, *swim_velocity;
  my_double *particle_radius, *particle_density;
  my_double *gravity_coeff;
  my_double *critical_shear_rate, *jump_time;
  vector vec;
  double *type_counter_local,*type_counter_all;

    rcounts = (int*)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);

    fprintf(stderr,"me : %d , name_offset %d\n", me , name_offset);


 /* restart from file */
  sprintf(fnamein,"part.h5");  
  fin = fopen("part.h5","r");  
  
 if(restart && fin != NULL){
  
     read_point_particle_h5();
     sendrecv_particles();
      if(ROOT) fprintf(stderr,"The %s file is present!\n Particles will be initialized from file.\n",fnamein);
      fclose(fin);

    }else{
    
      if(ROOT) fprintf(stderr,"Warning message -> %s file is missing!\n Particles will be initialized from memory.\n",fnamein);
     

/* Restart from memory */

      repetitions = 1.0;

#ifdef LAGRANGE_RADIUSandDENSITY
/* if the particle_radius and the particle_density is specified instead of tau_drag and beta_coeff */
      /* RADIUS */
 particle_radius = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
 cycles = (property.particle_types-property.fluid_tracers)/property.particle_radius_types;
 /* default: linear increment */
 if(property.particle_radius_types > 1 )  step = (property.particle_radius_max - property.particle_radius_min)/(property.particle_radius_types-1.0); else step = 0;
 #ifdef LAGRANGE_RADIUSandDENSITY_INCREMENT_LOG_RADIUS
  /* geometric increment */
  if(property.particle_radius_types > 1 )  step = pow(property.particle_radius_max/property.particle_radius_min, 1.0/(property.particle_radius_types-1.0)); else step = 0;
 #endif
 for (n=0; n<(int)property.fluid_tracers; n++) particle_radius[n] = 0.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions); k++){
   for (i=0; i<(int)property.particle_radius_types; i++){
     for (j=0; j<(int)repetitions; j++){
 /* default: linear increment */
   particle_radius[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.particle_radius_types ] = property.particle_radius_min + i*step;
 #ifdef LAGRANGE_RADIUSandDENSITY_INCREMENT_LOG_RADIUS
  /* geometric increment */
   particle_radius[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.particle_radius_types ] = property.particle_radius_min * pow(step,(double)i);
 #endif
     }
   }
 }
 repetitions *= property.particle_radius_types; 
 for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d particle_radius %g\n",i,particle_radius[i]);

      /* DENSITY */
 particle_density = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
 cycles = (property.particle_types-property.fluid_tracers)/property.particle_density_types;
 /* linear increment */
 if(property.particle_density_types > 1 )  step = (property.particle_density_max - property.particle_density_min)/(property.particle_density_types-1.0); else step = 0;
 #ifdef LAGRANGE_RADIUSandDENSITY_INCREMENT_LOG_DENSITY
 /* geometric increment */
 if(property.particle_density_types > 1 )  step = pow(property.particle_density_max/property.particle_density_min, 1.0/(property.particle_density_types-1.0)); else step = 0;
 #endif
 for (n=0; n<(int)property.fluid_tracers; n++) particle_density[n] = 1.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.particle_density_types; i++){
     for (j=0; j<(int)repetitions; j++){
 /* default: linear increment */
       particle_density[ (int)property.fluid_tracers +  j+i*(int)repetitions+k*(int)repetitions*(int)property.particle_density_types ] = property.particle_density_min + i*step;
 #ifdef LAGRANGE_RADIUSandDENSITY_INCREMENT_LOG_DENSITY
 /* geometric increment */
       particle_density[ (int)property.fluid_tracers +  j+i*(int)repetitions+k*(int)repetitions*(int)property.particle_density_types ] = property.particle_density_min * pow(step,(double)i);
 #endif
     }
   }
 }
 repetitions *= property.particle_density_types; 
 for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d particle_density %g\n",i,particle_density[i]);


 /* derive tau_drag from radius and density */
 tau_drag = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
 for(i=0;i<property.particle_types;i++){
   /* we use the definition \tau = (2/9)*(\rho_p / \rho_f)*r^2/\nu  , with assumption \rho_f = 1 */
   tau_drag[i] = (2./9.) * particle_density[i] * pow(particle_radius[i],2.0)/(property.nu);
     }

#ifdef LAGRANGE_GRADIENT
 #ifdef LAGRANGE_ADDEDMASS
 /* derive beta_coeff from density , assuming fluid density equal to 1*/
 beta_coeff = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
 for(i=0;i<property.particle_types;i++){
   beta_coeff[i] = 3.0/(1.0+2.0*particle_density[i]);
   /* adjust tau for added mass \tau = r^2/(3*\beta*\nu), with assumption \rho_f = 1 */
   tau_drag[i] = pow(particle_radius[i],2.0)/(3.0*beta_coeff[i]*property.nu);
     }
 for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d beta_coeff %g\n",i,beta_coeff[i]);
 #endif
#endif
 for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d tau_drag %g\n",i,tau_drag[i]);

#else /* else of LAGRANGE_RADIUSandDENSITY */

/* build particle property types arrays */
 tau_drag = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
 cycles = (property.particle_types-property.fluid_tracers)/property.tau_drag_types;
 if(property.tau_drag_types > 1 )  step = (property.tau_drag_max - property.tau_drag_min)/(property.tau_drag_types-1.0); else step = 0;
 for (n=0; n<(int)property.fluid_tracers; n++) tau_drag[n] = 0.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.tau_drag_types; i++){
     for (j=0; j<(int)repetitions; j++){
       tau_drag[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.tau_drag_types ] = property.tau_drag_min + i*step;
     }
   }
 }
 repetitions *= property.tau_drag_types; 
 for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d tau_drag %g\n",i,tau_drag[i]);
#ifdef LAGRANGE_GRADIENT
 #ifdef LAGRANGE_ADDEDMASS
 beta_coeff = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
 cycles = (property.particle_types-property.fluid_tracers)/property.beta_coeff_types;
 if(property.beta_coeff_types > 1 )  step = (property.beta_coeff_max - property.beta_coeff_min)/(property.beta_coeff_types-1.0); else step = 0;
 for (n=0; n<(int)property.fluid_tracers; n++) beta_coeff[n] = 1.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.beta_coeff_types; i++){
     for (j=0; j<(int)repetitions; j++){
       beta_coeff[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.beta_coeff_types ] = property.beta_coeff_min + i*step;
     }
   }
 }
 repetitions *= property.beta_coeff_types; 
 for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d beta_coeff %g\n",i,beta_coeff[i]);
 #endif
#endif
#endif /* end of LAGRANGE_RADIUSandDENSITY */

#ifdef LAGRANGE_GRAVITY
 #ifdef LAGRANGE_GRAVITY_VARIABLE
 gravity_coeff = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
 cycles = (property.particle_types-property.fluid_tracers)/property.gravity_coeff_types;
 /* default: linear increment */
 if(property.gravity_coeff_types > 1 )  step = (property.gravity_coeff_max - property.gravity_coeff_min)/(property.gravity_coeff_types-1.0); else step = 0;
  #ifdef LAGRANGE_GRAVITY_VARIABLE_INCREMENT_LOG
  /* geometric increment */
  if(property.gravity_coeff_types > 1 )  step = pow(property.gravity_coeff_max/property.gravity_coeff_min, 1.0/(property.gravity_coeff_types-1.0)); else step = 0;
  #endif
 for (n=0; n<(int)property.fluid_tracers; n++) gravity_coeff[n] = 0.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.gravity_coeff_types; i++){
     for (j=0; j<(int)repetitions; j++){
 /* default: linear increment */
       gravity_coeff[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.gravity_coeff_types ] = property.gravity_coeff_min + i*step;
  #ifdef LAGRANGE_GRAVITY_VARIABLE_INCREMENT_LOG
 /* geometric increment */
       gravity_coeff[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.gravity_coeff_types ] = property.gravity_coeff_min * pow(step,(double)i);
  #endif
     }
   }
 }
 repetitions *= property.gravity_coeff_types; 
 for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d gravity_coeff %g\n",i,gravity_coeff[i]);
 #endif
#endif

#ifdef LAGRANGE_GRADIENT
 #ifdef LAGRANGE_ORIENTATION
  #ifdef LAGRANGE_ORIENTATION_JEFFREY
 aspect_ratio = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
 cycles = (property.particle_types-property.fluid_tracers)/property.aspect_ratio_types;
 /* deafult: linear increment */
 if(property.aspect_ratio_types > 1 )  step = (property.aspect_ratio_max - property.aspect_ratio_min)/(property.aspect_ratio_types-1.0); else step = 0;
   #ifdef LAGRANGE_ORIENTATION_JEFFREY_INCREMENT_LOG
 /* geometric increment */
 if(property.aspect_ratio_types > 1 )  step = pow(property.aspect_ratio_max/property.aspect_ratio_min, 1.0/(property.aspect_ratio_types-1.0)); else step = 0;
   #endif
 for (n=0; n<(int)property.fluid_tracers; n++) aspect_ratio[n] = 1.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.aspect_ratio_types; i++){
     for (j=0; j<(int)repetitions; j++){
 /* default: linear increment */
       aspect_ratio[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.aspect_ratio_types ] = property.aspect_ratio_min + i*step;
   #ifdef LAGRANGE_ORIENTATION_JEFFREY_INCREMENT_LOG
 /* geometric increment */
       aspect_ratio[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.aspect_ratio_types ] = property.aspect_ratio_min * pow(step,(double)i);
   #endif
     }
   }
 }
 repetitions *= property.aspect_ratio_types; 
 for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d aspect_ratio %g\n",i,aspect_ratio[i]);

   #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
  gyrotaxis_velocity = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
  cycles = (property.particle_types-property.fluid_tracers)/property.gyrotaxis_velocity_types;
 /* linear increment */
 //if(property.gyrotaxis_velocity_types > 1 )  step = (property.gyrotaxis_velocity_max - property.gyrotaxis_velocity_min)/(property.gyrotaxis_velocity_types-1.0); else step = 0;
 /* geometric increment */
 if(property.gyrotaxis_velocity_types > 1 )  step = pow(property.gyrotaxis_velocity_max/property.gyrotaxis_velocity_min, 1.0/(property.gyrotaxis_velocity_types-1.0)); else step = 0;
 for (n=0; n<(int)property.fluid_tracers; n++) gyrotaxis_velocity[n] = 0.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.gyrotaxis_velocity_types; i++){
     for (j=0; j<(int)repetitions; j++){
 /* linear increment */
 //      gyrotaxis_velocity[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.gyrotaxis_velocity_types ] = property.gyrotaxis_velocity_min + i*step;
/* geometric increment */
       gyrotaxis_velocity[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.gyrotaxis_velocity_types ] = property.gyrotaxis_velocity_min * pow(step,(double)i);
     }
   }
 }
 repetitions *= property.gyrotaxis_velocity_types; 
  for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d gyrotaxis_velocity %g\n",i,gyrotaxis_velocity[i]);
   #endif
  #endif /* LAGRANGE_ORIENTATION_JEFFREY */

  #ifdef LAGRANGE_ORIENTATION_DIFFUSION
  rotational_diffusion = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
  cycles = (property.particle_types-property.fluid_tracers)/property.rotational_diffusion_types;
  if(property.rotational_diffusion_types > 1 )  step = (property.rotational_diffusion_max - property.rotational_diffusion_min)/(property.rotational_diffusion_types-1.0); else step = 0;
 for (n=0; n<(int)property.fluid_tracers; n++) rotational_diffusion[n]  = 0.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.rotational_diffusion_types; i++){
     for (j=0; j<(int)repetitions; j++){
       rotational_diffusion[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.rotational_diffusion_types ] = property.rotational_diffusion_min + i*step;
     }
   }
 }
 repetitions *= property.rotational_diffusion_types; 
  for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d rotational_diffusion %g\n",i,rotational_diffusion[i]);
  #endif /* LAGRANGE_ORIENTATION_DIFFUSION */

  #ifdef LAGRANGE_ORIENTATION_ACTIVE
  swim_velocity = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
  cycles = (property.particle_types-property.fluid_tracers)/property.swim_velocity_types;
  if(property.swim_velocity_types > 1 )  step = (property.swim_velocity_max - property.swim_velocity_min)/(property.swim_velocity_types-1.0); else step = 0;
 for (n=0; n<(int)property.fluid_tracers; n++) swim_velocity[n]  = 0.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.swim_velocity_types; i++){
     for (j=0; j<(int)repetitions; j++){
       swim_velocity[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.swim_velocity_types ] = property.swim_velocity_min + i*step;
     }
   }
 }
 repetitions *= property.swim_velocity_types; 
  for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d swim_velocity %g\n",i,swim_velocity[i]);

   #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
  /* critical shear rate */
  critical_shear_rate = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
  cycles = (property.particle_types-property.fluid_tracers)/property.critical_shear_rate_types;
  if(property.critical_shear_rate_types > 1 )  step = (property.critical_shear_rate_max - property.critical_shear_rate_min)/(property.critical_shear_rate_types-1.0); else step = 0;
 for (n=0; n<(int)property.fluid_tracers; n++) critical_shear_rate[n] = 1.e+10; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.critical_shear_rate_types; i++){
     for (j=0; j<(int)repetitions; j++){
       critical_shear_rate[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.critical_shear_rate_types ] = property.critical_shear_rate_min + i*step;
     }
   }
 }
 repetitions *= property.critical_shear_rate_types; 
  for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d critical_shear_rate %g\n",i,critical_shear_rate[i]);

  /* jump_time */
  jump_time = (my_double*) malloc(sizeof(my_double)*property.particle_types); 
  cycles = (property.particle_types-property.fluid_tracers)/property.jump_time_types;
  if(property.jump_time_types > 1 )  step = (property.jump_time_max - property.jump_time_min)/(property.jump_time_types-1.0); else step = 0;
 for (n=0; n<(int)property.fluid_tracers; n++) jump_time[n] = 0.0; /* this is for the tracer */
 for (k=0; k<(int)(cycles/repetitions)  ; k++){
   for (i=0; i<(int)property.jump_time_types; i++){
     for (j=0; j<(int)repetitions; j++){
       jump_time[ (int)property.fluid_tracers + j+i*(int)repetitions+k*(int)repetitions*(int)property.jump_time_types ] = property.jump_time_min + i*step;
     }
   }
 }
 repetitions *= property.jump_time_types; 
  for(i=0;i<property.particle_types;i++) if(ROOT) fprintf(stderr,"type %d jump_time %g\n",i,jump_time[i]);
   #endif /* LAGRANGE_ORIENTATION_ACTIVE_JUMP */
  #endif /* LAGRANGE_ORIENTATION_ACTIVE */
 #endif /* LAGRANGE_ORIENTATION */
#endif /* LAGRANGE_GRADIENT */

/* Here we ASSIGN IDENTIFICATION NUMBERS (NAMES) to particles, and
   at the same time count the number of particles per family */
  
  /* malloc counters */
  type_counter_local = (my_double*) malloc(sizeof(my_double)*(int)property.particle_types);
  type_counter_all = (my_double*) malloc(sizeof(my_double)*(int)property.particle_types);

  /* set counters to zero */
 for (type=0;type<(int)property.particle_types;type++)  type_counter_local[type] = type_counter_all[type] = 0.0;

 /* loop on particles  particles */
 for (i=0;i<npart;i++) {

   /* name */
   (tracer+i)->name = i+name_offset;

   /* family */
   type = ((int)(tracer+i)->name)%(int)property.particle_types;
   type_counter_local[type] += 1.0;
 }
 /* Sum */
 MPI_Allreduce(type_counter_local, type_counter_all, (int)property.particle_types, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
 

/* write on file particle properties by family */
if(ROOT){

  fin = fopen("particle_properties.dat","w");
  for(i=0;i<property.particle_types;i++){

fprintf(fin,"type %d counter_all %d ",i, (int)type_counter_all[i]);

#ifdef LAGRANGE_RADIUSandDENSITY
fprintf(fin,"particle_radius %e ",particle_radius[i]);
fprintf(fin,"particle_density %e ",particle_density[i]);
fprintf(fin,"tau_drag %e ",tau_drag[i]);
#else
fprintf(fin,"tau_drag %e ",i,tau_drag[i]);
#endif

#ifdef LAGRANGE_GRADIENT
 #ifdef LAGRANGE_ADDEDMASS
  fprintf(fin,"beta_coeff %e ",beta_coeff[i]);
 #endif /* LAGRANGE_ADDEDMASS */
#endif

#ifdef LAGRANGE_GRAVITY
 #ifdef LAGRANGE_GRAVITY_VARIABLE
 fprintf(fin,"gravity_coeff %e ",gravity_coeff[i]);
 #endif
#endif

#ifdef LAGRANGE_GRADIENT
 #ifdef LAGRANGE_ORIENTATION
  #ifdef LAGRANGE_ORIENTATION_JEFFREY
 fprintf(fin,"aspect_ratio %e ",aspect_ratio[i]);
   #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
   fprintf(fin,"gyrotaxis_velocity %e ",gyrotaxis_velocity[i]);
   #endif
  #endif
  #ifdef LAGRANGE_ORIENTATION_DIFFUSION
   fprintf(fin,"rotational_diffusion %e ",rotational_diffusion[i]);
  #endif
  #ifdef LAGRANGE_ORIENTATION_ACTIVE
   fprintf(fin,"swim_velocity %e ",swim_velocity[i]);
   #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
    fprintf(fin,"jump_time %e ",jump_time[i]);
    fprintf(fin,"critical_shear_rate %e ",critical_shear_rate[i]);
   #endif
  #endif
 #endif /* LAGRANGE_ORIENTATION */
#endif /* LAGRANGE_GRADIENT */
 fprintf(fin,"\n");
 } 
  fclose(fin);
 }/* end if ROOT */


/* assign param values to particles */
for (i=0;i<npart;i++) {

/* name */
//(tracer+i)->name = i+name_offset;

type = ((int)(tracer+i)->name)%(int)property.particle_types;

/* viscous drag */
(tracer+i)->tau_drag = tau_drag[type];
//(tracer+i)->tau_drag = 0.0;
#ifdef LAGRANGE_GRAVITY
 #ifdef LAGRANGE_GRAVITY_VARIABLE
/* gravity coefficient: modulates gravity per particle type */
(tracer+i)->gravity_coeff = gravity_coeff[type];
 #endif
#endif
#ifdef LAGRANGE_GRADIENT
 #ifdef LAGRANGE_ADDEDMASS
 /* added mass */
 (tracer+i)->beta_coeff = beta_coeff[type];
//(tracer+i)->beta_coeff = 1.0;
 #endif
 #ifdef LAGRANGE_ORIENTATION
  #ifdef LAGRANGE_ORIENTATION_JEFFREY
 /* particle aspect ratio (assuming axi-symmetry) */
 (tracer+i)->aspect_ratio = aspect_ratio[type];
 //(tracer+i)->aspect_ratio = 100.0;
   #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
  /* gyrotaxis rotational parameter: the velocity  v_0 parameter,  as in F.De Lillo, M. Cencini et al., PRL 112, 044502 (2014) */
  //(tracer+i)->gyrotaxis_velocity = 1.0;
  (tracer+i)->gyrotaxis_velocity = gyrotaxis_velocity[type];
   #endif
  #endif
  #ifdef LAGRANGE_ORIENTATION_DIFFUSION
  /* rotational diffusion , units ? [rad^2 /time] */
  //(tracer+i)->rotational_diffusion = 0.1;
  (tracer+i)->rotational_diffusion = rotational_diffusion[type];
  #endif 
  #ifdef LAGRANGE_ORIENTATION_ACTIVE 
   //(tracer+i)->swim_velocity = 0.01;
   (tracer+i)->swim_velocity = swim_velocity[type];
   #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
   //(tracer+i)->jump_time = 0.01;
   (tracer+i)->jump_time = jump_time[type];
   //(tracer+i)->critical_shear_rate = 0.01;
   (tracer+i)->critical_shear_rate = critical_shear_rate[type];
   #endif  
  #endif
 #endif
#endif

/* position: randomly distributed particles */
(tracer+i)->x = LNX_START + myrand()*LNX;
(tracer+i)->y = LNY_START + myrand()*LNY;
(tracer+i)->z = LNZ_START + myrand()*LNZ;
/* if in 2d */
#ifdef GRID_POP_D2Q9
(tracer+i)->z = LNZ_START + 0.5*LNZ;
#endif

#ifdef LB_FLUID_FORCING_LANDSCAPE
/* This part is to not to put particle in static (solid)  LANDSCAPE */
 solid = 1.0;
 while(solid==1.0){

(tracer+i)->x = LNX_START + myrand()*LNX;
(tracer+i)->y = LNY_START + myrand()*LNY;
(tracer+i)->z = LNZ_START + myrand()*LNZ;
/* if in 2d */
#ifdef GRID_POP_D2Q9
(tracer+i)->z = LNZ_START + 0.5*LNZ;
#endif

for (n=0; n<LNX+TWO_BRD-1; n++) if(center_V[IDX(n, BRD, BRD)].x <= (tracer+i)->x && (tracer+i)->x < center_V[IDX(n+1,BRD, BRD)].x) im = n; 
ip =  im + 1;
for (j=0; j<LNY+TWO_BRD-1; j++) if(center_V[IDX(BRD, j, BRD)].y <= (tracer+i)->y && (tracer+i)->y < center_V[IDX(BRD, j+1, BRD)].y) jm = j;
jp =  jm + 1;
for (k=0; k<LNZ+TWO_BRD-1; k++) if(center_V[IDX(BRD, BRD, k)].z <= (tracer+i)->z && (tracer+i)->z < center_V[IDX(BRD, BRD, k+1)].z) km = k;
kp =  km + 1;


  solid = landscape[IDX(im, jm, km)] * landscape[IDX(ip, jm, km)] * landscape[IDX(im, jp, km)]
        * landscape[IDX(im, jm, kp)] * landscape[IDX(ip, jp, km)] * landscape[IDX(im, jp, kp)]
        * landscape[IDX(ip, jm, kp)] * landscape[IDX(ip, jp, kp)];
 }
#endif

/* velocity: null speed */
(tracer+i)->vx = 0.0;
(tracer+i)->vy = 0.0;
(tracer+i)->vz = 0.0;

#ifdef LAGRANGE_ORIENTATION
/*
(tracer+i)->px = 0.0;
(tracer+i)->py = 1.0;
(tracer+i)->pz = 0.0;
*/
/* uniform random oriented vector */
/*
theta = acos(2.0*myrand()-1.0);
phi = two_pi*myrand();
(tracer+i)->px = sin(theta)*sin(phi);
(tracer+i)->py = sin(theta)*cos(phi);
(tracer+i)->pz = cos(theta);
*/
 vec = random_vector();
#ifdef GRID_POP_D2Q9
 /* generate a random vector in the x,y plane */
 vec = random_vector_2d();
#endif
 (tracer+i)->px = vec.x;
 (tracer+i)->py = vec.y;
 (tracer+i)->pz = vec.z;

(tracer+i)->dt_px = 0.0;
(tracer+i)->dt_py = 0.0;
(tracer+i)->dt_pz = 0.0;
 #ifdef LAGRANGE_ORIENTATION_ACTIVE
  #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
/* time from jump is initialized in a arbitrary range between 0 and 10 jump_times */
    (tracer+i)->time_from_jump = 10.*myrand()*(tracer+i)->jump_time;
  #endif
 #endif
#endif
 }/* end of  loop on particles */

/*
 free(tau_drag);
 free(beta_coeff);
 free(aspect_ratio);
 free(gyrotaxis_velocity);
 free(rotational_diffusion);
 free(swim_velocity);
*/

 }/* end of if /else on restart */


}/* end of function */



/* Interpolation for the moment only with regular grid */
void interpolate_vector_at_particles(vector *f,char which_vector){

 point_particle part;
 vector v;

 int ipart,im,jm,km,ip,jp,kp;
 double dxm,dxp,dym,dyp,dzm,dzp;
 double vol_ip_jp_kp,vol_im_jp_kp , vol_ip_jm_kp , vol_ip_jp_km , vol_im_jm_kp , vol_ip_jm_km , vol_im_jp_km , vol_im_jm_km; 

 int i,j,k;

#ifdef LAGRANGE_GRADIENT
 tensor grad, grad_im_jm_km , grad_ip_jm_km , grad_im_jp_km , grad_im_jm_kp , grad_ip_jp_km , grad_im_jp_kp , grad_ip_jm_kp , grad_ip_jp_kp; 
#endif


 for (ipart=0;ipart<npart;ipart++) {

   // fprintf(stderr,"vel %g %g\n", f[IDX(BRD-1,BRD,BRD)].x , f[IDX(BRD,BRD,BRD)].x );

  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

  //fprintf(stderr,"\n part %g %g %g\n",part.x,part.y,part.z);

for (i=0; i<LNX+TWO_BRD-1; i++) if(center_V[IDX(i, BRD, BRD)].x <= part.x && part.x < center_V[IDX(i+1,BRD, BRD)].x) im = i; 
ip =  im + 1;
for (j=0; j<LNY+TWO_BRD-1; j++) if(center_V[IDX(BRD, j, BRD)].y <= part.y && part.y < center_V[IDX(BRD, j+1, BRD)].y) jm = j;
jp =  jm + 1;
for (k=0; k<LNZ+TWO_BRD-1; k++) if(center_V[IDX(BRD, BRD, k)].z <= part.z && part.z < center_V[IDX(BRD, BRD, k+1)].z) km = k;
kp =  km + 1;

//fprintf(stderr,"index %d %d %d %d %d %d\n", im,ip,jm,jp,km,kp);

//for (j=0;j<10;j++) fprintf(stderr,"%d center_V %e\n",j,center_V[IDX(im, j, km)].y);

  dxm = part.x - center_V[IDX(im, BRD, BRD)].x;
  dxp = center_V[IDX(ip, BRD, BRD)].x - part.x;
  dym = part.y - center_V[IDX(BRD, jm, BRD)].y;
  dyp = center_V[IDX(BRD, jp, BRD)].y - part.y;
  dzm = part.z - center_V[IDX(BRD, BRD, km)].z;
  dzp = center_V[IDX(BRD, BRD, kp)].z - part.z;

  //fprintf(stderr,"distance %g %g %g %g %g %g\n", dxm,dxp,dym,dyp,dzm,dzp);

  vol_ip_jp_kp = dxp*dyp*dzp;
  vol_im_jp_kp = dxm*dyp*dzp;
  vol_ip_jm_kp = dxp*dym*dzp;
  vol_ip_jp_km = dxp*dyp*dzm;
  vol_im_jm_kp = dxm*dym*dzp;
  vol_ip_jm_km = dxp*dym*dzm;
  vol_im_jp_km = dxm*dyp*dzm;
  vol_im_jm_km = dxm*dym*dzm;


  v.x = f[IDX(im, jm, km)].x * vol_ip_jp_kp + 
        f[IDX(ip, jm, km)].x * vol_im_jp_kp +
        f[IDX(im, jp, km)].x * vol_ip_jm_kp +
        f[IDX(im, jm, kp)].x * vol_ip_jp_km +
        f[IDX(ip, jp, km)].x * vol_im_jm_kp +
        f[IDX(im, jp, kp)].x * vol_ip_jm_km +
        f[IDX(ip, jm, kp)].x * vol_im_jp_km +
        f[IDX(ip, jp, kp)].x * vol_im_jm_km ;

  v.y = f[IDX(im, jm, km)].y * vol_ip_jp_kp + 
        f[IDX(ip, jm, km)].y * vol_im_jp_kp +
        f[IDX(im, jp, km)].y * vol_ip_jm_kp +
        f[IDX(im, jm, kp)].y * vol_ip_jp_km +
        f[IDX(ip, jp, km)].y * vol_im_jm_kp +
        f[IDX(im, jp, kp)].y * vol_ip_jm_km +
        f[IDX(ip, jm, kp)].y * vol_im_jp_km +
        f[IDX(ip, jp, kp)].y * vol_im_jm_km ;

  v.z = f[IDX(im, jm, km)].z * vol_ip_jp_kp + 
        f[IDX(ip, jm, km)].z * vol_im_jp_kp +
        f[IDX(im, jp, km)].z * vol_ip_jm_kp +
        f[IDX(im, jm, kp)].z * vol_ip_jp_km +
        f[IDX(ip, jp, km)].z * vol_im_jm_kp +
        f[IDX(im, jp, kp)].z * vol_ip_jm_km +
        f[IDX(ip, jm, kp)].z * vol_im_jp_km +
        f[IDX(ip, jp, kp)].z * vol_im_jm_km ;

  /* if it is the velocity */
  if(which_vector == 'u'){
    (tracer+ipart)->ux = v.x;
    (tracer+ipart)->uy = v.y;
    (tracer+ipart)->uz = v.z;
    /* for two dimensions */
    //    if(NZ==1) (tracer+ipart)->uz = 0.0;
  }

#ifdef LAGRANGE_GRADIENT
  /* here we interpolate also the gradient of the same field */
  grad_im_jm_km = gradient_vector(f,im,jm,km);
  grad_ip_jm_km = gradient_vector(f,ip,jm,km);
  grad_im_jp_km = gradient_vector(f,im,jp,km);
  grad_im_jm_kp = gradient_vector(f,im,jm,kp);
  grad_ip_jp_km = gradient_vector(f,ip,jp,km);
  grad_im_jp_kp = gradient_vector(f,im,jp,kp);
  grad_ip_jm_kp = gradient_vector(f,ip,jm,kp);
  grad_ip_jp_kp = gradient_vector(f,ip,jp,kp);

 grad.xx =  grad_im_jm_km.xx * vol_ip_jp_kp + 
            grad_ip_jm_km.xx * vol_im_jp_kp +
            grad_im_jp_km.xx * vol_ip_jm_kp +
            grad_im_jm_kp.xx * vol_ip_jp_km +
            grad_ip_jp_km.xx * vol_im_jm_kp +
            grad_im_jp_kp.xx * vol_ip_jm_km +
            grad_ip_jm_kp.xx * vol_im_jp_km +
            grad_ip_jp_kp.xx * vol_im_jm_km ;

 grad.xy =  grad_im_jm_km.xy * vol_ip_jp_kp + 
            grad_ip_jm_km.xy * vol_im_jp_kp +
            grad_im_jp_km.xy * vol_ip_jm_kp +
            grad_im_jm_kp.xy * vol_ip_jp_km +
            grad_ip_jp_km.xy * vol_im_jm_kp +
            grad_im_jp_kp.xy * vol_ip_jm_km +
            grad_ip_jm_kp.xy * vol_im_jp_km +
            grad_ip_jp_kp.xy * vol_im_jm_km ;

 grad.xz =  grad_im_jm_km.xz * vol_ip_jp_kp + 
            grad_ip_jm_km.xz * vol_im_jp_kp +
            grad_im_jp_km.xz * vol_ip_jm_kp +
            grad_im_jm_kp.xz * vol_ip_jp_km +
            grad_ip_jp_km.xz * vol_im_jm_kp +
            grad_im_jp_kp.xz * vol_ip_jm_km +
            grad_ip_jm_kp.xz * vol_im_jp_km +
            grad_ip_jp_kp.xz * vol_im_jm_km ;

 grad.yx =  grad_im_jm_km.yx * vol_ip_jp_kp + 
            grad_ip_jm_km.yx * vol_im_jp_kp +
            grad_im_jp_km.yx * vol_ip_jm_kp +
            grad_im_jm_kp.yx * vol_ip_jp_km +
            grad_ip_jp_km.yx * vol_im_jm_kp +
            grad_im_jp_kp.yx * vol_ip_jm_km +
            grad_ip_jm_kp.yx * vol_im_jp_km +
            grad_ip_jp_kp.yx * vol_im_jm_km ;

 grad.yy =  grad_im_jm_km.yy * vol_ip_jp_kp + 
            grad_ip_jm_km.yy * vol_im_jp_kp +
            grad_im_jp_km.yy * vol_ip_jm_kp +
            grad_im_jm_kp.yy * vol_ip_jp_km +
            grad_ip_jp_km.yy * vol_im_jm_kp +
            grad_im_jp_kp.yy * vol_ip_jm_km +
            grad_ip_jm_kp.yy * vol_im_jp_km +
            grad_ip_jp_kp.yy * vol_im_jm_km ;

 grad.yz =  grad_im_jm_km.yz * vol_ip_jp_kp + 
            grad_ip_jm_km.yz * vol_im_jp_kp +
            grad_im_jp_km.yz * vol_ip_jm_kp +
            grad_im_jm_kp.yz * vol_ip_jp_km +
            grad_ip_jp_km.yz * vol_im_jm_kp +
            grad_im_jp_kp.yz * vol_ip_jm_km +
            grad_ip_jm_kp.yz * vol_im_jp_km +
            grad_ip_jp_kp.yz * vol_im_jm_km ;

 grad.zx =  grad_im_jm_km.zx * vol_ip_jp_kp + 
            grad_ip_jm_km.zx * vol_im_jp_kp +
            grad_im_jp_km.zx * vol_ip_jm_kp +
            grad_im_jm_kp.zx * vol_ip_jp_km +
            grad_ip_jp_km.zx * vol_im_jm_kp +
            grad_im_jp_kp.zx * vol_ip_jm_km +
            grad_ip_jm_kp.zx * vol_im_jp_km +
            grad_ip_jp_kp.zx * vol_im_jm_km ;

 grad.zy =  grad_im_jm_km.zy * vol_ip_jp_kp + 
            grad_ip_jm_km.zy * vol_im_jp_kp +
            grad_im_jp_km.zy * vol_ip_jm_kp +
            grad_im_jm_kp.zy * vol_ip_jp_km +
            grad_ip_jp_km.zy * vol_im_jm_kp +
            grad_im_jp_kp.zy * vol_ip_jm_km +
            grad_ip_jm_kp.zy * vol_im_jp_km +
            grad_ip_jp_kp.zy * vol_im_jm_km ;

grad.zz =   grad_im_jm_km.zz * vol_ip_jp_kp + 
            grad_ip_jm_km.zz * vol_im_jp_kp +
            grad_im_jp_km.zz * vol_ip_jm_kp +
            grad_im_jm_kp.zz * vol_ip_jp_km +
            grad_ip_jp_km.zz * vol_im_jm_kp +
            grad_im_jp_kp.zz * vol_ip_jm_km +
            grad_ip_jm_kp.zz * vol_im_jp_km +
            grad_ip_jp_kp.zz * vol_im_jm_km ;


  /* if it is the velocity */
  if(which_vector == 'u'){
    (tracer+ipart)->dx_ux = grad.xx;
    (tracer+ipart)->dy_ux = grad.xy;
    (tracer+ipart)->dz_ux = grad.xz;
    (tracer+ipart)->dx_uy = grad.yx;
    (tracer+ipart)->dy_uy = grad.yy;
    (tracer+ipart)->dz_uy = grad.yz;
    (tracer+ipart)->dx_uz = grad.zx;
    (tracer+ipart)->dy_uz = grad.zy;
    (tracer+ipart)->dz_uz = grad.zz;

  }
#endif

 }/* end of for on ipart */

}


/* Interpolation for the moment only with regular grid */
void interpolate_scalar_at_particles(my_double *f,char which_scalar){

 point_particle part;
 my_double s;

 int ipart,im,jm,km,ip,jp,kp;
 double dxm,dxp,dym,dyp,dzm,dzp;
 double vol_ip_jp_kp,vol_im_jp_kp , vol_ip_jm_kp , vol_ip_jp_km , vol_im_jm_kp , vol_ip_jm_km , vol_im_jp_km , vol_im_jm_km; 

int i,j,k;

#ifdef LAGRANGE_GRADIENT
 vector grad, grad_im_jm_km , grad_ip_jm_km , grad_im_jp_km , grad_im_jm_kp , grad_ip_jp_km , grad_im_jp_kp , grad_ip_jm_kp , grad_ip_jp_kp; 
#endif

 for (ipart=0;ipart<npart;ipart++) {

   // fprintf(stderr,"vel %g %g\n", f[IDX(BRD-1,BRD,BRD)].x , f[IDX(BRD,BRD,BRD)].x );

  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

  //fprintf(stderr,"\n part %g %g %g\n",part.x,part.y,part.z);

for (i=0; i<LNX+TWO_BRD-1; i++) if(center_V[IDX(i, BRD, BRD)].x <= part.x && part.x < center_V[IDX(i+1,BRD, BRD)].x) im = i; 
ip =  im + 1;
for (j=0; j<LNY+TWO_BRD-1; j++) if(center_V[IDX(BRD, j, BRD)].y <= part.y && part.y < center_V[IDX(BRD, j+1, BRD)].y) jm = j;
jp =  jm + 1;
for (k=0; k<LNZ+TWO_BRD-1; k++) if(center_V[IDX(BRD, BRD, k)].z <= part.z && part.z < center_V[IDX(BRD, BRD, k+1)].z) km = k;
kp =  km + 1;

//fprintf(stderr,"index %d %d %d %d %d %d\n", im,ip,jm,jp,km,kp);

//for (j=0;j<10;j++) fprintf(stderr,"%d center_V %e\n",j,center_V[IDX(im, j, km)].y);

  dxm = part.x - center_V[IDX(im, BRD, BRD)].x;
  dxp = center_V[IDX(ip, BRD, BRD)].x - part.x;
  dym = part.y - center_V[IDX(BRD, jm, BRD)].y;
  dyp = center_V[IDX(BRD, jp, BRD)].y - part.y;
  dzm = part.z - center_V[IDX(BRD, BRD, km)].z;
  dzp = center_V[IDX(BRD, BRD, kp)].z - part.z;

  //fprintf(stderr,"distance %g %g %g %g %g %g\n", dxm,dxp,dym,dyp,dzm,dzp);

  vol_ip_jp_kp = dxp*dyp*dzp;
  vol_im_jp_kp = dxm*dyp*dzp;
  vol_ip_jm_kp = dxp*dym*dzp;
  vol_ip_jp_km = dxp*dyp*dzm;
  vol_im_jm_kp = dxm*dym*dzp;
  vol_ip_jm_km = dxp*dym*dzm;
  vol_im_jp_km = dxm*dyp*dzm;
  vol_im_jm_km = dxm*dym*dzm;


   s  = f[IDX(im, jm, km)] * vol_ip_jp_kp + 
        f[IDX(ip, jm, km)] * vol_im_jp_kp +
        f[IDX(im, jp, km)] * vol_ip_jm_kp +
        f[IDX(im, jm, kp)] * vol_ip_jp_km +
        f[IDX(ip, jp, km)] * vol_im_jm_kp +
        f[IDX(im, jp, kp)] * vol_ip_jm_km +
        f[IDX(ip, jm, kp)] * vol_im_jp_km +
        f[IDX(ip, jp, kp)] * vol_im_jm_km ;

#ifdef LB_TEMPERATURE
   /* if it is temperature */
if(which_scalar == 't')  (tracer+ipart)->t = s;
#endif

#ifdef LB_SCALAR
  /* if it is a scalar */
if(which_scalar == 's')  (tracer+ipart)->s = s;
#endif


#ifdef LAGRANGE_GRADIENT
  /* here we interpolate also the gradient of the same field */
  grad_im_jm_km = gradient_scalar(f,im,jm,km);
  grad_ip_jm_km = gradient_scalar(f,ip,jm,km);
  grad_im_jp_km = gradient_scalar(f,im,jp,km);
  grad_im_jm_kp = gradient_scalar(f,im,jm,kp);
  grad_ip_jp_km = gradient_scalar(f,ip,jp,km);
  grad_im_jp_kp = gradient_scalar(f,im,jp,kp);
  grad_ip_jm_kp = gradient_scalar(f,ip,jm,kp);
  grad_ip_jp_kp = gradient_scalar(f,ip,jp,kp);

 grad.x =  grad_im_jm_km.x * vol_ip_jp_kp + 
             grad_ip_jm_km.x * vol_im_jp_kp +
             grad_im_jp_km.x * vol_ip_jm_kp +
             grad_im_jm_kp.x * vol_ip_jp_km +
             grad_ip_jp_km.x * vol_im_jm_kp +
             grad_im_jp_kp.x * vol_ip_jm_km +
             grad_ip_jm_kp.x * vol_im_jp_km +
             grad_ip_jp_kp.x * vol_im_jm_km ;

 grad.y =  grad_im_jm_km.y * vol_ip_jp_kp + 
             grad_ip_jm_km.y * vol_im_jp_kp +
             grad_im_jp_km.y * vol_ip_jm_kp +
             grad_im_jm_kp.y * vol_ip_jp_km +
             grad_ip_jp_km.y * vol_im_jm_kp +
             grad_im_jp_kp.y * vol_ip_jm_km +
             grad_ip_jm_kp.y * vol_im_jp_km +
             grad_ip_jp_kp.y * vol_im_jm_km ;

 grad.z =  grad_im_jm_km.z * vol_ip_jp_kp + 
             grad_ip_jm_km.z * vol_im_jp_kp +
             grad_im_jp_km.z * vol_ip_jm_kp +
             grad_im_jm_kp.z * vol_ip_jp_km +
             grad_ip_jp_km.z * vol_im_jm_kp +
             grad_im_jp_kp.z * vol_ip_jm_km +
             grad_ip_jm_kp.z * vol_im_jp_km +
             grad_ip_jp_kp.z * vol_im_jm_km ;

#ifdef LB_TEMPERATURE
  /* if it is the temperature */
  if(which_scalar == 't'){
    (tracer+ipart)->dx_t = grad.x;
    (tracer+ipart)->dy_t = grad.y;
    (tracer+ipart)->dz_t = grad.z;
  }
#endif

#ifdef LB_SCALAR
  /* if it is the scalar */
  if(which_scalar == 's'){
    (tracer+ipart)->dx_s = grad.x;
    (tracer+ipart)->dy_s = grad.y;
    (tracer+ipart)->dz_s = grad.z;
  }
#endif

#endif /* endif on lagrange_gradient */

 }/* end of loop on ipart */

}/* end of interp scalar*/


#define H5FILE_NAME_PARTICLE "particle.h5"

/* general output function for particles */
void output_particles(){
  int i,j;
  int np = (int)property.particle_number;
  FILE *fout;

    int *rcounts;
    int name_offset = 0;

    hid_t       file_id, dataset_id, dataspace_id , group ;  /* identifiers */
    hid_t	plist_id;                 /* property list identifier */
    hid_t hdf5_type;
    hid_t xfer_plist, ret, property_id;
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dims[1], offset[1], count[1];
    herr_t      hdf5_status;
    herr_t status;
    int size;    
    int RANK = 1;

    my_double *aux;

    char NEW_H5FILE_NAME[128];
    char XMF_FILE_NAME[128];

#ifdef LAGRANGE_OUTPUT_DEBUG
    if(ROOT){
     for (i=0;i<npart;i++) {
      fprintf(stdout,"%g %e %e %e %e %e %e\n",time_now, (tracer+i)->x,(tracer+i)->y,(tracer+i)->z,(tracer+i)->vx,(tracer+i)->vy,(tracer+i)->vz);
     }
    }
#endif

    /* check if we have to dump */
    if(itime%((int)(property.time_dump_lagr/property.time_dt))==0){

#ifdef OUTPUT_H5
    /* First check how many particles in each processor and compute offset */
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);



    /* then alloc space */
    aux  = (my_double*) malloc(sizeof(my_double)*npart); 
    if(aux == NULL){ fprintf(stderr,"Not enough memory to allocate aux field t\n"); exit(-1);}

    hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);

                /* Create a new file using default properties */
		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);   
		
		file_id = H5Fcreate(H5FILE_NAME_PARTICLE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		group   = H5Gcreate (file_id, "/lagrange", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

		H5Pclose(plist_id);

		property_id  = H5Pcreate(H5P_DATASET_CREATE);


                /* Create the data space for the dataset. */
                dims[0] = (int)property.particle_number;

		filespace = H5Screate_simple(RANK, dims, NULL);   
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
		count[0] = npart;
		offset[0] = name_offset;

	        memspace = H5Screate_simple(RANK, count, NULL); 	

    /*
     * Select hyperslab in the file.
     */
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

		xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);                       

		/* WRITE PARTICLE NAME */
		dataset_id = H5Dcreate(group, "name", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->name;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		/* WRITE PARTICLE DRAG RESPONSE TIME */
		dataset_id = H5Dcreate(group, "tau_drag", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->tau_drag;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
                
		/* WRITE PARTICLE POSITIONS */
		dataset_id = H5Dcreate(group, "x", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->x;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "y", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->y;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "z", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->z;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		/* WRITE PARTICLE VELOCITIES */
		dataset_id = H5Dcreate(group, "vx", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->vx;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "vy", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->vy;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "vz", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->vz;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		/* WRITE PARTICLE ACCELERATIONS */
		dataset_id = H5Dcreate(group, "ax", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->ax;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "ay", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->ay;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "az", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->az;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

#ifdef LAGRANGE_GRAVITY
 #ifdef LAGRANGE_GRAVITY_VARIABLE
		/* GRAVITY ADJUSTABLE COEFFICIENT */
		dataset_id = H5Dcreate(group, "gravity_coeff", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->gravity_coeff;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
 #endif
#endif

#ifdef LAGRANGE_GRADIENT
		/* FLUID VELOCITY GRADIENT */
		dataset_id = H5Dcreate(group, "dx_ux", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dx_ux;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dy_ux", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dy_ux;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dz_ux", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dz_ux;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dx_uy", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dx_uy;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dy_uy", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dy_uy;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dz_uy", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dz_uy;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dx_uz", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dx_uz;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dy_uz", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dy_uz;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dz_uz", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dz_uz;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

 #ifdef LAGRANGE_ADDEDMASS
		/* ADDED MASS BETA COEFFICIENT */
		dataset_id = H5Dcreate(group, "beta_coeff", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->beta_coeff;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
 #endif
 #ifdef LAGRANGE_ORIENTATION
		/* ORIENTATION VECTOR */
		dataset_id = H5Dcreate(group, "px", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->px;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "py", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->py;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "pz", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->pz;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
		/* VELOCITY ROTATION VECTOR */
		dataset_id = H5Dcreate(group, "dt_px", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dt_px;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dt_py", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dt_py;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dt_pz", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dt_pz;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
  #ifdef LAGRANGE_ORIENTATION_JEFFREY
		/* PARTICLE ASPECT RATIO */
		dataset_id = H5Dcreate(group, "aspect_ratio", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->aspect_ratio;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
   #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
		/* GYROTACTIC PARAMETER */
		dataset_id = H5Dcreate(group, "gyrotaxis_velocity", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->gyrotaxis_velocity;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);	
   #endif
  #endif		
  #ifdef LAGRANGE_ORIENTATION_ACTIVE
		dataset_id = H5Dcreate(group, "swim_velocity", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->swim_velocity;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
   #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
		dataset_id = H5Dcreate(group, "critical_shear_rate", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->critical_shear_rate;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "shear_rate", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->shear_rate;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "jump_time", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->jump_time;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
   #endif /* LAGRANGE_ORIENTATION_ACTIVE_JUMP */
  #endif /* LAGRANGE_ORIENTATION_ACTIVE */
 #endif
#endif

#ifdef LB_TEMPERATURE
		/* WRITE PARTICLE TEMPERATURE */		
		dataset_id = H5Dcreate(group, "temperature", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->t;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
 #ifdef LAGRANGE_GRADIENT
		/* TEMPERATURE GRADIENT */
		dataset_id = H5Dcreate(group, "dx_t", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dx_t;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dy_t", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dy_t;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dz_t", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dz_t;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
 #endif		
#endif

#ifdef LB_SCALAR
		/* WRITE PARTICLE SCALAR */
		dataset_id = H5Dcreate(group, "scalar", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->s;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
 #ifdef LAGRANGE_GRADIENT
		/* SCALAR GRADIENT */
		dataset_id = H5Dcreate(group, "dx_s", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dx_s;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dy_s", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dy_s;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);

		dataset_id = H5Dcreate(group, "dz_s", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
		for(i=0;i<npart;i++) aux[i]=(tracer + i)->dz_s;
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, aux);
                status = H5Dclose(dataset_id);
 #endif		
#endif

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

  /* free scalar auxiliary field */
  free(aux);

 /* create the file names */
  sprintf(NEW_H5FILE_NAME,"%s/particle_%d.h5",OutDir,itime);

  /* we rename the file */
  if(ROOT) rename(H5FILE_NAME_PARTICLE, NEW_H5FILE_NAME);


  /* Xml file */
  if(ROOT){
  sprintf(XMF_FILE_NAME,"%s/particle_%d.xmf" ,OutDir,itime);
  sprintf(NEW_H5FILE_NAME,"particle_%d.h5",itime);
  size=sizeof(my_double);
  

                fout = fopen(XMF_FILE_NAME,"w");
		              
                fprintf(fout,"<?xml version=\"1.0\" ?>\n");
                fprintf(fout,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                fprintf(fout,"<Xdmf Version=\"2.0\">\n");
                fprintf(fout,"<Domain>\n");
                fprintf(fout,"<Grid Name=\"Points\" GridType=\"Uniform\">\n");
                fprintf(fout,"<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n",np);
                fprintf(fout,"<Geometry GeometryType=\"X_Y_Z\">\n");
                
                fprintf(fout,"<DataItem Name=\"x\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/x\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                
                fprintf(fout,"<DataItem Name=\"y\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/y\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");

                fprintf(fout,"<DataItem Name=\"z\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/z\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");

                fprintf(fout,"</Geometry>\n");                

		/* name */
                fprintf(fout,"<Attribute Name=\"name\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/name\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  

		/* tau_drag */
                fprintf(fout,"<Attribute Name=\"tau_drag\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/tau_drag\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");   

		/* velocity as a vector */
                fprintf(fout,"<Attribute Name=\"velocity\" AttributeType=\"Vector\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d 3\" \n   Function=\"JOIN($0 , $1, $2)\">\n",np);
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/vx\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/vy\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/vz\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  

		/* acceleration as a vector */
                fprintf(fout,"<Attribute Name=\"acceleration\" AttributeType=\"Vector\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d 3\" \n   Function=\"JOIN($0 , $1, $2)\">\n",np);
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/ax\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/ay\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/az\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  

#ifdef LAGRANGE_GRAVITY
 #ifdef LAGRANGE_GRAVITY_VARIABLE
		/* gravity_coeff */
                fprintf(fout,"<Attribute Name=\"gravity_coeff\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/gravity_coeff\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n"); 
 #endif
#endif

#ifdef LAGRANGE_GRADIENT
		/* fluid velocity gradient */
                fprintf(fout,"<Attribute Name=\"fluid velocity gradient\" AttributeType=\"Tensor\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d 9\" \n   Function=\"JOIN($0,$1,$2,$3,$4,$5,$6,$7,$8)\">\n",np);
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dx_ux\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dy_ux\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dz_ux\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dx_uy\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dy_uy\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dz_uy\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dx_uz\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dy_uz\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dz_uz\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
 #ifdef LAGRANGE_ADDEDMASS
		/* beta_coeff */
                fprintf(fout,"<Attribute Name=\"beta_coeff\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/beta_coeff\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n"); 
 #endif
 #ifdef LAGRANGE_ORIENTATION 
		/* orientation vector */
                fprintf(fout,"<Attribute Name=\"orientation\" AttributeType=\"Vector\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d 3\" \n   Function=\"JOIN($0 , $1, $2)\">\n",np);
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/px\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/py\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/pz\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");
  		/* orientation velocity vector */
                fprintf(fout,"<Attribute Name=\"angular velocity\" AttributeType=\"Vector\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d 3\" \n   Function=\"JOIN($0 , $1, $2)\">\n",np);
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dt_px\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dt_py\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dt_pz\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");
  #ifdef LAGRANGE_ORIENTATION_JEFFREY
		/* aspect ratio */
                fprintf(fout,"<Attribute Name=\"aspect_ratio\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/aspect_ratio\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
   #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
		/* gyrotaxis velocity parameter */
                fprintf(fout,"<Attribute Name=\"gyrotaxis_velocity\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/gyrotaxis_velocity\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
   #endif
  #endif		
  #ifdef LAGRANGE_ORIENTATION_ACTIVE
		/* swim velocity */
                fprintf(fout,"<Attribute Name=\"swim_velocity\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/swim_velocity\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
   #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
		/* critical_shear_rate */
                fprintf(fout,"<Attribute Name=\"critical_shear_rate\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/critical_shear_rate\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
		/* shear_rate */
                fprintf(fout,"<Attribute Name=\"shear_rate\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/shear_rate\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
		/* jump_time */
                fprintf(fout,"<Attribute Name=\"jump_time\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/jump_time\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
   #endif
  #endif
 #endif       
#endif

#ifdef LB_TEMPERATURE
		/* temperature at particle position */
                fprintf(fout,"<Attribute Name=\"temperature\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/temperature\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");   
 #ifdef LAGRANGE_GRADIENT
		/* temperature gradient at particle position */
                fprintf(fout,"<Attribute Name=\"temperature gradient\" AttributeType=\"Vector\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d 3\" \n   Function=\"JOIN($0 , $1, $2)\">\n",np);
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dx_t\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dy_t\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dz_t\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
 #endif       
#endif

#ifdef LB_SCALAR
		/* scalar at particle position */
                fprintf(fout,"<Attribute Name=\"scalar\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/scalar\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n"); 
 #ifdef LAGRANGE_GRADIENT
		/* scalar gradient at particle position */
                fprintf(fout,"<Attribute Name=\"scalar gradient\" AttributeType=\"Vector\" Center=\"Node\"> \n");
                fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d 3\" \n   Function=\"JOIN($0 , $1, $2)\">\n",np);
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dx_s\n",NEW_H5FILE_NAME); 
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dy_s\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"<DataItem Dimensions=\"%d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", np, size);
                fprintf(fout,"%s:/lagrange/dz_s\n",NEW_H5FILE_NAME);
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</DataItem>\n");
                fprintf(fout,"</Attribute>\n");  
 #endif                
#endif
                
                fprintf(fout,"</Grid>\n");
                fprintf(fout,"</Domain>\n");
                fprintf(fout,"</Xdmf>\n");
		                
                fclose(fout);
  }/* end of if root */

#endif /* end of OUTPUT_H5 */

    }/* end of if on itime , check if we have to dump particle */

}/* end of function output_particles */


/* dump lagrangian averages */
void dump_particle_averages(){

  output_particle *out_particle_local, *out_particle_all;
  int ipart, type;
  double norm, *counter_local, *counter_all;
  FILE *fout;
  char fname[128];

if(itime%((int)(property.time_dump_diagn/property.time_dt))==0){

  /* malloc */
  out_particle_local  = (output_particle*) malloc(sizeof(output_particle)*(int)property.particle_types);
  out_particle_all  = (output_particle*) malloc(sizeof(output_particle)*(int)property.particle_types);
  counter_local = (my_double*) malloc(sizeof(my_double)*(int)property.particle_types);
  counter_all = (my_double*) malloc(sizeof(my_double)*(int)property.particle_types);

  /* set to zero */
 for (type=0;type<(int)property.particle_types;type++) {
   counter_local[type] = counter_all[type] = 0.0;

   out_particle_local[type].vx = out_particle_all[type].vx = 0.0;
   out_particle_local[type].vy = out_particle_all[type].vy = 0.0;
   out_particle_local[type].vz = out_particle_all[type].vz = 0.0;

   out_particle_local[type].vx2 = out_particle_all[type].vx2 = 0.0;
   out_particle_local[type].vy2 = out_particle_all[type].vy2 = 0.0;
   out_particle_local[type].vz2 = out_particle_all[type].vz2 = 0.0;

   out_particle_local[type].vx4 = out_particle_all[type].vx4 = 0.0;
   out_particle_local[type].vy4 = out_particle_all[type].vy4 = 0.0;
   out_particle_local[type].vz4 = out_particle_all[type].vz4 = 0.0;

   out_particle_local[type].ax = out_particle_all[type].ax = 0.0;
   out_particle_local[type].ay = out_particle_all[type].ay = 0.0;
   out_particle_local[type].az = out_particle_all[type].az = 0.0;

   out_particle_local[type].ax2 = out_particle_all[type].ax2 = 0.0;
   out_particle_local[type].ay2 = out_particle_all[type].ay2 = 0.0;
   out_particle_local[type].az2 = out_particle_all[type].az2 = 0.0;

   out_particle_local[type].ax4 = out_particle_all[type].ax4 = 0.0;
   out_particle_local[type].ay4 = out_particle_all[type].ay4 = 0.0;
   out_particle_local[type].az4 = out_particle_all[type].az4 = 0.0;
 #ifdef LAGRANGE_GRADIENT
  #ifdef LAGRANGE_ORIENTATION
   out_particle_local[type].dt_px = out_particle_all[type].dt_px = 0.0;
   out_particle_local[type].dt_py = out_particle_all[type].dt_py = 0.0;
   out_particle_local[type].dt_pz = out_particle_all[type].dt_pz = 0.0;

   out_particle_local[type].dt_px2 = out_particle_all[type].dt_px2 = 0.0;
   out_particle_local[type].dt_py2 = out_particle_all[type].dt_py2 = 0.0;
   out_particle_local[type].dt_pz2 = out_particle_all[type].dt_pz2 = 0.0;

   out_particle_local[type].dt_px4 = out_particle_all[type].dt_px4 = 0.0;
   out_particle_local[type].dt_py4 = out_particle_all[type].dt_py4 = 0.0;
   out_particle_local[type].dt_pz4 = out_particle_all[type].dt_pz4 = 0.0;
  #endif 
 #endif
 }   

/* Begin loop on particles */
 for (ipart=0;ipart<npart;ipart++) {

 type = ((int)(tracer+ipart)->name)%(int)property.particle_types;

 counter_local[type] += 1.0;

 out_particle_local[type].vx += (tracer+ipart)->vx;
 out_particle_local[type].vy += (tracer+ipart)->vy;
 out_particle_local[type].vz += (tracer+ipart)->vz;

 out_particle_local[type].vx2 += (tracer+ipart)->vx * (tracer+ipart)->vx;
 out_particle_local[type].vy2 += (tracer+ipart)->vy * (tracer+ipart)->vy;
 out_particle_local[type].vz2 += (tracer+ipart)->vz * (tracer+ipart)->vz;

 out_particle_local[type].vx4 += pow((tracer+ipart)->vx,4.0); 
 out_particle_local[type].vy4 += pow((tracer+ipart)->vy,4.0); 
 out_particle_local[type].vz4 += pow((tracer+ipart)->vz,4.0); 

 out_particle_local[type].ax += (tracer+ipart)->ax;
 out_particle_local[type].ay += (tracer+ipart)->ay;
 out_particle_local[type].az += (tracer+ipart)->az;

 out_particle_local[type].ax2 += (tracer+ipart)->ax * (tracer+ipart)->ax;
 out_particle_local[type].ay2 += (tracer+ipart)->ay * (tracer+ipart)->ay;
 out_particle_local[type].az2 += (tracer+ipart)->az * (tracer+ipart)->az;

 out_particle_local[type].ax4 += pow((tracer+ipart)->ax,4.0); 
 out_particle_local[type].ay4 += pow((tracer+ipart)->ay,4.0); 
 out_particle_local[type].az4 += pow((tracer+ipart)->az,4.0); 

 #ifdef LAGRANGE_GRADIENT
  #ifdef LAGRANGE_ORIENTATION
 out_particle_local[type].dt_px += (tracer+ipart)->dt_px;
 out_particle_local[type].dt_py += (tracer+ipart)->dt_py;
 out_particle_local[type].dt_pz += (tracer+ipart)->dt_pz;

 out_particle_local[type].dt_px2 += (tracer+ipart)->dt_px * (tracer+ipart)->dt_px;
 out_particle_local[type].dt_py2 += (tracer+ipart)->dt_py * (tracer+ipart)->dt_py;
 out_particle_local[type].dt_pz2 += (tracer+ipart)->dt_pz * (tracer+ipart)->dt_pz;

 out_particle_local[type].dt_px4 += pow((tracer+ipart)->dt_px,4.0);
 out_particle_local[type].dt_py4 += pow((tracer+ipart)->dt_py,4.0);
 out_particle_local[type].dt_pz4 += pow((tracer+ipart)->dt_pz,4.0);
  #endif
 #endif

 }/* end of for on ipart */

   /* Sum all */
 //if(ROOT) for(type=0;type<(int)property.particle_types;type++) fprintf(stdout,"before %e %e type %d\n", out_particle_local[type].vx2, out_particle_all[type].vx2, type);
 MPI_Allreduce(out_particle_local, out_particle_all, (int)property.particle_types, MPI_output_particle_type, MPI_SUM_output_particle, MPI_COMM_WORLD );
 //if(ROOT) for(type=0;type<(int)property.particle_types;type++) fprintf(stdout,"after %e %e type %d\n", out_particle_local[type].vx2, out_particle_all[type].vx2, type);
 MPI_Allreduce(counter_local, counter_all, (int)property.particle_types, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

 for (type=0;type<(int)property.particle_types;type++) {

   /* Normalization */
   norm = 1.0/counter_all[type];

   //if(ROOT)fprintf(stdout,"counter %f type %d\n", counter_all[type], type);
   
   out_particle_all[type].vx *= norm;
   out_particle_all[type].vy *= norm;
   out_particle_all[type].vz *= norm;

   out_particle_all[type].vx2 *= norm;
   out_particle_all[type].vy2 *= norm;
   out_particle_all[type].vz2 *= norm;

   out_particle_all[type].vx4 *= norm;
   out_particle_all[type].vy4 *= norm;
   out_particle_all[type].vz4 *= norm;

   out_particle_all[type].ax *= norm;
   out_particle_all[type].ay *= norm;
   out_particle_all[type].az *= norm;

   out_particle_all[type].ax2 *= norm;
   out_particle_all[type].ay2 *= norm;
   out_particle_all[type].az2 *= norm;

   out_particle_all[type].ax4 *= norm;
   out_particle_all[type].ay4 *= norm;
   out_particle_all[type].az4 *= norm;
 #ifdef LAGRANGE_GRADIENT
  #ifdef LAGRANGE_ORIENTATION
   out_particle_all[type].dt_px *= norm;
   out_particle_all[type].dt_py *= norm;
   out_particle_all[type].dt_pz *= norm;

   out_particle_all[type].dt_px2 *= norm;
   out_particle_all[type].dt_py2 *= norm;
   out_particle_all[type].dt_pz2 *= norm;

   out_particle_all[type].dt_px4 *= norm;
   out_particle_all[type].dt_py4 *= norm;
   out_particle_all[type].dt_pz4 *= norm;
  #endif 
 #endif
 }


if(ROOT){
    sprintf(fname,"particle_averages.dat");
    fout = fopen(fname,"a");    
 for (type=0;type<(int)property.particle_types;type++) {
   fprintf(fout,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ",time_now, (double)type,
	   (double)out_particle_all[type].vx ,(double)out_particle_all[type].vy ,(double)out_particle_all[type].vz,
	   (double)out_particle_all[type].vx2,(double)out_particle_all[type].vy2,(double)out_particle_all[type].vz2,
	   (double)out_particle_all[type].vx4,(double)out_particle_all[type].vy4,(double)out_particle_all[type].vz4,
	   (double)out_particle_all[type].ax ,(double)out_particle_all[type].ay ,(double)out_particle_all[type].az,
	   (double)out_particle_all[type].ax2,(double)out_particle_all[type].ay2,(double)out_particle_all[type].az2,
	   (double)out_particle_all[type].ax4,(double)out_particle_all[type].ay4,(double)out_particle_all[type].az4 );
 #ifdef LAGRANGE_GRADIENT
  #ifdef LAGRANGE_ORIENTATION
   fprintf(fout,"%e %e %e %e %e %e %e %e %e ",
	   (double)out_particle_all[type].dt_px ,(double)out_particle_all[type].dt_py  ,(double)out_particle_all[type].dt_pz,
	   (double)out_particle_all[type].dt_px2,(double)out_particle_all[type].dt_py2 ,(double)out_particle_all[type].dt_pz2,
	   (double)out_particle_all[type].dt_px4, (double)out_particle_all[type].dt_py4,(double)out_particle_all[type].dt_pz4);
  #endif 
 #endif
   fprintf(fout,"\n");
 }
    fclose(fout);
 }/* ROOT */

  free(out_particle_local);
  free(out_particle_all);
  free(counter_local);
  free(counter_all);

 }/* if on time */

}/* end of dump_particle_averages */


/* advance in time particles and assign them to the right processors */
void move_particles(){

  int ipart,i,j;
  int npart_here,npart_there,all_npart_there,all_npart;  
  point_particle part;
  int *displs,*rcounts;
  my_double invtau;
  vector Dt_u;
  vector vec, v_old;
#ifdef LAGRANGE_ORIENTATION
  my_double matA[3][3],matS[3][3],matW[3][3];
  my_double scalOSO, f_alpha, alpha,norm , gyro;
  my_double vecF[3],vecFold[3],vecP[3],vecTMP[3],vecA[3];
 #ifdef LAGRANGE_ORIENTATION_DIFFUSION
  my_double vec_xi[3];
  my_double two_d_r;
  my_double matI[3][3];
 #endif 
  #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
  my_double shear_rate,jump_time_duration ,velocity_amplitude;
  #endif
#endif

#ifdef LAGRANGE_ORIENTATION
 #ifdef LAGRANGE_ORIENTATION_DIFFUSION
  /* Define the identity matrix */
	     matI[0][0] = matI[1][1] = matI[2][2] = 1.0;
	     matI[0][1] = matI[0][2] = matI[1][2] = 0.0;
	     matI[1][0] = matI[2][0] = matI[2][1] = 0.0;
 #endif 
#endif

  //fprintf(stderr,"me %d I am here, npart %d time %g\n",me, npart,time_now);

/* Begin loop on particles */
 for (ipart=0;ipart<npart;ipart++) {

   /* if we just have tracers */
   if((tracer+ipart)->tau_drag == 0.0){

   /* copy fluid velocity into particle velocity NOTE that this is true only for tracers */
   (tracer+ipart)->vx = (tracer+ipart)->ux;
   (tracer+ipart)->vy = (tracer+ipart)->uy;
   (tracer+ipart)->vz = (tracer+ipart)->uz;

#ifdef LAGRANGE_ORIENTATION 
 #ifdef LAGRANGE_ORIENTATION_ACTIVE 
   /* if the particle is alive there is an extra velocity to add */

  #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP

   #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP_TEMPERATURE
    #ifdef LB_TEMPERATURE
   /* the particle can react to the temperature value and decide if to jump */
   /* put temperature value in in the shear_rate */
   shear_rate = (tracer+ipart)->t;
    #endif
   #else
  /* the particle can react to flow gradients and decide if to jump */
  /* compute |s| */
  (tracer+ipart)->shear_rate =  shear_rate = sqrt(
  ((tracer+ipart)->dx_ux) *  ((tracer+ipart)->dx_ux) +
  ((tracer+ipart)->dy_uy) *  ((tracer+ipart)->dz_uy) +
  ((tracer+ipart)->dz_uz) *  ((tracer+ipart)->dz_uz) +
  0.5*( ((tracer+ipart)->dy_ux) + ((tracer+ipart)->dx_uy) )*( ((tracer+ipart)->dy_ux) + ((tracer+ipart)->dx_uy) ) +  
  0.5*( ((tracer+ipart)->dz_ux) + ((tracer+ipart)->dx_uz) )*( ((tracer+ipart)->dz_ux) + ((tracer+ipart)->dx_uz) ) +  
  0.5*( ((tracer+ipart)->dy_uz) + ((tracer+ipart)->dz_uy) )*( ((tracer+ipart)->dy_uz) + ((tracer+ipart)->dz_uy) ) );
   #endif

 jump_time_duration = - ((tracer+ipart)->jump_time)*log(0.01);

    if((tracer+ipart)->time_from_jump > jump_time_duration){
      /* set to zero the velocity amplitude */
      velocity_amplitude = 0.0;
      /* but we are in the condition to jump */

#ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP_HIGHPASS
	/* avoids low shear rate values */
      if( shear_rate  < (tracer+ipart)->critical_shear_rate ){
#else 
	/* avoids high shear rate values */	
      if( shear_rate  > (tracer+ipart)->critical_shear_rate ){
#endif


	/* reset the time from jump */ 	
	(tracer+ipart)->time_from_jump = 0.0;
        /* set initial velocity amplitude */
        velocity_amplitude = (tracer+ipart)->swim_velocity; 

     /* old version 
        // generarate a new random vector    
        vec = random_vector();
        // retrive last particle velocity 
        v_old.x = (tracer+ipart)->vx_old;
        v_old.y = (tracer+ipart)->vy_old;
        v_old.z = (tracer+ipart)->vz_old;
        // choose vec direction in the same emisphere of the particle velocity 
        if( scalar_product(vec,v_old) < 0.0 ) vec = vector_scale( -1.0 , vec);
        (tracer+ipart)->px = vec.x;
        (tracer+ipart)->py = vec.y;
        (tracer+ipart)->pz = vec.z;
     */
	/* new version : take for jump the present orientation */  
        (tracer+ipart)->px_jump = (tracer+ipart)->px;
        (tracer+ipart)->py_jump = (tracer+ipart)->py;
        (tracer+ipart)->pz_jump = (tracer+ipart)->pz;	
      } /* end if on shear rate */

    }else{
      /* hence (tracer+ipart)->time_from_jump < jump_time_duration  : we are already jumping */
      velocity_amplitude = (tracer+ipart)->swim_velocity * exp( -((tracer+ipart)->time_from_jump) / ((tracer+ipart)->jump_time) ); 
    }

      /* add jump part to the particle velocity */
      /* old version 
      (tracer+ipart)->vx += velocity_amplitude*((tracer+ipart)->px);
      (tracer+ipart)->vy += velocity_amplitude*((tracer+ipart)->py);
      (tracer+ipart)->vz += velocity_amplitude*((tracer+ipart)->pz);
      */
      /* new version : we jump with jeffrey */
      (tracer+ipart)->vx += velocity_amplitude*((tracer+ipart)->px_jump);
      (tracer+ipart)->vy += velocity_amplitude*((tracer+ipart)->py_jump);
      (tracer+ipart)->vz += velocity_amplitude*((tracer+ipart)->pz_jump);

      /* increase time from jump */
      (tracer+ipart)->time_from_jump += property.time_dt;

   #else /* NOT LAGRANGE_ORIENTATION_ACTIVE_JUMP */
  
      /* it just swim with direction p evolved by jeffrey,random, etc. (if they are activated) */
   (tracer+ipart)->vx += ((tracer+ipart)->swim_velocity)*((tracer+ipart)->px);
   (tracer+ipart)->vy += ((tracer+ipart)->swim_velocity)*((tracer+ipart)->py);
   (tracer+ipart)->vz += ((tracer+ipart)->swim_velocity)*((tracer+ipart)->pz);
   // fprintf(stderr,"%g %g %g %g \n",(tracer+ipart)->name , (tracer+ipart)->vx , (tracer+ipart)->swim_velocity , (tracer+ipart)->px );
  //  fprintf(stderr,"%g %g %g %g \n",(tracer+ipart)->name , (tracer+ipart)->vy , (tracer+ipart)->swim_velocity , (tracer+ipart)->py );
  #endif /* NOT LAGRANGE_ORIENTATION_ACTIVE_JUMP */
 #endif /* LAGRANGE_ORIENTATION_ACTIVE */
#endif /* LAGRANGE_ORIENTATION */

  if(itime==0 && resume==0){ 
    (tracer+ipart)->vx_old = (tracer+ipart)->vx;
    (tracer+ipart)->vy_old = (tracer+ipart)->vy;
    (tracer+ipart)->vz_old = (tracer+ipart)->vz;
  }
  /* Compute tracer acceleration : if v=u as here, then a = D_t u */
  (tracer+ipart)->ax = ((tracer+ipart)->vx - (tracer+ipart)->vx_old )/property.time_dt;
  (tracer+ipart)->ay = ((tracer+ipart)->vy - (tracer+ipart)->vy_old )/property.time_dt;
  (tracer+ipart)->az = ((tracer+ipart)->vz - (tracer+ipart)->vz_old )/property.time_dt;


   if(itime==0 && resume==0){
   /* Explicit Euler 1st order */
   (tracer+ipart)->x += property.time_dt*(tracer+ipart)->vx;
   (tracer+ipart)->y += property.time_dt*(tracer+ipart)->vy;
   (tracer+ipart)->z += property.time_dt*(tracer+ipart)->vz;
   }else{
   /* Adams-Bashforth 2nd order */
   (tracer+ipart)->x += property.time_dt*0.5*(3.0*(tracer+ipart)->vx - (tracer+ipart)->vx_old);
   (tracer+ipart)->y += property.time_dt*0.5*(3.0*(tracer+ipart)->vy - (tracer+ipart)->vy_old);
   (tracer+ipart)->z += property.time_dt*0.5*(3.0*(tracer+ipart)->vz - (tracer+ipart)->vz_old);
   }
   /* copy particle velocity in old */
   (tracer+ipart)->vx_old = (tracer+ipart)->vx; 
   (tracer+ipart)->vy_old = (tracer+ipart)->vy; 
   (tracer+ipart)->vz_old = (tracer+ipart)->vz;  

   /* copy fluid velocity in old */
   (tracer+ipart)->ux_old = (tracer+ipart)->ux; 
   (tracer+ipart)->uy_old = (tracer+ipart)->uy; 
   (tracer+ipart)->uz_old = (tracer+ipart)->uz; 

   }/* end if on fluid tracer */

   /* With drag force */ 
   if((tracer+ipart)->tau_drag != 0.0){
  
   invtau = 1.0 / (tracer+ipart)->tau_drag;
   (tracer+ipart)->ax = ((tracer+ipart)->ux - (tracer+ipart)->vx)*invtau;
   (tracer+ipart)->ay = ((tracer+ipart)->uy - (tracer+ipart)->vy)*invtau;
   (tracer+ipart)->az = ((tracer+ipart)->uz - (tracer+ipart)->vz)*invtau;


#ifdef LAGRANGE_GRAVITY /* note that works only if LB_FORCING_GRAVITY is defined */
   /*  add: -g to acceleration */ 
 #ifdef LAGRANGE_GRAVITY_VARIABLE
     (tracer+ipart)->ax -= (tracer+ipart)->gravity_coeff*property.gravity_x;
     (tracer+ipart)->ay -= (tracer+ipart)->gravity_coeff*property.gravity_y;
     (tracer+ipart)->az -= (tracer+ipart)->gravity_coeff*property.gravity_z; 
 #else   
     (tracer+ipart)->ax -= property.gravity_x;
     (tracer+ipart)->ay -= property.gravity_y;
     (tracer+ipart)->az -= property.gravity_z; 
 #endif
#endif

#ifdef LAGRANGE_ADDEDMASS
  /* With Added mass */ 
   if((tracer+ipart)->beta_coeff != 0.0){

#ifdef LAGRANGE_GRAVITY
  /*  add also: -\beta*g to acceleration */ 
 #ifdef LAGRANGE_GRAVITY_VARIABLE
     (tracer+ipart)->ax -= (  - (tracer+ipart)->beta_coeff )*(tracer+ipart)->gravity_coeff*property.gravity_x;
     (tracer+ipart)->ay -= (  - (tracer+ipart)->beta_coeff )*(tracer+ipart)->gravity_coeff*property.gravity_y;
     (tracer+ipart)->az -= (  - (tracer+ipart)->beta_coeff )*(tracer+ipart)->gravity_coeff*property.gravity_z; 
 #else    
     (tracer+ipart)->ax -= (  - (tracer+ipart)->beta_coeff )*property.gravity_x;
     (tracer+ipart)->ay -= (  - (tracer+ipart)->beta_coeff )*property.gravity_y;
     (tracer+ipart)->az -= (  - (tracer+ipart)->beta_coeff )*property.gravity_z; 
 #endif    
#endif

  if(itime==0 && resume==0){ 
    (tracer+ipart)->ux_old = (tracer+ipart)->ux;
    (tracer+ipart)->uy_old = (tracer+ipart)->uy;
    (tracer+ipart)->uz_old = (tracer+ipart)->uz;
  }

   /* Here I will write the computation of the fluid material derivative */
   Dt_u.x = ((tracer+ipart)->ux - (tracer+ipart)->ux_old )/property.time_dt
          + ((tracer+ipart)->ux - (tracer+ipart)->vx)*(tracer+ipart)->dx_ux 
          + ((tracer+ipart)->uy - (tracer+ipart)->vy)*(tracer+ipart)->dy_ux 
          + ((tracer+ipart)->uz - (tracer+ipart)->vz)*(tracer+ipart)->dz_ux;

   Dt_u.y = ((tracer+ipart)->uy - (tracer+ipart)->uy_old )/property.time_dt
          + ((tracer+ipart)->ux - (tracer+ipart)->vx)*(tracer+ipart)->dx_uy 
          + ((tracer+ipart)->uy - (tracer+ipart)->vy)*(tracer+ipart)->dy_uy 
          + ((tracer+ipart)->uz - (tracer+ipart)->vz)*(tracer+ipart)->dz_uy;

   Dt_u.z = ((tracer+ipart)->uz - (tracer+ipart)->uz_old )/property.time_dt         
	  + ((tracer+ipart)->ux - (tracer+ipart)->vx)*(tracer+ipart)->dx_uz 
          + ((tracer+ipart)->uy - (tracer+ipart)->vy)*(tracer+ipart)->dy_uz 
          + ((tracer+ipart)->uz - (tracer+ipart)->vz)*(tracer+ipart)->dz_uz;

   (tracer+ipart)->ax += Dt_u.x*(tracer+ipart)->beta_coeff;
   (tracer+ipart)->ay += Dt_u.y*(tracer+ipart)->beta_coeff;
   (tracer+ipart)->az += Dt_u.z*(tracer+ipart)->beta_coeff;

   }/* end of if on addedd mass */
#endif

   if(itime==0 && resume==0){
     (tracer+ipart)->vx += property.time_dt*(tracer+ipart)->ax;
     (tracer+ipart)->vy += property.time_dt*(tracer+ipart)->ay;
     (tracer+ipart)->vz += property.time_dt*(tracer+ipart)->az;
   }else{
     (tracer+ipart)->vx += property.time_dt*0.5*(3.0*(tracer+ipart)->ax - (tracer+ipart)->ax_old);
     (tracer+ipart)->vy += property.time_dt*0.5*(3.0*(tracer+ipart)->ay - (tracer+ipart)->ay_old);
     (tracer+ipart)->vz += property.time_dt*0.5*(3.0*(tracer+ipart)->az - (tracer+ipart)->az_old);
   }

   (tracer+ipart)->ax_old = (tracer+ipart)->ax; 
   (tracer+ipart)->ay_old = (tracer+ipart)->ay; 
   (tracer+ipart)->az_old = (tracer+ipart)->az; 

   if(itime==0 && resume==0){
     (tracer+ipart)->x += property.time_dt*(tracer+ipart)->vx; 
     (tracer+ipart)->y += property.time_dt*(tracer+ipart)->vy;
     (tracer+ipart)->z += property.time_dt*(tracer+ipart)->vz;
   }else{
     (tracer+ipart)->x += property.time_dt*0.5*(3.0*(tracer+ipart)->vx - (tracer+ipart)->vx_old);
     (tracer+ipart)->y += property.time_dt*0.5*(3.0*(tracer+ipart)->vy - (tracer+ipart)->vy_old);
     (tracer+ipart)->z += property.time_dt*0.5*(3.0*(tracer+ipart)->vz - (tracer+ipart)->vz_old);
   }
   (tracer+ipart)->vx_old = (tracer+ipart)->vx; 
   (tracer+ipart)->vy_old = (tracer+ipart)->vy; 
   (tracer+ipart)->vz_old = (tracer+ipart)->vz; 

   (tracer+ipart)->ux_old = (tracer+ipart)->ux; 
   (tracer+ipart)->uy_old = (tracer+ipart)->uy; 
   (tracer+ipart)->uz_old = (tracer+ipart)->uz; 

   }/* end of if on tau_drag different from zero */

#ifdef LAGRANGE_ORIENTATION

	      /* assign P vector */
     vecP[0] = (tracer+ipart)->px;
     vecP[1] = (tracer+ipart)->py;
     vecP[2] = (tracer+ipart)->pz;
              /* assign the last dP /dt  vector */
     vecFold[0] = (tracer+ipart)->dt_px;
     vecFold[1] = (tracer+ipart)->dt_py;
     vecFold[2] = (tracer+ipart)->dt_pz;

 #ifdef LAGRANGE_ORIENTATION_JEFFREY
 /* Here we implement Jeffrey equation */

              /* aspect ratio factor */
   alpha = (tracer+ipart)->aspect_ratio;
   f_alpha = (alpha*alpha-1.0)/(1.0+alpha*alpha);


	      /* velocity gradient matrix */
	      matA[0][0]=(tracer+ipart)->dx_ux ; matA[0][1]=(tracer+ipart)->dy_ux; matA[0][2]=(tracer+ipart)->dz_ux;
	      matA[1][0]=(tracer+ipart)->dx_uy ; matA[1][1]=(tracer+ipart)->dy_uy; matA[1][2]=(tracer+ipart)->dz_uy;
	      matA[2][0]=(tracer+ipart)->dx_uz ; matA[2][1]=(tracer+ipart)->dy_uz; matA[2][2]=(tracer+ipart)->dz_uz;	 

	      /* Compute Sij */
	      /*make simmetric*/
	      for (i=0; i<3; i++)
                for (j=0; j<3; j++){
                  matS[i][j] = 0.5*(matA[i][j]+matA[j][i]);
		}

	      /* Compute Wij */
	      /*make simmetric*/
	      for (i=0; i<3; i++)
                for (j=0; j<3; j++){
                  matW[i][j] = 0.5*(matA[i][j]-matA[j][i]);
		}

	      /* multiply S by the aspect ratio stretching factor */
	      for (i=0; i<3; i++)
                for (j=0; j<3; j++){
                  matS[i][j] *= f_alpha;
		}

  #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
	       /* gravitational gyrotaxis : the stretched S matrix has an extra term -1/(2*v0) * g_i p_j      */	      
	      if((tracer+ipart)->gyrotaxis_velocity !=0){
		gyro = - 0.5 / (tracer+ipart)->gyrotaxis_velocity;
   #ifdef LAGRANGE_GRAVITY	      
		/* the full term */	      
		vecA[0] = gyro * (- property.gravity_x  + (tracer+ipart)->ax);
		vecA[1] = gyro * (- property.gravity_y  + (tracer+ipart)->ay);
		vecA[2] = gyro * (- property.gravity_z  + (tracer+ipart)->az);
	      
   #else
		/* just for a test: fixed vector along z. Like g_x=0, g_y=0, g_z=-1.0 */
		vecA[0] = gyro * 0.0;  
		vecA[1] = gyro * 0.0;
		vecA[2] = gyro * 1.0;
		/* only acceleration */
		/*
		  vecA[0] = gyro * (tracer+ipart)->ax;
		  vecA[1] = gyro * (tracer+ipart)->ay;
		  vecA[2] = gyro * (tracer+ipart)->az;
		*/
   #endif	      
		matS[0][0]+= vecA[0]*vecP[0]; matA[0][1] += vecA[0]*vecP[1]; matA[0][2] += vecA[0]*vecP[2];
		matS[1][0]+= vecA[1]*vecP[0]; matA[1][1] += vecA[1]*vecP[1]; matA[1][2] += vecA[1]*vecP[2];
		matS[2][0]+= vecA[2]*vecP[0]; matA[2][1] += vecA[2]*vecP[1]; matA[2][2] += vecA[2]*vecP[2];	      
	      }/* end if on gyrotaxis_velocity !=0 , to avoid nan */
  #endif /* LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS */

	      /* Now we compute RHS of the Jeffrey equation */

	      /* first product TMP[i] = S[i][j]*P[j] */
	      vecTMP[0]=vecTMP[1]=vecTMP[2]=0;
	      for (i=0; i<3; i++)
		for (j=0; j<3; j++){
		  vecTMP[i] += matS[i][j]*vecP[j];
		}
	      /* then the product of the two vectors P[i]*TMP[i] to form a scalar */
	      scalOSO = 0.0;	      
	      for (i=0; i<3; i++){
		  scalOSO += vecP[i]*vecTMP[i];
		}

	    /* here I add on all the contributions */
	      for (i=0; i<3; i++){
		  vecF[i] = 0.0; 
		for (j=0; j<3; j++){	
		  vecF[i] += matW[i][j]*vecP[j] + ( matS[i][j]*vecP[j]);
		}
		  vecF[i] -=  vecP[i]*scalOSO;
	      }
	     	      
   /* if restart Euler 1st order G = G0 + (DT)*F  */
 if(itime==0 && resume==0){
		  for (i=0; i<3; i++){	
		    vecP[i] =  vecP[i] + property.time_dt*vecF[i];	  
		    /* copy the old term */
		    vecFold[i]  = vecF[i];
		}
 }else{
   /* AB 2nd order G = G0 + (DT/2)*(3*F - Fold)  */
	      for (i=0; i<3; i++){		 
		  vecP[i] =  vecP[i] + 0.5*property.time_dt*( 3.*vecF[i] -  vecFold[i] );
		  /* copy the old term */
		  vecFold[i]  = vecF[i];
		}
 }

 #endif /*LAGRANGE_ORIENTATION_JEFFREY */
 #ifdef LAGRANGE_ORIENTATION_DIFFUSION
	     vec_xi[0] =  random_gauss(0.0,1.0);
	     vec_xi[1] =  random_gauss(0.0,1.0);
	     vec_xi[2] =  random_gauss(0.0,1.0); 	     

	     vecTMP[0]=vecTMP[1]=vecTMP[2]=0;
	      for (i=0; i<3; i++)
		for (j=0; j<3; j++){
		  vecTMP[i] += ( matI[i][j]- vecP[i]*vecP[j] )*vec_xi[j];
		}

	      two_d_r = 2.0*(tracer+ipart)->rotational_diffusion;
	      /* stochastic part of the rotation equation */
	      for (i=0; i<3; i++)
 	      vecP[i] +=  sqrt(two_d_r*property.time_dt)*vecTMP[i];
 #endif /* LAGRANGE_ORIENTATION_DIFFUSION */

 #ifdef LAGRANGE_ORIENTATION_RANDOM
 /* randomly oriented vector */
 vec = random_vector();
#ifdef GRID_POP_D2Q9
 /* generate a random vector in the x,y plane */
 vec = random_vector_2d();
#endif
 vecP[0] = vec.x;
 vecP[1] = vec.y;
 vecP[2] = vec.z;

 /* compute acceleration */
 vecFold[0] = (vecP[0] - (tracer+ipart)->px)/property.time_dt;
 vecFold[1] = (vecP[1] - (tracer+ipart)->py)/property.time_dt;
 vecFold[2] = (vecP[2] - (tracer+ipart)->pz)/property.time_dt;
 #endif /* LAGRANGE_ORIENTATION_RANDOM */

	      /* normalize P vector */
	      norm=0.0;
	      for (i=0; i<3; i++) norm += vecP[i]*vecP[i];
	      for (i=0; i<3; i++) vecP[i]/=sqrt(norm);

	      /* assign P vector */
     (tracer+ipart)->px = vecP[0];
     (tracer+ipart)->py = vecP[1];
     (tracer+ipart)->pz = vecP[2];

              /* assign the just computed dP /dt  vector */     
     (tracer+ipart)->dt_px = vecFold[0];
     (tracer+ipart)->dt_py = vecFold[1];
     (tracer+ipart)->dt_pz = vecFold[2];
     
#endif /* end of lagrange orientation */


   /* In case of BC we use elastic bouncing rule for the particle */
  
#ifdef LB_FLUID_BC
 #ifdef LB_FLUID_BC_Y

   if( (tracer+ipart)->y < 0.0 ){
     (tracer+ipart)->y *= -1.0; 
     //(tracer+ipart)->vx *= -1.0;
     if( (tracer+ipart)->vy < 0.0 ) (tracer+ipart)->vy *= -1.0;
     //(tracer+ipart)->vz *= -1.0;
  #ifdef LAGRANGE_ORIENTATION_BC
     if( (tracer+ipart)->py < 0.0 ) (tracer+ipart)->py *= -1.0;
  #endif 
   }

   if( (tracer+ipart)->y >= property.SY ){
     (tracer+ipart)->y = property.SY- ( (tracer+ipart)->y-property.SY ); 
     //(tracer+ipart)->vx *= -1.0;
     if( (tracer+ipart)->vy > 0.0 ) (tracer+ipart)->vy *= -1.0;
     //(tracer+ipart)->vz *= -1.0;
  #ifdef LAGRANGE_ORIENTATION_BC
     if( (tracer+ipart)->py > 0.0 ) (tracer+ipart)->py *= -1.0;
  #endif 
   }
 #endif

 #ifdef LB_FLUID_BC_X

   if( (tracer+ipart)->x < 0.0 ){
     (tracer+ipart)->x *= -1.0; 
     (tracer+ipart)->vx *= -1.0;
     //(tracer+ipart)->vy *= -1.0;
     //(tracer+ipart)->vz *= -1.0;
   }

   if( (tracer+ipart)->x >= property.SX ){
     (tracer+ipart)->x = property.SX- ( (tracer+ipart)->x-property.SX ); 
     (tracer+ipart)->vx *= -1.0;
     //(tracer+ipart)->vy *= -1.0;
     //(tracer+ipart)->vz *= -1.0;
   }
 #endif

 #ifdef LB_FLUID_BC_Z

   if( (tracer+ipart)->z < 0.0 ){
     (tracer+ipart)->z *= -1.0; 
     //(tracer+ipart)->vx *= -1.0;
     //(tracer+ipart)->vy *= -1.0;
     (tracer+ipart)->vz *= -1.0;
   }

   if( (tracer+ipart)->z >= property.SZ ){
     (tracer+ipart)->z = property.SZ- ( (tracer+ipart)->z-property.SZ ); 
     //(tracer+ipart)->vx *= -1.0;
     //(tracer+ipart)->vy *= -1.0;
     (tracer+ipart)->vz *= -1.0;
   }
 #endif

#endif
   

}/* end of loop on particles */

  sendrecv_particles();

}/* end of move_particles */


/* rearrange particles on the right processors */
void sendrecv_particles(){

  int ipart,i;
  int npart_here,npart_there,all_npart_there,all_npart;  
  point_particle part;
  int *displs,*rcounts;

/* Here we perform the particle re-arrangement between processors */

 npart_here = 0;
 npart_there= 0;
 all_npart_there = 0;
 all_npart = 0;

for (ipart=0;ipart<npart;ipart++) {

  part.x = wrap( (tracer+ipart)->x ,  property.SX);
  part.y = wrap( (tracer+ipart)->y ,  property.SY);
  part.z = wrap( (tracer+ipart)->z ,  property.SZ); 

#ifdef LAGRANGE_WRAP /* this really makes particles to stay in the box */
  (tracer+ipart)->x = wrap( (tracer+ipart)->x ,  property.SX);
  (tracer+ipart)->y = wrap( (tracer+ipart)->y ,  property.SY);
  (tracer+ipart)->z = wrap( (tracer+ipart)->z ,  property.SZ); 
#endif


  //fprintf(stderr,"\n part %g %g %g\n",part.x,part.y,part.z);

  /* check how many particles are still in the local domain (here) and how many have to go (there) */

  /*
 if(  part.x >= center_V[IDX(0, BRD, BRD)].x && part.x < center_V[IDX(LNX+TWO_BRD-1,BRD, BRD)].x &&
      part.y >= center_V[IDX(BRD, 0, BRD)].y && part.y < center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y &&
      part.z >= center_V[IDX(BRD, BRD, 0)].z && part.z < center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z ){
  */
if(  part.x >= mesh[IDXG(BRD, BRD, BRD)].x && part.x < mesh[IDXG(LNXG+BRD-1,BRD, BRD)].x &&
     part.y >= mesh[IDXG(BRD, BRD, BRD)].y && part.y < mesh[IDXG(BRD,LNYG+BRD-1, BRD)].y &&
     part.z >= mesh[IDXG(BRD, BRD, BRD)].z && part.z < mesh[IDXG(BRD, BRD,LNZG+BRD-1)].z ){

      npart_here += 1;

      tracer_here  = (point_particle*) realloc( tracer_here, sizeof(point_particle)*npart_here);
      tracer_here[npart_here-1] = tracer[ipart];

      //fprintf(stderr,"Ehi! ipart %d\n",ipart);
  }else{

      npart_there += 1;   
      tracer_there  = (point_particle*) realloc(tracer_there,sizeof(point_particle)*npart_there);
      tracer_there[npart_there-1] = tracer[ipart];

      //fprintf(stderr,"AAA me %d , tracer_there[%d] %g,  tracer[%d] %g \n",me, npart_there-1, tracer_there[npart_there-1].x , ipart, tracer[ipart].x);

  }/* end of if else */

 }/* for loop on ipart */

//fprintf(stderr,"me %d : npart_here %d , npart_there %d\n",me,npart_here, npart_there);

/* first we communicate to other procs how many particles we have to give away and we sum up them among all the processors */

 MPI_Allreduce(&npart_there, &all_npart_there, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

#ifdef LAGRANGE_DEBUG
 if(ROOT)fprintf(stderr,"me %d : all_npart_there %d\n",me, all_npart_there);
#endif

 if(all_npart_there != 0){

    displs  = (int *)malloc(nprocs*sizeof(int)); 
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart_there, 1 , MPI_INT, rcounts, 1 , MPI_INT,MPI_COMM_WORLD);


    //displs[0]=0;
    //for (i=1;i<nprocs;i++) displs[i] += rcounts[i];
    //for (i=0;i<nprocs;i++) displs[i] = rcounts[i];
    displs[0] = 0;
    for (i=1; i<nprocs; i++) displs[i] = displs[i-1] + rcounts[i-1];


    //for (i=0;i<nprocs;i++) fprintf(stderr,"me %d : rcounts[%d] = %d\n",me,i, rcounts[i]);    
 
  /* space is allocate for coming particles */
  all_tracer_there = (point_particle*) realloc(all_tracer_there,sizeof(point_particle)*all_npart_there);
 
 /* Allgather to get all the migrant particles */
 MPI_Allgatherv(tracer_there, npart_there, MPI_point_particle_type, all_tracer_there, rcounts, displs, MPI_point_particle_type,MPI_COMM_WORLD);
 

 /* Begin loop on particles which just arrived */
 for (ipart=0;ipart<all_npart_there;ipart++) {

   //fprintf(stderr,"BBB me %d , part %g\n",me,(all_tracer_there+ipart)->x);

  part.x = wrap( (all_tracer_there+ipart)->x ,  property.SX);
  part.y = wrap( (all_tracer_there+ipart)->y ,  property.SY);
  part.z = wrap( (all_tracer_there+ipart)->z ,  property.SZ); 

  /* check how many particles are still in the local domain (here) and how many have to go (there) */
  /*
  if( part.x >= center_V[IDX(0, BRD, BRD)].x && part.x < center_V[IDX(LNX+TWO_BRD-1, BRD, BRD)].x &&
      part.y >= center_V[IDX(BRD, 0, BRD)].y && part.y < center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y &&
      part.z >= center_V[IDX(BRD, BRD, 0)].z && part.z < center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z ){
  */
  if(  part.x >= mesh[IDXG(BRD, BRD, BRD)].x && part.x < mesh[IDXG(LNXG+BRD-1,BRD, BRD)].x &&
       part.y >= mesh[IDXG(BRD, BRD, BRD)].y && part.y < mesh[IDXG(BRD,LNYG+BRD-1, BRD)].y &&
       part.z >= mesh[IDXG(BRD, BRD, BRD)].z && part.z < mesh[IDXG(BRD, BRD,LNZG+BRD-1)].z ){

      npart_here += 1;

      tracer_here  = (point_particle*) realloc(tracer_here,sizeof(point_particle)*npart_here);
      tracer_here[npart_here-1] = all_tracer_there[ipart];
  }/* end if */

 }/* for on ipart till all_npart_there */


 /* here I realloc tracer and copy the new data into it */

      tracer = (point_particle*) realloc(tracer,sizeof(point_particle)*npart_here);   

      for (ipart=0;ipart<npart_here;ipart++){  tracer[ipart] = tracer_here[ipart];}

       free(rcounts);
       free(displs);
 }/* end on if on all_npart_there */

      npart = npart_here;

      //fprintf(stderr,"me %d : new npart is %d\n",me, npart);


      /* for debug */
      /*
      fprintf(stderr," %g %g\n %g %g\n %g %g\n", center_V[IDX(0, BRD, BRD)].x , center_V[IDX(LNX+TWO_BRD-1,BRD, BRD)].x, 
	                                         center_V[IDX(BRD, 0 , BRD)].y , center_V[IDX(BRD,LNY+TWO_BRD-1, BRD)].y,
	                                         center_V[IDX(BRD, BRD, 0)].z , center_V[IDX(BRD, BRD,LNZ+TWO_BRD-1)].z );
      */
      /*
     fprintf(stderr," %g %g\n %g %g\n %g %g\n", mesh[IDXG(BRD, BRD, BRD)].x , mesh[IDXG(LNXG+BRD-1,BRD, BRD)].x, 
	                                         mesh[IDXG(BRD, BRD , BRD)].y , mesh[IDXG(BRD,LNYG+BRD-1, BRD)].y,
	                                         mesh[IDXG(BRD, BRD, BRD)].z , mesh[IDXG(BRD, BRD,LNZG+BRD-1)].z );
      */
      /* final check */
       MPI_Allreduce(&npart, &all_npart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

#ifdef LAGRANGE_DEBUG
       if(itime%10==0 && ROOT)fprintf(stderr,"------------------ Check npart = %d\n",all_npart);
#endif
       if(all_npart != (int)property.particle_number){
         if(ROOT) fprintf(stderr,"Total number of particles has changed during run!!!!\n Was %d now is %d\n Exit.\n",(int)property.particle_number, all_npart);
         MPI_Finalize();
         exit(0);
       }

 

}/* end of sendrecv particles */



/* here we have the function to write the full point_particle structure and to read it */

/* general output function for particles */

#define H5FILE_NAME_PART "part.h5"

#ifdef OUTPUT_H5

void write_point_particle_h5(){
  int i,j;
  int np = (int)property.particle_number;
  FILE *fout;

    int *rcounts;
    int name_offset = 0;

    /* First check how many particles in each processor and compute offset */
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);


    hid_t       file_id, dataset_id, dataspace_id , group ;  /* identifiers */
    hid_t	plist_id;                 /* property list identifier */
    hid_t hdf5_type;
    hid_t xfer_plist, ret, property_id;
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dims[1], offset[1], count[1];
    herr_t      hdf5_status;
    herr_t status;
    int size;    
    int RANK = 1;

    my_double *aux;

    char NEW_H5FILE_NAME[128];
    char XMF_FILE_NAME[128];
 
    char label[128]; 

    /* create point particle compound */
    hdf5_type = H5Tcreate (H5T_COMPOUND, sizeof(point_particle));

    /* define its offsets */
    sprintf(label,"name");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, name), H5T_NATIVE_DOUBLE);

    sprintf(label,"x");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, x), H5T_NATIVE_DOUBLE);
    sprintf(label,"y");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, y), H5T_NATIVE_DOUBLE);
    sprintf(label,"z");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, z), H5T_NATIVE_DOUBLE);

    sprintf(label,"vx");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vx), H5T_NATIVE_DOUBLE);
    sprintf(label,"vy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vy), H5T_NATIVE_DOUBLE);
    sprintf(label,"vz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vz), H5T_NATIVE_DOUBLE);

    sprintf(label,"vx_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vx_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"vy_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vy_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"vz_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vz_old), H5T_NATIVE_DOUBLE);

    sprintf(label,"ax");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ax), H5T_NATIVE_DOUBLE);
    sprintf(label,"ay");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ay), H5T_NATIVE_DOUBLE);
    sprintf(label,"az");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, az), H5T_NATIVE_DOUBLE);

    sprintf(label,"ax_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ax_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"ay_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ay_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"az_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, az_old), H5T_NATIVE_DOUBLE);

    sprintf(label,"tau_drag");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, tau_drag), H5T_NATIVE_DOUBLE);

#ifdef LB_FLUID
    sprintf(label,"ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uz), H5T_NATIVE_DOUBLE);

    sprintf(label,"ux_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ux_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"uy_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uy_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"uz_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uz_old), H5T_NATIVE_DOUBLE);
    //#endif
 #ifdef LAGRANGE_GRAVITY
  #ifdef LAGRANGE_GRAVITY_VARIABLE
     sprintf(label,"gravity_coeff");
     H5Tinsert(hdf5_type, label, HOFFSET(point_particle, gravity_coeff), H5T_NATIVE_DOUBLE);
  #endif
 #endif
 #ifdef LAGRANGE_GRADIENT
    sprintf(label,"dx_ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"dx_uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"dx_uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_uz), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_uz), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_uz), H5T_NATIVE_DOUBLE);
  #ifdef LAGRANGE_ADDEDMASS
     sprintf(label,"beta_coeff");
     H5Tinsert(hdf5_type, label, HOFFSET(point_particle, beta_coeff), H5T_NATIVE_DOUBLE);
  #endif
  #ifdef LAGRANGE_ORIENTATION
    sprintf(label,"px");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, px), H5T_NATIVE_DOUBLE);
    sprintf(label,"py");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, py), H5T_NATIVE_DOUBLE);
    sprintf(label,"pz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, pz), H5T_NATIVE_DOUBLE);
    sprintf(label,"dt_px");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dt_px), H5T_NATIVE_DOUBLE);
    sprintf(label,"dt_py");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dt_py), H5T_NATIVE_DOUBLE);
    sprintf(label,"dt_pz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dt_pz), H5T_NATIVE_DOUBLE);
   #ifdef LAGRANGE_ORIENTATION_JEFFREY
    sprintf(label,"aspect_ratio");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, aspect_ratio), H5T_NATIVE_DOUBLE);
    #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
    sprintf(label,"gyrotaxis_velocity");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, gyrotaxis_velocity), H5T_NATIVE_DOUBLE);
    #endif
   #endif  
   #ifdef LAGRANGE_ORIENTATION_DIFFUSION
    sprintf(label,"rotational_diffusion");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, rotational_diffusion), H5T_NATIVE_DOUBLE);
   #endif 
   #ifdef LAGRANGE_ORIENTATION_ACTIVE
    sprintf(label,"swim_velocity");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, swim_velocity), H5T_NATIVE_DOUBLE);
    #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
    sprintf(label,"critical_shear_rate");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, critical_shear_rate), H5T_NATIVE_DOUBLE);
    sprintf(label,"shear_rate");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, shear_rate), H5T_NATIVE_DOUBLE);
    sprintf(label,"jump_time");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, jump_time), H5T_NATIVE_DOUBLE);
    /* jump orientation */
    sprintf(label,"px_jump");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, px_jump), H5T_NATIVE_DOUBLE);
    sprintf(label,"py_jump");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, py_jump), H5T_NATIVE_DOUBLE);
    sprintf(label,"pz_jump");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, pz_jump), H5T_NATIVE_DOUBLE);
    #endif
   #endif 
  #endif /* LAGRANGE_ORIENTATION */
 #endif /* LAGRANGE_GRADIENT */
#endif /* LB_FLUID */

#ifdef LB_TEMPERATURE
    sprintf(label,"t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, t), H5T_NATIVE_DOUBLE);
 #ifdef LAGRANGE_GRADIENT
    sprintf(label,"dx_t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_t), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_t), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_t), H5T_NATIVE_DOUBLE);
 #endif
#endif

#ifdef LB_SCALAR
    sprintf(label,"s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, s), H5T_NATIVE_DOUBLE);
 #ifdef LAGRANGE_GRADIENT
    sprintf(label,"dx_s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_s), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_s), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_s), H5T_NATIVE_DOUBLE);
 #endif
#endif

    /*************************************************************/
  
                /* Create a new file using default properties */
		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);   
		
		file_id = H5Fcreate(H5FILE_NAME_PART, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		group   = H5Gcreate (file_id, "/particles", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

		H5Pclose(plist_id);

		property_id  = H5Pcreate(H5P_DATASET_CREATE);


                /* Create the data space for the dataset. */
                dims[0] = (int)property.particle_number;

		filespace = H5Screate_simple(RANK, dims, NULL);   
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
		count[0] = npart;
		offset[0] = name_offset;

	        memspace = H5Screate_simple(RANK, count, NULL); 	

    /*
     * Select hyperslab in the file.
     */
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

		xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);                       

		/* WRITE POINT_PARTICLE STRUCTURE */
		dataset_id = H5Dcreate(group, "point_particle", hdf5_type, filespace,H5P_DEFAULT, H5P_DEFAULT ,H5P_DEFAULT);
                ret = H5Dwrite(dataset_id, hdf5_type, memspace, filespace, xfer_plist, tracer);
                status = H5Dclose(dataset_id);                

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

 /* create the file names */
 // sprintf(NEW_H5FILE_NAME,"part.h5");

  /* we rename the file */
  //if(ROOT) rename(H5FILE_NAME_PARTICLE, NEW_H5FILE_NAME);

}/* end of write point particle */


/* Now the function to read point particles */
void read_point_particle_h5(){
  int i,j;
  int np = (int)property.particle_number;
  FILE *fout;

    int *rcounts;
    int name_offset = 0;

    hid_t       file_id, dataset_id, dataspace_id , group ;  /* identifiers */
    hid_t	plist_id;                 /* property list identifier */
    hid_t hdf5_type;
    hid_t xfer_plist, ret, property_id;
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dims[1], offset[1], count[1];
    herr_t      hdf5_status;
    herr_t status;
    int size;    
    int RANK = 1;

    my_double *aux;

    char NEW_H5FILE_NAME[128];
    char XMF_FILE_NAME[128];
 
    char label[128]; 


    /* First check how many particles in each processor and compute offset */
    rcounts = (int *)malloc(nprocs*sizeof(int)); 

    MPI_Allgather(&npart, 1 , MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);

    for (i=0;i<me;i++) name_offset += rcounts[i];

    free(rcounts);

    /* create point particle compound */
    hdf5_type = H5Tcreate (H5T_COMPOUND, sizeof(point_particle));

    /* define its offsets */
    sprintf(label,"name");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, name), H5T_NATIVE_DOUBLE);

    sprintf(label,"x");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, x), H5T_NATIVE_DOUBLE);
    sprintf(label,"y");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, y), H5T_NATIVE_DOUBLE);
    sprintf(label,"z");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, z), H5T_NATIVE_DOUBLE);

    sprintf(label,"vx");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vx), H5T_NATIVE_DOUBLE);
    sprintf(label,"vy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vy), H5T_NATIVE_DOUBLE);
    sprintf(label,"vz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vz), H5T_NATIVE_DOUBLE);

    sprintf(label,"vx_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vx_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"vy_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vy_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"vz_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, vz_old), H5T_NATIVE_DOUBLE);

    sprintf(label,"ax");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ax), H5T_NATIVE_DOUBLE);
    sprintf(label,"ay");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ay), H5T_NATIVE_DOUBLE);
    sprintf(label,"az");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, az), H5T_NATIVE_DOUBLE);

    sprintf(label,"ax_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ax_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"ay_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ay_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"az_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, az_old), H5T_NATIVE_DOUBLE);

    sprintf(label,"tau_drag");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, tau_drag), H5T_NATIVE_DOUBLE);

#ifdef LB_FLUID
    sprintf(label,"ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uz), H5T_NATIVE_DOUBLE);

    sprintf(label,"ux_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, ux_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"uy_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uy_old), H5T_NATIVE_DOUBLE);
    sprintf(label,"uz_old");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, uz_old), H5T_NATIVE_DOUBLE);
    //#endif

 #ifdef LAGRANGE_GRADIENT
    sprintf(label,"dx_ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_ux");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_ux), H5T_NATIVE_DOUBLE);
    sprintf(label,"dx_uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_uy");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_uy), H5T_NATIVE_DOUBLE);
    sprintf(label,"dx_uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_uz), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_uz), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_uz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_uz), H5T_NATIVE_DOUBLE);
  #ifdef LAGRANGE_ADDEDMASS
     sprintf(label,"beta_coeff");
     H5Tinsert(hdf5_type, label, HOFFSET(point_particle, beta_coeff), H5T_NATIVE_DOUBLE);
  #endif
  #ifdef LAGRANGE_ORIENTATION
    sprintf(label,"px");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, px), H5T_NATIVE_DOUBLE);
    sprintf(label,"py");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, py), H5T_NATIVE_DOUBLE);
    sprintf(label,"pz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, pz), H5T_NATIVE_DOUBLE);
    sprintf(label,"dt_px");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dt_px), H5T_NATIVE_DOUBLE);
    sprintf(label,"dt_py");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dt_py), H5T_NATIVE_DOUBLE);
    sprintf(label,"dt_pz");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dt_pz), H5T_NATIVE_DOUBLE);
   #ifdef LAGRANGE_ORIENTATION_JEFFREY
    sprintf(label,"aspect_ratio");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, aspect_ratio), H5T_NATIVE_DOUBLE);
    #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
    sprintf(label,"gyrotaxis_velocity");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, gyrotaxis_velocity), H5T_NATIVE_DOUBLE); 
    #endif
   #endif
   #ifdef LAGRANGE_ORIENTATION_DIFFUSION
    sprintf(label,"rotational_diffusion");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, rotational_diffusion), H5T_NATIVE_DOUBLE); 
   #endif
   #ifdef LAGRANGE_ORIENTATION_ACTIVE
    sprintf(label,"swim_velocity");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, swim_velocity), H5T_NATIVE_DOUBLE);
    #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
    sprintf(label,"critical_shear_rate");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, critical_shear_rate), H5T_NATIVE_DOUBLE);
    sprintf(label,"shear_rate");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, shear_rate), H5T_NATIVE_DOUBLE);
    sprintf(label,"jump_time");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, jump_time), H5T_NATIVE_DOUBLE);
    /* jump orientation */
    sprintf(label,"px_jump");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, px_jump), H5T_NATIVE_DOUBLE);
    sprintf(label,"py_jump");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, py_jump), H5T_NATIVE_DOUBLE);
    sprintf(label,"pz_jump");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, pz_jump), H5T_NATIVE_DOUBLE);
    #endif
   #endif
  #endif
 #endif
#endif

#ifdef LB_TEMPERATURE
    sprintf(label,"t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, t), H5T_NATIVE_DOUBLE);
#ifdef LAGRANGE_GRADIENT
    sprintf(label,"dx_t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_t), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_t), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_t");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_t), H5T_NATIVE_DOUBLE);
#endif
#endif

#ifdef LB_SCALAR
    sprintf(label,"s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, s), H5T_NATIVE_DOUBLE);
#ifdef LAGRANGE_GRADIENT
    sprintf(label,"dx_s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dx_s), H5T_NATIVE_DOUBLE);
    sprintf(label,"dy_s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dy_s), H5T_NATIVE_DOUBLE);
    sprintf(label,"dz_s");
    H5Tinsert(hdf5_type, label, HOFFSET(point_particle, dz_s), H5T_NATIVE_DOUBLE);
#endif
#endif

    /*************************************************************/
  
                /* Create a new file using default properties */
		plist_id = H5Pcreate(H5P_FILE_ACCESS);
		hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);   
		
		file_id = H5Fopen(H5FILE_NAME_PART, H5F_ACC_RDONLY, H5P_DEFAULT); 
		group   = H5Gopen(file_id, "/particles", H5P_DEFAULT);

		H5Pclose(plist_id);

		property_id  = H5Pcreate(H5P_DATASET_CREATE);


                /* Create the data space for the dataset. */
                dims[0] = (int)property.particle_number;

		filespace = H5Screate_simple(RANK, dims, NULL);   
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
		count[0] = npart;
		offset[0] = name_offset;

	        memspace = H5Screate_simple(RANK, count, NULL); 	

    /*
     * Select hyperslab in the file.
     */
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

		xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);                       

		/* READ POINT_PARTICLE STRUCTURE */
		dataset_id = H5Dopen(group, "point_particle", H5P_DEFAULT);
                ret = H5Dread(dataset_id, hdf5_type, memspace, filespace,  H5P_DEFAULT, tracer);
                status = H5Dclose(dataset_id);                

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

  /* here we rearrange particles */
  sendrecv_particles();

}/* end of read point particle */


#endif


#endif
