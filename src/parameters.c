#include "common_object.h"

my_double read_parameter(char * variable){

  char fnamein[256], fnameout[256];
  char name[256] = "NULL";
  double val;
  FILE *fin, *fout;
  int i;
  int cmp=1;

  sprintf(fnamein,"param.in");
  sprintf(fnameout,"param.out");
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

  if(ROOT){
    sprintf(OutDir,"RUN");
    mkdir(OutDir, S_IWUSR|S_IXUSR|S_IRUSR);
    fprintf(stderr,"OutDir is %s\n",OutDir);   

    remove("param.out");
    /* read parameters from file */


    /* resume */
    sprintf(name,"resume");
    resume = (int)read_parameter(name);

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

    /* time stepping parameters */
    sprintf(name,"time_dt");
    property.time_dt = (double)read_parameter(name); 
    fprintf(stderr,"time step: %g\n",property.time_dt);

    sprintf(name,"time_start");
    property.time_start = (double)read_parameter(name); 
    fprintf(stderr,"Time start: %g\n",property.time_start);

    sprintf(name,"time_end");
    property.time_end = (double)read_parameter(name); 
    fprintf(stderr,"Time end: %g\n",property.time_end);
  
    sprintf(name,"time_dump_field");
    property.time_dump_field = (double)read_parameter(name);
    fprintf(stderr,"Time dump fields: %g\n",property.time_dump_field);
   
    sprintf(name,"time_dump_diagn");
    property.time_dump_diagn = (double)read_parameter(name);
    fprintf(stderr,"Time dump fields: %g\n",property.time_dump_diagn);


#ifdef LB_FLUID
    /* relaxation time and viscosity for fluid */
  fprintf(stderr,"YES <- LB_FLUID\n");
  sprintf(name,"tau_u");
  property.tau_u = read_parameter(name);
  fprintf(stderr,"Properties:\ntau_u %g\n",(double)property.tau_u);
  property.nu = property.tau_u/3.0;
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
  fprintf(stderr,"Properties:\Amp_x %g\n",(double)property.Amp_x);
  fprintf(stderr,"Properties:\Amp_y %g\n",(double)property.Amp_y);
  fprintf(stderr,"Properties:\Amp_z %g\n",(double)property.Amp_z);
#endif
#endif

#ifdef LB_TEMPERATURE
    /* relaxation time and viscosity for temperature */
  fprintf(stderr,"YES <- LB_TEMPERATURE\n");
  sprintf(name,"tau_t");
  property.tau_t = read_parameter(name);
  fprintf(stderr,"Properties:\ntau_t %g\n",(double)property.tau_t);
  property.kappa = property.tau_t/3.0;
  fprintf(stderr,"thermal diffusivity %g\n",(double)property.kappa);
  sprintf(name,"T_bot");
  property.T_bot = read_parameter(name);
  sprintf(name,"T_top");
  property.T_top = read_parameter(name);
  sprintf(name,"T_ref");
  property.T_ref = read_parameter(name);
  property.deltaT = property.T_bot-property.T_top;
  fprintf(stderr,"T_bot %g , T_top %g , deltaT %g\n",(double)property.T_bot, (double)property.T_top, (double)property.deltaT);
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

  fprintf(stderr,"And the Rayleigh Number is -> Ra = %e\n", property.beta_t*property.gravity_y*property.deltaT*pow(property.SY,3.0)/(property.nu*property.kappa) );

#endif 
#endif

#ifdef LB_SCALAR
    /* relaxation time and viscosity for temperature */
  fprintf(stderr,"YES <- LB_SCALAR\n");
  sprintf(name,"tau_s");
  property.tau_s = read_parameter(name);
  fprintf(stderr,"Properties:\ntau_s %g\n",(double)property.tau_s);
  property.chi = property.tau_s/3.0;
  fprintf(stderr,"mass diffusivity %g\n",(double)property.chi);
#endif


  /* size of types, just for a check */
    fprintf(stderr,"Size of float %d\n",sizeof(float));
    fprintf(stderr,"Size of double %d\n",sizeof(double));
    fprintf(stderr,"Size of long double %d\n",sizeof(long double));
    fprintf(stderr,"Size of my_double %d\n",sizeof(my_double));
    fprintf(stderr,"\n");
  }/* if ROOT*/

 /* Now broadcast all properties */
 MPI_Bcast(&resume, 1, MPI_INT, 0, MPI_COMM_WORLD);

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

}


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


void allocate_fields(){
 mesh  = (vector*) malloc(sizeof(vector)*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 if(mesh == NULL){ fprintf(stderr,"Not enough memory to allocate mesh\n"); exit(-1);}
 set_to_zero_vector( mesh,(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD));

 mesh_flag  = (int*) malloc(sizeof(int)*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 if(mesh_flag == NULL){ fprintf(stderr,"Not enough memory to allocate mesh_flag\n"); exit(-1);}
 set_to_zero_int( mesh_flag,(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD));

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

 
 //#ifdef METHOD_CENTERED
#if (defined METHOD_CENTERED || defined METHOD_MYQUICK)
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
 
#ifdef METHOD_MYQUICK
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



#ifdef LB_FLUID
 p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(p == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 set_to_zero_pop( p,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 rhs_p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(rhs_p == NULL){ fprintf(stderr,"Not enough memory to allocate rhs_p\n"); exit(-1);}
 set_to_zero_pop( rhs_p,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 old_rhs_p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_rhs_p == NULL){ fprintf(stderr,"Not enough memory to allocate old_rhs_p\n"); exit(-1);}
 set_to_zero_pop( old_rhs_p,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 u = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(u == NULL){ fprintf(stderr,"Not enough memory to allocate u\n"); exit(-1);}
 set_to_zero_vector( u,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 dens  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(dens == NULL){ fprintf(stderr,"Not enough memory to allocate dens\n"); exit(-1);}
 set_to_zero_my_double( dens,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));


#ifdef LB_FLUID_FORCING
 force  = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(force == NULL){ fprintf(stderr,"Not enough memory to allocate force\n"); exit(-1);}
 set_to_zero_vector( force,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
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

#endif

#ifdef LB_FLUID
 ruler_x_local  = (output*) malloc(sizeof(output)*NX);
 ruler_y_local  = (output*) malloc(sizeof(output)*NY);
 ruler_z_local  = (output*) malloc(sizeof(output)*NZ);
 ruler_x  = (output*) malloc(sizeof(output)*NX);
 ruler_y  = (output*) malloc(sizeof(output)*NY);
 ruler_z  = (output*) malloc(sizeof(output)*NZ);
#endif

#ifdef LB_TEMPERATURE
 g = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(g == NULL){ fprintf(stderr,"Not enough memory to allocate g\n"); exit(-1);}
 set_to_zero_pop( g,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 rhs_g  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(rhs_g == NULL){ fprintf(stderr,"Not enough memory to allocate rhs_g\n"); exit(-1);}
 set_to_zero_pop( rhs_g,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 old_rhs_g  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_rhs_g == NULL){ fprintf(stderr,"Not enough memory to allocate old_rhs_g\n"); exit(-1);}
 set_to_zero_pop( old_rhs_g,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 t  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(t == NULL){ fprintf(stderr,"Not enough memory to allocate t\n"); exit(-1);}
 set_to_zero_my_double( t,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

#ifdef LB_TEMPERATURE_FORCING
 t_source  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(t_source == NULL){ fprintf(stderr,"Not enough memory to allocate source_t\n"); exit(-1);}
 set_to_zero_scalar( t_source,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif
#endif

#ifdef LB_SCALAR
 h = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(h == NULL){ fprintf(stderr,"Not enough memory to allocate h\n"); exit(-1);}
 set_to_zero_pop( h,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 rhs_h  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(rhs_h == NULL){ fprintf(stderr,"Not enough memory to allocate rhs_h\n"); exit(-1);}
 set_to_zero_pop( rhs_h,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 old_rhs_h  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(old_rhs_h == NULL){ fprintf(stderr,"Not enough memory to allocate old_rhs_h\n"); exit(-1);}
 set_to_zero_pop( old_rhs_h,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));

 s  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(s == NULL){ fprintf(stderr,"Not enough memory to allocate s\n"); exit(-1);}
 set_to_zero_my_double( s,(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
#endif

}
