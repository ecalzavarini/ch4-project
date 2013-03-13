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
  /* pressure gradient */
  fprintf(stderr,"YES <- LB_FLUID_FORCING\n");
  sprintf(name,"gradP");
  property.gradP = read_parameter(name);
  fprintf(stderr,"Properties:\ngradP %g\n",(double)property.gradP);
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
