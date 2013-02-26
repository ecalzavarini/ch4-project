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

    sprintf(name,"time_max");
    property.time_max = (double)read_parameter(name); 
    fprintf(stderr,"Total time steps: %g\n",property.time_max);
  
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
#ifdef LB_FLUID_FORCING_POISEUILLE
  /* pressure gradient */
  fprintf(stderr,"YES <- LB_FLUID_FORCING_POISEUILLE\n");
  sprintf(name,"gradP");
  property.gradP = read_parameter(name);
  fprintf(stderr,"Properties:\ngradP %g\n",(double)property.gradP);
#endif
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
 MPI_Bcast(&property, 1,MPI_Property_type, 0, MPI_COMM_WORLD);
 

#ifdef DEBUG
 for(i=0;i<nprocs;i++){
   if(i==me){ fprintf(stderr,"me %d , System Size: NX %d NY %d NZ %d\n", me , NX , NY, NZ);
             fprintf(stderr,"me %d , time_dt %g\n", me , property.time_dt);
   }
 }
#endif

}



void allocate_fields(){
 mesh  = (vector*) malloc(sizeof(vector)*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 if(mesh == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

 mesh_flag  = (int*) malloc(sizeof(int)*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 if(mesh_flag == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

 center_V = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(center_V == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

 /* borders mesh */
 xp_mesh  = (vector*) malloc(sizeof(vector)*BRD*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 xm_mesh  = (vector*) malloc(sizeof(vector)*BRD*(LNYG+TWO_BRD)*(LNZG+TWO_BRD));
 xp_flag  = (int*) malloc(sizeof(int)*BRD*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 xm_flag  = (int*) malloc(sizeof(int)*BRD*(LNYG+TWO_BRD)*(LNZG+TWO_BRD)); 
 if(xp_mesh == NULL || xm_mesh == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 if(xp_flag == NULL || xm_flag == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 yp_mesh  = (vector*) malloc(sizeof(vector)*BRD*(LNXG+TWO_BRD)*(LNZG+TWO_BRD)); 
 ym_mesh  = (vector*) malloc(sizeof(vector)*BRD*(LNXG+TWO_BRD)*(LNZG+TWO_BRD)); 
 yp_flag  = (int*) malloc(sizeof(int)*BRD*(LNXG+TWO_BRD)*(LNZG+TWO_BRD)); 
 ym_flag  = (int*) malloc(sizeof(int)*BRD*(LNXG+TWO_BRD)*(LNZG+TWO_BRD));
 if(yp_mesh == NULL || ym_mesh == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 if(yp_flag == NULL || ym_flag == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 zp_mesh  = (vector*) malloc(sizeof(vector)*BRD*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)); 
 zm_mesh  = (vector*) malloc(sizeof(vector)*BRD*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)); 
 zp_flag  = (int*) malloc(sizeof(int)*BRD*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)); 
 zm_flag  = (int*) malloc(sizeof(int)*BRD*(LNXG+TWO_BRD)*(LNYG+TWO_BRD)); 
 if(zp_mesh == NULL || zm_mesh == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 if(zp_flag == NULL || zm_flag == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}


#ifdef LB
 coeff_xp = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_xp == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_xm = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_xm == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_yp = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_yp == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_ym = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_ym == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_zp = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_zp == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_zm = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(coeff_zm == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
#endif


#ifdef LB_FLUID
 p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(p == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

 rhs_p  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(rhs_p == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

 u = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(u == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

 dens  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(dens == NULL){ fprintf(stderr,"Not enough memory to allocate dens\n"); exit(-1);}

#ifdef LB_FLUID_FORCING
 force  = (vector*) malloc(sizeof(vector)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(force == NULL){ fprintf(stderr,"Not enough memory to allocate force\n"); exit(-1);}
#endif


 /* borders pop */
 xp_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
 xm_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
 if(xp_pop == NULL || xm_pop == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 yp_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
 ym_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
 if(yp_pop == NULL || ym_mesh == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 zp_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
 zm_pop  = (pop*) malloc(sizeof(pop)*BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
 if(zp_pop == NULL || zm_pop == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

#endif


#ifdef LB_FLUID
 ruler_x  = (vector*) malloc(sizeof(vector)*NX);
 ruler_y  = (vector*) malloc(sizeof(vector)*NY);
 ruler_z  = (vector*) malloc(sizeof(vector)*NZ);
#endif

}
