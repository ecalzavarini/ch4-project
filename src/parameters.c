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
    sprintf(name,"NX");
    property.NX = NX = (int)read_parameter(name);
    sprintf(name,"NY");
    property.NY = NY = (int)read_parameter(name);
    sprintf(name,"NZ");
    property.NZ = NZ = (int)read_parameter(name);
    fprintf(stderr,"System Size:\nNX %d \nNY %d \nNZ %d\n", NX , NY, NZ);

    /*
      sprintf(name,"max_step");
      max_step = (int)read_parameter(name); 
      fprintf(stderr,"Total time steps: %d\n",max_step);
  
      sprintf(name,"time_dump_field");
      time_dump_field = (int)read_parameter(name);
      fprintf(stderr,"Time dump fields: %d\n",time_dump_field);
    */

    fprintf(stderr,"Size of float %d\n",sizeof(float));
    fprintf(stderr,"Size of double %d\n",sizeof(double));
    fprintf(stderr,"Size of long double %d\n",sizeof(long double));
    fprintf(stderr,"Size of my_double %d\n",sizeof(my_double));
    fprintf(stderr,"\n");
  }/* if ROOT*/

 /* Now broadcast */
 MPI_Bcast(&NX, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&NY, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef DEBUG
 for(i=0;i<nprocs;i++){
   if(i==me) fprintf(stderr,"me %d , System Size: NX %d NY %d NZ %d\n", me , NX , NY, NZ);
 }
#endif

}



void allocate_fields(){
 mesh  = (mesh_type*) malloc(sizeof(mesh_type)*(LNX+BX)*(LNY+BY)*(LNZ+BZ)); 
 if(mesh == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

#ifdef LB
 coeff_xp = (pop*) malloc(sizeof(pop)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(coeff_xp == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_xm = (pop*) malloc(sizeof(pop)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(coeff_xm == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_yp = (pop*) malloc(sizeof(pop)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(coeff_yp == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_ym = (pop*) malloc(sizeof(pop)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(coeff_ym == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_zp = (pop*) malloc(sizeof(pop)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(coeff_zp == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
 coeff_zm = (pop*) malloc(sizeof(pop)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(coeff_zm == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
#endif

 center_V = (vector*) malloc(sizeof(vector)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(center_V == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

#ifdef LB_FLUID
 p  = (pop*) malloc(sizeof(pop)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(p == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}

 u = (vector*) malloc(sizeof(vector)*(LNX+BX-1)*(LNY+BY-1)*(LNZ+BZ-1)); 
 if(u == NULL){ fprintf(stderr,"Not enough memory to allocate p\n"); exit(-1);}
#endif

}
