#include "common_object.h"

void dump_averages(int itime){
  int i,j,k;
  FILE *fout;
  char fname[128];
  my_double ux2,uy2,uz2;


#ifdef LB_FLUID
    out_local.ux = out_local.uy = out_local.uz = 0.0; 
    out_local.ux2 = out_local.uy2 = out_local.uz2 = 0.0; 
    out_local.rho = 0.0;
    out_local.ene = 0.0;
    out_local.eps = 0.0;
#endif


	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {

#ifdef LB_FLUID

			  out_local.ux+=u[IDX(i, j, k)].x;
			  out_local.uy+=u[IDX(i, j, k)].y;
			  out_local.uz+=u[IDX(i, j, k)].z;

			  ux2=u[IDX(i, j, k)].x*u[IDX(i, j, k)].x;
			  out_local.ux2+=ux2;
			    
			  uy2=u[IDX(i, j, k)].y*u[IDX(i, j, k)].y;
			  out_local.uy2+=uy2;

			  uz2=u[IDX(i, j, k)].z*u[IDX(i, j, k)].z;
			  out_local.uz2+=uz2;

			  out_local.rho += dens[IDX(i, j, k)];
			  out_local.ene += 0.5*(ux2+uy2+uz2);
			  
			  out_local.eps = 0.0;

#endif
			} /* for i j k */      
  
#ifdef LB_FLUID
  MPI_Allreduce(&out_local, &out_all, 1, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );

  if(ROOT){
    sprintf(fname,"velocity_averages.dat");
    fout = fopen(fname,"a");
    fprintf(fout,"%d %e %e %e %e %e %e %e %e\n",itime, (double)out_all.ene, (double)out_all.rho, (double)out_all.ux, (double)out_all.uy, (double)out_all.uz, (double)out_all.ux2 , (double)out_all.uy2, (double)out_all.uz2);
    fclose(fout);
  }
#endif
  

}

