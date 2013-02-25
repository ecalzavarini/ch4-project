#include "common_object.h"

void dump_averages(int itime){
  int i,j,k;
  FILE *fout;
  char fname[128];
  my_double ux2,uy2,uz2;

 


#ifdef LB_FLUID
    out.ux = out.uy = out.uz = 0.0; 
    out.ux2 = out.uy2 = out.uz2 = 0.0; 
    out.rho = 0.0;
    out.ene = 0.0;
    out.eps = 0.0;
#endif


	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {

#ifdef LB_FLUID

			  out.ux+=u[IDX(i, j, k)].x;
			  out.uy+=u[IDX(i, j, k)].y;
			  out.uz+=u[IDX(i, j, k)].z;

			  ux2=u[IDX(i, j, k)].x*u[IDX(i, j, k)].x;
			  out.ux2+=ux2;
			    
			  uy2=u[IDX(i, j, k)].y*u[IDX(i, j, k)].y;
			  out.uy2+=uy2;

			  uz2=u[IDX(i, j, k)].z*u[IDX(i, j, k)].z;
			  out.uz2+=uz2;

			  out.rho += dens[IDX(i, j, k)];
			  out.ene += 0.5*(ux2+uy2+uz2);
			  
			  out.eps = 0.0;

#endif
			} /* for i j k */      
  

#ifdef LB_FLUID
    sprintf(fname,"velocity.dat");
    fout = fopen(fname,"a");
    fprintf(fout,"%d %e %e %e %e %e %e %e %e\n",itime, (double)out.ene, (double)out.rho, (double)out.ux, (double)out.uy, (double)out.uz, (double)out.ux2 , (double)out.uy2, (double)out.uz2);
    fclose(fout);

#endif
  

}

