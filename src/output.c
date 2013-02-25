#include "common_object.h"

void dump_averages(int tstep){
  FILE *fout;
  char fname[128];
  my_double vx2,vy2,vz2;

 


#ifdef LB_FLUID
    out.ux = rbout.uy = out.uz = 0.0; 
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

			  out.ux2+=ux2=u[IDX(i, j, k)].x*v[IDX(i, j, k)].x;
			  out.uy2+=uy2=u[IDX(i, j, k)].y*v[IDX(i, j, k)].y;
			  out.uz2+=uz2=u[IDX(i, j, k)].z*v[IDX(i, j, k)].z;

			  out.rho += dens[IDX(i, j, k)];
			  out.ene += 0.5*(ux2+uy2+uz2);
			  
			  out.eps = 0.0;

#endif
			} /* for i j k */      
  

#ifdef LB_FLUID
    sprintf(fname,"velocity.dat");
    fout = fopen(fname,"a");
    fprintf(fout,"%d %e %e %e %e %e\n",tstep, (double)out.ux, (double)out.uy,(double)rbout.rho, (double)rbout.ux2 , (double)rbout.uy2);
    fclose(fout);

    sprintf(fname,"velocity_y.dat");
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) fprintf(fout,"%d %e %e %e %e %e\n",y, (double)vx_y[y-1], (double)vy_y[y-1], (double)rho_y[y-1],  (double)vx2_y[y-1], (double)vy2_y[y-1]);
    fclose(fout);
#endif
  }

}

