#include "common_object.h"

void dump_averages(int itime){
  int i,j,k;
  FILE *fout;
  char fname[128];
  my_double x,y,z,ux,uy,uz,ux2,uy2,uz2,ene,rho,eps;


#ifdef LB_FLUID
    out_local.ux = out_local.uy = out_local.uz = 0.0; 
    out_local.ux2 = out_local.uy2 = out_local.uz2 = 0.0; 
    out_local.rho = 0.0;
    out_local.ene = 0.0;
    out_local.eps = 0.0;

    for (i = 0; i < NX; i++){ 
      ruler_x_local[i].x   = ruler_x_local[i].y   = ruler_x_local[i].z  = 0.0;
      ruler_x_local[i].ux  = ruler_x_local[i].uy  = ruler_x_local[i].uz  = 0.0;
      ruler_x_local[i].ux2 = ruler_x_local[i].uy2 = ruler_x_local[i].uz2 = 0.0;
      ruler_x_local[i].rho = ruler_x_local[i].ene = ruler_x_local[i].eps = 0.0;

      ruler_x[i].x   = ruler_x[i].y   = ruler_x[i].z  = 0.0;
      ruler_x[i].ux  = ruler_x[i].uy  = ruler_x[i].uz  = 0.0;
      ruler_x[i].ux2 = ruler_x[i].uy2 = ruler_x[i].uz2 = 0.0;
      ruler_x[i].rho = ruler_x[i].ene = ruler_x[i].eps = 0.0;
    }
    for (i = 0; i < NY; i++){ 
      ruler_y_local[i].x   = ruler_y_local[i].y   = ruler_y_local[i].z  = 0.0;
      ruler_y_local[i].ux  = ruler_y_local[i].uy  = ruler_y_local[i].uz  = 0.0;
      ruler_y_local[i].ux2 = ruler_y_local[i].uy2 = ruler_y_local[i].uz2 = 0.0;
      ruler_y_local[i].rho = ruler_y_local[i].ene = ruler_y_local[i].eps = 0.0;

      ruler_y[i].x   = ruler_y[i].y   = ruler_y[i].z  = 0.0;
      ruler_y[i].ux  = ruler_y[i].uy  = ruler_y[i].uz  = 0.0;
      ruler_y[i].ux2 = ruler_y[i].uy2 = ruler_y[i].uz2 = 0.0;
      ruler_y[i].rho = ruler_y[i].ene = ruler_y[i].eps = 0.0;
    } 
    for (i = 0; i < NZ; i++){ 
      ruler_z_local[i].x   = ruler_z_local[i].y   = ruler_z_local[i].z  = 0.0;
      ruler_z_local[i].ux  = ruler_z_local[i].uy  = ruler_z_local[i].uz  = 0.0;
      ruler_z_local[i].ux2 = ruler_z_local[i].uy2 = ruler_z_local[i].uz2 = 0.0;
      ruler_z_local[i].rho = ruler_z_local[i].ene = ruler_z_local[i].eps = 0.0;

      ruler_z[i].x   = ruler_z[i].y   = ruler_z[i].z  = 0.0;
      ruler_z[i].ux  = ruler_z[i].uy  = ruler_z[i].uz  = 0.0;
      ruler_z[i].ux2 = ruler_z[i].uy2 = ruler_z[i].uz2 = 0.0;
      ruler_z[i].rho = ruler_z[i].ene = ruler_z[i].eps = 0.0;
    }
#endif


	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {

#ifdef LB_FLUID
			  out_local.x += x = center_V[IDX(i, j, k)].x;
			  ruler_x_local[i -BRD + LNX_START].x += x;
			  ruler_y_local[j -BRD + LNY_START].x += x;
			  ruler_z_local[k -BRD + LNZ_START].x += x;

			  out_local.y += y = center_V[IDX(i, j, k)].y;
			  ruler_x_local[i -BRD + LNX_START].y += y;
			  ruler_y_local[j -BRD + LNY_START].y += y;
			  ruler_z_local[k -BRD + LNZ_START].y += y;

			  out_local.z += z = center_V[IDX(i, j, k)].z;
			  ruler_x_local[i -BRD + LNX_START].z += z;
			  ruler_y_local[j -BRD + LNY_START].z += z;
			  ruler_z_local[k -BRD + LNZ_START].z += z;

			  out_local.ux += ux = u[IDX(i, j, k)].x;
			  ruler_x_local[i -BRD + LNX_START].ux += ux;
			  ruler_y_local[j -BRD + LNY_START].ux += ux;
			  ruler_z_local[k -BRD + LNZ_START].ux += ux;

			  out_local.uy += uy = u[IDX(i, j, k)].y;
			  ruler_x_local[i -BRD + LNX_START].uy += uy;
			  ruler_y_local[j -BRD + LNY_START].uy += uy;
			  ruler_z_local[k -BRD + LNZ_START].uy += uy;

			  out_local.uz += uz = u[IDX(i, j, k)].z;
			  ruler_x_local[i -BRD + LNX_START].uz += uz;
			  ruler_y_local[j -BRD + LNY_START].uz += uz;
			  ruler_z_local[k -BRD + LNZ_START].uz += uz;

			  ux2 = u[IDX(i, j, k)].x*u[IDX(i, j, k)].x;
			  out_local.ux2 += ux2;
			  ruler_x_local[i -BRD + LNX_START].ux2 += ux2;
			  ruler_y_local[j -BRD + LNY_START].ux2 += ux2;
			  ruler_z_local[k -BRD + LNZ_START].ux2 += ux2;

			  uy2 = u[IDX(i, j, k)].y*u[IDX(i, j, k)].y;
			  out_local.uy2 += uy2;
			  ruler_x_local[i -BRD + LNX_START].uy2 += uy2;
			  ruler_y_local[j -BRD + LNY_START].uy2 += uy2;
			  ruler_z_local[k -BRD + LNZ_START].uy2 += uy2;

			  uz2 = u[IDX(i, j, k)].z*u[IDX(i, j, k)].z;
			  out_local.uz2 += uz2;
			  ruler_x_local[i -BRD + LNX_START].uz2 += uz2;
			  ruler_y_local[j -BRD + LNY_START].uz2 += uz2;
			  ruler_z_local[k -BRD + LNZ_START].uz2 += uz2;

			  out_local.rho += rho = dens[IDX(i, j, k)];
			  ruler_x_local[i -BRD + LNX_START].rho += rho;
			  ruler_y_local[j -BRD + LNY_START].rho += rho;
			  ruler_z_local[k -BRD + LNZ_START].rho += rho;

			  out_local.ene += ene = 0.5*(ux2+uy2+uz2);
			  ruler_x_local[i -BRD + LNX_START].ene += ene;
			  ruler_y_local[j -BRD + LNY_START].ene += ene;
			  ruler_z_local[k -BRD + LNZ_START].ene += ene;
			  
			  out_local.eps += eps = 0.0;
			  ruler_x_local[i -BRD + LNX_START].eps += eps;
			  ruler_y_local[j -BRD + LNY_START].eps += eps;
			  ruler_z_local[k -BRD + LNZ_START].eps += eps;
#endif
			} /* for i j k */      



#ifdef LB_FLUID
  MPI_Allreduce(&out_local, &out_all, 1, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );

  MPI_Allreduce(ruler_x_local, ruler_x, NX, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );
  MPI_Allreduce(ruler_y_local, ruler_y, NY, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );
  MPI_Allreduce(ruler_z_local, ruler_z, NZ, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );

  if(ROOT){
    sprintf(fname,"velocity_averages.dat");
    fout = fopen(fname,"a");
    fprintf(fout,"%d %e %e %e %e %e %e %e %e\n",itime, (double)out_all.ene, (double)out_all.rho, (double)out_all.ux, (double)out_all.uy, (double)out_all.uz, (double)out_all.ux2 , (double)out_all.uy2, (double)out_all.uz2);
    fclose(fout);

    sprintf(fname,"velocity_averages_x.dat");
    fout = fopen(fname,"w");
    for (i = 0; i < NX; i++) fprintf(fout,"%e %e %e %e %e %e %e %e %e\n",(double)ruler_x[i].x, (double)ruler_x[i].ene, (double)ruler_x[i].rho, (double)ruler_x[i].ux, (double)ruler_x[i].uy, (double)ruler_x[i].uz, (double)ruler_x[i].ux2 , (double)ruler_x[i].uy2, (double)ruler_x[i].uz2);
    fclose(fout);

    sprintf(fname,"velocity_averages_y.dat");
    fout = fopen(fname,"w");
    for (j = 0; j < NY; j++) fprintf(fout,"%e %e %e %e %e %e %e %e %e\n",(double)ruler_y[j].y, (double)ruler_y[j].ene, (double)ruler_y[j].rho, (double)ruler_y[j].ux, (double)ruler_y[j].uy, (double)ruler_y[j].uz, (double)ruler_y[j].ux2 , (double)ruler_y[j].uy2, (double)ruler_y[j].uz2);
    fclose(fout);

    sprintf(fname,"velocity_averages_z.dat");
    fout = fopen(fname,"w");
    for (k = 0; k < NZ; k++) fprintf(fout,"%e %e %e %e %e %e %e %e %e\n",(double)ruler_z[k].z, (double)ruler_z[k].ene, (double)ruler_z[k].rho, (double)ruler_z[k].ux, (double)ruler_z[k].uy, (double)ruler_z[k].uz, (double)ruler_z[k].ux2 , (double)ruler_z[k].uy2, (double)ruler_z[k].uz2);
    fclose(fout);

  }
#endif
  

}

