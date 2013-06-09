#include "common_object.h"


void dump_averages(){
  int i,j,k;
  FILE *fout;
  char fname[128];
  my_double x,y,z,ux,uy,uz,ux2,uy2,uz2,ene,rho,eps;
  my_double norm;
  tensor S;
#ifdef LB_TEMPERATURE
  my_double temp,t2,epst,dxt,dyt,dzt,uxt,uyt,uzt,nux,nuy,nuz;
  vector grad;
#ifdef LB_TEMPERATURE_MELTING   
  my_double lf, dtlf, enth;
#endif 
#endif


  x=y=z=ux=uy=uz=ux2=uy2=uz2=ene=rho=eps=0.0;

#ifdef LB_TEMPERATURE
temp = t2 = epst = dxt = dyt = dzt = uxt= uyt = uzt = nux = nuy = nuz= 0.0;
#endif

#ifdef LB_FLUID
    out_local.x = out_local.y = out_local.z = 0.0; 
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


#ifdef LB_TEMPERATURE
    out_local.dxt = out_local.dyt = out_local.dzt = 0.0; 
    out_local.uxt = out_local.uyt = out_local.uzt = 0.0; 
    out_local.nux = out_local.nuy = out_local.nuz = 0.0; 
    out_local.t = 0.0;
    out_local.t2 = 0.0;
    out_local.epst = 0.0;

    for (i = 0; i < NX; i++){ 
      ruler_x_local[i].dxt   = ruler_x_local[i].dyt   = ruler_x_local[i].dzt  = 0.0;
      ruler_x_local[i].uxt  = ruler_x_local[i].uyt  = ruler_x_local[i].uzt  = 0.0;
      ruler_x_local[i].nux = ruler_x_local[i].nuy = ruler_x_local[i].nuz = 0.0;
      ruler_x_local[i].t = ruler_x_local[i].t2 = ruler_x_local[i].epst = 0.0;

      ruler_x[i].dxt   = ruler_x[i].dyt   = ruler_x[i].dzt  = 0.0;
      ruler_x[i].uxt  = ruler_x[i].uyt  = ruler_x[i].uzt  = 0.0;
      ruler_x[i].nux = ruler_x[i].nuy = ruler_x[i].nuz = 0.0;
      ruler_x[i].t = ruler_x[i].t2 = ruler_x[i].epst = 0.0;
    }
    for (i = 0; i < NY; i++){ 
      ruler_y_local[i].dxt   = ruler_y_local[i].dyt   = ruler_y_local[i].dzt  = 0.0;
      ruler_y_local[i].uxt  = ruler_y_local[i].uyt  = ruler_y_local[i].uzt  = 0.0;
      ruler_y_local[i].nux = ruler_y_local[i].nuy = ruler_y_local[i].nuz = 0.0;
      ruler_y_local[i].t = ruler_y_local[i].t2 = ruler_y_local[i].epst = 0.0;

      ruler_y[i].dxt   = ruler_y[i].dyt   = ruler_y[i].dzt  = 0.0;
      ruler_y[i].uxt  = ruler_y[i].uyt  = ruler_y[i].uzt  = 0.0;
      ruler_y[i].nux = ruler_y[i].nuy = ruler_y[i].nuz = 0.0;
      ruler_y[i].t = ruler_y[i].t2 = ruler_y[i].epst = 0.0;
    } 
    for (i = 0; i < NZ; i++){ 
      ruler_z_local[i].dxt   = ruler_z_local[i].dyt   = ruler_z_local[i].dzt  = 0.0;
      ruler_z_local[i].uxt  = ruler_z_local[i].uyt  = ruler_z_local[i].uzt  = 0.0;
      ruler_z_local[i].nux = ruler_z_local[i].nuy = ruler_z_local[i].nuz = 0.0;
      ruler_z_local[i].t = ruler_z_local[i].t2 = ruler_z_local[i].epst = 0.0;

      ruler_z[i].dxt   = ruler_z[i].dyt   = ruler_z[i].dzt  = 0.0;
      ruler_z[i].uxt  = ruler_z[i].uyt  = ruler_z[i].uzt  = 0.0;
      ruler_z[i].nux = ruler_z[i].nuy = ruler_z[i].nuz = 0.0;
      ruler_z[i].t = ruler_z[i].t2 = ruler_z[i].epst = 0.0;
    }
#ifdef LB_TEMPERATURE_MELTING
    out_local.lf = out_local.dtlf = out_local.enth = 0.0; 

    for (i = 0; i < NX; i++){ 
      ruler_x_local[i].lf   = ruler_x_local[i].dtlf   = ruler_x_local[i].enth  = 0.0;
      ruler_x[i].lf   = ruler_x[i].dtlf   = ruler_x[i].enth  = 0.0;
    }

    for (i = 0; i < NY; i++){ 
      ruler_y_local[i].lf   = ruler_y_local[i].dtlf   = ruler_y_local[i].enth  = 0.0;
      ruler_y[i].lf   = ruler_y[i].dtlf   = ruler_y[i].enth  = 0.0;
    } 

    for (i = 0; i < NZ; i++){ 
      ruler_z_local[i].lf   = ruler_z_local[i].dtlf   = ruler_z_local[i].enth  = 0.0;
      ruler_z[i].lf   = ruler_z[i].dtlf   = ruler_z[i].enth  = 0.0;
    }
#endif
#endif


    /* Here we send recv to compute the gradients */
#ifdef LB_FLUID
    sendrecv_borders_vector(u);
#endif
#ifdef LB_TEMPERATURE
    sendrecv_borders_scalar(t);
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
			  
			  S=strain_tensor(p,i, j, k);
			  // fprintf(stderr,"SXY %e\n",S.xy);
			  out_local.eps += eps = S.xy;
			  ruler_x_local[i -BRD + LNX_START].eps += eps;
			  ruler_y_local[j -BRD + LNY_START].eps += eps;
			  ruler_z_local[k -BRD + LNZ_START].eps += eps;
#endif


#ifdef LB_TEMPERATURE
			  /* NB: this gradient definition assumes a (refined or not) cartesian grid */
			  /* note thet t field must be communicate first -> it has been done above */
			  grad=gradient_scalar(t,i,j,k);

			  dxt = grad.x;
			  out_local.dxt += dxt;
			  ruler_x_local[i -BRD + LNX_START].dxt += dxt;
			  ruler_y_local[j -BRD + LNY_START].dxt += dxt;
			  ruler_z_local[k -BRD + LNZ_START].dxt += dxt;

			  dyt = grad.y;
			  out_local.dyt += dyt;
			  ruler_x_local[i -BRD + LNX_START].dyt += dyt;
			  ruler_y_local[j -BRD + LNY_START].dyt += dyt;
			  ruler_z_local[k -BRD + LNZ_START].dyt += dyt;

			  dzt = grad.z;
			  out_local.dzt += dzt;
			  ruler_x_local[i -BRD + LNX_START].dzt += dzt;
			  ruler_y_local[j -BRD + LNY_START].dzt += dzt;
			  ruler_z_local[k -BRD + LNZ_START].dzt += dzt;
			  

			  out_local.uxt += uxt = u[IDX(i, j, k)].x*t[IDX(i, j, k)];
			  ruler_x_local[i -BRD + LNX_START].uxt += uxt;
			  ruler_y_local[j -BRD + LNY_START].uxt += uxt;
			  ruler_z_local[k -BRD + LNZ_START].uxt += uxt;

			  out_local.uyt += uyt = u[IDX(i, j, k)].y*t[IDX(i, j, k)];
			  ruler_x_local[i -BRD + LNX_START].uyt += uyt;
			  ruler_y_local[j -BRD + LNY_START].uyt += uyt;
			  ruler_z_local[k -BRD + LNZ_START].uyt += uyt;

			  out_local.uzt += uzt = u[IDX(i, j, k)].z*t[IDX(i, j, k)];
			  ruler_x_local[i -BRD + LNX_START].uzt += uzt;
			  ruler_y_local[j -BRD + LNY_START].uzt += uzt;
			  ruler_z_local[k -BRD + LNZ_START].uzt += uzt;

			  nux = ( - property.kappa*dxt + uxt );
			  out_local.nux += nux;
			  ruler_x_local[i -BRD + LNX_START].nux += nux;
			  ruler_y_local[j -BRD + LNY_START].nux += nux;
			  ruler_z_local[k -BRD + LNZ_START].nux += nux;

			  nuy = ( - property.kappa*dyt + uyt );
			  out_local.nuy += nuy;
			  ruler_x_local[i -BRD + LNX_START].nuy += nuy;
			  ruler_y_local[j -BRD + LNY_START].nuy += nuy;
			  ruler_z_local[k -BRD + LNZ_START].nuy += nuy;

			  nuz = ( - property.kappa*dzt + uzt );
			  out_local.nuz += nuz;
			  ruler_x_local[i -BRD + LNX_START].nuz += nuz;
			  ruler_y_local[j -BRD + LNY_START].nuz += nuz;
			  ruler_z_local[k -BRD + LNZ_START].nuz += nuz;

			  out_local.t += temp = t[IDX(i, j, k)];
			  ruler_x_local[i -BRD + LNX_START].t += temp;
			  ruler_y_local[j -BRD + LNY_START].t += temp;
			  ruler_z_local[k -BRD + LNZ_START].t += temp;

			  t2=t[IDX(i, j, k)]*t[IDX(i, j, k)];
			  out_local.t2 += t2;
			  ruler_x_local[i -BRD + LNX_START].t2 += t2;
			  ruler_y_local[j -BRD + LNY_START].t2 += t2;
			  ruler_z_local[k -BRD + LNZ_START].t2 += t2;
			  
			  epst = property.kappa*(dxt*dxt + dyt*dyt + dzt*dzt);
			  out_local.epst += epst; 
			  ruler_x_local[i -BRD + LNX_START].epst += epst;
			  ruler_y_local[j -BRD + LNY_START].epst += epst;
			  ruler_z_local[k -BRD + LNZ_START].epst += epst;

#ifdef LB_TEMPERATURE_MELTING   
			  lf = liquid_frac[IDX(i, j, k)];
			  out_local.lf += lf;
			  ruler_x_local[i -BRD + LNX_START].lf += lf;
			  ruler_y_local[j -BRD + LNY_START].lf += lf;
			  ruler_z_local[k -BRD + LNZ_START].lf += lf;

			  dtlf = ( liquid_frac[IDX(i, j, k)]-liquid_frac_old[IDX(i, j, k)] )/property.time_dt;
			  out_local.dtlf += dtlf;
			  ruler_x_local[i -BRD + LNX_START].dtlf += dtlf;
			  ruler_y_local[j -BRD + LNY_START].dtlf += dtlf;
			  ruler_z_local[k -BRD + LNZ_START].dtlf += dtlf;

			  enth = property.specific_heat*t[IDX(i,j,k)] + property.latent_heat*liquid_frac[IDX(i,j,k)];
			  out_local.enth += enth;
			  ruler_x_local[i -BRD + LNX_START].enth += enth;
			  ruler_y_local[j -BRD + LNY_START].enth += enth;
			  ruler_z_local[k -BRD + LNZ_START].enth += enth;
#endif 

#endif

			} /* for i j k */      




	/* Sum all */
  MPI_Allreduce(&out_local, &out_all, 1, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );

  MPI_Allreduce(ruler_x_local, ruler_x, NX, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );
  MPI_Allreduce(ruler_y_local, ruler_y, NY, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );
  MPI_Allreduce(ruler_z_local, ruler_z, NZ, MPI_output_type, MPI_SUM_output, MPI_COMM_WORLD );

#ifdef LB_FLUID

#define OUTPUT_NORM
#ifdef OUTPUT_NORM
  /* normalization */
  norm = 1.0/(my_double)(NX*NY*NZ);
  out_all.ene *= norm;
  out_all.rho *= norm;
  out_all.ux  *= norm;
  out_all.uy  *= norm;
  out_all.uz  *= norm;
  out_all.ux2 *= norm;
  out_all.uy2 *= norm;
  out_all.uz2 *= norm;
  out_all.eps *= norm;


  norm = 1.0/(my_double)(NY*NZ);
  for (i = 0; i < NX; i++){
    ruler_x[i].x *= norm;
    ruler_x[i].ene *= norm;
    ruler_x[i].rho *= norm;
    ruler_x[i].ux *= norm; 
    ruler_x[i].uy *= norm;
    ruler_x[i].uz *= norm;
    ruler_x[i].ux2 *= norm;
    ruler_x[i].uy2 *= norm;
    ruler_x[i].uz2 *= norm;
    ruler_x[i].eps *= norm;
  }

  norm = 1.0/(my_double)(NX*NZ);
  for (i = 0; i < NY; i++){
    ruler_y[i].y *= norm;
    ruler_y[i].ene *= norm;
    ruler_y[i].rho *= norm;
    ruler_y[i].ux *= norm; 
    ruler_y[i].uy *= norm;
    ruler_y[i].uz *= norm;
    ruler_y[i].ux2 *= norm;
    ruler_y[i].uy2 *= norm;
    ruler_y[i].uz2 *= norm;
    ruler_y[i].eps *= norm;
  }

  norm = 1.0/(my_double)(NX*NY);
  for (i = 0; i < NZ; i++){
    ruler_z[i].z *= norm;
    ruler_z[i].ene *= norm;
    ruler_z[i].rho *= norm;
    ruler_z[i].ux *= norm; 
    ruler_z[i].uy *= norm;
    ruler_z[i].uz *= norm;
    ruler_z[i].ux2 *= norm;
    ruler_z[i].uy2 *= norm;
    ruler_z[i].uz2 *= norm;
    ruler_z[i].eps *= norm;
  }
#endif
#endif 

#ifdef LB_TEMPERATURE

#ifdef OUTPUT_NORM
  /* normalization */
  norm = 1.0/(my_double)(NX*NY*NZ);
  out_all.t *= norm;
  out_all.t2 *= norm;
  out_all.epst  *= norm;
  out_all.dxt  *= norm;
  out_all.dyt  *= norm;
  out_all.dzt *= norm;
  out_all.uxt *= norm;
  out_all.uyt *= norm;
  out_all.uzt *= norm;
  out_all.nux *= norm/( property.kappa*property.deltaT/property.SY );
  out_all.nuy *= norm/( property.kappa*property.deltaT/property.SY );
  out_all.nuz *= norm/( property.kappa*property.deltaT/property.SY );

  norm = 1.0/(my_double)(NY*NZ);
  for (i = 0; i < NX; i++){
    ruler_x[i].t *= norm;
    ruler_x[i].t2 *= norm;
    ruler_x[i].epst *= norm;
    ruler_x[i].dxt *= norm; 
    ruler_x[i].dyt *= norm;
    ruler_x[i].dzt *= norm;
    ruler_x[i].uxt *= norm;
    ruler_x[i].uyt *= norm;
    ruler_x[i].uzt *= norm;
    ruler_x[i].nux *= norm/( property.kappa*property.deltaT/property.SY );
    ruler_x[i].nuy *= norm/( property.kappa*property.deltaT/property.SY );
    ruler_x[i].nuz *= norm/( property.kappa*property.deltaT/property.SY );
  }

  norm = 1.0/(my_double)(NX*NZ);
  for (i = 0; i < NY; i++){
    ruler_y[i].t *= norm;
    ruler_y[i].t2 *= norm;
    ruler_y[i].epst *= norm;
    ruler_y[i].dxt *= norm; 
    ruler_y[i].dyt *= norm;
    ruler_y[i].dzt *= norm;
    ruler_y[i].uxt *= norm;
    ruler_y[i].uyt *= norm;
    ruler_y[i].uzt *= norm;
    ruler_y[i].nux *= norm/( property.kappa*property.deltaT/property.SY );
    ruler_y[i].nuy *= norm/( property.kappa*property.deltaT/property.SY );
    ruler_y[i].nuz *= norm/( property.kappa*property.deltaT/property.SY );
  }

  norm = 1.0/(my_double)(NX*NY);
  for (i = 0; i < NZ; i++){
    ruler_z[i].t *= norm;
    ruler_z[i].t2 *= norm;
    ruler_z[i].epst *= norm;
    ruler_z[i].dxt *= norm; 
    ruler_z[i].dyt *= norm;
    ruler_z[i].dzt *= norm;
    ruler_z[i].uxt *= norm;
    ruler_z[i].uyt *= norm;
    ruler_z[i].uzt *= norm;
    ruler_z[i].nux *= norm/( property.kappa*property.deltaT/property.SY );
    ruler_z[i].nuy *= norm/( property.kappa*property.deltaT/property.SY );
    ruler_z[i].nuz *= norm/( property.kappa*property.deltaT/property.SY );
  }
#endif

#ifdef LB_TEMPERATURE_MELTING
  /* normalization */
  norm = 1.0/(my_double)(NX*NY*NZ);
  out_all.lf *= norm;
  out_all.dtlf *= norm;
  out_all.enth  *= norm;

  norm = 1.0/(my_double)(NY*NZ);
  for (i = 0; i < NX; i++){
    ruler_x[i].lf *= norm;
    ruler_x[i].dtlf *= norm;
    ruler_x[i].enth *= norm;
  }

  norm = 1.0/(my_double)(NX*NZ);
  for (i = 0; i < NY; i++){
    ruler_y[i].lf *= norm;
    ruler_y[i].dtlf *= norm;
    ruler_y[i].enth *= norm;
  }

  norm = 1.0/(my_double)(NX*NY);
  for (i = 0; i < NZ; i++){
    ruler_z[i].lf *= norm;
    ruler_z[i].dtlf *= norm;
    ruler_z[i].enth *= norm;
  }
#endif
#endif 



#ifdef LB_FLUID
  if(ROOT){
    sprintf(fname,"velocity_averages.dat");
    fout = fopen(fname,"a");    
    fprintf(fout,"%e %e %e %e %e %e %e %e %e %e\n",time_now, (double)out_all.ene,(double)out_all.rho, (double)out_all.ux, (double)out_all.uy, (double)out_all.uz, (double)out_all.ux2 , (double)out_all.uy2, (double)out_all.uz2, (double)out_all.eps);
    fclose(fout);

    sprintf(fname,"velocity_averages_x.dat");
    fout = fopen(fname,"w");
    for (i = 0; i < NX; i++) fprintf(fout,"%e %e %e %e %e %e %e %e %e %e\n",(double)ruler_x[i].x, (double)ruler_x[i].ene, (double)ruler_x[i].rho, (double)ruler_x[i].ux, (double)ruler_x[i].uy, (double)ruler_x[i].uz, (double)ruler_x[i].ux2 , (double)ruler_x[i].uy2, (double)ruler_x[i].uz2, (double)ruler_x[i].eps);
    fclose(fout);

    sprintf(fname,"velocity_averages_y.dat");
    fout = fopen(fname,"w");
    for (j = 0; j < NY; j++) fprintf(fout,"%e %e %e %e %e %e %e %e %e %e\n",(double)ruler_y[j].y, (double)ruler_y[j].ene, (double)ruler_y[j].rho, (double)ruler_y[j].ux, (double)ruler_y[j].uy, (double)ruler_y[j].uz, (double)ruler_y[j].ux2 , (double)ruler_y[j].uy2, (double)ruler_y[j].uz2,(double)ruler_y[i].eps);
    fclose(fout);

    sprintf(fname,"velocity_averages_z.dat");
    fout = fopen(fname,"w");
    for (k = 0; k < NZ; k++) fprintf(fout,"%e %e %e %e %e %e %e %e %e %e\n",(double)ruler_z[k].z, (double)ruler_z[k].ene, (double)ruler_z[k].rho, (double)ruler_z[k].ux, (double)ruler_z[k].uy, (double)ruler_z[k].uz, (double)ruler_z[k].ux2 , (double)ruler_z[k].uy2, (double)ruler_z[k].uz2, (double)ruler_z[i].eps);
    fclose(fout);

  }
#endif
  
#ifdef LB_TEMPERATURE
  if(ROOT){
    sprintf(fname,"temperature_averages.dat");
    fout = fopen(fname,"a");    
    fprintf(fout,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",time_now, (double)out_all.t,(double)out_all.t2, (double)out_all.epst, (double)out_all.dxt, (double)out_all.dyt, (double)out_all.dzt , (double)out_all.uxt, (double)out_all.uyt, (double)out_all.uzt,(double)out_all.nux, (double)out_all.nuy, (double)out_all.nuz );
    fclose(fout);

    sprintf(fname,"temperature_averages_x.dat");
    fout = fopen(fname,"w");
    for (i = 0; i < NX; i++) fprintf(fout,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",(double)ruler_x[i].x, (double)ruler_x[i].t, (double)ruler_x[i].t2, (double)ruler_x[i].epst, (double)ruler_x[i].dxt, (double)ruler_x[i].dyt, (double)ruler_x[i].dzt , (double)ruler_x[i].uxt, (double)ruler_x[i].uyt, (double)ruler_x[i].uzt, (double)ruler_x[i].nux, (double)ruler_x[i].nuy, (double)ruler_x[i].nuz);
    fclose(fout);

    sprintf(fname,"temperature_averages_y.dat");
    fout = fopen(fname,"w");
    for (j = 0; j < NY; j++) fprintf(fout,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",(double)ruler_y[j].y, (double)ruler_y[j].t, (double)ruler_y[j].t2, (double)ruler_y[j].epst, (double)ruler_y[j].dxt, (double)ruler_y[j].dyt, (double)ruler_y[j].dzt , (double)ruler_y[j].uxt, (double)ruler_y[j].uyt,(double)ruler_y[i].uzt, (double)ruler_y[j].nux, (double)ruler_y[j].nuy, (double)ruler_y[j].nuz);
    fclose(fout);

    sprintf(fname,"temperature_averages_z.dat");
    fout = fopen(fname,"w");
    for (k = 0; k < NZ; k++) fprintf(fout,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",(double)ruler_z[k].z, (double)ruler_z[k].t, (double)ruler_z[k].t2, (double)ruler_z[k].epst, (double)ruler_z[k].dxt, (double)ruler_z[k].dyt, (double)ruler_z[k].dzt , (double)ruler_z[k].uxt, (double)ruler_z[k].uyt, (double)ruler_z[i].uzt, (double)ruler_z[k].nux, (double)ruler_z[k].nuy, (double)ruler_z[k].nuz );
    fclose(fout);

  }

#ifdef LB_TEMPERATURE_MELTING
  if(ROOT){
    sprintf(fname,"melting_averages.dat");
    fout = fopen(fname,"a");    
    fprintf(fout,"%e %e %e %e %e\n",time_now, (double)out_all.lf,(double)out_all.dtlf, (double)out_all.enth);
    fclose(fout);

    sprintf(fname,"melting_averages_x.dat");
    fout = fopen(fname,"w");
    for (i = 0; i < NX; i++) fprintf(fout,"%e %e %e %e %e\n",(double)ruler_x[i].x, (double)ruler_x[i].lf, (double)ruler_x[i].dtlf, (double)ruler_x[i].enth);
    fclose(fout);

    sprintf(fname,"melting_averages_y.dat");
    fout = fopen(fname,"w");
    for (j = 0; j < NY; j++) fprintf(fout,"%e %e %e %e %e\n",(double)ruler_y[j].y, (double)ruler_y[j].lf, (double)ruler_y[j].dtlf, (double)ruler_y[j].enth);
    fclose(fout);

    sprintf(fname,"melting_averages_z.dat");
    fout = fopen(fname,"w");
    for (k = 0; k < NZ; k++) fprintf(fout,"%e %e %e %e %e\n",(double)ruler_z[k].z, (double)ruler_z[k].lf, (double)ruler_z[k].dtlf, (double)ruler_z[k].enth);
    fclose(fout);
  }
#endif
#endif

#ifdef OUTPUT_H5
    if(itime%((int)(property.time_dump_field/property.time_dt))==0) output_h5();
#endif


#ifdef OUTPUT_ASCII
#ifdef LB_FLUID

  if(ROOT && itime%((int)(property.time_dump_field/property.time_dt))==0 ){
  /* Here dumps the velocity field */
    sprintf(fname,"%s/vel.%d",OutDir,itime);
  fout = fopen(fname,"w");

  for(k=BRD;k<LNZ+BRD;k++){
    for(j=BRD;j<LNY+BRD;j++){
      for(i=BRD;i<LNX+BRD;i++){ 

        fprintf(fout,"%e %e %e %e %e %e %e\n", 
		(double)center_V[IDX(i, j, k)].x, (double)center_V[IDX(i, j, k)].y, (double)center_V[IDX(i, j, k)].z, 
                (double)u[IDX(i,j,k)].x, (double)u[IDX(i,j,k)].y, (double)u[IDX(i,j,k)].z,
		(double)dens[IDX(i,j,k)] );  
      } 
      fprintf(fout,"\n");
    }
    fprintf(fout,"\n");
  }
    fclose(fout);
  }
#endif

#ifdef LB_TEMPERATURE
  if(ROOT  && itime%((int)(property.time_dump_field/property.time_dt))==0){
  /* Here dumps the temperature field */
    sprintf(fname,"%s/temp.%d",OutDir,itime);
  fout = fopen(fname,"w");

  for(k=BRD;k<LNZ+BRD;k++){
    for(j=BRD;j<LNY+BRD;j++){
      for(i=BRD;i<LNX+BRD;i++){ 

        fprintf(fout,"%e %e %e %e\n", 
		(double)center_V[IDX(i, j, k)].x, (double)center_V[IDX(i, j, k)].y, (double)center_V[IDX(i, j, k)].z, 
		(double)t[IDX(i,j,k)] );  
      } 
      fprintf(fout,"\n");
    }
    fprintf(fout,"\n");
  }
    fclose(fout);
  }
#endif

#endif
}

