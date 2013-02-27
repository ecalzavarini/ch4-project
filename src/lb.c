#include "common_object.h"

void design_lb(){

	NPOP = 19;

	/* Settings for D3Q19 */
	wgt[0] = 1. / 3.;
	wgt[1] = 1. / 18.;
	wgt[2] = 1. / 18.;
	wgt[3] = 1. / 18.;
	wgt[4] = 1. / 18.;
	wgt[5] = 1. / 18.;
	wgt[6] = 1. / 18.;
	wgt[7] = 1. / 36.;
	wgt[8] = 1. / 36.;
	wgt[9] = 1. / 36.;
	wgt[10] = 1. / 36.;
	wgt[11] = 1. / 36.;
	wgt[12] = 1. / 36.;
	wgt[13] = 1. / 36.;
	wgt[14] = 1. / 36.;
	wgt[15] = 1. / 36.;
	wgt[16] = 1. / 36.;
	wgt[17] = 1. / 36.;
	wgt[18] = 1. / 36.;


	/* Lattice speeds, D3Q19 */
	c[0].x = 0.;	c[0].y = 0.;	c[0].z = 0.;
	c[1].x = 1.;	c[1].y = 0.;	c[1].z = 0.;
	c[2].x = -1.;	c[2].y = 0.;	c[2].z = 0.;
	c[3].x = 0.;	c[3].y = 1.;	c[3].z = 0.;
	c[4].x = 0.;	c[4].y = -1.;	c[4].z = 0.;
	c[5].x = 0.;	c[5].y = 0.;	c[5].z = 1.;
	c[6].x = 0.;	c[6].y = 0.;	c[6].z = -1.;
	c[7].x = 1.;	c[7].y = 1.;	c[7].z = 0.;
	c[8].x = 1.;	c[8].y = -1.;	c[8].z = 0.;
	c[9].x = -1.;	c[9].y = 1.;	c[9].z = 0.;
	c[10].x = -1.;	c[10].y = -1.;	c[10].z = 0.;
	c[11].x = 1.;	c[11].y = 0.;	c[11].z = 1.;
	c[12].x = -1.;	c[12].y = 0.;	c[12].z = 1.;
	c[13].x = 1.;	c[13].y = 0.;	c[13].z = -1.;
	c[14].x = -1.;	c[14].y = 0.;	c[14].z = -1.;
	c[15].x = 0.;	c[15].y = 1.;	c[15].z = 1.;
	c[16].x = 0.;	c[16].y = 1.;	c[16].z = -1.;
	c[17].x = 0.;	c[17].y = -1.;	c[17].z = 1.;
	c[18].x = 0.;	c[18].y = -1.;	c[18].z = -1.;


	/* in typedef.h  
#define mvx(a) (a.p[1] -a.p[2] +a.p[7]  +a.p[8] -a.p[9]  -a.p[10] +a.p[11] -a.p[12]+a.p[13]-a.p[14]) 
#define mvy(a) (a.p[3] -a.p[4] +a.p[7]  -a.p[8]  +a.p[9]  -a.p[10] +a.p[15]+a.p[16]-a.p[17]-a.p[18]) 
#define mvz(a) (a.p[5] -a.p[6] +a.p[11] +a.p[12] -a.p[13] -a.p[14] +a.p[15] -a.p[16]+a.p[17]-a.p[18])
#define m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8]+a.p[9]+a.p[10]+a.p[11]+a.p[12]+a.p[13]+a.p[14]+a.p[15]+a.p[16]+a.p[17]+a.p[18])
	 */

	/*
	 * This function provides the opposite velocity c_j=inv[i] for
	 * velocity c_i
	 */
	inv[0] = 0;
	inv[1] = 2;
	inv[2] = 1;
	inv[3] = 4;
	inv[4] = 3;
	inv[5] = 6;
	inv[6] = 5;

	inv[7] = 10;
	inv[8] = 9;
	inv[9] = 8;
	inv[10] = 7;
	inv[11] = 14;
	inv[12] = 13;
	inv[13] = 12;
	inv[14] = 11;
	inv[15] = 18;
	inv[16] = 17;
	inv[17] = 16;
	inv[18] = 15;

	/* speed of sound constants */
	cs = 1.0 / sqrt(3.0);
	invcs = 1.0 / cs;
	cs2 = (1.0 / 3.0);
	invcs2 = 1.0 / cs2;
	cs4 = (1.0 / 9.0);
	invcs4 = 1.0 / cs4;
	twocs2 = 2.0 * cs2;
	invtwocs2 = 1.0 / twocs2;
	twocs4 = 2.0 * cs4;
	invtwocs4 = 1.0 / twocs4;

}


pop equilibrium(pop * f, int i, int j, int k){
	int             pp;
	my_double       ux, uy, uz;
	my_double       rhof;
	my_double       cu, u2;
	pop             f_eq;

	rhof = m(f[IDX(i, j, k)]);

	ux = u[IDX(i, j, k)].x;
	uy = u[IDX(i, j, k)].y;
	uy = u[IDX(i, j, k)].z;

	u2 = (ux * ux + uy * uy + uz * uz);

	/* equilibrium distribution */
	for (pp = 0; pp < 9; pp++) {
		cu = (c[pp].x * ux + c[pp].y * uy + c[pp].z * uz);
		f_eq.p[pp] = rhof * wgt[pp] * (1.0 + invcs2 * cu + invtwocs4 * cu * cu - invtwocs2 * u2);
	}

	return f_eq;
}

/**************************************************/
void hydro_fields(){
	int i, j, k;
#ifdef DEBUG
	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;
#endif

	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {

#ifdef LB_FLUID
				dens[IDX(i, j, k)] = m(p[IDX(i, j, k)]);
#ifdef METHOD_FORCING_GUO
				u[IDX(i, j, k)].x = (mvx(p[IDX(i, j, k)]) + 0.5 * force[IDX(i, j, k)].x) / dens[IDX(i, j, k)];
				u[IDX(i, j, k)].y = (mvy(p[IDX(i, j, k)]) + 0.5 * force[IDX(i, j, k)].y) / dens[IDX(i, j, k)];
				u[IDX(i, j, k)].z = (mvz(p[IDX(i, j, k)]) + 0.5 * force[IDX(i, j, k)].z) / dens[IDX(i, j, k)];
				/* set to zero after computing velocity */
				force[IDX(i, j, k)].x = force[IDX(i, j, k)].y = force[IDX(i, j, k)].y = 0.0;
#else
				u[IDX(i, j, k)].x = mvx(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
				u[IDX(i, j, k)].y = mvy(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
				u[IDX(i, j, k)].z = mvz(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
#endif
#endif
			}/* for i,j,k */

#ifdef DEBUG
	/* Each processor prints its mesh */
	sprintf(fnamein, "velocity.%d.out", me);
	fout = fopen(fnamein, "w");
	
	for (k = BRD; k < LNZ + BRD; k++)
		for (j = BRD; j < LNY + BRD; j++)
		for (i = BRD; i < LNX + BRD; i++)
		/*
	for (k = 0; k < LNZ + TWO_BRD; k++)
	     for (j = 0; j < LNY + TWO_BRD; j++)
		  for (i = 0; i < LNX + TWO_BRD; i++)
		    */
		  fprintf(fout, "%d %d %d %e %e %e %e\n", i, j, k, center_V[IDX(i, j, k)].x, center_V[IDX(i, j, k)].y, center_V[IDX(i, j, k)].z, 
			                                           u[IDX(i, j, k)].x, u[IDX(i, j, k)].y, u[IDX(i, j, k)].z , dens[IDX(i, j, k)]);
	fclose(fout);
#endif

}

/********************************************/
void time_stepping(pop *f, pop *rhs_f){
  int i,j,k,pp;

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

	for(pp=0;pp<NPOP;pp++){

	  /* Euler first order */
	  f[IDX(i,j,k)].p[pp] += property.time_dt*rhs_f[IDX(i,j,k)].p[pp];

	}

      }/* for i,j,k */
}



/********************************************/ 
void sendrecv_borders_pop(pop *f){
  int i,j,k,brd_size;
  MPI_Status status1;

  /*     BRD|LNX|BRD     */
  /* Copy borders along x */
  brd_size = BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_pop[IDX_XBRD(i,j,k)] = f[IDX(i+LNX,j,k)];
      }

  MPI_Sendrecv( xp_pop, brd_size, MPI_vector_type, me_xp, 10,
                xm_pop, brd_size, MPI_vector_type, me_xm, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_pop[IDX_XBRD(i,j,k)];
        xm_pop[IDX_XBRD(i,j,k)] = f[IDX(i+BRD,j,k)];
      }
 MPI_Sendrecv( xm_pop, brd_size, MPI_vector_type, me_xm, 10,
               xp_pop, brd_size, MPI_vector_type, me_xp, 10, MPI_COMM_WORLD, &status1);

 for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k)] = xp_pop[IDX_XBRD(i,j,k)];
      }



  /* Copy borders along x */
  brd_size = BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){ 
        yp_pop[IDX_YBRD(i,j,k)] = f[IDX(i,j+LNY,k)];
      }

  MPI_Sendrecv( yp_pop, brd_size, MPI_vector_type, me_yp, 10,
                ym_pop, brd_size, MPI_vector_type, me_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++) {
        f[IDX(i,j,k)] = ym_pop[IDX_YBRD(i,j,k)];
        ym_pop[IDX_YBRD(i,j,k)] = f[IDX(i,j+BRD,k)];
      }
 MPI_Sendrecv( ym_pop, brd_size, MPI_vector_type, me_ym, 10,
               yp_pop, brd_size, MPI_vector_type, me_yp, 10, MPI_COMM_WORLD, &status1);

 for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){ 
	f[IDX(i,j+LNY+BRD,k)] = yp_pop[IDX_YBRD(i,j,k)];
      }


  /* Copy borders along x */
  brd_size = BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){ 
        zp_pop[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+LNZ)];
      }

  MPI_Sendrecv( zp_pop, brd_size, MPI_vector_type, me_zp, 10,
                zm_pop, brd_size, MPI_vector_type, me_zm, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++) {
        f[IDX(i,j,k)] = zm_pop[IDX_ZBRD(i,j,k)];
        zm_pop[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+BRD)];
      }
 MPI_Sendrecv( zm_pop, brd_size, MPI_vector_type, me_zm, 10,
               zp_pop, brd_size, MPI_vector_type, me_zp, 10, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){ 
	f[IDX(i,j,k+LNZ+BRD)] = zp_pop[IDX_ZBRD(i,j,k)];
      }

}/* end send rcv */
