#include "common_object.h"

void design_lb(){

  /* 4 lines below, to be removed */
  //NPOP = 19;
  /* commit pop type */
  //MPI_Type_contiguous(NPOP, MPI_DOUBLE, &MPI_pop_type);
  //MPI_Type_commit(&MPI_pop_type);

#ifdef GRID_POP_D2Q9
	wgt[0] = 4. / 9.;
	wgt[1] = 1. / 36.;
	wgt[2] = 1. / 9.;
	wgt[3] = 1. / 36.;
	wgt[4] = 1. / 9.;
	wgt[5] = 1. / 36.;
	wgt[6] = 1. / 9.;
	wgt[7] = 1. / 36.;
	wgt[8] = 1. / 9.;

	c[0].x = 0.;	c[0].y = 0.;   c[0].z=0.;
	c[1].x = -1.;	c[1].y = 1.;   c[1].z=0.;
	c[2].x = -1.;	c[2].y = 0.;   c[2].z=0.;
	c[3].x = -1.;	c[3].y = -1.;  c[3].z=0.;
	c[4].x = 0.;	c[4].y = -1.;  c[4].z=0.;
	c[5].x = 1.;	c[5].y = -1.;  c[5].z=0.;
	c[6].x = 1.;	c[6].y = 0.;   c[6].z=0.;
	c[7].x = 1.;	c[7].y = 1.;   c[7].z=0.;
	c[8].x = 0.;	c[8].y = 1.;   c[8].z=0.;

	inv[0] = 0;
	inv[1] = 5;
	inv[2] = 6;
	inv[3] = 7;
	inv[4] = 8;
	inv[5] = 1;
	inv[6] = 2;
	inv[7] = 3;
	inv[8] = 4;

	inv_x[0] = 0;	       inv_y[0] = 0;   inv_z[0] = 0;
	inv_x[1] = 7;	       inv_y[1] = 3;   inv_z[1] = 1;
	inv_x[2] = 6;	       inv_y[2] = 2;   inv_z[2] = 2;
	inv_x[3] = 5;	       inv_y[3] = 1;   inv_z[3] = 3;
	inv_x[4] = 4;	       inv_y[4] = 8;   inv_z[4] = 4;
	inv_x[5] = 3;	       inv_y[5] = 7;   inv_z[5] = 5;
	inv_x[6] = 2;	       inv_y[6] = 6;   inv_z[6] = 6;
	inv_x[7] = 1;	       inv_y[7] = 5;   inv_z[7] = 7;
	inv_x[8] = 8;	       inv_y[8] = 4;   inv_z[8] = 8;

#endif

#ifdef GRID_POP_D3Q15
	wgt[0] = 2. / 9.;
	wgt[1] = 1. / 9.;
	wgt[2] = 1. / 9.;
	wgt[3] = 1. / 9.;
	wgt[4] = 1. / 72.;
	wgt[5] = 1. / 72.;
	wgt[6] = 1. / 72.;
	wgt[7] = 1. / 72.;
	wgt[8] = 1. / 9.;
	wgt[9] = 1. / 9.;
	wgt[10] = 1. / 9.;
	wgt[11] = 1. / 72.;
	wgt[12] = 1. / 72.;
	wgt[13] = 1. / 72.;
	wgt[14] = 1. / 72.;


	/* Lattice speeds, D3Q15 */
	c[0].x = 0.;	c[0].y = 0.;	c[0].z = 0.;
	c[1].x = -1.;	c[1].y = 0.;	c[1].z = 0.;
	c[2].x = 0.;	c[2].y = -1.;	c[2].z = 0.;
	c[3].x = 0.;	c[3].y = 0.;	c[3].z = -1.;
	c[4].x = -1.;	c[4].y = -1.;	c[4].z = -1.;
	c[5].x = -1.;	c[5].y = -1.;	c[5].z = 1.;
	c[6].x = -1.;	c[6].y = 1.;	c[6].z = -1.;
	c[7].x = -1.;	c[7].y = 1.;	c[7].z = 1.;
	c[8].x = 1.;	c[8].y = 0.;	c[8].z = 0.;
	c[9].x = 0.;	c[9].y = 1.;	c[9].z = 0.;
	c[10].x = 0.;	c[10].y = 0.;	c[10].z = 1.;
	c[11].x = 1.;	c[11].y = 1.;	c[11].z = 1.;
	c[12].x = 1.;	c[12].y = 1.;	c[12].z = -1.;
	c[13].x = 1.;	c[13].y = -1.;	c[13].z = 1.;
	c[14].x = 1.;	c[14].y = -1.;	c[14].z = -1.;

	/*
	 * This function provides the opposite velocity c_j=inv[i] for
	 * velocity c_i
	 */
	inv[0] = 0;
	inv[1] = 8;
	inv[2] = 9;
	inv[3] = 10;
	inv[4] = 11;
	inv[5] = 12;
	inv[6] = 13;
	inv[7] = 14;
	inv[8] = 1;
	inv[9] = 2;
	inv[10] = 3;
	inv[11] = 4;
	inv[12] = 5;
	inv[13] = 6;
	inv[14] = 7;


	inv_x[0] = 0;	       inv_y[0] = 0;	     inv_z[0] = 0;
	inv_x[1] = 8;	       inv_y[1] = 1;	     inv_z[1] = 1;
	inv_x[2] = 2;	       inv_y[2] = 9;	     inv_z[2] = 2;
	inv_x[3] = 3;	       inv_y[3] = 3;	     inv_z[3] = 10;
	inv_x[4] = 14;	       inv_y[4] = 6;	     inv_z[4] = 5;
	inv_x[5] = 13;	       inv_y[5] = 7;	     inv_z[5] = 4;
	inv_x[6] = 12;	       inv_y[6] = 4;	     inv_z[6] = 7;
	inv_x[7] = 11;	       inv_y[7] = 5;	     inv_z[7] = 6;
	inv_x[8] = 1;	       inv_y[8] = 8;	     inv_z[8] = 8;
	inv_x[9] = 9;	       inv_y[9] = 2;	     inv_z[9] = 9;
	inv_x[10] = 10;        inv_y[10] = 10;	     inv_z[10] = 3;
	inv_x[11] = 7;	       inv_y[11] = 13;	     inv_z[11] = 12;
	inv_x[12] = 6;	       inv_y[12] = 14;	     inv_z[12] = 11;
	inv_x[13] = 5;	       inv_y[13] = 11;	     inv_z[13] = 14;  
	inv_x[14] = 4;	       inv_y[14] = 12;	     inv_z[14] = 13;
#endif

#ifdef GRID_POP_D3Q19
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



	inv_x[0] = 0;	       inv_y[0] = 0;	     inv_z[0] = 0;
	inv_x[1] = 2;	       inv_y[1] = 1;	     inv_z[1] = 1;
	inv_x[2] = 1;	       inv_y[2] = 2;	     inv_z[2] = 2;
	inv_x[3] = 3;	       inv_y[3] = 4;	     inv_z[3] = 3;
	inv_x[4] = 4;	       inv_y[4] = 3;	     inv_z[4] = 4;
	inv_x[5] = 5;	       inv_y[5] = 5;	     inv_z[5] = 6;
	inv_x[6] = 6;	       inv_y[6] = 6;	     inv_z[6] = 5;
	inv_x[7] = 9;	       inv_y[7] = 8;	     inv_z[7] = 7;
	inv_x[8] = 10;	       inv_y[8] = 7;	     inv_z[8] = 8;
	inv_x[9] = 7;	       inv_y[9] = 10;	     inv_z[9] = 9;
	inv_x[10] = 8;        inv_y[10] = 9;	     inv_z[10] = 10;
	inv_x[11] = 12;       inv_y[11] = 11;	     inv_z[11] = 13;
	inv_x[12] = 11;       inv_y[12] = 12;	     inv_z[12] = 14;
	inv_x[13] = 14;       inv_y[13] = 13;	     inv_z[13] = 11;  
	inv_x[14] = 13;       inv_y[14] = 14;	     inv_z[14] = 12;
	inv_x[15] = 15;       inv_y[15] = 17;	     inv_z[15] = 16;
	inv_x[16] = 16;       inv_y[16] = 18;	     inv_z[16] = 15;
	inv_x[17] = 17;       inv_y[17] = 15;	     inv_z[17] = 18;  
	inv_x[18] = 18;       inv_y[18] = 16;	     inv_z[18] = 17;

#endif

#ifdef GRID_POP_D3Q27
	/* D3Q27 */
	wgt[0] = 8. / 27.;
	wgt[1] = 2. / 27.;
	wgt[2] = 2. / 27.;
	wgt[3] = 2. / 27.;
	wgt[4] = 1. / 54.;
	wgt[5] = 1. / 54.;
	wgt[6] = 1. / 54.;
	wgt[7] = 1. / 54.;
	wgt[8] = 1. / 54.;
	wgt[9] = 1. / 54.;
	wgt[10] = 1. / 216.;
	wgt[11] = 1. / 216.;
	wgt[12] = 1. / 216.;
	wgt[13] = 1. / 216.;
	wgt[14] = 2. / 27.;
	wgt[15] = 2. / 27.;
	wgt[16] = 2. / 27.;
	wgt[17] = 1. / 54.;
	wgt[18] = 1. / 54.;
	wgt[19] = 1. / 54.;
	wgt[20] = 1. / 54.;
	wgt[21] = 1. / 54.;
	wgt[22] = 1. / 54.;
	wgt[23] = 1. / 216.;
	wgt[24] = 1. / 216.;
	wgt[25] = 1. / 216.;
	wgt[26] = 1. / 216.;

	/* Lattice speeds, D3Q27 */
	c[0].x = 0.;	c[0].y = 0.;	c[0].z = 0.;
	c[1].x = -1.;	c[1].y = 0.;	c[1].z = 0.;
	c[2].x = 0.;	c[2].y = -1.;	c[2].z = 0.;
	c[3].x = 0.;	c[3].y = 0.;	c[3].z = -1.;
	c[4].x = -1.;	c[4].y = -1.;	c[4].z = 0.;
	c[5].x = -1.;	c[5].y = 1.;	c[5].z = 0.;
	c[6].x = -1.;	c[6].y = 0.;	c[6].z = -1.;
	c[7].x = -1.;	c[7].y = 0.;	c[7].z = 1.;
	c[8].x = 0.;	c[8].y = -1.;	c[8].z = -1.;
	c[9].x = 0.;	c[9].y = -1.;	c[9].z = 1.;
	c[10].x = -1.;	c[10].y = -1.;	c[10].z = -1.;
	c[11].x = -1.;	c[11].y = -1.;	c[11].z = 1.;
	c[12].x = -1.;	c[12].y = 1.;	c[12].z = -1.;
	c[13].x = -1.;	c[13].y = 1.;	c[13].z = 1.;
	c[14].x = 1.;	c[14].y = 0.;	c[14].z = 0.;
	c[15].x = 0.;	c[15].y = 1.;	c[15].z = 0.;
	c[16].x = 0.;	c[16].y = 0.;	c[16].z = 1.;
	c[17].x = 1.;	c[17].y = 1.;	c[17].z = 0.;
	c[18].x = 1.;	c[18].y = -1.;	c[18].z = 0.;
	c[19].x = 1.;	c[19].y = 0.;	c[19].z = 1.;
	c[20].x = 1.;	c[20].y = 0.;	c[20].z = -1.;
	c[21].x = 0.;	c[21].y = 1.;	c[21].z = 1.;
	c[22].x = 0.;	c[22].y = 1.;	c[22].z = -1.;
	c[23].x = 1.;	c[23].y = 1.;	c[23].z = 1.;
	c[24].x = 1.;	c[24].y = 1.;	c[24].z = -1.;
	c[25].x = 1.;	c[25].y = -1.;	c[25].z = 1.;
	c[26].x = 1.;	c[26].y = -1.;	c[26].z = -1.;



	/*
	 * This function provides the opposite velocity c_j=inv[i] for
	 * velocity c_i
	 */
	inv[0] = 0;
	inv[1] = 14;
	inv[2] = 15;
	inv[3] = 16;
	inv[4] = 17;
	inv[5] = 18;
	inv[6] = 19;
	inv[7] = 20;
	inv[8] = 21;
	inv[9] = 22;
	inv[10] = 23;
	inv[11] = 24;
	inv[12] = 25;
	inv[13] = 26;
	inv[14] = 1;
	inv[15] = 2;
	inv[16] = 3;
	inv[17] = 4;
	inv[18] = 5;
	inv[19] = 6;
	inv[20] = 7;
	inv[21] = 8;
	inv[22] = 9;
	inv[23] = 10;
	inv[24] = 11;
	inv[25] = 12;
	inv[26] = 13;



	inv_x[0] = 0;	       inv_y[0] = 0;	     inv_z[0] = 0;
	inv_x[1] = 14;	       inv_y[1] = 1;	     inv_z[1] = 1;
	inv_x[2] = 2;	       inv_y[2] = 15;	     inv_z[2] = 2;
	inv_x[3] = 3;	       inv_y[3] = 3;	     inv_z[3] = 16;
	inv_x[4] = 18;	       inv_y[4] = 5;	     inv_z[4] = 4;
	inv_x[5] = 17;	       inv_y[5] = 4;	     inv_z[5] = 5;
	inv_x[6] = 20;	       inv_y[6] = 6;	     inv_z[6] = 7;
	inv_x[7] = 19;	       inv_y[7] = 7;	     inv_z[7] = 6;
	inv_x[8] = 8;	       inv_y[8] = 22;	     inv_z[8] = 9;
	inv_x[9] = 9;	       inv_y[9] = 21;	     inv_z[9] = 8;
	inv_x[10] = 26;       inv_y[10] = 12;	     inv_z[10] = 11;
	inv_x[11] = 25;       inv_y[11] = 13;	     inv_z[11] = 10;
	inv_x[12] = 24;       inv_y[12] = 10;	     inv_z[12] = 13;
	inv_x[13] = 23;       inv_y[13] = 11;	     inv_z[13] = 12;  
	inv_x[14] = 1;	      inv_y[14] = 14;	     inv_z[14) = 14;
	inv_x[15] = 15;       inv_y[15] = 2;	     inv_z[15) = 15;
	inv_x[16] = 16;       inv_y[16] = 16;	     inv_z[16] = 3;
	inv_x[17] = 5;        inv_y[17] = 18;	     inv_z[17]= 17;
	inv_x[18] = 4;        inv_y[18] = 17;	     inv_z[18] = 18;
	inv_x[19] = 7;        inv_y[19] = 19;	     inv_z[19] = 20;
	inv_x[20] = 6;	      inv_y[20] = 20;	     inv_z[20] = 19;
	inv_x[21] = 21;       inv_y[21] = 9;	     inv_z[21] = 22;
	inv_x[22] = 22;       inv_y[22] = 8;	     inv_z[22] = 21;
	inv_x[23] = 13;       inv_y[23] = 25;	     inv_z[23] = 24;
	inv_x[24] = 12;       inv_y[24] = 26;	     inv_z[24] = 23;
	inv_x[25] = 11;       inv_y[25] = 23;	     inv_z[25] = 26;  
	inv_x[26] = 10;       inv_y[26] = 24;	     inv_z[26] = 25;

#endif



	/* Speed of sound constants */

	/* just for a check: it shall be = 1 */
	/*
	int pp;
	cs2 = 0.0;
	for(pp=0;pp<NPOP;pp++) cs2 += wgt[pp]*(c[pp].x*c[pp].x + c[pp].y*c[pp].y + c[pp].z*c[pp].z);
	fprintf(stderr,"cs2 %e\n",cs2);
	*/
	cs  = 1.0 /sqrt(3.0);  /* Why this value? Is this a property of the D3Q19 lattice or just a choice of convenience?*/
	cs2 = pow((double)cs,2.0); //1.0 / 3.0;
	cs4 = pow((double)cs,4.0); //1.0 / 9.0;
	invcs  = 1.0 / cs;
	invcs2 = 1.0 / cs2;
	invcs4 = 1.0 / cs4;
	twocs2 = 2.0 * cs2;
	twocs4 = 2.0 * cs4;
	invtwocs2 = 1.0 / twocs2;
	invtwocs4 = 1.0 / twocs4;

}


pop equilibrium(pop * f, int i, int j, int k){
	int             pp;
	my_double       ux, uy, uz;
	my_double       rhof;
	my_double       cu, u2;
	pop             f_eq;

#ifdef METHOD_LOG
	rhof = m(f[IDX(i, j, k)],1./property.tau_u);
#else
	rhof = m(f[IDX(i, j, k)]);
#endif
	ux = u[IDX(i, j, k)].x;
	uy = u[IDX(i, j, k)].y;
	uz = u[IDX(i, j, k)].z;				

	u2 = (ux * ux + uy * uy + uz * uz);

	/* equilibrium distribution */
	for (pp = 0; pp < NPOP; pp++) {
		cu = (c[pp].x * ux + c[pp].y * uy + c[pp].z * uz);
		f_eq.p[pp] = rhof * wgt[pp] * (1.0 + invcs2 * cu + invtwocs4 * cu * cu - invtwocs2 * u2);

		
#ifdef METHOD_LOG
		f_eq.p[pp] = property.tau_u*log(f_eq.p[pp]);
#endif
		
	}

	//fprintf(stderr,"i %d j %d k %d\n",i,j,k);

	return f_eq;
}



/**************************************************/
void hydro_fields(char which_pop){
	int i, j, k;
#ifdef DEBUG
	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;
#endif
	vector u0, u0_all;
	my_double norm;

	for (i = BRD; i < LNX+BRD; i++)
		for (j = BRD; j < LNY+BRD; j++)
			for (k = BRD; k < LNZ+BRD; k++) {


#ifdef LB_FLUID
			  if(which_pop == 'p'){
 #ifdef LB_FLUID_PAST
			    /* NOTE: copy element by element is expensive, should be optimized in the future */
			   old_u[IDX(i, j, k)].x = u[IDX(i, j, k)].x;
			   old_u[IDX(i, j, k)].y = u[IDX(i, j, k)].y;
			   old_u[IDX(i, j, k)].z = u[IDX(i, j, k)].z;
			   old_dens[IDX(i, j, k)] = dens[IDX(i, j, k)];
 #endif
#ifdef METHOD_LOG
			  dens[IDX(i, j, k)] = m(p[IDX(i, j, k)],1./property.tau_u);			       
			  u[IDX(i, j, k)].x = mvx(p[IDX(i, j, k)],1./property.tau_u) / dens[IDX(i, j, k)];
			  u[IDX(i, j, k)].y = mvy(p[IDX(i, j, k)],1./property.tau_u) / dens[IDX(i, j, k)];
			  u[IDX(i, j, k)].z = mvz(p[IDX(i, j, k)],1./property.tau_u) / dens[IDX(i, j, k)];
#else
				dens[IDX(i, j, k)] = m(p[IDX(i, j, k)]);

				#ifndef METHOD_FORCING_GUO
				u[IDX(i, j, k)].x = mvx(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
				u[IDX(i, j, k)].y = mvy(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
				u[IDX(i, j, k)].z = mvz(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
				#else
				 #ifdef LB_FLUID_FORCING
				/* to be used only with METHOD_STREAMING */
				/* TO BE USED if the force is a force */
				//u[IDX(i, j, k)].x= ( mvx(p[IDX(i, j, k)]) + 0.5*property.time_dt*force[IDX(i, j, k)].x )/dens[IDX(i, j, k)];
   				//u[IDX(i, j, k)].y= ( mvy(p[IDX(i, j, k)]) + 0.5*property.time_dt*force[IDX(i, j, k)].y )/dens[IDX(i, j, k)];
				//u[IDX(i, j, k)].z= ( mvz(p[IDX(i, j, k)]) + 0.5*property.time_dt*force[IDX(i, j, k)].z )/dens[IDX(i, j, k)];
				/* TO BE USED if the force is a force per unit mass i.e. acceleration */
				u[IDX(i, j, k)].x= mvx(p[IDX(i, j, k)])/dens[IDX(i, j, k)] + 0.5*property.time_dt*force[IDX(i, j, k)].x ;
   				u[IDX(i, j, k)].y= mvy(p[IDX(i, j, k)])/dens[IDX(i, j, k)] + 0.5*property.time_dt*force[IDX(i, j, k)].y ;
				u[IDX(i, j, k)].z= mvz(p[IDX(i, j, k)])/dens[IDX(i, j, k)] + 0.5*property.time_dt*force[IDX(i, j, k)].z ;	
				 #else
				u[IDX(i, j, k)].x = mvx(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
				u[IDX(i, j, k)].y = mvy(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
				u[IDX(i, j, k)].z = mvz(p[IDX(i, j, k)]) / dens[IDX(i, j, k)];
				 #endif
				#endif
#endif

			  }
#endif

#ifdef LB_TEMPERATURE 
			  if(which_pop == 'g'){
 #ifdef LB_TEMPERATURE_PAST
			    /* NOTE: copy element by element is expensive, should be optimized in the future */
			   old_t[IDX(i, j, k)] = t[IDX(i, j, k)];
 #endif
				/* compute temperature field */
                               #ifndef METHOD_FORCING_MALASPINAS
				t[IDX(i, j, k)] = m(g[IDX(i, j, k)]);
                               #else
                       	/* to be used only with METHOD_STREAMING */
                                #ifndef LB_TEMPERATURE_FORCING
	                         t[IDX(i, j, k)] = m(g[IDX(i, j, k)]);
                                #else	                         
				 /* not as in eq (8.40) , pp. 310 of Timm Krueger at al. book ---> there is an error */
				 /* now as in Takeshi Seta, PHYSICAL REVIEW E 87, 063304 (2013) eq (22)  */ 
	                         t[IDX(i, j, k)] = m(g[IDX(i, j, k)]) + 0.5*property.time_dt*t_source[IDX(i, j, k)];
                                #endif
                               #endif
			  }
#endif

#ifdef LB_SCALAR
                          if(which_pop == 'h'){                  
				    /* compute scalar field */
                               #ifndef METHOD_FORCING_MALASPINAS	    
				s[IDX(i, j, k)] = m(h[IDX(i, j, k)]); 
                               #else
                       	/* to be used only with METHOD_STREAMING */
                                #ifndef LB_SCALAR_FORCING
	                         s[IDX(i, j, k)] = m(h[IDX(i, j, k)]);
                                #else	              
				 /* not as in eq (8.40) , pp. 310 of Timm Krueger at al. book ---> there is an error */
				 /* now as in Takeshi Seta, PHYSICAL REVIEW E 87, 063304 (2013) eq (22)  */ 			        
				 s[IDX(i, j, k)] = m(h[IDX(i, j, k)]) + 0.5*property.time_dt*s_source[IDX(i, j, k)]; 
                                #endif
                               #endif
			  }                           
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
		   fprintf(fout, "%d %d %d %e %e %e %e %e %e %e\n", i, j, k, center_V[IDX(i, j, k)].x, center_V[IDX(i, j, k)].y, center_V[IDX(i, j, k)].z, 
		  	                                           u[IDX(i, j, k)].x, u[IDX(i, j, k)].y, u[IDX(i, j, k)].z , dens[IDX(i, j, k)]);
	fclose(fout);
#endif

}
/* end of hydro_fields */


#ifdef METHOD_STREAMING
void streaming(pop *f, pop *rhs_f,int i,int j,int k){
//void streaming(pop *f, pop *rhs_f){
// int i, j, k;
 int pp;
 int ii, jj, kk;

 
 //  for(k=BRD;k<LNZ+BRD;k++){
 //    for(j=BRD;j<LNY+BRD;j++){
 //       for(i=BRD;i<LNX+BRD;i++){ 

     for(pp=0;pp<NPOP;pp++){
       
       ii = i-(int)c[pp].x;
       jj = j-(int)c[pp].y;
       kk = k-(int)c[pp].z;
       //f[IDX(i,j,k)].p[pp]  = rhs_f[IDX(ii,jj,kk)].p[pp];            
#ifdef LB_FLUID_FORCING_LANDSCAPE
       /* Here we avoid streaming for solid regions, i.e., for landscape = 1 */
       if(landscape[IDX(i, j, k)] == 0.0) f[IDX(i,j,k)].p[pp]  = rhs_f[IDX(ii,jj,kk)].p[pp];
#else
       /* the normal streaming case : pull method */ 
       if(property.time_dt == 1.0){ 
	 /* when the time step is = 1 */
       f[IDX(i,j,k)].p[pp]  = rhs_f[IDX(ii,jj,kk)].p[pp];
       }else{       
	 /* when the time step is not 1 */
       f[IDX(i,j,k)].p[pp]  = property.time_dt*rhs_f[IDX(ii,jj,kk)].p[pp];
       }
#endif

// if(j==50) fprintf(stderr,"%d %d %d pp %d %e %e\n", ii , jj , kk, pp, f[IDX(i,j,k)].p[pp], m(f[IDX(i,j,k)]) );


     }/* pp */

     //           }/* i */
     //        }/* j */
     //      }/* k */
}
#endif


/********************************************/
void time_stepping(pop *f, pop *rhs_f, pop *old_rhs_f, pop *old_old_rhs_f,my_double tau,pop *f_eq,char which_pop){
  int i,j,k,pp;
  //pop f_eq;
  my_double dt_over_tau,fac1,fac2;
  my_double rho1,rho2;
  my_double A1, A2, A3,dt,invtau;

  /* useful constants for time advancing */
#ifdef METHOD_EXPONENTIAL
  dt = property.time_dt;
  invtau = 1.0/tau;
  dt_over_tau = property.time_dt/tau;
  fac1 = exp(-dt_over_tau);
  fac2 =1.0/(1.0 + dt_over_tau); 
#endif

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

#ifdef DEBUG_HARD
	//f_eq=equilibrium(f,i,j,k);
#endif
      

	for(pp=0;pp<NPOP;pp++){
#ifdef METHOD_STEPPING_EULER
	  /* Euler first order */
#ifdef METHOD_EXPONENTIAL
	  /* my recipe a.k.a.  Integrating factor method (Lawson 1967) */
	  // f[IDX(i,j,k)].p[pp] =  fac1*(f[IDX(i,j,k)].p[pp] + property.time_dt*rhs_f[IDX(i,j,k)].p[pp]);

	  /* Exponential Euler Method a.k.a. Exponential time differencing (Certaine 1960)*/
	     f[IDX(i,j,k)].p[pp] =  fac1*f[IDX(i,j,k)].p[pp] + (1.0-fac1)*tau*rhs_f[IDX(i,j,k)].p[pp];

	  /* IMEX’ (Implicit-Explicit) or semi-implicit  Ascher, Ruuth & Wetton (1995) */
	  //f[IDX(i,j,k)].p[pp] =  fac2*(f[IDX(i,j,k)].p[pp] + property.time_dt*rhs_f[IDX(i,j,k)].p[pp]);
#else
          f[IDX(i,j,k)].p[pp] += property.time_dt*rhs_f[IDX(i,j,k)].p[pp];	
#endif

#ifdef DEBUG_HARD
	  //fprintf(stderr,"f_new %e f_eq %e diff %e\n",f[IDX(i,j,k)].p[pp], f_eq.p[pp], f[IDX(i,j,k)].p[pp]-f_eq.p[pp]);
#endif 
#endif

#ifdef METHOD_STEPPING_AB2
	  /* Adams Bashforth 2nd order */
	  //fprintf(stderr,"itime = %d\n",itime);
#ifdef METHOD_EXPONENTIAL
	  /* my AB2 exponential */
	  //if(itime==1)  
	  //f[IDX(i,j,k)].p[pp] =  fac1*(f[IDX(i,j,k)].p[pp] + property.time_dt*rhs_f[IDX(i,j,k)].p[pp]);
	  //else
	  //f[IDX(i,j,k)].p[pp] = fac1*(f[IDX(i,j,k)].p[pp] +  property.time_dt*(1.5*rhs_f[IDX(i,j,k)].p[pp]-0.5*fac1*old_rhs_f[IDX(i,j,k)].p[pp])); 

	  /* S. M. Cox and P. C. Matthews. Exponential Time Differencing for Stiff Systems.
	     J. Comput. Phys., 176:430{455, 2002. */	  
	  if(itime==1)   f[IDX(i,j,k)].p[pp] =  fac1*f[IDX(i,j,k)].p[pp] + (1.0-fac1)*tau*rhs_f[IDX(i,j,k)].p[pp];
	  else{
	    f[IDX(i,j,k)].p[pp] *= fac1;
	    A1 = ((fac1-1.0)*(1.0-dt_over_tau)+dt_over_tau)*(tau/dt_over_tau);
	    A2 = (-fac1+1.0-dt_over_tau)*(tau/dt_over_tau);  
	    f[IDX(i,j,k)].p[pp] += ( rhs_f[IDX(i,j,k)].p[pp]*A1 + old_rhs_f[IDX(i,j,k)].p[pp]*A2 );
	      // was working fine:
	      //f[IDX(i,j,k)].p[pp] = f[IDX(i,j,k)].p[pp]*fac1 + ( rhs_f[IDX(i,j,k)].p[pp]*((1.0-dt_over_tau)*fac1-1.0+2.0*dt_over_tau) 
	      //+ old_rhs_f[IDX(i,j,k)].p[pp]*(-fac1+1.0-dt_over_tau) )*(tau/dt_over_tau);
	  }
#else
	 if(itime==1) f[IDX(i,j,k)].p[pp] += property.time_dt*rhs_f[IDX(i,j,k)].p[pp];
	 else
	 f[IDX(i,j,k)].p[pp] += property.time_dt*(1.5*rhs_f[IDX(i,j,k)].p[pp]-0.5*old_rhs_f[IDX(i,j,k)].p[pp]); 
#endif

	 old_rhs_f[IDX(i,j,k)].p[pp] = rhs_f[IDX(i,j,k)].p[pp];
#endif

	 /* 3rd order */
#ifdef METHOD_STEPPING_AB3
#ifdef METHOD_EXPONENTIAL
	  /* S. M. Cox and P. C. Matthews. Exponential Time Differencing for Stiff Systems J. Comput. Phys., 176:430{455, 2002. */	  
	  if(itime==1) f[IDX(i,j,k)].p[pp] =  fac1*f[IDX(i,j,k)].p[pp] + (1.0-fac1)*tau*rhs_f[IDX(i,j,k)].p[pp];
	  if(itime==2)
	  f[IDX(i,j,k)].p[pp] = f[IDX(i,j,k)].p[pp]*fac1 + ( rhs_f[IDX(i,j,k)].p[pp]*((1.0-dt_over_tau)*fac1-1.0+2.0*dt_over_tau) 
	  						     + old_rhs_f[IDX(i,j,k)].p[pp]*(-fac1+1.0-dt_over_tau) )*(tau/dt_over_tau);	 
	  if(itime > 2){
	    A1 = -( (-2.0 + (2.0 - dt_over_tau*(3.0 - 2.0*dt_over_tau) )*exp(-dt_over_tau) + dt_over_tau*(5.0 - 6.0*dt_over_tau) )*pow(tau,3.0))/(2.0*pow(dt,2.0));
	    A2 = -( ( 4.0 - (4.0 * (1.0 - dt_over_tau) ) / exp(dt_over_tau) -  dt_over_tau*( 8.0 - 6.0*dt_over_tau)  )*pow(tau,3.0))/(2.0*pow(dt,2.0));
	    A3 = -( (-2.0 + (2.0 - dt_over_tau)/exp(dt_over_tau) + dt_over_tau*(3.0 - 2.0*dt_over_tau) )*pow(tau,3.0))/ (2.0*pow(dt,2.0));

          f[IDX(i,j,k)].p[pp] = f[IDX(i,j,k)].p[pp]*fac1 + A1*rhs_f[IDX(i,j,k)].p[pp] + A2*old_rhs_f[IDX(i,j,k)].p[pp] + A3*old_old_rhs_f[IDX(i,j,k)].p[pp];
	  }
#else 
	  /* here normal 3rd order AB */
	 if(itime==1) f[IDX(i,j,k)].p[pp] += property.time_dt*rhs_f[IDX(i,j,k)].p[pp];
	 if(itime==2) f[IDX(i,j,k)].p[pp] += property.time_dt*(1.5*rhs_f[IDX(i,j,k)].p[pp]-0.5*old_rhs_f[IDX(i,j,k)].p[pp]);
	 if(itime>2)  f[IDX(i,j,k)].p[pp] += (property.time_dt/12.)*(23.*rhs_f[IDX(i,j,k)].p[pp]-16.*old_rhs_f[IDX(i,j,k)].p[pp] + 5.0*old_old_rhs_f[IDX(i,j,k)].p[pp]);
#endif
	 old_old_rhs_f[IDX(i,j,k)].p[pp] = old_rhs_f[IDX(i,j,k)].p[pp];
	 old_rhs_f[IDX(i,j,k)].p[pp] = rhs_f[IDX(i,j,k)].p[pp];
#endif

	}/* pp */

#ifdef METHOD_STREAMING
	streaming(f, rhs_f, i , j , k);
#endif
      }/* for i,j,k */


#ifdef METHOD_HEUN
  /* this is Heun corrector scheme method, it depends on method Euler (METHOD_EULER must be defined) */

  /* 0) adjust bc for the new f field*/
 sendrecv_borders_pop(f);
#if (defined LB_FLUID_BC || defined LB_TEMPERATURE_BC || defined LB_SCALAR_BC)
  boundary_conditions();
#endif
  /* 1) we compute the advection with the new f field*/
     /* 1a) first we need the new hydro fields */
     //if(which_pop == 'p') hydro_fields('p');
     /* 1b) then we need to build the forcing */ 	
         //build_forcing();
     /* 1c) finally the new the advection */ 	
  compute_advection(f,rhs_f,tau,f_eq,which_pop);
  /* 2) perform the corrector step */
  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++)
	for(pp=0;pp<NPOP;pp++){

	  f[IDX(i,j,k)].p[pp] += 0.5*property.time_dt*( rhs_f[IDX(i,j,k)].p[pp] - old_rhs_f[IDX(i,j,k)].p[pp] );	  

	}
#endif
  /*
  #ifdef METHOD_STREAMING
      streaming(f, rhs_f);
  #endif
  */
}



/********************************************/ 
void sendrecv_borders_pop(pop *f){
  int i,j,k,brd_size;
  MPI_Status status1;

#ifdef DEBUG_HARD
 fprintf(stderr,"before\n");
 for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++)
	fprintf(stderr,"%d %d %d %e\n", i , j , k, m(f[IDX(i,j,k)]));
#endif 

#ifdef METHOD_MYQUICK
  /* if the method is Quick, then  BRD=2 : we have to be careful  when the NX,NY or NZ values are =1 */
  /* in that case we copy the same value everywhere */
  if(NX==1){
  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
     f[IDX(0,j,k)] = f[IDX(1,j,k)] = f[IDX(3,j,k)] = f[IDX(4,j,k)] =  f[IDX(2,j,k)];
  } 

  if(NY==1){
  for(k=BRD;k<LNZ+BRD;k++)
    for(i=BRD;i<LNY+BRD;i++)
     f[IDX(i,0,k)] = f[IDX(i,1,k)] = f[IDX(i,3,k)] = f[IDX(i,4,k)] =  f[IDX(i,2,k)];
  } 

  if(NZ==1){
  for(i=BRD;i<LNX+BRD;i++)
    for(j=BRD;j<LNY+BRD;j++)
     f[IDX(i,j,0)] = f[IDX(i,j,1)] = f[IDX(i,j,3)] = f[IDX(i,j,4)] =  f[IDX(i,j,2)];
  } 
#endif


#ifdef NEW_SENDRECV
  //IDX(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )
     
  MPI_Sendrecv( f + IDX(LNX,0,0)   , 1, MPI_pop_plane_x, me_xp, 14,
                f + IDX(0,0,0)     , 1, MPI_pop_plane_x, me_xm, 14, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(BRD,0,0)    , 1, MPI_pop_plane_x, me_xm, 15,
                f + IDX(LNX+BRD,0,0), 1, MPI_pop_plane_x, me_xp, 15, MPI_COMM_WORLD, &status1);
  
  MPI_Sendrecv( f + IDX(0,LNY,0)  , 1, MPI_pop_plane_y, me_yp, 12,
                f + IDX(0,0,0)  , 1, MPI_pop_plane_y, me_ym, 12, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(0,BRD,0)  , 1, MPI_pop_plane_y, me_ym, 13,
                f + IDX(0,LNY+BRD,0)  , 1, MPI_pop_plane_y, me_yp, 13, MPI_COMM_WORLD, &status1);
  
  MPI_Sendrecv( f + IDX(0,0,LNZ)  , 1, MPI_pop_plane_z, me_zp, 10,
                f + IDX(0,0,0)  , 1, MPI_pop_plane_z, me_zm, 10, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(0,0,BRD)  , 1, MPI_pop_plane_z, me_zm, 11,
                f + IDX(0,0,LNZ+BRD)  , 1, MPI_pop_plane_z, me_zp, 11, MPI_COMM_WORLD, &status1);
  
  #else /* = NEW_SENDRECV is not defined */


  /*     BRD|LNX|BRD     */
  /* Copy borders along x */
   
  brd_size = BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_pop[IDX_XBRD(i,j,k)] = f[IDX(i+LNX,j,k)];
      }

  MPI_Sendrecv( xp_pop, brd_size, MPI_pop_type, me_xp, 10,
                xm_pop, brd_size, MPI_pop_type, me_xm, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_pop[IDX_XBRD(i,j,k)];
        xm_pop[IDX_XBRD(i,j,k)] = f[IDX(i+BRD,j,k)];
      }
 MPI_Sendrecv( xm_pop, brd_size, MPI_pop_type, me_xm, 11,
               xp_pop, brd_size, MPI_pop_type, me_xp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k)] = xp_pop[IDX_XBRD(i,j,k)];
      }
  


  /* Copy borders along y */
 
  brd_size = BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_pop[IDX_YBRD(i,j,k)] = f[IDX(i,j+LNY,k)];
	//	fprintf(stderr,"\n SEND %d  %d %d %e %d %e\n",i,j,k,yp_pop[IDX_YBRD(i,j,k)].p[7],j+LNY,f[IDX(i,j+LNY,k)].p[7] );
      }

  MPI_Sendrecv( yp_pop, brd_size, MPI_pop_type, me_yp, 10,
                ym_pop, brd_size, MPI_pop_type, me_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = ym_pop[IDX_YBRD(i,j,k)];
	//fprintf(stderr,"\n RECV %d  %d %d %e %d %e\n",i,j,k,ym_pop[IDX_YBRD(i,j,k)].p[7],j+LNY,f[IDX(i,j,k)].p[7] );
	ym_pop[IDX_YBRD(i,j,k)] = f[IDX(i,j+BRD,k)];
      }
  
 MPI_Sendrecv( ym_pop, brd_size, MPI_pop_type, me_ym, 13,
               yp_pop, brd_size, MPI_pop_type, me_yp, 13, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
       	f[IDX(i,j+LNY+BRD,k)] = yp_pop[IDX_YBRD(i,j,k)];
      }
  

  /* Copy borders along z */
 
  brd_size = BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        zp_pop[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+LNZ)];
      }

  MPI_Sendrecv( zp_pop, brd_size, MPI_pop_type, me_zp, 14,
                zm_pop, brd_size, MPI_pop_type, me_zm, 14, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = zm_pop[IDX_ZBRD(i,j,k)];
        zm_pop[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+BRD)];
      }
 MPI_Sendrecv( zm_pop, brd_size, MPI_pop_type, me_zm, 15,
               zp_pop, brd_size, MPI_pop_type, me_zp, 15, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j,k+LNZ+BRD)] = zp_pop[IDX_ZBRD(i,j,k)];
      }
 
#endif /* NEW_SENDRECV */

#ifdef METHOD_EDGES_AND_CORNERS

  /* First we communicate the 8 corner cubes (they are either 1x1x1 or 2x2x2 depending on BRD) */ 
  
  brd_size = BRD*BRD*BRD;

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_zp_corner_pop[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+LNZ)];
        xp_yp_zm_corner_pop[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+BRD)];
        xp_ym_zp_corner_pop[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+LNZ)];
        xm_yp_zp_corner_pop[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+LNZ)];
      }

  MPI_Sendrecv( xp_yp_zp_corner_pop, brd_size, MPI_pop_type, me_xp_yp_zp, 10,
                xm_ym_zm_corner_pop, brd_size, MPI_pop_type, me_xm_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_yp_zm_corner_pop, brd_size, MPI_pop_type, me_xp_yp_zm, 10,
                xm_ym_zp_corner_pop, brd_size, MPI_pop_type, me_xm_ym_zp, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_ym_zp_corner_pop, brd_size, MPI_pop_type, me_xp_ym_zp, 10,
                xm_yp_zm_corner_pop, brd_size, MPI_pop_type, me_xm_yp_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_zp_corner_pop, brd_size, MPI_pop_type, me_xm_yp_zp, 10,
                xp_ym_zm_corner_pop, brd_size, MPI_pop_type, me_xp_ym_zm, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_ym_zm_corner_pop[IDX_CORNER(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = xm_ym_zp_corner_pop[IDX_CORNER(i,j,k)];
        f[IDX(i,j+LNY+BRD,k)] = xm_yp_zm_corner_pop[IDX_CORNER(i,j,k)];
        f[IDX(i+LNX+BRD,j,k)] = xp_ym_zm_corner_pop[IDX_CORNER(i,j,k)];
        xm_ym_zm_corner_pop[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+BRD)];
        xm_ym_zp_corner_pop[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+LNZ)];
        xm_yp_zm_corner_pop[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+BRD)];
        xp_ym_zm_corner_pop[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+BRD)];
      }
 MPI_Sendrecv( xm_ym_zm_corner_pop, brd_size, MPI_pop_type, me_xm_ym_zm, 11,
               xp_yp_zp_corner_pop, brd_size, MPI_pop_type, me_xp_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_ym_zp_corner_pop, brd_size, MPI_pop_type, me_xm_ym_zp, 11,
               xp_yp_zm_corner_pop, brd_size, MPI_pop_type, me_xp_yp_zm, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_yp_zm_corner_pop, brd_size, MPI_pop_type, me_xm_yp_zm, 11,
               xp_ym_zp_corner_pop, brd_size, MPI_pop_type, me_xp_ym_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_zm_corner_pop, brd_size, MPI_pop_type, me_xp_ym_zm, 11,
               xm_yp_zp_corner_pop, brd_size, MPI_pop_type, me_xm_yp_zp, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k+LNZ+BRD)] = xp_yp_zp_corner_pop[IDX_CORNER(i,j,k)];
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)] = xp_yp_zm_corner_pop[IDX_CORNER(i,j,k)];
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)] = xp_ym_zp_corner_pop[IDX_CORNER(i,j,k)];
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)] = xm_yp_zp_corner_pop[IDX_CORNER(i,j,k)];
      }
 

 /* Then we communicate the 12 edges  */
 
 /* along x */
 
 brd_size = BRD*BRD*(LNX+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_zp_edge_pop[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+LNZ)];
        yp_zm_edge_pop[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+BRD)];
      }

  MPI_Sendrecv( yp_zp_edge_pop, brd_size, MPI_pop_type, me_yp_zp, 10,
                ym_zm_edge_pop, brd_size, MPI_pop_type, me_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( yp_zm_edge_pop, brd_size, MPI_pop_type, me_yp_zm, 10,
                ym_zp_edge_pop, brd_size, MPI_pop_type, me_ym_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = ym_zm_edge_pop[IDX_EDGE_X(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = ym_zp_edge_pop[IDX_EDGE_X(i,j,k)];
        ym_zm_edge_pop[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+BRD)];
        ym_zp_edge_pop[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+LNZ)];
      }
 MPI_Sendrecv( ym_zm_edge_pop, brd_size, MPI_pop_type, me_ym_zm, 11,
               yp_zp_edge_pop, brd_size, MPI_pop_type, me_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( ym_zp_edge_pop, brd_size, MPI_pop_type, me_ym_zp, 11,
               yp_zm_edge_pop, brd_size, MPI_pop_type, me_yp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)] = yp_zp_edge_pop[IDX_EDGE_X(i,j,k)];
	f[IDX(i,j+LNY+BRD,k)] = yp_zm_edge_pop[IDX_EDGE_X(i,j,k)];
      }
 
 /* along y */
 
 brd_size = BRD*BRD*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_zp_edge_pop[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+LNZ)];
        xp_zm_edge_pop[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+BRD)];
      }

  MPI_Sendrecv( xp_zp_edge_pop, brd_size, MPI_pop_type, me_xp_zp, 10,
                xm_zm_edge_pop, brd_size, MPI_pop_type, me_xm_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_zm_edge_pop, brd_size, MPI_pop_type, me_xp_zm, 10,
                xm_zp_edge_pop, brd_size, MPI_pop_type, me_xm_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_zm_edge_pop[IDX_EDGE_Y(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = xm_zp_edge_pop[IDX_EDGE_Y(i,j,k)];
        xm_zm_edge_pop[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+BRD)];
        xm_zp_edge_pop[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+LNZ)];
      }
 MPI_Sendrecv( xm_zm_edge_pop, brd_size, MPI_pop_type, me_xm_zm, 11,
               xp_zp_edge_pop, brd_size, MPI_pop_type, me_xp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_zp_edge_pop, brd_size, MPI_pop_type, me_xm_zp, 11,
               xp_zm_edge_pop, brd_size, MPI_pop_type, me_xp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)] = xp_zp_edge_pop[IDX_EDGE_Y(i,j,k)];
	f[IDX(i+LNX+BRD,j,k)] = xp_zm_edge_pop[IDX_EDGE_Y(i,j,k)];
      }
 
 
 /* along z */
 
 brd_size = BRD*BRD*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_edge_pop[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+LNY,k)];
        xm_yp_edge_pop[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+LNY,k)];
      }

  MPI_Sendrecv( xp_yp_edge_pop, brd_size, MPI_pop_type, me_xp_yp, 10,
                xm_ym_edge_pop, brd_size, MPI_pop_type, me_xm_ym, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_edge_pop, brd_size, MPI_pop_type, me_xm_yp, 10,
                xp_ym_edge_pop, brd_size, MPI_pop_type, me_xp_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_ym_edge_pop[IDX_EDGE_Z(i,j,k)];
        f[IDX(i+LNX+BRD,j,k)] = xp_ym_edge_pop[IDX_EDGE_Z(i,j,k)];
        xm_ym_edge_pop[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+BRD,k)];
        xp_ym_edge_pop[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+BRD,k)];
      }
 MPI_Sendrecv( xm_ym_edge_pop, brd_size, MPI_pop_type, me_xm_ym, 11,
               xp_yp_edge_pop, brd_size, MPI_pop_type, me_xp_yp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_edge_pop, brd_size, MPI_pop_type, me_xp_ym, 11,
               xm_yp_edge_pop, brd_size, MPI_pop_type, me_xm_yp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)] = xp_yp_edge_pop[IDX_EDGE_Z(i,j,k)];
	f[IDX(i,j+LNY+BRD,k)] = xm_yp_edge_pop[IDX_EDGE_Z(i,j,k)];
      }
 
#endif

#ifdef DEBUG_HARD
 fprintf(stderr,"after\n");
 for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++)
	fprintf(stderr,"%d %d %d %e\n", i , j , k, m(f[IDX(i,j,k)]));
#endif
       
}/* end send rcv */

/*************************************************************/


pop equilibrium_given_velocity(vector v , my_double rho){
	int             pp;
	my_double       ux, uy, uz;
	my_double       rhof;
	my_double       cu, u2;
	pop             f_eq;

	ux = v.x;
	uy = v.y;
	uz = v.z;
	u2 = (ux * ux + uy * uy + uz * uz);

	/* equilibrium distribution */
	for (pp = 0; pp < NPOP; pp++) {
		cu = (c[pp].x * ux + c[pp].y * uy + c[pp].z * uz);
		f_eq.p[pp] = rho * wgt[pp] * (1.0 + invcs2 * cu + invtwocs4 * cu * cu - invtwocs2 * u2);		 
	}

	return f_eq;
}


/* function to copy a pop field */
void copy_pop(pop *f, pop *f_copy){
  int i,j,k;
  for(k=0;k<LNZ+TWO_BRD;k++)
    for(j=0;j<LNY+TWO_BRD;j++)
      for(i=0;i<LNX+TWO_BRD;i++){

	f_copy[IDX(i,j,k)] = f[IDX(i,j,k)]; 	

	}
}

/*******************/
/*  END of lb.c    */
/*******************/
