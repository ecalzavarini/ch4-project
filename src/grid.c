#include "common_object.h"

void read_mesh(){

	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;
	int             i, j, k;

	sprintf(fnamein, "mesh.in");
	fin = fopen(fnamein, "r");
	if (fin != NULL) {
		fprintf(stderr, "Mesh file %s has been found!\n", fnamein);
	} else {

		fprintf(stderr, "Warning message -> %s file is missing!\n Starting from grid generated on the fly\n ", fnamein);


		/* set field to zero */
		for (k =0; k < LNZG+TWO_BRD; k++)
			for (j =0; j < LNYG+TWO_BRD; j++)
				for (i = 0; i < LNXG+TWO_BRD; i++) {
				  mesh[IDXG(i, j, k)].x = 0.0;
				  mesh[IDXG(i, j, k)].y = 0.0;
				  mesh[IDXG(i, j, k)].z = 0.0;			 
				  mesh_flag[IDXG(i, j, k)] = 0;
				}

		/* moving on the bulk only */
		for (k = BRD; k < LNZG+BRD; k++)
			for (j = BRD; j < LNYG+BRD; j++)
				for (i = BRD; i < LNXG+BRD; i++) {
					mesh[IDXG(i, j, k)].x = (my_double) (i + LNXG_START-BRD);
					mesh[IDXG(i, j, k)].y = (my_double) (j + LNYG_START-BRD);
					mesh[IDXG(i, j, k)].z = (my_double) (k + LNZG_START-BRD);
					/*
					 * flag: 1 is bulk , 0 is wall , -1
					 * is dormient
					 */
					mesh_flag[IDXG(i, j, k)] = 1;
								
#ifdef GRID_RANDOM
					  if(i<LNXG+BRD-1) mesh[IDXG(i, j, k)].x += 0.25*(my_double)(2.0*drand48()-1.0);
					  if(j<LNYG+BRD-1) mesh[IDXG(i, j, k)].y += 0.25*(my_double)(2.0*drand48()-1.0);
					  if(j<LNZG+BRD-1) mesh[IDXG(i, j, k)].z += 0.25*(my_double)(2.0*drand48()-1.0);
#endif

				} /* for ijk */

	} /* end else */


	/* here we copy the borders from the neighbors */
	//	sendrecv_borders_mesh(mesh , mesh_flag);

#ifdef DEBUG
	/* Each processor prints its mesh */
	sprintf(fnamein, "mesh.%d.out", me);
	fout = fopen(fnamein, "w");
	
	for (k = BRD; k < LNZG+BRD; k++){
	  for (j = BRD; j < LNYG+BRD; j++){
	    for (i = BRD; i < LNXG+BRD; i++){ 
				  /*			  
	for (k = 0; k < LNZG + TWO_BRD; k++)
	  for (j = 0; j < LNYG + TWO_BRD; j++)
	    for (i = 0; i < LNXG + TWO_BRD; i++)
				  */
       fprintf(fout, "%d %d %d %e %e %e %d\n", i, j, k, mesh[IDXG(i, j, k)].x, mesh[IDXG(i, j, k)].y, mesh[IDXG(i, j, k)].z , mesh_flag[IDXG(i, j, k)]);
	    }
	    fprintf(fout,"\n");
	  }
	  fprintf(fout,"\n");
	} /* ijk */
	fclose(fout);

	/*  I rewrite it swapping x with y for gnuplot use */
	/* Each processor prints its mesh */
	sprintf(fnamein, "mesh.%d.swap.out", me);
	fout = fopen(fnamein, "w");
	
	for (k = BRD; k < LNZG+BRD; k++){
	    for (i = BRD; i < LNXG+BRD; i++){ 
	      for (j = BRD; j < LNYG+BRD; j++){
				  /*			  
	for (k = 0; k < LNZG + TWO_BRD; k++)
	  for (j = 0; j < LNYG + TWO_BRD; j++)
	    for (i = 0; i < LNXG + TWO_BRD; i++)
				  */
       fprintf(fout, "%d %d %d %e %e %e %d\n", i, j, k, mesh[IDXG(i, j, k)].x, mesh[IDXG(i, j, k)].y, mesh[IDXG(i, j, k)].z , mesh_flag[IDXG(i, j, k)]);
	    }
	    fprintf(fout,"\n");
	  }
	  fprintf(fout,"\n");
	} /* ijk */
	fclose(fout);

#endif

}
/*****************************/

vector vector_product(vector a , vector b){
  vector c;
  c.x =  a.y * b.z - a.z * b.y;
  c.y = -a.x * b.z + a.z * b.x;
  c.z =  a.x * b.y - a.y * b.x;
  return c;
}

my_double scalar_product(vector a , vector b){
  my_double c;
  c =  a.x * b.x + a.y * b.y + a.z * b.z;
  return c;
}

vector vector_difference(vector a , vector b){
  vector c;
  c.x =  a.x - b.x;
  c.y =  a.y - b.y;
  c.z =  a.z - b.z;
  return c;
}

vector vector_sum(vector a , vector b){
  vector c;
  c.x =  a.x + b.x;
  c.y =  a.y + b.y;
  c.z =  a.z + b.z;
  return c;
}


/**************************************************************************************************/

void compute_volumes(){
  int i,j,k,n,pp;

	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;

	/* Allocating array for storing information of control volume */
	vector       P0, P1, P2, P3, P4, P5, P6, P7, P8;
	vector       D17, D35, D03, D12, D05, D14,D06, D24, D27, D36, D47, D56;
	my_double       S1357, S0145, S0246, S0123, S4567, S2367;
	vector       N1357, N0145, N0246, N0123, N4567,N2367;
	vector       M0167, M0257, M0347;
	my_double       V, V1, V2, V3;


	for (i = BRD; i < LNXG + BRD - 1 ; i++) 
		for (j = BRD; j < LNYG + BRD - 1; j++) 
			for (k = BRD; k < LNZG + BRD - 1; k++) {


/* points definition */
/*    
          z   
          ^ 
          |                         4 -----6
          |                        /|    / |
          | ------>y              / |   /  | 
         /                       5-----7   |
        /                        |  0--|---2 
       /                         | /   |  /
      /                          |/    |/
     x                           1-----3

 */

				P0 = mesh[IDXG(i, j, k)];
				P1 = mesh[IDXG(i + 1, j, k)];
				P2 = mesh[IDXG(i, j + 1, k)];
				P3 = mesh[IDXG(i + 1, j + 1, k)];
				P4 = mesh[IDXG(i, j, k + 1)];
				P5 = mesh[IDXG(i + 1, j, k + 1)];
				P6 = mesh[IDXG(i, j + 1, k + 1)];
				P7 = mesh[IDXG(i + 1, j + 1, k + 1)];
				P8.x = (P0.x + P1.x + P2.x + P3.x + P4.x + P5.x + P6.x + P7.x) / 8;
				P8.y = (P0.y + P1.y + P2.y + P3.y + P4.y + P5.y + P6.y + P7.y) / 8;
				P8.z = (P0.z + P1.z + P2.z + P3.z + P4.z + P5.z + P6.z + P7.z) / 8;
#ifdef DEBUG
	fprintf(stdout, "%e %e %e %e %e %e %e %e %e \n", P0.x, P1.x, P2.x, P3.x, P4.x, P5.x, P6.x, P7.x, P8.x);
	fprintf(stdout, "%e %e %e %e %e %e %e %e %e \n", P0.y, P1.y, P2.y, P3.y, P4.y, P5.y, P6.y, P7.y, P8.y);
	fprintf(stdout, "%e %e %e %e %e %e %e %e %e \n", P0.z, P1.z, P2.z, P3.z, P4.z, P5.z, P6.z, P7.z, P8.z);
#endif
				/*
				 * diagonals definition of each surface of a
				 * polyhedron treating as a vector
				 */
	D17 = vector_difference(P7 , P1);
	D35 = vector_difference(P5 , P3);
	D03 = vector_difference(P3 , P0);
	D12 = vector_difference(P2 , P1);
	D05 = vector_difference(P5 , P0);
	D14 = vector_difference(P4 , P1);
	D06 = vector_difference(P6 , P0);
	D24 = vector_difference(P4 , P2);
	D27 = vector_difference(P7 , P2);
	D36 = vector_difference(P6 , P3);
	D47 = vector_difference(P7 , P4);
	D56 = vector_difference(P6 , P5);
#ifdef DEBUG
	    	fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e \n", D17.x, D35.x, D03.x, D12.x, D05.x, D14.x, D06.x, D24.x, D27.x, D36.x, D47.x, D56.x);
	    	fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e \n", D17.y, D35.y, D03.y, D12.y, D05.y, D14.y, D06.y, D24.y, D27.y, D36.y, D47.y, D56.y);
	    	fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e \n", D17.z, D35.z, D03.z, D12.z, D05.z, D14.z, D06.z, D24.z, D27.z, D36.z, D47.z, D56.z);
#endif
				


				/*
				 * SURFACE AREA definition  Vector formulas
				 * for a quadrilateral: The area of a
				 * quadrilateral ABCD can be calculated using
				 * vectors. Let vectors AC and BD form the
				 * diagonals from A to C and from B to D. The
				 * area of the quadrilateral is then S = 0.5 *
				 * mod(AC X BD) and the magnitude of AC X BD
				 * is calculated by finding the determinant
				 * of matrix formed
				 */
				S1357 = 0.5 * fabs(sqrt((pow((D17.y * D35.z - D17.z * D35.y), 2.0)) + (pow((D17.z * D35.x - D17.x * D35.z), 2.0)) + (pow((D17.x * D35.y - D17.y * D35.x), 2.0))));
				S0145 = 0.5 * fabs(sqrt((pow((D05.y * D14.z - D05.z * D14.y), 2.0)) + (pow((D05.z * D14.x - D05.x * D14.z), 2.0)) + (pow((D05.x * D14.y - D05.y * D14.x), 2.0))));
				S0246 = 0.5 * fabs(sqrt((pow((D06.y * D24.z - D06.z * D24.y), 2.0)) + (pow((D06.z * D24.x - D06.x * D24.z), 2.0)) + (pow((D06.x * D24.y - D06.y * D24.x), 2.0))));
				S2367 = 0.5 * fabs(sqrt((pow((D27.y * D36.z - D27.z * D36.y), 2.0)) + (pow((D27.z * D36.x - D27.x * D36.z), 2.0)) + (pow((D27.x * D36.y - D27.y * D36.x), 2.0))));
				S0123 = 0.5 * fabs(sqrt((pow((D03.y * D12.z - D03.z * D12.y), 2.0)) + (pow((D03.z * D12.x - D03.x * D12.z), 2.0)) + (pow((D03.x * D12.y - D03.y * D12.x), 2.0))));
				S4567 = 0.5 * fabs(sqrt((pow((D47.y * D56.z - D47.z * D56.y), 2.0)) + (pow((D47.z * D56.x - D47.x * D56.z), 2.0)) + (pow((D47.x * D56.y - D47.y * D56.x), 2.0))));
				fprintf(stdout, "\n%e %e %e %e %e %e \n", S1357, S0145, S0246, S2367, S0123, S4567);

				/*
				 * NORMAL VECTOR definition Normal vector of
				 * a plane with known two vectors can be
				 * calculated by finding the cross product of
				 * these two vectors. Here, AC X BD will give
				 * the normal vector to the plane where AC
				 * and BD lie. Also, X[0] = i component X[1]
				 * = j component X[2] = k component
				 * 
				 * These vector are constructed in a way to have
				 * their norms equal to the surface S to
				 * which they are perpendicular e.g.  M = 0.5
				 * ( AC X BD )
				 */

				N1357 =  vector_product(D17 , D35);
				N0145 =  vector_product(D05 , D14);
				N0246 =  vector_product(D24 , D06);
				N2367 =  vector_product(D36 , D27);
				N0123 =  vector_product(D12 , D03);
				N4567 =  vector_product(D47 , D56);			       	       

 #ifdef DEBUG_HARD
	/* Each processor prints its mesh */

	  if(ROOT){ 
	sprintf(fnamein, "normal.%d.out", me);
	fout = fopen(fnamein, "a");
	fprintf(fout, "N1357 %e %e %e\n", N1357.x , N1357.y , N1357.z); fflush(fout);
	fprintf(fout, "N0145 %e %e %e\n", N0145.x , N0145.y , N0145.z); fflush(fout);
	fprintf(fout, "N0246 %e %e %e\n", N0246.x , N0246.y , N0246.z); fflush(fout);
	fprintf(fout, "N2367 %e %e %e\n", N2367.x , N2367.y , N2367.z); fflush(fout);
	fprintf(fout, "N0123 %e %e %e\n", N0123.x , N0123.y , N0123.z); fflush(fout);
	fprintf(fout, "N4567 %e %e %e\n", N4567.x , N4567.y , N4567.z); fflush(fout);
		    fclose(fout);
	  }
 #endif 				

				/*
				 * volume definition Ref:  Jeffrey Grandy,
				 * Efficient Computation of Volume of
				 * Hexahedral Cells, Lawrence Livermore
				 * National Laboratory
				 */
				
				M0167 = vector_sum(D17, D06);				
				M0257 = vector_sum(D05, D27);				
				M0347 = vector_sum(D03, D47);

			
	V1 = fabs(M0167.x * (D27.y * D03.z - D27.z * D03.y) - D27.x * (M0167.y * D03.z - M0167.z * D03.y) + D03.x * (M0167.y * D27.z - M0167.z * D27.y));			
	V2 = fabs(D06.x * (M0257.y * D47.z - M0257.z * D47.y) - M0257.x * (D06.y * D47.z - D06.z * D47.y) + D47.x * (D06.y * M0257.z - D06.z * M0257.y));			
	V3 = fabs(D17.x * (D05.y * M0347.z - D05.z * M0347.y) - D05.x * (D17.y * M0347.z - D17.z * M0347.y) + M0347.x * (D17.y * D05.z - D17.z * D05.y));			
	V = (V1 + V2 + V3)/12.0;

				/* testing volume */
				fprintf(stdout, "%e \n", V);

				center_V[IDX(i, j, k)] = P8;


	for(pp=0;pp<NPOP;pp++){	
	  coeff_xp[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N1357, c[pp])/V;
	  coeff_xm[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N0246, c[pp])/V;
	  coeff_yp[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N2367, c[pp])/V;
          coeff_ym[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N0145, c[pp])/V; 
	  coeff_zp[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N4567, c[pp])/V; 
	  coeff_zm[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N0123, c[pp])/V; 

 #ifdef DEBUG_HARD
	/* Each processor prints its mesh */

	  if(ROOT){ 
	sprintf(fnamein, "coeff.%d.out", me);
	fout = fopen(fnamein, "a");
	fprintf(fout, "%d %d %d %d coeff_xp %e coeff_xm %e\n",i,j,k,pp,coeff_xp[IDX(i, j, k)].p[pp], coeff_xm[IDX(i, j, k)].p[pp] ); fflush(fout);
	fprintf(fout, "%d %d %d %d coeff_yp %e coeff_ym %e\n",i,j,k,pp,coeff_yp[IDX(i, j, k)].p[pp], coeff_ym[IDX(i, j, k)].p[pp] ); fflush(fout);
	fprintf(fout, "%d %d %d %d coeff_zp %e coeff_zm %e\n",i,j,k,pp,coeff_zp[IDX(i, j, k)].p[pp], coeff_zm[IDX(i, j, k)].p[pp] ); fflush(fout);
		    fclose(fout);
	  }
 #endif 
 			}/* for pp */	
	}  /* for i , j , k */

}



