#include "common_object.h"

#ifdef GRID_REFINED
void make_grid_rulers(my_double fac){
  int i;
  my_double dx;

  dx=2./(my_double)NX;
  for(i=0; i< NXG; i++){	
    /* 1) build an array with NGX points, uniformly spaced points in the interval -1,1 */ 
    grid_ruler_x[i] = -1.0+dx*((double)i);
    /* 2) build an array with NGX points, clustered to the walls in the interval -1,1 */ 
    grid_ruler_x[i] = (1./fac)*tanh(grid_ruler_x[i]*atanh(fac));
    /* 3) rescale the -1,1 interval to the 0,NXG interval*/
    grid_ruler_x[i] = ((my_double)NX/2.0)*(grid_ruler_x[i]+1.0);					  
 }
		
  dx=2./(my_double)NY;			 
  for(i=0; i< NYG; i++){					    
    grid_ruler_y[i] = -1.0+dx*((double)i);
    grid_ruler_y[i] = (1./fac)*tanh(grid_ruler_y[i]*atanh(fac));
    grid_ruler_y[i] = ((my_double)NY/2.0)*(grid_ruler_y[i]+1.0);
    // if(ROOT) fprintf(stderr,"grid_ruler_y[%d] = %e\n",i,grid_ruler_y[i]);
 }

  dx=2./(my_double)NZ;
  for(i=0; i< NZG; i++){					    
    grid_ruler_z[i] = -1.0+dx*((double)i);
    grid_ruler_z[i] = (1./fac)*tanh(grid_ruler_z[i]*atanh(fac));
    grid_ruler_z[i] = ((my_double)NZ/2.0)*(grid_ruler_z[i]+1.0);					  
 }

}
#endif

void read_mesh(){

	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;
	int             i, j, k;

#ifdef GRID_REFINED
	my_double stretch=0.9;
#endif
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


#ifdef GRID_REFINED					  
		 make_grid_rulers(stretch);
#endif

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



#ifdef GRID_REFINED					  			       
					  //mesh[IDXG(i, j, k)].x = grid_ruler_x[i+LNXG_START-BRD];					  
					  mesh[IDXG(i, j, k)].y = grid_ruler_y[j+LNYG_START-BRD];
					  //mesh[IDXG(i, j, k)].z = grid_ruler_z[k+LNZG_START-BRD];
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
/**************************************/
/* Useful functions on vectors        */

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

vector vector_norm(vector a){
  my_double norm;
  norm = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  a.x /= norm; 
  a.y /= norm; 
  a.z /= norm; 
  return a;
}


vector average_4vectors(vector a , vector b, vector c, vector d){
  vector e;
  e.x =  0.25*(a.x + b.x + c.x + d.x);
  e.y =  0.25*(a.y + b.y + c.y + d.y);
  e.z =  0.25*(a.z + b.z + c.z + d.z);
  return e;
}

my_double vector_length(vector a){
  my_double norm;
  norm = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  return norm;
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
	  /*
	  if(LNY_START == 0 && j == LNY-TWO_BRD-1) coeff_xp[IDX(i, j, k)].p[pp]=0.0;
	  if(LNY_END == NY  && j == BRD) coeff_xm[IDX(i, j, k)].p[pp]=0.0;
	  */
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


//#ifdef METHOD_CENTERED

#if (defined METHOD_CENTERED || defined METHOD_QUICK)
/****************************************************************************************************/ 
void sendrecv_borders_scalar(my_double *f){
  int i,j,k,brd_size;
  MPI_Status status1;

  /*     BRD|LNX|BRD     */
  /* Copy borders along x */
  brd_size = BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_scalar[IDX_XBRD(i,j,k)] = f[IDX(i+LNX,j,k)];
      }

  MPI_Sendrecv( xp_scalar, brd_size, MPI_my_double_type, me_xp, 10,
                xm_scalar, brd_size, MPI_my_double_type, me_xm, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_scalar[IDX_XBRD(i,j,k)];
        xm_scalar[IDX_XBRD(i,j,k)] = f[IDX(i+BRD,j,k)];
      }
 MPI_Sendrecv( xm_scalar, brd_size, MPI_my_double_type, me_xm, 11,
               xp_scalar, brd_size, MPI_my_double_type, me_xp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k)] = xp_scalar[IDX_XBRD(i,j,k)];
      }



  /* Copy borders along y */
  brd_size = BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_scalar[IDX_YBRD(i,j,k)] = f[IDX(i,j+LNY,k)];
      }

  MPI_Sendrecv( yp_scalar, brd_size, MPI_my_double_type, me_yp, 10,
                ym_scalar, brd_size, MPI_my_double_type, me_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = ym_scalar[IDX_YBRD(i,j,k)];
	ym_scalar[IDX_YBRD(i,j,k)] = f[IDX(i,j+BRD,k)];
      }
  
 MPI_Sendrecv( ym_scalar, brd_size, MPI_my_double_type, me_ym, 13,
               yp_scalar, brd_size, MPI_my_double_type, me_yp, 13, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
       	f[IDX(i,j+LNY+BRD,k)] = yp_scalar[IDX_YBRD(i,j,k)];
      }
  

  /* Copy borders along z */
  brd_size = BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        zp_scalar[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+LNZ)];
      }

  MPI_Sendrecv( zp_scalar, brd_size, MPI_my_double_type, me_zp, 14,
                zm_scalar, brd_size, MPI_my_double_type, me_zm, 14, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = zm_scalar[IDX_ZBRD(i,j,k)];
        zm_scalar[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+BRD)];
      }
 MPI_Sendrecv( zm_scalar, brd_size, MPI_my_double_type, me_zm, 15,
               zp_scalar, brd_size, MPI_my_double_type, me_zp, 15, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j,k+LNZ+BRD)] = zp_scalar[IDX_ZBRD(i,j,k)];
      }

}/* end send rcv */

#endif


/**************************************************************************************************/
/* Here begin the function to compute the interpolation coefficients */

void compute_interpolation_coefficients(){

  int i,j,k;
  vector xp,xm,yp,ym,zp,zm;
  vector P0, P1, P2, P3, P4, P5, P6, P7, P8;

  /* for quick */
  vector XPF1, XPF2, XPB1, XPB2, XMF1, XMF2, XMB1, XMB2;
  vector YPF1, YPF2, YPB1, YPB2, YMF1, YMF2, YMB1, YMB2;   
  vector ZPF1, ZPF2, ZPB1, ZPB2, ZMF1, ZMF2, ZMB1, ZMB2;


#ifdef METHOD_CENTERED
	for (i = BRD; i < LNXG + BRD - 1 ; i++) 
		for (j = BRD; j < LNYG + BRD - 1; j++) 
			for (k = BRD; k < LNZG + BRD - 1; k++) {
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

/* we compute here the position of the centers of the 6 surfaces */
xp = average_4vectors(P1, P3, P5, P7);
xm = average_4vectors(P0, P2, P4, P6); 
yp = average_4vectors(P2, P3, P6, P7); 
ym = average_4vectors(P0, P1, P4, P5);
zp = average_4vectors(P4, P5, P6, P7);
zm = average_4vectors(P0, P1, P2, P3); 

/* we find the distence of the centers from such 6 surfece centers */

xp = vector_difference( center_V[IDX(i, j, k)], xp);
xm = vector_difference( center_V[IDX(i, j, k)], xm);
yp = vector_difference( center_V[IDX(i, j, k)], yp);
ym = vector_difference( center_V[IDX(i, j, k)], ym);
zp = vector_difference( center_V[IDX(i, j, k)], zp);
zm = vector_difference( center_V[IDX(i, j, k)], zm);

/* we compute the vector length (norm) */
interp_xp[IDX(i,j,k)] = vector_length(xp);
interp_xm[IDX(i,j,k)] = vector_length(xm);
interp_yp[IDX(i,j,k)] = vector_length(yp);
interp_ym[IDX(i,j,k)] = vector_length(ym);
interp_zp[IDX(i,j,k)] = vector_length(zp);
interp_zm[IDX(i,j,k)] = vector_length(zm);

			}/*i,j,k*/

	/*  send receive the borders */
	sendrecv_borders_scalar(interp_xp);
	sendrecv_borders_scalar(interp_xm);
	sendrecv_borders_scalar(interp_yp);
	sendrecv_borders_scalar(interp_ym);
	sendrecv_borders_scalar(interp_zp);
	sendrecv_borders_scalar(interp_zm);

	/* now compute the coefficients */
	for (i = BRD; i < LNX + BRD; i++) 
		for (j = BRD; j < LNY + BRD; j++) 
			for (k = BRD; k < LNZ + BRD; k++) {

interp_xp[IDX(i,j,k)] =  interp_xp[IDX(i,j,k)] / (interp_xm[IDX(i+1,j,k)]+interp_xp[IDX(i,j,k)]);
interp_xm[IDX(i,j,k)] =  interp_xm[IDX(i,j,k)] / (interp_xp[IDX(i-1,j,k)]+interp_xm[IDX(i,j,k)]);
interp_yp[IDX(i,j,k)] =  interp_yp[IDX(i,j,k)] / (interp_ym[IDX(i,j+1,k)]+interp_yp[IDX(i,j,k)]);
interp_ym[IDX(i,j,k)] =  interp_ym[IDX(i,j,k)] / (interp_yp[IDX(i,j-1,k)]+interp_ym[IDX(i,j,k)]);
interp_zp[IDX(i,j,k)] =  interp_zp[IDX(i,j,k)] / (interp_zm[IDX(i,j,k+1)]+interp_zp[IDX(i,j,k)]);
interp_zm[IDX(i,j,k)] =  interp_zm[IDX(i,j,k)] / (interp_zp[IDX(i,j,k-1)]+interp_zm[IDX(i,j,k)]);

			}/* i,j,k*/

#endif

#ifdef METHOD_QUICK
	for (i = BRD; i < LNXG + BRD - 1 ; i++) 
		for (j = BRD; j < LNYG + BRD - 1; j++) 
			for (k = BRD; k < LNZG + BRD - 1; k++) {
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

/* we compute here the position of the centers of the 6 surfaces */
xp = average_4vectors(P1, P3, P5, P7);
xm = average_4vectors(P0, P2, P4, P6); 
yp = average_4vectors(P2, P3, P6, P7); 
ym = average_4vectors(P0, P1, P4, P5);
zp = average_4vectors(P4, P5, P6, P7);
zm = average_4vectors(P0, P1, P2, P3); 


/* we find the coefficients for XP  */

 XPF1 =  ((xp -  center_V[IDX(i, j, k)])*(xp -  center_V[IDX(i-1, j, k)]))/(( center_V[IDX(i+1, j, k)] -  center_V[IDX(i, j, k)]) * ( center_V[IDX(i+1, j, k)] -  center_V[IDX(i-1, j, k)]));
 XPF2 =  ((xp -  center_V[IDX(i, j, k)])*( center_V[IDX(i+1, j, k)] - xp ))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i-1, j, k)]) * ( center_V[IDX(i+1, j, k)] -  center_V[IDX(i-1, j, k)]));

XPB1 = ((xp -  center_V[IDX(i+1, j, k)])*(xp -  center_V[IDX(i+2, j, k)]))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i+1, j, k)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i+2, j, k)]));
XPB2 = ((xp -  center_V[IDX(i+1, j, k)])*(center_V[IDX(i, j, k)] - xp))/(( center_V[IDX(i+1, j, k)] -  center_V[IDX(i+2, j, k)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i+2, j, k)]));

/* we find the coefficients for XM  */

XMF1 =  ((xm -  center_V[IDX(i, j, k)])*(xm -  center_V[IDX(i+1, j, k)]))/(( center_V[IDX(i-1, j, k)] -  center_V[IDX(i, j, k)]) * ( center_V[IDX(i-1, j, k)] -  center_V[IDX(i+1, j, k)]));
XMF2 =  ((xm -  center_V[IDX(i, j, k)])*(center_V[IDX(i-1, j, k)] - xm))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i+1, j, k)]) * ( center_V[IDX(i-1, j, k)] -  center_V[IDX(i+1, j, k)]));

XMB1 =  ((xm -  center_V[IDX(i-1, j, k)])*(xm -  center_V[IDX(i-2, j, k)]))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i-1, j, k)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i-2, j, k)]));
XMB2 =  ((xm -  center_V[IDX(i-1, j, k)])*(center_V[IDX(i, j, k)] - xm))/(( center_V[IDX(i-1, j, k)] -  center_V[IDX(i-2, j, k)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i-2, j, k)]));

/* we find the coefficients for YP  */

 YPF1 =  ((yp -  center_V[IDX(i, j, k)])*(yp -  center_V[IDX(i, j-1, k)]))/(( center_V[IDX(i, j+1, k)] -  center_V[IDX(i, j, k)]) * ( center_V[IDX(i, j+1, k)] -  center_V[IDX(i, j-1, k)]));
 YPF2 =  ((yp -  center_V[IDX(i, j, k)])*( center_V[IDX(i, j+1, k)] - yp ))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i, j-1, k)]) * ( center_V[IDX(i, j+1, k)] -  center_V[IDX(i, j-1, k)]));

YPB1 = ((yp -  center_V[IDX(i, j+1, k)])*(yp -  center_V[IDX(i, j+2, k)]))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i, j+1, k)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i, j+2, k)]));
YPB2 = ((yp -  center_V[IDX(i, j+1, k)])*(center_V[IDX(i, j, k)] - yp))/(( center_V[IDX(i, j+1, k)] -  center_V[IDX(i, j+2, k)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i, j+2, k)]));

/* we find the coefficients for YM  */
YMF1 =  ((ym -  center_V[IDX(i, j, k)])*(ym -  center_V[IDX(i, j+1, k)]))/(( center_V[IDX(i, j-1, k)] -  center_V[IDX(i, j, k)]) * ( center_V[IDX(i, j-1, k)] -  center_V[IDX(i, j+1, k)]));
YMF2 =  ((ym -  center_V[IDX(i, j, k)])*(center_V[IDX(i, j-1, k)] - ym))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i, j+1, k)]) * ( center_V[IDX(i, j-1, k)] -  center_V[IDX(i, j+1, k)]));

YMB1 =  ((ym -  center_V[IDX(i, j-1, k)])*(ym -  center_V[IDX(i, j-2, k)]))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i, j-1, k)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i, j-2, k)]));
YMB2 =  ((ym -  center_V[IDX(i, j-1, k)])*(center_V[IDX(i, j, k)] - ym))/(( center_V[IDX(i, j-1, k)] -  center_V[IDX(i, j-2, k)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i, j-2, k)]));


/* we find the coefficients for ZP  */

 ZPF1 =  ((zp -  center_V[IDX(i, j, k)])*(zp -  center_V[IDX(i, j, k-1)]))/(( center_V[IDX(i, j, k+1)] -  center_V[IDX(i, j, k)]) * ( center_V[IDX(i, j, k+1)] -  center_V[IDX(i, j, k-1)]));
 ZPF2 =  ((zp -  center_V[IDX(i, j, k)])*( center_V[IDX(i, j, k+1)] - zp ))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i, j, k-1)]) * ( center_V[IDX(i, j, k+1)] -  center_V[IDX(i, j, k-1)]));

ZPB1 = ((zp -  center_V[IDX(i, j, k+1)])*(zp -  center_V[IDX(i, j, k+2)]))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i, j, k+1)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i, j, k+2)]));
ZPB2 = ((zp -  center_V[IDX(i, j, k+1)])*(center_V[IDX(i, j, k)] - zp))/(( center_V[IDX(i, j, k+1)] -  center_V[IDX(i, j, k+2)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i, j, k+2)]));

/* we find the coefficients for ZM  */

ZMF1 =  ((zm -  center_V[IDX(i, j, k)])*(zm -  center_V[IDX(i, j, k+1)]))/(( center_V[IDX(i, j, k-1)] -  center_V[IDX(i, j, k)]) * ( center_V[IDX(i, j, k-1)] -  center_V[IDX(i, j, k+1)]));
ZMF2 =  ((zm -  center_V[IDX(i, j, k)])*(center_V[IDX(i, j, k-1)] - zm))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i, j, k+1)]) * ( center_V[IDX(i, j, k-1)] -  center_V[IDX(i, j, k+1)]));

ZMB1 =  ((zm -  center_V[IDX(i, j, k-1)])*(zm -  center_V[IDX(i, j, k-2)]))/(( center_V[IDX(i, j, k)] -  center_V[IDX(i, j, k-1)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i, j, k-2)]));
ZMB2 =  ((zm -  center_V[IDX(i, j, k-1)])*(center_V[IDX(i, j, k)] - zm))/(( center_V[IDX(i, j, k-1)] -  center_V[IDX(i, j, k-2)]) * ( center_V[IDX(i, j, k)] -  center_V[IDX(i, j, k-2)]));


/* we compute the vector length (norm) of the coefficients G1 and G2 for QUICK */
quick_xpf1[IDX(i,j,k)] = vector_length(XPF1);
quick_xpf2[IDX(i,j,k)] = vector_length(XPF2);
quick_xpb1[IDX(i,j,k)] = vector_length(XPB1);
quick_xpb2[IDX(i,j,k)] = vector_length(XPB2);
quick_xmf1[IDX(i,j,k)] = vector_length(XMF1);
quick_xmf2[IDX(i,j,k)] = vector_length(XMF2);
quick_xmb1[IDX(i,j,k)] = vector_length(XMB1);
quick_xmb2[IDX(i,j,k)] = vector_length(XMB2);
quick_ypf1[IDX(i,j,k)] = vector_length(YPF1);
quick_ypf2[IDX(i,j,k)] = vector_length(YPF2);
quick_ypb1[IDX(i,j,k)] = vector_length(YPB1);
quick_ypb2[IDX(i,j,k)] = vector_length(YPB2);
quick_ymf1[IDX(i,j,k)] = vector_length(YMF1);
quick_ymf2[IDX(i,j,k)] = vector_length(YMF2);
quick_ymb1[IDX(i,j,k)] = vector_length(YMB1);
quick_ymb2[IDX(i,j,k)] = vector_length(YMB2);
quick_zpf1[IDX(i,j,k)] = vector_length(ZPF1);
quick_zpf2[IDX(i,j,k)] = vector_length(ZPF2);
quick_zpb1[IDX(i,j,k)] = vector_length(ZPB1);
quick_zpb2[IDX(i,j,k)] = vector_length(ZPB2);
quick_zmf1[IDX(i,j,k)] = vector_length(ZMF1);
quick_zmf2[IDX(i,j,k)] = vector_length(ZMF2);
quick_zmb1[IDX(i,j,k)] = vector_length(ZMB1);
quick_zmb2[IDX(i,j,k)] = vector_length(ZMB2);


			}/*i,j,k*/

	/*  send receive the borders */
	sendrecv_borders_scalar(quick_xpf1);
	sendrecv_borders_scalar(quick_xpf2);
	sendrecv_borders_scalar(quick_xpb1);
	sendrecv_borders_scalar(quick_xpb2);
	sendrecv_borders_scalar(quick_xmf1);
	sendrecv_borders_scalar(quick_xmf2);
	sendrecv_borders_scalar(quick_xmb1);
	sendrecv_borders_scalar(quick_xmb2);
	sendrecv_borders_scalar(quick_ypf1);
	sendrecv_borders_scalar(quick_ypf2);
	sendrecv_borders_scalar(quick_ypb1);
	sendrecv_borders_scalar(quick_ypb2);
	sendrecv_borders_scalar(quick_ymf1);
	sendrecv_borders_scalar(quick_ymf2);
	sendrecv_borders_scalar(quick_ymb1);
	sendrecv_borders_scalar(quick_ymb2);
	sendrecv_borders_scalar(quick_zpf1);
	sendrecv_borders_scalar(quick_zpf2);
	sendrecv_borders_scalar(quick_zpb1);
	sendrecv_borders_scalar(quick_zpb2);
	sendrecv_borders_scalar(quick_zmf1);
	sendrecv_borders_scalar(quick_zmf2);
	sendrecv_borders_scalar(quick_zmb1);
	sendrecv_borders_scalar(quick_zmb2);

#endif


}/* end func interpolation coefficients */




/**************************************************************************************************/
#ifdef LB_BC
void prepare_boundary_conditions(){
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
        vector cn[NPOP];


	if(LNX_END == NX){
	  norm_xp_pop=(pop*) malloc(sizeof(pop)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
          if(norm_xp_pop == NULL){fprintf(stderr,"Not enough memory to allocate norm_xp_pop\n"); exit(-1);}
	}
	if(LNX_START == 0){
	  norm_xm_pop=(pop*) malloc(sizeof(pop)*(LNY+TWO_BRD)*(LNZ+TWO_BRD));
         if(norm_xm_pop == NULL){fprintf(stderr,"Not enough memory to allocate norm_xm_pop\n"); exit(-1);}
	}

	if(LNY_END == NY){
	  norm_yp_pop  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
         if(norm_yp_pop == NULL){fprintf(stderr,"Not enough memory to allocate norm_yp_pop\n"); exit(-1);}
	}
	if(LNY_START == 0){
	  norm_ym_pop  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNZ+TWO_BRD)); 
         if(norm_ym_pop == NULL){fprintf(stderr,"Not enough memory to allocate norm_ym_pop\n"); exit(-1);}
	}

	if(LNZ_END == NZ){
	  norm_zp_pop  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD)); 
         if(norm_zp_pop == NULL){fprintf(stderr,"Not enough memory to allocate norm_zp_pop\n"); exit(-1);}
	}

	if(LNZ_START == 0){
	  norm_zm_pop  = (pop*) malloc(sizeof(pop)*(LNX+TWO_BRD)*(LNY+TWO_BRD));
         if(norm_zm_pop == NULL){fprintf(stderr,"Not enough memory to allocate norm_zm_pop\n"); exit(-1);}
	}
 

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

				/* diagonals */
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
				
	/* get the outgoing normal vectors from the vector product of the diagonals */
				N1357 =  vector_product(D17 , D35);
				N0145 =  vector_product(D05 , D14);
				N0246 =  vector_product(D24 , D06);
				N2367 =  vector_product(D36 , D27);
				N0123 =  vector_product(D12 , D03);
				N4567 =  vector_product(D47 , D56);

				/* Normalization */
			       	N1357 =  vector_norm(N1357);
				N0145 =  vector_norm(N0145);
				N0246 =  vector_norm(N0246);
				N2367 =  vector_norm(N2367);
				N0123 =  vector_norm(N0123);
				N4567 =  vector_norm(N4567);				

				for(pp=0;pp<NPOP;pp++) cn[pp] =  vector_norm(c[pp]);
				cn[0].x=cn[0].y=cn[0].z=0.0;
		
	for(pp=0;pp<NPOP;pp++){	
	if(LNX_END == NX)   norm_xp_pop[IDX_X(j, k)].p[pp] = scalar_product(N1357, cn[pp]);

	if(LNX_START == 0)  norm_xm_pop[IDX_X(j, k)].p[pp] = scalar_product(N0246, cn[pp]);

	if(LNY_END == NY)   norm_yp_pop[IDX_Y(i, k)].p[pp] = scalar_product(N2367, cn[pp]);

	if(LNY_START == 0)  norm_ym_pop[IDX_Y(i, k)].p[pp] = scalar_product(N0145, cn[pp]); 

	if(LNZ_END == NZ)   norm_zp_pop[IDX_Z(i, j)].p[pp] = scalar_product(N4567, cn[pp]); 

	if(LNZ_START == 0)  norm_zm_pop[IDX_Z(i, j)].p[pp] = scalar_product(N0123, cn[pp]); 
 			}/* for pp */

	
	}  /* for i , j , k */

}
#endif
