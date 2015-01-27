#include "common_object.h"



#ifdef GRID_REFINED
void make_grid_rulers(my_double fac){
  int i;
  my_double dx,dy,dz;
  my_double fac2 = 2.0;
#ifdef GRID_REFINED_ONESIDED
  fac2 = 1.0;
#endif

  //dx=2./(my_double)NX;
  dx=fac2/(my_double)NX;
  for(i=0; i< NXG; i++){	
    /* 1) build an array with NGX points, uniformly spaced points in the interval -1,1 */ 
    grid_ruler_x[i] = -1.0+dx*((double)i);
    /* 2) build an array with NGX points, clustered to the walls in the interval -1,1 */ 
    grid_ruler_x[i] = (1./fac)*tanh(grid_ruler_x[i]*atanh(fac));
    /* 3) rescale the -1,1 interval to the 0,SX interval*/
    grid_ruler_x[i] = (property.SX/fac2)*(grid_ruler_x[i]+1.0);
 }
		
  //dy=2./(my_double)NY;			 
  dy=fac2/(my_double)NY;
  for(i=0; i< NYG; i++){					    
    grid_ruler_y[i] = -1.0+dy*((double)i);
    grid_ruler_y[i] = (1./fac)*tanh(grid_ruler_y[i]*atanh(fac));
    grid_ruler_y[i] = (property.SY/fac2)*(grid_ruler_y[i]+1.0);
    // if(ROOT) fprintf(stderr,"grid_ruler_y[%d] = %e\n",i,grid_ruler_y[i]);
 }

  //dz=2./(my_double)NZ;
  dz=fac2/(my_double)NZ;
  for(i=0; i< NZG; i++){					    
    grid_ruler_z[i] = -1.0+dz*((double)i);
    grid_ruler_z[i] = (1./fac)*tanh(grid_ruler_z[i]*atanh(fac));
    grid_ruler_z[i] = (property.SZ/fac2)*(grid_ruler_z[i]+1.0);					  
 }

#ifdef GRID_REFINED_BULK
  dx=fac2/(my_double)NX;
  for(i=0; i< NXG; i++){	
    /* 1) build an array with NGX points, uniformly spaced points in the interval -1,1 */ 
    grid_ruler_x[i] = -1.0+dx*((double)i);
    /* 2) build an array with NGX points, clustered to the walls in the interval -1,1 */ 
    grid_ruler_x[i] = (1.-fac)*sinh(grid_ruler_x[i]*asinh(1./(1.-fac)));
    /* 3) rescale the -1,1 interval to the 0,SX interval*/
    grid_ruler_x[i] = (property.SX/fac2)*(grid_ruler_x[i]+1.0);
 }

 
  dy=fac2/(my_double)NY;
  for(i=0; i< NYG; i++){					    
    grid_ruler_y[i] = -1.0+dy*((double)i);
    grid_ruler_y[i] = (1.-fac)*sinh(grid_ruler_y[i]*asinh(1./(1.-fac)));
    grid_ruler_y[i] = (property.SY/fac2)*(grid_ruler_y[i]+1.0);
    // if(ROOT) fprintf(stderr,"grid_ruler_y[%d] = %e\n",i,grid_ruler_y[i]);
 }

  dz=fac2/(my_double)NZ;
  for(i=0; i< NZG; i++){					    
    grid_ruler_z[i] = -1.0+dz*((double)i);
    grid_ruler_z[i] = (1.-fac)*tanh(grid_ruler_z[i]*atanh(1./(1.-fac)));
    grid_ruler_z[i] = (property.SZ/fac2)*(grid_ruler_z[i]+1.0);					  
 }
#endif

#ifdef GRID_REFINED_CHEBYSHEV /* Chebyshev polynomial nodes , from http://en.wikipedia.org/wiki/Chebyshev_nodes */
  /* x */
  for(i=0; i< NXG; i++){	
    grid_ruler_x[NX-i] = 0.5*property.SX*(1.0 + cos((i+0.5)*(one_pi/NXG)) );
 }
	
  /* y */	
  for(i=0; i< NYG; i++){					    
     grid_ruler_y[NY-i] = 0.5*property.SY*(1.0 + cos((i+0.5)*(one_pi/NYG)) );
 }

  /* z */	
  for(i=0; i< NZG; i++){					    
     grid_ruler_z[NZ-i] = 0.5*property.SZ*(1.0 + cos((i+0.5)*(one_pi/NZG)) );				  
 }
#endif

#ifdef GRID_REFINED_SINH
/*  Rudy Kunnen spacing 
     ly = 1.0d0
     ny = 64

     ny2 = ny/2
     by=6.5d0
     ay=0.5d0*ly/SINH(0.5d0*by)

     do 120 j=0,ny2 
     y(j) = ay*SINH(by*j/(1.0d0*ny))
120  continue

     do 130 j=ny2+1,ny     
     y(j) = ly - y(ny-j) 
130  continue      
*/
  fac2=6.5; /* stretching factor */

  /* x */
  dx=0.5*property.SX/sinh(0.5*fac2);
  for(i=0; i< ceil(NXG/2)+1; i++) grid_ruler_x[i+1] = dx*sinh(fac2*(i+0.5)/NX);
  grid_ruler_x[0]=0.0;
  for(i=1; i< ceil(NXG/2)+1; i++) grid_ruler_x[i] = 0.5*(grid_ruler_x[i+1] + grid_ruler_x[i]);
  for(i=1; i< ceil(NXG/2)+1; i++) grid_ruler_x[NX-i] = property.SX -  grid_ruler_x[i];

  /* y */		
  dy=0.5*property.SY/sinh(0.5*fac2);
  for(i=0; i< ceil(NYG/2)+1; i++) grid_ruler_y[i+1] = dy*sinh(fac2*(i+0.5)/NY);
  grid_ruler_y[0]=0.0;
  for(i=1; i< ceil(NYG/2)+1; i++) grid_ruler_y[i] = 0.5*(grid_ruler_y[i+1] + grid_ruler_y[i]);
  for(i=1; i< ceil(NYG/2)+1; i++) grid_ruler_y[NY-i] = property.SY -  grid_ruler_y[i];

  /* z */		
  dz=0.5*property.SZ/sinh(0.5*fac2);
  for(i=0; i< ceil(NZG/2)+1; i++) grid_ruler_z[i+1] = dz*sinh(fac2*(i+0.5)/NZ);
  grid_ruler_z[0]=0.0;
  for(i=1; i< ceil(NZG/2)+1; i++) grid_ruler_z[i] = 0.5*(grid_ruler_z[i+1] + grid_ruler_z[i]);
  for(i=1; i< ceil(NZG/2)+1; i++) grid_ruler_z[NZ-i] = property.SZ -  grid_ruler_z[i];

#endif


}
#endif




void read_mesh(){

	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;
	int             i, j, k;

#ifdef GRID_REFINED
	my_double stretch=0.98;
#endif
	sprintf(fnamein, "mesh.in");
	fin = fopen(fnamein, "r");
	if (fin != NULL) {
	  if(ROOT) fprintf(stderr, "Mesh file %s has been found!\n", fnamein);
	} else {

	  if(ROOT) fprintf(stderr, "Warning message -> %s file is missing!\n Starting from grid generated on the fly\n ", fnamein);


		/* set field to zero */
		for (k =0; k < LNZG+TWO_BRD; k++)
			for (j =0; j < LNYG+TWO_BRD; j++)
				for (i = 0; i < LNXG+TWO_BRD; i++) {
				  mesh[IDXG(i, j, k)].x = 0.0;
				  mesh[IDXG(i, j, k)].y = 0.0;
				  mesh[IDXG(i, j, k)].z = 0.0;			 
				  //mesh_flag[IDXG(i, j, k)] = 0;
				}


#ifdef GRID_REFINED					  
		 make_grid_rulers(stretch);
#endif

		/* moving on the bulk only */
		for (k = BRD; k < LNZG+BRD; k++)
			for (j = BRD; j < LNYG+BRD; j++)
				for (i = BRD; i < LNXG+BRD; i++) {
				        mesh[IDXG(i, j, k)].x = (my_double) (i + LNXG_START-BRD)*property.SX/property.NX;
					mesh[IDXG(i, j, k)].y = (my_double) (j + LNYG_START-BRD)*property.SY/property.NY;
					mesh[IDXG(i, j, k)].z = (my_double) (k + LNZG_START-BRD)*property.SZ/property.NZ;
					/*
					 * flag: 1 is bulk , 0 is wall , -1
					 * is dormient
					 */
					//mesh_flag[IDXG(i, j, k)] = 1;
								
#ifdef GRID_RANDOM
					  if(i<LNXG+BRD-1) mesh[IDXG(i, j, k)].x += 0.25*(my_double)(2.0*myrand()-1.0);
					  if(j<LNYG+BRD-1) mesh[IDXG(i, j, k)].y += 0.25*(my_double)(2.0*myrand()-1.0);
					  if(j<LNZG+BRD-1) mesh[IDXG(i, j, k)].z += 0.25*(my_double)(2.0*myrand()-1.0);
#endif



#ifdef GRID_REFINED
#ifdef GRID_REFINED_X					  			       
					  mesh[IDXG(i, j, k)].x = grid_ruler_x[i+LNXG_START-BRD];
#endif
#ifdef GRID_REFINED_Y					  
					  mesh[IDXG(i, j, k)].y = grid_ruler_y[j+LNYG_START-BRD];
#endif
#ifdef GRID_REFINED_Z
					  mesh[IDXG(i, j, k)].z = grid_ruler_z[k+LNZG_START-BRD];
#endif
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
     //fprintf(fout, "%d %d %d %e %e %e %d\n", i, j, k, mesh[IDXG(i, j, k)].x, mesh[IDXG(i, j, k)].y, mesh[IDXG(i, j, k)].z , mesh_flag[IDXG(i, j, k)]);
       fprintf(fout, "%d %d %d %e %e %e\n", i, j, k, mesh[IDXG(i, j, k)].x, mesh[IDXG(i, j, k)].y, mesh[IDXG(i, j, k)].z );
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
      //fprintf(fout, "%d %d %d %e %e %e %d\n", i, j, k, mesh[IDXG(i, j, k)].x, mesh[IDXG(i, j, k)].y, mesh[IDXG(i, j, k)].z , mesh_flag[IDXG(i, j, k)]);
	fprintf(fout, "%d %d %d %e %e %e\n", i, j, k, mesh[IDXG(i, j, k)].x, mesh[IDXG(i, j, k)].y, mesh[IDXG(i, j, k)].z);
	    }
	    fprintf(fout,"\n");
	  }
	  fprintf(fout,"\n");
	} /* ijk */
	fclose(fout);

#endif
}/* end of read mesh */



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

vector vector_scale(my_double a , vector b){
  vector c;
  c.x =  a*b.x;
  c.y =  a*b.y;
  c.z =  a*b.z;
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

#ifdef DEBUG			
	fprintf(stdout, "\n%e %e %e %e %e %e \n", S1357, S0145, S0246, S2367, S0123, S4567);
#endif
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

#ifdef DEBUG
				/* testing volume */
				fprintf(stdout, "Volume : %e \n", V);
#endif

				center_V[IDX(i, j, k)] = P8;
			    //if(P8.x==0.5)fprintf(stdout, "WWW me %d , %d %d %d , %g %g %g\n",me, i,j,k,center_V[IDX(i, j, k)].x,center_V[IDX(i, j, k)].y,center_V[IDX(i, j, k)].z);


	for(pp=0;pp<NPOP;pp++){	
	  coeff_xp[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N1357, c[pp])/V;
	  coeff_xm[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N0246, c[pp])/V;
	  coeff_yp[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N2367, c[pp])/V;
          coeff_ym[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N0145, c[pp])/V; 
	  coeff_zp[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N4567, c[pp])/V; 
	  coeff_zm[IDX(i, j, k)].p[pp] = 0.5*scalar_product(N0123, c[pp])/V; 
	  /*
	  if(LNY_START == 0 && j == BRD ) coeff_ym[IDX(i, j, k)].p[pp]=0.0;
	  if(LNY_END == NY  && j == LNY-TWO_BRD-1 ) coeff_yp[IDX(i, j, k)].p[pp]=0.0;
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

#ifdef DEBUG_HARD
	  /* check that the sum of the six coefficients sum up to zero :  \Sum_j=1,6 c_i*n_j*S_j/V = 0  */
	  my_double  value;
	  value = coeff_xp[IDX(i, j, k)].p[pp] + coeff_xm[IDX(i, j, k)].p[pp] + coeff_yp[IDX(i, j, k)].p[pp] 
	        + coeff_ym[IDX(i, j, k)].p[pp] + coeff_zp[IDX(i, j, k)].p[pp] + coeff_zm[IDX(i, j, k)].p[pp];
	  fprintf(stdout, "p %d value %e\n",pp,value);
#endif


	  /* If the simulation is 2D we do not want any flux in those directions, therefore we set them to zero */
	  /* Enrico: A test showed that it has no importance */
	  /*
	  if(NX==1){ coeff_xp[IDX(i, j, k)].p[pp] = coeff_xm[IDX(i, j, k)].p[pp] = 0.0;}
	  if(NY==1){ coeff_yp[IDX(i, j, k)].p[pp] = coeff_ym[IDX(i, j, k)].p[pp] = 0.0;}
	  if(NZ==1){ coeff_zp[IDX(i, j, k)].p[pp] = coeff_zm[IDX(i, j, k)].p[pp] = 0.0;}
	  */
 			}/* for pp */	
	}  /* for i , j , k */

}




#if (defined METHOD_CENTERED || defined METHOD_MYQUICK || defined METHOD_STREAMING || defined METHOD_UPWIND)
/****************************************************************************************************/ 
void sendrecv_borders_scalar(my_double *f){
  int i,j,k,brd_size;
  MPI_Status status1;

#ifdef NEW_SENDRECV
  //IDX(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )
     
  MPI_Sendrecv( f + IDX(LNX,0,0)   , 1, MPI_my_double_plane_x, me_xp, 14,
                f + IDX(0,0,0)     , 1, MPI_my_double_plane_x, me_xm, 14, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(BRD,0,0)    , 1, MPI_my_double_plane_x, me_xm, 15,
                f + IDX(LNX+BRD,0,0), 1, MPI_my_double_plane_x, me_xp, 15, MPI_COMM_WORLD, &status1);
  
  MPI_Sendrecv( f + IDX(0,LNY,0)  , 1, MPI_my_double_plane_y, me_yp, 12,
                f + IDX(0,0,0)  , 1, MPI_my_double_plane_y, me_ym, 12, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(0,BRD,0)  , 1, MPI_my_double_plane_y, me_ym, 13,
                f + IDX(0,LNY+BRD,0)  , 1, MPI_my_double_plane_y, me_yp, 13, MPI_COMM_WORLD, &status1);
  
  MPI_Sendrecv( f + IDX(0,0,LNZ)  , 1, MPI_my_double_plane_z, me_zp, 10,
                f + IDX(0,0,0)  , 1, MPI_my_double_plane_z, me_zm, 10, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(0,0,BRD)  , 1, MPI_my_double_plane_z, me_zm, 11,
                f + IDX(0,0,LNZ+BRD)  , 1, MPI_my_double_plane_z, me_zp, 11, MPI_COMM_WORLD, &status1);
  
#else /* = NEW_SENDRECV is not defined */


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

#endif /* NEW_SENDRECV */

#ifdef METHOD_EDGES_AND_CORNERS

  /* First we communicate the 8 corner cubes (they are either 1x1x1 or 2x2x2 depending on BRD) */ 
  
  brd_size = BRD*BRD*BRD;

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_zp_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+LNZ)];
        xp_yp_zm_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+BRD)];
        xp_ym_zp_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+LNZ)];
        xm_yp_zp_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+LNZ)];
      }

  MPI_Sendrecv( xp_yp_zp_corner_scalar, brd_size, MPI_my_double_type, me_xp_yp_zp, 10,
                xm_ym_zm_corner_scalar, brd_size, MPI_my_double_type, me_xm_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_yp_zm_corner_scalar, brd_size, MPI_my_double_type, me_xp_yp_zm, 10,
                xm_ym_zp_corner_scalar, brd_size, MPI_my_double_type, me_xm_ym_zp, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_ym_zp_corner_scalar, brd_size, MPI_my_double_type, me_xp_ym_zp, 10,
                xm_yp_zm_corner_scalar, brd_size, MPI_my_double_type, me_xm_yp_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_zp_corner_scalar, brd_size, MPI_my_double_type, me_xm_yp_zp, 10,
                xp_ym_zm_corner_scalar, brd_size, MPI_my_double_type, me_xp_ym_zm, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_ym_zm_corner_scalar[IDX_CORNER(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = xm_ym_zp_corner_scalar[IDX_CORNER(i,j,k)];
        f[IDX(i,j+LNY+BRD,k)] = xm_yp_zm_corner_scalar[IDX_CORNER(i,j,k)];
        f[IDX(i+LNX+BRD,j,k)] = xp_ym_zm_corner_scalar[IDX_CORNER(i,j,k)];
        xm_ym_zm_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+BRD)];
        xm_ym_zp_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+LNZ)];
        xm_yp_zm_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+BRD)];
        xp_ym_zm_corner_scalar[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+BRD)];
      }
 MPI_Sendrecv( xm_ym_zm_corner_scalar, brd_size, MPI_my_double_type, me_xm_ym_zm, 11,
               xp_yp_zp_corner_scalar, brd_size, MPI_my_double_type, me_xp_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_ym_zp_corner_scalar, brd_size, MPI_my_double_type, me_xm_ym_zp, 11,
               xp_yp_zm_corner_scalar, brd_size, MPI_my_double_type, me_xp_yp_zm, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_yp_zm_corner_scalar, brd_size, MPI_my_double_type, me_xm_yp_zm, 11,
               xp_ym_zp_corner_scalar, brd_size, MPI_my_double_type, me_xp_ym_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_zm_corner_scalar, brd_size, MPI_my_double_type, me_xp_ym_zm, 11,
               xm_yp_zp_corner_scalar, brd_size, MPI_my_double_type, me_xm_yp_zp, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k+LNZ+BRD)] = xp_yp_zp_corner_scalar[IDX_CORNER(i,j,k)];
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)] = xp_yp_zm_corner_scalar[IDX_CORNER(i,j,k)];
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)] = xp_ym_zp_corner_scalar[IDX_CORNER(i,j,k)];
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)] = xm_yp_zp_corner_scalar[IDX_CORNER(i,j,k)];
      }
 

 /* Then we communicate the 12 edges  */
 
 /* along x */
 
 brd_size = BRD*BRD*(LNX+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_zp_edge_scalar[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+LNZ)];
        yp_zm_edge_scalar[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+BRD)];
      }

  MPI_Sendrecv( yp_zp_edge_scalar, brd_size, MPI_my_double_type, me_yp_zp, 10,
                ym_zm_edge_scalar, brd_size, MPI_my_double_type, me_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( yp_zm_edge_scalar, brd_size, MPI_my_double_type, me_yp_zm, 10,
                ym_zp_edge_scalar, brd_size, MPI_my_double_type, me_ym_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = ym_zm_edge_scalar[IDX_EDGE_X(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = ym_zp_edge_scalar[IDX_EDGE_X(i,j,k)];
        ym_zm_edge_scalar[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+BRD)];
        ym_zp_edge_scalar[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+LNZ)];
      }
 MPI_Sendrecv( ym_zm_edge_scalar, brd_size, MPI_my_double_type, me_ym_zm, 11,
               yp_zp_edge_scalar, brd_size, MPI_my_double_type, me_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( ym_zp_edge_scalar, brd_size, MPI_my_double_type, me_ym_zp, 11,
               yp_zm_edge_scalar, brd_size, MPI_my_double_type, me_yp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)] = yp_zp_edge_scalar[IDX_EDGE_X(i,j,k)];
	f[IDX(i,j+LNY+BRD,k)] = yp_zm_edge_scalar[IDX_EDGE_X(i,j,k)];
      }
 
 /* along y */
 
 brd_size = BRD*BRD*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_zp_edge_scalar[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+LNZ)];
        xp_zm_edge_scalar[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+BRD)];
      }

  MPI_Sendrecv( xp_zp_edge_scalar, brd_size, MPI_my_double_type, me_xp_zp, 10,
                xm_zm_edge_scalar, brd_size, MPI_my_double_type, me_xm_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_zm_edge_scalar, brd_size, MPI_my_double_type, me_xp_zm, 10,
                xm_zp_edge_scalar, brd_size, MPI_my_double_type, me_xm_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_zm_edge_scalar[IDX_EDGE_Y(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = xm_zp_edge_scalar[IDX_EDGE_Y(i,j,k)];
        xm_zm_edge_scalar[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+BRD)];
        xm_zp_edge_scalar[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+LNZ)];
      }
 MPI_Sendrecv( xm_zm_edge_scalar, brd_size, MPI_my_double_type, me_xm_zm, 11,
               xp_zp_edge_scalar, brd_size, MPI_my_double_type, me_xp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_zp_edge_scalar, brd_size, MPI_my_double_type, me_xm_zp, 11,
               xp_zm_edge_scalar, brd_size, MPI_my_double_type, me_xp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)] = xp_zp_edge_scalar[IDX_EDGE_Y(i,j,k)];
	f[IDX(i+LNX+BRD,j,k)] = xp_zm_edge_scalar[IDX_EDGE_Y(i,j,k)];
      }
 
 
 /* along z */
 
 brd_size = BRD*BRD*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_edge_scalar[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+LNY,k)];
        xm_yp_edge_scalar[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+LNY,k)];
      }

  MPI_Sendrecv( xp_yp_edge_scalar, brd_size, MPI_my_double_type, me_xp_yp, 10,
                xm_ym_edge_scalar, brd_size, MPI_my_double_type, me_xm_ym, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_edge_scalar, brd_size, MPI_my_double_type, me_xm_yp, 10,
                xp_ym_edge_scalar, brd_size, MPI_my_double_type, me_xp_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_ym_edge_scalar[IDX_EDGE_Z(i,j,k)];
        f[IDX(i+LNX+BRD,j,k)] = xp_ym_edge_scalar[IDX_EDGE_Z(i,j,k)];
        xm_ym_edge_scalar[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+BRD,k)];
        xp_ym_edge_scalar[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+BRD,k)];
      }
 MPI_Sendrecv( xm_ym_edge_scalar, brd_size, MPI_my_double_type, me_xm_ym, 11,
               xp_yp_edge_scalar, brd_size, MPI_my_double_type, me_xp_yp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_edge_scalar, brd_size, MPI_my_double_type, me_xp_ym, 11,
               xm_yp_edge_scalar, brd_size, MPI_my_double_type, me_xm_yp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)] = xp_yp_edge_scalar[IDX_EDGE_Z(i,j,k)];
	f[IDX(i,j+LNY+BRD,k)] = xm_yp_edge_scalar[IDX_EDGE_Z(i,j,k)];
      }
 
#endif /* METHOD_EDGES_AND_CORNERS */

}/* end scalar send rcv function */

/****************************************************************************************************/ 
void sendrecv_borders_vector(vector *f){
  int i,j,k,brd_size;
  MPI_Status status1;

#ifdef NEW_SENDRECV
  //IDX(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )
     
  MPI_Sendrecv( f + IDX(LNX,0,0)   , 1, MPI_vector_plane_x, me_xp, 14,
                f + IDX(0,0,0)     , 1, MPI_vector_plane_x, me_xm, 14, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(BRD,0,0)    , 1, MPI_vector_plane_x, me_xm, 15,
                f + IDX(LNX+BRD,0,0), 1, MPI_vector_plane_x, me_xp, 15, MPI_COMM_WORLD, &status1);
  
  MPI_Sendrecv( f + IDX(0,LNY,0)  , 1, MPI_vector_plane_y, me_yp, 12,
                f + IDX(0,0,0)  , 1, MPI_vector_plane_y, me_ym, 12, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(0,BRD,0)  , 1, MPI_vector_plane_y, me_ym, 13,
                f + IDX(0,LNY+BRD,0)  , 1, MPI_vector_plane_y, me_yp, 13, MPI_COMM_WORLD, &status1);
  
  MPI_Sendrecv( f + IDX(0,0,LNZ)  , 1, MPI_vector_plane_z, me_zp, 10,
                f + IDX(0,0,0)  , 1, MPI_vector_plane_z, me_zm, 10, MPI_COMM_WORLD, &status1);
  MPI_Sendrecv( f + IDX(0,0,BRD)  , 1, MPI_vector_plane_z, me_zm, 11,
                f + IDX(0,0,LNZ+BRD)  , 1, MPI_vector_plane_z, me_zp, 11, MPI_COMM_WORLD, &status1);
  
#else /* = NEW_SENDRECV is not defined */

  /*     BRD|LNX|BRD     */
  /* Copy borders along x */
  brd_size = BRD*(LNY+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_vector[IDX_XBRD(i,j,k)] = f[IDX(i+LNX,j,k)];
      }

  MPI_Sendrecv( xp_vector, brd_size, MPI_vector_type, me_xp, 10,
                xm_vector, brd_size, MPI_vector_type, me_xm, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_vector[IDX_XBRD(i,j,k)];
        xm_vector[IDX_XBRD(i,j,k)] = f[IDX(i+BRD,j,k)];
      }
 MPI_Sendrecv( xm_vector, brd_size, MPI_vector_type, me_xm, 11,
               xp_vector, brd_size, MPI_vector_type, me_xp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k)] = xp_vector[IDX_XBRD(i,j,k)];
      }



  /* Copy borders along y */
  brd_size = BRD*(LNX+TWO_BRD)*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_vector[IDX_YBRD(i,j,k)] = f[IDX(i,j+LNY,k)];
      }

  MPI_Sendrecv( yp_vector, brd_size, MPI_vector_type, me_yp, 10,
                ym_vector, brd_size, MPI_vector_type, me_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = ym_vector[IDX_YBRD(i,j,k)];
	ym_vector[IDX_YBRD(i,j,k)] = f[IDX(i,j+BRD,k)];
      }
  
 MPI_Sendrecv( ym_vector, brd_size, MPI_vector_type, me_ym, 13,
               yp_vector, brd_size, MPI_vector_type, me_yp, 13, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
       	f[IDX(i,j+LNY+BRD,k)] = yp_vector[IDX_YBRD(i,j,k)];
      }
  

  /* Copy borders along z */
  brd_size = BRD*(LNX+TWO_BRD)*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        zp_vector[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+LNZ)];
      }

  MPI_Sendrecv( zp_vector, brd_size, MPI_vector_type, me_zp, 14,
                zm_vector, brd_size, MPI_vector_type, me_zm, 14, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = zm_vector[IDX_ZBRD(i,j,k)];
        zm_vector[IDX_ZBRD(i,j,k)] = f[IDX(i,j,k+BRD)];
      }
 MPI_Sendrecv( zm_vector, brd_size, MPI_vector_type, me_zm, 15,
               zp_vector, brd_size, MPI_vector_type, me_zp, 15, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j,k+LNZ+BRD)] = zp_vector[IDX_ZBRD(i,j,k)];
      }


#endif /* NEW_SENDRECV */

#ifdef METHOD_EDGES_AND_CORNERS

  /* First we communicate the 8 corner cubes (they are either 1x1x1 or 2x2x2 depending on BRD) */ 
  
  brd_size = BRD*BRD*BRD;

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_zp_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+LNZ)];
        xp_yp_zm_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+LNY,k+BRD)];
        xp_ym_zp_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+LNZ)];
        xm_yp_zp_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+LNZ)];
      }

  MPI_Sendrecv( xp_yp_zp_corner_vector, brd_size, MPI_vector_type, me_xp_yp_zp, 10,
                xm_ym_zm_corner_vector, brd_size, MPI_vector_type, me_xm_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_yp_zm_corner_vector, brd_size, MPI_vector_type, me_xp_yp_zm, 10,
                xm_ym_zp_corner_vector, brd_size, MPI_vector_type, me_xm_ym_zp, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_ym_zp_corner_vector, brd_size, MPI_vector_type, me_xp_ym_zp, 10,
                xm_yp_zm_corner_vector, brd_size, MPI_vector_type, me_xm_yp_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_zp_corner_vector, brd_size, MPI_vector_type, me_xm_yp_zp, 10,
                xp_ym_zm_corner_vector, brd_size, MPI_vector_type, me_xp_ym_zm, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_ym_zm_corner_vector[IDX_CORNER(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = xm_ym_zp_corner_vector[IDX_CORNER(i,j,k)];
        f[IDX(i,j+LNY+BRD,k)] = xm_yp_zm_corner_vector[IDX_CORNER(i,j,k)];
        f[IDX(i+LNX+BRD,j,k)] = xp_ym_zm_corner_vector[IDX_CORNER(i,j,k)];
        xm_ym_zm_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+BRD)];
        xm_ym_zp_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+BRD,k+LNZ)];
        xm_yp_zm_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+BRD,j+LNY,k+BRD)];
        xp_ym_zm_corner_vector[IDX_CORNER(i,j,k)] = f[IDX(i+LNX,j+BRD,k+BRD)];
      }
 MPI_Sendrecv( xm_ym_zm_corner_vector, brd_size, MPI_vector_type, me_xm_ym_zm, 11,
               xp_yp_zp_corner_vector, brd_size, MPI_vector_type, me_xp_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_ym_zp_corner_vector, brd_size, MPI_vector_type, me_xm_ym_zp, 11,
               xp_yp_zm_corner_vector, brd_size, MPI_vector_type, me_xp_yp_zm, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_yp_zm_corner_vector, brd_size, MPI_vector_type, me_xm_yp_zm, 11,
               xp_ym_zp_corner_vector, brd_size, MPI_vector_type, me_xp_ym_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_zm_corner_vector, brd_size, MPI_vector_type, me_xp_ym_zm, 11,
               xm_yp_zp_corner_vector, brd_size, MPI_vector_type, me_xm_yp_zp, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k+LNZ+BRD)] = xp_yp_zp_corner_vector[IDX_CORNER(i,j,k)];
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)] = xp_yp_zm_corner_vector[IDX_CORNER(i,j,k)];
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)] = xp_ym_zp_corner_vector[IDX_CORNER(i,j,k)];
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)] = xm_yp_zp_corner_vector[IDX_CORNER(i,j,k)];
      }
 

 /* Then we communicate the 12 edges  */
 
 /* along x */
 
 brd_size = BRD*BRD*(LNX+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
        yp_zp_edge_vector[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+LNZ)];
        yp_zm_edge_vector[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+LNY,k+BRD)];
      }

  MPI_Sendrecv( yp_zp_edge_vector, brd_size, MPI_vector_type, me_yp_zp, 10,
                ym_zm_edge_vector, brd_size, MPI_vector_type, me_ym_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( yp_zm_edge_vector, brd_size, MPI_vector_type, me_yp_zm, 10,
                ym_zp_edge_vector, brd_size, MPI_vector_type, me_ym_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++) {
        f[IDX(i,j,k)] = ym_zm_edge_vector[IDX_EDGE_X(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = ym_zp_edge_vector[IDX_EDGE_X(i,j,k)];
        ym_zm_edge_vector[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+BRD)];
        ym_zp_edge_vector[IDX_EDGE_X(i,j,k)] = f[IDX(i,j+BRD,k+LNZ)];
      }
 MPI_Sendrecv( ym_zm_edge_vector, brd_size, MPI_vector_type, me_ym_zm, 11,
               yp_zp_edge_vector, brd_size, MPI_vector_type, me_yp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( ym_zp_edge_vector, brd_size, MPI_vector_type, me_ym_zp, 11,
               yp_zm_edge_vector, brd_size, MPI_vector_type, me_yp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
	f[IDX(i,j+LNY+BRD,k+LNZ+BRD)] = yp_zp_edge_vector[IDX_EDGE_X(i,j,k)];
	f[IDX(i,j+LNY+BRD,k)] = yp_zm_edge_vector[IDX_EDGE_X(i,j,k)];
      }
 
 /* along y */
 
 brd_size = BRD*BRD*(LNY+TWO_BRD);

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_zp_edge_vector[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+LNZ)];
        xp_zm_edge_vector[IDX_EDGE_Y(i,j,k)] = f[IDX(i+LNX,j,k+BRD)];
      }

  MPI_Sendrecv( xp_zp_edge_vector, brd_size, MPI_vector_type, me_xp_zp, 10,
                xm_zm_edge_vector, brd_size, MPI_vector_type, me_xm_zm, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xp_zm_edge_vector, brd_size, MPI_vector_type, me_xp_zm, 10,
                xm_zp_edge_vector, brd_size, MPI_vector_type, me_xm_zp, 10, MPI_COMM_WORLD, &status1); 

  for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_zm_edge_vector[IDX_EDGE_Y(i,j,k)];
        f[IDX(i,j,k+LNZ+BRD)] = xm_zp_edge_vector[IDX_EDGE_Y(i,j,k)];
        xm_zm_edge_vector[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+BRD)];
        xm_zp_edge_vector[IDX_EDGE_Y(i,j,k)] = f[IDX(i+BRD,j,k+LNZ)];
      }
 MPI_Sendrecv( xm_zm_edge_vector, brd_size, MPI_vector_type, me_xm_zm, 11,
               xp_zp_edge_vector, brd_size, MPI_vector_type, me_xp_zp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xm_zp_edge_vector, brd_size, MPI_vector_type, me_xm_zp, 11,
               xp_zm_edge_vector, brd_size, MPI_vector_type, me_xp_zm, 11, MPI_COMM_WORLD, &status1);

 for(k=0;k<BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j,k+LNZ+BRD)] = xp_zp_edge_vector[IDX_EDGE_Y(i,j,k)];
	f[IDX(i+LNX+BRD,j,k)] = xp_zm_edge_vector[IDX_EDGE_Y(i,j,k)];
      }
 
 
 /* along z */
 
 brd_size = BRD*BRD*(LNZ+TWO_BRD);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
        xp_yp_edge_vector[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+LNY,k)];
        xm_yp_edge_vector[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+LNY,k)];
      }

  MPI_Sendrecv( xp_yp_edge_vector, brd_size, MPI_vector_type, me_xp_yp, 10,
                xm_ym_edge_vector, brd_size, MPI_vector_type, me_xm_ym, 10, MPI_COMM_WORLD, &status1); 
  MPI_Sendrecv( xm_yp_edge_vector, brd_size, MPI_vector_type, me_xm_yp, 10,
                xp_ym_edge_vector, brd_size, MPI_vector_type, me_xp_ym, 10, MPI_COMM_WORLD, &status1); 

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++) {
        f[IDX(i,j,k)] = xm_ym_edge_vector[IDX_EDGE_Z(i,j,k)];
        f[IDX(i+LNX+BRD,j,k)] = xp_ym_edge_vector[IDX_EDGE_Z(i,j,k)];
        xm_ym_edge_vector[IDX_EDGE_Z(i,j,k)] = f[IDX(i+BRD,j+BRD,k)];
        xp_ym_edge_vector[IDX_EDGE_Z(i,j,k)] = f[IDX(i+LNX,j+BRD,k)];
      }
 MPI_Sendrecv( xm_ym_edge_vector, brd_size, MPI_vector_type, me_xm_ym, 11,
               xp_yp_edge_vector, brd_size, MPI_vector_type, me_xp_yp, 11, MPI_COMM_WORLD, &status1);
 MPI_Sendrecv( xp_ym_edge_vector, brd_size, MPI_vector_type, me_xp_ym, 11,
               xm_yp_edge_vector, brd_size, MPI_vector_type, me_xm_yp, 11, MPI_COMM_WORLD, &status1);

 for(k=BRD;k<LNZ+BRD;k++)
    for(j=0;j<BRD;j++)
      for(i=0;i<BRD;i++){ 
	f[IDX(i+LNX+BRD,j+LNY+BRD,k)] = xp_yp_edge_vector[IDX_EDGE_Z(i,j,k)];
	f[IDX(i,j+LNY+BRD,k)] = xm_yp_edge_vector[IDX_EDGE_Z(i,j,k)];
      }
 
#endif /* METHOD_EDGES_AND_CORNERS */


}/* end vector send rcv function */

#endif


/**************************************************************************************************/
/* Here begin the function to compute the interpolation coefficients */

void compute_interpolation_coefficients(){

  int i,j,k;
  vector xp,xm,yp,ym,zp,zm;
  vector w_xp,w_xm,w_yp,w_ym,w_zp,w_zm;
  vector dxp,dxm,dyp,dym,dzp,dzm;
  vector P0, P1, P2, P3, P4, P5, P6, P7, P8;
  vector edge;

  /* for quick */
  vector XPF1, XPF2, XPB1, XPB2, XMF1, XMF2, XMB1, XMB2;
  vector YPF1, YPF2, YPB1, YPB2, YMF1, YMF2, YMB1, YMB2;   
  vector ZPF1, ZPF2, ZPB1, ZPB2, ZMF1, ZMF2, ZMB1, ZMB2;

#if (defined METHOD_STREAMING || defined METHOD_UPWIND)
  //#ifdef METHOD_STREAMING
  //#ifdef METHOD_UPWIND
	/* exchange centers */
	sendrecv_borders_vector(center_V);

        /*
	for (j = BRD; j < LNY + BRD ; j++) 
	for (k = BRD; k < LNZ + BRD ; k++){
        */
	for (j = 0; j < LNY + TWO_BRD ; j++) 
	for (k = 0; k < LNZ + TWO_BRD ; k++){
	if(LNX_START == 0){
	  center_V[IDX(0, j, k)].x = - center_V[IDX(1, j , k)].x;
	}	
	if(LNX_END == NX){
	  center_V[IDX(LNX+TWO_BRD-1, j, k)].x =  (property.SX - center_V[IDX(LNX+BRD-1, j, k)].x)+property.SX;
	}
	}
	
	/*
       	for (i = BRD; i < LNX + BRD; i++) 
	for (k = BRD; k < LNZ + BRD; k++) {
	*/
       	for (i = 0; i < LNX + TWO_BRD; i++) 
	for (k = 0; k < LNZ + TWO_BRD; k++) {
	if(LNY_START == 0){
	  center_V[IDX(i, 0, k)].y = - center_V[IDX(i, 1 , k)].y;
	}	
	if(LNY_END == NY){
	  center_V[IDX(i, LNY+TWO_BRD-1, k)].y =  (property.SY - center_V[IDX(i, LNY+BRD-1, k)].y)+property.SY;
	}
	}

	/*
	for (i = BRD; i < LNX + BRD ; i++) 
	for (j = BRD; j < LNY + BRD ; j++){ 
	*/
	for (i = 0; i < LNX + TWO_BRD ; i++) 
	for (j = 0; j < LNY + TWO_BRD ; j++){ 
	if(LNZ_START == 0){
	  center_V[IDX(i, j, 0)].z = - center_V[IDX(i, j, 1)].z;
	}	
	if(LNZ_END == NZ){
	  center_V[IDX(i, j, LNZ+TWO_BRD-1)].z =  (property.SZ - center_V[IDX(i, j, LNZ+BRD-1)].z)+property.SZ;
	}
	}

#ifdef DEBUG_HARD
	
	fprintf(stderr,"\n TWO_BRD %d LNX %d , LNY %d LNZ %d\n", TWO_BRD, LNX+TWO_BRD, LNY+TWO_BRD, LNZ+TWO_BRD);
	/*	
	for (i = 0; i < LNX + TWO_BRD; i++) 
	  for (j = 0; j < LNY + TWO_BRD; j++) 
	    for (k = 0; k < LNZ + TWO_BRD; k++) {
		
	      //fprintf(stderr,"i j k %ld %ld %ld\n",i,j,k);
	        fprintf(stderr,"%lf %lf %lf\n", center_V[IDX(i,j,k)].x , center_V[IDX(i,j,k)].y , center_V[IDX(i,j,k)].z );

			}	
	*/	
#endif
	
#endif

#ifdef METHOD_UPWIND_SKEW
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
w_xp = average_4vectors(P1, P3, P5, P7);
w_xm = average_4vectors(P0, P2, P4, P6); 
w_yp = average_4vectors(P2, P3, P6, P7); 
w_ym = average_4vectors(P0, P1, P4, P5);
w_zp = average_4vectors(P4, P5, P6, P7);
w_zm = average_4vectors(P0, P1, P2, P3); 

/* from inside to outside */
/* we find the distance of the centers from such 6 surface centers */
 xp = vector_difference(w_xp, center_V[IDX(i, j, k)]);
 xm = vector_difference(w_xm, center_V[IDX(i, j, k)]);
 yp = vector_difference(w_yp, center_V[IDX(i, j, k)]);
 ym = vector_difference(w_ym, center_V[IDX(i, j, k)]);
 zp = vector_difference(w_zp, center_V[IDX(i, j, k)]);
 zm = vector_difference(w_zm, center_V[IDX(i, j, k)]);
 /* we find the distance between the two neighboring nodes */ 
 dxp = vector_difference(center_V[IDX(i+1, j, k)], center_V[IDX(i, j, k)]);
 dxm = vector_difference(center_V[IDX(i-1, j, k)], center_V[IDX(i, j, k)]);
 dyp = vector_difference(center_V[IDX(i, j+1, k)], center_V[IDX(i, j, k)]);
 dym = vector_difference(center_V[IDX(i, j-1, k)], center_V[IDX(i, j, k)]);
 dzp = vector_difference(center_V[IDX(i, j, k+1)], center_V[IDX(i, j, k)]);
 dzm = vector_difference(center_V[IDX(i, j, k-1)], center_V[IDX(i, j, k)]);

 interp_xp[IDX(i,j,k)] = fabs(xp.x);
 interp_xm[IDX(i,j,k)] = fabs(xm.x);			      
 interp_yp[IDX(i,j,k)] = fabs(yp.y);
 interp_ym[IDX(i,j,k)] = fabs(ym.y);			      
 interp_zp[IDX(i,j,k)] = fabs(zp.z);
 interp_zm[IDX(i,j,k)] = fabs(zm.z);

 interp3_xp[IDX(i,j,k)] = fabs(dxp.x);
 interp3_xm[IDX(i,j,k)] = fabs(dxm.x);
 interp3_yp[IDX(i,j,k)] = fabs(dyp.y);
 interp3_ym[IDX(i,j,k)] = fabs(dym.y);
 interp3_zp[IDX(i,j,k)] = fabs(dzp.z);
 interp3_zm[IDX(i,j,k)] = fabs(dzm.z);

			}/*i,j,k*/

	/*  send receive the coefficients */	
	sendrecv_borders_scalar(interp_xp);
	sendrecv_borders_scalar(interp_xm);
	sendrecv_borders_scalar(interp_yp);
	sendrecv_borders_scalar(interp_ym);
	sendrecv_borders_scalar(interp_zp);
	sendrecv_borders_scalar(interp_zm);

	sendrecv_borders_scalar(interp3_xp);
	sendrecv_borders_scalar(interp3_xm);
	sendrecv_borders_scalar(interp3_yp);
	sendrecv_borders_scalar(interp3_ym);
	sendrecv_borders_scalar(interp3_zp);
	sendrecv_borders_scalar(interp3_zm);
#endif


#ifdef METHOD_CENTERED
	/* exchange centers */
	sendrecv_borders_vector(center_V);

	for (j = BRD; j < LNY + BRD ; j++) 
	for (k = BRD; k < LNZ + BRD ; k++){
	if(LNX_START == 0){
	  center_V[IDX(0, j, k)].x = - center_V[IDX(1, j , k)].x;
	}	
	if(LNX_END == NX){
	  center_V[IDX(LNX+TWO_BRD-1, j, k)].x =  (property.SX - center_V[IDX(LNX+BRD-1, j, k)].x)+property.SX;
	}
	}
	

       	for (i = BRD; i < LNX + BRD; i++) 
	for (k = BRD; k < LNZ + BRD; k++) {
	if(LNY_START == 0){
	  center_V[IDX(i, 0, k)].y = - center_V[IDX(i, 1 , k)].y;
	}	
	if(LNY_END == NY){
	  center_V[IDX(i, LNY+TWO_BRD-1, k)].y =  (property.SY - center_V[IDX(i, LNY+BRD-1, k)].y)+property.SY;
	}
	}

	
	for (i = BRD; i < LNX + BRD ; i++) 
	for (j = BRD; j < LNY + BRD ; j++){ 
	if(LNZ_START == 0){
	  center_V[IDX(i, j, 0)].z = - center_V[IDX(i, j, 1)].z;
	}	
	if(LNZ_END == NZ){
	  center_V[IDX(i, j, LNZ+TWO_BRD-1)].z =  (property.SZ - center_V[IDX(i, j, LNZ+BRD-1)].z)+property.SZ;
	}
	}


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

/* we find the distance of the centers from such 6 surface centers */
 xp = vector_difference(xp, center_V[IDX(i, j, k)]);
 xm = vector_difference(center_V[IDX(i, j, k)],xm);
 yp = vector_difference(yp, center_V[IDX(i, j, k)]);
 ym = vector_difference(center_V[IDX(i, j, k)],ym);
 zp = vector_difference(zp, center_V[IDX(i, j, k)]);
 zm = vector_difference(center_V[IDX(i, j, k)],zm);
 /* we find the distance between the two neighboring nodes */ 
 dxp = vector_difference(center_V[IDX(i+1, j, k)], center_V[IDX(i, j, k)]);
 dxm = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i-1, j, k)]);
 dyp = vector_difference(center_V[IDX(i, j+1, k)], center_V[IDX(i, j, k)]);
 dym = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j-1, k)]);
 dzp = vector_difference(center_V[IDX(i, j, k+1)], center_V[IDX(i, j, k)]);
 dzm = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j, k-1)]);

/* now correct for the boundary nodes */
/*
 if(i == LNX+BRD-1 && LNX_END == NX){ 
   edge.x=property.SX;edge.y=0;edge.z=0;
   dxp = vector_difference(edge, center_V[IDX(i, j, k)]); 
   dxp=vector_scale(2.0,dxp);
 }

 if(i == BRD && LNX_START == 0){
  edge.x=0;edge.y=0;edge.z=0;
  dxm = vector_difference(center_V[IDX(i, j, k)],edge);
  dxm=vector_scale(2.0,dxm);
 }
 
 if(j == LNY+BRD-1 && LNY_END == NY){ 
   edge.x=0;edge.y=property.SY;edge.z=0;
   dyp = vector_difference(edge, center_V[IDX(i, j, k)]); 
   dyp=vector_scale(2.0,dyp);
 }

 if(j == BRD && LNY_START == 0){
  edge.x=0;edge.y=0;edge.z=0;
  dym = vector_difference(center_V[IDX(i, j, k)],edge);
  dym=vector_scale(2.0,dym);
 }

 if(k == LNZ+BRD-1 && LNZ_END == NZ){ 
   edge.x=0;edge.y=0;edge.z=property.SZ;
   dzp = vector_difference(edge, center_V[IDX(i, j, k)]); 
   dzp=vector_scale(2.0,dzp);
 }

 if(k == BRD && LNZ_START == 0){
  edge.x=0;edge.y=0;edge.z=0;
  dzm = vector_difference(center_V[IDX(i, j, k)],edge);
  dzm=vector_scale(2.0,dzm);
 }
*/
 interp_xp[IDX(i,j,k)] = xp.x/dxp.x;
 interp_xm[IDX(i,j,k)] = xm.x/dxm.x;
 interp_yp[IDX(i,j,k)] = yp.y/dyp.y;
 interp_ym[IDX(i,j,k)] = ym.y/dym.y;
 interp_zp[IDX(i,j,k)] = zp.z/dzp.z;
 interp_zm[IDX(i,j,k)] = zm.z/dzm.z;

			}/*i,j,k*/

	/*  send receive the coefficients */	
	sendrecv_borders_scalar(interp_xp);
	sendrecv_borders_scalar(interp_xm);
	sendrecv_borders_scalar(interp_yp);
	sendrecv_borders_scalar(interp_ym);
	sendrecv_borders_scalar(interp_zp);
	sendrecv_borders_scalar(interp_zm);
	
	
#ifdef BEBUG_HARD
	for (i = BRD; i < LNX + BRD; i++) 
		for (j = BRD; j < LNY + BRD; j++) 
			for (k = BRD; k < LNZ + BRD; k++) {

			  fprintf(stderr,"ZZZZ me %d i j k %d %d %d | %e %e\n", me, i,j,k,interp_yp[IDX(i,j,k)], interp_ym[IDX(i,j,k)]);
}/* i,j,k*/
#endif

#endif


#ifdef METHOD_MYQUICK

	//vector w_xp,w_xm,w_yp,w_ym,w_zp,w_zm;
  vector xp2,xm2,yp2,ym2,zp2,zm2;
  vector dxp2,dxm2,dyp2,dym2,dzp2,dzm2;
  vector xp3,xm3,yp3,ym3,zp3,zm3;
  vector dxp3,dxm3,dyp3,dym3,dzp3,dzm3;

	/* exchange centers */
  	sendrecv_borders_vector(center_V);
	
	for (j = BRD; j < LNY + BRD ; j++) 
	for (k = BRD; k < LNZ + BRD ; k++){
	if(LNX_START == 0){
	  center_V[IDX(1, j, k)].x = - center_V[IDX(2, j , k)].x;
	  center_V[IDX(0, j, k)].x = - center_V[IDX(3, j , k)].x;
	  if(NX==1) center_V[IDX(0, j, k)].x = center_V[IDX(1, j, k)].x - 2.0*center_V[IDX(2, j , k)].x;
	}	
	if(LNX_END == NX){
	  center_V[IDX(LNX+TWO_BRD-2, j, k)].x =  (property.SX - center_V[IDX(LNX+BRD-1, j, k)].x)+property.SX;
	  center_V[IDX(LNX+TWO_BRD-1, j, k)].x =  (property.SX - center_V[IDX(LNX+BRD-2, j, k)].x)+property.SX;
	}
	}
	

       	for (i = BRD; i < LNX + BRD; i++) 
	for (k = BRD; k < LNZ + BRD; k++) {
	if(LNY_START == 0){
	  center_V[IDX(i, 1, k)].y = - center_V[IDX(i, 2 , k)].y;
	  center_V[IDX(i, 0, k)].y = - center_V[IDX(i, 3 , k)].y;
	  if(NY==1) center_V[IDX(i, 0, k)].y = center_V[IDX(i, 1, k)].y - 2.0*center_V[IDX(i, 2 , k)].y;
	}	
	if(LNY_END == NY){
	  center_V[IDX(i, LNY+TWO_BRD-2, k)].y =  (property.SY - center_V[IDX(i, LNY+BRD-1, k)].y)+property.SY;
	  center_V[IDX(i, LNY+TWO_BRD-1, k)].y =  (property.SY - center_V[IDX(i, LNY+BRD-2, k)].y)+property.SY;
	}
	}

	
	for (i = BRD; i < LNX + BRD ; i++) 
	for (j = BRD; j < LNY + BRD ; j++){ 
	if(LNZ_START == 0){
	  center_V[IDX(i, j, 1)].z = - center_V[IDX(i, j, 2)].z;
	  center_V[IDX(i, j, 0)].z = - center_V[IDX(i, j, 3)].z;
	  if(NZ==1) center_V[IDX(i, j, 0)].z = center_V[IDX(i, j, 1)].z - 2.0*center_V[IDX(i, j , 2)].z;
	}	
	if(LNZ_END == NZ){
	  center_V[IDX(i, j, LNZ+TWO_BRD-2)].z =  (property.SZ - center_V[IDX(i, j, LNZ+BRD-1)].z)+property.SZ;
	  center_V[IDX(i, j, LNZ+TWO_BRD-1)].z =  (property.SZ - center_V[IDX(i, j, LNZ+BRD-2)].z)+property.SZ;
	}
	}
	

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
w_xp = average_4vectors(P1, P3, P5, P7);
w_xm = average_4vectors(P0, P2, P4, P6); 
w_yp = average_4vectors(P2, P3, P6, P7); 
w_ym = average_4vectors(P0, P1, P4, P5);
w_zp = average_4vectors(P4, P5, P6, P7);
w_zm = average_4vectors(P0, P1, P2, P3); 


/*******************************************************************/
/*  Part 1: from inside to the outside */ 

/* we find the distance of the centers from such 6 surface centers */
/*  x - xu */
 xp = vector_difference(w_xp, center_V[IDX(i, j, k)]);
 xm = vector_difference(w_xm, center_V[IDX(i, j, k)]);
 yp = vector_difference(w_yp, center_V[IDX(i, j, k)]);
 ym = vector_difference(w_ym, center_V[IDX(i, j, k)]);
 zp = vector_difference(w_zp, center_V[IDX(i, j, k)]);
 zm = vector_difference(w_zm, center_V[IDX(i, j, k)]);
/* we find the distance of the centers from such 6 surface centers */
/* x - xuu */
 xp2 = vector_difference(w_xp, center_V[IDX(i-1, j, k)]);
 xm2 = vector_difference(w_xm, center_V[IDX(i+1, j, k)]);
 yp2 = vector_difference(w_yp, center_V[IDX(i, j-1, k)]);
 ym2 = vector_difference(w_ym, center_V[IDX(i, j+1, k)]);
 zp2 = vector_difference(w_zp, center_V[IDX(i, j, k-1)]);
 zm2 = vector_difference(w_zm, center_V[IDX(i, j, k+1)]);
 /* we find the distance between the two neighboring nodes */ 
 /* xd - xu */
 dxp = vector_difference(center_V[IDX(i+1, j, k)], center_V[IDX(i, j, k)]);
 dxm = vector_difference(center_V[IDX(i-1, j, k)], center_V[IDX(i, j, k)]);
 dyp = vector_difference(center_V[IDX(i, j+1, k)], center_V[IDX(i, j, k)]);
 dym = vector_difference(center_V[IDX(i, j-1, k)], center_V[IDX(i, j, k)]);
 dzp = vector_difference(center_V[IDX(i, j, k+1)], center_V[IDX(i, j, k)]);
 dzm = vector_difference(center_V[IDX(i, j, k-1)], center_V[IDX(i, j, k)]);
/* we find the distance between the two next neighboring nodes */ 
/* xd - xuu */
 dxp2 = vector_difference(center_V[IDX(i+1, j, k)], center_V[IDX(i-1, j, k)]);
 dxm2 = vector_difference(center_V[IDX(i+1, j, k)], center_V[IDX(i-1, j, k)]);
 dyp2 = vector_difference(center_V[IDX(i, j+1, k)], center_V[IDX(i, j-1, k)]);
 dym2 = vector_difference(center_V[IDX(i, j+1, k)], center_V[IDX(i, j-1, k)]);
 dzp2 = vector_difference(center_V[IDX(i, j, k+1)], center_V[IDX(i, j, k-1)]);
 dzm2 = vector_difference(center_V[IDX(i, j, k+1)], center_V[IDX(i, j, k-1)]);
 /* xd - x */
 xp3 = vector_difference(center_V[IDX(i+1, j, k)],w_xp);
 xm3 = vector_difference(center_V[IDX(i-1, j, k)],w_xm);
 yp3 = vector_difference(center_V[IDX(i, j+1, k)],w_yp);
 ym3 = vector_difference(center_V[IDX(i, j-1, k)],w_ym);
 zp3 = vector_difference(center_V[IDX(i, j, k+1)],w_zp);
 zm3 = vector_difference(center_V[IDX(i, j, k-1)],w_zm);
 /* xu - xuu */
 dxp3 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i-1, j, k)]);
 dxm3 = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i+1, j, k)]);
 dyp3 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j-1, k)]);
 dym3 = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j+1, k)]);
 dzp3 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j, k-1)]);
 dzm3 = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j, k+1)]);

 /* This is UPWIND g1 = (x-xu)*(x-xuu) / (xd-xu)*(xd-xuu)         */
 interp_xp[IDX(i,j,k)] = fabs( (xp.x*xp2.x)/(dxp.x*dxp2.x) );
 interp_xm[IDX(i,j,k)] = fabs( (xm.x*xm2.x)/(dxm.x*dxm2.x) );
 interp_yp[IDX(i,j,k)] = fabs( (yp.y*yp2.y)/(dyp.y*dyp2.y) );
 interp_ym[IDX(i,j,k)] = fabs( (ym.y*ym2.y)/(dym.y*dym2.y) );
 interp_zp[IDX(i,j,k)] = fabs( (zp.z*zp2.z)/(dzp.z*dzp2.z) );
 interp_zm[IDX(i,j,k)] = fabs( (zm.z*zm2.z)/(dzm.z*dzm2.z) );

 /* This is UPWIND g2 = (x-xu)*(xd-x) / (xu-xuu)*(xd-xuu)         */
 interp2_xp[IDX(i,j,k)] = fabs( (xp.x*xp3.x)/(dxp3.x*dxp2.x) );
 interp2_xm[IDX(i,j,k)] = fabs( (xm.x*xm3.x)/(dxm3.x*dxm2.x) );
 interp2_yp[IDX(i,j,k)] = fabs( (yp.y*yp3.y)/(dyp3.y*dyp2.y) );
 interp2_ym[IDX(i,j,k)] = fabs( (ym.y*ym3.y)/(dym3.y*dym2.y) );
 interp2_zp[IDX(i,j,k)] = fabs( (zp.z*zp3.z)/(dzp3.z*dzp2.z) );
 interp2_zm[IDX(i,j,k)] = fabs( (zm.z*zm3.z)/(dzm3.z*dzm2.z) );


 /*******************************************************************/
/*  Part 2: from outside to the inside */ 

/* we find the distance of the centers from such 6 surface centers */
/*  x - xu */
 xp = vector_difference(w_xp, center_V[IDX(i+1, j, k)]);
 xm = vector_difference(w_xm, center_V[IDX(i-1, j, k)]);
 yp = vector_difference(w_yp, center_V[IDX(i, j+1, k)]);
 ym = vector_difference(w_ym, center_V[IDX(i, j-1, k)]);
 zp = vector_difference(w_zp, center_V[IDX(i, j, k+1)]);
 zm = vector_difference(w_zm, center_V[IDX(i, j, k-1)]);
/* we find the distance of the centers from such 6 surface centers */
/* x - xuu */
 xp2 = vector_difference(w_xp, center_V[IDX(i+2, j, k)]);
 xm2 = vector_difference(w_xm, center_V[IDX(i-2, j, k)]);
 yp2 = vector_difference(w_yp, center_V[IDX(i, j+2, k)]);
 ym2 = vector_difference(w_ym, center_V[IDX(i, j-2, k)]);
 zp2 = vector_difference(w_zp, center_V[IDX(i, j, k+2)]);
 zm2 = vector_difference(w_zm, center_V[IDX(i, j, k-2)]);
 /* we find the distance between the two neighboring nodes */ 
 /* xd - xu */
 dxp = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i+1, j, k)]);
 dxm = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i-1, j, k)]);
 dyp = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j+1, k)]);
 dym = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j-1, k)]);
 dzp = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j, k+1)]);
 dzm = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j, k-1)]);
/* we find the distance between the two next neighboring nodes */ 
/* xd - xuu */
 dxp2 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i+2, j, k)]);
 dxm2 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i-2, j, k)]);
 dyp2 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j+2, k)]);
 dym2 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j-2, k)]);
 dzp2 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j, k+2)]);
 dzm2 = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j, k-2)]);
 /* xd - x */
 xp3 = vector_difference(center_V[IDX(i, j, k)],w_xp);
 xm3 = vector_difference(center_V[IDX(i, j, k)],w_xm);
 yp3 = vector_difference(center_V[IDX(i, j, k)],w_yp);
 ym3 = vector_difference(center_V[IDX(i, j, k)],w_ym);
 zp3 = vector_difference(center_V[IDX(i, j, k)],w_zp);
 zm3 = vector_difference(center_V[IDX(i, j, k)],w_zm);
 /* xu - xuu */
 dxp3 = vector_difference(center_V[IDX(i+1, j, k)], center_V[IDX(i+2, j, k)]);
 dxm3 = vector_difference(center_V[IDX(i-1, j, k)],center_V[IDX(i-2, j, k)]);
 dyp3 = vector_difference(center_V[IDX(i, j+1, k)], center_V[IDX(i, j+2, k)]);
 dym3 = vector_difference(center_V[IDX(i, j-1, k)],center_V[IDX(i, j-2, k)]);
 dzp3 = vector_difference(center_V[IDX(i, j, k+1)], center_V[IDX(i, j, k+2)]);
 dzm3 = vector_difference(center_V[IDX(i, j, k-1)],center_V[IDX(i, j, k-2)]);

 /* This is DOWNWIND g3 = (x-xu)*(x-xuu) / (xd-xu)*(xd-xuu)        */
 interp3_xp[IDX(i,j,k)] = (xp.x*xp2.x)/(dxp.x*dxp2.x);
 interp3_xm[IDX(i,j,k)] = (xm.x*xm2.x)/(dxm.x*dxm2.x);
 interp3_yp[IDX(i,j,k)] = (yp.y*yp2.y)/(dyp.y*dyp2.y);
 interp3_ym[IDX(i,j,k)] = (ym.y*ym2.y)/(dym.y*dym2.y);
 interp3_zp[IDX(i,j,k)] = (zp.z*zp2.z)/(dzp.z*dzp2.z);
 interp3_zm[IDX(i,j,k)] = (zm.z*zm2.z)/(dzm.z*dzm2.z);

 /* This is DOWNWIND g4 = (x-xu)*(xd-x) / (xu-xuu)*(xd-xuu)        */   
 interp4_xp[IDX(i,j,k)] = (xp.x*xp3.x)/(dxp3.x*dxp2.x);
 interp4_xm[IDX(i,j,k)] = (xm.x*xm3.x)/(dxm3.x*dxm2.x);
 interp4_yp[IDX(i,j,k)] = (yp.y*yp3.y)/(dyp3.y*dyp2.y);
 interp4_ym[IDX(i,j,k)] = (ym.y*ym3.y)/(dym3.y*dym2.y);
 interp4_zp[IDX(i,j,k)] = (zp.z*zp3.z)/(dzp3.z*dzp2.z);
 interp4_zm[IDX(i,j,k)] = (zm.z*zm3.z)/(dzm3.z*dzm2.z);


			}/*i,j,k*/

	/*  send receive the coefficients */	
	sendrecv_borders_scalar(interp_xp);
	sendrecv_borders_scalar(interp_xm);
	sendrecv_borders_scalar(interp_yp);
	sendrecv_borders_scalar(interp_ym);
	sendrecv_borders_scalar(interp_zp);
	sendrecv_borders_scalar(interp_zm);

	sendrecv_borders_scalar(interp2_xp);
	sendrecv_borders_scalar(interp2_xm);
	sendrecv_borders_scalar(interp2_yp);
	sendrecv_borders_scalar(interp2_ym);
	sendrecv_borders_scalar(interp2_zp);
	sendrecv_borders_scalar(interp2_zm);

	sendrecv_borders_scalar(interp3_xp);
	sendrecv_borders_scalar(interp3_xm);
	sendrecv_borders_scalar(interp3_yp);
	sendrecv_borders_scalar(interp3_ym);
	sendrecv_borders_scalar(interp3_zp);
	sendrecv_borders_scalar(interp3_zm);

	sendrecv_borders_scalar(interp4_xp);
	sendrecv_borders_scalar(interp4_xm);
	sendrecv_borders_scalar(interp4_yp);
	sendrecv_borders_scalar(interp4_ym);
	sendrecv_borders_scalar(interp4_zp);
	sendrecv_borders_scalar(interp4_zm);
	
	
	#ifdef BEBUG_HARD
	for (i = BRD; i < LNX + BRD; i++) 
		for (j = BRD; j < LNY + BRD; j++) 
			for (k = BRD; k < LNZ + BRD; k++) {

fprintf(stderr,"ZZZZ me %d i j k %d %d %d | %f %f %f %f\n", me, i,j,k,
	//center_V[IDX(i,j,k)].x , center_V[IDX(i,j,k)].y , center_V[IDX(i,j,k)].z);
	//interp_xp[IDX(i,j,k)], interp2_xp[IDX(i,j,k)], interp3_xp[IDX(i,j,k)], interp4_xp[IDX(i,j,k)]);
	//interp_xm[IDX(i,j,k)], interp2_xm[IDX(i,j,k)], interp3_xm[IDX(i,j,k)], interp4_xm[IDX(i,j,k)]);
       	 interp_yp[IDX(i,j,k)], interp2_yp[IDX(i,j,k)], interp3_yp[IDX(i,j,k)], interp4_yp[IDX(i,j,k)]);
	//interp_ym[IDX(i,j,k)], interp2_ym[IDX(i,j,k)], interp3_ym[IDX(i,j,k)], interp4_ym[IDX(i,j,k)]);
        //interp_zp[IDX(i,j,k)], interp2_zp[IDX(i,j,k)], interp3_zp[IDX(i,j,k)], interp4_zp[IDX(i,j,k)]);
        //interp_zm[IDX(i,j,k)], interp2_zm[IDX(i,j,k)], interp3_zm[IDX(i,j,k)], interp4_zm[IDX(i,j,k)]);
}/* i,j,k*/
	//exit(1);
	#endif

#endif

	/***** Linear Upwind ****/
#ifdef METHOD_UPWIND_LINEAR

  //vector w_xp,w_xm,w_yp,w_ym,w_zp,w_zm;
  //vector xp2,xm2,yp2,ym2,zp2,zm2;
  //vector dxp2,dxm2,dyp2,dym2,dzp2,dzm2;
  //vector xp3,xm3,yp3,ym3,zp3,zm3;
  //vector dxp3,dxm3,dyp3,dym3,dzp3,dzm3;

	/* exchange centers */
  	sendrecv_borders_vector(center_V);
	
	for (j = BRD; j < LNY + BRD ; j++) 
	for (k = BRD; k < LNZ + BRD ; k++){
	if(LNX_START == 0){
	  center_V[IDX(1, j, k)].x = - center_V[IDX(2, j , k)].x;
	  center_V[IDX(0, j, k)].x = - center_V[IDX(3, j , k)].x;
	  if(NX==1) center_V[IDX(0, j, k)].x = center_V[IDX(1, j, k)].x - 2.0*center_V[IDX(2, j , k)].x;
	}	
	if(LNX_END == NX){
	  center_V[IDX(LNX+TWO_BRD-2, j, k)].x =  (property.SX - center_V[IDX(LNX+BRD-1, j, k)].x)+property.SX;
	  center_V[IDX(LNX+TWO_BRD-1, j, k)].x =  (property.SX - center_V[IDX(LNX+BRD-2, j, k)].x)+property.SX;
	}
	}
	

       	for (i = BRD; i < LNX + BRD; i++) 
	for (k = BRD; k < LNZ + BRD; k++) {
	if(LNY_START == 0){
	  center_V[IDX(i, 1, k)].y = - center_V[IDX(i, 2 , k)].y;
	  center_V[IDX(i, 0, k)].y = - center_V[IDX(i, 3 , k)].y;
	  if(NY==1) center_V[IDX(i, 0, k)].y = center_V[IDX(i, 1, k)].y - 2.0*center_V[IDX(i, 2 , k)].y;
	}	
	if(LNY_END == NY){
	  center_V[IDX(i, LNY+TWO_BRD-2, k)].y =  (property.SY - center_V[IDX(i, LNY+BRD-1, k)].y)+property.SY;
	  center_V[IDX(i, LNY+TWO_BRD-1, k)].y =  (property.SY - center_V[IDX(i, LNY+BRD-2, k)].y)+property.SY;
	}
	}

	
	for (i = BRD; i < LNX + BRD ; i++) 
	for (j = BRD; j < LNY + BRD ; j++){ 
	if(LNZ_START == 0){
	  center_V[IDX(i, j, 1)].z = - center_V[IDX(i, j, 2)].z;
	  center_V[IDX(i, j, 0)].z = - center_V[IDX(i, j, 3)].z;
	  if(NZ==1) center_V[IDX(i, j, 0)].z = center_V[IDX(i, j, 1)].z - 2.0*center_V[IDX(i, j , 2)].z;
	}	
	if(LNZ_END == NZ){
	  center_V[IDX(i, j, LNZ+TWO_BRD-2)].z =  (property.SZ - center_V[IDX(i, j, LNZ+BRD-1)].z)+property.SZ;
	  center_V[IDX(i, j, LNZ+TWO_BRD-1)].z =  (property.SZ - center_V[IDX(i, j, LNZ+BRD-2)].z)+property.SZ;
	}
	}
	

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
w_xp = average_4vectors(P1, P3, P5, P7);
w_xm = average_4vectors(P0, P2, P4, P6); 
w_yp = average_4vectors(P2, P3, P6, P7); 
w_ym = average_4vectors(P0, P1, P4, P5);
w_zp = average_4vectors(P4, P5, P6, P7);
w_zm = average_4vectors(P0, P1, P2, P3); 

#ifndef METHOD_UPWIND_LINEAR_IMPROVED
/*******************************************************************/
/*  Part 1: from inside to the outside */ 

/* we find the distance of the centers from such 6 surface centers */
/* x - xuu */
 xp = vector_difference(w_xp, center_V[IDX(i-1, j, k)]);
 xm = vector_difference(w_xm, center_V[IDX(i+1, j, k)]);
 yp = vector_difference(w_yp, center_V[IDX(i, j-1, k)]);
 ym = vector_difference(w_ym, center_V[IDX(i, j+1, k)]);
 zp = vector_difference(w_zp, center_V[IDX(i, j, k-1)]);
 zm = vector_difference(w_zm, center_V[IDX(i, j, k+1)]);
 /* we find the distance between the two neighboring nodes */ 
 /* xu - xuu */
 dxp = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i-1, j, k)]);
 dxm = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i+1, j, k)]);
 dyp = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j-1, k)]);
 dym = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j+1, k)]);
 dzp = vector_difference(center_V[IDX(i, j, k)], center_V[IDX(i, j, k-1)]);
 dzm = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j, k+1)]);

 /* This is UPWIND g1 = (x-xuu)/(xu-xuu)  */
 interp5_xp[IDX(i,j,k)] = fabs( xp.x/dxp.x );
 interp5_xm[IDX(i,j,k)] = fabs( xm.x/dxm.x );
 interp5_yp[IDX(i,j,k)] = fabs( yp.y/dyp.y );
 interp5_ym[IDX(i,j,k)] = fabs( ym.y/dym.y );
 interp5_zp[IDX(i,j,k)] = fabs( zp.z/dzp.z );
 interp5_zm[IDX(i,j,k)] = fabs( zm.z/dzm.z );

 /*******************************************************************/
/*  Part 2: from outside to the inside */ 

/* we find the distance of the centers from such 6 surface centers */
/* x - xuu */
 xp = vector_difference(w_xp, center_V[IDX(i+2, j, k)]);
 xm = vector_difference(w_xm, center_V[IDX(i-2, j, k)]);
 yp = vector_difference(w_yp, center_V[IDX(i, j+2, k)]);
 ym = vector_difference(w_ym, center_V[IDX(i, j-2, k)]);
 zp = vector_difference(w_zp, center_V[IDX(i, j, k+2)]);
 zm = vector_difference(w_zm, center_V[IDX(i, j, k-2)]);
 /* we find the distance between the two neighboring nodes */ 
 /* xu - xuu */
 dxp = vector_difference(center_V[IDX(i+1, j, k)],center_V[IDX(i+2, j, k)]);
 dxm = vector_difference(center_V[IDX(i-1, j, k)],center_V[IDX(i-2, j, k)]);
 dyp = vector_difference(center_V[IDX(i, j+1, k)],center_V[IDX(i, j+2, k)]);
 dym = vector_difference(center_V[IDX(i, j-1, k)],center_V[IDX(i, j-2, k)]);
 dzp = vector_difference(center_V[IDX(i, j, k+1)],center_V[IDX(i, j, k+2)]);
 dzm = vector_difference(center_V[IDX(i, j, k-1)],center_V[IDX(i, j, k-2)]);

 /* This is DOWNWIND g3 = (x-xuu)/(xu-xuu)    */
 interp6_xp[IDX(i,j,k)] =  fabs(xp.x/dxp.x);
 interp6_xm[IDX(i,j,k)] =  fabs(xm.x/dxm.x);
 interp6_yp[IDX(i,j,k)] =  fabs(yp.y/dyp.y);
 interp6_ym[IDX(i,j,k)] =  fabs(ym.y/dym.y);
 interp6_zp[IDX(i,j,k)] =  fabs(zp.z/dzp.z);
 interp6_zm[IDX(i,j,k)] =  fabs(zm.z/dzm.z);


#else
 /* In this case the coefficient are computed according to Kalyan improved method , comments here below shall be updated*/

/*  Part 1: from inside to the outside */ 

/* we find the distance of the centers from such 6 surface centers */
/* x - xuu */
 xp = vector_difference(w_xp, center_V[IDX(i, j, k)]);
 xm = vector_difference(w_xm, center_V[IDX(i, j, k)]);
 yp = vector_difference(w_yp, center_V[IDX(i, j, k)]);
 ym = vector_difference(w_ym, center_V[IDX(i, j, k)]);
 zp = vector_difference(w_zp, center_V[IDX(i, j, k)]);
 zm = vector_difference(w_zm, center_V[IDX(i, j, k)]);
 /* we find the distance between the two neighboring nodes */ 
 /* xu - xuu */
 dxp = vector_difference(center_V[IDX(i+1, j, k)], center_V[IDX(i-1, j, k)]);
 dxm = vector_difference(center_V[IDX(i+1, j, k)],center_V[IDX(i-1, j, k)]);
 dyp = vector_difference(center_V[IDX(i, j+1, k)], center_V[IDX(i, j-1, k)]);
 dym = vector_difference(center_V[IDX(i, j+1, k)],center_V[IDX(i, j-1, k)]);
 dzp = vector_difference(center_V[IDX(i, j, k+1)], center_V[IDX(i, j, k-1)]);
 dzm = vector_difference(center_V[IDX(i, j, k+1)],center_V[IDX(i, j, k-1)]);

 /* This is UPWIND g1 = (x-xuu)/(xu-xuu)  */
 interp5_xp[IDX(i,j,k)] = fabs( xp.x/dxp.x );
 interp5_xm[IDX(i,j,k)] = fabs( xm.x/dxm.x );
 interp5_yp[IDX(i,j,k)] = fabs( yp.y/dyp.y );
 interp5_ym[IDX(i,j,k)] = fabs( ym.y/dym.y );
 interp5_zp[IDX(i,j,k)] = fabs( zp.z/dzp.z );
 interp5_zm[IDX(i,j,k)] = fabs( zm.z/dzm.z );

 /*******************************************************************/
/*  Part 2: from outside to the inside */ 

/* we find the distance of the centers from such 6 surface centers */
/* x - xuu */
 xp = vector_difference(w_xp, center_V[IDX(i+1, j, k)]);
 xm = vector_difference(w_xm, center_V[IDX(i-1, j, k)]);
 yp = vector_difference(w_yp, center_V[IDX(i, j+1, k)]);
 ym = vector_difference(w_ym, center_V[IDX(i, j-1, k)]);
 zp = vector_difference(w_zp, center_V[IDX(i, j, k+1)]);
 zm = vector_difference(w_zm, center_V[IDX(i, j, k-1)]);
 /* we find the distance between the two neighboring nodes */ 
 /* xu - xuu */
 dxp = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i+2, j, k)]);
 dxm = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i-2, j, k)]);
 dyp = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j+2, k)]);
 dym = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j-2, k)]);
 dzp = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j, k+2)]);
 dzm = vector_difference(center_V[IDX(i, j, k)],center_V[IDX(i, j, k-2)]);

 /* This is DOWNWIND g3 = (x-xuu)/(xu-xuu)    */
 interp6_xp[IDX(i,j,k)] =  fabs(xp.x/dxp.x);
 interp6_xm[IDX(i,j,k)] =  fabs(xm.x/dxm.x);
 interp6_yp[IDX(i,j,k)] =  fabs(yp.y/dyp.y);
 interp6_ym[IDX(i,j,k)] =  fabs(ym.y/dym.y);
 interp6_zp[IDX(i,j,k)] =  fabs(zp.z/dzp.z);
 interp6_zm[IDX(i,j,k)] =  fabs(zm.z/dzm.z);
#endif

			}/*i,j,k*/

	/*  send receive the coefficients */	
	sendrecv_borders_scalar(interp5_xp);
	sendrecv_borders_scalar(interp5_xm);
	sendrecv_borders_scalar(interp5_yp);
	sendrecv_borders_scalar(interp5_ym);
	sendrecv_borders_scalar(interp5_zp);
	sendrecv_borders_scalar(interp5_zm);

	sendrecv_borders_scalar(interp6_xp);
	sendrecv_borders_scalar(interp6_xm);
	sendrecv_borders_scalar(interp6_yp);
	sendrecv_borders_scalar(interp6_ym);
	sendrecv_borders_scalar(interp6_zp);
	sendrecv_borders_scalar(interp6_zm);  
#endif


}/* end func interpolation coefficients */


/**************************************************************************************************/
#ifdef LB_FLUID_BC
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



#ifdef LB_FLUID_FORCING_LANDSCAPE
/* this makes a flag for a cylinder with axis along y */
my_double cylinder(int i , int j , int k, vector center, my_double radius, my_double height){
  my_double dum;

  if(center_V[IDX(i,j,k)].z <= height && sqrt(pow(center_V[IDX(i,j,k)].x-center.x, 2.0)+pow(center_V[IDX(i,j,k)].y-center.y, 2.0)) <= radius )
    dum = 1.0;
  else
    dum =0.0;

    return dum; 
}

/* this makes a flag for a cubick block */
my_double cubic_block(int i , int j , int k, vector center, vector size){
  my_double dum;

  if(fabs(center_V[IDX(i,j,k)].x-center.x) <= size.x/2.0 && fabs(center_V[IDX(i,j,k)].y) <= size.y && fabs(center_V[IDX(i,j,k)].z-center.z) <= size.z/2.0 )
    dum = 1.0;
  else
    dum =0.0;

    return dum; 
}

void read_landscape(){


	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;
	int             i, j, k;
	my_double  fac;
	int count;
	vector center;
	my_double radius;
	my_double height;
	vector size;


	sprintf(fnamein, "landscape.in");
	fin = fopen(fnamein, "r");
	if (fin != NULL) {
	  if(ROOT) fprintf(stderr, "Landscape file %s has been found!\n", fnamein);
	} else {

	  if(ROOT) fprintf(stderr, "Warning message -> %s file is missing!\n Starting from grid generated on the fly\n ", fnamein);


	  /*
		for (k =0; k < LNZ+TWO_BRD; k++)
			for (j =0; j < LNY+TWO_BRD; j++)
				for (i = 0; i < LNX+TWO_BRD; i++) {
				  if( sqrt(pow(center_V[IDX(i,j,k)].x-property.SX/2.0, 2.0)+pow(center_V[IDX(i,j,k)].y-property.SY/2.0, 2.0)) < 1)
				  landscape[IDX(i, j, k)] = 1.0;
				}
	  */

	  // defining different types of bluff bodies : cyinder, cube, buildings 
#ifdef  LB_FLUID_FORCING_LANDSCAPE_CYLINDER
			    for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){
				  center.x = (property.SX/5.0)+1.0;     //property.SX/2.0;
				  center.y = (property.SY/2.0)+3.0;     //property.SY/2.0;  
				  center.z = property.SZ/2.0;
				  radius = (property.SY/10.0)+1.0;   //200.0;
				  height = 300.0;
				  landscape[IDX(i, j, k)] += cylinder(i,j,k, center, radius,  height);
				}
#endif

#ifdef LB_FLUID_FORCING_LANDSCAPE_CUBE
			    for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){
				  center.x = property.SX/2.0 ;
    				  center.z = property.SZ/2.0 ;
				  size.x = 50.0;        //80.0;
				  size.y = 50.0;           //200.0;
				  size.z =  80.0;
				  landscape[IDX(i, j, k)] += cubic_block(i,j,k, center, size);
				}
#endif

#ifdef LB_FLUID_FORCING_LANDSCAPE_BUILDINGS
			    for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){
				  center.x = property.SX/2.58 ;
    				  center.z = property.SZ/2.0 ;
				  size.x = 10.0;        //80.0;
				  size.y = 20.0;           //200.0;
				  size.z =  10.0;
				  landscape[IDX(i, j, k)] += cubic_block(i,j,k, center, size);
				}

			    for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){
				  center.x = property.SX/2.05 ;
    				  center.z = property.SZ/2.0 ;
				  size.x = 10.0;        //80.0;
				  size.y = 50.0;           //200.0;
				  size.z =  10.0;
				  landscape[IDX(i, j, k)] += cubic_block(i,j,k, center, size);
				}

			    for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){
				  center.x = property.SX/1.57 ;
    				  center.z = property.SZ/2.0 ;
				  size.x = 10.0;        //80.0;
				  size.y = 30.0;           //200.0;
				  size.z =  10.0;
				  landscape[IDX(i, j, k)] += cubic_block(i,j,k, center, size);
				}

			    for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){
				  center.x = property.SX/1.27 ;
    				  center.z = property.SZ/2.0 ;
				  size.x = 10.0;        //80.0;
				  size.y = 60.0;           //200.0;
				  size.z =  10.0;
				  landscape[IDX(i, j, k)] += cubic_block(i,j,k, center, size);
				}

#endif

#ifdef LB_FLUID_FORCING_LANDSCAPE_BUILDINGS2
                            for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){
				  center.x = property.SX/2.034 ;
    				  center.z = property.SZ/2.0 ;
				  size.x = 10.0;        //80.0;
				  size.y = 50.0;           //200.0;
				  size.z =  10.0;
				  landscape[IDX(i, j, k)] += cubic_block(i,j,k, center, size);
				}

			    for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){
				  center.x = property.SX/1.558 ;
    				  center.z = property.SZ/2.0 ;
				  size.x = 10.0;        //80.0;
				  size.y = 70.0;           //200.0;
				  size.z =  10.0;
				  landscape[IDX(i, j, k)] += cubic_block(i,j,k, center, size);
				}

#endif

			    /*

	               //  for setting the walls to zero 

			    for (i = 0; i < LNX+TWO_BRD; i++) 
			      for (k =0; k < LNZ+TWO_BRD; k++)
			        for (j =0; j < LNY+TWO_BRD; j++){

				  if((LNY_END == NY && j==LNY+BRD-1 ) || (LNY_START == 0 && j==BRD )){
				  //if((LNY_END == NY && j>=LNY+BRD-2  ) || (LNY_START == 0 &&  j<=BRD+1 )){
				      landscape[IDX(i, j, k)] = 1;
				      }

                             }
*/
}
}
#endif
