#include "common_object.h"


//void sum_output(output *a, output *b, int *length, MPI_Datatype *dtype);
//void sum_vector(vector *a, vector *b, int *length, MPI_Datatype *dtype);

void sum_my_double(my_double *a, my_double *b, int *length, MPI_Datatype *dtype){
  int i;

  for (i = 0; i < *length; i++) {

    b[i]  += a[i];

  }
}

 
void sum_output(output *a, output *b, int *length, MPI_Datatype *dtype){
  int i;

for (i = 0; i < *length; i++) {

#ifdef LB_FLUID
  (b+i)->x   += (a+i)->x;
  (b+i)->y   += (a+i)->y;
  (b+i)->z   += (a+i)->z;
  (b+i)->ux  += (a+i)->ux;
  (b+i)->uy  += (a+i)->uy;
  (b+i)->uz  += (a+i)->uz;
  (b+i)->ux2 += (a+i)->ux2;
  (b+i)->uy2 += (a+i)->uy2;
  (b+i)->uz2 += (a+i)->uz2;
  (b+i)->rho += (a+i)->rho;
  (b+i)->ene += (a+i)->ene;
  (b+i)->eps += (a+i)->eps;
  (b+i)->rho2 += (a+i)->rho2;
#endif
#ifdef LB_TEMPERATURE
  (b+i)->dxt  += (a+i)->dxt;
  (b+i)->dyt  += (a+i)->dyt;
  (b+i)->dzt  += (a+i)->dzt;
  (b+i)->uxt  += (a+i)->uxt;
  (b+i)->uyt  += (a+i)->uyt;
  (b+i)->uzt  += (a+i)->uzt;
  (b+i)->nux  += (a+i)->nux;
  (b+i)->nuy  += (a+i)->nuy;
  (b+i)->nuz  += (a+i)->nuz;
  (b+i)->t    += (a+i)->t;
  (b+i)->t2   += (a+i)->t2;
  (b+i)->epst += (a+i)->epst;
  //  (b+i)->lb += (a+i)->lb;
 #ifdef LB_TEMPERATURE_MELTING
  (b+i)->lf   += (a+i)->lf;
  (b+i)->dtlf += (a+i)->dtlf;
  (b+i)->enth += (a+i)->enth;
  (b+i)->lf2   += (a+i)->lf2;
 #endif
#endif
#ifdef LB_SCALAR
  (b+i)->dxs  += (a+i)->dxs;
  (b+i)->dys  += (a+i)->dys;
  (b+i)->dzs  += (a+i)->dzs;
  (b+i)->uxs  += (a+i)->uxs;
  (b+i)->uys  += (a+i)->uys;
  (b+i)->uzs  += (a+i)->uzs;
  (b+i)->nusx  += (a+i)->nusx;
  (b+i)->nusy  += (a+i)->nusy;
  (b+i)->nusz  += (a+i)->nusz;
  (b+i)->s    += (a+i)->s;
  (b+i)->s2   += (a+i)->s2;
  (b+i)->epss += (a+i)->epss;
#endif
 }

}


void sum_vector(vector *a, vector *b, int *length, MPI_Datatype *dtype){
  int i;

  for (i = 0; i < *length; i++) {

    (b+i)->x   += (a+i)->x;
    (b+i)->y   += (a+i)->y;
    (b+i)->z   += (a+i)->z;
  }
}

#ifdef LAGRANGE
void sum_output_particle(output_particle *a, output_particle *b,  int *length, MPI_Datatype *dtype){
  int i;
  for(i=0;i<*length;i++){
    (b+i)->vx += (a+i)->vx;
    (b+i)->vy += (a+i)->vy;
    (b+i)->vz += (a+i)->vz;

    (b+i)->vx2 += (a+i)->vx2;
    (b+i)->vy2 += (a+i)->vy2;
    (b+i)->vz2 += (a+i)->vz2;

    (b+i)->vx4 += (a+i)->vx4;
    (b+i)->vy4 += (a+i)->vy4;
    (b+i)->vz4 += (a+i)->vz4;

    (b+i)->ax += (a+i)->ax;
    (b+i)->ay += (a+i)->ay;
    (b+i)->az += (a+i)->az;

    (b+i)->ax2 += (a+i)->ax2;
    (b+i)->ay2 += (a+i)->ay2;
    (b+i)->az2 += (a+i)->az2;

    (b+i)->ax4 += (a+i)->ax4;
    (b+i)->ay4 += (a+i)->ay4;
    (b+i)->az4 += (a+i)->az4;

 #ifdef LAGRANGE_GRADIENT
  #ifdef LAGRANGE_ORIENTATION
    (b+i)->dt_px += (a+i)->dt_px;
    (b+i)->dt_py += (a+i)->dt_py;
    (b+i)->dt_pz += (a+i)->dt_pz;

    (b+i)->dt_px2 += (a+i)->dt_px2;
    (b+i)->dt_py2 += (a+i)->dt_py2;
    (b+i)->dt_pz2 += (a+i)->dt_pz2;

    (b+i)->dt_px4 += (a+i)->dt_px4;
    (b+i)->dt_py4 += (a+i)->dt_py4;
    (b+i)->dt_pz4 += (a+i)->dt_pz4;
   #ifdef LAGRANGE_ORIENTATION_SECONDORIENTATION
    (b+i)->dt_nx += (a+i)->dt_nx;
    (b+i)->dt_ny += (a+i)->dt_ny;
    (b+i)->dt_nz += (a+i)->dt_nz;

    (b+i)->dt_nx2 += (a+i)->dt_nx2;
    (b+i)->dt_ny2 += (a+i)->dt_ny2;
    (b+i)->dt_nz2 += (a+i)->dt_nz2;

    (b+i)->dt_nx4 += (a+i)->dt_nx4;
    (b+i)->dt_ny4 += (a+i)->dt_ny4;
    (b+i)->dt_nz4 += (a+i)->dt_nz4;
   #endif
  #endif
 #endif
#ifdef LB_TEMPERATURE
    (b+i)->t  +=  (a+i)->t;
    (b+i)->t2 += (a+i)->t2;
    (b+i)->t4 += (a+i)->t4;
    (b+i)->dt_t += (a+i)->dt_t;
    (b+i)->dt_t2 += (a+i)->dt_t2;
	#ifdef LAGRANGE_TEMPERATURE 
		(b+i)->t_p += (a+i)->t_p;
	#endif	
#endif
 #ifdef LAGRANGE_OUTPUT_FLUID_AVERAGES
    (b+i)->ux += (a+i)->ux;
    (b+i)->uy += (a+i)->uy;
    (b+i)->uz += (a+i)->uz;

    (b+i)->ux2 += (a+i)->ux2;
    (b+i)->uy2 += (a+i)->uy2;
    (b+i)->uz2 += (a+i)->uz2;

    (b+i)->ux4 += (a+i)->ux4;
    (b+i)->uy4 += (a+i)->uy4;
    (b+i)->uz4 += (a+i)->uz4;
 #endif
  }
}
#endif

/* The following is the first function called in main */
void initialization_MPI(int argc, char **argv){
	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (ROOT)
		fprintf(stderr, "Hello world! I am ROOT: processor %d of %d\n", me, nprocs);


/* define and commit types */
 MPI_Type_contiguous(1, MPI_DOUBLE , &MPI_MY_DOUBLE); /* this is the default : double precision */
#ifdef PRECISION_FLOAT
 MPI_Type_contiguous(1, MPI_FLOAT , &MPI_MY_DOUBLE);
#endif
#ifdef PRECISION_DOUBLE_LONG
 MPI_Type_contiguous(1, MPI_LONG_DOUBLE , &MPI_MY_DOUBLE);
#endif

 MPI_Type_commit(&MPI_MY_DOUBLE);
 MPI_Op_create( (MPI_User_function *)sum_my_double, 1, &MPI_SUM_my_double );

#ifdef DEBUG
	/* check  of the sum between my_double operation */
	my_double prova, prova_all;
	prova = 1;
	prova_all = 2;
	MPI_Allreduce(&prova, &prova_all, 1, MPI_MY_DOUBLE, MPI_SUM_my_double, MPI_COMM_WORLD );
	if (ROOT) fprintf(stderr, "MPI_SUM_my_double check %e x %d = %e\n",prova,nprocs,prova_all);
/* end of check */
#endif 

 //MPI_Type_contiguous(1, MPI_MY_DOUBLE, &MPI_my_double_type); /* duplicated but still used in the code e.g. for SENDRECV */
 //MPI_Type_commit(&MPI_my_double_type);

 MPI_Type_contiguous(NPROP, MPI_MY_DOUBLE , &MPI_property_type);
 MPI_Type_commit(&MPI_property_type);

 MPI_Type_contiguous(NOUT, MPI_MY_DOUBLE , &MPI_output_type);
 MPI_Type_commit(&MPI_output_type);
 MPI_Op_create( (MPI_User_function *)sum_output, 1, &MPI_SUM_output );
 if(ROOT) fprintf(stderr,"--------> NOUT size %d\n",(int)NOUT);

 MPI_Type_contiguous(3, MPI_MY_DOUBLE, &MPI_vector_type);
 MPI_Type_commit(&MPI_vector_type);
 MPI_Op_create( (MPI_User_function *)sum_vector, 1, &MPI_SUM_vector );

#ifdef LB
 //NPOP = 19;
 /* commit pop type */
 MPI_Type_contiguous(NPOP, MPI_MY_DOUBLE, &MPI_pop_type);
 MPI_Type_commit(&MPI_pop_type);
#endif

#ifdef LAGRANGE 
 MPI_Type_contiguous(SIZE_OF_POINT_PARTICLE, MPI_MY_DOUBLE, &MPI_point_particle_type);
 MPI_Type_commit(&MPI_point_particle_type);

 MPI_Type_contiguous(SIZE_OF_OUTPUT_PARTICLE, MPI_MY_DOUBLE , &MPI_output_particle_type);
 MPI_Type_commit(&MPI_output_particle_type);
 MPI_Op_create( (MPI_User_function *)sum_output_particle, 1, &MPI_SUM_output_particle);
 if(ROOT) fprintf(stderr,"--------> SIZE_OF_OUTPUT_PARTICLE size %d\n",(int)SIZE_OF_OUTPUT_PARTICLE);
#endif

 /* Initialize random seeds */
#ifdef SYSTEM_RANDOM48
  seed = time(NULL);
  srand48(me+seed);
#else
  initdum = -time(NULL); /* the minus sign is very important */
  initdum -= me;
  idum = &initdum;
#endif

#ifdef VERBOSE
fprintf(stderr,"RANDOM me %d myrand(): %e\n",me, myrand());
#endif


}


/**********************************************************************************/

int compare(const void *a, const void *b)
{
	return (*(int *) a - *(int *) b);
}


int check_prime(int x){
  int i;
  int isprime = 0;
    for(i=2;i<x;i++)
        if(x%i == 0) {
            isprime++;
            break;
        }
    return isprime;
}


void processor_splitting()
{

	int             tnprocs, nleft, ncubic;
	int             sort3d[3], sort2d[2];
	int             i, j, k, ami, next, next_x, next_y, next_z;
	int error=0;
	int             prime;

	/* grid processor splitting */
	LNX = NX;
	LNY = NY;
	LNZ = NZ;
	tnprocs = nprocs;
	nxprocs = nyprocs = nzprocs = 1;
	prime = 2;

	if (ROOT) {

	  /* Try divisibility by prime factors */
		while (nxprocs * nyprocs * nzprocs < nprocs) {
			sort3d[0] = (LNX % prime == 0) ? LNX : 0;
			sort3d[1] = (LNY % prime == 0) ? LNY : 0;
			sort3d[2] = (LNZ % prime == 0) ? LNZ : 0;
			qsort(sort3d, 3, sizeof(int), compare);
			if (sort3d[2] == LNX) {
				LNX /= prime;
				nxprocs *= prime;
			} else {
				if (sort3d[2] == LNY) {
					LNY /= prime;
					nyprocs *= prime;
				} else {
					if (sort3d[2] == LNZ) {
						LNZ /= prime;
						nzprocs *= prime;
					}
				}
			}
			//fprintf(stderr,"LNX %d , LNY %d , LNZ %d\n",LNX,LNY,LNZ);
			//fprintf(stderr,"nxprocs %d , nyprocs %d , nzprocs %d\n",nxprocs,nyprocs,nzprocs);
			if (LNX % prime != 0 && LNY % prime != 0 && LNZ % prime != 0){
			  /************ take the next prime factor *****************/
			  for(i=prime+1;;i++){
			    if(check_prime(i)){
			      prime=i;
			      break;
			    }
			  }
			  /*****************************/
			}
			
		}/* end of while loop */

		if (nxprocs * nyprocs * nzprocs == nprocs) {
			fprintf(stderr, "good splitting!\n");
			fprintf(stderr, "LNX %d LNY %d LNZ %d\n", LNX, LNY, LNZ);
			fprintf(stderr, "nxprocs %d nyprocs %d nzprocs %d total %d\n", nxprocs, nyprocs, nzprocs, nxprocs * nyprocs * nzprocs);
		} else {
				fprintf(stderr, "bad splitting: expected %d procs but available %d\nAborting.\n", nxprocs * nyprocs * nzprocs, nprocs);
				error=1;
		}

	} /* if ROOT */

	/* check for error */
	MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(error){
	  MPI_Finalize();
	  exit(-1);
	}

	/* Now broadcast */
	MPI_Bcast(&LNX, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&LNY, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&LNZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nxprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nyprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nzprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef DEBUG
	for (i = 0; i < nprocs; i++) {
		if (i == me)
			fprintf(stderr, "me %d , System Size: LNX %d LNY %d LNZ %d\n   nxprocs %d nyprocs %d nzprocs %d\n", me, LNX, LNY, LNZ, nxprocs, nyprocs, nzprocs);
	}
#endif

	/* here we find the coordinate of the processors */
	for (k = 0; k < nzprocs; k++)
		for (j = 0; j < nyprocs; j++)
			for (i = 0; i < nxprocs; i++) {
				ami = k * (nyprocs * nxprocs) + j * nxprocs + i;
				if (ami == me) {
					mex = i;
					mey = j;
					mez = k;
				}
			}

#ifdef DEBUG
	fprintf(stderr, "me %d , mex %d  mey %d mez %d\n", me, mex, mey, mez);
#endif

	 	
	/* every procs finds its neighbors */
	/*
	next = (mex+1+nxprocs) % nxprocs; me_xp = mez * (nyprocs * nxprocs) + mey * nxprocs + next;
	next = (mex-1+nxprocs) % nxprocs; me_xm = mez * (nyprocs * nxprocs) + mey * nxprocs + next;
	next = (mey+1+nyprocs) % nyprocs; me_yp = mez * (nyprocs * nxprocs) + next * nxprocs + mex; 
	next = (mey-1+nyprocs) % nyprocs; me_ym = mez * (nyprocs * nxprocs) + next * nxprocs + mex;
	next = (mez+1+nzprocs) % nzprocs; me_zp = next * (nyprocs * nxprocs) + mey * nxprocs + mex; 
	next = (mez-1+nzprocs) % nzprocs; me_zm = next * (nyprocs * nxprocs) + mey * nxprocs + mex;

#ifdef DEBUG
  fprintf(stderr, "me  %d old way: me_xp %d, me_xm %d, me_yp %d, me_ym %d, me_zp %d, me_zm %d \n",me, me_xp, me_xm, me_yp, me_ym, me_zp, me_zm);
#endif
	*/

	/* Map of of next neighbor processors  */ 	
        me_next  = (int*) malloc(sizeof(int)*27);
	for (k = -1; k <=1; k++)
		for (j = -1; j <=1; j++)
			for (i = -1; i <=1; i++) {

			  next_x = (mex+i+nxprocs) % nxprocs;
			  next_y = (mey+j+nyprocs) % nyprocs;
			  next_z = (mez+k+nzprocs) % nzprocs;
			  me_next[IDX_NEXT(i,j,k)] = next_z * (nyprocs * nxprocs) + next_y * nxprocs + next_x;

			  // fprintf(stderr,"me %d, i=%d j=%d k=%d, me_next=%d\n", me,i,j,k,me_next[IDX_NEXT(i,j,k)]);
			}
	
 /* for the faces */
 me_xp = me_next[IDX_NEXT(  1,  0,  0)];
 me_xm = me_next[IDX_NEXT( -1,  0,  0)];
 me_yp = me_next[IDX_NEXT(  0,  1,  0)];
 me_ym = me_next[IDX_NEXT(  0, -1,  0)];
 me_zp = me_next[IDX_NEXT(  0,  0,  1)]; 
 me_zm = me_next[IDX_NEXT(  0,  0, -1)]; 	

#ifdef DEBUG
 fprintf(stderr, "me %d new way: me_xp %d, me_xm %d, me_yp %d, me_ym %d, me_zp %d, me_zm %d \n",me, me_xp, me_xm, me_yp, me_ym, me_zp, me_zm);
#endif
 
#ifdef METHOD_EDGES_AND_CORNERS

 /* for the corners */
 me_xp_yp_zp = me_next[IDX_NEXT(  1,  1,  1)];
 me_xm_ym_zm = me_next[IDX_NEXT( -1, -1, -1)];
 me_xp_yp_zm = me_next[IDX_NEXT(  1,  1, -1)];
 me_xm_ym_zp = me_next[IDX_NEXT( -1, -1,  1)];
 me_xp_ym_zp = me_next[IDX_NEXT(  1, -1,  1)]; 
 me_xm_yp_zm = me_next[IDX_NEXT( -1,  1, -1)]; 
 me_xm_yp_zp = me_next[IDX_NEXT( -1,  1,  1)]; 
 me_xp_ym_zm = me_next[IDX_NEXT(  1, -1, -1)]; 
	
 /* for the edges */
 me_xp_yp = me_next[IDX_NEXT(  1,  1,  0)]; 
 me_xm_ym = me_next[IDX_NEXT( -1, -1,  0)]; 
 me_xp_ym = me_next[IDX_NEXT(  1, -1,  0)]; 
 me_xm_yp = me_next[IDX_NEXT( -1,  1,  0)]; 
 me_yp_zp = me_next[IDX_NEXT(  0,  1,  1)]; 
 me_ym_zm = me_next[IDX_NEXT(  0, -1, -1)]; 
 me_yp_zm = me_next[IDX_NEXT(  0,  1, -1)]; 
 me_ym_zp = me_next[IDX_NEXT(  0, -1,  1)]; 
 me_xp_zp = me_next[IDX_NEXT(  1,  0,  1)]; 
 me_xm_zm = me_next[IDX_NEXT( -1,  0, -1)]; 
 me_xp_zm = me_next[IDX_NEXT(  1,  0, -1)]; 
 me_xm_zp = me_next[IDX_NEXT( -1,  0,  1)]; 

#endif


#ifdef PARALLEL_MPICART
 /* Alternative way to split the computational domain on a grid of processor, it is LESS performant than the previous mode */
 MPI_Comm Comm_3d_grid;
 int ndim = 3;
 int dimensions[3];
 int periodic[3];
 int reorder = 1;
 int coord[3];
 int me_rank;


 /* Give partial freedom to MPI to chose thes sizes of the grid of processors */
 if(NX==1){dimensions[0] = 1;}else{dimensions[0] = 0;}
 if(NY==1){dimensions[1] = 1;}else{dimensions[1] = 0;}
 if(NZ==1){dimensions[2] = 1;}else{dimensions[2] = 0;}
 fprintf(stderr,"Restrictions (%d, %d, %d) give decomposition ", dimensions[0], dimensions[1], dimensions[2]);
 MPI_Dims_create(nprocs, 3, dimensions);
 fprintf(stderr,"(%d, %d, %d).\n", dimensions[0], dimensions[1], dimensions[2]);
 
 nxprocs = dimensions[0]; 
 nyprocs = dimensions[1];
 nzprocs = dimensions[2];
	
 error = 0;
 if(NX%nxprocs==0){ LNX=NX/nxprocs; }else{error=1;}
 if(NY%nyprocs==0){ LNY=NY/nyprocs; }else{error=1;}
 if(NZ%nzprocs==0){ LNZ=NZ/nzprocs; }else{error=1;}
 if(error){
    MPI_Finalize();
    exit(-1);
  }
	
 periodic[0] = periodic[1] = periodic[2] =  1;
 MPI_Cart_create(MPI_COMM_WORLD,ndim,dimensions,periodic,reorder,&Comm_3d_grid);
 MPI_Comm_rank(Comm_3d_grid,&me_rank);
 MPI_Cart_coords(Comm_3d_grid,me_rank,ndim,coord);

 fprintf(stderr,"I am %d: (%d,%d,%d); originally %d\n",me_rank,coord[0],coord[1],coord[2], me);

 mex = coord[0];
 mey = coord[1];
 mez = coord[2];
 me = me_rank;       

 /* 6 faces */
 coord[0]=mex+1; coord[1]=mey  ; coord[2]=mez  ; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp);
 coord[0]=mex-1; coord[1]=mey  ; coord[2]=mez  ; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm);
 coord[0]=mex  ; coord[1]=mey+1; coord[2]=mez  ; MPI_Cart_rank(Comm_3d_grid, coord, &me_yp);
 coord[0]=mex  ; coord[1]=mey-1; coord[2]=mez  ; MPI_Cart_rank(Comm_3d_grid, coord, &me_ym);
 coord[0]=mex  ; coord[1]=mey  ; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_zp);
 coord[0]=mex  ; coord[1]=mey  ; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_zm);
 #ifdef METHOD_EDGES_AND_CORNERS
  /* 8 corners */
 coord[0]=mex+1; coord[1]=mey+1; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp_yp_zp);
 coord[0]=mex-1; coord[1]=mey-1; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm_ym_zm);
 coord[0]=mex+1; coord[1]=mey+1; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp_yp_zm);
 coord[0]=mex-1; coord[1]=mey-1; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm_ym_zp);
 coord[0]=mex+1; coord[1]=mey-1; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp_ym_zp);
 coord[0]=mex-1; coord[1]=mey+1; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm_yp_zm);
 coord[0]=mex-1; coord[1]=mey+1; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm_yp_zp);
 coord[0]=mex+1; coord[1]=mey-1; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp_ym_zm);
 /* 12 edges */
 coord[0]=mex+1; coord[1]=mey+1; coord[2]=mez  ; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp_yp);
 coord[0]=mex-1; coord[1]=mey-1; coord[2]=mez  ; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm_ym);
 coord[0]=mex+1; coord[1]=mey-1; coord[2]=mez  ; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp_ym);
 coord[0]=mex-1; coord[1]=mey+1; coord[2]=mez  ; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm_yp);
 coord[0]=mex  ; coord[1]=mey+1; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_yp_zp);
 coord[0]=mex  ; coord[1]=mey-1; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_ym_zm);
 coord[0]=mex  ; coord[1]=mey+1; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_yp_zm);
 coord[0]=mex  ; coord[1]=mey-1; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_ym_zp);
 coord[0]=mex+1; coord[1]=mey  ; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp_zp);
 coord[0]=mex-1; coord[1]=mey  ; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm_zm);
 coord[0]=mex+1; coord[1]=mey  ; coord[2]=mez-1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xp_zm);
 coord[0]=mex-1; coord[1]=mey  ; coord[2]=mez+1; MPI_Cart_rank(Comm_3d_grid, coord, &me_xm_zp);  
 #endif
#endif /* end of PARALLEL_MPICART */
 

#ifdef DEBUG
	fprintf(stderr, "me %d , me_xp %d  me_xm %d me_yp %d me_ym %d me_zp %d me_zm %d\n", me, me_xp, me_xm, me_yp, me_ym, me_zp, me_zm);
#endif




	/* processor rulers for vertices (computational node points) */
	LNX_START = LNX * mex;
	LNX_END = LNX * (mex + 1);
	LNY_START = LNY * mey;
	LNY_END = LNY * (mey + 1);
	LNZ_START = LNZ * mez;
	LNZ_END = LNZ * (mez + 1);

#ifdef DEBUG
	fprintf(stderr, "me %d LNX_START %d , LNY_START %d  LNZ_START %d\n", me, LNX_START ,LNY_START, LNZ_START);
	fprintf(stderr, "me %d LNX_END %d , LNY_END %d  LNZ_END %d\n", me, LNX_END ,LNY_END, LNZ_END);
#endif

	/* processor rulers for the grid */
	NXG=NX+1;
	NYG=NY+1;
	NZG=NZ+1;
	LNXG=LNX+1;
	LNYG=LNY+1;
	LNZG=LNZ+1;
	/*
	LNXG_START = LNXG * mex;
	LNXG_END = LNXG * (mex + 1);
	LNYG_START = LNYG * mey;
	LNYG_END = LNYG * (mey + 1);
	LNZG_START = LNZG * mez;
	LNZG_END = LNZG * (mez + 1);
	*/	
	LNXG_START = LNX * mex;
	LNXG_END = LNXG * (mex + 1);
	LNYG_START = LNY * mey;
	LNYG_END = LNYG * (mey + 1);
	LNZG_START = LNZ * mez;
	LNZG_END = LNZG * (mez + 1);
	

#ifdef DEBUG
	fprintf(stderr, "me %d LNXG_START %d , LNYG_START %d  LNZG_START %d\n", me, LNXG_START ,LNYG_START, LNZG_START);
	fprintf(stderr, "me %d LNXG_END %d , LNYG_END %d  LNZG_END %d\n", me, LNXG_END ,LNYG_END, LNZG_END);
#endif

	

#ifdef NEW_SENDRECV
	/* Creates a vector (strided) datatype */
	/* Input: 
	   count -> number of blocks (nonnegative integer)
	   blocklength -> number of elements in each block (nonnegative integer)
	   stride -> number of elements between start of each block (integer)
	   oldtype -> old datatype (handle)
	   IDX(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )
	   IDX(i,j,k) (k*(LNY+TWO_BRD)+j)*(LNX+TWO_BRD)+i
	*/
	/* pop */
    MPI_Type_vector((LNZ+TWO_BRD)*(LNY+TWO_BRD), BRD , (LNX+TWO_BRD) , MPI_pop_type, &MPI_pop_plane_x);
    MPI_Type_vector((LNZ+TWO_BRD), (LNX+TWO_BRD)*BRD , (LNX+TWO_BRD)*(LNY+TWO_BRD) , MPI_pop_type, &MPI_pop_plane_y);
    MPI_Type_vector(BRD, (LNX+TWO_BRD)*(LNY+TWO_BRD), (LNX+TWO_BRD)*(LNY+TWO_BRD), MPI_pop_type, &MPI_pop_plane_z);

    MPI_Type_commit(&MPI_pop_plane_x);
    MPI_Type_commit(&MPI_pop_plane_y);
    MPI_Type_commit(&MPI_pop_plane_z);

      /* vector */ 
    MPI_Type_vector((LNZ+TWO_BRD)*(LNY+TWO_BRD), BRD , (LNX+TWO_BRD) , MPI_vector_type, &MPI_vector_plane_x);
    MPI_Type_vector((LNZ+TWO_BRD), (LNX+TWO_BRD)*BRD , (LNX+TWO_BRD)*(LNY+TWO_BRD) , MPI_vector_type, &MPI_vector_plane_y);
    MPI_Type_vector(BRD, (LNX+TWO_BRD)*(LNY+TWO_BRD), (LNX+TWO_BRD)*(LNY+TWO_BRD), MPI_vector_type, &MPI_vector_plane_z);
  
    MPI_Type_commit(&MPI_vector_plane_x);
    MPI_Type_commit(&MPI_vector_plane_y);
    MPI_Type_commit(&MPI_vector_plane_z);

     /* scalar */
    MPI_Type_vector((LNZ+TWO_BRD)*(LNY+TWO_BRD), BRD , (LNX+TWO_BRD) , MPI_MY_DOUBLE, &MPI_my_double_plane_x);
    MPI_Type_vector((LNZ+TWO_BRD), (LNX+TWO_BRD)*BRD , (LNX+TWO_BRD)*(LNY+TWO_BRD) , MPI_MY_DOUBLE, &MPI_my_double_plane_y);
    MPI_Type_vector(BRD, (LNX+TWO_BRD)*(LNY+TWO_BRD), (LNX+TWO_BRD)*(LNY+TWO_BRD), MPI_MY_DOUBLE, &MPI_my_double_plane_z);

    MPI_Type_commit(&MPI_my_double_plane_x);
    MPI_Type_commit(&MPI_my_double_plane_y);
    MPI_Type_commit(&MPI_my_double_plane_z);
#endif

}



void measure_time(){

  FILE* fout;
  /* t1,t2 and tick are global variables */
	tick = MPI_Wtick();
	if(ROOT){ 
	 remove("run_time.dat");
         fout = fopen("run_time.dat","w");
	 fprintf(fout,"Total execution time %e\n",t2-t1);
	 fprintf(stdout,"Total execution time %e\n",t2-t1);
	 fprintf(fout,"Time steps %d\n",itime);	 
	 fprintf(fout,"Execution time per time step %e\n",(t2-t1)/(double)itime);
	 fprintf(fout,"Execution time per grid point %e\n",(t2-t1)/(double)(NX*NY*NZ));
	 fprintf(fout,"Execution time per time step and grid point %e\n",(t2-t1)/(double)(itime*NX*NY*NZ));
	 fprintf(fout,"Time tick on this machine %e\n", tick);
	 fprintf(fout,"Number of processes for this run %d\n",nprocs);
         fprintf(fout,"Execution time per process %e\n",(t2-t1)/(double)nprocs);
	 fclose(fout);	 
	}

}
