#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <mpi.h>
#include "define.h"
#include "typedef.h"
#include "functions.h"


/* MPI variables */
int me;
int mex,mey,mez;
int me_xp, me_xm, me_yp,me_ym, me_zp,me_zm;
int nprocs;
int nxprocs, nyprocs, nzprocs;
MPI_Datatype MPI_property_type , MPI_pop_type , MPI_vector_type, MPI_output_type;
MPI_Op MPI_SUM_output;

/* random seed */
unsigned int seed;

/* System size , for the nodes*/
int NX , NY , NZ;
int LNX , LNY , LNZ;
int LNX_END , LNY_END , LNZ_END;
int LNX_START , LNY_START , LNZ_START;
/* for the grid */
int NXG , NYG , NZG;
int LNXG , LNYG , LNZG;
int LNXG_END , LNYG_END , LNZG_END;
int LNXG_START , LNYG_START , LNZG_START;


/* Array of properties*/
prop property;

/* Array of outputs */
output out_local,out_all;
output *ruler_x_local, *ruler_y_local, *ruler_z_local;
output *ruler_x, *ruler_y, *ruler_z;

/* Mesh matrix */
vector *mesh, *center_V;
int *mesh_flag;
vector *xp_mesh,*xm_mesh,*yp_mesh,*ym_mesh,*zp_mesh,*zm_mesh;
int *xp_flag,*xm_flag,*yp_flag,*ym_flag,*zp_flag,*zm_flag;
pop *coeff_xp, *coeff_xm, *coeff_yp, *coeff_ym, *coeff_zp, *coeff_zm;

pop *xp_pop,*xm_pop,*yp_pop,*ym_pop,*zp_pop,*zm_pop;
#ifdef LB_BC
pop *norm_xp_pop,*norm_xm_pop,*norm_yp_pop,*norm_ym_pop,*norm_zp_pop,*norm_zm_pop;
#endif

/* Populations */
int NPOP;

/* LB speeds & weights */
my_double wgt[19]; 
vector c[19];
int dirp[19] ,inv[19];
my_double twocs2 , twocs4;
my_double invtwocs2 , invtwocs4;
my_double cs, cs2 , cs4 , cs22 , cssq;
my_double invcs, invcs2, invcs4;


#ifdef LB_FLUID
pop *p, *rhs_p, *old_rhs_p;
my_double *dens;
vector *u;
#ifdef LB_FLUID_FORCING
vector *force;
#endif
#endif



/* time */
int itime;
my_double time_now;
char OutDir[256];























