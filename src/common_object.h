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
extern int me;
extern int mex,mey,mez;
extern int me_xp, me_xm, me_yp,me_ym, me_zp,me_zm;
extern int nprocs;
extern int nxprocs, nyprocs, nzprocs;
extern MPI_Datatype MPI_property_type , MPI_pop_type , MPI_vector_type, MPI_output_type;
extern MPI_Op MPI_SUM_output;

/* random seed */
extern unsigned int seed;

/* System size , center nodes */
extern int NX , NY , NZ;
extern int LNX , LNY , LNZ;
extern int LNX_END , LNY_END , LNZ_END;
extern int LNX_START , LNY_START , LNZ_START;
/* for the grid */
extern int NXG , NYG , NZG;
extern int LNXG , LNYG , LNZG;
extern int LNXG_END , LNYG_END , LNZG_END;
extern int LNXG_START , LNYG_START , LNZG_START;


/* Array of properties*/
extern prop property;

/* Array of outputs */
extern output out_local, out_all;
extern output *ruler_x_local, *ruler_y_local, *ruler_z_local;
extern output *ruler_x, *ruler_y, *ruler_z;

/* mesh */
extern vector *mesh, *center_V;
extern int *mesh_flag;
extern vector *xp_mesh,*xm_mesh,*yp_mesh,*ym_mesh,*zp_mesh,*zm_mesh;
extern int *xp_flag,*xm_flag,*yp_flag,*ym_flag,*zp_flag,*zm_flag;
extern pop *coeff_xp, *coeff_xm, *coeff_yp, *coeff_ym, *coeff_zp, *coeff_zm;

extern pop *xp_pop,*xm_pop,*yp_pop,*ym_pop,*zp_pop,*zm_pop;
#ifdef LB_BC
extern pop *norm_xp_pop,*norm_xm_pop,*norm_yp_pop,*norm_ym_pop,*norm_zp_pop,*norm_zm_pop;
#endif

/* Populations */
extern int NPOP;

/* LB speeds & weights */
extern my_double wgt[19]; 
extern vector c[19];
extern my_double dirp[19] ,inv[19];
extern my_double twocs2 , twocs4;
extern my_double invtwocs2 , invtwocs4;
extern my_double cs, cs2 , cs4 , cs22 , cssq;
extern my_double invcs, invcs2, invcs4;

#ifdef LB_FLUID
extern pop *p, *rhs_p, *old_rhs_p;
extern my_double *dens;
extern vector *u;
#ifdef LB_FLUID_FORCING
extern vector *force;
#endif
#endif

/* time */
extern int itime;
extern my_double time_now;
extern char OutDir[256];























