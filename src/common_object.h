#include <stdio.h>
#include <math.h>
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
extern int nprocs;
extern int nxprocs, nyprocs, nzprocs;

/* System size */
extern int NX , NY , NZ;
extern int LNX , LNY , LNZ;
extern int LNX_END , LNY_END , LNZ_END;
extern int LNX_START , LNY_START , LNZ_START;

/* Array of properties*/
extern prop property;

/* mesh */
extern mesh_type *mesh;
extern surf_type *nS_over_V;
extern vector *center_V;



/* Populations */
extern int NPOP;
extern pop *p;
extern pop *p_old;
extern vector *v,*vold;
extern my_double *dens;
extern vector *force;

/* LB speeds & weights */
extern my_double wgt[19]; 
extern vector c[19];
extern my_double dirp[19] ,inv[19];
extern my_double cs , cs2 , cs4 , twocs2 , twocs4;
extern my_double invcs , invcs2 ,  invcs4 , invtwocs2 , invtwocs4;

/* time */
extern int itime;
extern int max_step;
extern int time_dump_field;
extern char OutDir[256];























