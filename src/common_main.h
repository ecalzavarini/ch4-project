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
int me;
int mex,mey,mez;
int nprocs;
int nxprocs, nyprocs, nzprocs;

/* System size */
int NX , NY , NZ;
int LNX , LNY , LNZ;
int LNX_END , LNY_END , LNZ_END;
int LNX_START , LNY_START , LNZ_START;

/* Array of properties*/
prop property;

/* Mesh matrix */
mesh_type *mesh;
surf_type *nS_over_V;
vector *center_V;

/* Populations */
int NPOP;
pop *p;
pop *p_old;
vector *v,*vold;
my_double *dens;
vector *force;

/* LB speeds & weights */
my_double wgt[19]; 
vector c[19];
my_double dirp[19] ,inv[19];
my_double cs , cs2 , cs4 , twocs2 , twocs4;
my_double invcs , invcs2 ,  invcs4 , invtwocs2 , invtwocs4;

/* time */
int itime;
int max_step;
int time_dump_field;
char OutDir[256];























