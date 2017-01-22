#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <mpi.h>
#ifdef NO_XMLHEADERS
 #include "define.h"
#endif
#include "typedef.h"
#include "functions.h"

#include <hdf5.h>

/* MPI variables */
int me;
int mex,mey,mez;
int me_xp, me_xm, me_yp,me_ym, me_zp,me_zm;
int  me_xp_yp_zp, me_xm_ym_zm , me_xp_yp_zm , me_xm_ym_zp , me_xp_ym_zp , me_xm_yp_zm , me_xm_yp_zp , me_xp_ym_zm; 
int  me_xp_yp, me_xm_ym, me_xp_ym , me_xm_yp , me_yp_zp , me_ym_zm , me_yp_zm , me_ym_zp , me_xp_zp , me_xm_zm , me_xp_zm , me_xm_zp ;
int *me_next;
int nprocs;
int nxprocs, nyprocs, nzprocs;
MPI_Datatype MPI_property_type , MPI_pop_type , MPI_vector_type, MPI_output_type, MPI_my_double_type;
MPI_Op MPI_SUM_output;
MPI_Op MPI_SUM_vector;

MPI_Datatype MPI_pop_plane_x,MPI_pop_plane_y,MPI_pop_plane_z;
MPI_Datatype MPI_my_double_plane_x,MPI_my_double_plane_y,MPI_my_double_plane_z;
MPI_Datatype MPI_vector_plane_x,MPI_vector_plane_y,MPI_vector_plane_z;

#ifdef LAGRANGE
MPI_Datatype MPI_point_particle_type;
MPI_Datatype MPI_output_particle_type;
MPI_Op MPI_SUM_output_particle;
#endif

/* random seed */
#ifdef RANDOM48
unsigned int seed;
#else
long * idum;
long initdum;
#endif

/* memory allocation size */
my_double memory_local;
my_double memory_all;
my_double memory_max;

/* resume */
int resume;
#ifdef LB_FLUID
int resume_u;
#endif
#ifdef LB_TEMPERATURE
int resume_t;
#endif
#ifdef LB_SCALAR
int resume_s;
#endif


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

/* for intial conditions */
my_double *channel_u, *channel_y;

/* Array of outputs */
output out_local,out_all;
output *ruler_x_local, *ruler_y_local, *ruler_z_local;
output *ruler_x, *ruler_y, *ruler_z;
output *ruler_x_running, *ruler_y_running, *ruler_z_running;

/* Mesh matrix */
vector *mesh, *center_V;
my_double *grid_ruler_x , *grid_ruler_y ,*grid_ruler_z;
int *mesh_flag;
vector *xp_mesh,*xm_mesh,*yp_mesh,*ym_mesh,*zp_mesh,*zm_mesh;
int *xp_flag,*xm_flag,*yp_flag,*ym_flag,*zp_flag,*zm_flag;

#if (defined METHOD_CENTERED || defined METHOD_MYQUICK || defined METHOD_STREAMING || defined METHOD_UPWIND)
my_double *interp_xp,*interp_xm,*interp_yp,*interp_ym,*interp_zp,*interp_zm;
my_double *xp_scalar,*xm_scalar,*yp_scalar,*ym_scalar,*zp_scalar,*zm_scalar;
vector *xp_vector,*xm_vector,*yp_vector,*ym_vector,*zp_vector,*zm_vector;
#endif
#if (defined METHOD_MYQUICK  || defined METHOD_UPWIND_SKEW)
my_double *interp2_xp,*interp2_xm,*interp2_yp,*interp2_ym,*interp2_zp,*interp2_zm;
my_double *interp3_xp,*interp3_xm,*interp3_yp,*interp3_ym,*interp3_zp,*interp3_zm;
my_double *interp4_xp,*interp4_xm,*interp4_yp,*interp4_ym,*interp4_zp,*interp4_zm;
#endif
#ifdef METHOD_UPWIND_LINEAR
my_double *interp5_xp,*interp5_xm,*interp5_yp,*interp5_ym,*interp5_zp,*interp5_zm;
my_double *interp6_xp,*interp6_xm,*interp6_yp,*interp6_ym,*interp6_zp,*interp6_zm;
#endif

pop *coeff_xp, *coeff_xm, *coeff_yp, *coeff_ym, *coeff_zp, *coeff_zm;
pop *xp_pop,*xm_pop,*yp_pop,*ym_pop,*zp_pop,*zm_pop;
#ifdef LB_FLUID_BC
pop *norm_xp_pop,*norm_xm_pop,*norm_yp_pop,*norm_ym_pop,*norm_zp_pop,*norm_zm_pop;
#endif

#ifdef METHOD_EDGES_AND_CORNERS
/* 8 corners */
pop *xp_yp_zp_corner_pop;
pop *xp_yp_zm_corner_pop;
pop *xp_ym_zp_corner_pop;
pop *xp_ym_zm_corner_pop;
pop *xm_yp_zp_corner_pop;
pop *xm_yp_zm_corner_pop;
pop *xm_ym_zp_corner_pop;
pop *xm_ym_zm_corner_pop;

/* 12 edges */
pop *xp_yp_edge_pop;
pop *xp_ym_edge_pop;
pop *xm_yp_edge_pop;
pop *xm_ym_edge_pop;

pop *xp_zp_edge_pop;
pop *xp_zm_edge_pop;
pop *xm_zp_edge_pop;
pop *xm_zm_edge_pop;

pop *yp_zp_edge_pop;
pop *yp_zm_edge_pop;
pop *ym_zp_edge_pop;
pop *ym_zm_edge_pop;

/* And now the same for vector */
/* 8 corners */
vector *xp_yp_zp_corner_vector;
vector *xp_yp_zm_corner_vector;
vector *xp_ym_zp_corner_vector;
vector *xp_ym_zm_corner_vector;
vector *xm_yp_zp_corner_vector;
vector *xm_yp_zm_corner_vector;
vector *xm_ym_zp_corner_vector;
vector *xm_ym_zm_corner_vector;

/* 12 edges */
vector *xp_yp_edge_vector;
vector *xp_ym_edge_vector;
vector *xm_yp_edge_vector;
vector *xm_ym_edge_vector;

vector *xp_zp_edge_vector;
vector *xp_zm_edge_vector;
vector *xm_zp_edge_vector;
vector *xm_zm_edge_vector;

vector *yp_zp_edge_vector;
vector *yp_zm_edge_vector;
vector *ym_zp_edge_vector;
vector *ym_zm_edge_vector;

/* And now the same for scalar */
/* 8 corners */
my_double *xp_yp_zp_corner_scalar;
my_double *xp_yp_zm_corner_scalar;
my_double *xp_ym_zp_corner_scalar;
my_double *xp_ym_zm_corner_scalar;
my_double *xm_yp_zp_corner_scalar;
my_double *xm_yp_zm_corner_scalar;
my_double *xm_ym_zp_corner_scalar;
my_double *xm_ym_zm_corner_scalar;

/* 12 edges */
my_double *xp_yp_edge_scalar;
my_double *xp_ym_edge_scalar;
my_double *xm_yp_edge_scalar;
my_double *xm_ym_edge_scalar;

my_double *xp_zp_edge_scalar;
my_double *xp_zm_edge_scalar;
my_double *xm_zp_edge_scalar;
my_double *xm_zm_edge_scalar;

my_double *yp_zp_edge_scalar;
my_double *yp_zm_edge_scalar;
my_double *ym_zp_edge_scalar;
my_double *ym_zm_edge_scalar;


#endif



/* Populations */
//int NPOP;

/* LB speeds & weights */
my_double wgt[NPOP]; 
vector c[NPOP];
int dirp[NPOP] ,inv[NPOP] , inv_x[NPOP], inv_y[NPOP], inv_z[NPOP];
my_double twocs2 , twocs4;
my_double invtwocs2 , invtwocs4;
my_double cs, cs2 , cs4 , cs22 , cssq;
my_double invcs, invcs2, invcs4;


#ifdef LB_FLUID
 pop *p, *rhs_p, *old_rhs_p, *old_old_rhs_p, *p_eq;
 my_double *dens;
 vector *u;
 #ifdef LB_FLUID_PAST
  my_double *old_dens;
  vector *old_u;
 #endif
 #ifdef METHOD_REDEFINED_POP
 pop *f_aux;
 #endif
 #ifdef LB_FLUID_FORCING
 vector *force;
  #ifdef LB_FLUID_FORCING_LANDSCAPE
  my_double *landscape;
  #endif
  #ifdef LB_FLUID_FORCING_HIT
  int nk;
  vector *vk, *phi_u;
  my_double *vk2;
  int randomization_itime;
  #endif
 #endif

 #ifdef LB_FLUID_LES
  #ifdef LB_FLUID_LES_SISM
  vector *u_mean;
   #ifdef LB_FLUID_LES_SISM_KALMAN
   vector *u_mean_kalman_pre;
   vector *sqr_var;
   vector *sqr_var_kalman;
   vector *K_kalman;
   vector *P_kalman;
   vector *P_kalman_pre;
   #endif
  #endif
 #endif
#endif


#ifdef LB_TEMPERATURE
 pop *g, *rhs_g, *old_rhs_g, *old_old_rhs_g, *g_eq;
 my_double *t;
 #ifdef LB_TEMPERATURE_PAST
  my_double *old_t;
 #endif
 #ifdef LB_TEMPERATURE_FORCING
  my_double *t_source;
   #ifdef LB_TEMPERATURE_FORCING_PAST
    my_double *old_t_source;
   #endif
  #ifdef LB_TEMPERATURE_FORCING_HIT
  int nk_t;
  vector *vk_t, *phi_t;
  my_double *vk2_t;
  int randomization_itime_t;
  #endif
 #endif
 #ifdef LB_TEMPERATURE_MELTING
 my_double *liquid_frac, *liquid_frac_old;
 //my_double *enthalpy;
 #endif
#endif

#ifdef LB_SCALAR
 pop *h, *rhs_h, *old_rhs_h, *old_old_rhs_h, *h_eq;
 my_double *s;
 #ifdef LB_SCALAR_FORCING
 my_double *s_source;
  #ifdef LB_SCALAR_FORCING_HIT
  int nk_s;
  vector *vk_s, *phi_s;
  my_double *vk2_s;
  int randomization_itime_s;
  #endif
 #endif
#endif


/* time */
int itime;
my_double time_now;
char OutDir[256];

#ifdef TIMING
my_double t1,t2,tick;
#endif


/* From here Lagrangian definitions */
#ifdef LAGRANGE
point_particle *tracer , *tracer_here, *tracer_there, *all_tracer_there;
int npart;
#endif


















