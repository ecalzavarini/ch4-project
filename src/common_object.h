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
extern int me;
extern int mex,mey,mez;
extern int me_xp, me_xm, me_yp,me_ym, me_zp,me_zm;
extern int  me_xp_yp_zp, me_xm_ym_zm , me_xp_yp_zm , me_xm_ym_zp , me_xp_ym_zp , me_xm_yp_zm , me_xm_yp_zp , me_xp_ym_zm; 
extern int  me_xp_yp, me_xm_ym, me_xp_ym , me_xm_yp , me_yp_zp , me_ym_zm , me_yp_zm , me_ym_zp , me_xp_zp , me_xm_zm , me_xp_zm , me_xm_zp ;
extern int *me_next;
extern int nprocs;
extern int nxprocs, nyprocs, nzprocs;
extern MPI_Datatype MPI_property_type , MPI_pop_type , MPI_vector_type, MPI_output_type,MPI_my_double_type;
extern MPI_Op MPI_SUM_output;
extern MPI_Op MPI_SUM_vector;

extern MPI_Datatype MPI_pop_plane_x,MPI_pop_plane_y,MPI_pop_plane_z;
extern MPI_Datatype MPI_my_double_plane_x,MPI_my_double_plane_y,MPI_my_double_plane_z;
extern MPI_Datatype MPI_vector_plane_x,MPI_vector_plane_y,MPI_vector_plane_z;

#ifdef LAGRANGE
extern MPI_Datatype MPI_point_particle_type;
extern MPI_Datatype MPI_output_particle_type;
extern MPI_Op MPI_SUM_output_particle;
#endif

/* random seed */
#ifdef RANDOM48
extern unsigned int seed;
#else
extern long * idum;
extern long initdum;
#endif

/* memory allocation size */
extern my_double memory_local;
extern my_double memory_all;
extern my_double memory_max;

/* resume */
extern int resume;
#ifdef LB_FLUID
extern int resume_u;
#endif
#ifdef LB_TEMPERATURE
extern int resume_t;
#endif
#ifdef LB_SCALAR
extern int resume_s;
#endif


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

/* for intial conditions */
extern my_double *channel_u, *channel_y;

/* Array of outputs */
extern output out_local, out_all;
extern output *ruler_x_local, *ruler_y_local, *ruler_z_local;
extern output *ruler_x, *ruler_y, *ruler_z;
extern output *ruler_x_running, *ruler_y_running, *ruler_z_running;

/* mesh */
extern vector *mesh, *center_V;
extern my_double *grid_ruler_x , *grid_ruler_y ,*grid_ruler_z;
extern int *mesh_flag;
extern vector *xp_mesh,*xm_mesh,*yp_mesh,*ym_mesh,*zp_mesh,*zm_mesh;
extern int *xp_flag,*xm_flag,*yp_flag,*ym_flag,*zp_flag,*zm_flag;

#if (defined METHOD_CENTERED || defined METHOD_MYQUICK || defined METHOD_STREAMING || defined METHOD_UPWIND)
extern my_double *interp_xp,*interp_xm,*interp_yp,*interp_ym,*interp_zp,*interp_zm; 
extern my_double *xp_scalar,*xm_scalar,*yp_scalar,*ym_scalar,*zp_scalar,*zm_scalar;
extern vector *xp_vector,*xm_vector,*yp_vector,*ym_vector,*zp_vector,*zm_vector;
#endif
#if (defined METHOD_MYQUICK  || defined METHOD_UPWIND_SKEW)
extern my_double *interp2_xp,*interp2_xm,*interp2_yp,*interp2_ym,*interp2_zp,*interp2_zm;
extern my_double *interp3_xp,*interp3_xm,*interp3_yp,*interp3_ym,*interp3_zp,*interp3_zm;
extern my_double *interp4_xp,*interp4_xm,*interp4_yp,*interp4_ym,*interp4_zp,*interp4_zm;
#endif
#ifdef METHOD_UPWIND_LINEAR 
extern my_double *interp5_xp,*interp5_xm,*interp5_yp,*interp5_ym,*interp5_zp,*interp5_zm;
extern my_double *interp6_xp,*interp6_xm,*interp6_yp,*interp6_ym,*interp6_zp,*interp6_zm;
#endif


extern pop *coeff_xp, *coeff_xm, *coeff_yp, *coeff_ym, *coeff_zp, *coeff_zm;
extern pop *xp_pop,*xm_pop,*yp_pop,*ym_pop,*zp_pop,*zm_pop;


#ifdef LB_FLUID_BC
extern pop *norm_xp_pop,*norm_xm_pop,*norm_yp_pop,*norm_ym_pop,*norm_zp_pop,*norm_zm_pop;
#endif

#ifdef METHOD_EDGES_AND_CORNERS
/* 8 corners */
extern pop *xp_yp_zp_corner_pop;
extern pop *xp_yp_zm_corner_pop;
extern pop *xp_ym_zp_corner_pop;
extern pop *xp_ym_zm_corner_pop;
extern pop *xm_yp_zp_corner_pop;
extern pop *xm_yp_zm_corner_pop;
extern pop *xm_ym_zp_corner_pop;
extern pop *xm_ym_zm_corner_pop;

/* 12 edges */
extern pop *xp_yp_edge_pop;
extern pop *xp_ym_edge_pop;
extern pop *xm_yp_edge_pop;
extern pop *xm_ym_edge_pop;

extern pop *xp_zp_edge_pop;
extern pop *xp_zm_edge_pop;
extern pop *xm_zp_edge_pop;
extern pop *xm_zm_edge_pop;

extern pop *yp_zp_edge_pop;
extern pop *yp_zm_edge_pop;
extern pop *ym_zp_edge_pop;
extern pop *ym_zm_edge_pop;

/* And now the same for vector */
/* 8 corners */
extern vector *xp_yp_zp_corner_vector;
extern vector *xp_yp_zm_corner_vector;
extern vector *xp_ym_zp_corner_vector;
extern vector *xp_ym_zm_corner_vector;
extern vector *xm_yp_zp_corner_vector;
extern vector *xm_yp_zm_corner_vector;
extern vector *xm_ym_zp_corner_vector;
extern vector *xm_ym_zm_corner_vector;

/* 12 edges */
extern vector *xp_yp_edge_vector;
extern vector *xp_ym_edge_vector;
extern vector *xm_yp_edge_vector;
extern vector *xm_ym_edge_vector;

extern vector *xp_zp_edge_vector;
extern vector *xp_zm_edge_vector;
extern vector *xm_zp_edge_vector;
extern vector *xm_zm_edge_vector;

extern vector *yp_zp_edge_vector;
extern vector *yp_zm_edge_vector;
extern vector *ym_zp_edge_vector;
extern vector *ym_zm_edge_vector;

/* And now the same for scalar */
/* 8 corners */
extern my_double *xp_yp_zp_corner_scalar;
extern my_double *xp_yp_zm_corner_scalar;
extern my_double *xp_ym_zp_corner_scalar;
extern my_double *xp_ym_zm_corner_scalar;
extern my_double *xm_yp_zp_corner_scalar;
extern my_double *xm_yp_zm_corner_scalar;
extern my_double *xm_ym_zp_corner_scalar;
extern my_double *xm_ym_zm_corner_scalar;

/* 12 edges */
extern my_double *xp_yp_edge_scalar;
extern my_double *xp_ym_edge_scalar;
extern my_double *xm_yp_edge_scalar;
extern my_double *xm_ym_edge_scalar;

extern my_double *xp_zp_edge_scalar;
extern my_double *xp_zm_edge_scalar;
extern my_double *xm_zp_edge_scalar;
extern my_double *xm_zm_edge_scalar;

extern my_double *yp_zp_edge_scalar;
extern my_double *yp_zm_edge_scalar;
extern my_double *ym_zp_edge_scalar;
extern my_double *ym_zm_edge_scalar;

#endif


/* Populations */
//extern int NPOP;

/* LB speeds & weights */
extern my_double wgt[NPOP]; 
extern vector c[NPOP];
extern int dirp[NPOP] ,inv[NPOP] , inv_x[NPOP], inv_y[NPOP], inv_z[NPOP];
extern my_double twocs2 , twocs4;
extern my_double invtwocs2 , invtwocs4;
extern my_double cs, cs2 , cs4 , cs22 , cssq;
extern my_double invcs, invcs2, invcs4;

#ifdef LB_FLUID
extern pop *p, *rhs_p, *old_rhs_p, *old_old_rhs_p, *p_eq;
extern my_double *dens;
extern vector *u;
 #ifdef LB_FLUID_PAST
 extern my_double *old_dens;
 extern vector *old_u;
 #endif
 #ifdef METHOD_REDEFINED_POP
 extern pop *f_aux;
 #endif
 #ifdef LB_FLUID_FORCING
 extern vector *force;
  #ifdef LB_FLUID_FORCING_LANDSCAPE
  extern my_double *landscape;
  #endif
  #ifdef LB_FLUID_FORCING_HIT
  extern int nk; 
  extern vector *vk, *phi_u;
  extern my_double *vk2;
  extern int randomization_itime;
  #endif
 #endif

 #ifdef LB_FLUID_LES
  #ifdef LB_FLUID_LES_SISM
  extern vector *u_mean;
   #ifdef LB_FLUID_LES_SISM_KALMAN
   extern vector *u_mean_kalman_pre;
   extern vector *sqr_var;
   extern vector *sqr_var_kalman;
   extern vector *K_kalman;
   extern vector *P_kalman;
   extern vector *P_kalman_pre;
   #endif
  #endif
 #endif
#endif

#ifdef LB_TEMPERATURE
extern pop *g, *rhs_g, *old_rhs_g, *old_old_rhs_g , *g_eq;
extern my_double *t;
 #ifdef LB_TEMPERATURE_PAST
 extern my_double *old_t;
 #endif
 #ifdef LB_TEMPERATURE_FORCING
 extern my_double *t_source;
   #ifdef LB_TEMPERATURE_FORCING_PAST
 extern my_double *old_t_source;
   #endif
  #ifdef LB_TEMPERATURE_FORCING_HIT
  extern int nk_t;
  extern vector *vk_t, *phi_t;
  extern my_double *vk2_t;
  extern int randomization_itime_t;
  #endif
  #ifdef LB_TEMPERATURE_MELTING
  extern my_double *liquid_frac, *liquid_frac_old;
  //extern my_double *enthalpy;
  #endif
 #endif
#endif

#ifdef LB_SCALAR
extern pop *h, *rhs_h, *old_rhs_h, *old_old_rhs_h, *h_eq;
extern my_double *s;
 #ifdef LB_SCALAR_FORCING
 extern my_double *s_source;
  #ifdef LB_SCALAR_FORCING_HIT
  extern int nk_s;
  extern vector *vk_s, *phi_s;
  extern my_double *vk2_s;
  extern int randomization_itime_s;
  #endif
 #endif
#endif



/* time */
extern int itime;
extern my_double time_now;
extern char OutDir[256];

#ifdef TIMING
extern my_double t1,t2,tick;
#endif

/* From here Lagrangian definitions */
#ifdef LAGRANGE
extern point_particle *tracer , *tracer_here, *tracer_there, *all_tracer_there;
extern int npart;
#endif























