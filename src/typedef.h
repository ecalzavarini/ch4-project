#ifdef NO_XMLHEADERS
 #include "define.h"
#endif

//#define my_double long double 
#define my_double double

/* Useful structures */
typedef struct {
  my_double x;
  my_double y;
  my_double z;
} vector;

/* Useful structures */
typedef struct {
  my_double xx,xy,xz;
  my_double yx,yy,yz;
  my_double zx,zy,zz;
} tensor;

/* Useful structures */
typedef struct {
  my_double xx,xy,xz;
  my_double yy,yz;
  my_double zz;
} sym_tensor;

#ifdef GRID_POP_D2Q9
#define NPOP 9
#endif
#ifdef GRID_POP_D3Q15
#define NPOP 15
#endif
#ifdef GRID_POP_D3Q19
#define NPOP 19
#endif
#ifdef GRID_POP_D3Q27
#define NPOP 27
#endif


typedef struct {
  my_double p[NPOP];
} pop;

typedef struct {
  my_double val;
  my_double name;
} param;

typedef struct {
  my_double NX , NY , NZ;
  my_double SX , SY , SZ; 
  my_double time_dt, time_start, time_end, time_dump_field, time_dump_diagn;
#ifdef LB_FLUID
  my_double tau_u , nu , nu_add;
 #ifdef LB_FLUID_BC
  #ifdef LB_FLUID_BC_Y
   #ifdef LB_FLUID_BC_Y_M_VELOCITY
    my_double ym_wall_velocity_x , ym_wall_velocity_z ;
   #endif
   #ifdef LB_FLUID_BC_Y_P_VELOCITY
    my_double yp_wall_velocity_x , yp_wall_velocity_z ;
   #endif
  #endif
 #endif
 #ifdef LB_FLUID_INITIAL_ADD_NOISE
  my_double noise_u;
 #endif
 #ifdef LB_FLUID_FORCING
  my_double Amp_x,Amp_y,Amp_z;
  #ifdef LB_FLUID_FORCING_GRAVITY
   my_double gravity_x,gravity_y,gravity_z;
  #endif
 #endif
 #ifdef LB_TEMPERATURE
  my_double tau_t, kappa , kappa_add;
  my_double T_bot,T_top,T_ref,T_ref2,deltaT;
  my_double grad_T_bot,grad_T_top;
  #ifdef LB_TEMPERATURE_BUOYANCY
  my_double beta_t,beta2_t;
  //  my_double gravity_x,gravity_y,gravity_z;
  #endif
  #ifdef LB_TEMPERATURE_INITIAL_ADD_NOISE
  my_double noise_t;
  #endif
  #ifdef LB_TEMPERATURE_FORCING
  my_double Amp_t;
   #ifdef LB_TEMPERATURE_FORCING_RADIATION
  my_double attenuation; 
   #endif
  #endif
  #ifdef LB_TEMPERATURE_MELTING
  my_double T_solid,latent_heat,specific_heat;
   #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID
    #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_SEMISPHERE
     my_double cavity_radius; 
    #endif
    #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_CAVITY
     my_double cavity_SX,cavity_SY,cavity_SZ; 
    #endif
    #ifdef LB_TEMPERATURE_MELTING_INITIAL_LIQUID_LAYER
     my_double layer_height;
    #endif
   #endif
   #ifdef LB_TEMPERATURE_MELTING_SOLID_DIFFUSIVITY
    my_double tau_solid, kappa_solid;
   #endif
  #endif
 #endif
 #ifdef LB_SCALAR
  my_double tau_s, chi, chi_add;
  my_double S_bot,S_top,S_ref,deltaS;
  #ifdef LB_SCALAR_BUOYANCY
  my_double beta_s;
  #endif
  #ifdef LB_SCALAR_FORCING
  my_double Amp_s;
  #endif
 #endif
#endif

#ifdef LAGRANGE
  my_double time_dump_lagr;
  my_double particle_number;
  my_double particle_types; 
  my_double fluid_tracers;
  my_double tau_drag_types, tau_drag_min , tau_drag_max;
  #ifdef LAGRANGE_RADIUSandDENSITY
   my_double particle_radius_types, particle_radius_min , particle_radius_max;
   my_double particle_density_types, particle_density_min , particle_density_max;
  #endif
  #ifdef LAGRANGE_GRAVITY
   #ifdef LAGRANGE_GRAVITY_VARIABLE
   my_double gravity_coeff_types, gravity_coeff_min , gravity_coeff_max;
   #endif
  #endif
 #ifdef LAGRANGE_GRADIENT
  #ifdef LAGRANGE_ADDEDMASS
  my_double beta_coeff_types, beta_coeff_min , beta_coeff_max;
  #endif
  #ifdef LAGRANGE_ORIENTATION
   #ifdef LAGRANGE_ORIENTATION_JEFFREY
   my_double aspect_ratio_types, aspect_ratio_min , aspect_ratio_max;
    #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
    my_double gyrotaxis_velocity_types, gyrotaxis_velocity_min, gyrotaxis_velocity_max;
    #endif
   #endif
   #ifdef LAGRANGE_ORIENTATION_DIFFUSION
   my_double rotational_diffusion_types, rotational_diffusion_min, rotational_diffusion_max;
   #endif
   #ifdef LAGRANGE_ORIENTATION_ACTIVE
   my_double swim_velocity_types, swim_velocity_min, swim_velocity_max;
    #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
    my_double jump_time_types, jump_time_min, jump_time_max;
    my_double critical_shear_rate_types, critical_shear_rate_min, critical_shear_rate_max;
    #endif
   #endif
  #endif /* LAGRANGE_ORIENTATION */ 
  #ifdef LAGRANGE_POLYMER
   my_double tau_polymer, nu_polymer;
  #endif /* LAGRANGE_POLYMER */
 #endif /* LAGRANGE_GRADIENT */
#endif /* LAGRANGE */
} prop;


/* useful macro */
#define SIZE_STRUCT(x) (sizeof(x) / sizeof(my_double))

#define NPROP SIZE_STRUCT(prop)


typedef struct {
#ifdef LB_FLUID
  my_double x,y,z;
  my_double ux,uy,uz;
  my_double ux2,uy2,uz2;
  my_double rho,ene,eps;
  my_double rho2;
#endif
#ifdef LB_TEMPERATURE
  my_double t,t2,epst;
  my_double dxt,dyt,dzt;
  my_double uxt,uyt,uzt;
  my_double nux,nuy,nuz;
  my_double lb;
#ifdef LB_TEMPERATURE_MELTING	
  my_double lf, lf2, dtlf, enth;
#endif		
#endif
#ifdef LB_SCALAR
  my_double s,s2,epss;
  my_double dxs,dys,dzs;
  my_double uxs,uys,uzs;
  my_double nusx,nusy,nusz;
#endif	

} output;

#define NOUT SIZE_STRUCT(output)

#ifndef METHOD_LOG
/* WARNING mvx means rho*vx momentum ,  m means rho density*/

#ifdef GRID_POP_D2Q9
#define mvx(a) (a.p[7] +a.p[6] +a.p[5] -a.p[1] -a.p[2] -a.p[3])
#define mvy(a) (a.p[1] +a.p[8] +a.p[7] -a.p[3] -a.p[4] -a.p[5])
#define mvz(a) 0.0
#define m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8])
#endif

#ifdef GRID_POP_D3Q15
#define mvx(a) (a.p[11] +a.p[12] +a.p[8]  +a.p[13]  +a.p[14]  -a.p[7] -a.p[5] -a.p[4] -a.p[6] -a.p[1])
#define mvy(a) (a.p[7] +a.p[11] +a.p[12]  +a.p[6]  +a.p[9]  -a.p[5] -a.p[4] -a.p[14] -a.p[13] -a.p[2])
#define mvz(a) (a.p[7] +a.p[11] +a.p[5]  +a.p[13]  +a.p[10] -a.p[4] -a.p[14] -a.p[6] -a.p[12] -a.p[3])
#define m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8]+a.p[9]+a.p[10]+a.p[11]+a.p[12]+a.p[13]+a.p[14])
#endif

#ifdef GRID_POP_D3Q19 
#define mvx(a) (a.p[1] -a.p[2] +a.p[7]  +a.p[8]  -a.p[9]  -a.p[10] +a.p[11] -a.p[12]+a.p[13]-a.p[14])
#define mvy(a) (a.p[3] -a.p[4] +a.p[7]  -a.p[8]  +a.p[9]  -a.p[10] +a.p[15] +a.p[16]-a.p[17]-a.p[18])
#define mvz(a) (a.p[5] -a.p[6] +a.p[11] +a.p[12] -a.p[13] -a.p[14] +a.p[15] -a.p[16]+a.p[17]-a.p[18])
#define m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8]+a.p[9]+a.p[10]+a.p[11]+a.p[12]+a.p[13]+a.p[14]+a.p[15]+a.p[16]+a.p[17]+a.p[18])
#endif

#ifdef GRID_POP_D3Q27
#define mvx(a) (a.p[23] +a.p[17] +a.p[24]  +a.p[20]  +a.p[26] +a.p[18] +a.p[25] +a.p[19] +a.p[14] -a.p[13] -a.p[7] -a.p[11] -a.p[4] -a.p[10] -a.p[6] -a.p[12] -a.p[5] -a.p[1])
#define mvy(a) (a.p[13] +a.p[21] +a.p[23]  +a.p[17]  +a.p[24] +a.p[22] +a.p[12] +a.p[5] +a.p[15] -a.p[11] -a.p[4] -a.p[10] -a.p[8] -a.p[26] -a.p[18] -a.p[25] -a.p[9] -a.p[2])
#define mvz(a) (a.p[13] +a.p[21] +a.p[23]  +a.p[19]  +a.p[25] +a.p[9] +a.p[11] +a.p[7] +a.p[16] -a.p[12] -a.p[22] -a.p[24] -a.p[20] -a.p[26] -a.p[8] -a.p[10] -a.p[6] -a.p[3])
#define m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8]+a.p[9]+a.p[10]+a.p[11]+a.p[12]+a.p[13]+a.p[14]+a.p[15]+a.p[16]+a.p[17]+a.p[18]+ a.p[19]+ a.p[20]+ a.p[21]+ a.p[22]+ a.p[23]+ a.p[24]+ a.p[25]+ a.p[26])
#endif

#else
/* WARNING mvx means rho*vx momentum ,  m means rho density*/
#define mvx(a,b) (exp(a.p[1]*b)-exp(a.p[2]*b)+exp(a.p[7]*b)+exp(a.p[8]*b)-exp(a.p[9]*b)-exp(a.p[10]*b)+exp(a.p[11]*b)-exp(a.p[12]*b)+exp(a.p[13]*b)-exp(a.p[14]*b))
#define mvy(a,b) (exp(a.p[3]*b)-exp(a.p[4]*b)+exp(a.p[7]*b)-exp(a.p[8]*b)+exp(a.p[9]*b)-exp(a.p[10]*b)+exp(a.p[15]*b)+exp(a.p[16]*b)-exp(a.p[17]*b)-exp(a.p[18]*b))
#define mvz(a,b) (exp(a.p[5]*b)-exp(a.p[6]*b)+exp(a.p[11]*b)+exp(a.p[12]*b)-exp(a.p[13]*b)-exp(a.p[14]*b)+exp(a.p[15]*b)-exp(a.p[16]*b)+exp(a.p[17]*b)-exp(a.p[18]*b))
#define m(a,b) (exp(a.p[0]*b)+exp(a.p[1]*b)+exp(a.p[2]*b)+exp(a.p[3]*b)+exp(a.p[4]*b)+exp(a.p[5]*b)+exp(a.p[6]*b)+exp(a.p[7]*b)+exp(a.p[8]*b)+exp(a.p[9]*b)+exp(a.p[10]*b)+exp(a.p[11]*b)+exp(a.p[12]*b)+exp(a.p[13]*b)+exp(a.p[14]*b)+exp(a.p[15]*b)+exp(a.p[16]*b)+exp(a.p[17]*b)+exp(a.p[18]*b))
#endif


#define ROOT (!me)
#define AMIROOT (!me)

/* total frame size */
#ifdef METHOD_MYQUICK
#define BRD 2
#else 
#define BRD 1
#endif
#define TWO_BRD 2*BRD

/* index on the grid */
#define IDXG(i,j,k) ( (int)(k)*(LNYG+TWO_BRD)*(LNXG+TWO_BRD)+(int)(j)*(LNXG+TWO_BRD)+(int)(i) )
/* and index for grid planes */
#define IDXG_XBRD(i,j,k) ( (int)(k)*(LNYG+TWO_BRD)*(BRD)+(int)(j)*(BRD)+(int)(i) )
#define IDXG_YBRD(i,j,k) ( (int)(k)*(BRD)*(LNXG+TWO_BRD)+(int)(j)*(LNXG+TWO_BRD)+(int)(i) )
#define IDXG_ZBRD(i,j,k) ( (int)(k)*(LNYG+TWO_BRD)*(LNXG+TWO_BRD)+(int)(j)*(LNXG+TWO_BRD)+(int)(i) )


/* index on the vertices */
#define IDX(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )

/* and index for vertices planes for BRD thickness slices (faces)*/
#define IDX_XBRD(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(BRD)+(int)(j)*(BRD)+(int)(i) )
#define IDX_YBRD(i,j,k) ( (int)(k)*(BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )
#define IDX_ZBRD(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )

/* and index for vertices planes for unit thickness slice*/
#define IDX_X(j,k) ( (int)(k)*(LNY+TWO_BRD)+(int)(j))
#define IDX_Y(i,k) ( (int)(k)*(LNX+TWO_BRD)+(int)(i))
#define IDX_Z(i,j) ( (int)(j)*(LNX+TWO_BRD)+(int)(i))

/* index for corner cubelets */
#define IDX_CORNER(i,j,k) ( (int)(k)*BRD*BRD+(int)(j)*BRD+(int)(i) )

/* index for corner edges */
#define IDX_EDGE_X(i,j,k) ( (int)(k)*(BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )
#define IDX_EDGE_Y(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(BRD)+(int)(j)*(BRD)+(int)(i) )
#define IDX_EDGE_Z(i,j,k) ( (int)(k)*(BRD)*(BRD)+(int)(j)*(BRD)+(int)(i) )


/* neighbouring index for processors*/
#define IDX_NEXT(i,j,k) ( (int)(k+1)*9+(int)(j+1)*3+(int)(i+1) )


#define two_pi 2*3.14159265359
#define one_pi 3.14159265359

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))



/*  Lagrangian part */
#ifdef LAGRANGE

#define wrap(x,s) ( (x) - floor((x) / s) * s )

/*
typedef struct {
  my_double x,y,z,vx,vy,vz,ax,ay,az;
} point;

typedef struct {
  my_double name, kind;
  point now,old;
} point_particle;
*/

typedef struct {
  my_double name;
  my_double tau_drag;
  //  my_double kind;
  my_double x,y,z;
  my_double vx,vy,vz;
  my_double vx_old,vy_old,vz_old;
  my_double ax,ay,az;
  my_double ax_old,ay_old,az_old;
#ifdef LB_FLUID
  my_double ux, uy , uz; /* fluid velocity */
  my_double ux_old, uy_old , uz_old; 
 #ifdef LAGRANGE_GRAVITY
  #ifdef LAGRANGE_GRAVITY_VARIABLE
  my_double gravity_coeff;
  #endif
 #endif
 #ifdef LAGRANGE_GRADIENT
  my_double dx_ux,dy_ux,dz_ux;
  my_double dx_uy,dy_uy,dz_uy;
  my_double dx_uz,dy_uz,dz_uz;
  #ifdef LAGRANGE_ADDEDMASS
   my_double beta_coeff; 
  #endif
  #ifdef LAGRANGE_ORIENTATION
  my_double px,py,pz;
  my_double dt_px,dt_py,dt_pz;
   #ifdef LAGRANGE_ORIENTATION_JEFFREY
   my_double aspect_ratio;
    #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
    my_double gyrotaxis_velocity;
    #endif
   #endif 
   #ifdef LAGRANGE_ORIENTATION_DIFFUSION
   my_double rotational_diffusion;
   #endif
   #ifdef LAGRANGE_ORIENTATION_ACTIVE
   my_double swim_velocity;
    #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
    my_double shear_rate,critical_shear_rate;
    my_double time_from_jump,jump_time;
    my_double px_jump,py_jump,pz_jump;
    #endif
   #endif
  #endif
  #ifdef LAGRANGE_POLYMER
  my_double cxx,cyy,czz,cxy,cyz,cxz;
  my_double dt_cxx,dt_cyy,dt_czz,dt_cxy,dt_cyz,dt_cxz;
  #endif
 #endif
#endif
#ifdef LB_TEMPERATURE
  my_double t;  /* temperature value */
  my_double t_old;
  my_double dt_t;
 #ifdef LAGRANGE_GRADIENT
  my_double dx_t,dy_t,dz_t;
 #endif
#endif
#ifdef LB_SCALAR
  my_double s;  /*scalar value */
  my_double s_old;
  my_double dt_s;
 #ifdef LAGRANGE_GRADIENT
  my_double dx_s,dy_s,dz_s;
 #endif
#endif

} point_particle;

#define SIZE_OF_POINT_PARTICLE SIZE_STRUCT(point_particle)

typedef struct {
  my_double vx,vy,vz;
  my_double vx2,vy2,vz2;
  my_double vx4,vy4,vz4;
  my_double ax,ay,az;
  my_double ax2,ay2,az2;
  my_double ax4,ay4,az4;
 #ifdef LAGRANGE_GRADIENT
  #ifdef LAGRANGE_ORIENTATION
  my_double dt_px,dt_py,dt_pz;
  my_double dt_px2,dt_py2,dt_pz2;
  my_double dt_px4,dt_py4,dt_pz4;
  #endif
 #endif
 #ifdef LB_TEMPERATURE
  my_double t,t2,t4; 
  my_double dt_t,dt_t2;
 #endif
} output_particle;

#define SIZE_OF_OUTPUT_PARTICLE SIZE_STRUCT(output_particle)

#endif /* LAGRANGE */
