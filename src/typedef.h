#include "define.h"

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
  my_double tau_u , nu;
#ifdef LB_FLUID_FORCING
  my_double Amp_x,Amp_y,Amp_z;
#endif
#ifdef LB_TEMPERATURE
  my_double tau_t, kappa;
  my_double T_bot,T_top,T_ref,deltaT;
#ifdef LB_TEMPERATURE_BUOYANCY
  my_double beta_t,beta2_t;
  my_double gravity_x,gravity_y,gravity_z;
#endif
#ifdef LB_TEMPERATURE_FORCING
  my_double Amp_t;
#ifdef LB_TEMPERATURE_FORCING_RADIATION
  my_double attenuation; 
#endif
#endif
#ifdef LB_TEMPERATURE_MELTING
  my_double T_solid,latent_heat,specific_heat;
#endif
#endif
#ifdef LB_SCALAR
  my_double tau_s, chi;
#endif
#endif
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
#endif
#ifdef LB_TEMPERATURE
  my_double t,t2,epst;
  my_double dxt,dyt,dzt;
  my_double uxt,uyt,uzt;
  my_double nux,nuy,nuz;
  my_double lb;
#ifdef LB_TEMPERATURE_MELTING	
  my_double lf, dtlf, enth;
#endif		
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
  my_double name, kind;
  my_double x,y,z,vx,vy,vz,ax,ay,az;
  my_double x_old,y_old,z_old,vx_old,vy_old,vz_old;
#ifdef LB_TEMPERATURE
  my_double t;
#endif

} point_particle;

#define SIZE_OF_POINT_PARTICLE SIZE_STRUCT(point_particle)

#endif
