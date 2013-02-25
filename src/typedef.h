#include "define.h"

//#define my_double long double 
#define my_double double

/* Useful structures */
typedef struct {
  my_double x;
  my_double y;
  my_double z;
} vector;

typedef struct {
  my_double p[19];
} pop;

typedef struct {
  my_double val;
  my_double name;
} param;

typedef struct {
  double NX , NY , NZ;
  double time_dt, time_max, time_dump_field, time_dump_diagn;
#ifdef LB_FLUID
  double tau_u , nu;
#ifdef LB_FLUID_FORCING_POISEUILLE
  double gradP;
#endif
#endif
} prop;

typedef struct {
#ifdef LB_FLUID
  my_double ux,uy,uz,ux2,uy2,uz2;
  my_double rho,ene,eps;
#endif
} output;



#define NPROP 7
#ifdef LB_FLUID
#define NPROP 9
#ifdef LB_FLUID_FORCING_POISEUILLE
#define NPROP 10
#endif
#endif




/* WARNING mvx means rho*vx momentum ,  m means rho density*/
#define mvx(a) (a.p[1] -a.p[2] +a.p[7]  +a.p[8]  -a.p[9]  -a.p[10] +a.p[11] -a.p[12]+a.p[13]-a.p[14])
#define mvy(a) (a.p[3] -a.p[4] +a.p[7]  -a.p[8]  +a.p[9]  -a.p[10] +a.p[15] +a.p[16]-a.p[17]-a.p[18])
#define mvz(a) (a.p[5] -a.p[6] +a.p[11] +a.p[12] -a.p[13] -a.p[14] +a.p[15] -a.p[16]+a.p[17]-a.p[18])
#define m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8]+a.p[9]+a.p[10]+a.p[11]+a.p[12]+a.p[13]+a.p[14]+a.p[15]+a.p[16]+a.p[17]+a.p[18])

#define ROOT (!me)

/* total frame size */
#define BRD 1
#define TWO_BRD 2*BRD

/* index on the grid */
#define IDXG(i,j,k) ( (int)(k)*(LNYG+TWO_BRD)*(LNXG+TWO_BRD)+(int)(j)*(LNXG+TWO_BRD)+(int)(i) )
/* and index for grid planes */
#define IDXG_XBRD(i,j,k) ( (int)(k)*(LNYG+TWO_BRD)*(BRD)+(int)(j)*(BRD)+(int)(i) )
#define IDXG_YBRD(i,j,k) ( (int)(k)*(BRD)*(LNXG+TWO_BRD)+(int)(j)*(LNXG+TWO_BRD)+(int)(i) )
#define IDXG_ZBRD(i,j,k) ( (int)(k)*(LNYG+TWO_BRD)*(LNXG+TWO_BRD)+(int)(j)*(LNXG+TWO_BRD)+(int)(i) )


/* index on the vertices */
#define IDX(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )
/* and index for vertices planes */
#define IDX_XBRD(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(BRD)+(int)(j)*(BRD)+(int)(i) )
#define IDX_YBRD(i,j,k) ( (int)(k)*(BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )
#define IDX_ZBRD(i,j,k) ( (int)(k)*(LNY+TWO_BRD)*(LNX+TWO_BRD)+(int)(j)*(LNX+TWO_BRD)+(int)(i) )


#define two_pi 2*3.14159265359
