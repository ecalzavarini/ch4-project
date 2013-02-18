#include "define.h"

//#define my_double long double 
#define my_double double

/* Useful structures */
typedef struct {
  my_double x , y , z;
  int flag;
} mesh_type;


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


#define NPROP 7
#ifdef LB_FLUID
#define NPROP 9
#ifdef LB_FLUID_FORCING_POISEUILLE
#define NPROP 10
#endif
#endif

typedef struct {
  my_double x;
  my_double y;
  my_double z;
} vector;


/* WARNING mvx means rho*vx momentum ,  m means rho density*/
#define mvx(a) (a.p[1] -a.p[2] +a.p[7]  +a.p[8]  -a.p[9]  -a.p[10] +a.p[11] -a.p[12]+a.p[13]-a.p[14])
#define mvy(a) (a.p[3] -a.p[4] +a.p[7]  -a.p[8]  +a.p[9]  -a.p[10] +a.p[15] +a.p[16]-a.p[17]-a.p[18])
#define mvz(a) (a.p[5] -a.p[6] +a.p[11] +a.p[12] -a.p[13] -a.p[14] +a.p[15] -a.p[16]+a.p[17]-a.p[18])
#define m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8]+a.p[9]+a.p[10]+a.p[11]+a.p[12]+a.p[13]+a.p[14]+a.p[15]+a.p[16]+a.p[17]+a.p[18])

/* index on the grid */
#define IDXG(i,j,k) ( (int)(k)*(LNY+1+BY)*(LNX+1+BX)+(int)(j)*(LNX+1+BX)+(int)(i) )
/* index on the vertices */
#define IDX(i,j,k) ( (int)(k)*(LNY+BY)*(LNX+BX)+(int)(j)*(LNX+BX)+(int)(i) )

#define ROOT (!me)

/* total frame size */
#define BX 2
#define BY 2
#define BZ 2




