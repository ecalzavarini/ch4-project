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
  int NX , NY , NZ;
} prop;

typedef struct {
  my_double x;
  my_double y;
  my_double z;
} vector;


/* WARNING vx means rho*vx */
#define vx(a) (a.p[1]+a.p[5]+a.p[8]-a.p[3]-a.p[6]-a.p[7])
#define vy(a) (a.p[2]+a.p[5]+a.p[6]-a.p[4]-a.p[7]-a.p[8])
#define  m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8])


#define IDX(i,j,k) ( (int)(k)*(LNY+BY)*(LNX+BX)+(int)(j)*(LNX+BX)+(int)(i) )

#define ROOT (!me)

/* total frame size */
#define BX 2
#define BY 2
#define BZ 2




