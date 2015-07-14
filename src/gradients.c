#include "common_object.h"

/* Contains several functions to compute the gradients of the hydrodynamic fields */

#ifdef LB_FLUID
/* be careful still not working */
tensor strain_tensor(pop *f,int i, int j, int k){
  int  pp;
  pop f_eq;
  tensor S;

  S.xx = S.xy = S.xz = 0.0;
  S.yx = S.yy = S.yz = 0.0;
  S.zx = S.zy = S.zz = 0.0;

      /* equilibrium distribution */
   f_eq=equilibrium(f,i,j,k);

      for (pp=0; pp<NPOP; pp++){
	S.xx += c[pp].x*c[pp].x*(f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
	S.yy += c[pp].y*c[pp].y*(f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
	S.zz += c[pp].z*c[pp].z*(f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
	S.xy += c[pp].x*c[pp].y*(f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
	S.xz += c[pp].x*c[pp].z*(f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
	S.yz += c[pp].y*c[pp].z*(f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
      }
      S.yx = S.xy;
      S.zz = S.xz;
      S.zy = S.yz;

      //fprintf(stderr,"SXY %e\n",f_eq.p[0]);
      return S;
}
#endif

/* to define the order of the next gradients */
#define FIRST_ORDER_GRAD
//#define SECOND_ORDER_GRAD
#define LAGRANGE_INTERPOLATION_FIX

#ifdef LB_FLUID
tensor gradient_vector(vector *t, int i, int j, int k){

  tensor tens;
  my_double a0,a1,a2,h1,h2;

/* in the bulk , centered 2nd order finite difference */
   tens.xx = ( t[IDX(i+1, j, k)].x - t[IDX(i-1, j, k)].x )/( center_V[IDX(i+1, j, k)].x - center_V[IDX(i-1, j, k)].x );
   tens.xy = ( t[IDX(i, j+1, k)].x - t[IDX(i, j-1, k)].x )/( center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j-1, k)].y );
   tens.xz = ( t[IDX(i, j, k+1)].x - t[IDX(i, j, k-1)].x )/( center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k-1)].z );

   tens.yx = ( t[IDX(i+1, j, k)].y - t[IDX(i-1, j, k)].y )/( center_V[IDX(i+1, j, k)].x - center_V[IDX(i-1, j, k)].x );
   tens.yy = ( t[IDX(i, j+1, k)].y - t[IDX(i, j-1, k)].y )/( center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j-1, k)].y );
   tens.yz = ( t[IDX(i, j, k+1)].y - t[IDX(i, j, k-1)].y )/( center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k-1)].z );

   tens.zx = ( t[IDX(i+1, j, k)].z - t[IDX(i-1, j, k)].z )/( center_V[IDX(i+1, j, k)].x - center_V[IDX(i-1, j, k)].x );
   tens.zy = ( t[IDX(i, j+1, k)].z - t[IDX(i, j-1, k)].z )/( center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j-1, k)].y );
   tens.zz = ( t[IDX(i, j, k+1)].z - t[IDX(i, j, k-1)].z )/( center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k-1)].z );


#ifdef FIRST_ORDER_GRAD
   /* at the x boundaries , one sided 1st order difference*/
if(LNX_START == 0 && i == BRD){
   tens.xx = ( t[IDX(i+1, j, k)].x - t[IDX(i, j, k)].x )/( center_V[IDX(i+1, j, k)].x - center_V[IDX(i, j, k)].x );
   tens.yx = ( t[IDX(i+1, j, k)].y - t[IDX(i, j, k)].y )/( center_V[IDX(i+1, j, k)].x - center_V[IDX(i, j, k)].x );
   tens.zx = ( t[IDX(i+1, j, k)].z - t[IDX(i, j, k)].z )/( center_V[IDX(i+1, j, k)].x - center_V[IDX(i, j, k)].x );
 }
 if(LNX_END == NX && i == LNX+BRD-1){ 
   tens.xx = ( t[IDX(i, j, k)].x - t[IDX(i-1, j, k)].x )/( center_V[IDX(i, j, k)].x - center_V[IDX(i-1, j, k)].x );
   tens.yx = ( t[IDX(i, j, k)].y - t[IDX(i-1, j, k)].y )/( center_V[IDX(i, j, k)].x - center_V[IDX(i-1, j, k)].x );
   tens.zx = ( t[IDX(i, j, k)].z - t[IDX(i-1, j, k)].z )/( center_V[IDX(i, j, k)].x - center_V[IDX(i-1, j, k)].x );
 }

   /* at the y boundaries */
if(LNY_START == 0 && j == BRD){
   tens.xy = ( t[IDX(i, j+1, k)].x - t[IDX(i, j, k)].x )/( center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j, k)].y );
   tens.yy = ( t[IDX(i, j+1, k)].y - t[IDX(i, j, k)].y )/( center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j, k)].y );
   tens.zy = ( t[IDX(i, j+1, k)].z - t[IDX(i, j, k)].z )/( center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j, k)].y ); 
 }
 if(LNY_END == NY && j == LNY+BRD-1){ 
   tens.xy = ( t[IDX(i, j, k)].x - t[IDX(i, j-1, k)].x )/( center_V[IDX(i, j, k)].y - center_V[IDX(i, j-1, k)].y );
   tens.yy = ( t[IDX(i, j, k)].y - t[IDX(i, j-1, k)].y )/( center_V[IDX(i, j, k)].y - center_V[IDX(i, j-1, k)].y );
   tens.zy = ( t[IDX(i, j, k)].z - t[IDX(i, j-1, k)].z )/( center_V[IDX(i, j, k)].y - center_V[IDX(i, j-1, k)].y );
 }

   /* at the z boundaries */
if(LNZ_START == 0 && k == BRD){
   tens.xz = ( t[IDX(i, j, k+1)].x - t[IDX(i, j, k)].x )/( center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k)].z );
   tens.yz = ( t[IDX(i, j, k+1)].y - t[IDX(i, j, k)].y )/( center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k)].z );
   tens.zz = ( t[IDX(i, j, k+1)].z - t[IDX(i, j, k)].z )/( center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k)].z );
 }
 if(LNZ_END == NZ && k == LNZ+BRD-1){ 
   tens.xz = ( t[IDX(i, j, k)].x - t[IDX(i, j, k-1)].x )/( center_V[IDX(i, j, k)].z - center_V[IDX(i, j, k-1)].z );
   tens.yz = ( t[IDX(i, j, k)].y - t[IDX(i, j, k-1)].y )/( center_V[IDX(i, j, k)].z - center_V[IDX(i, j, k-1)].z );
   tens.zz = ( t[IDX(i, j, k)].z - t[IDX(i, j, k-1)].z )/( center_V[IDX(i, j, k)].z - center_V[IDX(i, j, k-1)].z );
 }
#endif

#ifdef SECOND_ORDER_GRAD
   /* at the x boundaries , one sided 2nd order difference*/
   /* but it seems to be less precise than the first order */
if(LNX_START == 0 && i == BRD){
    h1 =  center_V[IDX(i+1, j, k)].x - center_V[IDX(i, j, k)].x;
    h2 =  center_V[IDX(i+2, j, k)].x - center_V[IDX(i+1, j, k)].x; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    tens.xx = a0*t[IDX(i, j, k)].x + a1*t[IDX(i+1, j, k)].x + a2*t[IDX(i+2, j, k)].x;
    tens.yx = a0*t[IDX(i, j, k)].y + a1*t[IDX(i+1, j, k)].y + a2*t[IDX(i+2, j, k)].y;
    tens.zx = a0*t[IDX(i, j, k)].z + a1*t[IDX(i+1, j, k)].z + a2*t[IDX(i+2, j, k)].z;
 }
 if(LNX_END == NX && i == LNX+BRD-1){ 
    h1 =  center_V[IDX(i, j, k)].x - center_V[IDX(i-1, j, k)].x;
    h2 =  center_V[IDX(i-1, j, k)].x - center_V[IDX(i-2, j, k)].x; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    tens.xx = a0*t[IDX(i, j, k)].x + a1*t[IDX(i-1, j, k)].x + a2*t[IDX(i-2, j, k)].x;  
    tens.yx = a0*t[IDX(i, j, k)].y + a1*t[IDX(i-1, j, k)].y + a2*t[IDX(i-2, j, k)].y;  
    tens.zx = a0*t[IDX(i, j, k)].z + a1*t[IDX(i-1, j, k)].z + a2*t[IDX(i-2, j, k)].z;   
 }

   /* at the y boundaries */
if(LNY_START == 0 && j == BRD){
    h1 =  center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j, k)].y;
    h2 =  center_V[IDX(i, j+2, k)].y - center_V[IDX(i, j+1, k)].y; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    tens.xy = a0*t[IDX(i, j, k)].x + a1*t[IDX(i, j+1, k)].x + a2*t[IDX(i, j+2, k)].x;
    tens.yy = a0*t[IDX(i, j, k)].y + a1*t[IDX(i, j+1, k)].y + a2*t[IDX(i, j+2, k)].y;
    tens.zy = a0*t[IDX(i, j, k)].z + a1*t[IDX(i, j+1, k)].z + a2*t[IDX(i, j+2, k)].z;

    //fprintf(stderr,"0 grad.y %e \n",grad.y);
 }
 if(LNY_END == NY && j == LNY+BRD-1){ 
    h1 =  center_V[IDX(i, j, k)].y - center_V[IDX(i, j-1, k)].y;
    h2 =  center_V[IDX(i, j-1, k)].y - center_V[IDX(i, j-2, k)].y; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    tens.xy = a0*t[IDX(i, j, k)].x + a1*t[IDX(i, j-1, k)].x + a2*t[IDX(i, j-2, k)].x;
    tens.yy = a0*t[IDX(i, j, k)].y + a1*t[IDX(i, j-1, k)].y + a2*t[IDX(i, j-2, k)].y;
    tens.zy = a0*t[IDX(i, j, k)].z + a1*t[IDX(i, j-1, k)].z + a2*t[IDX(i, j-2, k)].z; 
    //fprintf(stderr,"50 grad.y %e \n",grad.y);
 }

   /* at the z boundaries */
if(LNZ_START == 0 && k == BRD){
  h1 =  center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k)].z;
  h2 =  center_V[IDX(i, j, k+2)].z - center_V[IDX(i, j, k+1)].z; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    tens.xz = a0*t[IDX(i, j, k)].x + a1*t[IDX(i, j, k+1)].x + a2*t[IDX(i, j, k+2)].x;
    tens.yz = a0*t[IDX(i, j, k)].y + a1*t[IDX(i, j, k+1)].y + a2*t[IDX(i, j, k+2)].y;
    tens.zz = a0*t[IDX(i, j, k)].z + a1*t[IDX(i, j, k+1)].z + a2*t[IDX(i, j, k+2)].z;

 }
 if(LNZ_END == NZ && k == LNZ+BRD-1){ 
   h1 =  center_V[IDX(i, j, k)].z - center_V[IDX(i, j, k-1)].z;
   h2 =  center_V[IDX(i, j, k-1)].z - center_V[IDX(i, j, k-2)].z; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    tens.xz = a0*t[IDX(i, j, k)].x + a1*t[IDX(i, j, k-1)].x + a2*t[IDX(i, j, k-2)].x; 
    tens.yz = a0*t[IDX(i, j, k)].y + a1*t[IDX(i, j, k-1)].y + a2*t[IDX(i, j, k-2)].y; 
    tens.zz = a0*t[IDX(i, j, k)].z + a1*t[IDX(i, j, k-1)].z + a2*t[IDX(i, j, k-2)].z; 
 }
#endif


#ifdef LAGRANGE_INTERPOLATION_FIX
   /* at the x boundaries , one sided 1st order difference*/
 //if(i == 0){
if(i == BRD-1){
   tens.xx = ( t[IDX(i+2, j, k)].x - t[IDX(i+1, j, k)].x )/( center_V[IDX(i+2, j, k)].x - center_V[IDX(i+1, j, k)].x );
   tens.yx = ( t[IDX(i+2, j, k)].y - t[IDX(i+1, j, k)].y )/( center_V[IDX(i+2, j, k)].x - center_V[IDX(i+1, j, k)].x );
   tens.zx = ( t[IDX(i+2, j, k)].z - t[IDX(i+1, j, k)].z )/( center_V[IDX(i+2, j, k)].x - center_V[IDX(i+1, j, k)].x );
 }
// if(i == LNX+TWO_BRD-1){ 
if(i == LNX+BRD){ 
   tens.xx = ( t[IDX(i-1, j, k)].x - t[IDX(i-2, j, k)].x )/( center_V[IDX(i-1, j, k)].x - center_V[IDX(i-2, j, k)].x );
   tens.yx = ( t[IDX(i-1, j, k)].y - t[IDX(i-2, j, k)].y )/( center_V[IDX(i-1, j, k)].x - center_V[IDX(i-2, j, k)].x );
   tens.zx = ( t[IDX(i-1, j, k)].z - t[IDX(i-2, j, k)].z )/( center_V[IDX(i-1, j, k)].x - center_V[IDX(i-2, j, k)].x );
 }

   /* at the y boundaries */
//if(j == 0){
if(j == BRD-1){
   tens.xy = ( t[IDX(i, j+2, k)].x - t[IDX(i, j+1, k)].x )/( center_V[IDX(i, j+2, k)].y - center_V[IDX(i, j+1, k)].y );
   tens.yy = ( t[IDX(i, j+2, k)].y - t[IDX(i, j+1, k)].y )/( center_V[IDX(i, j+2, k)].y - center_V[IDX(i, j+1, k)].y );
   tens.zy = ( t[IDX(i, j+2, k)].z - t[IDX(i, j+1, k)].z )/( center_V[IDX(i, j+2, k)].y - center_V[IDX(i, j+1, k)].y ); 
 }
//if(j == LNY+TWO_BRD-1){
if(j == LNY+BRD){ 
   tens.xy = ( t[IDX(i, j-1, k)].x - t[IDX(i, j-2, k)].x )/( center_V[IDX(i, j-1, k)].y - center_V[IDX(i, j-2, k)].y );
   tens.yy = ( t[IDX(i, j-1, k)].y - t[IDX(i, j-2, k)].y )/( center_V[IDX(i, j-1, k)].y - center_V[IDX(i, j-2, k)].y );
   tens.zy = ( t[IDX(i, j-1, k)].z - t[IDX(i, j-2, k)].z )/( center_V[IDX(i, j-1, k)].y - center_V[IDX(i, j-2, k)].y );
 }

   /* at the z boundaries */
//if(k == 0){
if(k == BRD-1){
   tens.xz = ( t[IDX(i, j, k+2)].x - t[IDX(i, j, k+1)].x )/( center_V[IDX(i, j, k+2)].z - center_V[IDX(i, j, k+1)].z );
   tens.yz = ( t[IDX(i, j, k+2)].y - t[IDX(i, j, k+1)].y )/( center_V[IDX(i, j, k+2)].z - center_V[IDX(i, j, k+1)].z );
   tens.zz = ( t[IDX(i, j, k+2)].z - t[IDX(i, j, k+1)].z )/( center_V[IDX(i, j, k+2)].z - center_V[IDX(i, j, k+1)].z );
 }
// if(k == LNZ+TWO_BRD-1){ 
if(k == LNZ+BRD){ 
   tens.xz = ( t[IDX(i, j, k-1)].x - t[IDX(i, j, k-2)].x )/( center_V[IDX(i, j, k-1)].z - center_V[IDX(i, j, k-2)].z );
   tens.yz = ( t[IDX(i, j, k-1)].y - t[IDX(i, j, k-2)].y )/( center_V[IDX(i, j, k-1)].z - center_V[IDX(i, j, k-2)].z );
   tens.zz = ( t[IDX(i, j, k-1)].z - t[IDX(i, j, k-2)].z )/( center_V[IDX(i, j, k-1)].z - center_V[IDX(i, j, k-2)].z );
 }
#endif

      return tens;
}
#endif


//#ifdef LB_TEMPERATURE
vector gradient_scalar(my_double *t, int i, int j, int k){

  vector grad;
  my_double a0,a1,a2,h1,h2;

  /* in the bulk , centered 2nd order finite difference */
   grad.x = ( t[IDX(i+1, j, k)] - t[IDX(i-1, j, k)] )/( center_V[IDX(i+1, j, k)].x - center_V[IDX(i-1, j, k)].x );
   grad.y = ( t[IDX(i, j+1, k)] - t[IDX(i, j-1, k)] )/( center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j-1, k)].y );
   grad.z = ( t[IDX(i, j, k+1)] - t[IDX(i, j, k-1)] )/( center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k-1)].z );


#ifdef FIRST_ORDER_GRAD
   /* at the x boundaries , one sided 1st order difference*/
if(LNX_START == 0 && i == BRD){
    grad.x = ( t[IDX(i+1, j, k)] - t[IDX(i, j, k)] )/( center_V[IDX(i+1, j, k)].x - center_V[IDX(i, j, k)].x );
 }
 if(LNX_END == NX && i == LNX+BRD-1){ 
   grad.x = ( t[IDX(i, j, k)] - t[IDX(i-1, j, k)] )/( center_V[IDX(i, j, k)].x - center_V[IDX(i-1, j, k)].x );
 }

   /* at the y boundaries */
if(LNY_START == 0 && j == BRD){
   grad.y = ( t[IDX(i, j+1, k)] - t[IDX(i, j, k)] )/( center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j, k)].y );
 }
 if(LNY_END == NY && j == LNY+BRD-1){ 
   grad.y = ( t[IDX(i, j, k)] - t[IDX(i, j-1, k)] )/( center_V[IDX(i, j, k)].y - center_V[IDX(i, j-1, k)].y );
 }

   /* at the z boundaries */
if(LNZ_START == 0 && k == BRD){
   grad.z = ( t[IDX(i, j, k+1)] - t[IDX(i, j, k)] )/( center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k)].z ); 
 }
 if(LNZ_END == NZ && k == LNZ+BRD-1){ 
   grad.z = ( t[IDX(i, j, k)] - t[IDX(i, j, k-1)] )/( center_V[IDX(i, j, k)].z - center_V[IDX(i, j, k-1)].z );
 }
#endif

#ifdef SECOND_ORDER_GRAD
   /* at the x boundaries , one sided 2nd order difference*/
   /* but it seems to be less precise than the first order */
if(LNX_START == 0 && i == BRD){

  h1 =  center_V[IDX(i+1, j, k)].x - center_V[IDX(i, j, k)].x;
  h2 =  center_V[IDX(i+2, j, k)].x - center_V[IDX(i+1, j, k)].x; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    grad.x = a0*t[IDX(i, j, k)] + a1*t[IDX(i+1, j, k)] + a2*t[IDX(i+2, j, k)];

 }
 if(LNX_END == NX && i == LNX+BRD-1){ 

   h1 =  center_V[IDX(i, j, k)].x - center_V[IDX(i-1, j, k)].x;
   h2 =  center_V[IDX(i-1, j, k)].x - center_V[IDX(i-2, j, k)].x; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    grad.x = a0*t[IDX(i, j, k)] + a1*t[IDX(i-1, j, k)] + a2*t[IDX(i-2, j, k)];   
 }

   /* at the y boundaries */
if(LNY_START == 0 && j == BRD){
    h1 =  center_V[IDX(i, j+1, k)].y - center_V[IDX(i, j, k)].y;
    h2 =  center_V[IDX(i, j+2, k)].y - center_V[IDX(i, j+1, k)].y; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    grad.y = a0*t[IDX(i, j, k)] + a1*t[IDX(i, j+1, k)] + a2*t[IDX(i, j+2, k)];
    //fprintf(stderr,"0 grad.y %e \n",grad.y);
 }
 if(LNY_END == NY && j == LNY+BRD-1){ 
    h1 =  center_V[IDX(i, j, k)].y - center_V[IDX(i, j-1, k)].y;
    h2 =  center_V[IDX(i, j-1, k)].y - center_V[IDX(i, j-2, k)].y; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    grad.y = a0*t[IDX(i, j, k)] + a1*t[IDX(i, j-1, k)] + a2*t[IDX(i, j-2, k)]; 
    //fprintf(stderr,"50 grad.y %e \n",grad.y);
 }

   /* at the z boundaries */
if(LNZ_START == 0 && k == BRD){
  h1 =  center_V[IDX(i, j, k+1)].z - center_V[IDX(i, j, k)].z;
  h2 =  center_V[IDX(i, j, k+2)].z - center_V[IDX(i, j, k+1)].z; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    grad.z = a0*t[IDX(i, j, k)] + a1*t[IDX(i, j, k+1)] + a2*t[IDX(i, j, k+2)];
 }
 if(LNZ_END == NZ && k == LNZ+BRD-1){ 
   h1 =  center_V[IDX(i, j, k)].z - center_V[IDX(i, j, k-1)].z;
   h2 =  center_V[IDX(i, j, k-1)].z - center_V[IDX(i, j, k-2)].z; 
    a0 =  (2.0*h1+h2 )/( h1*(h1+h2) );
    a1 = -( h1+h2 )/( h1*h2 );
    a2 =  ( h1 )/( h2*(h1+h2) ); 
    grad.z = a0*t[IDX(i, j, k)] + a1*t[IDX(i, j, k-1)] + a2*t[IDX(i, j, k-2)]; 
 }
#endif

#ifdef LAGRANGE_INTERPOLATION_FIX
   /* at the x boundaries , one sided 1st order difference*/
 //if(i == 0){
if(i == BRD-1){
   grad.x = ( t[IDX(i+2, j, k)] - t[IDX(i+1, j, k)] )/( center_V[IDX(i+2, j, k)].x - center_V[IDX(i+1, j, k)].x );
 }
// if(i == LNX+TWO_BRD-1){ 
 if(i == LNX+BRD){ 
   grad.x = ( t[IDX(i-1, j, k)] - t[IDX(i-2, j, k)] )/( center_V[IDX(i-1, j, k)].x - center_V[IDX(i-2, j, k)].x );
 }

   /* at the y boundaries */
 //if(j == 0){
if(j == BRD-1){
   grad.y = ( t[IDX(i, j+2, k)] - t[IDX(i, j+1, k)] )/( center_V[IDX(i, j+2, k)].y - center_V[IDX(i, j+1, k)].y );
 }
// if(j == LNY+TWO_BRD-1){ 
 if(j == LNY+BRD){ 
   grad.y = ( t[IDX(i, j-1, k)] - t[IDX(i, j-2, k)] )/( center_V[IDX(i, j-1, k)].y - center_V[IDX(i, j-2, k)].y );
 }

   /* at the z boundaries */
 //if(k == 0){
if(k == BRD-1){
   grad.z = ( t[IDX(i, j, k+2)] - t[IDX(i, j, k+1)] )/( center_V[IDX(i, j, k+2)].z - center_V[IDX(i, j, k+1)].z ); 
 }
//if(k == LNZ+TWO_BRD-1){ 
 if(k == LNZ+BRD){ 
   grad.z = ( t[IDX(i, j, k-1)] - t[IDX(i, j, k-2)] )/( center_V[IDX(i, j, k-1)].z - center_V[IDX(i, j, k-2)].z );
 }
#endif

    return grad;
}
//#endif
