#include "common_object.h"


void compute_advection(pop * f){

  int i,j,k,pp;
  pop adv;

  /* check this index */
  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 

 	for(pp=0;pp<NPOP;pp++){

#ifdef METHOD_CENTERED
 /* centered difference scheme */
	  adv.p[pp] = coeff_xp[IDX(i,j,k)].p[pp]*0.5*(p[IDX(i+1,j,k)].p[pp] + p[IDX(i,j,k)].p[pp])+ 
                      coeff_xm[IDX(i,j,k)].p[pp]*0.5*(p[IDX(i-1,j,k)].p[pp] + p[IDX(i,j,k)].p[pp])+
                      coeff_yp[IDX(i,j,k)].p[pp]*0.5*(p[IDX(i,j+1,k)].p[pp] + p[IDX(i,j,k)].p[pp])+
                      coeff_ym[IDX(i,j,k)].p[pp]*0.5*(p[IDX(i,j-1,k)].p[pp] + p[IDX(i,j,k)].p[pp])+
                      coeff_zp[IDX(i,j,k)].p[pp]*0.5*(p[IDX(i,j,k+1)].p[pp] + p[IDX(i,j,k)].p[pp])+
                      coeff_zm[IDX(i,j,k)].p[pp]*0.5*(p[IDX(i,j,k-1)].p[pp] + p[IDX(i,j,k)].p[pp]);
#endif
#ifdef METHOD_UPWIND
 /* first order upwind scheme */
 if(coeff_xp[IDX(i,j,k)].p[pp] >= 0.0)
   adv.p[pp] += coeff_xp[IDX(i,j,k)].p[pp]*p[IDX(i,j,k)].p[pp];
 else
   adv.p[pp] += coeff_xp[IDX(i,j,k)].p[pp]*p[IDX(i+1,j,k)].p[pp];

 if(coeff_xm[IDX(i,j,k)].p[pp] >= 0.0)
   adv.p[pp] += coeff_xm[IDX(i,j,k)].p[pp]*p[IDX(i,j,k)].p[pp];
 else
   adv.p[pp] += coeff_xm[IDX(i,j,k)].p[pp]*p[IDX(i-1,j,k)].p[pp];

 if(coeff_yp[IDX(i,j,k)].p[pp] >= 0.0)
   adv.p[pp] += coeff_yp[IDX(i,j,k)].p[pp]*p[IDX(i,j,k)].p[pp];
 else
   adv.p[pp] += coeff_yp[IDX(i,j,k)].p[pp]*p[IDX(i,j+1,k)].p[pp];

 if(coeff_ym[IDX(i,j,k)].p[pp] >= 0.0)
   adv.p[pp] += coeff_ym[IDX(i,j,k)].p[pp]*p[IDX(i,j,k)].p[pp];
 else
   adv.p[pp] += coeff_ym[IDX(i,j,k)].p[pp]*p[IDX(i,j-1,k)].p[pp];

 if(coeff_zp[IDX(i,j,k)].p[pp] >= 0.0)
   adv.p[pp] += coeff_zp[IDX(i,j,k)].p[pp]*p[IDX(i,j,k)].p[pp];
 else
   adv.p[pp] += coeff_zp[IDX(i,j,k)].p[pp]*p[IDX(i,j,k+1)].p[pp];

 if(coeff_zm[IDX(i,j,k)].p[pp] >= 0.0)
   adv.p[pp] += coeff_zm[IDX(i,j,k)].p[pp]*p[IDX(i,j,k)].p[pp];
 else
   adv.p[pp] += coeff_zm[IDX(i,j,k)].p[pp]*p[IDX(i,j,k-1)].p[pp];

#endif

 rhs_p[IDX(i,j,k)].p[pp] = adv.p[pp];

	}/* for pp */
      }/* for i, j , k */
}


void compute_collision(pop * f){
  int i, j, k, pp;
  my_double invtau , one_minus_invtau;
  pop f_eq;

  one_minus_invtau = (1.0 - invtau);

  for(k=BRD;k<LNZ+BRD;k++)
    for(j=BRD;j<LNY+BRD;j++)
      for(i=BRD;i<LNX+BRD;i++){ 
      
	f_eq=equilibrium(f,i,j,k);
      /* collision */
	for (pp=0; pp<NPOP; pp++){
      /* original */
       //	f[IDX(i,j,k)].p[pp] = f[IDX(i,j,k)].p[pp] - invtau * (f[IDX(i,j,k)].p[pp] - f_eq.p[pp]);
      /* compact */
	f[IDX(i,j,k)].p[pp] =  one_minus_invtau * f[IDX(i,j,k)].p[pp] + invtau * f_eq.p[pp];
	}
      }/* i,j,k*/


}




