#include "common_object.h"

void design_lb(){

  int i,j,k, pp;

  NPOP = 19;

  /* Settings for D3Q19 */
  wgt[ 0]=1./3.;
  wgt[ 1]=1./18.;  wgt[ 2]=1./18.;  wgt[ 3]=1./18.;
  wgt[ 4]=1./18.;  wgt[ 5]=1./18.;  wgt[ 6]=1./18.;
  wgt[ 7]=1./36.;  wgt[ 8]=1./36.;  wgt[ 9]=1./36.;
  wgt[10]=1./36.;  wgt[11]=1./36.;  wgt[12]=1./36.;
  wgt[13]=1./36.;  wgt[14]=1./36.;  wgt[15]=1./36.;
  wgt[16]=1./36.;  wgt[17]=1./36.;  wgt[18]=1./36.;


  /* Lattice speeds, D3Q19 */
  c[0].x = 0.;  c[0].y = 0.;  c[0].z = 0.;
  c[1].x = 1.;  c[1].y = 0.;  c[1].z = 0.;
  c[2].x =-1.;  c[2].y = 0.;  c[2].z = 0.;
  c[3].x = 0.;  c[3].y = 1.;  c[3].z = 0.;
  c[4].x = 0.;  c[4].y =-1.;  c[4].z = 0.;
  c[5].x = 0.;  c[5].y = 0.;  c[5].z = 1.;
  c[6].x = 0.;  c[6].y = 0.;  c[6].z =-1.;
  c[7].x = 1.;  c[7].y = 1.;  c[7].z = 0.;
  c[8].x = 1.;  c[8].y =-1.;  c[8].z = 0.;
  c[9].x =-1.;  c[9].y = 1.;  c[9].z = 0.;
  c[10].x=-1.;  c[10].y=-1.;  c[10].z= 0.;
  c[11].x= 1.;  c[11].y= 0.;  c[11].z= 1.;
  c[12].x=-1.;  c[12].y= 0.;  c[12].z= 1.;
  c[13].x= 1.;  c[13].y= 0.;  c[13].z=-1.;
  c[14].x=-1.;  c[14].y= 0.;  c[14].z=-1.;
  c[15].x= 0.;  c[15].y= 1.;  c[15].z= 1.;
  c[16].x= 0.;  c[16].y= 1.;  c[16].z=-1.;
  c[17].x= 0.;  c[17].y=-1.;  c[17].z= 1.;
  c[18].x= 0.;  c[18].y=-1.;  c[18].z=-1.;


  /* This function provides the opposite velocity
     c_j=inv[i] for velocity c_i */
  inv[ 0] = 0;
  inv[ 1] = 2;
  inv[ 2] = 1;
  inv[ 3] = 4;
  inv[ 4] = 3;
  inv[ 5] = 6;
  inv[ 6] = 5;

  inv[ 7] = 10;
  inv[ 8] = 9;
  inv[ 9] = 8;
  inv[10] = 7;
  inv[11] = 14;
  inv[12] = 13;
  inv[13] = 12;
  inv[14] = 11;
  inv[15] = 18;
  inv[16] = 17;
  inv[17] = 16;
  inv[18] = 15;

  /* speed of sound constants */
  cs=1.0/sqrt(3.0);   invcs = 1.0/cs;
  cs2=(1.0 / 3.0);    invcs2 = 1.0/cs2;
  cs4=(1.0 / 9.0);    invcs4 = 1.0/cs4;
  twocs2=2.0*cs2; invtwocs2 = 1.0/twocs2;
  twocs4=2.0*cs4; invtwocs4 = 1.0/twocs4;


}


pop equilibrium(pop *f, int i, int j, int k) {
  int pp;
  my_double ux, uy,uz;
  my_double rhof;
  my_double cu, u2;
  pop f_eq;

  rhof = m(f[IDX(i,j,k)]);

      ux = u[IDX(i,j,k)].x;
      uy = u[IDX(i,j,k)].y;
      uy = u[IDX(i,j,k)].z;

      u2 = (ux*ux +  uy*uy + uz*uz);
  
      /* equilibrium distribution */
      for (pp=0; pp<9; pp++){
	cu = (c[pp].x*ux + c[pp].y*uy + c[pp].z*uz);
	f_eq.p[pp] = rhof * wgt[pp] * (1.0 + invcs2*cu  + invtwocs4*cu*cu - invtwocs2*u2 );
      }

      return f_eq;    
}


#ifdef TOBEDONE
/**************************************************/
void hydro_fields(int i){
  int x,y;
  my_double tmpx, tmpy;
  FILE *ferr;
  my_double error;

  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {

#ifdef FLUID
      dens[IDX(y,x)] = m(p[IDX(y,x)]);
      #ifdef METHOD_FORCING_GUO
      v[IDX(y,x)].vx = ( vx(p[IDX(y,x)]) + 0.5*force[IDX(y,x)].x )/dens[IDX(y,x)];
      v[IDX(y,x)].vy = ( vy(p[IDX(y,x)]) + 0.5*force[IDX(y,x)].y )/dens[IDX(y,x)];
     /* set to zero after computing velocity */
      force[IDX(y,x)].x = force[IDX(y,x)].y = 0.0;
      #else
      v[IDX(y,x)].vx = vx(p[IDX(y,x)])/dens[IDX(y,x)];
      v[IDX(y,x)].vy = vy(p[IDX(y,x)])/dens[IDX(y,x)];
      #endif
#endif

#ifdef FLUID_POROSITY
      /* following Guo & Zhao , PHYSICAL REVIEW E 66, 036304 (2002) */
#ifdef METHOD_FORCING_GUO
      eps=porosity[IDX(y,x)];
      vx_temp = ( vx(p[IDX(y,x)]) + 0.5*eps*force[IDX(y,x)].x )/dens[IDX(y,x)];
      vy_temp = ( vy(p[IDX(y,x)]) + 0.5*eps*force[IDX(y,x)].y )/dens[IDX(y,x)];
      v_amp = sqrt(vx_temp*vx_temp + vy_temp*vy_temp); 
      c0 = 0.5 + 0.25*porperty.nu*(1-eps*eps)/(eps*eps);
      c1 = 0.0; 	
      v[IDX(y,x)].vx = vx_temp/(c0+sqrt(c0*c0 + c1*v_amp));
      v[IDX(y,x)].vy = vy_temp/(c0+sqrt(c0*c0 + c1*v_amp));
      /* set to zero after computing velocity */
      force[IDX(y,x)].x = force[IDX(y,x)].y = 0.0;
#endif
#endif

#ifdef TEMPERATURE
      tt[IDX(y,x)]    =  t(g[IDX(y,x)]);
#endif

#ifdef SALT
      ss[IDX(y,x)]    =  t(s[IDX(y,x)]);
#endif
    }

  /* check thermalization */
  if (i%500 == 0) {
    ferr  = fopen("error.dat","a");
    error = 0.0;
    for (y=1; y<NY+1; y++) 
      for (x=1; x<NX+1; x++) {
	tmpx = v[IDX(y,x)].vx - vold[IDX(y,x)].vx;
	tmpy = v[IDX(y,x)].vy - vold[IDX(y,x)].vy;  
	error += (tmpx * tmpx + tmpy * tmpy);
      }
    fprintf(ferr,"%d %g\n",i,error);
    fflush(ferr);

    if ((error < 10e-11) && (i!=0)) {
      fprintf(stderr,"Run termalized\n");
      fclose(ferr);
    }
  }/* if */

  for (y=1; y<NY+1; y++) 
    for (x=1; x<NX+1; x++) {
      vold[IDX(y,x)].vx = v[IDX(y,x)].vx;
      vold[IDX(y,x)].vy = v[IDX(y,x)].vy;
#ifdef TEMPERATURE
      ttold[IDX(y,x)] = tt[IDX(y,x)];
#endif
    }
}
#endif
