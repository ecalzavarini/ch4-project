#include "common_object.h"

void output_averages(int tstep){
  FILE *fout;
  char fname[128];
  int x, y, yp,ym;
  my_double Nusselt; 
  my_double nu = 0.0;
  my_double kappa = 0.0;  
  my_double tmp1,tmp2,tmp3, tmp4, tmp5, vt, v2, t1, t2 , dyt;
  my_double vx_y[NY] ,  vy_y[NY];
  my_double vx2_y[NY], vy2_y[NY];
  my_double rho_y[NY];
  my_double nusselt_y[NY], vyt_y[NY] , dyt_y[NY];
  my_double t_y[NY], t2_y[NY];
  my_double  s_y[NY] , s2_y[NY];
  my_double lf_y[NY];

#ifdef FLUID
  nu = property.nu;
#endif 

  if(flag==0){
#ifdef FLUID
    rbout.vx = rbout.vy = 0.0; 
    rbout.vx2 = rbout.vy2 =0.0; 
    rbout.rho = 0.0;
#endif
    for (y=0; y<NY; y++){
      nusselt_y[y] = 0.0;
      t_y[y] = vx_y[y] = vy_y[y] = 0.0;
      t2_y[y] = vx2_y[y] = vy2_y[y] = 0.0;
      rho_y[y] = 0.0;
    }
  }

  for (y=1; y<NY+1; y++){
    yp=y+1;
    ym=y-1;
    for (x=1; x<NX+1; x++){

#ifdef FLUID      
#ifdef METHOD_FORCING_GUO      
      tmp1 = (vx(p[IDX(y,x)]) + 0.5*force[IDX(y,x)].x)/m(p[IDX(y,x)]); 
      tmp2 = (vy(p[IDX(y,x)]) + 0.5*force[IDX(y,x)].y)/m(p[IDX(y,x)]); 
#else
      tmp1 = vx(p[IDX(y,x)])/m(p[IDX(y,x)]);
      tmp2 = vy(p[IDX(y,x)])/m(p[IDX(y,x)]);     
#endif
#endif


#ifdef FLUID
      rbout.vx2 +=  tmp1*tmp1;  vx2_y[y-1] +=  tmp1*tmp1;
      rbout.vy2 +=  tmp2*tmp2;  vy2_y[y-1] +=  tmp2*tmp2;
      rbout.vx +=  tmp1;        vx_y[y-1] +=  tmp1;
      rbout.vy +=  tmp2;        vy_y[y-1] +=  tmp2;
      rbout.rho += m(p[IDX(y,x)]) ;        rho_y[y-1]  +=  m(p[IDX(y,x)]);   
#endif

    }       
  }

  if(flag==1){
    // fprintf(stderr, "flag %d, tstep %d\n", flag, tstep);
#ifdef FLUID
    rbout.vx2 /= 2.0*NX*NY;
    rbout.vy2 /= 2.0*NX*NY;
    rbout.rho /= 2.0*NX*NY;
    rbout.vx /= 2.0*NX*NY;
    rbout.vy /= 2.0*NX*NY;
#endif

    for (y=1; y<NY+1; y++){
#ifdef FLUID
      vx_y[y-1] /= 2.0*NX;
      vy_y[y-1] /= 2.0*NX;
      rho_y[y-1] /= 2.0*NX;
      vx2_y[y-1] /= 2.0*NX;
      vy2_y[y-1] /= 2.0*NX;
#endif
    }


#ifdef LB_FLUID
    sprintf(fname,"velocity.dat");
    fout = fopen(fname,"a");
    fprintf(fout,"%d %e %e %e %e %e\n",tstep, (double)rbout.vx, (double)rbout.vy,(double)rbout.rho, (double)rbout.vx2 , (double)rbout.vy2);
    fclose(fout);

    sprintf(fname,"velocity_y.dat");
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) fprintf(fout,"%d %e %e %e %e %e\n",y, (double)vx_y[y-1], (double)vy_y[y-1], (double)rho_y[y-1],  (double)vx2_y[y-1], (double)vy2_y[y-1]);
    fclose(fout);
#endif
  }

}

