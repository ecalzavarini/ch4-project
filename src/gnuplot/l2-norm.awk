#grid=128, Umax =0.05
#th= -(4*Umax/(L^2))*y*(y-L)
BEGIN{
    sum=0.0;
    SY=128.;
}
{
    ux=$4;
    th=-0.00001220703125*($1)*(($1)-SY);
    #printf "%d %e %e\n",NR,ux,th;
    diff=(ux-th)**2.0;
    sum += diff;
    #printf "%d %e\n",NR,sum;
}END{
    printf "grid with %d points, l2 norm is %e\n", NR, sqrt(sum)/NR;
 }