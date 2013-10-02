BEGIN{
    a = -9.75657e-06;
    b = 0.000624421;
    c = 1.36916e-06; 
    #f(x)=a*x*x + b*x+c;
    sum=0.0;
}
{
    ux=$4;
    th=a*$1*$1 + b*$1+c;
    #printf "%d %e %e\n",NR,ux,th;
    diff=(ux-th)**2.0;
    sum += diff;
    #printf "%d %e\n",NR,sum;
}END{
    printf "grid with %d points, l2 norm is %e\n", NR, sqrt(sum)/NR;
}

