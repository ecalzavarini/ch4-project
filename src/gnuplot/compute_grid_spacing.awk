 BEGIN{
 i=0;
 }
 {
 x[i]=$1; 
 i++;
 }
END{
    for( j=0; j<i; j++){
	    if(j!=0 && j != i) a = 0.5*(x[j+1]+x[j]) - 0.5*(x[j]+x[j-1]);
	    if(j==0 ) a = 0.5*(x[j+1]+x[j]) - 0.0;
	    if(j==i-1) a = SY - 0.5*(x[j]+x[j-1]);
# 	    printf "%d %e %e\n", j, x[j], a;
	    printf "%e\n",a;
    }
 }
