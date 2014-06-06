system("make")
system("make clean-dat")
system("./ch4")

#definition of derivative
d2(x,y) = ($0 == 0) ? (x1 = x, y1 = y, 1/0) : (x2 = x1, x1 = x, y2 = y1, y1 = y, (y1-y2)/(x1-x2))

#read parameters from param.in
tau_u = system("grep tau_u param.in|awk '{print $2}'")
time_dt = system("grep time_dt param.in|awk '{print $2}'")
SY = system("grep SY param.in|awk '{print $2}'")
NY = system("grep NY param.in|awk '{print $2}'")

#parameters
tau_real = tau_u-0.5*time_dt

ENE0=0.5*(0.01**2.)  # value of initial enery u_x^2

TIME_IN = 10*time_dt   # we fit from 10 time steps
#Fit over 12 order of magnitudes

MIN = system("tail -1 velocity_averages.dat|awk '{print $7}'")
if(MIN < 0.5*(0.01**2.)*1.e-12){ 
TIME_END  =  system("cat velocity_averages.dat|awk '{ if( $7< 0.5*(0.01**2.)*1.e-12 ){print $1;}}'>aaa ; head -1 aaa")
}else{
TIME_END  = system("grep time_end param.in|awk '{print $2}'")
}
print TIME_END
####
cs2 = 3.0
nu=tau_real/cs2


p [TIME_IN:TIME_END] 'velocity_averages.dat' u 1:(d2($1,log($7)) / (-2.0*((2.*pi)/SY)**2.)  ) w l , nu

f(x) = a 
fit [TIME_IN:TIME_END] f(x) 'velocity_averages.dat' u 1:(d2($1,log($7)) / (-2.0*((2.*pi)/SY)**2.)  ) via a

unset log
p [TIME_IN:TIME_END] 'velocity_averages.dat' u 1:(d2($1,log($7)) / (-2.0*((2.*pi)/SY)**2.)  ) w l lw 3, nu , a lw 1



print "#tau_real, nu_LB , nu_TH , 100.*(nu-a)/nu , time_dt , tau_u , NY , SY , TIME_END"
print  tau_real," ", a," ",nu," ",100.*(a-nu)/nu," ",time_dt," ",tau_u," ",NY," ",SY," ",TIME_END 

set print "aaa"
print "#tau_real, nu_LB , nu_TH , 100.*(nu-a)/nu , time_dt , tau_u , NY , SY , TIME_END"
print  tau_real," ", a," ",nu," ",100.*(a-nu)/nu," ",time_dt," ",tau_u," ",NY," ",SY," ",TIME_END

system("cat aaa >> measure_viscosity.txt")