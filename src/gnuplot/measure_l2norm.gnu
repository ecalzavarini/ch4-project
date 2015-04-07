reset
system("make -f Makefile.ch4")
system("make -f Makefile.ch4 clean-dat")
system("./ch4")

d(x) = ($0 == 0) ? (x1 = x, x1) : (x2 = x1, x1 = x, x1-x2 )

#read parameters from param.in
tau_u = system("grep tau_u param.in|awk '{print $2}'")
time_dt = system("grep time_dt param.in|awk '{print $2}'")
SY = system("grep SY param.in|awk '{print $2}'")
NY = system("grep NY param.in|awk '{print $2}'")
Umax = system("grep Amp_x param.in|awk '{print $2}'")


#parameters
tau_real = tau_u-0.5*time_dt
cs2 = 3.0
nu=tau_real/cs2
Re = Umax*SY/nu
########
# expected u
U(x) = -(4.0*Umax/(SY**2.0))*x*(x-SY)

p [0:SY] 'velocity_averages_y.dat' u 1:4 w lp, U(x)

#pause -1 

######################################


!cat velocity_averages_y.dat| awk -v SY=`grep SY param.in|awk '{print $2}'` -f compute_grid_spacing.awk > bbb



p [0:SY] '<paste bbb velocity_averages_y.dat' u 2:( $1*(($5-U($2))**2.)  )  w lp



stats '<paste bbb velocity_averages_y.dat' u ( $1*(($5-U($2))**2.)  )  
a = STATS_sum

stats '<paste bbb velocity_averages_y.dat' u ( $1*((U($2))**2.)  ) 
b = STATS_sum

err = sqrt(a/b)

print "Re ",Re 
print  "err ", err
dx = (SY*1.0)/(NY*1.0)

print "# dx, err, Re, tau_real, nu, time_dt , tau_u , NY , SY "
print  dx," ",err," ",Re," ",tau_real," ",nu," ",time_dt," ",tau_u," ",NY," ",SY

set print "aaa"
print "# dx, err, Re, tau_real, nu, time_dt , tau_u , NY , SY "
print  dx," ",err," ",Re," ",tau_real," ",nu," ",time_dt," ",tau_u," ",NY," ",SY

system("cat aaa >> measure_l2norm.txt")