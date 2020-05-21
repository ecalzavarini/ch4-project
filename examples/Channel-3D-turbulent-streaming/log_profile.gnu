GNUTERM = "x11"

nu = ( 0.504762-0.5)/3.  #0.125/3.*4. #0.06/3.
Lhalf= 64. #1344.0/2.0*2.5*4.
NY = 128.
Amp= 0.011165/2.5  #0.015 
#F=2.0*Amp*((4.*nu)*(Lhalf*2.0)**(-2.));
F=(Amp**2.)/Lhalf

b=0.03  #015  #  0.015
f(x)=a*x
fit [0:b*Lhalf] f(x) 'velocity_averages_y.dat' u 1:4 via a
a1=a
fit [0:b*Lhalf] f(x) 'velocity_averages_y.dat' u (2.*Lhalf-$1):4 via a
a2=a
a= 0.5*(a1+a2)
set log x


set grid 
set key bot
p [1:Lhalf*sqrt(a/nu)][0:20] 'velocity_averages_y.dat' u ($1*sqrt(a/nu)):($4/sqrt(nu*a)) w lp , x , (1./0.41)*log(x)+5.8,\
'velocity_averages_y_run.dat' u ($1*sqrt(a/nu)):($4/sqrt(nu*a)) w l lc 2 lw 2,\
'velocity_averages_y.dat' u ((2.*Lhalf-$1)*sqrt(a/nu)):($4/sqrt(nu*a)) w lp,\
'velocity_averages_y_run.dat' u ((2.*Lhalf-$1)*sqrt(a/nu)):($4/sqrt(nu*a)) w l lw 2

print "Re_tau (derivative)", Lhalf*sqrt(a/nu)

print "Re_tau (large scale)", Lhalf*sqrt(F*Lhalf)/nu

print "press return"
pause -1 

unset log x
#p [0:2.*Lhalf*sqrt(F*Lhalf)/nu][:3.5]
p [0:2.*Lhalf*sqrt(a/nu)] 'velocity_averages_y.dat' u ($1*sqrt(a/nu)):(sqrt($7-$4**2.)/sqrt(nu*a)) w lp, '' u ($1*sqrt(a/nu)):(sqrt($8-$5**2.)/sqrt(nu*a)) w lp, '' u ($1*sqrt(a/nu)):(sqrt($9-$6**2.)/sqrt(nu*a)) w lp
rep 'velocity_averages_y_run.dat' u ($1*sqrt(a/nu)):(sqrt($7-$4**2.)/sqrt(nu*a)) w lp, '' u ($1*sqrt(a/nu)):(sqrt($8-$5**2.)/sqrt(nu*a)) w lp, '' u ($1*sqrt(a/nu)):(sqrt($9-$6**2.)/sqrt(nu*a)) w lp





reset

print "press return"
pause -1

#p 'DAT1/velocity_averages.dat' u 1:2 w lp, '' u 1:7 w lp, '' u 1:8 w lp, '' u 1:9 w lp
p 'velocity_averages.dat' u 1:2 ev 10 w l
#, '' u 1:7 w lp, '' u 1:8 w lp, '' u 1:9 w lp
