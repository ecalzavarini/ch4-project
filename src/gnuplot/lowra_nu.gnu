nu = 0.5/3.0
H = 50.0
T=H**2./nu

set xl "t/t_D"
set yl "Nu(t)"
set key bot right

nu=1.56*(17080./1708.)**0.296

p [:6000./T] 'temperature_averages.dat.st' u ($1/T):12 t 'STREAMING' w l , 'temperature_averages.dat.fv' u ($1/T):12 ev 6 t 'FINITE VOLUME' w p lc 3 pt 6 , nu