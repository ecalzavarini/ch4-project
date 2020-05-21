set  format y "%e"
p [] 'velocity_averages_y.dat' u 1:(abs($4-( -4.e-6*($1-50.)**2. + 0.01 ))) w lp