set  format y "%e"
p [] 'velocity_averages_y.dat' u 1:(abs($4-( -4.*0.01/(64.**2.)*($1-32.)**2. + 0.01 ))) w lp