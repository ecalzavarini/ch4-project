reset 
a=200.
set xr[0:1280]
set yr[0:640]
set zr[0:1]
set cbr[-1.0:1.0]
set view map
set palette @MATLAB
unset key 
set xlabel '{/Helvetica-Oblique x}'
set ylabel '{/Helvetica-Oblique y}'
set cblabel '{/Helvetica-Oblique T}'
set title '{/Helvetica-Oblique 2D Rayleigh Benard Cell, Ra = 2.5 10^{6} , Pr =1}'

set pm3d interp 4,4

sp '<paste RUN/vel.14000 RUN/temp.14000' u 1:2:(0):11 w pm3d ,\
'' u 1:2:(0.5):(a*$4):(a*$5):(0) ev 2 w vect  lc "black" head size 1.0,20,60 

set terminal postscript eps size 15.0,7.5 enhanced color font 'Helvetica,20' linewidth 1
set output 'vector_field.eps'
rep
!open vector_field.eps

set term x11
rep
