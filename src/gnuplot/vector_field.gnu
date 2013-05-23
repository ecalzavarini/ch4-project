reset 
set xr[0:128]
set yr[0:64]
set zr[0:1]
set cbr[-1:1]
set view map
set palette @MATLAB
unset key 
set xlabel '{/Helvetica-Oblique x}'
set ylabel '{/Helvetica-Oblique y}'
set cblabel '{/Helvetica-Oblique T}'
set title '{/Helvetica-Oblique 2D Rayleigh Benard Cell, Ra = 2.5 10^{6} , Pr =1}'

#set pm3d interp 10,10
sp '<paste RUN/vel.200 RUN/temp.200' u 1:2:(0):11 w pm3d ,\
'' u 1:2:(0.5):(20.*$4):(20.*$5):(0) ev 2 w vect  lc "black" head size 1.0,20,60 

set terminal postscript eps size 15.0,7.5 enhanced color font 'Helvetica,20' linewidth 1
set output 'vector_field.eps'
rep
!open vector_field.eps

set term x11
rep
