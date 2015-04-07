#!/bin/bash

rm -f measure_viscosity.txt



for i in `seq 1 30`;
do
echo "i = " $i

TAU=0.005
echo $TAU
TAU_R=`echo "0.05 + "$TAU"*"$i | bc`
TAU_OLD=`cat param.in |grep "tau_u"`
TAU_NEW="tau_u "$TAU_R
echo "OLD "$TAU_OLD
echo "NEW "$TAU_NEW
cat param.in | sed s/"$TAU_OLD"/"$TAU_NEW"/  > aaa
cp aaa param.in
gnuplot measure_viscosity.gnu

done 
