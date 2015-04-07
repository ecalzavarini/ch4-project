#!/bin/bash

make -f Makefile.ch4
rm -f measure_viscosity.txt



for i in `seq 1 10`;
do
echo "i = " $i


NY=`echo "2^"$i | bc`
TAU_OLD=`cat param.in |grep "NY "`
TAU_NEW="NY "$NY
echo "OLD "$TAU_OLD
echo "NEW "$TAU_NEW
cat param.in | sed s/"$TAU_OLD"/"$TAU_NEW"/  > aaa
cp aaa param.in

TAU_OLD=`cat param.in |grep "SY"`
TAU_NEW="SY "$NY
echo "OLD "$TAU_OLD
echo "NEW "$TAU_NEW
cat param.in | sed s/"$TAU_OLD"/"$TAU_NEW"/  > aaa
cp aaa param.in

gnuplot measure_viscosity_fixedRe.gnu

done 
