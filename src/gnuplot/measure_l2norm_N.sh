#!/bin/bash

rm -f measure_l2norm.txt

for i in `seq 1 128`;
do
echo "i = " $i

NY=`cat param.in |grep "NY"|awk '{print $2}'`
NY_OLD="NY "$NY



NY_TOGO=`echo $NY $i |awk '{print $1+1.}'`
NY_NEW="NY "$NY_TOGO


echo "NY OLD "$NY_OLD
echo "NY NEW "$NY_NEW

cat param.in | sed s/"$NY_OLD"/"$NY_NEW"/ > aaa
cp aaa param.in


gnuplot measure_l2norm.gnu

done 
