#!/bin/bash

rm -f measure_l2norm.txt

for i in `seq 1 10`;
do
echo "i = " $i


NY=`cat param.in |grep "NY"|awk '{print $2}'`
SY=`cat param.in |grep "SY"|awk '{print $2}'`
AMP=`cat param.in|grep "Amp_x" |awk '{print $2}'`
NY_OLD="NY "$NY
SY_OLD="SY "$SY
AMP_OLD="Amp_x "$AMP


NY_TOGO=`echo $NY $i |awk '{print $1*2.0}'`
SY_TOGO=`echo $SY $i |awk '{print $1*2.0}'`
AMP_TOGO=`echo $AMP $i |awk '{print $1/2.0}'`
NY_NEW="NY "$NY_TOGO
SY_NEW="SY "$SY_TOGO
AMP_NEW="Amp_x "$AMP_TOGO

echo "NY OLD "$NY_OLD
echo "NY NEW "$NY_NEW
echo "SY OLD "$SY_OLD
echo "SY NEW "$SY_NEW
echo "AMP OLD "$AMP_OLD
echo "AMP NEW "$AMP_NEW

cat param.in | sed s/"$SY_OLD"/"$SY_NEW"/ > aaa
cp aaa param.in
cat param.in | sed s/"$NY_OLD"/"$NY_NEW"/ > aaa
cp aaa param.in
cat param.in     | sed s/"$AMP_OLD"/"$AMP_NEW"/ > bbb
cp bbb param.in

gnuplot measure_l2norm.gnu

done 
