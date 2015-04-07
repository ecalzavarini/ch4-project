#!/bin/bash

rm -f measure_l2norm.txt

for i in `seq 0 50`;
do
echo "i = " $i

SY=`cat param.in |grep "SY"|awk '{print $2}'`
AMP=`cat param.in|grep "Amp_x" |awk '{print $2}'`
SY_OLD="SY "$SY
AMP_OLD="Amp_x "$AMP


SY_TOGO=`echo $SY $i |awk '{if($2>0) print $1*1.2; else  print $1;}'`
AMP_TOGO=`echo $AMP $i |awk '{if($2>0)print $1/1.2; else print $1;}'`
SY_NEW="SY "$SY_TOGO
AMP_NEW="Amp_x "$AMP_TOGO


SX=`cat param.in |grep "SX"|awk '{print $2}'`
SX_OLD="SX "$SX
SX_TOGO=`echo $SX $i |awk '{if($2>0) print $1*1.2; else  print $1;}'`
SX_NEW="SX "$SX_TOGO


SZ=`cat param.in |grep "SZ"|awk '{print $2}'`
SZ_TOGO=`echo $SZ $i |awk '{if($2>0) print $1*1.2; else  print $1;}'`
SZ_OLD="SZ "$SZ
SZ_NEW="SZ "$SZ_TOGO

time_dt=`cat param.in |grep "time_dt"|awk '{print $2}'`
time_dt_OLD="time_dt "$time_dt

tau_u=`cat param.in |grep "tau_u"|awk '{print $2}'`
tau_u_OLD="tau_u "$tau_u

time_dt_TOGO=$SX_TOGO
time_dt_NEW="time_dt "$time_dt_TOGO

tau_u_TOGO=`echo $time_dt_TOGO |awk '{printf "%lf", 0.5+$1/2.0}'`
tau_u_NEW="tau_u "$tau_u_TOGO


echo "SY OLD "$SY_OLD
echo "SY NEW "$SY_NEW
echo "AMP OLD "$AMP_OLD
echo "AMP NEW "$AMP_NEW

echo "SX OLD "$SX_OLD
echo "SX NEW "$SX_NEW

echo "SZ OLD "$SZ_OLD
echo "SZ NEW "$SZ_NEW

echo "time_dt OLD "$time_dt_OLD
echo "time_dt NEW "$time_dt_NEW

echo "tau_u OLD "$tau_u_OLD
echo "tau_u NEW "$tau_u_NEW

cat param.in | sed s/"$SY_OLD"/"$SY_NEW"/ > aaa
cp aaa param.in
cat param.in | sed s/"$AMP_OLD"/"$AMP_NEW"/ > aaa
cp aaa param.in
cat param.in | sed s/"$SX_OLD"/"$SX_NEW"/ > aaa
cp aaa param.in
cat param.in | sed s/"$SZ_OLD"/"$SZ_NEW"/ > aaa
cp aaa param.in
cat param.in | sed s/"$tau_u_OLD"/"$tau_u_NEW"/ > aaa
cp aaa param.in
cat param.in | sed s/"$time_dt_OLD"/"$time_dt_NEW"/ > aaa
cp aaa param.in

gnuplot measure_l2norm.gnu

done 
