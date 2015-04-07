#!/bin/bash

rm -f measure_l2norm.txt

for i in `seq 1 50`;
do
echo "i = " $i


# take the values from file

NY=`cat param.in |grep "NY"|awk '{print $2}'`
NY_OLD="NY "$NY

SY=`cat param.in |grep "SY"|awk '{print $2}'`

SX=`cat param.in |grep "SX"|awk '{print $2}'`
SX_OLD="SX "$SX

SZ=`cat param.in |grep "SZ"|awk '{print $2}'`
SZ_OLD="SZ "$SZ

time_dt=`cat param.in |grep "time_dt"|awk '{print $2}'`
time_dt_OLD="time_dt "$time_dt

tau_u=`cat param.in |grep "tau_u"|awk '{print $2}'`
tau_u_OLD="tau_u "$tau_u

# modify the values

#NY_TOGO=`echo $NY $i |awk '{if($1<20) print $1+1.; else print $1+10.;}'`
NY_TOGO=`echo $NY $i |awk '{if($1<=20) printf "%d", $1+1.; else printf "%d", $1+4.;}'`
NY_NEW="NY "$NY_TOGO

SX_TOGO=`echo $NY_TOGO $SY |awk '{printf "%e", $2/$1;}'`
SX_NEW="SX "$SX_TOGO

SZ_TOGO=`echo $NY_TOGO $SY |awk '{printf "%e", $2/$1;}'`
SZ_NEW="SZ "$SZ_TOGO

time_dt_TOGO=$SX_TOGO
time_dt_NEW="time_dt "$time_dt_TOGO

tau_u_TOGO=`echo $time_dt_TOGO |awk '{printf "%lf", 0.5+$1/2.0}'`
tau_u_NEW="tau_u "$tau_u_TOGO


echo "NY OLD "$NY_OLD
echo "NY NEW "$NY_NEW

echo "SX OLD "$SX_OLD
echo "SX NEW "$SX_NEW

echo "SZ OLD "$SZ_OLD
echo "SZ NEW "$SZ_NEW

echo "time_dt OLD "$time_dt_OLD
echo "time_dt NEW "$time_dt_NEW

echo "tau_u OLD "$tau_u_OLD
echo "tau_u NEW "$tau_u_NEW


cat param.in | sed s/"$NY_OLD"/"$NY_NEW"/ > aaa
cp aaa param.in

cat param.in | sed s/"$SX_OLD"/"$SX_NEW"/ > aaa
cp aaa param.in

cat param.in | sed s/"$SZ_OLD"/"$SZ_NEW"/ > aaa
cp aaa param.in

cat param.in | sed s/"$time_dt_OLD"/"$time_dt_NEW"/ > aaa
cp aaa param.in

cat param.in | sed s/"$tau_u_OLD"/"$tau_u_NEW"/ > aaa
cp aaa param.in



gnuplot measure_l2norm.gnu

done 
