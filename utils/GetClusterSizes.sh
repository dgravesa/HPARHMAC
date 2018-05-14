#!/bin/sh
# Daniel Graves

if [ $# -ne 1 ]
then
	echo usage: "./GetClusterSizes hmac_output"
	exit
fi

nclusters=`head -1 $1 | cut -f 1 -d ' '`
nc_m1=`expr $nclusters - 1`
is=`seq 0 $nc_m1`

for i in $is
do
	echo -e "$i:\t\c"
	grep "^$i$" $1 | wc -l
done

