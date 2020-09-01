#! /bin/bash

rlzStart=$1
rlzEnd=$2

for i in `seq $rlzStart $rlzEnd` 
	do
	before=$(date +"%T")

	#-----run realization
	nice -n 10 ./runOneRlz.sh &> log_$i &
	wait

	#-----copy data
	cp -r ../data/test ../data/data_$i
	cp log_$i ../data/data_$i/.
	#-----clean data
	rm log_$i
	rm -r ../data/test/
	wait

	#-----
	after=$(date +"%T")
	echo "Simulation No. $i finished! Start: $before End: $after"
	done

