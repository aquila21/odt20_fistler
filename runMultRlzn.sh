#! /bin/bash
nCores=1
nRlzn=1
runCase=$1
 
rm -r $1
cd ./data
rm -r gathered_data/
mkdir gathered_data/
cd ..

echo "Start decomposition (nCores: $nCores, nRlzn: $nRlzn)"
let "r=$nRlzn % $nCores" 
x=0
mkdir $runCase
cd ./$runCase
for i in `seq 1 $nCores`
	do
	#-----Create CPU directories
	mkdir CPU_$i
	cp -r ../input ../source ../run ./CPU_$i/.
	mkdir CPU_$i/data
	cd ./CPU_$i/run/

	#-----Compute realizations per CPU
	let "rlzStart=($i-1)*(($nRlzn-$r)/$nCores)+$x + 1 "
    if [ "$i" -le "$r" ]
    then
        let "x=$((x+1))"
    fi
    let "rlzEnd=$i*((($nRlzn-$r)/$nCores))+$x "
	echo "Core $i: $rlzStart till $rlzEnd"
	#-----Start run-scripts with realization indixes
	./run.sh  $rlzStart $rlzEnd  $&> log_run &
	cd ../../.
	done
 


