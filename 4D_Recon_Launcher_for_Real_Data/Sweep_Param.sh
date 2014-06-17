#!/bin/sh -l

param_idx=`seq 0 0`
#$msg="data_20_lowres"
export MSG="data_19"
for i in $param_idx; do
	export PARAM_INDEX=$i
	#qsub -N K_16_r_16_$i -v PARAM_INDEX Conte_SimData_Submit.sh
	qsub -I -N $MSG$i -v PARAM_INDEX Conte_RealData_Submit.sh 
	sleep 0.2
done 
