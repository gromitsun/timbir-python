#!/bin/sh -l

param_idx=`seq 0 8`
#$msg="data_20_lowres"
export MSG="data_14_"
for i in $param_idx; do
	export PARAM_INDEX=$i
	qsub -N $MSG$i -v PARAM_INDEX Conte_RealData_Submit.sh
	#qsub -N $MSG$i -v PARAM_INDEX Conte_SimData_Submit.sh 
	sleep 0.2
done 
