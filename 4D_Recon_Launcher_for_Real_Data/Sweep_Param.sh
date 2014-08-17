#!/bin/sh -l

<<<<<<< HEAD
param_idx=`seq 0 11`
#$msg="data_20_lowres"
export MSG="alsi-60-"
=======
param_idx=`seq 0 15`
#$msg="data_20_lowres"
export MSG="zero_mean_"
>>>>>>> 232251ff45fae046b45696cf0aa5c8ee1e8cdcd5
for i in $param_idx; do
	export PARAM_INDEX=$i
	qsub -N $MSG$i -v PARAM_INDEX Conte_RealData_Submit.sh
#	qsub -N $MSG$i -v PARAM_INDEX Conte_SimData_Submit.sh 
	sleep 0.2
done 
