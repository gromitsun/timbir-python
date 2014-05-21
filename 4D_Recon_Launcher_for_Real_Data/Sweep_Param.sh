param_idx=`seq 0 0`
for i in $param_idx; do
	export PARAM_INDEX=$i
	qsub -N rmse_$i -v PARAM_INDEX Conte_SimData_Submit.sh
#        qsub -N dataset_17_$i -v PARAM_INDEX Conte_RealData_Submit.sh 
	sleep 0.2
done 
