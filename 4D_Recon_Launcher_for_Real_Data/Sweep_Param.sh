param_idx=`seq 1 1`
for i in $param_idx; do
	export PARAM_INDEX=$i
	qsub -N param_$i -I -v PARAM_INDEX Single_Node_Job.sub 
	sleep 0.2
done 


