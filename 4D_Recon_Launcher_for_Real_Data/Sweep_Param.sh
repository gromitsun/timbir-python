param_idx=`seq 0 3`
for i in $param_idx; do
	export PARAM_INDEX=$i
	qsub -N sim_$i -v PARAM_INDEX Single_Node_Job.sub 
#        qsub -N z_270 -v PARAM_INDEX Multi_Node_Job.sub 
	sleep 0.2
done 
