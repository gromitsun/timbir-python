param_idx=`seq 1 25`
for i in $param_idx; do
	export PARAM_INDEX=$i
	qsub -N r_1_K_1_sim -v PARAM_INDEX Single_Node_Job.sub 
#        qsub -N z_270 -v PARAM_INDEX Multi_Node_Job.sub 
	sleep 0.2
done 
