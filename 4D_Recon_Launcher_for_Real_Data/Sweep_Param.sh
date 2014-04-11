param_idx=`seq 1 12`
for i in $param_idx; do
	export PARAM_INDEX=$i
	qsub -N K_16_r_16_Sim -v PARAM_INDEX Single_Node_Job.sub 
   #    qsub -N K_32_R_32_X800_ProjStart24Kset14 -v PARAM_INDEX Multi_Node_Job.sub 
	sleep 0.2
done 


