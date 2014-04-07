param_idx=`seq 1 1`
for i in $param_idx; do
	export PARAM_INDEX=$i
	qsub -N Dataset_03_Res1040 -v PARAM_INDEX Single_Node_Job.sub 
   #    qsub -N K_32_R_32_X800_ProjStart24Kset14 -v PARAM_INDEX Multi_Node_Job.sub 
	sleep 0.2
done 


