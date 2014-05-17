param_idx=`seq 0 4`
for i in $param_idx; do
	export PARAM_INDEX=$i
#	qsub -N K_16_r_16_$i -v PARAM_INDEX Conte_SimData_Submit.sh
        qsub -N sls_att_$i -v PARAM_INDEX Conte_RealData_Submit.sh 
	sleep 0.2
done 
