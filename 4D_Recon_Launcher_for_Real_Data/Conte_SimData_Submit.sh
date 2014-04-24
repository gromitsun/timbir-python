#!/bin/sh -l
# FILENAME:  Conte_SimData_Submit.sh

#PBS -V
#PBS -q standby
#PBS -l nodes=1:ppn=16
#PBS -l walltime=4:00:00

# start top in batch mode, 30 sec interval between reports, in background

#top -b -d 30 -u mohank >& $HOME/top.out &
#tpid=$!          # Save top's process ID

#module load intel
module load python
module load openmpi

cd $PBS_O_WORKDIR
export PARALLEL=1
export OMP_NUM_THREADS=32
#uniq < $PBS_NODEFILE > nodefile

declare -a sigma_s=(1000000 2000000)
declare -a sigma_t=(8000 16000)
s_idx=`expr $PARAM_INDEX / 2`
t_idx=`expr $PARAM_INDEX % 2`

python XT_Main.py --run_setup --run_recon --MBIR_ATT_SIM --num_nodes 1 --Purdue --Path2Data /blah/ --Path2WhitesDarks /blah/ --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/Att_Recon_Float_Debug/ --rot_center 132 --vox_size 0.65 --proj_start 0 --proj_num 1024 --x_width 1024 --recon_x_width 256 --z_start 0 --z_width 4 --recon_z_width 4 --vox_stop_thresh 0.05 --cost_stop_thresh 1 --sigma_s ${sigma_s[$s_idx]} --sigma_t ${sigma_t[$t_idx]} --K 16 --N_theta 256 --r 16 --multres_xy 4 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string proj_start_0 --min_time_btw_views 0 --rotation_speed 100 --ZingerT 4 --maxHU 40000 --minHU 0           

checkjob -v $PBS_JOBID

#kill $tpid
