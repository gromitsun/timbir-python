#!/bin/sh -l
# FILENAME:  Conte_RealData_Submit.sh

#PBS -V
#PBS -q standby
#PBS -l nodes=32:ppn=16
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
uniq < $PBS_NODEFILE > nodefile

#declare -a sigma_s=(1000000 2000000)
#declare -a sigma_t=(8000 16000)
#s_idx=`expr $PARAM_INDEX / 2`
#t_idx=`expr $PARAM_INDEX % 2`

python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 32 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_041024/APS14_AlCu_21.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_041024/APS14_AlCu_19.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_032014_AlCu/ --rot_center 256 --vox_size 1.3 --proj_start 13288 --proj_num 3072 --x_width 1600 --recon_x_width 512 --z_start 0 --z_width 1024 --recon_z_width 256 --vox_stop_thresh 0.5 --cost_stop_thresh 1 --sigma_s 500000 --sigma_t 20000 --K 32 --N_theta 1536 --r 32 --multres_xy 4 --multres_z 2 --do_VarEstimate 1 --MaxIter 1000 --msg_string _proj_start_13288_FastCode --min_time_btw_views 0.0047 --rotation_speed 720 --ZingerT 4 --maxHU 10000 --minHU 0           

#python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 1 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_041024/APS14_AlCu_60.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_041024/APS14_AlCu_60.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_032014_AlCu/ --rot_center 529.95 --vox_size 0.65 --proj_start 118000 --proj_num 3072 --x_width 2080 --recon_x_width 1040 --z_start 300 --z_width 16 --recon_z_width 4 --vox_stop_thresh 0.5 --cost_stop_thresh 1 --sigma_s 500000 --sigma_t 20000 --K 16 --N_theta 3072 --r 16 --multres_xy 4 --multres_z 1 --do_VarEstimate 1 --MaxIter 1000 --msg_string _Dataset_60_proj_start_118000 --min_time_btw_views 0.00868 --rotation_speed 108 --ZingerT 4 --maxHU 10000 --minHU 0           

checkjob -v $PBS_JOBID

#kill $tpid
