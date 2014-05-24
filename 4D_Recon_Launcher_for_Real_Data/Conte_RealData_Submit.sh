#!/bin/sh -l
# FILENAME:  Conte_RealData_Submit.sh

#PBS -V
#PBS -q standby
#PBS -l nodes=70:ppn=16
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

declare -a sigma_s=(100000 200000 400000)
declare -a sigma_t=(10000)
s_idx=`expr $PARAM_INDEX / 1`
t_idx=`expr $PARAM_INDEX % 1`

#python XT_Main.py --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 70 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_14.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_032014_AlCu/ --rot_center 256 --vox_size 1.3 --proj_start 31856 --proj_num 6144 --x_width 1600 --recon_x_width 512 --z_start 0 --z_width 1050 --recon_z_width 350 --vox_stop_thresh 0.5 --cost_stop_thresh 1 --sigma_s 500000 --sigma_t 20000 --K 32 --N_theta 1536 --r 32 --multres_xy 4 --multres_z 1 --do_VarEstimate 1 --MaxIter 1000 --msg_string _dataset_14 --min_time_btw_views 0.0047 --rotation_speed 720 --ZingerT 4 --maxHU 10000 --minHU 0 --BH_Quad_Coef 0.5           

python XT_Main.py --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 70 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_17.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_032014_AlCu/ --rot_center 256 --vox_size 1.3 --proj_start 47500 --proj_num 6144 --x_width 1600 --recon_x_width 512 --z_start 0 --z_width 1050 --recon_z_width 350 --vox_stop_thresh 0.5 --cost_stop_thresh 1 --sigma_s 500000 --sigma_t 20000 --K 32 --N_theta 1536 --r 32 --multres_xy 4 --multres_z 1 --do_VarEstimate 1 --MaxIter 1000 --msg_string _dataset_17_proj_start_47500 --min_time_btw_views 0.0047 --rotation_speed 360 --ZingerT 4 --maxHU 50000 --minHU 0 --BH_Quad_Coef 0.5           

python XT_Main.py --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 70 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_17.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_032014_AlCu/ --rot_center 256 --vox_size 1.3 --proj_start 47500 --proj_num 6144 --x_width 1600 --recon_x_width 512 --z_start 0 --z_width 1050 --recon_z_width 350 --vox_stop_thresh 0.5 --cost_stop_thresh 1 --sigma_s 500000 --sigma_t 20000 --K 32 --N_theta 1536 --r 32 --multres_xy 4 --multres_z 1 --do_VarEstimate 1 --MaxIter 1000 --msg_string _dataset_17_proj_start_47500 --min_time_btw_views 0.0047 --rotation_speed 360 --ZingerT 4 --maxHU 50000 --minHU 0 --BH_Quad_Coef 0.5           

#### DATASET 21 #####
#python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 3 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_21.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_032014_AlCu/ --rot_center 256 --vox_size 1.3 --proj_start 13288 --proj_num 3072 --x_width 1600 --recon_x_width 512 --z_start 512 --z_width 24 --recon_z_width 12 --vox_stop_thresh 0.5 --cost_stop_thresh 1 --sigma_s 500000 --sigma_t 20000 --K 32 --N_theta 1536 --r 32 --multres_xy 4 --multres_z 2 --do_VarEstimate 1 --MaxIter 1000 --msg_string _proj_start_13288_OldCode --min_time_btw_views 0.0047 --rotation_speed 720 --ZingerT 4 --maxHU 10000 --minHU 0 --BH_Quad_Coef 0.5           

#python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 32 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_21.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_032014_AlCu/ --rot_center 800 --vox_size 1.3 --proj_start 10000 --proj_num 10752 --x_width 1600 --recon_x_width 1600 --z_start 0 --z_width 128 --recon_z_width 128 --vox_stop_thresh 0.5 --cost_stop_thresh 1 --sigma_s 500000 --sigma_t 5000 --K 32 --N_theta 1536 --r 32 --multres_xy 4 --multres_z 1 --do_VarEstimate 1 --MaxIter 1000 --msg_string _z_0 --min_time_btw_views 0.0047 --rotation_speed 720 --ZingerT 4 --maxHU 10000 --minHU 0 --BH_Quad_Coef 0.5 

#### DATASET 60 ####
#python XT_Main.py --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 1 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_60_AlSi.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_60_AlSi.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_Dataset_60/ --rot_center 529.95 --vox_size 0.65 --proj_start 150000 --proj_num 3072 --x_width 2080 --recon_x_width 1040 --z_start 300 --z_width 8 --recon_z_width 4 --vox_stop_thresh 1 --cost_stop_thresh 1 --sigma_s ${sigma_s[$s_idx]} --sigma_t ${sigma_t[$t_idx]} --K 16 --N_theta 3072 --r 16 --multres_xy 4 --multres_z 1 --do_VarEstimate 1 --MaxIter 1000 --msg_string _sweep --min_time_btw_views 0.00868 --rotation_speed 108 --ZingerT 4 --maxHU 10000 --minHU 0 --BH_Quad_Coef 0.0

#python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 1 --Purdue --Path2Data $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_60_AlSi.hdf --Path2WhitesDarks $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_60_AlSi.hdf --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/APS_Dataset_60/ --rot_center 529.95 --vox_size 0.65 --proj_start 150000 --proj_num 3072 --x_width 2080 --recon_x_width 1040 --z_start 300 --z_width 8 --recon_z_width 4 --vox_stop_thresh 1 --cost_stop_thresh 1 --sigma_s ${sigma_s[$s_idx]} --sigma_t ${sigma_t[$t_idx]} --K 16 --N_theta 3072 --r 16 --multres_xy 4 --multres_z 1 --do_VarEstimate 1 --MaxIter 1000 --msg_string _sweep --min_time_btw_views 0.00868 --rotation_speed 108 --ZingerT 4 --maxHU 10000 --minHU 0 --BH_Quad_Coef 0.0

### SLS Dataset ###
#python XT_Main.py --gen_outfiles --MBIR --PHCON --REAL_DATA --num_nodes 1 --Purdue --Path2Data $RCAC_SCRATCH/SLS_Datasets/SLS_Dataset.h5 --Path2WhitesDarks $RCAC_SCRATCH/SLS_Datasets/SLS_Dataset.h5 --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/PhCon_Real/ --rot_center 636 --vox_size 0.74 --proj_start 0 --proj_num 1000 --x_width 2560 --recon_x_width 1280 --z_start 0 --z_width 8 --recon_z_width 4 --vox_stop_thresh 2 --cost_stop_thresh 10 --sigma_s 100000 --sigma_t 100000 --K 1 --N_theta 1000 --r 1 --multres_xy 6 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string _N_r_1280 --min_time_btw_views 0.0 --rotation_speed 108 --ZingerT 4 --maxHU 10000 --minHU 0 --BH_Quad_Coef 0.0

#python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --PHCON --REAL_DATA --num_nodes 1 --Purdue --Path2Data $RCAC_SCRATCH/SLS_Datasets/SLS_Dataset_sub4.h5 --Path2WhitesDarks $RCAC_SCRATCH/SLS_Datasets/SLS_Dataset_sub4.h5 --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $RCAC_SCRATCH/Recon_Runs/PhCon_Real/ --rot_center 636 --vox_size 0.74 --proj_start 0 --proj_num 250 --x_width 2560 --recon_x_width 1280 --z_start 0 --z_width 8 --recon_z_width 4 --vox_stop_thresh 0.5 --cost_stop_thresh 10  --sigma_s ${sigma_s[$s_idx]} --sigma_t ${sigma_t[$t_idx]} --K 1 --N_theta 250 --r 1 --multres_xy 6 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string _PhCon --min_time_btw_views 0.0 --rotation_speed 108 --ZingerT 4 --maxHU 10000 --minHU 0 --BH_Quad_Coef 0.0

checkjob -v $PBS_JOBID

#kill $tpid
