#!/bin/tcsh

#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=04:00:00
#PBS -N dataset_14


cd $PBS_O_WORKDIR
setenv  OMP_NUM_THREADS DSL
setenv OMP_NUM_THREADS 32
module load PrgEnv-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py
module load cray-hdf5

python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --REAL_DATA --num_nodes 1 --Edison --Path2Data $SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_14.hdf --Path2WhitesDarks $SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf --Path2Phantom $SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard.bin --Path2Mask $SCRATCH/Sim_Datasets/phantom_Cahn_Hilliard_mask.bin --run_folder $SCRATCH/Recon_Runs/APS_032014_AlCu/ --rot_center 256 --vox_size 1.3 --proj_start 31856 --proj_num 6144 --x_width 1600 --recon_x_width 512 --z_start 0 --z_width 12 --recon_z_width 4 --vox_stop_thresh 0.5 --cost_stop_thresh 1 --sigma_s 500000 --sigma_t 20000 --K 32 --N_theta 1536 --r 32 --multres_xy 4 --multres_z 1 --do_VarEstimate 1 --MaxIter 1000 --msg_string _dataset_14 --min_time_btw_views 0.0047 --rotation_speed 720 --ZingerT 4 --maxHU 10000 --minHU 0 --BH_Quad_Coef 0.5           

