#!/bin/sh -l
# FILENAME:  Conte_SimData_Submit.sh

#PBS -V
#PBS -q bouman
#PBS -l nodes=1:ppn=16
#PBS -l walltime=8:00:00

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

#declare -a sigma_s=(50000 100000 250000 500000)
#declare -a sigma_t=(125 250 500 1000 2000)
declare -a sigma_s=(100000 200000 400000)
declare -a sigma_t=(1000, 2000, 4000, 8000)

s_idx=`expr $PARAM_INDEX / 4`
t_idx=`expr $PARAM_INDEX % 4`

#python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --SIM_DATA --num_nodes 1 --Purdue --Path2Data /blah/ --Path2WhitesDarks /blah/ --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/CH_Phantom_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/CH_Mask_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --run_folder $RCAC_SCRATCH/Recon_Runs/Att_Recon_Sim_256/ --rot_center 132 --vox_size 0.65 --proj_start 0 --proj_num 1024 --x_width 1024 --recon_x_width 256 --z_start 0 --z_width 4 --recon_z_width 4 --vox_stop_thresh 0.1 --cost_stop_thresh 0.1 --sigma_s ${sigma_s[$s_idx]} --sigma_t ${sigma_t[$t_idx]} --K 16 --N_theta 256 --r 16 --multres_xy 3 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string _param_sweep --min_time_btw_views 0 --rotation_speed 100 --ZingerT 4 --maxHU 40000 --minHU 0 --BH_Quad_Coef 0.0 --converged_object_file $RCAC_SCRATCH/Recon_Runs/Att_Recon_Sim_256/XT_Results/MBIR_sigs_1000000.0_sigt_16000.0_r_16_K_16_N_theta_256_N_p_1024_converged/object.hdf5

#python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --SIM_DATA --num_nodes 1 --Purdue --Path2Data /blah/ --Path2WhitesDarks /blah/ --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/CH_Phantom_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/CH_Mask_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --run_folder $RCAC_SCRATCH/Recon_Runs/Att_Recon_Sim_256/ --rot_center 132 --vox_size 0.65 --proj_start 0 --proj_num 1024 --x_width 1024 --recon_x_width 256 --z_start 0 --z_width 4 --recon_z_width 4 --vox_stop_thresh 0.1 --cost_stop_thresh 0.1 --sigma_s ${sigma_s[$s_idx]} --sigma_t ${sigma_t[$t_idx]} --K 1 --N_theta 256 --r 1 --multres_xy 3 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string _param_sweep --min_time_btw_views 0 --rotation_speed 100 --ZingerT 4 --maxHU 40000 --minHU 0 --BH_Quad_Coef 0.0 --converged_object_file $RCAC_SCRATCH/Recon_Runs/Att_Recon_Sim_256/XT_Results/MBIR_sigs_1000000.0_sigt_16000.0_r_16_K_16_N_theta_256_N_p_1024_converged/object.hdf5

python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --SIM_DATA --num_nodes 1 --Purdue --Path2Data /blah/ --Path2WhitesDarks /blah/ --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/CH_Phantom_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/CH_Mask_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --run_folder $RCAC_SCRATCH/Recon_Runs/Att_Recon_Sim_256/ --rot_center 132 --vox_size 0.65 --proj_start 0 --proj_num 1024 --x_width 1024 --recon_x_width 256 --z_start 0 --z_width 4 --recon_z_width 4 --vox_stop_thresh 0.1 --cost_stop_thresh 0.1 --sigma_s ${sigma_s[$s_idx]} --sigma_t ${sigma_t[$t_idx]} --K 1 --N_theta 16 --r 1 --multres_xy 3 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string _param_sweep --min_time_btw_views 0 --rotation_speed 100 --ZingerT 4 --maxHU 40000 --minHU 0 --BH_Quad_Coef 0.0 --converged_object_file $RCAC_SCRATCH/Recon_Runs/Att_Recon_Sim_256/XT_Results/MBIR_sigs_1000000.0_sigt_16000.0_r_16_K_16_N_theta_256_N_p_1024_converged/object.hdf5

#------- Usman's Phantom Recon -------#
#python XT_Main.py --run_setup --run_recon --gen_outfiles --MBIR --ATT --SIM_DATA --num_nodes 1 --Purdue --Path2Data /blah/ --Path2WhitesDarks /blah/ --Path2Phantom $RCAC_SCRATCH/Sim_Datasets/CH_Phantom_Usman_2.bin --Path2Mask $RCAC_SCRATCH/Sim_Datasets/CH_Mask_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --run_folder $RCAC_SCRATCH/Recon_Runs/Att_Recon_Sim_256/ --rot_center 128 --vox_size 0.65 --proj_start 0 --proj_num 420 --x_width 256 --recon_x_width 256 --z_start 0 --z_width 4 --recon_z_width 4 --vox_stop_thresh 0.05 --cost_stop_thresh 0.1 --sigma_s ${sigma_s[$s_idx]} --sigma_t ${sigma_t[$t_idx]} --K 1 --N_theta 30 --r 1 --multres_xy 3 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string _Usman --min_time_btw_views 0 --rotation_speed 100 --ZingerT 4 --maxHU 40000 --minHU 0 --BH_Quad_Coef 0.0 --converged_object_file /dummy_path/ --phantom_xy_width 256 --phantom_z_width 4 --proj_start_4_RMSE 30 --proj_num_4_RMSE 360

checkjob -v $PBS_JOBID

#kill $tpid
