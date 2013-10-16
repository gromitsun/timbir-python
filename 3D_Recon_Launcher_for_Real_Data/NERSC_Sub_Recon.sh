#!/bin/tcsh

#PBS -q debug
#PBS -l mppwidth=48
#PBS -l walltime=00:30:00
#PBS -N mpi_test
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

set mpi_tasks=1
cd $PBS_O_WORKDIR
setenv OMP_NUM_THREADS 32
setenv CRAY_ROOTFS DSL

module load PrgEnv-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py

#python XT_MBIR_3D.py --setup_launch_folder --run_reconstruction --NERSC --input_hdf5 $GSCRATCH/LBNL_Datasets/20130118_150718_chevron_orig_again_WLa.h5 --output_hdf5 $GSCRATCH/Recon_Runs/LBNL_Recons/XT_run/20130118_150718_chevron_orig_again_WLa/ --pix_size 1.33 --num_views 2048 --x_width 1024 --z_start 130 --z_numElts 12 --num_nodes 1 --num_threads 32 --rot_center 462 --final_res_multiple 1 --smoothness 1 --zinger_thresh 50 --stop_threshold 35 --view_subsmpl_fact 8 --num_res 8

python XT_MBIR_3D.py --setup_launch_folder --run_reconstruction --NERSC --input_hdf5 $GSCRATCH/LBNL_Datasets/20131004_013841_parikh_soil_microaggregate_1-2_0-15.h5 --code_launch_folder $GSCRATCH/Recon_Runs/LBNL_Recons/XT_run/ --output_hdf5 $GSCRATCH/Results/Edison/ --pix_size 0.65 --num_views 1024 --x_width 2560 --recon_x_width 2560 --z_start 1000 --z_numElts 64 --num_nodes 2 --num_threads 32 --rot_center 1276 --final_res_multiple 1 --smoothness .15 --zinger_thresh 50 --stop_threshold 35 --view_subsmpl_fact 4 --num_res 4
