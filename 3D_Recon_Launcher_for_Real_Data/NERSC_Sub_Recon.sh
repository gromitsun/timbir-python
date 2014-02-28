#!/bin/tcsh

#PBS -q regular
#PBS -l mppwidth=48
#PBS -l walltime=12:00:00
#PBS -N Ceramic_Data
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR
setenv OMP_NUM_THREADS 32
setenv CRAY_ROOTFS DSL

module load PrgEnv-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py

python XT_MBIR_3D.py --setup_launch_folder --run_reconstruction --Edison --input_hdf5 $GSCRATCH/LBNL_Datasets/20130719_222857_HN2_2011-RT-20N-scan1_x0y0.h5 --group_hdf5 /20130719_222857_HN2_2011-RT-20N-scan1_x0y0 --code_launch_folder $GSCRATCH/Recon_Runs/LBNL_Recons/XT_run/ --output_hdf5 $GSCRATCH/Results/Ceramic/ --pix_size 0.65 --num_views 1024 --x_width 2560 --recon_x_width 2560 --z_start 1000 --z_numElts 32 --num_nodes 1 --num_threads 32 --rot_center 1333 --smoothness 0.50 --view_subsmpl_fact 4 --num_bright 30 --num_dark 30 --zinger_thresh 3.5 --Variance_Est 1 --stop_threshold 10

