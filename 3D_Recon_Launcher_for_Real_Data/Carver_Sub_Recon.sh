#!/bin/tcsh

#PBS -q debug
#PBS -l nodes=1:ppn=1
#PBS -l pvmem=20GB
#PBS -l walltime=00:30:00

#PBS -N mpi_test
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR

module swap pgi intel
module swap openmpi openmpi-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py

setenv OMP_NUM_THREADS 8

python XT_MBIR_3D.py --setup_launch_folder --run_reconstruction --Carver --input_hdf5 $GSCRATCH/LBNL_Datasets/20130110_235013_km3_peeled1_x0y0.h5 --code_launch_folder $GSCRATCH/Recon_Runs/LBNL_Recons/XT_run/ --output_hdf5 $GSCRATCH/Results/Carver/ --pix_size 3.33 --num_views 2048 --x_width 2560 --recon_x_width 2560 --z_start 200 --z_numElts 3 --num_nodes 1 --num_threads 8 --rot_center 1279 --smoothness 1.1 --view_subsmpl_fact 8 --num_bright_dark 20 --multires_2D --num_res 5


