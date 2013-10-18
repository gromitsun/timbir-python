#!/bin/tcsh

#PBS -q regular
#PBS -l nodes=1:ppn=1
#PBS -l pvmem=20GB
#PBS -l walltime=02:00:00
#PBS -N mpi_test
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR

#module swap pgi intel
#module swap openmpi openmpi-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py

setenv OMP_NUM_THREADS 8

python XT_MBIR_3D.py --setup_launch_folder --run_reconstruction --Carver --input_hdf5 $GSCRATCH/LBNL_Datasets/20131004_013841_parikh_soil_microaggregate_1-2_0-15.h5 --code_launch_folder $GSCRATCH/Recon_Runs/LBNL_Recons/XT_run/ --output_hdf5 $GSCRATCH/Results/Carver/ --pix_size 0.65 --num_views 1024 --x_width 2560 --recon_x_width 1024 --z_start 1000 --z_numElts 32 --num_nodes 1 --num_threads 8 --rot_center 1276 --smoothness .024 --view_subsmpl_fact 4 --num_bright_dark 30
