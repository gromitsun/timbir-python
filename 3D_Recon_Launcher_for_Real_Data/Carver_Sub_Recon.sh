#!/bin/tcsh

#PBS -q regular
#PBS -l nodes=4:ppn=1
#PBS -l pvmem=20GB
#PBS -l walltime=04:00:00
#PBS -N mpi_test
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

set mpi_tasks=4
cd $PBS_O_WORKDIR

#module swap pgi intel
#module swap openmpi openmpi-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py

setenv OMP_NUM_THREADS 8

python XT_MBIR_3D.py --setup_launch_folder --run_reconstruction --Carver --input_hdf5 $SCRATCH/LBNL_Datasets/20131004_013841_parikh_soil_microaggregate_1-2_0-15.h5 --output_hdf5 $SCRATCH/Recon_Runs/LBNL_Recons/XT_run/ --pix_size 0.65 --num_views 1024 --x_width 2560 --z_start 1000 --z_numElts 256 --num_nodes 4 --num_threads 8 --rot_center 1276 --final_res_multiple 1 --smoothness .25 --zinger_thresh 40 --stop_threshold 35 --view_subsmpl_fact 2 
