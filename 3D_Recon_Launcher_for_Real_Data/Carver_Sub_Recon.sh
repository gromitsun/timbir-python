#!/bin/tcsh

#PBS -l nodes=2:ppn=1
#PBS -l pvmem=10GB
#PBS -l walltime=00:30:00
#PBS -N mpi_test
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

set mpi_tasks=2
cd $PBS_O_WORKDIR

module load python/2.7.3
module load h5py
module load pil
module load mpi4py

setenv OMP_NUM_THREADS 8

#python XT_Main.py --run_reconstruction --node_num $mpi_tasks --Carver
python XT_MBIR_3D.py --setup_launch_folders --run_reconstruction --Carver --input_hdf5 $SCRATCH/LBNL_Datasets/20131004_013841_parikh_soil_microaggregate_1-2_0-15.h5 --output_hdf5 $SCRATCH/Recon_Runs/LBNL_Recons/XT_run/ --pix_size 0.65 --num_views 1024 --x_width 1024 --z_start 1000 --z_numElts 64 --num_nodes 2 --num_threads 8 --rot_center 512 --final_res_multiple 1 --smoothness .25 --zinger_thresh 500 --stop_threshold 35
