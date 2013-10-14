#!/bin/tcsh

#PBS -q regular
#PBS -l mppwidth=48
#PBS -l walltime=5:00:00
#PBS -N mpi_test
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

set mpi_tasks=2
cd $PBS_O_WORKDIR
setenv OMP_NUM_THREADS 32
setenv CRAY_ROOTFS DSL

module load PrgEnv-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py

# for Intel compiled programs
# the "-cc numa_node" option should be used if the number of threads is less than or equal 24 with Hyperthreading
# (note: use "-cc none" instead for other number of threads)
# aprun -j 2 -n 128 -N 2 -S 1 -d 24 -cc numa_node ./my_executable      
# for Cray or GNU compiled programs
#python XT_Main.py --run_reconstruction --node_num $mpi_tasks --NERSC
python XT_MBIR_3D.py --setup_launch_folders --run_reconstruction --NERSC --input_hdf5 /global/scratch2/sd/svenkata/LBNL_Datasets/20131004_013841_parikh_soil_microaggregate_1-2_0-15.h5 --output_hdf5 /scratch1/scratchdirs/svenkata/Recon_Runs/LBNL_Recons/XT_run/ --pix_size 0.65 --num_views 1024 --x_width 1024 --z_start 1000 --z_numElts 64 --num_nodes 2 --num_threads 32 --rot_center 512 --final_res_multiple 1 --smoothness .25 --zinger_thresh 500 --stop_threshold 35
