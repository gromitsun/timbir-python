#!/bin/tcsh

#PBS -q regular
#PBS -l mppwidth=48
#PBS -l walltime=24:00:00
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
python XT_Main_3D.py --run_reconstruction --node_num $mpi_tasks --NERSC
