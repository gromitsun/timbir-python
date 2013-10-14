#!/bin/tcsh

#PBS -l nodes=2:ppn=1
#PBS -l pvmem=10GB
#PBS -l walltime=5:00:00
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

python XT_Main.py --run_reconstruction --node_num $mpi_tasks --Carver
