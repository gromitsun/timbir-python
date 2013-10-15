#! bin/tcsh 

#PBS -l nodes=2:ppn=1
#PBS -l pvmem=10GB
#PBS -l walltime=00:20:00
#PBS -N setupFolders
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

set mpi_tasks=2 #Total tasks for C-Code 
cd $PBS_O_WORKDIR

#module swap PrgEnv-pgi PrgEnv-intel
#module load PrgEnv-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py

setenv OMP_NUM_THREADS 8
mpirun -np 1 python XT_Main.py --setup_launch_folder --node_num $mpi_tasks --Carver
