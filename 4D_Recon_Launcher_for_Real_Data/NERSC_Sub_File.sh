#PBS -q debug
#PBS -l mppwidth=48
#PBS -l walltime=6:00:00
#PBS -N mpi_test
#PBS -j oe 
#PBS -V

mpi_tasks=2
cd $PBS_O_WORKDIR
setenv OMP_NUM_THREADS 16
module load PrgEnv-gnu
module load python
module load h5py

# for Intel compiled programs
# the "-cc numa_node" option should be used if the number of threads is less than or equal 24 with Hyperthreading
# (note: use "-cc none" instead for other number of threads)
# aprun -j 2 -n 128 -N 2 -S 1 -d 24 -cc numa_node ./my_executable      
# for Cray or GNU compiled programs
aprun -j 2 -n $mpi_tasks -N 2 -S 1 -d 16 python XT_Main.py --setup_launch_folder --node_num $mpi_tasks
python XT_Main.py --run_reconstruction --node_num $mpi_tasks
