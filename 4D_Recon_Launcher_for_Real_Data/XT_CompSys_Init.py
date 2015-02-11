

import os


def CompSys_Init(args, recon, files):
	if (args.PC):
		recon['num_threads'] = 8
		files['scratch'] = '../..'
		files['data_scratch'] = '/Users/aditya/Academics/Graduate_Courses/ECE699/Time_Varying_XRay_Tomography/C_Code/Workspace_Argonne'
		recon['run_command'] = 'mpiexec -n ' + str(recon['node_num'])
		recon['compile_command'] = 'mpicc -ansi -Wall -fopenmp '
		#recon['compile_command'] = 'mpicc -fopenmp '
		recon['HPC'] = 'PC' 
		recon['rank'] = 0
	elif (args.Quest):
		recon['num_threads'] = 20
		files['scratch'] = '/tmp'
		files['data_scratch'] = '/tmp'
		if (recon['node_num'] > 1):
			recon['run_command'] = 'mpiexec -n ' + str(recon['node_num']) + ' -machinefile nodefile '
		else:
			recon['run_command'] = 'mpiexec -n ' + str(recon['node_num'])
		recon['compile_command'] = 'mpicc -ansi -Wall -fopenmp -O3 -I/software/hdf5/1.8.12/include -L/software/hdf5/1.8.12/lib -lhdf5 ' # gcc
		#recon['compile_command'] = 'mpiicc -ansi -Wall -openmp -O3 -I/software/hdf5/1.8.10/include -L/software/hdf5/1.8.10/lib -lhdf5 '  # intel
		#recon['compile_command'] = 'mpiicc -ansi -Wall -openmp -O3 -I/software/hdf5/1.8.12/include -L/software/hdf5/1.8.12/lib -lhdf5 -L/software/mpi/openmpi-1.6.3-gcc-4.6.3-trq4/lib64 -lmpi'  # intel
		recon['HPC'] = 'Quest' 
		#recon['rank'] = MPI.COMM_WORLD.rank
		recon['rank'] = 0
	elif (args.Purdue):
		recon['num_threads'] = 32
		files['scratch'] = os.environ['RCAC_SCRATCH']
		files['data_scratch'] = os.environ['RCAC_SCRATCH']
		if (recon['node_num'] > 1):
			recon['run_command'] = 'mpiexec -n ' + str(recon['node_num']) + ' -machinefile nodefile '
		else:
			recon['run_command'] = 'mpiexec -n ' + str(recon['node_num'])
		recon['compile_command'] = 'mpicc -ansi -Wall -openmp -O3 -lhdf5 '
		recon['HPC'] = 'Purdue' 
		#recon['rank'] = MPI.COMM_WORLD.rank
		recon['rank'] = 0
	elif (args.Edison or args.Hopper):
		files['scratch'] = os.environ['SCRATCH']
		files['data_scratch'] = os.environ['GSCRATCH']
		if (args.Edison):
			recon['num_threads'] = 32
			recon['run_command'] = 'aprun -j 2 -n ' + str(recon['node_num']) + ' -N 1 -d ' + str(recon['num_threads']) + ' -cc none '
			recon['compile_command'] = 'cc -Wall -ansi -openmp -O3 '
		else:
			recon['num_threads'] = 16
			#recon['run_command'] = 'aprun -n ' + str(recon['node_num']) + ' -N 1 -d ' + str(recon['num_threads']) + ' -cc none valgrind --leak-check=yes --num-callers=500 --track-origins=yes '
			recon['run_command'] = 'aprun -n ' + str(recon['node_num']) + ' -N 1 -d ' + str(recon['num_threads']) + ' -cc none '
			#recon['compile_command'] = 'cc -g -O0 -openmp '
			recon['compile_command'] = 'cc -fopenmp '
		recon['HPC'] = 'NERSC'
		recon['rank'] = 0
	elif (args.linuxPC):
		recon['num_threads'] = 4
		files['scratch'] = '../..'
		files['data_scratch'] = '/home/y/temp/'
		recon['run_command'] = 'mpiexec -n ' + str(recon['node_num'])
		recon['compile_command'] = 'mpicc -ansi -Wall -fopenmp '
		#recon['compile_command'] = 'mpicc -fopenmp '
		recon['HPC'] = 'linuxPC' 
		recon['rank'] = 0		
	else:
		error_by_flag(1, 'HPC system not recognized')
	return (recon, files)
