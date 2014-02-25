

import os


def CompSys_Init(args, recon, files):
	if (args.PC):
		recon['num_threads'] = 1
		files['scratch'] = '../..'
		files['data_scratch'] = '/Users/aditya/Academics/Graduate_Courses/ECE699/Time_Varying_XRay_Tomography/C_Code/Workspace_Argonne'
		recon['run_command'] = 'mpiexec -n ' + str(recon['node_num'])
		recon['compile_command'] = 'mpicc -ansi -Wall -fopenmp '
		recon['HPC'] = 'PC' 
		recon['rank'] = 0
	elif (args.Purdue):
		recon['num_threads'] = 32
		files['scratch'] = os.environ['RCAC_SCRATCH']
		files['data_scratch'] = os.environ['RCAC_SCRATCH']
		if (recon['node_num'] > 1):
			recon['run_command'] = 'mpiexec -n ' + str(recon['node_num']) + ' -machinefile nodefile '
		else:
			recon['run_command'] = 'mpiexec -n ' + str(recon['node_num'])
		recon['compile_command'] = 'mpicc -ansi -Wall -openmp -lhdf5 '
		recon['HPC'] = 'Purdue' 
		#recon['rank'] = MPI.COMM_WORLD.rank
		recon['rank'] = 0
	elif (args.Edison or args.Hopper):
		files['scratch'] = os.environ['SCRATCH']
		files['data_scratch'] = os.environ['GSCRATCH']
		if (args.Edison):
			recon['num_threads'] = 32
			recon['run_command'] = 'aprun -j 2 -n ' + str(recon['node_num']) + ' -N 1 -d ' + str(recon['num_threads']) + ' -cc none '
			recon['compile_command'] = 'cc -Wall -ansi -openmp '
		else:
			recon['num_threads'] = 16
			#recon['run_command'] = 'aprun -n ' + str(recon['node_num']) + ' -N 1 -d ' + str(recon['num_threads']) + ' -cc none valgrind --leak-check=yes --num-callers=500 --track-origins=yes '
			recon['run_command'] = 'aprun -n ' + str(recon['node_num']) + ' -N 1 -d ' + str(recon['num_threads']) + ' -cc none '
			#recon['compile_command'] = 'cc -g -O0 -openmp '
			recon['compile_command'] = 'cc -fopenmp '
		recon['HPC'] = 'NERSC'
		recon['rank'] = 0		
	else:
		error_by_flag(1, 'HPC system not recognized')
	return (recon, files)
