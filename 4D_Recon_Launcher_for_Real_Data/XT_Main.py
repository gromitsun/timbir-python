#! /usr/bin/python

from XT_Initialize import proj_init, recon_init, files_init 
#from XT_Projections import generate_projections
from XT_MBIR_Reconstruction import do_MBIR_reconstruction
#from XT_FBP_Reconstruction import do_FBP_reconstruction
from XT_IOMisc import error_by_flag
from XT_IOMisc import write_object2HDF
from XT_ObjectHDFIO import writepar_object2HDF
from XT_IOMisc import write_tiff_from_object_bin_file
import argparse
import time
import os
from mpi4py import MPI

def main():
	start_time = time.time()

	parser = argparse.ArgumentParser()
	parser.add_argument("--setup_launch_folder", help="Specify whether you want to setup the launch folder", action="store_true")
	parser.add_argument("--run_reconstruction", help="Run reconstruction code", action="store_true")
	parser.add_argument("--NERSC", help="Use NERSC when running on NERSC systems", action="store_true")
	parser.add_argument("--Purdue", help="Use Purdue when running on Conte or Carter", action="store_true")
	parser.add_argument("-n", "--node_num", type=int, help="Specifies number of nodes")
	args = parser.parse_args()

	recon = {}
	files = {}
	recon['node_num'] = args.node_num
	if (args.Purdue):
		recon['num_threads'] = 32
		files['scratch'] = os.environ['RCAC_SCRATCH']
		recon['run_command'] = 'mpiexec -n ' + str(recon['node_num']) + ' -machinefile ' + os.getcwd() + '/nodefile '
		recon['compile_command'] = 'mpicc '
	elif (args.NERSC):
		recon['num_threads'] = 16
		files['scratch'] = os.environ['SCRATCH']
		recon['run_command'] = 'aprun -j 2 -n ' + str(recon['node_num']) + ' -N 2 -S 1 -d ' + str(recon['num_threads'])
		recon['compile_command'] = 'cc '
	else:
		error_by_flag(1, 'HPC system not recognized')	
	
	proj = proj_init(files)
	recon = recon_init(proj, recon)
	files = files_init(files)

	recon['rank'] = MPI.COMM_WORLD.rank	
	recon['set_up_launch_folder'] = 0	
	if args.setup_launch_folder:
		recon['set_up_launch_folder'] = 1	
	
	recon['reconstruct'] = 0
	if (args.run_reconstruction):
		recon['reconstruct'] = 1	

	if (recon['recon_type'] == 'MBIR'):
		print 'main: Will do MBIR reconstruction'
		do_MBIR_reconstruction(proj, recon, files)
		if (args.run_reconstruction):	
			write_tiff_from_object_bin_file (proj, recon, files)
			write_object2HDF (proj, recon, files)
	elif (recon['recon_type'] == 'FBP'):
		print 'main: Will do FBP reconstruction'
		do_FBP_reconstruction(proj, recon, files)	
	else:
		print 'ERROR: main: Reconstruction type not recognized'		

	print 'main: Done!'
	print 'main: Total time taken by program - ' + str((time.time() - start_time)/60.0) + 'mins'


def bin2tiff_main ():
	proj = proj_init()
	recon = recon_init(proj)
	files = files_init()

	parser = argparse.ArgumentParser()
	parser.add_argument("-n", "--node_num", type=int, help="Specifies number of nodes")
	args = parser.parse_args()
	recon['node_num'] = args.node_num
	write_object2HDF(proj,recon,files)
	#write_tiff_from_object_bin_file (proj, recon, files)
	

main()
#bin2tiff_main()
