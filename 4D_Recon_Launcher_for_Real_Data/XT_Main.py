#! /usr/bin/python

from XT_Initialize import proj_init, recon_init, files_init 
#from XT_Projections import generate_projections
from XT_MBIR_Reconstruction import do_MBIR_reconstruction
#from XT_FBP_Reconstruction import do_FBP_reconstruction
from XT_IOMisc import error_by_flag
from XT_IOMisc import write_object2HDF
from XT_ObjectHDFIO import writepar_object2HDF
from XT_ObjectHDFIO import writepar_tiff_from_object_bin_file
from XT_IOMisc import write_tiff_from_object_bin_file
from XT_Interpolate import compute_RMSE_of_recon
import argparse
import time
from XT_Att_Tomo_Sim import attenuation_tomo_sim_init
from XT_Att_Tomo_Real import attenuation_tomo_real_init
from XT_PhCon_Tomo_Real import phase_contrast_tomo_real_init
from XT_CompSys_Init import CompSys_Init
#from mpi4py import MPI
from XT_ArgParser import ArgParser

def main():
	start_time = time.time()

	proj = {}
	recon = {}
	files = {}
	args = ArgParser()
	recon['node_num'] = args.num_nodes
	recon, files = CompSys_Init(args, recon, files)
	
	recon['modality'] = 'ATT'
	if (args.MBIR_PHCON_REAL):
		recon['modality'] = 'PHCON'

	proj = proj_init(proj, args)
	recon = recon_init(proj, recon, args)
	files = files_init(files, args)
	
	recon['set_up_launch_folder'] = 0	
	if args.run_setup:
		recon['set_up_launch_folder'] = 1

	recon['reconstruct'] = 0
	if (args.run_recon):
		recon['reconstruct'] = 1

	if (recon['recon_type'] == 'MBIR'):
		print 'main: Will do MBIR reconstruction'
		do_MBIR_reconstruction(proj, recon, files)
	elif (recon['recon_type'] == 'FBP'):
		do_MBIR_reconstruction(proj, recon, files)
		print 'main: Will do FBP reconstruction'
#		do_FBP_reconstruction(proj, recon, files)
	else:
		print 'ERROR: main: Reconstruction type not recognized'

	if (args.MBIR_ATT_SIM == 1):
		compute_RMSE_of_recon (proj, recon, files)
			
#	if (args.create_objectHDFtiff):
#		if (recon['HPC'] == 'Purdue'):
#			#recon['size'] = MPI.COMM_WORLD.size
#			recon['size'] = 1
#			writepar_object2HDF (proj, recon, files)
#			writepar_tiff_from_object_bin_file (proj, recon, files)
	if (args.gen_outfiles):
		write_tiff_from_object_bin_file (proj, recon, files)
		write_object2HDF (proj, recon, files)
	
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
	write_tiff_from_object_bin_file (proj, recon, files)
	

main()
#bin2tiff_main()
