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

def main():
	start_time = time.time()

	parser = argparse.ArgumentParser()
	parser.add_argument("--create_objectHDFtiff", help="specify whether you want to create HDF and tiff output files", action="store_true")
	parser.add_argument("--setup_launch_folder", help="Specify whether you want to setup the launch folder", action="store_true")
	parser.add_argument("--run_reconstruction", help="Run reconstruction code", action="store_true")
	
	parser.add_argument("--MBIR_ATT_REAL", help="Recon of attenuation modality from real data", action="store_true")
	parser.add_argument("--MBIR_PHCON_REAL", help="Recon of phase constrast modality from real data", action="store_true")
	parser.add_argument("--MBIR_ATT_SIM", help="Recon of attenuation modality from simulated phantom", action="store_true")

	parser.add_argument("--Edison", help="Use Edison when running on Edison", action="store_true")
	parser.add_argument("--Hopper", help="Use Hopper when running on Hopper", action="store_true")
	parser.add_argument("--Purdue", help="Use Purdue when running on Conte or Carter", action="store_true")
	parser.add_argument("--PC", help="Use PC when running on PC", action="store_true")
	parser.add_argument("-n", "--num_MPI_process", type=int, help="Specifies number of nodes")
	args = parser.parse_args()

	proj = {}
	recon = {}
	files = {}
	recon['node_num'] = args.num_MPI_process
	recon, files = CompSys_Init(args, recon, files)
	
	if (args.MBIR_ATT_SIM):
		recon['modality'] = 'ATT'
		proj, recon, files = attenuation_tomo_sim_init (proj, recon, files)
	elif (args.MBIR_ATT_REAL):
		recon['modality'] = 'ATT'
		proj, recon, files = attenuation_tomo_real_init (proj, recon, files)
	elif (args.MBIR_PHCON_REAL):
		recon['modality'] = 'PHCON'
		proj, recon, files = phase_contrast_tomo_real_init (proj, recon, files)

	proj = proj_init(proj, files)
	recon = recon_init(proj, recon)
	files = files_init(files)
	
	recon['set_up_launch_folder'] = 0	
	if args.setup_launch_folder:
		recon['set_up_launch_folder'] = 1

	recon['reconstruct'] = 0
	if (args.run_reconstruction):
		recon['reconstruct'] = 1

	if (recon['recon_type'] == 'MBIR'):
		print 'main: Will do MBIR reconstruction'
		do_MBIR_reconstruction(proj, recon, files)
	elif (recon['recon_type'] == 'FBP'):
		print 'main: Will do FBP reconstruction'
		do_FBP_reconstruction(proj, recon, files)
	else:
		print 'ERROR: main: Reconstruction type not recognized'

	if (args.MBIR_ATT_SIM == 1):
		compute_RMSE_of_recon (proj, recon, files)
			
	if (args.create_objectHDFtiff):
		if (recon['HPC'] == 'Purdue'):
			recon['size'] = MPI.COMM_WORLD.size
			writepar_object2HDF (proj, recon, files)
			writepar_tiff_from_object_bin_file (proj, recon, files)
		else:
		#	write_tiff_from_object_bin_file (proj, recon, files)
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
