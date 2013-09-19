#! /usr/bin/python

from XT_Initialize import proj_init, recon_init, files_init 
from XT_Projections import generate_projections
from XT_MBIR_Reconstruction import do_MBIR_reconstruction
#from XT_FBP_Reconstruction import do_FBP_reconstruction
from XT_IOMisc import write_object2HDF
from XT_IOMisc import write_tiff_from_object_bin_file
import time

def main():
	start_time = time.time()

	proj = proj_init()
	recon = recon_init(proj)
	files = files_init()

	proj = generate_projections(proj, recon, files)
	print 'main: Time since start of program - ' + str((time.time() - start_time)/60.0) + 'mins'
	
	if (recon['recon_type'] == 'MBIR'):
		print 'main: Will do MBIR reconstruction'
		do_MBIR_reconstruction(proj, recon, files)	
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

	write_object2HDF(proj,recon,files)
	write_tiff_from_object_bin_file (proj, recon, files)
	

main()
#bin2tiff_main()
