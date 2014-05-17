import argparse
from XT_IOMisc import error_by_flag

def ArgParser ():
	parser = argparse.ArgumentParser()
	parser.add_argument("--gen_outfiles", help="specify whether you want to create HDF and tiff output files", action="store_true")
	parser.add_argument("--run_setup", help="Specify whether you want to setup the launch folder", action="store_true")
	parser.add_argument("--run_recon", help="Run reconstruction code", action="store_true")
	parser.add_argument("--restart", help="Restart the reconstruction", action="store_true")
	
	parser.add_argument("--ATT", help="Imaging modality is attenuation contrast tomography", action="store_true")
	parser.add_argument("--PHCON", help="Imaging modality is phase constrast tomography", action="store_true")
	parser.add_argument("--REAL_DATA", help="Recon of real data", action="store_true")
	parser.add_argument("--SIM_DATA", help="Recon of simulated data", action="store_true")
	parser.add_argument("--MBIR", help="Do MBIR reconstruction", action="store_true")
	parser.add_argument("--FBP", help="Do FBP reconstruction", action="store_true")
	
	parser.add_argument("--Edison", help="Use Edison when running on Edison", action="store_true")
	parser.add_argument("--Hopper", help="Use Hopper when running on Hopper", action="store_true")
	parser.add_argument("--Purdue", help="Use Purdue when running on Conte or Carter", action="store_true")
	parser.add_argument("--PC", help="Use PC when running on PC", action="store_true")
	parser.add_argument("-n", "--num_nodes", type=int, help="Specifies number of nodes")

        parser.add_argument("--Path2Data", help="Full path of the input hdf5 file containing the data")
        parser.add_argument("--Path2WhitesDarks", help="Full path of the input hdf5 file containing the white field and dark field images")
        parser.add_argument("--Path2Phantom", help="Full path of the input hdf5 file containing the data")
        parser.add_argument("--Path2Mask", help="Full path of the input hdf5 file containing the white field and dark field images")
        parser.add_argument("--run_folder", help="Full path where the code executable is going to reside")
        parser.add_argument("--rot_center", help="Center of rotation of the object in units of detector pixels",type=float,default=-1)
        parser.add_argument("--vox_size", help="Size of the detector pixels in micro meter",type=float)
	parser.add_argument("--proj_num", help="Number of projections used for reconstruction",type=int)
	parser.add_argument("--proj_start", help="first projection used for reconstruction",type=int)
        parser.add_argument("--x_width", help="Number of detector elements along x-direction",type=int)
        parser.add_argument("--recon_x_width", help="Number of detector elements along x-direction",type=int,default=-1)
        parser.add_argument("--z_start", help="Starting detector pixel along z-direction",type=int)
        parser.add_argument("--z_width", help="Number of detector pixels to use along z-direction",type=int)
        parser.add_argument("--recon_z_width", help="Number of detector pixels to use along z-direction",type=int)
        parser.add_argument("--vox_stop_thresh",help="Stopping thresold as a percentage of average change in pixel values in percentage",type=float,default=1) 
        parser.add_argument("--cost_stop_thresh",help="Stopping thresold on cost",type=float,default=10) 
        parser.add_argument("--sigma_s",help="spatial regularization parameter",type=float) 
        parser.add_argument("--sigma_t",help="temporal regularization parameter",type=float) 
        
	parser.add_argument("--K",help="K is the number of sub-frames",type=int)
        parser.add_argument("--N_theta",help="N_theta is the number of frames",type=int)       
        parser.add_argument("--r",help="r is the number of reconstructions per frame",type=int)       
        parser.add_argument("--multres_xy",help="The number of multiresolution stages in x-y plane",type=int)       
        parser.add_argument("--multres_z",help="The number of multiresolution stages along z-axis",type=int)       
        parser.add_argument("--do_VarEstimate",help="Do variance estimate",type=int)       
        parser.add_argument("--MaxIter",help="Maximum number of iterations",type=int)       
        parser.add_argument("--msg_string", help="message string to be appended to folder name at the end")
        parser.add_argument("--min_time_btw_views",help="minimum time between views",type=float,default=0) 
        parser.add_argument("--rotation_speed",help="rotation speed in deg/s",type=float,default=100) 
        parser.add_argument("--ZingerT",help="Zinger threshold beyond which error sinogram measurement is classified as anamalous",type=float,default=4.0) 
        parser.add_argument("--maxHU",help="Maximum value of attenuation",type=float,default=40000) 
        parser.add_argument("--minHU",help="Minimum value of attenuation",type=float,default=0.0) 
        parser.add_argument("--BH_Quad_Coef",help="Beamhardening quadratic coefficient",type=float,default=0.0) 
	args = parser.parse_args()
	
	if (args.multres_xy < args.multres_z):
		error_by_flag(1, "Number of multi-resolution stages in x-y should be greater than along z axis")
	if (args.multres_xy < 2):
		error_by_flag(1, "Must have atleast two multiresolution stages")
	return args