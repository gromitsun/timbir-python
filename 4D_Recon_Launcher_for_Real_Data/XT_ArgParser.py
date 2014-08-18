import argparse
from XT_IOMisc import error_by_flag

def ArgParser ():
	parser = argparse.ArgumentParser()
	parser.add_argument("--gen_outfiles", help="Specify whether you want to generate output files (i.e. of the reconstruction) in HDF and Tiff formats", action="store_true")
	parser.add_argument("--run_setup", help="Specify whether you want to setup the launch folder where the reconstruction will be run. You may not want to set up the launch folder if it already exists and you don't want to overwrite it.", action="store_true")
	parser.add_argument("--run_recon", help="Run the reconstruction", action="store_true")
        parser.add_argument("--restart_stage", help="Multi-resolution stage number from where the reconstruction is to be re-started",type=int)
	parser.add_argument("--restart", help="Restart the reconstruction from the mult-resolution stage specified by the input '--restart_stage'", action="store_true")
	parser.add_argument("--same_stage", help="Restart the reconstruction from the last running stage if output files were updated every iteration by setting write_per_iter = 1", action="store_true")
	
	parser.add_argument("--ATT", help="Imaging modality is attenuation contrast tomography", action="store_true")
	parser.add_argument("--PHCON", help="Imaging modality is phase constrast tomography", action="store_true")
	parser.add_argument("--REAL_DATA", help="Reconstruction of real data", action="store_true")
	parser.add_argument("--SIM_DATA", help="Reconstruction of simulated data", action="store_true")
	parser.add_argument("--MBIR", help="Do MBIR reconstruction", action="store_true")
	parser.add_argument("--FBP", help="Do FBP reconstruction", action="store_true")
	
	parser.add_argument("--Edison", help="Cluster used is Edison at NERSC", action="store_true")
	parser.add_argument("--Hopper", help="Cluster used is Hopper at NERSC", action="store_true")
	parser.add_argument("--Purdue", help="Cluster used is Conte or Carter at Purdue", action="store_true")
	parser.add_argument("--PC", help="Use 'PC' when running on a mac", action="store_true")
	parser.add_argument("-n", "--num_nodes", type=int, help="Specifies number of nodes/processors ( = number of MPI processes)")

        parser.add_argument("--Path2Data", help="Full path of the input hdf5 file containing the data", default="/dummy_path/")
        parser.add_argument("--Path2Whites", help="Full path of the input hdf5 file containing the white field images", default="/dummy_path/")
        parser.add_argument("--Path2Darks", help="Full path of the input hdf5 file containing the dark field images", default="/dummy_path/")
        parser.add_argument("--Path2Phantom", help="Path to the binary file containing the 4D phantom", default="/dummy_path/")
        parser.add_argument("--Path2Mask", help="Path to the binary file containing the mask for the 4D phantom over which the RMSE is computed", default="/dummy_path/")
        parser.add_argument("--run_folder", help="Full path to the folder where the reconstruction is run and results are stored. Inside this folder, there should be two other folders called 'XT_run' (inside this folder a new folder is automatically created and the reconstruction is run from within this folder) and 'XT_Results' (inside this folder a new folder is automatically created and the output fildes are stored within this folder)")
	#The automatically created folders have several parameters like sigma_s, sigma_t, etc., in their names. This prevents overwriting of folders during parameter sweeps, etc. 
        parser.add_argument("--rot_center", help="Center of rotation of the object in units of detector pixels (in the same units as the value specified in '--x_width')",type=float)
        parser.add_argument("--vox_size", help="Size of the detector pixels in micro meter",type=float)
	parser.add_argument("--proj_num", help="Number of projections used for reconstruction",type=int)
	parser.add_argument("--proj_start", help="first projection used for reconstruction",type=int)
        parser.add_argument("--x_width", help="Actual number of detector pixels along x-direction",type=int)
        parser.add_argument("--recon_x_width", help="Downsampled number of detector pixels along x-direction",type=int)
	#downsampling is from length 'x_width' to 'recon_x_width'
        parser.add_argument("--z_start", help="Starting detector pixel along z-direction",type=int)
        parser.add_argument("--z_width", help="Number of detector pixels to use along z-direction",type=int)
        parser.add_argument("--recon_z_width", help="Downsampled number of detector pixels to use along z-direction",type=int)
	#downsampling is from length 'z_width' to 'recon_z_width'
        parser.add_argument("--vox_stop_thresh",help="Stopping thresold expressed as a percentage of average magnitude of voxel updates",type=float) 
        parser.add_argument("--cost_stop_thresh",help="Stopping thresold on cost",type=float,default=10) #Not used if cost calculation is switched off 
        parser.add_argument("--sigma_s",help="spatial regularization parameter",type=float) 
        parser.add_argument("--sigma_t",help="temporal regularization parameter",type=float) 
        
	parser.add_argument("--K",help="K is the number of sub-frames",type=int)
        parser.add_argument("--N_theta",help="N_theta is the number of frames",type=int)       
        parser.add_argument("--r",help="r is the number of reconstructions per frame",type=int)       
	parser.add_argument("--num_cycles", help="Number of full frames of data (or cycles)",type=int)
        parser.add_argument("--multres_xy",help="The number of multiresolution stages in x-y plane",type=int)       
        parser.add_argument("--multres_z",help="The number of multiresolution stages along z-axis",type=int)       
        parser.add_argument("--do_VarEstimate",help="Estimate the variance parameter",type=int)       
        parser.add_argument("--MaxIter",help="Maximum number of iterations in each multi-resolution stage",type=int, default=1000)       
        parser.add_argument("--msg_string", help="message string to be appended to folder name at the end")
        parser.add_argument("--min_time_btw_views",help="minimum time between views",type=float,default=0) 
        parser.add_argument("--rotation_speed",help="rotation speed in deg/s",type=float,default=100) 
        parser.add_argument("--ZingerT",help="Zinger threshold T beyond which error sinogram measurement is classified as anamalous",type=float,default=4.0) 
        parser.add_argument("--maxHU",help="Maximum value of attenuation (used to scale the output tiff files)",type=float,default=40000) 
        parser.add_argument("--minHU",help="Minimum value of attenuation (used to scale the output tiff files)",type=float,default=0.0) 
        parser.add_argument("--BH_Quad_Coef",help="Leading coefficient of quadratic polynomial used for beamhardening correction",type=float,default=0.0) 
        parser.add_argument("--RMSE_converged",help="Compute RMSE between the reconstruction and the converged result at the finest resolution", action="store_true") 
        parser.add_argument("--converged_object_file", help="Full path of the input hdf5 file containing the converged object with which RMSE comparison should be done", default="/dummy_path/")
        
	parser.add_argument("--phantom_xy_width", help="Number of voxels of the phantom along x-direction",type=int, default=1024)
        parser.add_argument("--phantom_z_width", help="Number of voxels of the phantom along z-direction",type=int, default=4)
        parser.add_argument("--proj_start_4_RMSE", help="RMSE is computed from the time the specified projection was acquired",type=int, default=256)
        parser.add_argument("--proj_num_4_RMSE", help="RMSE is computed over the time duration required to acquire the specified number of projections",type=int, default=512)
	parser.add_argument("--real_is_double", help="Do all reconstruction using double type as 'real' variable. Default is 'float'", action="store_true")
	parser.add_argument("--no_offset_est", help="If set, will not estimate the offset error, d, used to correct ring artifact", action="store_true")
	parser.add_argument("--no_zero_mean_offset", help="If set, will not enforce zero mean constraint on offset error, d, used to correct ring artifacts", action="store_true")
	args = parser.parse_args()
	
	if (args.multres_xy < args.multres_z):
		error_by_flag(1, "Number of multi-resolution stages in x-y should be greater than along z axis")
	if (args.multres_xy < 2):
		error_by_flag(1, "Must have atleast two multiresolution stages")
	return args
