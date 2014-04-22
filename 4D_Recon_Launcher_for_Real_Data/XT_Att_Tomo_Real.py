
import numpy as np

def attenuation_tomo_real_init (proj, recon, files):
	#Dataset paths (Typically, stored in scratch pointed to by the environment variable $RCAC_SCRATCH)
	#Path2Dataset - File path to the HDF file containing the dataset
	#Path2WhiteDark - File path to the HDF file containing the white and dark images
	proj['Path2Dataset'] = files['data_scratch'] + "/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_21.hdf"
	proj['Path2WhiteDark'] = files['data_scratch'] + "/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf"
#	proj['Path2Dataset'] = "../../Datasets/K_16_N_theta_2000_RotSpeed_100_Exp_4_ROI_1000x2080_Ramp_2/k-16-4ms-last_22.hdf"
#	proj['Path2WhiteDark'] = "../../Datasets/K_16_N_theta_2000_RotSpeed_100_Exp_4_ROI_1000x2080_Ramp_2/k-16-4ms-last_31.hdf"
	#voxel_size is the side length of each voxel (in micrometer(um))
	files['Result_Folder'] = files['scratch'] + "/Recon_Runs/APS_032014_AlCu/XT_Result_Repository/" #results will be stored here
	files['Launch_Folder'] = files['scratch'] + "/Recon_Runs/APS_032014_AlCu/XT_run/" #location where reconstruction will be run

	proj['use_slice_white'] = -1 
	proj['recon_N_r'] = 1024/2 # recon_N_r is detector resolution along r-axis used in reconstruction. Subsampled from the actual detector resolution of N_r)
	proj['N_t'] = 1024 # N_t is number of detector slices used in recon
	proj['recon_N_t'] = 1024/2 # Subsamples N_t to recon_N_t before doing reconstruction
	proj['K'] = 32 #Number of interlaced sub-frames in a frame
        recon['r'] = 32 #Number of reconstruction per frame
        recon['msg_string'] = "_ROI_512x512x512_proj_start_16360"
        
	proj['voxel_size'] = 1.3
	proj['slice_t_start'] = 0 # slice_t_start is first detector slice used in recon
	proj['rotation_center_r'] = 512.0/2 # Same units as recon_N_r
	proj['proj_start'] = 16360 # Index of the first view used for reconstruction 
	proj['proj_num'] = 1536*2 # Total number of views used for reconstruction
	proj['N_p'] = 1536*25 # N_p is just used for angle generation in python. Should be greater than proj_start + proj_num
	proj['N_theta'] = 1536 #number of views in a frame
	proj['N_r'] = 1600 #detector resolution along r-axis
	proj['min_time_btw_views'] = 0.0047 #a view will be deleted if the time between two consecutive views is less than min_time_btw_views
	proj['rotation_speed'] = 720 #rotation speed in degrees per second
	
        #recon['sigma_s'] = [100*(10**5)] #spatial regularization. lesser the value smoother the reconstruction
	#recon['sigma_t'] = [(10**4)] #temporal regularization. lesser the value smoother the reconstruction
        recon['sigma_s'] = [5*(10**5)] #spatial regularization. lesser the value smoother the reconstructioni
	recon['sigma_t'] = [2*(10**4)] #temporal regularization. lesser the value smoother the reconstruction
	recon['c_s'] = 10**-6 #parameter of spatial qGGMRF prior controlling the tradeoff between gaussian and generalized gaussian models
	recon['c_t'] = 10**-6 #same but for temporal qGGMRF term
	recon['ZingerT'] = 3 #threshold on error sinogram value above which measurement is classified as zinger
	recon['ZingerDel'] = 0.1 #generalized huber function parameter
	recon['maxHU'] = 10000 #maximum value of reconstructed attenuation coefficient 
	recon['minHU'] = 0 #minimum value of reconstructed attenuation coefficient
	recon['voxel_thresh'] = [0.5, 0.5, 0.5, 0.5] #convergence threshold (percentage change in average magnitude of updates
        recon['cost_thresh'] = [10, 10, 10, 10] #convergence threshold on cost (percentage of change in cost normalized with the change in 1st iteration) 
        recon['delta_xy'] = [8, 4, 2, 1] #voxel size as a multiple of detector pixel size in x-y plane 
        recon['delta_z'] = [4, 4, 2, 1] #voxel size as a multiple of detector pixel size in z plane
        recon['initICD'] = [0, 2, 3, 3]
	# 0 - initialize all voxels to 0
	# 1 - initialize at same resolution as previous stage in multiresolution
	# 2 - initialize by upsampling by a factor of 2 in x-y plane (2D multiresolution)
	# 3 - initialize by upsampling by a factor of 2 in x-y-z (3D multiresolution)
        recon['WritePerIter'] = [0]*len(recon['delta_xy']) #1 - write tiff and binary files of object every iteratio. 0 - Don't write
	recon['WritePerIter'][-1] = 1
        recon['updateProjOffset'] = [0, 2, 3, 3] 
#        recon['updateProjOffset'] = [0, 0, 0, 0] 
	# 0 - Don't update the additive projection offset error term, d.
	# 1 - Initialize d from binary file but don't update
	# 2 - Update d after initializing it with 0
	# 3 - Update d after initializing it from binary file
	recon['readSino4mHDF'] = [0]*len(recon['delta_xy']) # 1 - read projection data from HDF file in C. 0 - read from python if required
	recon['readSino4mHDF'][0] = 1
	recon['iterations'] = [500]*len(recon['delta_xy']) #number of iterations
	recon['do_VarEstimate'] = [1]*len(recon['delta_xy']) #1 - estimate variance term parameter
	recon['Estimate_of_Var'] = 1; #initial estimate for variance parameter	

	recon['recon_type'] = 'MBIR' #MBIR - Does MBIR reconstruction. FBP - Does FBP reconstruction
        recon['sinobin'] = 1 #1 - read projection data from binary file (use this in this configuration. don't change.)
        recon['initMagUpMap'] = [1]*len(recon['delta_xy']) 
	recon['initMagUpMap'][0] = 0 #1 - initialize magnitude update map from binary file. 0 - Don't initialize.
        recon['only_Edge_Updates'] = [0]*len(recon['delta_xy']) # Keep all at 0
        recon['writeTiff'] = [1]*len(recon['delta_xy']) #1 - write reconstructed object to tiff files. 0 - Don't write.
	recon['BH_Quad_Coef'] = 0.0; #0.5 is the quadratic term coefficient used in the beam hardening correction

	# The variables below will be ignored by C-code in this mode
	proj['Expected_Counts'] = 29473 
	proj['phantom_N_xy'] = 2048
	# phantom_N_z is the resolution of phantom along z
	proj['phantom_N_z'] = 4
	proj['Path2Phantom'] = files['data_scratch'] + "/Sim_Datasets/blah.bin"

	return proj, recon, files
