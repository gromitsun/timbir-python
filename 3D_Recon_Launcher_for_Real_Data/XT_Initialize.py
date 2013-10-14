#! /usr/bin/python

import numpy as np 
from XT_Interlaced_Angles import gen_interlaced_views_0_to_Inf,clip_list_of_views
from math import ceil
from XT_IOMisc import error_by_flag

"""The variable 'proj' is a dictionary datatype containing information about the projections and the detector.
------- Dictionary Key - Explantion -----------------
	Path2dataset - Path to the dataset from which to read the count (or projection) data
	N_r - Actual resolution of the detector along the r - axis which is perpendicular to the axis of rotation
	slice_t - Slice of the actual detector which should be reconstructed
	recon_N_r - Resolution of detector at which reconstruction will be done
	rotation_center_r - Center of rotation along the r - axis assuming recon_N_r is the resolution
	proj_start - First projection to be read from the dataset
	proj_num - Total number of projections to be read from the dataset
	N_p - Total number of projections in the dataset
	K - Number of interlaced sub-frames
	N_theta - Number of projections in a full frame (or total number of angles) 
	length_r - Length of the actual detector along r-axis (in units of micro meter (um))
	L - Number of projections in a sub-frame
	min_time_btw_views - The minimum time between views required by the trigger generation circuit (mostly is 5ms)
	fps_of_camera - Frame rate of the camera
	angles - Angles of projections used for reconstruction
	times - Times at which these projections were acquired """

def proj_init (files):
	proj = {}

	#proj['Path2Dataset'] = files['data_scratch'] + "/LBNL_Datasets/20131004_013841_parikh_soil_microaggregate_1-2_0-15.h5"
	#proj['Path2WhiteDark'] = files['data_scratch'] + "/LBNL_Datasets/20131004_013841_parikh_soil_microaggregate_1-2_0-15.h5"
	#proj['Dataset_Name'] = '20131004_013841_parikh_soil_microaggregate_1-2_0-15'
	
	proj['Path2Dataset'] = files['data_scratch'] + "/LBNL_Datasets/20130118_150718_chevron_orig_again_WLa.h5"
	proj['Path2WhiteDark'] = files['data_scratch'] + "/LBNL_Datasets/20130118_150718_chevron_orig_again_WLa.h5"
	proj['Dataset_Name'] = '20130118_150718_chevron_orig_again_WLa'
	proj['Num_Bright_Dark'] = 30
	#proj['Path2Dataset'] = "/Volumes/Stack-1/APS_Datasets/Solidification_Small_Datasets/K_32_N_theta_1984_RotSpeed_100_Exp_2_ROI_2000x2080_Ramp_5/k-32-02ms_1.hdf"
	#proj['Path2WhiteDark'] = "/Volumes/Stack-1/APS_Datasets/Solidification_Small_Datasets/K_32_N_theta_1984_RotSpeed_100_Exp_2_ROI_2000x2080_Ramp_5/k-32-02ms_1.hdf"
	
	proj['view_subsmpl_fact'] = 4
	proj['recon_N_r'] = 2560/4 #Total number of detector elements to be used (crops to nearest power of 2 and then down samples to specified number 
	proj['slice_t_start'] = 100 #parallel to z
	proj['N_t'] = 4*4 #Number of slices 
	proj['recon_N_t'] = 4 #Downsampled to N_t
	#proj['rotation_center_r'] = 1276.0 #detector pixels from left; To Do  
	proj['rotation_center_r'] = 1179.0/4 #detector pixels from left; To Do  
	proj['proj_start'] = 0 #view index start
	proj['proj_num'] = 2048 #num of views to use 
	proj['N_p'] = 2048 #total number of supposed to be taken. For 3D take equal to proj_num
	proj['K'] = 1 #Set to 1
	proj['N_theta'] = 2048 #Equal to proj num

	proj['N_r'] = 2560 #Total number of detector pixels
	proj['length_r'] = 0.65*proj['N_r'] #0.65 is pixel size in micro meter	
	proj['length_t'] = 0.65*proj['N_t']	
	proj['L'] = proj['N_theta']/proj['K'] 
	
	min_time_btw_views = 0 #default to 0; in units of sec
	rotation_speed = 100 #rotation set to any number for 3d; degrees per sec

	fps_of_camera = proj['L']/(180.0/rotation_speed) 
	angles, times = gen_interlaced_views_0_to_Inf(proj['K'], proj['N_theta'], proj['N_p']) #Generate the list
	angles_clip, times_clip, angles_del, times_del = clip_list_of_views (angles, times, min_time_btw_views, rotation_speed)
	print 'clip_list_of_views: Number of views deleted are ', angles_del.size, '. Deleted views are ' + str(angles_del) 
	print 'clip_list_of_views: Deleted views\' times are ' + str(times_del)

	proj['angles'] = angles_clip[proj['proj_start'] : proj['proj_start'] + proj['proj_num'] : proj['view_subsmpl_fact']]
	proj['times'] = times_clip[proj['proj_start'] : proj['proj_start'] + proj['proj_num'] : proj['view_subsmpl_fact']]
	
	proj['recon_N_p'] = len(proj['angles']) #total number of angles to be used	
	print 'proj_init: Total number of projections used for reconstruction is ' + str(proj['recon_N_p'])

	return proj

""" Variable 'recon' is a dictionary datatype containing information about the object (or sample) and the reconstruction algorithm.
	Dictionary Key - Explanation.
	recon_type - MBIR or FBP
	num_threads - Number of threads to be used the reconstruction code (its a C code using OpenMP).
	voxel_thresh, cost_thresh, delta_xy, initICD, writeTiff are lists, each element of which corresponds to one stage in multiresolution reconstruction.
	voxel_thresh - Threshold in units of HU at which the algorithm converges. The average value of absolute change in voxel value is compared to threshold.
	cost_thresh - Threshold on a certain function of the cost to better enforce convergence. Its a percentage (range is 0 - 100). Not used if cost is not calculated
	delta_xy - Voxel size as a multiple of detector bin size at which reconstruction should be done.
	initICD - Value of '0' means initialize the object to '0'. If '1', reads in object.bin, interpolates by a factor of 2 by pixel replication and initializes the object. If '2' interpolates by a factor of 2 by bilinear interpolation and initializes the object.
	writeTiff - If '0', reconstruction code does not write the reconstructions to tiff files. If '1', it does.
	radius_obj - Radius from the center at which voxels are updated (time saving measure if we know the field of view of object)
	r, c_s, c_t, sigma_s, sigma_t are lists. Each corresponing value launches a separate reconstruction with those parameters. All should have same length
	r - Number of reconstructions in one full frame.
	c_s - Spatial qGGMRF prior parameter (to start with just keep it atleast 100 times less than abs(average_voxel_value/sigma_s)
	c_t - Temporal qGGMRD prior parameter (same as for c_s)
	sigma_s - Spatial qGGMRF regularization parameter. Lower the value higher the regularization (or higher spatial blur)
	sigma_t - Temporal qGGMRF regularization parameter. Lower the value higher the regularization (or higher temporal blur)
	p - qGGMRF parameter which controls sharpness (Keep p = 1.2)
	iterations - Maximum limit on total number of iterations in each multiresolution stage (Expectation is that algorithm will converge well before this)  
	alpha - Overrelaxation parameter which speeds up convergence (Should be less than 2)
	time_reg - If '1' does time regularization (Keep it at 1)
	Rtime0 - Starting time of first time slice of object.
	Rtime_num - Number of time slices of object
	Rtime_delta - Time between time slices of object
	multi_res_stages - Number of multi-resolution stages
	N_xy - Resolution of object along x and y axis (two axis perpendicular to axis of rotation)
	maxHU - Maximum value of attenuation coefficient displayed in a image
	minHU - Minimum value of attenuation coefficient displayed in a image
	reconstruct - If '1', will do reconstruction (obviously, keep it at 1)"""

""" r, c_s, c_t, sigma_s, sigma_t are lists. Each corresponding item in the lists will be used to run a instance of reconstruction. """

def recon_init (proj, recon):
	recon['recon_type'] = 'MBIR'
	
	recon['r'] = [1] #recon per frame 
	recon['c_s'] = [10**-6] #10^-6
	recon['c_t'] = [10**-4]

	recon['sigma_s'] = [8*(10**5)] #need to automatically set. To Do
	recon['sigma_t'] = [(10**2)] #Ignored for 3d recon
	
	recon['ZingerT'] = [100000] #Need to set automatically
	recon['ZingerDel'] = [0.1]

	recon['init_object4mHDF'] = 0
	
	recon['maxHU'] = 20000 
	recon['minHU'] = 0
	
	recon['radius_obj'] = 0.65*proj['N_r']/2 #Used or not used?
	recon['BH_Quad_Coef'] = 0;#need to make zero
	
	#recon['voxel_thresh'] = [10]
        #recon['cost_thresh'] = [10]
        #recon['delta_xy'] = [1]
        #recon['delta_z'] = [1]
        #recon['initICD'] = [2]
        #recon['sinobin'] = 1 
        #recon['writeTiff'] = [1]
        #recon['WritePerIter'] = [1]
        #recon['updateProjOffset'] = [3]
        #recon['iterations'] = [5]
        #recon['only_Edge_Updates'] = [0]
        #recon['initMagUpMap'] = [1]
	
	recon['voxel_thresh'] = [10, 20, 35, 35] #4 stage multi-resolution, with stopping in HU
        recon['cost_thresh'] = [10, 10, 10, 10] #percentage change presnt-prev / present - initial
        recon['delta_xy'] = [8, 4, 2, 1] #Multi-res multiple of det pixel size
        recon['delta_z'] = [8, 4, 2, 1] #Mutli-res multi-resolution 
        recon['initICD'] = [0, 3, 3, 3] #upsampling factor; 0 - no umpsample , 2 - xy upsampling, 3 - xy,z upsampling
        recon['sinobin'] = 1 #1-multi-res or 3-mult-grid  
        recon['writeTiff'] = [1, 1, 1, 1] #1 writes upon termination
        recon['WritePerIter'] = [0, 0, 0, 1] #Writes after each iteration
        recon['updateProjOffset'] = [0, 2, 3, 3] #update gain fluction 0 - no estimation, 1 - initialize and not estimated, 2 - not read but estimated , 3 initialized and estimated
        recon['iterations'] = [30, 20, 10, 10] #max iter
        recon['only_Edge_Updates'] = [0, 0, 0, 0] #DO NOT USE
        recon['initMagUpMap'] = [0, 1, 1, 1] #Update Mag Map 
	
	recon['init_with_FBP'] = 0 #DO NOT USE
	#recon['num_threads'] = 1
	recon['positivity_constraint'] = 0;#0 means no positivity ; 1 positivty used
	
	recon['p'] = 1.2
	recon['alpha'] = 1.5
	recon['time_reg'] = 1 #0 disable, 1 enable

	recon['Rtime0'] = proj['times'][0] #always 1 for 3-D recon;
	recon['Rtime_num'] = [ceil(recon['r'][i]*float(proj['recon_N_p'])/proj['N_theta']) for i in range(len(recon['r']))]
#	recon['Rtime_num'] = proj['r']*proj['N_p']/proj['N_theta']
	recon['Rtime_delta'] = [(proj['times'][-1]-proj['times'][0])/recon['Rtime_num'][i] for i in range(len(recon['r']))]
#	recon['Rtime_delta'] = proj['N_p']*proj['delta_time']/recon['Rtime_num']

	recon['multi_res_stages'] = len(recon['delta_xy'])
	recon['N_xy'] = proj['recon_N_r']/recon['delta_xy'][-1] #-1 means last elemen in the list
	recon['N_z'] = proj['recon_N_t']/recon['delta_z'][-1]
	
	recon['calculate_cost'] = 0 #0 for no 1 for yes
	recon['set_up_launch_folder'] = 0
	recon['NHICD'] = 1 #Enable or disable NHICD algorithm

	recon['FBP_N_xy'] = proj['recon_N_r']
	if (recon['init_with_FBP'] == 1):
		recon['FBP_N_xy'] = recon['FBP_N_xy']/recon['delta_xy'][0]

	if (proj['N_t'] % proj['recon_N_t'] != 0):
		error_by_flag(1, 'ERROR: recon_init: recon_N_t must divide N_t')

	recon['updateProjOffset'] = np.asarray(recon['updateProjOffset'])	
	if (len(recon['r']) != len(recon['c_s']) or len(recon['r']) != len(recon['c_t']) or len(recon['r']) != len(recon['sigma_s']) or len(recon['r']) != len(recon['sigma_t'])):
		error_by_flag (1, 'ERROR: recon_init: Lengths of r, c_t, c_s, sigma_s, sigma_t does not match')

	if (len(recon['delta_xy']) != len(recon['voxel_thresh']) or len(recon['delta_xy']) != len(recon['cost_thresh']) or len(recon['delta_xy']) != len(recon['initICD']) or len(recon['delta_xy']) != len(recon['writeTiff']) or len(recon['delta_xy']) != len(recon['delta_z'])):
		error_by_flag (1, 'ERROR: recon_init: Lengths of delta_xy, voxel_thresh, cost_thresh, initICD, writeTiff, sigma_t does not match')
	
	return recon

""" Variable 'files' is dictionary datatype containing the path to different folders 
	C_Source_Folder - Path to the folder containing the source code for reconstruction
	Launch_Folder - Path to the folder where the reconstruction code will be run
	Results_Folder - Folder to which the results will be copied
	copy_executables - If '1', copies the executables from the source code folder and hence is not compiled again (Keep it at 0)
	copy_projections - If '1', copies projection.bin and weight.bin from source code folder. If '0' reads the projection data from HDF files, as is described in XT_Projections.py"""

def files_init (files):
	files['C_Source_Folder'] = "../Source_Code_3D/"
	#files['Result_Folder'] = "../XT_Result_Repository/"
	files['Result_Folder'] = files['scratch'] + "/Recon_Runs/LBNL_Recons/XT_Result_Repository/" #Unncessary?
#	files['Result_Folder'] = files['scratch'] + "/Results/"

	files['Proj_Offset_File'] = "../Source_Code_3D/proj_offset.bin" #Not used if 0 is mult-res gain parameter estimation
	files['Launch_Folder'] = files['scratch'] + "/Recon_Runs/LBNL_Recons/XT_run/" #input by programmers
	#files['Launch_Folder'] = "../XT_run/"
#        files['Launch_Folder'] = files['scratch'] + "/LaunchFolder/"
	files['copy_executables'] = 0 #0 - dont exec, copy code + compile; TO DO : Test if we can only copy this 
	files['copy_projections'] = 0 #0 always for 3D

	return files
