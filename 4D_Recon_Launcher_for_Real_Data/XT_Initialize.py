#! /usr/bin/python

import numpy as np
from XT_Interlaced_Angles import gen_interlaced_views_0_to_Inf,clip_list_of_views
from math import ceil
from XT_IOMisc import error_by_flag
import os

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

def proj_init (proj, files):

#	proj['Path2Dataset'] = files['data_scratch'] + "/Argonne_Datasets/K_16_N_theta_2000_RotSpeed_100_Exp_4_ROI_1000x2080_Ramp_2/k-16-4ms-last_22.hdf"
#	proj['Path2WhiteDark'] = files['data_scratch'] + "/Argonne_Datasets/K_16_N_theta_2000_RotSpeed_100_Exp_4_ROI_1000x2080_Ramp_2/k-16-4ms-last_31.hdf"
#	proj['Path2Dataset'] = files['data_scratch'] + "/Argonne_Datasets/K_32_N_theta_1984_RotSpeed_100_Exp_8_ROI_2000x2080_Ramp_5/k-32-08ms_1.hdf"
#	proj['Path2WhiteDark'] = files['data_scratch'] + "/Argonne_Datasets/K_32_N_theta_1984_RotSpeed_100_Exp_8_ROI_2000x2080_Ramp_5/k-32-08ms_1.hdf"
	
	proj['length_r'] = proj['voxel_size']*proj['N_r']	
	proj['length_t'] = proj['voxel_size']*proj['N_t']	
	proj['L'] = proj['N_theta']/proj['K']
	
	fps_of_camera = proj['L']/(180.0/proj['rotation_speed'])
	angles, times = gen_interlaced_views_0_to_Inf(proj['K'], proj['N_theta'], proj['N_p'])
	angles_clip, times_clip, angles_del, times_del = clip_list_of_views (angles, times, proj['min_time_btw_views'], proj['rotation_speed'])
	print 'clip_list_of_views: Number of views deleted are ', angles_del.size, '. Deleted views are ' + str(angles_del)
	print 'clip_list_of_views: Deleted views\' times are ' + str(times_del)

	proj['angles'] = angles_clip[proj['proj_start'] : proj['proj_start'] + proj['proj_num']]
	proj['times'] = times_clip[proj['proj_start'] : proj['proj_start'] + proj['proj_num']]
	
	proj['recon_N_p'] = len(proj['angles']) 	
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
	recon['sigma_s'],recon['sigma_t'] = np.meshgrid(recon['sigma_s'], recon['sigma_t'])

	recon['sigma_s'] = recon['sigma_s'].flatten()
	recon['sigma_t'] = recon['sigma_t'].flatten()
	recon['r'] = [recon['r']]*len(recon['sigma_s'])
	recon['c_s'] = [recon['c_s']]*len(recon['r'])
	recon['c_t'] = [recon['c_t']]*len(recon['r'])
	recon['ZingerT'] = [recon['ZingerT']]*len(recon['r'])
	recon['ZingerDel'] = [recon['ZingerDel']]*len(recon['r'])
	
	param_idx = int(os.environ['PARAM_INDEX'])
	# if param_idx is 0, the one python process will run recon for all params
	if (param_idx > 0):
		recon['r'] = [recon['r'][param_idx-1]]
		recon['c_s'] = [recon['c_s'][param_idx-1]]
		recon['c_t'] = [recon['c_t'][param_idx-1]]
		recon['sigma_s'] = [recon['sigma_s'][param_idx-1]]
		recon['sigma_t'] = [recon['sigma_t'][param_idx-1]]
		recon['ZingerT'] = [recon['ZingerT'][param_idx-1]]
		recon['ZingerDel'] = [recon['ZingerDel'][param_idx-1]]
	
#	if (proj['recon_N_p']/proj['N_theta'] > 2):
#		recon['Proj0RMSE'] = proj['N_theta']
#		recon['ProjNumRMSE'] = (proj['recon_N_p']/proj['N_theta']-2)*proj['N_theta']
#	else:
#		recon['Proj0RMSE'] = proj['N_theta']/proj['K'] 
#		recon['ProjNumRMSE'] = proj['recon_N_p'] - 2*proj['N_theta']/proj['K']

	recon['init_object4mHDF'] = 0
	
	recon['radius_obj'] = proj['voxel_size']*proj['N_r']
	
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
	
	recon['init_with_FBP'] = 0
	#recon['num_threads'] = 1
	recon['positivity_constraint'] = 0;
	
	recon['p'] = 1.2
	recon['alpha'] = 1.5
	recon['time_reg'] = 1
	recon['Rtime0'] = proj['times'][0]

	recon['Rtime_num'] = [ceil(recon['r'][i]*float(proj['recon_N_p'])/proj['N_theta']) for i in range(len(recon['r']))]
#	recon['Rtime_num'] = proj['r']*proj['N_p']/proj['N_theta']
	recon['Rtime_delta'] = [(proj['times'][-1]-proj['times'][0])/recon['Rtime_num'][i] for i in range(len(recon['r']))]
#	recon['Rtime_delta'] = proj['N_p']*proj['delta_time']/recon['Rtime_num']

	recon['multi_res_stages'] = len(recon['delta_xy'])
	recon['N_xy'] = proj['recon_N_r']/recon['delta_xy'][-1]
	recon['N_z'] = proj['recon_N_t']/recon['delta_z'][-1]
	recon['zSlice4RMSE'] = recon['N_z']/2
	
	recon['calculate_cost'] = 0
	recon['set_up_launch_folder'] = 0
	recon['NHICD'] = 1


	recon['FBP_N_xy'] = proj['recon_N_r']
	if (recon['init_with_FBP'] == 1):
		recon['FBP_N_xy'] = recon['FBP_N_xy']/recon['delta_xy'][0]

	if (proj['N_t'] % proj['recon_N_t'] != 0):
		error_by_flag(1, 'ERROR: recon_init: recon_N_t must divide N_t')

	recon['updateProjOffset'] = np.asarray(recon['updateProjOffset'])	
	if (len(recon['r']) != len(recon['c_s']) or len(recon['r']) != len(recon['c_t']) or len(recon['r']) != len(recon['sigma_s']) or len(recon['r']) != len(recon['sigma_t'])):
		error_by_flag (1, 'ERROR: recon_init: Lengths of r, c_t, c_s, sigma_s, sigma_t does not match')

	if (len(recon['delta_xy']) != len(recon['voxel_thresh']) or len(recon['delta_xy']) != len(recon['cost_thresh']) or len(recon['delta_xy']) != len(recon['initICD']) or len(recon['delta_xy']) != len(recon['writeTiff']) or len(recon['delta_xy']) != len(recon['delta_z'])):
		error_by_flag (1, 'ERROR: recon_init: Lengths of delta_xy, voxel_thresh, cost_thresh, initICD, writeTiff does not match')
	
	return recon

""" Variable 'files' is dictionary datatype containing the path to different folders 
	C_Source_Folder - Path to the folder containing the source code for reconstruction
	Launch_Folder - Path to the folder where the reconstruction code will be run
	Results_Folder - Folder to which the results will be copied
	copy_executables - If '1', copies the executables from the source code folder and hence is not compiled again (Keep it at 0)
	copy_projections - If '1', copies projection.bin and weight.bin from source code folder. If '0' reads the projection data from HDF files, as is described in XT_Projections.py"""

def files_init (files):
	files['C_Source_Folder'] = "../Source_Code_4D_Wang/"
	files['Proj_Offset_File'] = "../Source_Code_4D_Wang/proj_offset.bin"
	files['copy_executables'] = 0
	files['copy_projections'] = 0

	return files
