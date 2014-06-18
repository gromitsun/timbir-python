
import numpy as np
from XT_IOMisc import error_by_flag


def attenuation_tomo_sim_init (proj, recon, files):
	proj['Expected_Counts'] = 29473 
	proj['phantom_N_xy'] = 1024
	# phantom_N_z is the resolution of phantom along z
	proj['phantom_N_z'] = 32

	if (proj['slice_t_start'] + proj['N_t'] > proj['phantom_N_z']):
		error_by_flag(1,"proj['slice_t_start'] + proj['N_t'] > proj['phantom_N_z']")
	# Subsamples N_t to recon_N_t when doing reconstruction
	recon['Proj0RMSE'] = 256
	recon['ProjNumRMSE'] = 256*2

	# reg params for N_theta = 16, K = 1, r = 1
	# Minimum for 16*(10**5) & 4*(10**5)
	# recon['sigma_s'] = [2*(10**5), 4*(10**5), 8*(10**5), 16*(10**5)]
	# recon['sigma_t'] = [5*(10**4), (10**5), 2*(10**5), 4*(10**5)]

	# reg params for N_theta = 256, K = 1, r = 1
	# Minimum for 8*(10**5) & 10**4
	#recon['sigma_s'] = [25*(10**4), 5*(10**5), (10**6), 2*(10**6), 4*(10**6)]
	#recon['sigma_t'] = [2000, 1000, 500, 250, 100, 50]
	
	# reg params for N_theta = 256, K = 16, r = 16
	# Minimum for 10**6 & 4*(10**4)
	#recon['sigma_s'] = [25*(10**4), 5*(10**5), (10**6), 2*(10**6), 4*(10**6)]
	#recon['sigma_t'] = [16000, 8000, 4000, 2000, 1000]

	recon['sigma_s'] = [25*(10**4), 5*(10**5), 10*(10**5), 2*(10**6), 4*(10**6)]
	recon['sigma_t'] = [125, 250, 500, 1000, 2*(10**3)]
	#recon['sigma_s'] = [25*(10**4), 5*(10**5), 10*(10**5), 2*(10**6)]
	#recon['sigma_t'] = [2*(10**3), 4*(10**3), 8*(10**3), 16*(10**3), 32*(10**3)]
	recon['c_s'] = 10**-6
	recon['c_t'] = 10**-6
	recon['ZingerT'] = 4
	recon['ZingerDel'] = 0.1
	recon['maxHU'] = 40000
	recon['minHU'] = 0
	
	recon['voxel_thresh'] = [0.05, 0.05, 0.05]
        recon['cost_thresh'] = [0.1, 0.1, 0.1]
        recon['delta_xy'] = [4, 2, 1]
        recon['delta_z'] = [1, 1, 1]
        recon['initICD'] = [0, 2, 2]
        recon['WritePerIter'] = [0, 0, 1]
        recon['updateProjOffset'] = [0, 2, 3]
        recon['iterations'] = [500, 500, 400]
	recon['do_VarEstimate'] = [0]*len(recon['voxel_thresh'])
	recon['Estimate_of_Var'] = 1.0;	
	
	if (proj['N_r'] == proj['phantom_N_xy']):
		error_by_flag(1,"ERROR: N_r is not equal to phantom_N_xy")
	proj['min_time_btw_views'] = 0
	proj['rotation_speed'] = 100
	recon['recon_type'] = 'MBIR'
        recon['sinobin'] = 2
        recon['initMagUpMap'] = [1]*len(recon['voxel_thresh'])
	recon['initMagUpMap'][0] = 0
        recon['only_Edge_Updates'] = [0]*len(recon['voxel_thresh'])
	recon['readSino4mHDF'] = [0]*len(recon['voxel_thresh'])	
        recon['writeTiff'] = [1]*len(recon['voxel_thresh'])
	recon['BH_Quad_Coef'] = 0;
	# Two variables below will be ignored by C-code in this mode
	proj['Path2Dataset'] = files['data_scratch'] + "/blah/"
	proj['Path2WhiteDark'] = files['data_scratch'] + "/blah/"
	
	return (proj, recon, files)	
