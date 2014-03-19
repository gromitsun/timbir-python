
import numpy as np


def attenuation_tomo_real_init (proj, recon, files):
	proj['Path2Dataset'] = files['data_scratch'] + "/Argonne_Datasets/K_16_N_theta_2000_RotSpeed_100_Exp_4_ROI_1000x2080_Ramp_2/k-16-4ms-last_22.hdf"
	proj['Path2WhiteDark'] = files['data_scratch'] + "/Argonne_Datasets/K_16_N_theta_2000_RotSpeed_100_Exp_4_ROI_1000x2080_Ramp_2/k-16-4ms-last_31.hdf"
	proj['voxel_size'] = 0.65
	# recon_N_r is detector resolution along r-axis used in reconstruction. Subsampled from a resolution of N_r (actual detector resolution)
	proj['recon_N_r'] = 512
	# slice_t_start is first detector slice used in recon
	proj['slice_t_start'] = 0
	# N_t is number of detector slices used in recon
	proj['N_t'] = 4*4
	# Subsamples N_t to recon_N_t when doing reconstruction
	proj['recon_N_t'] = 4
	# Same units as recon_N_r
	proj['rotation_center_r'] = 264.75
	proj['proj_start'] = 998
	proj['proj_num'] = 2000
	# N_p is just used for angle generation in python. Should be greater than proj_start + proj_num
	proj['N_p'] = 2000*4
	proj['K'] = 16
	proj['N_theta'] = 2000
	proj['N_r'] = 2080
	proj['min_time_btw_views'] = 0.014039
	proj['rotation_speed'] = 100

	recon['sigma_s'] = [150*(10**5)]
	recon['sigma_t'] = [2*(10**4)]
	recon['r'] = 16
	recon['c_s'] = 10**-6
	recon['c_t'] = 10**-6
	recon['ZingerT'] = 4
	recon['ZingerDel'] = 0.1
	recon['maxHU'] = 60000
	recon['minHU'] = 10000
	recon['voxel_thresh'] = [0.5, 0.5, 0.5, 0.5]
        recon['cost_thresh'] = [10, 10, 10, 10]
        recon['delta_xy'] = [8, 4, 2, 1]
        recon['delta_z'] = [1, 1, 1, 1]
        recon['initICD'] = [0, 2, 2, 2]
        recon['WritePerIter'] = [0, 0, 0, 1]
        recon['updateProjOffset'] = [0, 2, 3, 3]
	recon['readSino4mHDF'] = [1, 0, 0, 0]	
        recon['iterations'] = [500, 400, 400, 400]
	recon['do_VarEstimate'] = [1]*len(recon['voxel_thresh'])
	recon['Estimate_of_Var'] = 1;	
	
	files['Result_Folder'] = files['scratch'] + "/Recon_Runs/Recon_Att_Real_New/XT_Result_Repository/"
	files['Launch_Folder'] = files['scratch'] + "/Recon_Runs/Recon_Att_Real_New/XT_run/"

	recon['recon_type'] = 'MBIR'
        recon['sinobin'] = 1
        recon['initMagUpMap'] = [0, 1, 1, 1]
        recon['only_Edge_Updates'] = [0, 0, 0, 0]
        recon['writeTiff'] = [1, 1, 1, 1]
	recon['BH_Quad_Coef'] = 0.5;
	# Two variables below will be ignored by C-code in this mode
	
	proj['Expected_Counts'] = 29473 
	proj['phantom_N_xy'] = 2048
	# phantom_N_z is the resolution of phantom along z
	proj['phantom_N_z'] = 4
	proj['Path2Phantom'] = files['data_scratch'] + "/Sim_Datasets/blah.bin"

	return proj, recon, files
