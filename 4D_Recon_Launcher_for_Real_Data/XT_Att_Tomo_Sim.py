
import numpy as np


def attenuation_tomo_sim_init (proj, recon, files):
	proj['Path2Phantom'] = files['data_scratch'] + "/Sim_Datasets/phantom_Cahn_Hilliard.bin"
	proj['Expected_Counts'] = 29473 
	proj['phantom_N_xy'] = 2048
	# phantom_N_z is the resolution of phantom along z
	proj['phantom_N_z'] = 4
	proj['voxel_size'] = 0.65

	# recon_N_r is detector resolution along r-axis used in reconstruction
	proj['recon_N_r'] = 512
	# slice_t_start is first slice in phantom where recon starts
	proj['slice_t_start'] = 0
	# N_t is number of phantom slices extracted for reconstruction
	# Assumes same voxel size for phantom and detector (pixel size)
	proj['N_t'] = 4
	# Subsamples N_t to recon_N_t when doing reconstruction
	proj['recon_N_t'] = 4
	# Same units as recon_N_r
	proj['rotation_center_r'] = 256
	proj['proj_start'] = 0
	proj['proj_num'] = 512
	# N_p is just used for angle generation in python. Should be greater than proj_start + proj_num
	proj['N_p'] = 512
	proj['K'] = 16
	proj['N_theta'] = 512

	recon['sigma_s'] = [10*(10**5)]
	recon['sigma_t'] = [45*(10**4)]
	recon['r'] = 16
	recon['c_s'] = 10**-6
	recon['c_t'] = 10**-6
	recon['ZingerT'] = 100000
	recon['ZingerDel'] = 0.1
	recon['maxHU'] = 43000
	recon['minHU'] = 5000
	recon['voxel_thresh'] = [0.5, 0.5, 0.5, 0.5]
        recon['cost_thresh'] = [10, 10, 10, 10]
        recon['delta_xy'] = [8, 4, 2, 1]
        recon['delta_z'] = [1, 1, 1, 1]
        recon['initICD'] = [0, 2, 2, 2]
        recon['WritePerIter'] = [0, 0, 0, 1]
        recon['updateProjOffset'] = [0, 2, 3, 3]
        recon['iterations'] = [500, 400, 400, 400]
	recon['do_VarEstimate'] = [1]*len(recon['voxel_thresh'])
	recon['Estimate_of_Var'] = 0.124;	
	
	files['Result_Folder'] = files['scratch'] + "/Recon_Runs/Recon_New_Sim_Data/XT_Result_Repository/"
	files['Launch_Folder'] = files['scratch'] + "/Recon_Runs/Recon_New_Sim_Data/XT_run/"

	proj['N_r'] = proj['phantom_N_xy']
	proj['min_time_btw_views'] = 0
	proj['rotation_speed'] = 100
	recon['recon_type'] = 'MBIR'
        recon['sinobin'] = 2
        recon['initMagUpMap'] = [0, 1, 1, 1]
        recon['only_Edge_Updates'] = [0, 0, 0, 0]
	recon['readSino4mHDF'] = [0, 0, 0, 0]	
        recon['writeTiff'] = [1, 1, 1, 1]
	recon['BH_Quad_Coef'] = 0;
	# Two variables below will be ignored by C-code in this mode
	proj['Path2Dataset'] = files['data_scratch'] + "/blah/"
	proj['Path2WhiteDark'] = files['data_scratch'] + "/blah/"
	
	return (proj, recon, files)	
