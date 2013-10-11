import numpy as np
from XT_IOMisc import create_folder, convert_um2HU, convert_HU2uint8, write_array2tif
from skimage.transform import iradon
from math import floor

def gen_centered_projections (projections, rotation_center):
	centered_proj_sz = 2*rotation_center+1
	centered_projections = np.zeros((projections.shape[0], centered_proj_sz, projections.shape[2]), dtype = np.float64, order = 'C') 
	centered_projections[:, 0:min(centered_proj_sz,projections.shape[1]), :] = projections[:, 0:min(centered_proj_sz,projections.shape[1]), :]
	return centered_projections

def do_FBP_reconstruction (proj, recon, files):

########################### CHANGE RESULT FOLDER BELOW #########################
	result_folder = files['Result_Folder']
################################################################################

	centered_projections = gen_centered_projections(proj['projections'],floor(proj['rotation_center_r']))

	for k in range(len(recon['Rtime_num'])):		
		path2results = result_folder + 'FBP_' + 'r_' + str(recon['r'][k]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '/'
		create_folder(path2results)
		
		for i in range(int(recon['Rtime_num'][k])):
			time_start = recon['Rtime0'] + i*recon['Rtime_delta'][k]
			time_stop = recon['Rtime0'] + (i + 1)*recon['Rtime_delta'][k]
			indices = np.bitwise_and(proj['times'] >= time_start, proj['times'] < time_stop)
			for j in range(proj['recon_N_t']):
				Object = iradon (np.transpose(centered_projections[indices,:, j]), 90-np.mod(proj['angles'][indices],360), output_size=proj['recon_N_r'], filter="hamming", interpolation='linear')/(proj['length_r']/proj['recon_N_r'])
				#Object = iradon (np.transpose(centered_projections[indices,:]), np.mod(proj['angles'][indices],180), output_size=centered_projections.shape[1], filter="hamming", interpolation='linear')/(proj['length_r']/proj['N_r'])
				Object = convert_um2HU(np.transpose(Object))
				Object_uint8 = convert_HU2uint8(Object, recon['maxHU'], recon['minHU'])
				write_array2tif(path2results + 'object_time' + str(i) + '_z_' + str(j) + '.tif', Object_uint8)
