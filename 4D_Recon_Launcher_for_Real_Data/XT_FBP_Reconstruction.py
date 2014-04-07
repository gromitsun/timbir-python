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

	for k in range(len(recon['Rtime_num'])):		
		path2launch = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][k]) + '_sigt_' + str(recon['sigma_t'][k]) + '_r_' + str(recon['r'][k]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
		projections = np.zeros((proj['recon_N_p'], proj['recon_N_r'], proj['recon_N_t']), dtype = np.float64, order = 'C')
		ratio = proj['recon_N_t']/recon['node_num']
		for rank in range(recon['node_num']):
			projections[:,:,rank*ratio:(rank+1)*ratio] = np.fromfile(path2launch + 'projection_n' + str(rank) + '.bin', dtype=np.float64, count=-1).reshape((proj['recon_N_p'],proj['recon_N_r'],ratio),order='C') 	
		centered_projections = gen_centered_projections(projections,floor(proj['rotation_center_r']))
		path2results = result_folder + 'FBP_' + 'r_' + str(recon['r'][k]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
		create_folder(path2results)
		
		for i in range(int(recon['Rtime_num'][k])):
			time_start = recon['Rtime0'] + i*recon['Rtime_delta'][k]
			time_stop = recon['Rtime0'] + (i + 1)*recon['Rtime_delta'][k]
			indices = np.bitwise_and(proj['times'] >= time_start, proj['times'] < time_stop)
			for k in range(recon['node_num']):
				ratio = proj['recon_N_t']/recon['node_num']
				idx = 0
				Obj2Bin = np.zeros((ratio, proj['recon_N_r'], proj['recon_N_r']), dtype=np.float64, order='C')
				for j in range(k*ratio, (k+1)*ratio):
					Object = iradon (np.transpose(centered_projections[indices,:, j]), 90-np.mod(proj['angles'][indices],360), output_size=proj['recon_N_r'], filter="hamming", interpolation='linear')/(proj['length_r']/proj['recon_N_r'])
					#Object = iradon (np.transpose(centered_projections[indices,:]), np.mod(proj['angles'][indices],180), output_size=centered_projections.shape[1], filter="hamming", interpolation='linear')/(proj['length_r']/proj['N_r'])
					Object = np.transpose(Object)
					Obj2Bin[idx,:,:] = Object
					idx = idx + 1
					Object = convert_um2HU(Object)
					Object_uint8 = convert_HU2uint8(Object, recon['maxHU'], recon['minHU'])
					write_array2tif(path2launch + 'object_time' + str(i) + '_z_' + str(j) + '.tif', Object_uint8)
				fid = open(path2launch + 'object_n' + str(k) + '_time_' + str(i)  + '.bin','wb')
				Obj2Bin.tofile(fid)
				fid.close()		
