
from scipy import interpolate
from scipy.misc import imresize
import scipy
import numpy as np
from XT_IOMisc import error_by_flag
import struct
from XT_Projections import decimate_count_data_in_r
from XT_Projections import decimate_count_data_in_t
#from XT_Conversions import convert_mu2HU
from XT_IOMisc import *
import scipy.io as sio
import os

#import cv2,cv

#def write_array2video (proj,recon,files):
#	for i in range(len(recon['r'])):
#		path2results = files['Results_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_L_' + str(proj['L']) + '_K_' + str(proj['K']) + '_N_p_' + str(proj['N_p']) + '/'
#		cv2.VideoWriter.open(path2results+'object.avi', cv2.CV_FOURCC('M','J','P','G'), files['video_framerate'], (recon['N_xy'], recon['N_xy']), isColor=0) 
#			error_by_flag(1, 'ERROR: write_array2video: Cannot open VideoWriter object') 
#		for k in range(recon['N_p']):
#			frame = convert_mu2HU2uint8(recon['object'][k::])
#			cv2.VideoWriter.write(frame)
	
def compute_RMSE_of_recon (proj, recon, files):
	if (proj['phantom_N_xy'] % recon['N_xy'] != 0):
		error_by_flag(1,"ERROR: N_xy should divide phantom_N_xy")
	RMSE_EDGE = np.zeros((len(recon['r'])), dtype=np.float64, order='C')
	RMSE_FULL = np.zeros((len(recon['r'])), dtype=np.float64, order='C')
	for i in range(len(recon['r'])):
		start = recon['Rtime0'] + recon['Rtime_delta'][i]/2
		stop = start + recon['Rtime_delta'][i]*(recon['Rtime_num'][i]-1)
		obj_times = np.linspace(start, stop, recon['Rtime_num'][i], endpoint=True)

		#path2results = files['Results_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_L_' + str(proj['L']) + '_K_' + str(proj['K']) + '_N_p_' + str(proj['N_p']) + '/'
		if (recon['recon_type'] == 'MBIR'):
			path2launch = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
			path2results = files['Result_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
		elif (recon['recon_type'] == 'FBP'):	
			path2launch = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
			path2results = files['Result_Folder'] + 'FBP' + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'

		node = recon['zSlice4RMSE']/(recon['N_z']/recon['node_num'])
		slice = recon['zSlice4RMSE'] % (recon['N_z']/recon['node_num'])
		print 'compute_RMSE_of_recon: RMSE will be computed from slice ' + str(slice) + ' from node ' + str(node)
	
		obj = np.zeros((recon['Rtime_num'][i], recon['N_xy'], recon['N_xy']), dtype=np.float64, order='C');	
		for j in range(int(recon['Rtime_num'][i])):	
			fid = open(path2launch + 'object_n' + str(node) + '_time_' + str(j) + '.bin', 'rb')
			data = RealData4mBin(path2launch + 'object_n' + str(node) + '_time_' + str(j) + '.bin', slice*recon['N_xy']*recon['N_xy'], recon['N_xy']*recon['N_xy'], recon['real_var_type'])
			obj[j,:,:] = np.reshape(data, (recon['N_xy'], recon['N_xy']), order='C')

		f_interp = interpolate.interp1d(obj_times, obj, kind='cubic', axis=0, copy=False)
		proj_times = proj['times'][recon['Proj0RMSE']:recon['Proj0RMSE']+recon['ProjNumRMSE']]
		interped_object = f_interp(proj_times)
	
		if (proj['slice_t_start'] + proj['N_t'] > proj['phantom_N_z']):
			error_by_flag(1, "ERROR: slice_t_start + N_t should be less than phantom_N_z")
	
		phantom = np.zeros((recon['ProjNumRMSE'], recon['N_xy'], recon['N_xy']), dtype=np.float32, order='C');	
		Object2Bin = np.zeros((recon['ProjNumRMSE'], recon['N_xy'], recon['N_xy']), dtype=np.float32, order='C');	
		Phantom2Bin = np.zeros((recon['ProjNumRMSE'], recon['N_xy'], recon['N_xy']), dtype=np.float32, order='C');	
		mask = np.zeros((recon['ProjNumRMSE'], recon['N_xy'], recon['N_xy']), dtype=np.uint8, order='C');	
		fid_phtm = open(proj['Path2Phantom'],'rb')
		fid_mask = open(proj['Path2Mask'],'rb')
		t_ratio = proj['N_t']/proj['recon_N_t']
		z_ratio = proj['N_t']/recon['N_z']
		index = 0
		rmse_edge = 0
		rmse_full = 0
		count_mask = 0
		for j in range(recon['Proj0RMSE'],recon['Proj0RMSE']+recon['ProjNumRMSE']):
			offset = j*proj['phantom_N_z']*proj['phantom_N_xy']*proj['phantom_N_xy']
			offset = offset + recon['zSlice4RMSE']*proj['phantom_N_xy']*proj['phantom_N_xy']*z_ratio
			offset = offset + proj['slice_t_start']*proj['phantom_N_xy']*proj['phantom_N_xy']
			
			fid_phtm.seek(offset*4, 0)
			data = fid_phtm.read(proj['phantom_N_xy']*proj['phantom_N_xy']*z_ratio*4)
			data = struct.unpack(str(proj['phantom_N_xy']*proj['phantom_N_xy']*z_ratio)+'f',data)
			data = np.reshape(data, (z_ratio, proj['phantom_N_xy'], proj['phantom_N_xy']), order='C')
			data = np.mean(data, axis=0, dtype=np.float64)
			data = decimate_count_data_in_r(data, proj['phantom_N_xy'], recon['N_xy'])/proj['phantom_N_xy']*recon['N_xy']
			data = decimate_count_data_in_t(data, proj['phantom_N_xy'], recon['N_xy'])/proj['phantom_N_xy']*recon['N_xy']
			phantom = data
			
			fid_mask.seek(offset, 0)
			data = fid_mask.read(proj['phantom_N_xy']*proj['phantom_N_xy']*z_ratio)
			data = struct.unpack(str(proj['phantom_N_xy']*proj['phantom_N_xy']*z_ratio)+'B',data)
			data = np.reshape(data, (z_ratio, proj['phantom_N_xy'], proj['phantom_N_xy']), order='C')
			data = np.mean(data, axis=0, dtype=np.float64)
			data = decimate_count_data_in_r(data, proj['phantom_N_xy'], recon['N_xy'])/proj['phantom_N_xy']*recon['N_xy']
			data = decimate_count_data_in_t(data, proj['phantom_N_xy'], recon['N_xy'])/proj['phantom_N_xy']*recon['N_xy']
			mask = data > 0
			
			arr2tif = np.hstack((phantom, interped_object[index,:,:]))
			arr2tif = convert_um2HU (arr2tif)
			arr2tif = np.hstack((arr2tif, recon['minHU']+(recon['maxHU']-recon['minHU'])*(mask.astype(np.float64))))
			arr2tif = convert_HU2uint8(arr2tif, recon['maxHU'], recon['minHU'])
			write_array2tif (path2results + 'ObjVsPhantom_' + str(index) + '.tif', arr2tif)
			
			obj_temp = interped_object[index,:,:]
			rmse_edge = rmse_edge + np.sum((phantom[mask]-obj_temp[mask])*(phantom[mask]-obj_temp[mask]),dtype=np.float64)
			rmse_full = rmse_full + np.sum((phantom-obj_temp)*(phantom-obj_temp),dtype=np.float64)
			count_mask = count_mask + np.sum(mask.astype(np.float64),dtype=np.float64)
			Object2Bin[index,:,:] = obj_temp.astype(np.float32)
			Phantom2Bin[index,:,:] = phantom.astype(np.float32)
			index = index + 1			

		fid_phtm.close()
		fid_mask.close()
#		write_Video4mArray(path2launch+'object_interpolated.mp4', recon['N_xy'], recon['N_xy'], recon['r'][i], interped_object)
#		write_Video4mArray(path2launch+'phantom.mp4', recon['N_xy'], recon['N_xy'], recon['r'][i], phantom)
		RMSE_EDGE[i] = np.sqrt(rmse_edge/count_mask)
		RMSE_EDGE[i] = convert_um2HU(RMSE_EDGE[i])
		
		RMSE_FULL[i] = np.sqrt(rmse_full/(recon['N_xy']*recon['N_xy']*recon['ProjNumRMSE']))
		RMSE_FULL[i] = convert_um2HU(RMSE_FULL[i])
		
		text_file = open(path2results + 'RMSE.txt', 'w')
		text_file.write('RMSE_EDGE = ' + str(RMSE_EDGE[i]) + '\n')
		text_file.write('RMSE = ' + str(RMSE_FULL[i]) + '\n')
		text_file.write('z slice = ' + str(recon['zSlice4RMSE']) + '\n')
		text_file.write('K = ' + str(proj['K']) + '\n')
		text_file.write('N_theta = ' + str(proj['N_theta']) + '\n')
		text_file.write('r = ' + str(recon['r'][i]) + '\n')
		text_file.write('sigma_s = ' + str(recon['sigma_s'][i]) + '\n')
		text_file.write('sigma_t = ' + str(recon['sigma_t'][i]) + '\n')
		text_file.write('c_s = ' + str(recon['c_s'][i]) + '\n')
		text_file.write('c_t = ' + str(recon['c_t'][i]) + '\n')
		text_file.write('ZingerT = ' + str(recon['ZingerT'][i]) + '\n')
		text_file.write('ZingerDel = ' + str(recon['ZingerDel'][i]) + '\n')
		text_file.write('Proj0RMSE = ' + str(recon['Proj0RMSE']) + '\n')
		text_file.write('ProjNumRMSE = ' + str(recon['ProjNumRMSE']) + '\n')
		text_file.write('Real Variable Type = ' + str(recon['real_var_type']) + '\n')
		text_file.close()
		
		fid = open(path2results + 'Reconstruction.bin', 'wb')
		Object2Bin.tofile(fid)		
		fid.close()	
		
		fid = open(path2results + 'Phantom.bin', 'wb')
		Phantom2Bin.tofile(fid)		
		fid.close()	

		Data2Mat = {}	
		Data2Mat['RMSE_EDGE'] = RMSE_EDGE[i]
		Data2Mat['RMSE'] = RMSE_FULL[i]
		Data2Mat['recon_type'] = recon['recon_type']
		Data2Mat['Z_Slice'] = recon['zSlice4RMSE']
		Data2Mat['K'] = proj['K']
		Data2Mat['N_theta'] = proj['N_theta']
		Data2Mat['r'] = recon['r'][i]
		Data2Mat['sigma_s'] = recon['sigma_s'][i]
		Data2Mat['sigma_t'] = recon['sigma_t'][i]
		Data2Mat['c_s'] = recon['c_s'][i]
		Data2Mat['c_t'] = recon['c_t'][i]
		Data2Mat['ZingerT'] = recon['ZingerT'][i]
		Data2Mat['ZingerDel'] = recon['ZingerDel'][i]
		Data2Mat['Proj0RMSE'] = recon['Proj0RMSE']
		Data2Mat['ProjNumRMSE'] = recon['ProjNumRMSE']
		Data2Mat['real_var_type'] = recon['real_var_type']
		Data2Mat['fineres_voxthresh'] = recon['voxel_thresh'][-1]
		#create_param_sweep_file (Data2Mat,files)
		
		sio.savemat(path2results + 'Params.mat', Data2Mat)
	print RMSE_EDGE
	print RMSE_FULL
		#	a = f_interp(proj_times[512])
		#	print a[256,256]
		#	print interped_object	
	
def create_param_sweep_file (Data2Mat,files):
	path2results = files['Result_Folder']
	if (os.path.isfile(path2results + 'Sweep.mat')):
		sweep = sio.loadmat(path2results + 'Sweep.mat')
		for i in range(len(sweep['K'])):
			if (sweep['recon_type'] == Data2Mat['recon_type'] and sweep['Z_Slice'][i] == Data2Mat['Z_Slice'] and sweep['K'][i] == Data2Mat['K'] and sweep['N_theta'][i] == Data2Mat['N_theta'] and sweep['r'][i] == Data2Mat['r'] and sweep['sigma_s'][i] == Data2Mat['sigma_s'] and sweep['sigma_t'][i] == Data2Mat['sigma_t'] and sweep['c_s'][i] == Data2Mat['c_s'] and sweep['c_t'][i] == Data2Mat['c_t'] and sweep['ZingerT'][i] == Data2Mat['ZingerT'] and sweep['ZingerDel'][i] == Data2Mat['ZingerDel'] and sweep['Proj0RMSE'][i] == Data2Mat['Proj0RMSE'] and sweep['ProjNumRMSE'][i] == Data2Mat['ProjNumRMSE'] and sweep['fineres_voxthresh'][i] == Data2Mat['fineres_voxthresh']):
				copy_data2sweep(sweep,Data2Mat,i)
		copy_data2sweep(sweep,Data2Mat,-1)
	else:
		sweep = {'RMSE_EDGE':[],'RMSE':[],'recon_type':[],'Z_Slice':[],'K':[],'N_theta':[],'r':[],'sigma_s':[],'sigma_t':[],'c_s':[],'c_t':[],'ZingerT':[],'ZingerDel':[],'Proj0RMSE':[],'ProjNumRMSE':[],'fineres_voxthresh':[]}
		copy_data2sweep(sweep,Data2Mat,-1)
	sio.savemat(path2results + 'Sweep.mat', sweep)
	

def copy_data2sweep (sweep, Data2Mat, add_or_rem):
	if (add_or_rem == -1):	
		sweep['RMSE_EDGE'].append(Data2Mat['RMSE_EDGE'])
		sweep['RMSE'].append(Data2Mat['RMSE'])
		sweep['recon_type'].append(Data2Mat['recon_type'])
		sweep['Z_Slice'].append(Data2Mat['Z_Slice'])
		sweep['K'].append(Data2Mat['K'])
		sweep['N_theta'].append(Data2Mat['N_theta'])
		sweep['r'].append(Data2Mat['r'])
		sweep['sigma_s'].append(Data2Mat['sigma_s'])
		sweep['sigma_t'].append(Data2Mat['sigma_t'])
		sweep['c_s'].append(Data2Mat['c_s'])
		sweep['c_t'].append(Data2Mat['c_t'])
		sweep['ZingerT'].append(Data2Mat['ZingerT'])
		sweep['ZingerDel'].append(Data2Mat['ZingerDel'])
		sweep['Proj0RMSE'].append(Data2Mat['Proj0RMSE'])
		sweep['ProjNumRMSE'].append(Data2Mat['ProjNumRMSE'])
		sweep['fineres_voxthresh'].append(Data2Mat['fineres_voxthresh'])
	else:
		sweep['RMSE_EDGE'].pop(add_or_rem)
		sweep['RMSE'].pop(add_or_rem)
		sweep['recon_type'].pop(add_or_rem)
		sweep['Z_Slice'].pop(add_or_rem)
		sweep['K'].pop(add_or_rem)
		sweep['N_theta'].pop(add_or_rem)
		sweep['r'].pop(add_or_rem)
		sweep['sigma_s'].pop(add_or_rem)
		sweep['sigma_t'].pop(add_or_rem)
		sweep['c_s'].pop(add_or_rem)
		sweep['c_t'].pop(add_or_rem)
		sweep['ZingerT'].pop(add_or_rem)
		sweep['ZingerDel'].pop(add_or_rem)
		sweep['Proj0RMSE'].pop(add_or_rem)
		sweep['ProjNumRMSE'].pop(add_or_rem)
		sweep['fineres_voxthresh'].pop(add_or_rem)
			
			
