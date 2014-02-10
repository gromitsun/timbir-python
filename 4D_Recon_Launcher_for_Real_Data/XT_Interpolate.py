
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
	RMSE = np.zeros((len(recon['r'])), dtype=np.float64, order='C')
	for i in range(len(recon['r'])):
		start = recon['Rtime0'] + recon['Rtime_delta'][i]/2
		stop = start + recon['Rtime_delta'][i]*(recon['Rtime_num'][i]-1)
		obj_times = np.linspace(start, stop, recon['Rtime_num'][i], endpoint=True)

		#path2results = files['Results_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_L_' + str(proj['L']) + '_K_' + str(proj['K']) + '_N_p_' + str(proj['N_p']) + '/'
		path2launch = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + '/'
	
		node = recon['zSlice4RMSE']/(recon['N_z']/recon['node_num'])
		slice = recon['zSlice4RMSE'] % (recon['N_z']/recon['node_num'])
		print 'compute_RMSE_of_recon: RMSE will be computed from slice ' + str(slice) + ' from node ' + str(node)
	
		obj = np.zeros((recon['Rtime_num'][i], recon['N_xy'], recon['N_xy']), dtype=np.float64, order='C');	
		for j in range(int(recon['Rtime_num'][i])):	
			fid = open(path2launch + 'object_n' + str(node) + '_time_' + str(j) + '.bin', 'rb')
			fid.seek(slice*recon['N_xy']*recon['N_xy']*8, 0)
			# Multiply above offset by 8 which is size of double in bytes
			data = fid.read(recon['N_xy']*recon['N_xy']*8)
			data = struct.unpack(str(recon['N_xy']*recon['N_xy'])+'d',data)	
			obj[j,:,:] = np.reshape(data, (recon['N_xy'], recon['N_xy']), order='C')
			fid.close()

		f_interp = interpolate.interp1d(obj_times, obj, kind='cubic', axis=0, copy=False)
		proj_times = proj['times'][recon['Proj0RMSE']:recon['Proj0RMSE']+recon['ProjNumRMSE']]
		interped_object = f_interp(proj_times)
	
		if (proj['slice_t_start'] + proj['N_t'] > proj['phantom_N_z']):
			error_by_flag(1, "ERROR: slice_t_start + N_t should be less than phantom_N_z")
	
		phantom = np.zeros((recon['ProjNumRMSE'], recon['N_xy'], recon['N_xy']), dtype=np.float32, order='C');	
		fid = open(proj['Path2Phantom'],'rb')
		t_ratio = proj['N_t']/proj['recon_N_t']
		z_ratio = proj['N_t']/recon['N_z']
		index = 0
		for j in range(recon['Proj0RMSE'],recon['Proj0RMSE']+recon['ProjNumRMSE']):
			offset = j*proj['phantom_N_z']*proj['phantom_N_xy']*proj['phantom_N_xy']*4
			offset = offset + recon['zSlice4RMSE']*proj['phantom_N_xy']*proj['phantom_N_xy']*z_ratio*4
			offset = offset + proj['slice_t_start']*proj['phantom_N_xy']*proj['phantom_N_xy']*t_ratio*4
			fid.seek(offset, 0)
			data = fid.read(proj['phantom_N_xy']*proj['phantom_N_xy']*z_ratio*4)
		#	print j
			data = struct.unpack(str(proj['phantom_N_xy']*proj['phantom_N_xy']*z_ratio)+'f',data)
		#	print np.mean(data)
			data = np.reshape(data, (z_ratio, proj['phantom_N_xy'], proj['phantom_N_xy']), order='C')
		#	print np.mean(data)
			data = np.mean(data, axis=0, dtype=np.float32)
		#	print np.mean(data)
#			data = imresize(data, 0.5, interp='nearest', mode=None)
			data = decimate_count_data_in_r(data, proj['phantom_N_xy'], recon['N_xy'])/proj['phantom_N_xy']*recon['N_xy']
			data = decimate_count_data_in_t(data, proj['phantom_N_xy'], recon['N_xy'])/proj['phantom_N_xy']*recon['N_xy']
			phantom[index,:,:] = data
		#	print np.mean(phantom[index,:,:])
		#	print np.mean(interped_object)
			index = index + 1			

		fid.close()
#		write_Video4mArray(path2launch+'object_interpolated.mp4', recon['N_xy'], recon['N_xy'], recon['r'][i], interped_object)
#		write_Video4mArray(path2launch+'phantom.mp4', recon['N_xy'], recon['N_xy'], recon['r'][i], phantom)
		RMSE[i] = np.sqrt(np.mean((phantom-interped_object)*(phantom-interped_object)))
		RMSE[i] = convert_um2HU(RMSE[i])
		path2results = files['Result_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '_RMSE_' + str(int(RMSE[i])) + '/'
		create_folder(path2results)	
		for j in range(recon['ProjNumRMSE']):
			arr2tif = np.hstack((phantom[j,:,:], interped_object[j,:,:]))
			arr2tif = convert_um2HU (arr2tif)
			arr2tif = convert_HU2uint8(arr2tif, recon['maxHU'], recon['minHU'])
			write_array2tif (path2results + 'ObjVsPhantom_' + str(j) + '.tif', arr2tif)
		
		text_file = open(path2results + 'RMSE.txt', 'w')
		text_file.write('RMSE = ' + str(RMSE[i]) + '\n')
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
		text_file.close()
		fid = open(path2results + 'Interped_ObjSlice.bin', 'wb')
		interped_object.tofile(fid)		
		fid.close()	

	print RMSE
		#	a = f_interp(proj_times[512])
		#	print a[256,256]
		#	print interped_object	
	

