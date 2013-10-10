
from scipy import interpolate
import numpy as np
from XT_IOMisc import error_by_flag
#import cv2,cv

#def write_array2video (proj,recon,files):
#	for i in range(len(recon['r'])):
#		path2results = files['Results_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_L_' + str(proj['L']) + '_K_' + str(proj['K']) + '_N_p_' + str(proj['N_p']) + '/'
#		cv2.VideoWriter.open(path2results+'object.avi', cv2.CV_FOURCC('M','J','P','G'), files['video_framerate'], (recon['N_xy'], recon['N_xy']), isColor=0) 
#			error_by_flag(1, 'ERROR: write_array2video: Cannot open VideoWriter object') 
#		for k in range(recon['N_p']):
#			frame = convert_mu2HU2uint8(recon['object'][k::])
#			cv2.VideoWriter.write(frame)
	
def interpolate_recon_in_1D (proj, recon, files):
	for i in range(len(recon['r'])):
		start = recon['Rtime0'] + recon['Rtime_delta'][i]/2
		stop = start + recon['Rtime_delta'][i]*(recon['Rtime_num'][i]-1)
		obj_times = np.linspace(start, stop, recon['Rtime_num'][i], endpoint=True)

		path2results = files['Results_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_L_' + str(proj['L']) + '_K_' + str(proj['K']) + '_N_p_' + str(proj['N_p']) + '/'
		fid = open(path2results + 'object.bin')
		obj = np.fromfile(fid, dtype=np.float64, count=-1)
		obj = np.reshape(obj, (recon['Rtime_num'][i], recon['N_xy'], recon['N_xy']), order='C')

		f_interp = interpolate.interp1d(obj_times, obj, kind='linear', axis=0, copy=False)
		proj_times = np.linspace(start, stop, proj['N_p'], endpoint=False)
	#	interped_object = f_interp(proj_times)
		a = f_interp(proj_times[512])
		print a[256,256]
	#	print interped_object	
	
		
