
from PIL import Image
from os import system
import numpy as np
import h5py

def error_by_flag (flag, message):
	if (flag != 0):
		print 'ERROR_FLAG: ' + str(flag) 
		print message
		exit()

def create_folder (path2folder):
	print 'create_folder: Will create ' + path2folder 
	flag = system('mkdir ' + path2folder)
	if (flag != 0):
		print 'create_folder: Deleting ' + path2folder
		rm_flag = system('rm -r ' + path2folder)
		error_by_flag(rm_flag, 'ERROR: cannot remove ' + path2folder)

		print 'do_reconstruction: Creating ' + path2folder 
		flag = system('mkdir ' + path2folder)
		error_by_flag(flag, 'ERROR: cannot create ' + path2folder + '. Delete folder manually.')


def convert_um2HU (Object):
	density_scale = 0.001**4
	att_scale = 100**4

	AIR_MASS_ATT_COEFF = 0.496372335005353 * att_scale
	WATER_MASS_ATT_COEFF = 0.521225397034623 * att_scale

	WATER_DENSITY = 1.0 * density_scale
	AIR_DENSITY = 0.001205 * density_scale
	HOUNSFIELD_WATER_MAP = 1000
	HOUNSFIELD_AIR_MAP = 0
	slope = (HOUNSFIELD_WATER_MAP-HOUNSFIELD_AIR_MAP)/(WATER_MASS_ATT_COEFF*WATER_DENSITY-AIR_MASS_ATT_COEFF*AIR_DENSITY)
	c = -slope*(AIR_MASS_ATT_COEFF*AIR_DENSITY)
	Object = Object*slope + c
	return Object


def convert_HU2uint8 (Object, maxval, minval):
#	maxval = np.amax(Object)
#	minval = np.amin(Object)
	Object = (255.0/(maxval - minval))*(Object - minval);
	Object[Object <= 0] = 0
	Object[Object >= 255.0] = 255.0
	return Object.astype(np.uint8)

def write_array2tif (filename, array):
	print 'Average value of image is ' + str(np.mean(array));
	img = Image.fromarray(array)
	img.save(filename)

def read_slice_of_tiff_file (filename, slice_num):
	img_file = Image.open(filename)
	Nx,Ny = img_file.size
	img_crop = img_file.crop((0, slice_num, Nx, slice_num + 1))
	img_data = img_crop.getdata()
	img = np.asarray(img_data, dtype=np.uint16)
	return img

def write_tiff_from_object_bin_file (proj, recon, files):
	Object = np.zeros((recon['N_z'], recon['N_xy'], recon['N_xy']), dtype = np.float64, order = 'C')
	for i in range(len(recon['r'])):
                path2results = files['Result_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '_zinger_' + str(recon['ZingerT'][i]) + '_' + str(recon['ZingerDel'][i]) + '/'
		for j in range(int(recon['Rtime_num'][i])):
			for k in range(recon['node_num']):
				zpernode = recon['N_z']/recon['node_num']	
				Object[k*zpernode:(k+1)*zpernode,:,:] = np.fromfile(path2results + 'object_n' + str(k) + '_time_' + str(j) + '.bin', dtype = np.float64, count = -1).reshape((zpernode, recon['N_xy'], recon['N_xy']), order = 'C')
			for k in range(recon['N_z']):
				slice_object = convert_um2HU(Object[k,:,:])
		#		print 'write_tiff_from_object_bin_file: Average value of time slice ' + str(j) + ' in HU is ' + str(np.mean(slice_object))
				write_array2tif(path2results + 'object_z_' + str(k) + '_time_' + str(j) + '.tif', convert_HU2uint8(slice_object, recon['maxHU'], recon['minHU']))
			print 'write_tiff_from_object_bin_file: Average value of time slice ' + str(j) + ' in um^-1 is ' + str(np.mean(Object))

def write_object2HDF (proj, recon, files):
	Object = np.zeros((recon['N_z'], recon['N_xy'], recon['N_xy']), dtype = np.float64, order = 'C')
        for i in range(len(recon['r'])):
                path2results = files['Result_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '_zinger_' + str(recon['ZingerT'][i]) + '_' + str(recon['ZingerDel'][i]) + '/'
                file = h5py.File(path2results + 'object.hdf5', 'w');
#               dset = file.create_dataset('object', (recon['Rtime_num'][i], recon['N_z'], recon['N_xy'], recon['N_xy']), dtype=np.float32, chunks=True, compression='lzf');
                dset = file.create_dataset('object', (recon['Rtime_num'][i], recon['N_z'], recon['N_xy'], recon['N_xy']), dtype=np.float32, chunks=True);
                dset_off = file.create_dataset('proj_offset', (proj['recon_N_r'], proj['recon_N_t']), dtype=np.float32, chunks=True)
		for k in range(recon['node_num']):
                	proj_offset = np.fromfile(path2results + 'proj_offset_n' + str(k) + '.bin', dtype = np.float64, count = -1).reshape((proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), order = 'C')
                	dset_off[:,k*proj['recon_N_t']/recon['node_num']:(k+1)*proj['recon_N_t']/recon['node_num']] = proj_offset.astype(np.float32)
                for j in range(int(recon['Rtime_num'][i])):
			zpernode = recon['N_z']/recon['node_num']	
			for k in range(recon['node_num']):
                        	Object[k*zpernode:(k+1)*zpernode,:,:] = np.fromfile(path2results + 'object_n' + str(k) + '_time_' + str(j) + '.bin', dtype = np.float64, count = -1).reshape((zpernode, recon['N_xy'], recon['N_xy']), order = 'C')

                        Object[Object < 0] = 0
                        dset[j,:,:,:] = Object.astype(np.float32)
                file.close

