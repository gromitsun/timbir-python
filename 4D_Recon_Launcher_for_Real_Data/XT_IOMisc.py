
from PIL import Image
from os import system
import numpy as np
import h5py
import subprocess as sp
import struct

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
#	print 'Average value of image is ' + str(np.mean(array));
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
                path2launch = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
                path2results = files['Result_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
		for j in range(int(recon['Rtime_num'][i])):
			for k in range(recon['node_num']):
				zpernode = recon['N_z']/recon['node_num']	
				temp = RealData4mBin(path2launch + 'object_n' + str(k) + '_time_' + str(j) + '.bin', 0, recon['N_xy']*recon['N_xy']*zpernode, recon['real_var_type'])
				Object[k*zpernode:(k+1)*zpernode,:,:] = temp.reshape((zpernode, recon['N_xy'], recon['N_xy']), order = 'C')
			for k in range(recon['N_z']):
				slice_object = convert_um2HU(Object[k,:,:])
		#		print 'write_tiff_from_object_bin_file: Average value of time slice ' + str(j) + ' in HU is ' + str(np.mean(slice_object))
				write_array2tif(path2results + 'object_time_z_' + str(j) + '_' + str(k) + '.tif', convert_HU2uint8(slice_object, recon['maxHU'], recon['minHU']))
			print 'write_tiff_from_object_bin_file: Average value of time slice ' + str(j) + ' in um^-1 is ' + str(np.mean(Object))

def write_object2HDF (proj, recon, files):
	Object = np.zeros((recon['N_z'], recon['N_xy'], recon['N_xy']), dtype = np.float64, order = 'C')
        for i in range(len(recon['r'])):
                path2results = files['Result_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
                path2launch = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
		for k in range(recon['node_num']):
			temp = RealData4mBin(path2launch + 'proj_offset_n' + str(k) + '.bin', 0, proj['recon_N_r']*proj['recon_N_t']/recon['node_num'], recon['real_var_type'])
                	proj_offset = temp.reshape((proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), order = 'C')
                	offset[:,k*proj['recon_N_t']/recon['node_num']:(k+1)*proj['recon_N_t']/recon['node_num']] = proj_offset.astype(np.float32)
                for j in range(int(recon['Rtime_num'][i])):
                	file = h5py.File(path2results + 'object_time_' + str(j) + '.hdf5', 'w');
#               	dset = file.create_dataset('object', (recon['Rtime_num'][i], recon['N_z'], recon['N_xy'], recon['N_xy']), dtype=np.float32, chunks=True, compression='lzf');
                	dset = file.create_dataset('object', (recon['N_z'], recon['N_xy'], recon['N_xy']), dtype=np.float32, chunks=True);
			zpernode = recon['N_z']/recon['node_num']	
			for k in range(recon['node_num']):
				temp = RealData4mBin(path2launch + 'object_n' + str(k) + '_time_' + str(j) + '.bin', 0, recon['N_xy']*recon['N_xy']*zpernode, recon['real_var_type'])
                        	Object[k*zpernode:(k+1)*zpernode,:,:] = temp.reshape((zpernode, recon['N_xy'], recon['N_xy']), order = 'C')

                        Object[Object < 0] = 0
                        dset = Object.astype(np.float32)
                	if (j == 0):
				dset_off = file.create_dataset('proj_offset', (proj['recon_N_r'], proj['recon_N_t']), dtype=np.float32, chunks=True)
				dset_off = offset.astype(np.float32)
                	file.close

def write_Video4mArray(filename, rows, cols, fps, array):
	sz = array.shape
	pipe = sp.Popen(['ffmpeg',
        '-y', # (optional) overwrite the output file if it already exists
        '-f', 'rawvideo',
        '-vcodec','rawvideo',
        '-s', str(sz[1]) + 'x' + str(sz[2]), # size of one frame
        '-pix_fmt', 'rgb24',
        '-r', str(fps), # frames per second
        '-i', '-', # The input comes from a pipe
        '-an', # Tells FFMPEG not to expect any audio
        '-vcodec', 'mpeg',
        filename ],
        stdin=sp.PIPE,stdout=sp.PIPE, stderr=sp.PIPE)

	rgb_array[...,0] = array
	rgb_array[...,1] = array
	rgb_array[...,2] = array
	rgb_array.tofile(pipe.proc.stdin)


def RealData2Bin(filename, data, real_var_type):
	fid = open(filename, 'wb')
	if (real_var_type == 'double'):
		data = data.astype(np.float64)
	else:
		data = data.astype(np.float32)
	data.tofile(fid)
	fid.close()
		
def RealData4mBin(filename, offset, size, real_var_type):
	fid = open(filename, 'rb')
	if (real_var_type == 'double'):
		fid.seek(offset*8, 0)
		# Multiply above offset by 8 which is size of double in bytes
		data = fid.read(size*8)
		data = np.array(struct.unpack(str(size)+'d',data))	
	else:
		fid.seek(offset*4, 0)
		data = fid.read(size*4)
		data = np.array(struct.unpack(str(size)+'f',data))
	
	fid.close()
	return data
