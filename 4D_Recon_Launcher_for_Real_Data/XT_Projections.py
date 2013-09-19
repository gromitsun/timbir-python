
from XT_IOMisc import error_by_flag
import numpy as np
import h5py

def decimate_count_data_in_r (data, true_length, reduced_length):
	ratio = true_length/reduced_length
	shape = data.shape
	new_data = np.zeros(shape[0:len(shape)-1] + (reduced_length,), dtype = np.float64, order = 'C')

	for i in range(reduced_length):
		new_data[:,i] = np.sum(data[:,i*ratio:(i + 1)*ratio],1)
	return new_data

def decimate_count_data_in_t (data, true_length, reduced_length):
	ratio = true_length/reduced_length
	shape = data.shape
	new_data = np.zeros((reduced_length,)+shape[1:len(shape)], dtype = np.float64, order = 'C')

	for i in range(reduced_length):
		new_data[i,:] = np.sum(data[i*ratio:(i + 1)*ratio,:],0)
	return new_data
	


def generate_projections (proj, recon, files):
	if (recon['sinobin'] != 1 and recon['sinobin'] != 3):
		error_by_flag(1, "ERROR: sinobin should be either 1 or 3")

	if (recon['sinobin'] == 1):
		proj['projections'] = np.zeros((recon['node_num'], proj['recon_N_p'], proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), dtype = np.float64, order = 'C')
	else:	
		proj['bright'] = np.zeros((recon['node_num'], proj['recon_N_p'], proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), dtype = np.float64, order = 'C')

	proj['weight'] = np.zeros((recon['node_num'], proj['recon_N_p'], proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), dtype = np.float64, order = 'C')
#	proj['expected'] = np.zeros((proj['recon_N_p'], proj['recon_N_r']), dtype = np.float64, order = 'C')
	
	FILE = h5py.File(proj['Path2Dataset'], 'r')
	FILE_wd = h5py.File(proj['Path2WhiteDark'], 'r')
	dark_ptr = FILE_wd['/exchange/data_dark']
	white_ptr = FILE_wd['/exchange/data_white']
	data_ptr = FILE['/exchange/data']

	if (dark_ptr.shape[-1] != proj['N_r'] or white_ptr.shape[-1] != proj['N_r'] or data_ptr.shape[-1] != proj['N_r']):
		error_by_flag(1, 'ERROR: generate_projections: Number of detector bins does not match input, dark = ' + str(dark_ptr.shape[-1]) + ', white = ' + str(white_ptr.shape[-1]) + ', data = ' + str(data_ptr.shape[-1]))

	print 'generate_projections: Number of projections in the HDF5 file is ' + str(data_ptr.shape[0])
	
	extras_r = proj['N_r'] % proj['recon_N_r']
	true_length_r = proj['N_r'] - extras_r

	proj['length_r'] = proj['length_r']*true_length_r/proj['N_r']
	recon['radius_obj'] = recon['radius_obj']*true_length_r/proj['N_r']

	index_r = range(extras_r/2, extras_r/2 + true_length_r)
#	index_t_start = int(proj['N_t']/2) - int(proj['recon_N_t']/2)
	index_t_start = proj['slice_t_start']
	index_t_end = index_t_start + proj['N_t']
	white = white_ptr[3:, index_t_start:index_t_end, index_r].astype(np.uint16)
	white = np.mean(white.astype(np.float64), axis=0)
	dark = dark_ptr[3:, index_t_start:index_t_end, index_r].astype(np.uint16)
	dark = np.mean(dark.astype(np.float64), axis=0)
	
	count_expected = decimate_count_data_in_r((np.abs(white - dark)).astype(np.float64), true_length_r, proj['recon_N_r'])
	count_expected = decimate_count_data_in_t(count_expected, proj['N_t'], proj['recon_N_t'])
	#count_expected = decimate_count_data((np.abs(white)).astype(np.float64), true_length_r, proj['recon_N_r'])
	if (recon['use_same_folder_environ'] == 0):
		for i in range(proj['recon_N_p']):
			data = data_ptr[proj['proj_start'] + i, index_t_start:index_t_end, index_r].astype(np.uint16)
			count_data = decimate_count_data_in_r((np.abs(data - dark)).astype(np.float64), true_length_r, proj['recon_N_r'])	
			count_data = decimate_count_data_in_t(count_data, proj['N_t'], proj['recon_N_t'])
			#count_data = decimate_count_data((np.abs(data[i,...])).astype(np.float64), true_length_r, proj['recon_N_r'])	
			for j in range(recon['node_num']):
				zpernode = proj['recon_N_t']/recon['node_num']
				proj['weight'][j,i,...] = np.transpose(count_data[j*zpernode:(j+1)*zpernode,:])
				if (recon['sinobin'] == 1):
					proj['projections'][j,i,...] = np.transpose(np.log(count_expected[j*zpernode:(j+1)*zpernode,:]/count_data[j*zpernode:(j+1)*zpernode,:]))
				else:
					proj['bright'][j,i,...] = np.transpose(count_expected[j*zpernode:(j+1)*zpernode,:])
				
	if (recon['sinobin'] == 1):
		proj['projections'] = recon['BH_Quad_Coef']*(proj['projections']*proj['projections']) + proj['projections'];

	if (recon['sinobin'] == 1):
		print 'generate_projections: The average value of projection is ', np.mean(proj['projections'])
	else:	
		print 'generate_projections: The average value of bright field is ', np.mean(proj['bright'])

	print 'generate_projections: The average value of weight is ', np.mean(proj['weight'])

	FILE.close()
	FILE_wd.close()
	return proj


"""
def generate_projections (proj):
	error_by_flag(proj['tiff_cols'] % proj['N_r'], 'generate_projections: \'N_r\' should divide \'tiff_cols\'')	

	proj['projections'] = np.zeros((proj['N_p'], proj['N_r']), dtype = np.float64, order = 'C')
	proj['weight'] = np.zeros((proj['N_p'], proj['N_r']), dtype = np.float64, order = 'C')
	proj['expected'] = np.zeros((proj['N_p'], proj['N_r']), dtype = np.float64, order = 'C')
	idx = 0
	for i in range(proj['N_p']/proj['N_theta']):
		print 'generate_projections: Reading frame ' + str(i)
#		frame_number = 10 + i
#		frame_number = str(frame_number).zfill(4)

#		dark_1_filename = "/Volumes/LaCie/Data/26B_" + frame_number + "_/tif/26B_" + frame_number + "_" + '1'.zfill(4) + ".tif"
#		dark_1 = read_slice_of_tiff_file (dark_1_filename, proj['slice'])
#		dark_2_filename = "/Volumes/LaCie/Data/26B_" + frame_number + "_/tif/26B_" + frame_number + "_" + '2'.zfill(4) + ".tif"
#		dark_2 = read_slice_of_tiff_file (dark_2_filename, proj['slice'])
#		dark = (dark_1 + dark_2)/2
		
#		bright_1_filename = "/Volumes/LaCie/Data/26B_" + frame_number + "_/tif/26B_" + frame_number + "_" + '3'.zfill(4) + ".tif"
#		bright_1 = read_slice_of_tiff_file (bright_1_filename, proj['slice'])
#		bright_2_filename = "/Volumes/LaCie/Data/26B_" + frame_number + "_/tif/26B_" + frame_number + "_" + '4'.zfill(4) + ".tif"
#		bright_2 = read_slice_of_tiff_file (bright_2_filename, proj['slice'])
#		bright = (bright_1 + bright_2)/2

		bright = read_slice_of_tiff_file ("../../Datasets/Dataset_K_16_N_theta_256/Flat_Field.tif", proj['slice'])
		
		for j in range(proj['N_theta']):
#			view_number = 5 + j
#			view_number = str(view_number).zfill(4)
#			count_filename = "/Volumes/LaCie/Data/26B_" + frame_number + "_/tif/26B_" + frame_number + "_" + view_number + ".tif"
#			print count_filename
#			img = read_slice_of_tiff_file(count_filename, proj['slice'])

			img = read_slice_of_tiff_file("../../Datasets/Dataset_K_16_N_theta_256/View_" + str(i*proj['N_theta'] + j + 1) + ".tif", proj['slice'])
#			count_expected = np.abs(bright - dark).astype(np.float64)
#			count_data = np.abs(img - dark).astype(np.float64)
			count_expected = decimate_count_data((np.abs(bright)).astype(np.float64), proj['tiff_cols'], proj['N_r'])
			count_data = decimate_count_data((np.abs(img)).astype(np.float64), proj['tiff_cols'], proj['N_r'])
			if np.any(count_data <= 0.0) or np.any(count_expected <= 0.0):
				error_by_flag(1, "ERROR: generate_projections: count_data or count_expected has a zero element")
				
			proj['weight'][idx,:] = count_data
			proj['projections'][idx,:] = np.log(count_expected/count_data)
 			idx += 1

	return proj
"""
