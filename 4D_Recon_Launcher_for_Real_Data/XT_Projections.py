
from XT_IOMisc import error_by_flag
import numpy as np
import h5py
from XT_IOMisc import RealData4mBin
#from mpi4py import MPI

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
	


def generate_projections (proj, recon, files, path2launch):
	rank = recon['rank']
	if (recon['sinobin'] != 1 and recon['sinobin'] != 3):
		error_by_flag(1, "ERROR: sinobin should be either 1 or 3")

	if (proj['recon_N_t'] % recon['node_num'] != 0 or proj['N_t'] % recon['node_num'] != 0 or recon['N_z'] % recon['node_num'] != 0 ):
		error_by_flag(1, "ERROR: Number of nodes should divide N_t, recon_N_t and N_z")


	if (recon['sinobin'] == 1):
		projections = np.zeros((proj['recon_N_p'], proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), dtype = np.float64, order = 'C')
	else:	
		bright = np.zeros((proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), dtype = np.float64, order = 'C')

	weight = np.zeros((proj['recon_N_p'], proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), dtype = np.float64, order = 'C')
#	proj['expected'] = np.zeros((proj['recon_N_p'], proj['recon_N_r']), dtype = np.float64, order = 'C')
	
	print proj['Path2Dataset']
	FILE = h5py.File(proj['Path2Dataset'], 'r')
	
#	print proj['Path2WhiteDark']
#	FILE_wd = h5py.File(proj['Path2WhiteDark'], 'r')
#	dark_ptr = FILE_wd['/exchange/data_dark']
#	white_ptr = FILE_wd['/exchange/data_white']
	
	print proj['Path2Whites'] # Added by Yue
	print proj['Path2Darks'] # Added by Yue
	FILE_w = h5py.File(proj['Path2Whites'], 'r') # Added by Yue
	FILE_d = h5py.File(proj['Path2Darks'], 'r') # Added by Yue
	dark_ptr = FILE_d['/exchange/data_dark'] # Added by Yue
	white_ptr = FILE_w['/exchange/data_white'] # Added by Yue
	data_ptr = FILE['/exchange/data']

	if (dark_ptr.shape[-1] != proj['N_r'] or white_ptr.shape[-1] != proj['N_r'] or data_ptr.shape[-1] != proj['N_r']):
		error_by_flag(1, 'ERROR: generate_projections: Number of detector bins does not match input, dark = ' + str(dark_ptr.shape[-1]) + ', white = ' + str(white_ptr.shape[-1]) + ', data = ' + str(data_ptr.shape[-1]))
	
	if (rank == 0):
		print 'generate_projections: Number of projections in the HDF5 file is ' + str(data_ptr.shape[0])
	
	extras_r = proj['N_r'] % proj['recon_N_r']
	true_length_r = proj['N_r'] - extras_r
	x_start = proj['x_start'] # Added by Yue

	proj['length_r'] = proj['length_r']*true_length_r/proj['N_r']
	recon['radius_obj'] = recon['radius_obj']*true_length_r/proj['N_r']
	
#	index_r = range(extras_r/2, extras_r/2 + true_length_r)
	index_r = range(extras_r/2 + x_start, extras_r/2 + true_length_r + x_start) # Added by Yue
	index_t_start = proj['slice_t_start'] + rank*proj['N_t']/recon['node_num']
	index_t_end = proj['slice_t_start'] + (rank+1)*proj['N_t']/recon['node_num']
	max0idx = white_ptr.shape[0]
	if (proj['use_slice_white'] == -1):
		white = white_ptr[1:max0idx-1, index_t_start:index_t_end, index_r].astype(np.uint16)
		white = np.mean(white.astype(np.float64), axis=0)
	else:
		white = white_ptr[proj['use_slice_white'], index_t_start:index_t_end, index_r].astype(np.uint16)
	dark = dark_ptr[1:max0idx-1, index_t_start:index_t_end, index_r].astype(np.uint16)
	dark = np.mean(dark.astype(np.float64), axis=0)
	count_expected = decimate_count_data_in_r((np.abs(white - dark)).astype(np.float64), true_length_r, proj['recon_N_r'])
	count_expected = decimate_count_data_in_t(count_expected, proj['N_t']/recon['node_num'], proj['recon_N_t']/recon['node_num'])
	if (recon['sinobin'] != 1):
		bright = np.transpose(count_expected)
	for i in range(proj['recon_N_p']):
		data = data_ptr[proj['proj_start'] + i, index_t_start:index_t_end, index_r].astype(np.uint16)
		count_data = decimate_count_data_in_r((np.abs(data - dark)).astype(np.float64), true_length_r, proj['recon_N_r'])	
		count_data = decimate_count_data_in_t(count_data, proj['N_t']/recon['node_num'], proj['recon_N_t']/recon['node_num'])
		#count_data = decimate_count_data((np.abs(data[i,...])).astype(np.float64), true_length_r, proj['recon_N_r'])	
		weight[i,...] = np.transpose(count_data)
		if (recon['sinobin'] == 1):
			projections[i,...] = np.transpose(np.log(count_expected/count_data))
				
	if (recon['sinobin'] == 1):
		projections = recon['BH_Quad_Coef']*(projections*projections) + projections
		RealData2Bin(path2launch + 'projection_n' + str(rank) + '.bin', projections, recon['real_var_type'])
		print "Mean of projections from node " + str(rank) + " is ", np.mean(projections)
	else:
		RealData2Bin(path2launch + 'bright_n' + str(rank) + '.bin', bright, recon['real_var_type'])
		print "Mean of bright field data from node " + str(rank) + " is ", np.mean(bright)

	RealData2Bin(path2launch + 'weight_n' + str(rank) + '.bin', weight, recon['real_var_type'])
	print "Mean of weight data from node " + str(rank) + " is ", np.mean(weight)
	
	print "generate_projections: Generated projections for node ", rank
	
#	if (recon['sinobin'] == 1):
#		print 'generate_projections: The average value of projection is ', np.mean(proj['projections'])
#	else:	
#		print 'generate_projections: The average value of bright field is ', np.mean(proj['bright'])

#	print 'generate_projections: The average value of weight is ', np.mean(proj['weight'])

	FILE.close()
#	FILE_wd.close()
	FILE_w.close() # Added by Yue
	FILE_d.close() # Added by Yue
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
