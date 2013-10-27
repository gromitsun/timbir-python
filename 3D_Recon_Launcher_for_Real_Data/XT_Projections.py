""" ============================================================================
 * Copyright (c) 2013 K. Aditya Mohan (Purdue University)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of K. Aditya Mohan, Purdue
 * University, nor the names of its contributors may be used
 * to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """



from XT_IOMisc import error_by_flag
import numpy as np
import h5py
#from mpi4py import MPI

def calcAvgDeviationCount (weight, num_views, num_slices):
	dev = 0.0
	for i in range(num_views):
		for j in range(num_slices):
			data = weight[i,j,:]
			mean = np.mean(np.log(data[data > 0]))
			dev = dev + np.mean(np.absolute(np.log(data[data > 0]) - mean))
	dev = dev/(num_views*num_slices)
	return dev

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
	

def generate_projections_nersc (proj, recon, files, path2launch):
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
	
	FILE = h5py.File(proj['Path2Dataset'], 'r')
	FILE_wd = h5py.File(proj['Path2WhiteDark'], 'r')

	extras_r = proj['N_r'] % proj['recon_N_r']
	true_length_r = proj['N_r'] - extras_r

	proj['length_r'] = proj['length_r']*true_length_r/proj['N_r']
	recon['radius_obj'] = recon['radius_obj']*true_length_r/proj['N_r']
	
	index_r = range(extras_r/2, extras_r/2 + true_length_r)
	index_t_start = proj['slice_t_start'] + rank*proj['N_t']/recon['node_num']
	index_t_end = proj['slice_t_start'] + (rank+1)*proj['N_t']/recon['node_num']
	
        #bakdrk_idx = 0
	bakdrk_idx = proj['N_theta']
	white = FILE[proj['Dataset_Name'] + '/' + proj['Dataset_Name'] + 'bak_' + str(0).zfill(4) + '_' + str(bakdrk_idx).zfill(4)  + '.tif'][0, index_t_start:index_t_end, index_r].astype(np.uint16).astype(np.float64)
	dark = FILE[proj['Dataset_Name'] + '/' + proj['Dataset_Name'] + 'drk_' + str(0).zfill(4) + '_' + str(bakdrk_idx).zfill(4)  + '.tif'][0, index_t_start:index_t_end, index_r].astype(np.uint16).astype(np.float64)
	for i in range(1, proj['Num_Bright_Dark']):
		white = white + FILE[proj['Dataset_Name'] + '/' + proj['Dataset_Name'] + 'bak_' + str(i).zfill(4) + '_' + str(bakdrk_idx).zfill(4)  + '.tif'][0, index_t_start:index_t_end, index_r].astype(np.uint16).astype(np.float64)
		dark = dark + FILE[proj['Dataset_Name'] + '/' + proj['Dataset_Name'] + 'drk_' + str(i).zfill(4) + '_' + str(bakdrk_idx).zfill(4)  + '.tif'][0, index_t_start:index_t_end, index_r].astype(np.uint16).astype(np.float64)
	white = white/proj['Num_Bright_Dark']
	dark = dark/proj['Num_Bright_Dark']
	
	count_expected = decimate_count_data_in_r((np.abs(white - dark)).astype(np.float64), true_length_r, proj['recon_N_r'])
	count_expected = decimate_count_data_in_t(count_expected, proj['N_t']/recon['node_num'], proj['recon_N_t']/recon['node_num'])
	if (recon['sinobin'] != 1):
		bright = np.transpose(count_expected)
	for i in range(proj['recon_N_p']):
		data = FILE[proj['Dataset_Name'] + '/' + proj['Dataset_Name'] + '_0000_' + str(proj['proj_start'] + i*proj['view_subsmpl_fact']).zfill(4) + '.tif'][0, index_t_start:index_t_end, index_r].astype(np.uint16)
		count_data = decimate_count_data_in_r((np.abs(data - dark)).astype(np.float64), true_length_r, proj['recon_N_r'])	
		count_data = decimate_count_data_in_t(count_data, proj['N_t']/recon['node_num'], proj['recon_N_t']/recon['node_num'])
		#count_data = decimate_count_data((np.abs(data[i,...])).astype(np.float64), true_length_r, proj['recon_N_r'])
			
		weight[i,...] = np.transpose(count_data)
		if (recon['sinobin'] == 1):
			projections[i,...] = np.transpose(np.log(count_expected/count_data))
				
	if (recon['sinobin'] == 1):
		projections = recon['BH_Quad_Coef']*(projections*projections) + projections
		fid = open(path2launch + 'projection_n' + str(rank) + '.bin','wb')
		projections.tofile(fid)
		fid.close()
		print "Mean of projections from node " + str(rank) + " is ", np.mean(projections)
	else:
		fid = open(path2launch + 'bright_n' + str(rank) + '.bin','wb')
		bright.tofile(fid)
		fid.close()
		print "Mean of bright field data from node " + str(rank) + " is ", np.mean(bright)

	fid = open(path2launch + 'weight_n' + str(rank) + '.bin','wb')
	weight.tofile(fid)
	fid.close()
	print "Mean of weight data from node " + str(rank) + " is ", np.mean(weight)
	
	print "generate_projections: Generated projections for node ", rank
	dev = calcAvgDeviationCount (weight, proj['recon_N_p'], proj['recon_N_t']/recon['node_num'])/proj['length_r']
	
#	if (recon['sinobin'] == 1):
#		print 'generate_projections: The average value of projection is ', np.mean(proj['projections'])
#	else:	
#		print 'generate_projections: The average value of bright field is ', np.mean(proj['bright'])

#	print 'generate_projections: The average value of weight is ', np.mean(proj['weight'])

	FILE.close()
	FILE_wd.close()
	return dev


def generate_projections_purdue (proj, recon, files, path2launch):
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
	
	FILE = h5py.File(proj['Path2Dataset'], 'r')
	FILE_wd = h5py.File(proj['Path2WhiteDark'], 'r')
	dark_ptr = FILE_wd['/exchange/data_dark']
	white_ptr = FILE_wd['/exchange/data_white']
	data_ptr = FILE['/exchange/data']

	if (dark_ptr.shape[-1] != proj['N_r'] or white_ptr.shape[-1] != proj['N_r'] or data_ptr.shape[-1] != proj['N_r']):
		error_by_flag(1, 'ERROR: generate_projections: Number of detector bins does not match input, dark = ' + str(dark_ptr.shape[-1]) + ', white = ' + str(white_ptr.shape[-1]) + ', data = ' + str(data_ptr.shape[-1]))
	
	if (rank == 0):
		print 'generate_projections: Number of projections in the HDF5 file is ' + str(data_ptr.shape[0])
	
	extras_r = proj['N_r'] % proj['recon_N_r']
	true_length_r = proj['N_r'] - extras_r

	proj['length_r'] = proj['length_r']*true_length_r/proj['N_r']
	recon['radius_obj'] = recon['radius_obj']*true_length_r/proj['N_r']
	
	index_r = range(extras_r/2, extras_r/2 + true_length_r)
	index_t_start = proj['slice_t_start'] + rank*proj['N_t']/recon['node_num']
	index_t_end = proj['slice_t_start'] + (rank+1)*proj['N_t']/recon['node_num']
	max0idx = white_ptr.shape[0]
	white = white_ptr[1:max0idx-1, index_t_start:index_t_end, index_r].astype(np.uint16)
	white = np.mean(white.astype(np.float64), axis=0)
	dark = dark_ptr[1:max0idx-1:, index_t_start:index_t_end, index_r].astype(np.uint16)
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
		fid = open(path2launch + 'projection_n' + str(rank) + '.bin','wb')
		projections.tofile(fid)
		fid.close()
		print "Mean of projections from node " + str(rank) + " is ", np.mean(projections)
	else:
		fid = open(path2launch + 'bright_n' + str(rank) + '.bin','wb')
		bright.tofile(fid)
		fid.close()
		print "Mean of bright field data from node " + str(rank) + " is ", np.mean(bright)

	fid = open(path2launch + 'weight_n' + str(rank) + '.bin','wb')
	weight.tofile(fid)
	fid.close()
	
	print "Mean of weight data from node " + str(rank) + " is ", np.mean(weight)
	
	print "generate_projections: Generated projections for node ", rank
	
#	if (recon['sinobin'] == 1):
#		print 'generate_projections: The average value of projection is ', np.mean(proj['projections'])
#	else:	
#		print 'generate_projections: The average value of bright field is ', np.mean(proj['bright'])

#	print 'generate_projections: The average value of weight is ', np.mean(proj['weight'])

	FILE.close()
	FILE_wd.close()


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
