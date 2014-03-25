import numpy as np
import h5py


Path2Source = '../../Datasets/K_16_N_theta_2000_RotSpeed_100_Exp_4_ROI_1000x2080_Ramp_2/k-16-4ms-last_22.hdf'
Path2Dest = '../../Datasets/K_16_N_theta_2000_RotSpeed_100_Exp_4_ROI_1000x2080_Ramp_2/k-16-4ms-last_22_new.hdf'
first_slice_index = 500
num_slices_2_extract = 16
contain_bright_dark = 0

FILE_RD = h5py.File(Path2Source, 'r')
if (contain_bright_dark == 1):
	dark_rd = FILE_RD['/exchange/data_dark']
	white_rd = FILE_RD['/exchange/data_white']
	print dark_rd.shape
	print white_rd.shape
data_rd = FILE_RD['/exchange/data']
print data_rd.shape

FILE_WR = h5py.File(Path2Dest, 'w')
if (contain_bright_dark == 1):
	dark_wr = FILE_WR.create_dataset('/exchange/data_dark', (dark_rd.shape[0],num_slices_2_extract,dark_rd.shape[2]), dtype=np.uint16)
	white_wr = FILE_WR.create_dataset('/exchange/data_white', (white_rd.shape[0],num_slices_2_extract,white_rd.shape[2]), dtype=np.uint16)
data_wr = FILE_WR.create_dataset('/exchange/data', (data_rd.shape[0],num_slices_2_extract,data_rd.shape[2]), dtype=np.int16)

if (contain_bright_dark == 1):
	temp = white_rd[:,first_slice_index:first_slice_index + num_slices_2_extract,:]
	white_wr[:,:,:] = temp
	temp = dark_rd[:,first_slice_index:first_slice_index + num_slices_2_extract,:]
	dark_wr[:,:,:] = temp
	print np.mean(white_wr)
	print np.mean(dark_wr)

temp = data_rd[:,first_slice_index:first_slice_index + num_slices_2_extract,:]
data_wr[:,:,:] = temp
print np.mean(data_wr)

FILE_WR.close
FILE_RD.close


