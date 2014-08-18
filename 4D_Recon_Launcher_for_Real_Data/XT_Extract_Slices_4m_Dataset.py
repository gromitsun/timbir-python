import numpy as np
import argparse
import h5py


parser = argparse.ArgumentParser()
parser.add_argument("--Path2Source", help="Full path of the input hdf5 file from which the slices are extracted", default="/dummy_path/")
parser.add_argument("--Path2Dest", help="Full path of the output hdf5 file to which the extracted slices are written", default="/dummy_path/")
parser.add_argument("--z_start_extract",help="Starting slice (along z-axis) for extraction",type=int)
parser.add_argument("--z_num_extract",help="Number slices (along z-axis) extracted",type=int)
parser.add_argument("--contain_whites",help="Input dataset contains whites", action="store_true")
parser.add_argument("--contain_darks",help="Input dataset contains darks", action="store_true")
args = parser.parse_args()

Path2Source = args.Path2Source
Path2Dest = args.Path2Dest
z_start = args.z_start_extract
z_num = args.z_num_extract

if (Path2Source == Path2Dest):
	print "ERROR: Source and Destination files are same"
else:
	FILE_RD = h5py.File(Path2Source, 'r')
	#if (args.contain_whites):
	#	white_rd = FILE_RD['/exchange/data_white']
	#	print white_rd.shape
	#if (args.contain_darks):
	#	dark_rd = FILE_RD['/exchange/data_dark']
	#	print dark_rd.shape
	#data_rd = FILE_RD['/exchange/data']
	#print data_rd.shape

	FILE_WR = h5py.File(Path2Dest, 'w')
	if (args.contain_whites):
		shape = FILE_RD['/exchange/data_white'].shape
		white_wr = FILE_WR.create_dataset('/exchange/data_white', (shape[0],z_num,shape[2]), dtype=np.uint16)

	if (args.contain_darks):
		shape = FILE_RD['/exchange/data_dark'].shape
		dark_wr = FILE_WR.create_dataset('/exchange/data_dark', (shape[0],z_num,shape[2]), dtype=np.uint16)
	
	shape = FILE_RD['/exchange/data'].shape
	data_wr = FILE_WR.create_dataset('/exchange/data', (shape[0],z_num,shape[2]), dtype=np.uint16)

	if (args.contain_whites):
		white_wr = FILE_RD['/exchange/data_white'][:,z_start:z_start + z_num,:]
	#	white_wr[:,:,:] = temp
		print 'Average of whites is ', np.mean(white_wr)
		print 'Copied whites!'
	if (args.contain_darks):
		dark_wr = FILE_RD['/exchange/data_dark'][:,z_start:z_start + z_num,:]
	#	dark_wr[:,:,:] = temp
		print 'Average of darks is ', np.mean(dark_wr)
		print 'Copied darks!'

	data_wr = FILE_RD['/exchange/data'][:,z_start:z_start + z_num,:]
	#data_wr[:,:,:] = temp
	print 'Average of the data is ', np.mean(data_wr)
	print 'Copied data. Done!'

	#neg = data_wr[3,:,:] - np.mean(dark_wr,axis=0)
	#print np.min(neg)
	#print np.sum(neg<0)
	#print np.var(np.mean(dark_wr,axis=0))
	
	#FILE_WR.close
	#FILE_RD.close


