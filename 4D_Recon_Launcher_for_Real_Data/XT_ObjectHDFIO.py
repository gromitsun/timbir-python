import numpy as np
import h5py
from mpi4py import MPI
from XT_IOMisc import convert_um2HU
from XT_IOMisc import write_array2tif
from XT_IOMisc import convert_HU2uint8

def initpar_object4mHDF (proj, recon, files, i):
	rank = recon['rank']
	path2results = files['Result_Folder'] + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '_zinger_' + str(recon['ZingerT'][i]) + '_' + str(recon['ZingerDel'][i]) + '/'
	path2launch = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '/'
	file = h5py.File(path2results + 'object.hdf5', 'r', comm=MPI.COMM_WORLD);
	zpernode = proj['recon_N_t']/(recon['delta_z'][0]*recon['node_num'])
	
	proj_offset = file['/proj_offset'][:,rank*zpernode:(rank+1)*zpernode].astype(np.float64)
	fid = open(path2launch + 'proj_offset_n' + str(rank) + '.bin','wb')
	proj_offset.tofile(fid)
	fid.close()
	
	for j in range(int(recon['Rtime_num'][i])):
		object = file['/object'][j,rank*zpernode:(rank+1)*zpernode,:,:].astype(np.float64)
		fid = open(path2launch + 'object_n' + str(rank) + '_time_' + str(j) + '.bin','wb')
		object.tofile(fid)
		fid.close()
	file.close

	 		
def writepar_object2HDF (proj, recon, files):
	zpernode = recon['N_z']/recon['node_num']
	Object = np.zeros((zpernode, recon['N_xy'], recon['N_xy']), dtype = np.float64, order = 'C')
        for i in range(len(recon['r'])):
                path2results = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '/'
		recon_start = int(recon['rank']*recon['Rtime_num'][i]/recon['size'])
		recon_stop = int((recon['rank'] + 1)*recon['Rtime_num'][i]/recon['size'])
		for j in range(recon_start, recon_stop):
                	file = h5py.File(path2results + 'object_' + str(j) + '.hdf5', 'w')
#               	dset = file.create_dataset('object', (recon['Rtime_num'][i], recon['N_z'], recon['N_xy'], recon['N_xy']), dtype=np.float32, chunks=True, compression='lzf');
                	dset_obj = file.create_dataset('object', (recon['N_z'], recon['N_xy'], recon['N_xy']), dtype=np.float32, chunks=True)
                	if (j == 0):
				dset_off = file.create_dataset('proj_offset', (proj['recon_N_r'], proj['recon_N_t']), dtype=np.float32, chunks=True)
				for k in range(recon['node_num']):
                			proj_offset = np.fromfile(path2results + 'proj_offset_n' + str(k) + '.bin', dtype = np.float64, count = -1).reshape((proj['recon_N_r'], proj['recon_N_t']/recon['node_num']), order = 'C')
	                		dset_off[:,k*proj['recon_N_t']/recon['node_num']:(k+1)*proj['recon_N_t']/recon['node_num']] = proj_offset.astype(np.float32)
			for k in range(recon['node_num']):
                        	Object = np.fromfile(path2results + 'object_n' + str(k) + '_time_' + str(j) + '.bin', dtype = np.float64, count = -1).reshape((zpernode, recon['N_xy'], recon['N_xy']), order = 'C')

                        	Object[Object < 0] = 0
                        	dset_obj[k*zpernode:(k+1)*zpernode,:,:] = Object.astype(np.float32)

                	file.close


def writepar_tiff_from_object_bin_file (proj, recon, files):
	zpernode = recon['N_z']/recon['node_num']	
	Object = np.zeros((zpernode, recon['N_xy'], recon['N_xy']), dtype = np.float64, order = 'C')
	for i in range(len(recon['r'])):
                path2results = files['Launch_Folder'] + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '/'
		recon_start = int(recon['rank']*recon['Rtime_num'][i]/recon['size'])
		recon_stop = int((recon['rank'] + 1)*recon['Rtime_num'][i]/recon['size'])
		for j in range(recon_start, recon_stop):
			for k in range(recon['node_num']):
				Object = np.fromfile(path2results + 'object_n' + str(k) + '_time_' + str(j) + '.bin', dtype = np.float64, count = -1).reshape((zpernode, recon['N_xy'], recon['N_xy']), order = 'C')
				for l in range(zpernode):
					slice_object = convert_um2HU(Object[l,:,:])
		#			print 'write_tiff_from_object_bin_file: Average value of time slice ' + str(j) + ' in HU is ' + str(np.mean(slice_object))
					write_array2tif(path2results + 'object_z_' + str(k*zpernode+l) + '_time_' + str(j) + '.tif', convert_HU2uint8(slice_object, recon['maxHU'], recon['minHU']))

