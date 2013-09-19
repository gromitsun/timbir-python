
from os import system
from XT_IOMisc import create_folder, error_by_flag

def write_views2file (path2save, views, times):
	fid = open(path2save + 'view_info.txt', 'w')

	for i in range(len(views)):
		fid.write(str(times[i]) + ' - ' + str(views[i]) + ',\n')

	fid.close()

def do_MBIR_reconstruction(proj, recon, files):
	path2source = files['C_Source_Folder']

######## CHANGE LAUNCH AND RESULT FOLDER USING VARIABLE BELOW #################
	launch_folder = files['Launch_Folder']
	result_folder = files['Result_Folder']
###############################################################################

	for i in range(len(recon['sigma_s'])):
		if (recon['use_same_folder_environ'] == 0):
			path2launch = launch_folder + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + '/'
			create_folder(path2launch)

			print ('do_reconstruction: Copying source files')
			flag = system('cp ' + path2source + '*.cpp ' + path2source + '*.sh ' + path2source + '*.c ' + path2source + '*.h ' + path2launch + '.')
			error_by_flag(flag, 'ERROR: cannot copy source files')
			
			if (files['copy_executables'] == 1):
				print 'do_reconstructions: Copying executables'
				flag = system('cp ' + path2source + 'XT_Engine ' +  path2launch + '.')
				error_by_flag(flag, 'ERROR: cannot copy executables')
			
			if (recon['updateProjOffset'][0] == 1 or recon['updateProjOffset'][0] == 3):
				print ('do_reconstruction: Copying proj_offset.bin')
				flag = system('cp ' + path2source + 'proj_offset*.bin ' + path2launch + '.')

			if (files['copy_projections'] == 1):
				print ('do_reconstruction: Copying projection.bin and weight.bin')
				flag = system('cp ' + path2source + 'bright*.bin ' + path2source + 'projection*.bin ' + path2source + 'weight*.bin ' + path2launch + '.')
				error_by_flag(flag, 'ERROR: cannot copy projection.bin/weight.bin')
			else:
				if (recon['sinobin'] == 1):
					for node in range(recon['node_num']):
						fid = open(path2launch + 'projection_n' + str(node) + '.bin','wb')
						proj['projections'][node,...].tofile(fid)
						fid.close()
				else:
					for node in range(recon['node_num']):
						fid = open(path2launch + 'bright_n' + str(node) + '.bin','wb')
						proj['bright'][node,...].tofile(fid)
						fid.close()
		
				for node in range(recon['node_num']):
					fid = open(path2launch + 'weight_n' + str(node) + '.bin','wb')
					proj['weight'][node,...].tofile(fid)
					fid.close()

			macros = ' -fopenmp -DBH_QUAD_COEF="' + str(recon['BH_Quad_Coef']) + '" -DHOUNSFIELD_MAX="' + str(recon['maxHU']) + '" -DHOUNSFIELD_MIN="' + str(recon['minHU']) + '"' 
			if (recon['calculate_cost'] == 0):
				macros = macros + ' -DNO_COST_CALCULATE'

			if (recon['NHICD'] == 0):
				macros = macros + ' -DNO_NHICD'
			
			if (recon['positivity_constraint'] == 1):
				macros = macros + ' -DPOSITIVITY_CONSTRAINT'

			if (files['copy_executables'] == 0):
				print('do_reconstruction: Compiling C code')
#        		if (engine.use_tifflibrary == 1)
#            			flag = system(['cd ',path2launch,';g++ -Wall -ansi', macros, ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c allocate.c TiffUtilities.cpp randlib.c tiff.c XT_IOMisc.c -lm -L/usr/local/tiff/lib -I/usr/local/tiff/include -ltiff '],'-echo');
				flag = system('cd ' + path2launch + ';mpicc -Wall -ansi' + macros + ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c XT_NHICD.c allocate.c randlib.c tiff.c XT_IOMisc.c XT_ImageProc.c XT_MPI.c -lm ')
				error_by_flag(flag, 'ERROR: Not able to compile')
				print 'do_reconstruction: Compile successful!'			

			write_views2file(path2launch, proj['angles'], proj['times'])
		
		path2launch = launch_folder + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '/'
		flag = system('cp nodefile ' + path2launch + '.')
		error_by_flag(flag, 'ERROR: Was not able to copy nodefile')
		for multidx in range(len(recon['delta_xy'])):
			ZingerT = recon['ZingerT'][i]
			if (recon['ZingerDel'][i] == 0 and multidx == 0):
				ZingerT = 100000 
			command = 'mpiexec -n ' + str(recon['node_num']) + ' -machinefile nodefile ./XT_Engine --p ' + str(recon['p']) + ' --sigma_s ' + str(recon['sigma_s'][i]) + ' --sigma_t ' + str(recon['sigma_t'][i]) + ' --c_s ' + str(recon['c_s'][i]) + ' --c_t ' + str(recon['c_t'][i]) + ' --delta_xy ' + str(recon['delta_xy'][multidx]) + ' --delta_z ' + str(recon['delta_z'][multidx]) + ' --length_r ' + str(proj['length_r']) + ' --length_t ' + str(proj['length_t']) + ' --voxel_thresh ' + str(recon['voxel_thresh'][multidx]) + ' --cost_thresh ' + str(recon['cost_thresh'][multidx]) + ' --iter ' + str(recon['iterations'][multidx]) + ' --rotation_center ' + str(proj['rotation_center_r']) + ' --alpha ' + str(recon['alpha']) + ' --sinobin ' + str(recon['sinobin']) + ' --initICD ' + str(recon['initICD'][multidx]) + ' --Rtime0 ' + str(recon['Rtime0']) + ' --Rtime_delta ' + str(recon['Rtime_delta'][i]) + ' --Rtime_num ' + str(recon['Rtime_num'][i]) + ' --num_projections ' + str(proj['recon_N_p']) + ' --N_r ' + str(proj['recon_N_r']) + ' --N_t ' + str(proj['recon_N_t']) + ' --detector_slice_begin ' + '0' + ' --detector_slice_end ' + str(proj['recon_N_t']-1) + ' --num_threads ' + str(recon['num_threads']) + ' --radius_obj ' + str(recon['radius_obj']) + ' --updateProjOffset ' + str(recon['updateProjOffset'][multidx]) + ' --writeTiff ' + str(recon['writeTiff'][multidx]) + ' --zingerT ' + str(ZingerT) + ' --zingerDel ' + str(recon['ZingerDel'][i])
			
			if (recon['updateProjOffset'][multidx] == 1 or recon['updateProjOffset'][multidx] == 3):
				if (multidx == 0 and recon['sinobin'] == 3):
					error_by_flag(1, 'ERROR: Cannot initialize offset error at coarsest resolution')
	
			if (recon['time_reg'] == 1):
				command = command + ' --time_reg'

			if (recon['WritePerIter'][multidx] == 1):
				command = command + ' --WritePerIter'
			
			if (recon['reconstruct'] == 0):
				command = command + ' --dont_reconstruct'
		
			if (recon['NHICD'] == 0):
				command = command + ' --no_NHICD'		
			
			if (recon['only_Edge_Updates'][multidx] == 1):
				command = command + ' --only_Edge_Updates'
	
			print 'Will run reconstruction code for delta_xy = ' + str(recon['delta_xy'][multidx]) + ' and delta_z = ' + str(recon['delta_z'][multidx])
			flag = system('cd ' + path2launch + ';' + command)
			error_by_flag(flag, 'ERROR: Was not able to run - ' + command)

			if (recon['reconstruct'] == 0):
				break
				
		path2results = result_folder + 'MBIR_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + '_zinger_' + str(recon['ZingerT'][i]) + '_' + str(recon['ZingerDel'][i]) + '/'
		create_folder(path2results)	
	
		flag = system('cp ' + path2launch + 'object_*tif ' + path2results + '.')
		error_by_flag(flag, 'ERROR: Could not move object_*.tif files')
					
		flag = system('cp ' + path2launch + 'object*.bin ' + path2results + '.')
		error_by_flag(flag, 'ERROR: Could not move object.bin files')
					
		flag = system('cp ' + path2launch + 'view_info.txt ' + path2results + '.')
		error_by_flag(flag, 'ERROR: Could not move view_info.txt')
				
		flag = system('cp ' + path2launch + 'DEBUG*.log ' + path2results + '.')
		error_by_flag(flag, 'ERROR: Could not move debug.log')
			
		if (any(recon['updateProjOffset'] > 1)): 
			flag = system('cp ' + path2launch + 'proj_offset* ' + path2results + '.')
			error_by_flag(flag, 'ERROR: Could not move proj_offset*')

		if (recon['calculate_cost'] == 1):
			flag = system('cp ' + path2launch + 'cost* ' + path2results + '.')
			error_by_flag(flag, 'ERROR: Could not move cost.bin')
