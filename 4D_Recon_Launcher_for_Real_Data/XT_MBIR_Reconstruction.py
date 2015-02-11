
from os import system
from XT_IOMisc import create_folder, error_by_flag
from XT_Projections import generate_projections
from XT_ObjectHDFIO import initpar_object4mHDF

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
		path2launch = launch_folder + 'run_' + 'sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
		if (recon['set_up_launch_folder'] == 1):
			if (recon['rank'] == 0):
				print 'Setting up run folder by node with rank ', recon['rank']
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
				
				flag = system('cp ' + path2source + 'proj_offset_sim.bin ' + path2launch + '.')
			if (files['copy_projections'] == 1 and recon['rank'] == 0):
				if (recon['sinobin'] == 1):
					print ('do_reconstruction: Copying projection.bin and weight.bin')
					flag = system('cp ' + path2source + 'projection*.bin ' + path2source + 'weight*.bin ' + path2launch + '.')
				elif (recon['sinobin'] == 3):
					print ('do_reconstruction: Copying proj_offset_sim.bin, bright.bin and weight.bin')
					flag = system('cp ' + path2source + 'bright*.bin ' + path2source + 'weight*.bin ' + path2launch + '.')
				error_by_flag(flag, 'ERROR: cannot copy projection.bin/weight.bin')
			elif (recon['readSino4mHDF'][0] == 0 and recon['sinobin'] != 2):
				#if (recon['HPC'] == 'NERSC'):
				for j in range(recon['node_num']):
					recon['rank'] = j
					generate_projections(proj, recon, files, path2launch)
				recon['rank'] = 0
				#else:
				#	generate_projections (proj, recon, files, path2launch)
					

			macros = ' -DBH_QUAD_COEF="' + str(recon['BH_Quad_Coef']) + '" -DHOUNSFIELD_MAX="' + str(recon['maxHU']) + '" -DHOUNSFIELD_MIN="' + str(recon['minHU']) + '"' 
			#macros = ' -DBH_QUAD_COEF="' + str(recon['BH_Quad_Coef']) + '" -DHOUNSFIELD_MAX="' + str(recon['maxHU']) + '" -DHOUNSFIELD_MIN="' + str(recon['minHU']) + '"' 
	
			macros = macros + ' -DDATA_HDF_FILENAME="\\"' + proj['Path2Dataset'] + '\\""'
			macros = macros + ' -DWHITE_HDF_FILENAME="\\"' + proj['Path2Whites'] + '\\""'
			macros = macros + ' -DDARK_HDF_FILENAME="\\"' + proj['Path2Darks'] + '\\""'
			macros = macros + ' -DPROJECTION_HDF_START="' + str(proj['proj_start']) + '"'
			macros = macros + ' -DPATH_TO_PHANTOM="\\"' + proj['Path2Phantom'] + '\\""'
			macros = macros + ' -DEXPECTED_COUNTS_FOR_PHANTOM_DATA="' + str(proj['Expected_Counts']) + '"'
			macros = macros + ' -DCONVERGED_OBJECT_FILE="\\"' + files['Converged_Object'] + '\\""'
			if (recon['real_var_type'] == 'double'):
				macros = macros + ' -DREAL_IS_DOUBLE'

			if (recon['calculate_cost'] == 0):
				macros = macros + ' -DNO_COST_CALCULATE'

			if (recon['NHICD'] == 0):
				macros = macros + ' -DNO_NHICD'
			
			if (recon['positivity_constraint'] == 1):
				macros = macros + ' -DPOSITIVITY_CONSTRAINT'

#			if (any(recon['readSino4mHDF'])):
			macros = macros + ' -DREAD_PROJECTION_DATA_4M_HDF'
			if (recon['modality'] == 'PHCON'):
				macros = macros + ' -DPHASE_CONTRAST_TOMOGRAPHY'

			if (recon['rank'] == 0 and files['copy_executables'] == 0):
				print('do_reconstruction: Compiling C code')
#        		if (engine.use_tifflibrary == 1)
#            			flag = system(['cd ',path2launch,';g++ -Wall -ansi', macros, ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c allocate.c TiffUtilities.cpp randlib.c tiff.c XT_IOMisc.c -lm -L/usr/local/tiff/lib -I/usr/local/tiff/include -ltiff '],'-echo');
				print 'cd ' + path2launch + ';' + recon['compile_command'] + macros + ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c XT_NHICD.c allocate.c randlib.c tiff.c XT_IOMisc.c XT_ImageProc.c XT_MPI.c XT_HDFIO.c XT_VoxUpdate.c XT_ForwardProject.c XT_Filter.c -lm '
				flag = system('cd ' + path2launch + ';' + recon['compile_command'] + macros + ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c XT_NHICD.c allocate.c randlib.c tiff.c XT_IOMisc.c XT_ImageProc.c XT_MPI.c XT_HDFIO.c XT_VoxUpdate.c XT_ForwardProject.c XT_Filter.c -lm ')
				#print 'cd ' + path2launch + ';' + recon['compile_command'] + macros + ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c XT_NHICD.c allocate.c randlib.c tiff.c XT_IOMisc.c XT_ImageProc.c XT_MPI.c XT_HDFIO.c -lm '
				#flag = system('cd ' + path2launch + ';' + recon['compile_command'] + macros + ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c XT_NHICD.c allocate.c randlib.c tiff.c XT_IOMisc.c XT_ImageProc.c XT_MPI.c XT_HDFIO.c -lm ')
				error_by_flag(flag, 'ERROR: Not able to compile')
				print 'do_reconstruction: Compile successful!'			
			if (recon['rank'] == 0):
				write_views2file(path2launch, proj['angles'], proj['times'])
		
			if (recon['init_object4mHDF'] == 1):
				initpar_object4mHDF (proj, recon, files, i)

		if (recon['reconstruct'] == 1):
			if (recon['HPC'] == 'Purdue' and recon['node_num'] > 1):
				flag = system('cp nodefile ' + path2launch + '.')
				error_by_flag(flag, 'ERROR: Cannot copy nodefile to launch folder')
			elif (recon['HPC'] == 'Quest' and recon['node_num'] > 1): # Added by Yue
				flag = system('cp nodefile ' + path2launch + '.')
				error_by_flag(flag, 'ERROR: Cannot copy nodefile to launch folder')

			for multidx in range(recon['multstart'],len(recon['delta_xy'])):
				ZingerT = recon['ZingerT'][i]
				if (recon['ZingerDel'][i] == 0 and multidx == 0):
					ZingerT = 100000
		
				sinobin = recon['sinobin']
				if (multidx > 0 and recon['sinobin'] == 2):
					sinobin = 1
				 
#				command = recon['run_command'] + ' ./XT_Engine --p ' + str(recon['p']) + ' --sigma_s ' + str(recon['sigma_s'][i]) + ' --sigma_t ' + str(recon['sigma_t'][i]) + ' --c_s ' + str(recon['c_s'][i]) + ' --c_t ' + str(recon['c_t'][i]) + ' --delta_xy ' + str(recon['delta_xy'][multidx]) + ' --delta_z ' + str(recon['delta_z'][multidx]) + ' --length_r ' + str(proj['length_r']) + ' --length_t ' + str(proj['length_t']) + ' --voxel_thresh ' + str(recon['voxel_thresh'][multidx]) + ' --cost_thresh ' + str(recon['cost_thresh'][multidx]) + ' --iter ' + str(recon['iterations'][multidx]) + ' --rotation_center ' + str(proj['rotation_center_r']) + ' --alpha ' + str(recon['alpha']) + ' --sinobin ' + str(sinobin) + ' --initICD ' + str(recon['initICD'][multidx]) + ' --Rtime0 ' + str(recon['Rtime0']) + ' --Rtime_delta ' + str(recon['Rtime_delta'][i]) + ' --Rtime_num ' + str(recon['Rtime_num'][i]) + ' --num_projections ' + str(proj['recon_N_p']) + ' --N_r ' + str(proj['recon_N_r']) + ' --N_t ' + str(proj['recon_N_t']) + ' --detector_slice_begin ' + str(proj['slice_t_start']) + ' --detector_slice_num ' + str(proj['N_t']) + ' --num_threads ' + str(recon['num_threads']) + ' --radius_obj ' + str(recon['radius_obj']) + ' --updateProjOffset ' + str(recon['updateProjOffset'][multidx]) + ' --writeTiff ' + str(recon['writeTiff'][multidx]) + ' --zingerT ' + str(ZingerT) + ' --zingerDel ' + str(recon['ZingerDel'][i]) + ' --Est_of_Var ' + str(recon['Estimate_of_Var']) + ' --phantom_N_xy ' + str(proj['phantom_N_xy']) + ' --phantom_N_z ' + str(proj['phantom_N_z'])
				
				command = recon['run_command'] + ' ./XT_Engine --p ' + str(recon['p']) + ' --sigma_s ' + str(recon['sigma_s'][i]) + ' --sigma_t ' + str(recon['sigma_t'][i]) + ' --c_s ' + str(recon['c_s'][i]) + ' --c_t ' + str(recon['c_t'][i]) + ' --delta_xy ' + str(recon['delta_xy'][multidx]) + ' --delta_z ' + str(recon['delta_z'][multidx]) + ' --length_r ' + str(proj['length_r']) + ' --length_t ' + str(proj['length_t']) + ' --voxel_thresh ' + str(recon['voxel_thresh'][multidx]) + ' --cost_thresh ' + str(recon['cost_thresh'][multidx]) + ' --iter ' + str(recon['iterations'][multidx]) + ' --rotation_center ' + str(proj['rotation_center_r']) + ' --alpha ' + str(recon['alpha']) + ' --sinobin ' + str(sinobin) + ' --initICD ' + str(recon['initICD'][multidx]) + ' --Rtime0 ' + str(recon['Rtime0']) + ' --Rtime_delta ' + str(recon['Rtime_delta'][i]) + ' --Rtime_num ' + str(recon['Rtime_num'][i]) + ' --num_projections ' + str(proj['recon_N_p']) + ' --x_start ' + str(proj['x_start']) + ' --x_width ' + str(proj['x_width']) + ' --N_r ' + str(proj['recon_N_r']) + ' --N_t ' + str(proj['recon_N_t']) + ' --detector_slice_begin ' + str(proj['slice_t_start']) + ' --detector_slice_num ' + str(proj['N_t']) + ' --num_threads ' + str(recon['num_threads']) + ' --radius_obj ' + str(recon['radius_obj']) + ' --updateProjOffset ' + str(recon['updateProjOffset'][multidx]) + ' --writeTiff ' + str(recon['writeTiff'][multidx]) + ' --zingerT ' + str(ZingerT) + ' --zingerDel ' + str(recon['ZingerDel'][i]) + ' --Est_of_Var ' + str(recon['Estimate_of_Var']) + ' --phantom_N_xy ' + str(proj['phantom_N_xy']) + ' --phantom_N_z ' + str(proj['phantom_N_z']) # Added by Yue

				if (recon['zero_mean_offset'][multidx] == 1):
					command = command + ' --EnforceZeroMeanOffset'
					
				if (recon['time_reg'] == 1):
					command = command + ' --time_reg'

				if (recon['RMSE_converged'][multidx] == 1):
					command = command + ' --RMSE_converged'			
				if (recon['WritePerIter'][multidx] == 1):
					command = command + ' --WritePerIter'
				
				if (recon['NHICD'] == 0):
					command = command + ' --no_NHICD'		
			
				if (recon['only_Edge_Updates'][multidx] == 1):
					command = command + ' --only_Edge_Updates'
	
				if (recon['initMagUpMap'][multidx] == 1):
					command = command + ' --initMagUpMap'
				
				if (recon['do_VarEstimate'][multidx] == 1):
					command = command + ' --do_VarEst '
	
				if (recon['readSino4mHDF'][multidx] == 1):
					print 'Generating projection, weight, bright field data from HDF file'
					flag = system('cd ' + path2launch + ';' + command + ' --dont_reconstruct --readSino4mHDF')
					error_by_flag(flag, 'ERROR: Was not able to run - ' + command + ' --dont_reconstruct --readSino4mHDF')
				
				if (recon['recon_type'] == 'FBP'):
					command = command + ' --dont_reconstruct'
				print 'Will run reconstruction code for delta_xy = ' + str(recon['delta_xy'][multidx]) + ' and delta_z = ' + str(recon['delta_z'][multidx])
				flag = system('cd ' + path2launch + ';' + command)
				error_by_flag(flag, 'ERROR: Was not able to run - ' + command)

				if (recon['recon_type'] == 'FBP'):
					break
		
				
		path2results = result_folder + recon['recon_type'] + '_sigs_' + str(recon['sigma_s'][i]) + '_sigt_' + str(recon['sigma_t'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta']) + '_N_p_' + str(proj['recon_N_p']) + recon['msg_string'] + '/'
		create_folder(path2results)	
		
		if (recon['recon_type'] == 'FBP'):
			break
	
			#flag = system('cp ' + path2launch + 'object_*tif ' + path2results + '.')
			#error_by_flag(flag, 'ERROR: Could not move object_*.tif files')
					
		#flag = system('cp ' + path2launch + 'object*.bin ' + path2results + '.')
		#error_by_flag(flag, 'ERROR: Could not move object.bin files')
					
		flag = system('cp ' + path2launch + 'view_info.txt ' + path2results + '.')
		error_by_flag(flag, 'ERROR: Could not move view_info.txt')
				
		flag = system('cp ' + path2launch + 'DEBUG*.log ' + path2results + '.')
		error_by_flag(flag, 'ERROR: Could not move debug.log')
			
#		if (any(recon['updateProjOffset'] > 1)): 
#			flag = system('cp ' + path2launch + 'proj_offset* ' + path2results + '.')
#			error_by_flag(flag, 'ERROR: Could not move proj_offset*')

#		if (recon['calculate_cost'] == 1):
#			flag = system('cp ' + path2launch + 'cost* ' + path2results + '.')
#			error_by_flag(flag, 'ERROR: Could not move cost.bin')
