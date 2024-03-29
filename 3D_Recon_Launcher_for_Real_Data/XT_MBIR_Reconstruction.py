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


#! /usr/bin/python

from os import system
from XT_IOMisc import create_folder, error_by_flag
from XT_Projections import generate_projections_nersc
from XT_Projections import generate_projections_purdue
from XT_ObjectHDFIO import initpar_object4mHDF
from XT_IOMisc import convert_um2HU
import numpy as np

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
		path2launch = launch_folder + 'run_' + 'smooth_' + str(recon['smoothness'][i]) + '_slice_start_' + str(proj['slice_t_start']) + '_slice_num_' + str(proj['recon_N_t']) + '_T_' + str(recon['ZingerT'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + '/'
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

			if (files['copy_projections'] == 1 and recon['rank'] == 0):
				print ('do_reconstruction: Copying projection.bin and weight.bin')
				flag = system('cp ' + path2source + 'bright*.bin ' + path2source + 'projection*.bin ' + path2source + 'weight*.bin ' + path2launch + '.')
				error_by_flag(flag, 'ERROR: cannot copy projection.bin/weight.bin')
			else:
				dev = 0.0
			#	if (recon['HPC'] == 'NERSC'):
				for j in range(recon['node_num']):
					recon['rank'] = j
					dev = dev + generate_projections_nersc (proj, recon, files, path2launch)
				dev = dev/recon['node_num']			
				sigma_s = convert_um2HU(dev*np.array(recon['sigma_s']))
				#print convert_um2HU(dev)
				print recon['sigma_s']
				#print 'The regularization parameter is sigma_s = ' + sigma_s
				recon['rank'] = 0
			#	else:
			#		generate_projections_nersc(proj, recon, files, path2launch)
			#		sigma_s = recon['sigma_s']
					

			macros = ' -DBH_QUAD_COEF="' + str(recon['BH_Quad_Coef']) + '" -DHOUNSFIELD_MAX="' + str(recon['maxHU']) + '" -DHOUNSFIELD_MIN="' + str(recon['minHU']) + '"' 
                        macros = macros + ' -DDATA_HDF_FILENAME="\\"' + str(recon['BH_Quad_Coef']) + '\\""'
			macros = macros + ' -DWHITEDARK_HDF_FILENAME="\\"' + str(recon['BH_Quad_Coef']) + '\\""'
			macros = macros + ' -DPROJECTION_HDF_START="' + str(recon['BH_Quad_Coef']) + '"'
			macros = macros + ' -DPATH_TO_PHANTOM="\\"' + str(recon['BH_Quad_Coef']) + '\\""'
			macros = macros + ' -DEXPECTED_COUNTS_FOR_PHANTOM_DATA="' + str(recon['BH_Quad_Coef']) + '"'
			macros = macros + ' -DCONVERGED_OBJECT_FILE="\\"' + str(recon['BH_Quad_Coef']) + '\\""'
			if (recon['real_var_type'] == 'double'):
				macros = macros + ' -DREAL_IS_DOUBLE'
			
			if (recon['calculate_cost'] == 0):
				macros = macros + ' -DNO_COST_CALCULATE'

			if (recon['NHICD'] == 0):
				macros = macros + ' -DNO_NHICD'
			
			if (recon['positivity_constraint'] == 1):
				macros = macros + ' -DPOSITIVITY_CONSTRAINT'

			if (recon['rank'] == 0 and files['copy_executables'] == 0):
				print('do_reconstruction: Compiling C code')
#        		if (engine.use_tifflibrary == 1)
#            			flag = system(['cd ',path2launch,';g++ -Wall -ansi', macros, ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c allocate.c TiffUtilities.cpp randlib.c tiff.c XT_IOMisc.c -lm -L/usr/local/tiff/lib -I/usr/local/tiff/include -ltiff '],'-echo');
				flag = system('cd ' + path2launch + ';' + recon['compile_command'] + ' -Wall -ansi' + macros + ' -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c XT_NHICD.c allocate.c randlib.c tiff.c XT_IOMisc.c XT_ImageProc.c XT_MPI.c XT_VoxUpdate.c XT_ForwardProject.c XT_Filter.c -lm ')
				error_by_flag(flag, 'ERROR: Not able to compile')
				print 'do_reconstruction: Compile successful!'			
			if (recon['rank'] == 0):
				write_views2file(path2launch, proj['angles'], proj['times'])
		
			if (recon['init_object4mHDF'] == 1):
				initpar_object4mHDF (proj, recon, files, i)

		if (recon['reconstruct'] == 1):	
			if (recon['HPC'] == 'Purdue'):
				flag = system('cp nodefile ' + path2launch + '.')
				error_by_flag(flag, 'ERROR: Cannot copy nodefile to launch folder')
			for multidx in range(len(recon['delta_xy'])):
				ZingerT = recon['ZingerT'][i]
				if (recon['ZingerDel'][i] == 0 and multidx == 0):
					ZingerT = 100000 
				command = recon['run_command'] + ' ./XT_Engine --p ' + str(recon['p']) + ' --sigma_s ' + str(sigma_s[i]) + ' --sigma_t ' + str(recon['sigma_t'][i]) + ' --c_s ' + str(recon['c_s'][i]) + ' --c_t ' + str(recon['c_t'][i]) + ' --delta_xy ' + str(recon['delta_xy'][multidx]) + ' --delta_z ' + str(recon['delta_z'][multidx]) + ' --length_r ' + str(proj['length_r']) + ' --length_t ' + str(proj['length_t']) + ' --voxel_thresh ' + str(recon['voxel_thresh'][multidx]) + ' --cost_thresh ' + str(recon['cost_thresh'][multidx]) + ' --iter ' + str(recon['iterations'][multidx]) + ' --rotation_center ' + str(proj['rotation_center_r']) + ' --alpha ' + str(recon['alpha']) + ' --sinobin ' + str(recon['sinobin']) + ' --initICD ' + str(recon['initICD'][multidx]) + ' --Rtime0 ' + str(recon['Rtime0']) + ' --Rtime_delta ' + str(recon['Rtime_delta'][i]) + ' --Rtime_num ' + str(recon['Rtime_num'][i]) + ' --num_projections ' + str(proj['recon_N_p']) + ' --N_r ' + str(proj['recon_N_r']) + ' --N_t ' + str(proj['recon_N_t']) + ' --detector_slice_begin ' + '0' + ' --detector_slice_num ' + str(proj['recon_N_t']) + ' --num_threads ' + str(recon['num_threads']) + ' --radius_obj ' + str(recon['radius_obj']) + ' --updateProjOffset ' + str(recon['updateProjOffset'][multidx]) + ' --writeTiff ' + str(recon['writeTiff'][multidx]) + ' --zingerT ' + str(ZingerT) + ' --zingerDel ' + str(recon['ZingerDel'][i]) + ' --Est_of_Var ' + str(recon['Est_of_Var']) + ' --phantom_N_xy ' + str(0) + ' --phantom_N_z ' + str(0)	
			
				if (recon['time_reg'] == 1):
					command = command + ' --time_reg'
	
				if (recon['WritePerIter'][multidx] == 1):
					command = command + ' --WritePerIter'
				
				if (recon['reconstruct'] == 0):
					command = command + ' --dont_reconstruct'
			
				if (recon['NHICD'] == 0):
					command = command + ' --no_NHICD'		
				
				if (recon['do_VarEstimate'][multidx] == 1):
					command = command + ' --do_VarEst '
			
				if (recon['only_Edge_Updates'][multidx] == 1):
					command = command + ' --only_Edge_Updates'
	
				if (recon['initMagUpMap'][multidx] == 1):
					command = command + ' --initMagUpMap'
				
				print 'Will run reconstruction code for delta_xy = ' + str(recon['delta_xy'][multidx]) + ' and delta_z = ' + str(recon['delta_z'][multidx])
				flag = system('cd ' + path2launch + ';' + command)
				error_by_flag(flag, 'ERROR: Was not able to run - ' + command)

				if (recon['reconstruct'] == 0):
					break
				
		path2results = result_folder + 'MBIR_' + 'smooth_' + str(recon['smoothness'][i]) + '_slice_start_' + str(proj['slice_t_start']) + '_slice_num_' + str(proj['recon_N_t']) + '_T_' + str(recon['ZingerT'][i]) + '_r_' + str(recon['r'][i]) + '_K_' + str(proj['K']) + '_N_theta_' + str(proj['N_theta'])  + '_N_p_' + str(proj['recon_N_p']) + '/'
		create_folder(path2results)	

	        #Copy the output tiffs to the result folder
		flag = system('cp ' + path2launch + 'object_*tif ' + path2results + '.')
		error_by_flag(flag, 'ERROR: Could not move object_*.tif files')

		#Copy the output bins to the result folder
		flag = system('cp ' + path2launch + 'object*.bin ' + path2results + '.')
		error_by_flag(flag, 'ERROR: Could not move object.bin files')
					
		#flag = system('cp ' + path2launch + 'view_info.txt ' + path2results + '.')
		#error_by_flag(flag, 'ERROR: Could not move view_info.txt')
				
		#flag = system('cp ' + path2launch + 'DEBUG*.log ' + path2results + '.')
		#error_by_flag(flag, 'ERROR: Could not move debug.log')
			
#		if (any(recon['updateProjOffset'] > 1)): 
#			flag = system('cp ' + path2launch + 'proj_offset* ' + path2results + '.')
#			error_by_flag(flag, 'ERROR: Could not move proj_offset*')

#		if (recon['calculate_cost'] == 1):
#			flag = system('cp ' + path2launch + 'cost* ' + path2results + '.')
#			error_by_flag(flag, 'ERROR: Could not move cost.bin')
