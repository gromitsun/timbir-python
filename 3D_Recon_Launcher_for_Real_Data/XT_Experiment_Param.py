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



import numpy as np
from XT_Interlaced_Angles import gen_interlaced_views_0_to_Inf,clip_list_of_views

""" 	K - Number of sub-frames in a full frame
	N_theta - Number of projections in a full frame (Total number of discrete angles)
	fps - Frame rate of the camera (Number of projections per second)
	N_p - Total number of projections (used only to compute the deleted angles)
	min_time_btw_views - Minimum time between views. Limit imposed by the trigger system is 5ms  """

def initialize ():
	proj = {}

	proj['K'] = np.array([16])
	proj['N_theta'] = np.array([1600])
	proj['rotation_speed'] = np.array([60])
	N_p = 1600*3
	min_time_btw_views = 0.028
	
	proj['L'] = proj['N_theta']/proj['K']
	proj['fps'] = proj['L']/(180.0/proj['rotation_speed'])

	proj['num_angle_del'] = np.zeros(len(proj['K']), dtype = np.int32, order='C')
	for i in range(len(proj['K'])):
		angles, times = gen_interlaced_views_0_to_Inf(proj['K'][i], proj['N_theta'][i], N_p)
		angles_clip, times_clip, angles_del, times_del = clip_list_of_views (angles, times, min_time_btw_views, proj['rotation_speed'][i])
		proj['num_angle_del'][i] = angles_del.size 
	
	return proj

def gen_experiment_params (proj):
	proj['temporal_res'] = proj['L']/proj['fps']
	proj['delta_theta'] = 180.0/proj['N_theta']

	proj['max_exposure_time'] = 1.0/(proj['fps']*proj['K'])
	proj['time_delta_theta'] = proj['delta_theta']/proj['rotation_speed']
	
def write_params_to_console (proj):
	for i in range(len(proj['K'])):
		print 'K =', proj['K'][i], ', N_theta =', proj['N_theta'][i], ', L =', proj['L'][i], ', Temporal Resolution =', proj['temporal_res'][i], ', delta Theta =', proj['delta_theta'][i], ', Rotation Speed =', proj['rotation_speed'][i], ', Max Exposure Time (Time to traverse delta Theta) =', proj['max_exposure_time'][i], ', Number of angles deleted =', proj['num_angle_del'][i], ' \n '

def write_params_to_txt_file (proj):
	fid = open('experiment_parameters.txt','w')
	space = 15
	fid.write('K'.ljust(space) + 'N_theta'.ljust(space) + 'L'.ljust(space) + 'Time Res'.ljust(space) + 'del Theta'.ljust(space) + 'Rot Speed'.ljust(space) + 'Max Exp Time'.ljust(space) + 'Time del Theta'.ljust(space) + '\n')
	for i in range(len(proj['K'])):
		fid.write(str(proj['K'][i]).ljust(space) + str(proj['N_theta'][i]).ljust(space) + str(proj['L'][i]).ljust(space) + str(proj['temporal_res'][i]).ljust(space) + str(proj['delta_theta'][i]).ljust(space) + str(proj['rotation_speed'][i]).ljust(space) + str(proj['max_exposure_time'][i]).ljust(space) + str(proj['time_delta_theta'][i]).ljust(space) + '\n')

def main ():
	proj = initialize ()
	gen_experiment_params(proj)
#	write_params_to_txt_file(proj)
	write_params_to_console(proj)

main()
