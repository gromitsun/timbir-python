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




from math import log
import numpy as np

def clip_list_of_views (views, times, min_time_btw_views, rotation_speed):
	min_dTheta = rotation_speed*min_time_btw_views

	new_views = np.zeros(len(views), dtype = np.float64, order='C')
	new_times = np.zeros(len(views), dtype = np.float64, order='C')
	deleted_views = np.zeros(len(views), dtype = np.float64, order='C')
	deleted_times = np.zeros(len(views), dtype = np.float64, order='C')
	
	idx_keep = 1
	idx_rem = 0
	new_views[0] = views[0]
	for i in range(1,len(views)):
		if (views[i] - views[i-1] > min_dTheta):
			new_views[idx_keep] = views[i]
			new_times[idx_keep] = times[i]
			idx_keep = idx_keep + 1
		else:
			deleted_views[idx_rem] = views[i]
			deleted_times[idx_rem] = times[i]
			idx_rem = idx_rem + 1

	return (new_views[0:idx_keep], new_times[0:idx_keep], deleted_views[0:idx_rem], deleted_times[0:idx_rem])
	

def bit_reverse_decimal (val, maxbits):
	rev_val = np.zeros(len(val), dtype = np.float64, order='C')
	for i in range(len(val)):
		rev_val[i] = int(np.binary_repr(val[i].astype(np.int32), width=maxbits)[::-1], 2)
	return rev_val

def gen_interlaced_views_0_to_180 (K, N_theta, N_p):

	# Generates interlaced views as a function of projection time index 

	############ Input Variables ###############
	#
	# N_p - Total number of projections
	# K - Number of sub-frames in a full frame
	# N_theta - Total number of discrete projection angles
	#
	##########################################

	########### Output Variables #################
	#
	# views - Angles at which projections should be taken
	#         Every angle is between 0 and 180.
	#
	###############################################

	L = N_theta/K
	maxbits = int(log(K, 2))
	delta_theta = 180.0/N_theta

	i = np.arange(0, N_p, 1, dtype=np.float64)
	buf1 = np.mod(i, L)*K
	temp = np.mod(np.floor(i/L), K)
	buf2 = bit_reverse_decimal(temp, maxbits)

	views = (buf1 + buf2)
	times = views/K
	views = views*delta_theta

	return (views,times)
	


def gen_interlaced_views_0_to_360 (K, N_theta, N_p):

	# Generates interlaced views as a function of projection time index 

	############ Input Variables ###############
	#
	# N_p - Total number of projections
	# K - Number of sub-frames in a full frame
	# N_theta - Total number of discrete projection angles
	#
	##########################################

	########### Output Variables #################
	#
	# views - Angles at which projections should be taken
	#         Every angle is between 0 and 360.
	#
	###############################################

	L = N_theta/K
	maxbits = int(log(K, 2))
	delta_theta = 180.0/N_theta

	i = np.arange(0, N_p, 1, dtype=np.float64)
	buf1 = np.mod(i, 2*L)*K
	temp = np.mod(np.floor(i/L), K)
	buf2 = bit_reverse_decimal(temp, maxbits)

	views = (buf1 + buf2)
	times = views/K
	views = views*delta_theta

	return (views,times)


def gen_interlaced_views_0_to_Inf (K, N_theta, N_p):

	# Generates interlaced views as a function of projection time index 

	############ Input Variables ###############
	#
	# N_p - Total number of projections
	# K - Number of sub-frames in a full frame
	# N_theta - Total number of discrete projection angles
	#
	##########################################

	########### Output Variables #################
	#
	# views - Angles at which projections should be taken
	#         Every angle is between 0 and Infinity.
	#
	###############################################

	L = N_theta/K
	maxbits = int(log(K, 2))
	delta_theta = 180.0/N_theta

	i = np.arange(0, N_p, 1, dtype=np.float64)
	buf1 = i*K
	temp = np.mod(np.floor(i/L), K)
	buf2 = bit_reverse_decimal(temp, maxbits)

	views = (buf1 + buf2)
	times = views/K
	views = views*delta_theta

	return (views,times)



