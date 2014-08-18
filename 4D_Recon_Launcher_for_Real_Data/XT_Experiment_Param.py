
import argparse
import numpy as np
from XT_Interlaced_Angles import gen_interlaced_views_0_to_Inf,clip_list_of_views
import math

""" 	K - Number of sub-frames in a full frame
	N_theta - Number of projections in a full frame (Total number of discrete angles)
	fps - Frame rate of the camera (Number of projections per second)
	N_p - Total number of projections (used only to compute the deleted angles)
	min_time_btw_views - Minimum time between views. Limit imposed by the trigger system is 5ms  """

def initialize ():
	proj = {}

	detector_resolution = 2560 #Perpendicular to rotation axis
	framerate = 100 #projections per second
	max_time_btw_views = 0.0099 #maximum time between angles (minimum acquire time)
	temporal_resolution = 1 #expected temporal resolution of reconstructions in seconds
	num_of_cycles = 4

	proj['rotation_speed'] = 180/temporal_resolution #slew speed
	proj['L'] = int(math.floor((180.0*framerate)/proj['rotation_speed'])) 
	proj['N_theta'] = (2**math.ceil(math.log(float(detector_resolution)/proj['L'],2)))*proj['L']
	proj['K'] = proj['N_theta']/proj['L']
	proj['N_p'] = proj['N_theta']*num_of_cycles
	proj['fps'] = proj['L']/(180.0/proj['rotation_speed'])

	angles, times = gen_interlaced_views_0_to_Inf(proj['K'], proj['N_theta'], proj['N_p'])
	angles_clip, times_clip, angles_del, times_del = clip_list_of_views (angles, times, max_time_btw_views, proj['rotation_speed'])
	proj['num_angle_del'] = angles_del.size 
	
	return proj

def gen_experiment_params (proj):
	proj['temporal_res'] = proj['L']/proj['fps']
	proj['delta_theta'] = 180.0/proj['N_theta']

	proj['max_exposure_time'] = 1.0/(proj['fps']*proj['K'])
	proj['time_delta_theta'] = proj['delta_theta']/proj['rotation_speed']
	
def write_params_to_console (proj):
	print 'K =', proj['K'], '\nN_theta =', proj['N_theta'], '\nL =', proj['L'], '\nEffective Framerate =', proj['fps'], '\nTemporal Resolution =', proj['temporal_res'], '\ndelta Theta =', proj['delta_theta'], '\nRotation Speed =', proj['rotation_speed'], '\nMax Exposure Time (Time to traverse delta Theta) =', proj['max_exposure_time'], '\nNumber of angles deleted =', proj['num_angle_del'], '\nNumber of remaining angles =', proj['N_p']-proj['num_angle_del'], ' \n '

def write_params_to_txt_file (proj):
	fid = open('experiment_parameters.txt','w')
	space = 15
	fid.write('K'.ljust(space) + 'N_theta'.ljust(space) + 'L'.ljust(space) + 'Eff FPS'.ljust(space)  + 'Time Res'.ljust(space) + 'del Theta'.ljust(space) + 'Rot Speed'.ljust(space) + 'Max Exp Time'.ljust(space) + 'Time del Theta'.ljust(space) + '\n')
	fid.write(str(proj['K']).ljust(space) + str(proj['N_theta']).ljust(space) + str(proj['L']).ljust(space) + str(proj['fps']).ljust(space) + str(proj['temporal_res']).ljust(space) + str(proj['delta_theta']).ljust(space) + str(proj['rotation_speed']).ljust(space) + str(proj['max_exposure_time']).ljust(space) + str(proj['time_delta_theta']).ljust(space) + '\n')

def main ():
	proj = initialize ()
	gen_experiment_params(proj)
#	write_params_to_txt_file(proj)
	write_params_to_console(proj)

main()
