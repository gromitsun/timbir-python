
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
