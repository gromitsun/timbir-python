""" ============================================================================
 * Copyright (c) 2013 K. Aditya Mohan (Purdue University)
 * Copyright (c) 2013 Singanallur Venkatakrishnan (Purdue University)
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
 * Neither the name of K. Aditya Mohan, Singanallur Venkatakrishnan, Purdue
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

from XT_Initialize_v2 import proj_init, recon_init, files_init
from XT_MBIR_Reconstruction import do_MBIR_reconstruction
#from XT_FBP_Reconstruction import do_FBP_reconstruction
from XT_IOMisc import error_by_flag
from XT_IOMisc import write_object2HDF
#from XT_ObjectHDFIO import writepar_object2HDF
#from XT_ObjectHDFIO import writepar_tiff_from_object_bin_file
from XT_IOMisc import write_tiff_from_object_bin_file
import argparse
import time
import os
import math
#from mpi4py import MPI
def main():
	
	parser = argparse.ArgumentParser()

        #Parameters associated with the acquired data
        parser.add_argument("--input_hdf5", help="Full path of the input hdf5 file")
        parser.add_argument("--group_hdf5", help="Full path of the input group to reconstruct")
	parser.add_argument("--output_hdf5", help="Full path of the output hdf5 file")
        parser.add_argument("--code_launch_folder", help="Full path where the code executable is going to reside")
        parser.add_argument("--rot_center", help="Center of rotation of the object in units of detector pixels",type=float,default=-1)
        parser.add_argument("--pix_size", help="Size of the detector pixels in micro meter",type=float)
	parser.add_argument("--num_views", help="Number of views for the data",type=int)
        parser.add_argument("--num_bright",help="Number of bright fields acquired prior the the experiment",type=int)
        parser.add_argument("--num_dark",help="Number of dark fields acquired prior the the experiment",type=int)
        parser.add_argument("--view_subsmpl_fact", help="View skipping factor. This can be used to subset the data in terms of views",type=int,default=1)
    
        #Reconstruction parameters
        parser.add_argument("--x_width", help="Number of detector elements along x-direction",type=int)
        parser.add_argument("--recon_x_width", help="Number of detector elements along x-direction",type=int,default=-1)
        parser.add_argument("--z_start", help="Starting detector pixel along z-direction",type=int)
        parser.add_argument("--z_numElts", help="Number of detector pixels to use along z-direction",type=int)
        parser.add_argument("--p",help="qGGMRF prior model parameter use to control how smooth the edge smoothness",type=float,default=1.2)
        parser.add_argument("--final_res_multiple",help="Final desired resolution of the reconstruction as a multiple of the detector pixel size",type=int,default=1)
	parser.add_argument("--smoothness",help="A scaling parameter use to control the balance between resolution and noise in the reconstruction. This value is multiplied by a predtermined internal value",type=float,default=0.50)
        parser.add_argument("--zinger_thresh",help="Controls the rejection thresold of bad measurements that cause zingers. In a future version this will have proper units. At present try values in the range 1-50 to correct for zingers",type=float,default=10000)
       
        #Advanced parameters which the user need not worry about but can manipulate if necessary
        parser.add_argument("--stop_threshold",help="Stopping thresold as a percentage of average change in pixel values in percentage",type=float,default=10)       
        parser.add_argument("--max_iter",help="Maximum number of ICD iterations for the algorithm",type=int,default=30)
        parser.add_argument("--Variance_Est",help="Enter the an estimate for the variance",type=float,default=10)
        parser.add_argument("--num_res",help="Number of resolutions to be used",type=int,default=4)
        parser.add_argument("--multires_2D", help="specify whether to use only 2-D multiresolution. This will slow the 3-D reconstruction but is useful for reconstructing a few slices quickly", action="store_true")

        parser.add_argument("--num_nodes",help="Number of nodes to use",type=int,default=1)
	parser.add_argument("--num_threads",help="Number of threads per node",type=int,default=32)
	
        #Parameters from the old Main program which we need to talk about removing 
        parser.add_argument("--create_objectHDFtiffonly", help="specify whether you want to create HDF and tiff output files only", action="store_true")
	parser.add_argument("--setup_launch_folder", help="Specify whether you want to setup the launch folder", action="store_true")
	parser.add_argument("--run_reconstruction", help="Run reconstruction code", action="store_true")
	parser.add_argument("--Edison", help="Use Edison when running on Edison", action="store_true")
	parser.add_argument("--Purdue", help="Use Purdue when running on Conte or Carter", action="store_true")
        parser.add_argument("--Carver", help="Used to indicate HPC", action="store_true")
        parser.add_argument("--PC", help="Used to indicate code run is on a PC", action="store_true")

	parser.add_argument("--real_is_double", help="Do all reconstruction using double type as 'real' variable", action="store_true")
        args = parser.parse_args()
        
        inputs = {}
        inputs['input_hdf5'] = args.input_hdf5
        inputs['group_hdf5'] = args.group_hdf5
        #The group name has to start with a "/"
        if (inputs['group_hdf5'][0] != '/'):
            print 'Group name has to start with a /'
            return -1

        inputs['output_hdf5']= args.output_hdf5
	inputs['code_launch_folder']=args.code_launch_folder

	inputs['Variance_Est']=args.Variance_Est
        inputs['pix_size'] = args.pix_size
        inputs['num_bright'] = args.num_bright
        inputs['num_dark'] = args.num_dark
        inputs['num_views'] = args.num_views
	inputs['view_subsmpl_fact'] = args.view_subsmpl_fact
        inputs['x_width'] = args.x_width

        if(args.recon_x_width != -1):
            inputs['recon_x_width'] = args.recon_x_width
        else:
            inputs['recon_x_width'] = inputs['x_width']

        if(args.rot_center != -1): #This means the user has input a center of rotation
           if(args.x_width == args.recon_x_width): #If the recon width is same as original det, use same value
	       inputs['rot_center'] = args.rot_center        
           else: #else compute the center for this cropped/subsampled data
               temp = inputs['recon_x_width']/2.0 + (args.recon_x_width/(2.0**math.floor(math.log(args.x_width,2))))*(args.rot_center-args.x_width/2)
               inputs['rot_center'] = temp
               print 'New computed center ='
               print inputs['rot_center'] 
        else: #If nothing is entered simple default it to the center of the detector
	   inputs['rot_center'] = inputs['recon_x_width']/2

        inputs['z_start'] = args.z_start
        inputs['z_numElts'] = args.z_numElts
        inputs['p'] = args.p
        inputs['final_res_multiple'] = args.final_res_multiple
        inputs['smoothness'] = args.smoothness
        inputs['stop_threshold'] = args.stop_threshold
        inputs['max_iter'] = args.max_iter
        inputs['num_res'] = args.num_res
        

        if (args.multires_2D):
            inputs['multires_2D'] = 1 #A boolean variable
        else:
            inputs['multires_2D'] = 0
        print inputs['multires_2D']

        inputs['num_views'] = args.num_views
        inputs['zinger_thresh'] = args.zinger_thresh
        inputs['num_nodes'] = args.num_nodes   
	inputs['num_threads']= args.num_threads

	recon = {}
	files = {}
	recon['node_num'] = inputs['num_nodes']
	if (args.real_is_double):
		recon['real_var_type'] = 'double'
	else:
		recon['real_var_type'] = 'float'

        if (args.Purdue):
		recon['num_threads'] = inputs['num_threads']
		files['scratch'] = os.environ['RCAC_SCRATCH']
		files['data_scratch'] = os.environ['RCAC_SCRATCH']
		recon['run_command'] = 'mpiexec -n ' + str(recon['node_num']) + ' -machinefile nodefile '
		recon['compile_command'] = 'mpicc -openmp '
		recon['HPC'] = 'Purdue' 
		recon['rank'] = 0
	elif (args.Edison):
		recon['num_threads'] = inputs['num_threads']
		files['scratch'] = os.environ['SCRATCH']
		files['data_scratch'] = os.environ['GSCRATCH']
		recon['run_command'] = 'aprun -j 2 -n ' + str(recon['node_num']) + ' -N 1 -d ' + str(recon['num_threads']) + ' -cc none '
		recon['compile_command'] = 'cc -openmp '
		recon['HPC'] = 'NERSC'
		recon['rank'] = 0		
	elif (args.Carver):
		recon['num_threads'] = inputs['num_threads']
		files['scratch'] = os.environ['SCRATCH']
		files['data_scratch'] = os.environ['GSCRATCH']
		recon['run_command'] = 'mpirun -np ' + str(recon['node_num']) + ' -bynode '
		recon['compile_command'] = 'mpicc -openmp '
		recon['HPC'] = 'NERSC'
		recon['rank'] = 0
	elif (args.PC):
		recon['num_threads'] = inputs['num_threads']
		files['scratch'] = '../XT_run/'
		files['data_scratch'] = '../XT_run/'
		recon['run_command'] = 'mpiexec -n ' + str(recon['node_num']) 
		recon['compile_command'] = 'mpicc -fopenmp '
		recon['HPC'] = 'PC'
		recon['rank'] = 0
	else:
		error_by_flag(1, 'HPC system not recognized')
       
        proj = proj_init(inputs)
        recon = recon_init(proj, recon,inputs)
        files = files_init(files,inputs)

	recon['set_up_launch_folder'] = 0	
	if args.setup_launch_folder:
		recon['set_up_launch_folder'] = 1

	recon['reconstruct'] = 0
	if (args.run_reconstruction and args.create_objectHDFtiffonly == False):
		recon['reconstruct'] = 1

	if (recon['recon_type'] == 'MBIR'):
		print 'main: Will do MBIR reconstruction'
		do_MBIR_reconstruction(proj, recon, files)
		write_object2HDF (proj, recon, files)
	else:
		print 'ERROR: main: Reconstruction type not recognized'

	if (args.create_objectHDFtiffonly):
#		if (recon['HPC'] == 'Purdue'):
#			recon['size'] = MPI.COMM_WORLD.size
#			writepar_object2HDF (proj, recon, files)
#			writepar_tiff_from_object_bin_file (proj, recon, files)
#		else:
		write_tiff_from_object_bin_file (proj, recon, files)
    	
#        print proj
#        print recon
#        print files		
        
	print 'main: Done!'

main()
