#! /usr/bin/python

import argparse
import time
import os
from XT_Initialize_v2 import proj_init, recon_init, files_init

def main():
	
	parser = argparse.ArgumentParser()

        #Parameters associated with the acquired data
        parser.add_argument("--input_hdf5", help="Full path of the input hdf5 file")
	parser.add_argument("--output_hdf5", help="Full path of the output hdf5 file")
        parser.add_argument("--rot_center", help="Center of rotation of the object in units of detector pixels",type=float,default=1280)
        parser.add_argument("--pix_size", help="Size of the detector pixels in micro meter",type=float,default=.6)
	parser.add_argument("--num_views", help="Number of views for the data",type=int)

        #Reconstruction parameters
        parser.add_argument("--x_width", help="Number of detector elements to use along x-direction",type=int,default=2560)
        parser.add_argument("--z_start", help="Starting detector pixel along z-direction",type=int,default=0)
        parser.add_argument("--z_numElts", help="Number of detector pixels to use along z-direction",type=int,default=8)
        parser.add_argument("--p",help="qGGMRF prior model parameter use to control how smooth the edge smoothness",type=float,default=1.2)
        parser.add_argument("--final_res_multiple",help="Final desired resolution of the reconstruction as a multiple of the detector pixel size",type=int,default=1)
	parser.add_argument("--smoothness",help="A scaling parameter use to control the balance between resolution and noise in the reconstruction. This value is multiplied by a predtermined internal value",type=float,default=1)
        parser.add_argument("--zinger_thresh",help="Controls the rejection thresold of bad measurements that cause zingers. In a future version this will have proper units. At present try values in the range 1-50 to correct for zingers",type=float,default=10000)
       
        #Advanced parameters which the user need not worry about but can manipulate if necessary
        parser.add_argument("--stop_threshold",help="Stopping thresold as a percentage of average change in pixel values in units of HU",type=float,default=20)       
        parser.add_argument("--max_iter",help="Maximum number of ICD iterations for the algorithm",type=int,default=30)
        parser.add_argument("--num_res",help="Number of resolutions to be used",type=int,default=4)
        parser.add_argument("--num_nodes",help="Number of nodes to use",type=int,default=1)
	parser.add_argument("--num_threads",help="Number of threads per node",type=int,default=32)
	
        args = parser.parse_args()
        
        inputs = {}
        inputs['input_hdf5'] = args.input_hdf5
        inputs['output_hdf5']= args.output_hdf5
        inputs['rot_center'] = args.rot_center        
        inputs['pix_size'] = args.pix_size
        inputs['num_views'] = args.num_views
        inputs['x_width'] = args.x_width
        inputs['z_start'] = args.z_start
        inputs['z_numElts'] = args.z_numElts
        inputs['p'] = args.p
        inputs['final_res_multiple'] = args.final_res_multiple
        inputs['smoothness'] = args.smoothness
        inputs['stop_threshold'] = args.stop_threshold
        inputs['max_iter'] = args.max_iter
        inputs['num_res'] = args.num_res
        inputs['num_views'] = args.num_views
        inputs['zinger_thresh'] = args.zinger_thresh
        inputs['num_nodes'] = args.num_nodes   
	inputs['num_threads']= args.num_threads

	recon = {}
	files = {}
        
        proj = proj_init(inputs)
        recon = recon_init(proj, recon,inputs)
#        files = files_init(files,inputs)
    	
        print proj
        print recon
        print files		
        
	print 'main: Done!'

main()
