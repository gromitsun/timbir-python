#! /usr/bin/python

from SeparateBin2Tiff import Convert2Tiffs
import argparse
import time
import os
import math

def main():
	
	parser = argparse.ArgumentParser()

        #Parameters associated with the acquired data
        parser.add_argument("--input_path", help="Full path of the input bin files")
	parser.add_argument("--output_tiff_path", help="Full path of the output tiffs")
        parser.add_argument("--tiff_base_name", help="A name for files to which the node number will be appended and .tif added")
        parser.add_argument("--num_nodes",help="Number of nodes to use",type=int,default=1)
        parser.add_argument("--z_start", help="Starting detector pixel along z-direction",type=int)
        parser.add_argument("--z_numElts", help="Number of detector pixels to use along z-direction used for the reconstruction",type=int)
        parser.add_argument("--im_width", help="Width of the reconstructed image in pixels. The code assume the height is same as this",type=int)
        parser.add_argument("--pix_size", help="Size of the detector pixels in micro meter",type=float)   
        
        args = parser.parse_args()
        
        inputs = {}
        inputs['input_path'] = args.input_path
        inputs['output_tiff_path']= args.output_tiff_path
        inputs['tiff_base_name'] = args.tiff_base_name
        inputs['pix_size'] = args.pix_size
        inputs['im_width'] = args.im_width
        inputs['z_start'] = args.z_start
        inputs['z_numElts'] = args.z_numElts
        inputs['num_nodes'] = args.num_nodes

        print inputs
        
        exit_status = Convert2Tiffs(inputs)

	print exit_status    
	print 'main: Done!'

main()
