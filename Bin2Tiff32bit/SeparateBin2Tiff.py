#! /usr/bin/python

from os import system
import numpy as np 
from math import ceil

#Code to take in the inputs and generate the tiffs via calls to the Bin2Tiff executable produced from the C++ code.

#The code assume that bin files corresponding to node 0 say have the form
# object_n0_time_0.bin and so on. 

def Convert2Tiffs (inputs):

       slices_per_node = inputs['z_numElts']/inputs['num_nodes']

       for i in range(0,inputs['num_nodes']):
           temp_z_start = inputs['z_start']+ i*slices_per_node
           tmp_file = 'object_n'+str(i)+'_time_0.bin'
           flag = system('./Bin2Tiff ' + inputs['input_path'] + tmp_file + ' ' + inputs['output_tiff_path'] + ' ' + inputs['tiff_base_name'] + ' ' + str(inputs['im_width']) + ' ' + str(temp_z_start) + ' ' + str(slices_per_node) + ' ' + str(inputs['pix_size']))	
       	

       return flag
