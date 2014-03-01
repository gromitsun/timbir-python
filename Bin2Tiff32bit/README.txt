This code in this folder is used to convert the .bin files produced by each node of XT_MBIR_3D
and turn them into 32 bit tiff with values scaled in cm^-1. 

Usage:

1) Compile the C-Code on Carver: 
   module load tiff
   setenv GRIDREC "-L${TIFF_DIR}/lib -I${TIFF_DIR}/include -ltiff"
   g++ -o Bin2Tiff main.cpp grid_readwrite.cpp $GRIDREC 

2) Run the python script with the relevant inputs: 
Bin2Tiff.py [-h] [--input_path INPUT_PATH][--output_tiff_path OUTPUT_TIFF_PATH] [--tiff_base_name TIFF_BASE_NAME] [--num_nodes NUM_NODES][--z_start Z_START] [--z_numElts Z_NUMELTS][--im_width IM_WIDTH] [--pix_size PIX_SIZE]

z_start : starting slice used for the reconstruction of the volume
z_numElts : Number of z slices reconstructed in total 
im_width represents number of  rows of each slice 

For examples see main.cpp and Bin2Tiff.py 