1) Copy the entire 3DMBIR folder as is to home directory

2) cd to the 3D_Launcher_for_Real_Data

3) To set the necessary inputs for 

a) Edison  : Open NERSC_Sub_Recon.sh 

b) Carver : Open Carver_Sub_Recon.sh 

4) After setting the mppwidth to the number of nodes*24 go to the main call to change parameters 

Mandatory inputs: 

--setup_launch_folder : a flag to tell the code to create the necessary folder to copy the code into - always keep this
--run_reconstruction  : A flag to tell the code to reconstruct. Always keep this ON.  
--Edison : This is to set HPC specific settings; Options: Edison, Carver, Purdue and PC (for a regular desktop)  
--input_hdf5 : FULL path of data h5 file. Assume the h5 is in the format given to us by Dula
--group_hdf5 : The tile folder to be reconstructed  
--code_launch_folder : FULL path of a folder where the code is going to be copied and called from. Ideally this is in scratch
--output_hdf5 : Full path of a FOLDER where the outputs are going to be written. In a later version this will just be the output H5 file but for now its a folder. 
--pix_size : pixel size in micron 
----num_views : total number of views in the data. If the "0 degree projection"  is acquired two times just indicate the number of unique views. 
--x_width : number of detector pixels in x direction 
--num_bright : Number of brights acquired in the data at the end of the rotation
--num_dark : Number of darks acquired in the data at the end of the rotation
--z_start : the slice number along z in which to start the reconstruction
--z_numElts : number of z slices to reconstruct. For the default case this should be at least 32. If you want to reconstruct less see the advanced parameters section (so that the 3-D multi-resolution works)
--num_nodes : Number of nodes to be used
--num_threads : Number of threads per node
--rot_center : Center of rotation in pixels 

Optional parameters:

--view_subsmpl_fact : The value by which to subset the reconstruction. Can help speed up the code. A number like 2 means take every other view, 3 means take every 3rd view and so on.
 
--recon_x_width : The number of pixels in a subsampled version of the data. If the total number of pixels is 2500 and this number is set to 1024 the acquired data is subsampled to "1024" for each row of the detector. Strongly recommended while getting familiar with the code.  

--smoothness : The MBIR smoothness parameter . Internally a value is computed and then adjusted using this user input. A value of 1 means the code uses the internal value. 2 means a more smooth reconstruction. 1/2 means a less smooth reconstruction. Needs to be adjusted for the best recon. 

--num_res : Total number of resolutions to be used in the multi -resolution optimization. Defaulted to 4. 

----zinger_thresh : A threshold to remove zingers. A very large value means no zingers are corrected. A value like 2-3 may be used to remove zingers. Defaulted to 10000 (i.e. no zinger removal)

----stop_threshold : The stopping criteria - average change in voxel values in Hounsfeld Unit. Defaulted to 35. 

--final_res_multiple : The desired resolution to reconstruct at relative to the detector pixel size. Defaulted to 1 i.e. reconstructions are at the same size as detector. Can be made larger for coarse reconstructions. 

--multires_2D : Set this flag if you wish to do milt-resolution only in 2-D. Can be used to do for example a 3-slice reconstruction quickly. DONOT use for large 3D volumes. 

--Variance_Est : Variance scaling parameter. This is defaulted to 1.
See example of this file in NERSC_Sub_Recon.sh and Carver_Sub_Recon.sh



