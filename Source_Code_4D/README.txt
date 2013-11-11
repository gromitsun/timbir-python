************************** Required/Optional Input Files ****************************
projection_n<rank>.bin               - Projection data to be used for reconstruction by node <rank>
weight_n<rank>.bin                   - Weight data to be used for reconstruction by node <rank>
bright_n<rank>.bin                   - Bright field data to be used for reconstruction by node <rank>
NOTE                                 - Use either projection_n<rank> or bright_n<rank> and not both. Further, the 
                                       data should be arranged in row major order. The dimensions of the data should be
				       (N_p, N_r, N_t/(Number of nodes)) where N_p = total number of projections, N_r = 
				       number of detector elements along r-axis (parallel to x), N_t = number of detector elements
				       along t-axis (parallel to z).
proj_offset_n<rank>                  - Contains the additive offset error for each projection assigned to node <rank>. Dimension 
					is (N_r, N_t/(number of nodes)).
object_n<rank>_time_<time_slice>.bin - contains the attenuation coefficients of the object (in um^-1) 
				       reconstructed by node <rank> at time index <time_slice>. Dimension is (N_z, N_y, N_x).
mag_update_map_n<rank> 		     - contains the magnitude of update for each voxel in the last iteration
				       of the last coarse resolution reconstruction.
NOTE 				     - All the binary files listed are required or optional depending on whether the binary files
				       are required for initialization or not.



************************** Compiling C code ********************************
-------------- Source Files to be built ------------------
XT_Engine.c, XT_ICD_update.c, XT_Init.c, XT_genSinogram.c, XT_AMatrix.c,
XT_Profile.c, XT_NHICD.c, allocate.c, randlib.c, tiff.c, XT_IOMisc.c, XT_ImageProc.c, XT_MPI.c, XT_HDFIO.c

--------------- Macros Inputs ----------------
DATA_HDF_FILENAME - Path to HDF file including file name. 
The HDF file should contain the actual count data acquired during the experiment. 
WHITEDARK_HDF_FILENAME - Path to HDF file containing the white (bright field) and dark data (dark field).
PROJECTION_HDF_START - Starting index of projection from which data should be read from the HDF file.
BH_QUAD_COEF - Beamhardening correction parameter. If p is the value of projection, then the corrected projection 
value is computed as p = BH_QUAD_COEF*p*p + p
POSITIVITY_CONSTRAINT - If used, enforces positivity constraint
NO_COST_CALCULATE - If used, cost is not computed

------------ Arguments to compiler/builder ----------
Wall - Enforces some constraints on compilation. Helps avoid mistakes.
ansi - May be used to enforce ansi C code. 
fopenmp - GNU compiler argument which enables openmp based parallelization (REQUIRED ARGUMENT)

------------- Example ---------------------------------
mpicc -Wall -ansi -fopenmp -DDATA_HDF_FILENAME="\"data.hdf5\"" -DWHITEDARK_HDF_FILENAME="\"data.hdf5\"" -DPROJECTION_HDF_START="0" -DBH_QUAD_COEF="0.5" -DPOSITIVITY_CONSTRAINT -DNO_COST_CALCULATE -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c XT_NHICD.c allocate.c randlib.c tiff.c XT_IOMisc.c XT_ImageProc.c XT_MPI.c XT_HDFIO.c -lm 

************************** Running the code ******************************
--------------Input Arguments
p			- 'p' is the qGGMRF sharpness parameter. Typically p = 1.2
sigma_s			- qGGMRF spatial regularization parameter
sigma_t			- qGGMRF temporal regularization parameter
c_s, c_t		- Spatial and temporal qGGMRF parameter which dictates the transition from 
			  quadratic regularization to GGMRF regularization.
delta_xy		- The voxel size as a multiple of detector pixel size in the x-y plane (size is side length).
delta_z			- The voxel size as a multiple of detector pixel size along z axis. 

