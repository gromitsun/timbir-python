/* ============================================================================
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
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */



#ifndef XT_STRUCTURES_H
#define XT_STRUCTURES_H
#include <stdint.h>
#include "XT_Constants.h"
#include <stdbool.h>
 /* Axes conventions:

        . X
       .
      .
     .
    .
   .
   ---------------------> Y
   |
   |
   |
   |
   |
   |
   V
   Z
   */
typedef double Real_t;

/*Structure 'Sinogram' contains the sinogram itself and also other parameters related to the sinogram and the detector*/
  typedef struct
  {
   Real_t ***Projection; /*Stores the sinogram*/
   Real_t **ProjOffset; /*Projection offset*/
   bool ***ProjSelect;
    int32_t N_r;/*Number of measurements in r direction*/
    int32_t N_t;/*Number of measurements in t direction to be reconstructed*/
    int32_t N_p;/*Number of available angles*/
    int32_t total_t_slices;/*Number of slices in t-direction*/
    Real_t delta_r;/*Distance between successive measurements along r*/
    Real_t delta_t;/*Distance between successive measurements along t*/
    Real_t R0,RMax;/*location of leftmost corner of detector along r*/
    Real_t T0,TMax;/*location of rightmost corner of detector along r*/
    Real_t *cosine, *sine; /*stores the cosines and sines of the values in angles. angles consist of the angles of the different views of the object with respect to the detector*/
    Real_t Length_R; /*Length of the detector along the r-dimension*/
    Real_t Length_T; /*Length of the detector along the t-dimension*/
    Real_t OffsetR; /*increments of distance between the center of the voxel and the midpoint of the detector along r axis */
    Real_t OffsetT; /*increments of distance between the center of the voxel and the midpoint of the detector along t axis*/
    Real_t *ViewPtr;
    Real_t *TimePtr;

    int32_t slice_begin; /*Detector slice begin*/
    int32_t slice_num; /*Detector slice end*/
  } Sinogram;

  typedef struct
  {
    Real_t ****Object; /*Stores the reconstructed object*/
    Real_t Length_X;/*maximum possible length of the object along x*/
    Real_t Length_Y;/*max length of object along y*/
    Real_t Length_Z;/*max length of object along z*/
    int32_t N_x;/*Number of voxels in x direction*/
    int32_t N_z;/*Number of voxels in z direction*/
    int32_t N_y;/*Number of voxels in y direction*/
/*Now, we assume the voxel is a cube and also the resolution is equal in all the 3 dimensions*/
    int32_t N_time;/*Number of measurements in time*/
    /*Coordinates of the left corner of the object along x, y and z*/
    Real_t x0;
    Real_t z0;
    Real_t y0;
    Real_t delta_xy;/*Voxel size in the x-y direction*/
    Real_t delta_z;/*Voxel size in the z direction*/
    Real_t mult_xy;
    Real_t mult_z;

/*again, for the time being we assume, delta_x = delta_y = delta_z*/ 
    Real_t BeamWidth; /*Beamwidth of the detector response*/
    Real_t Sigma_S; /*regularization parameter over space (over x-y slice). Its a parameter of qGGMRF prior model*/ 
    Real_t Sigma_T; /*regularization parameter across time. Its a parameter of qGGMRF prior model*/ 

    Real_t MRF_P; /*parameter p of qGGMRF prior model*/
    Real_t C_S; /*parameter c of qGGMRF prior model*/
    Real_t C_T; /*parameter c of qGGMRF prior model*/

    int32_t **ProjIdxPtr;
    int32_t *ProjNum;

   Real_t gamma;
   Real_t Rtime0;
   Real_t delta_Rtime;
  } ScannedObject;

  /*Structure to store a single column(A_i) of the A-matrix*/
  typedef struct
  {
      Real_t* values; /*Store the non zero entries*/
      uint8_t count; /*The number of non zero values present in the column*/
      int32_t *index; /*This maps each value to its location in the column.*/
  } AMatrixCol;

typedef struct
  {
    int32_t NumIter; /*Maximum number of iterations that the ICD can be run. Normally, ICD converges before completing all the iterations and exits*/
    Real_t StopThreshold; /*ICD exits after the average update of the voxels becomes less than this threshold. Its specified in units of HU.*/
    Real_t RotCenter; /*Center of rotation of the object as measured on the detector in units of pixels*/ 
    int32_t CountAngles; 
    
    int32_t phantom_N_xy; 
    int32_t phantom_N_z; 
     
    Real_t radius_obj;	
    Real_t Sigma_S_Q;
    Real_t Sigma_T_Q;
    Real_t Sigma_S_Q_P;
    Real_t Sigma_T_Q_P;
    
    Real_t Spatial_Filter[NHOOD_Z_MAXDIM][NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM]; /*Filter is the weighting kernel used in the prior model*/
    Real_t Time_Filter[(NHOOD_TIME_MAXDIM-1)/2]; /*Filter is the weighting kernel used in the prior model*/
    
    Real_t*** Weight;
    Real_t var_est;
    
    Real_t alpha; /*Value of over-relaxation*/
    Real_t cost_thresh; 
    uint8_t time_reg; 
    uint8_t sinobin; /*If set, intializes the sinogram from bin file*/
    uint8_t initICD; /*used to specify the method of initializing the object before ICD*/
    uint8_t Write2Tiff; /*If set, tiff files are written*/
    uint8_t No_Projection_Noise; 
    uint8_t reconstruct; 
    int32_t num_threads;
    uint8_t updateProjOffset;
    uint8_t no_NHICD;
    uint8_t WritePerIter;
    int32_t num_z_blocks;	
    int32_t*** x_NHICD_select;
    int32_t*** y_NHICD_select;
    int32_t*** x_rand_select;
    int32_t*** y_rand_select;
    int32_t** UpdateSelectNum;
    int32_t** NHICDSelectNum;

    Real_t ErrorSinoThresh;
    Real_t ErrorSinoDelta;

    uint8_t only_Edge_Updates;
    int32_t max_HICD_iter;
    int32_t node_num;
    int32_t node_rank;

    uint8_t updateVar;
    uint8_t initMagUpMap;
    uint8_t readSino4mHDF;
    FILE *debug_file_ptr;
  } TomoInputs;


#endif /*#define XT_STRUCTURES_H*/
