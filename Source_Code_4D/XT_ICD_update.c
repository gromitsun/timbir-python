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



#include "XT_Constants.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "allocate.h"
#include "randlib.h"
#include "XT_genSinogram.h"
#include <time.h>
#include "XT_AMatrix.h"
#include "XT_Profile.h"
#include "XT_Structures.h"
#include "XT_IOMisc.h"
#include "XT_NHICD.h"
#include "omp.h"
#include "XT_ImageProc.h"
#include "XT_MPI.h"
#include <mpi.h>

/*computes the location of (i,j,k) th element in a 1D array*/
int32_t array_loc_1D (int32_t i, int32_t j, int32_t k, int32_t N_j, int32_t N_k)
{
	return (i*N_j*N_k + j*N_k + k);
}

/*finds the maximum in a array 'array_in' with number of elements being 'num'*/
int32_t find_max(int32_t* array_in, int32_t num)
{
	int32_t i, maxnum;
	maxnum = array_in[0];
	for (i=1; i<num; i++)
		if (array_in[i] > maxnum)
			maxnum = array_in[i];
	
	return(maxnum);
}

/*converts the value 'val' to hounsfield units and returns it*/
Real_t convert2Hounsfield (Real_t val)
{
	Real_t slope, c;
	
	slope=(HOUNSFIELD_WATER_MAP-HOUNSFIELD_AIR_MAP)/(WATER_MASS_ATT_COEFF*WATER_DENSITY-AIR_MASS_ATT_COEFF*AIR_DENSITY)/HFIELD_UNIT_CONV_CONST;
	c=-slope*(AIR_MASS_ATT_COEFF*AIR_DENSITY*HFIELD_UNIT_CONV_CONST);
	
	return (slope*val + c);
}

/*Computes the qGGMRF spatial prior cost value at delta = x_i - x_j. i & j being the voxel and its neighbor*/
Real_t CE_QGGMRF_Spatial_Value(Real_t delta, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  return ((pow(fabs(delta),MRF_Q)/TomoInputsPtr->Sigma_S_Q)/(ScannedObjectPtr->C_S + pow(fabs(delta),MRF_Q - ScannedObjectPtr->MRF_P)/TomoInputsPtr->Sigma_S_Q_P));
}

/*Computes the qGGMRF temporal prior cost value at delta = x_i - x_j. i & j being the voxel and its neighbor*/
Real_t CE_QGGMRF_Temporal_Value(Real_t delta, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  return ((pow(fabs(delta),MRF_Q)/TomoInputsPtr->Sigma_T_Q)/(ScannedObjectPtr->C_T + pow(fabs(delta),MRF_Q - ScannedObjectPtr->MRF_P)/TomoInputsPtr->Sigma_T_Q_P));
}

/*Computes the qGGMRF spatial prior derivative at delta = x_i - x_j. i & j being the voxel and its neighbor*/
Real_t CE_QGGMRF_Spatial_Derivative(Real_t delta, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  Real_t temp1,temp2,temp3;
  temp1=pow(fabs(delta),MRF_Q - ScannedObjectPtr->MRF_P)/(TomoInputsPtr->Sigma_S_Q_P);
  temp2=pow(fabs(delta),MRF_Q - 1);
  temp3 = ScannedObjectPtr->C_S + temp1;
  if(delta < 0)
    return ((-1*temp2/(temp3*TomoInputsPtr->Sigma_S_Q))*(MRF_Q - ((MRF_Q-ScannedObjectPtr->MRF_P)*temp1)/(temp3)));
  else
  {
    return ((temp2/(temp3*TomoInputsPtr->Sigma_S_Q))*(MRF_Q - ((MRF_Q-ScannedObjectPtr->MRF_P)*temp1)/(temp3)));
  }
}

/*Computes the qGGMRF temporal prior derivative at delta = x_i - x_j. i & j being the voxel and its neighbor*/
Real_t CE_QGGMRF_Temporal_Derivative(Real_t delta, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  Real_t temp1,temp2,temp3;
  temp1 = pow(fabs(delta),MRF_Q - ScannedObjectPtr->MRF_P)/(TomoInputsPtr->Sigma_T_Q_P);
  temp2 = pow(fabs(delta),MRF_Q - 1); 
  temp3 = ScannedObjectPtr->C_T + temp1;
  if(delta < 0)
    return ((-1*temp2/(temp3*TomoInputsPtr->Sigma_T_Q))*(MRF_Q - ((MRF_Q-ScannedObjectPtr->MRF_P)*temp1)/(temp3)));
  else
  {
    return ((temp2/(temp3*TomoInputsPtr->Sigma_T_Q))*(MRF_Q - ((MRF_Q-ScannedObjectPtr->MRF_P)*temp1)/(temp3)));
  }
}

/*Computes the qGGMRF spatial prior second derivative at delta = 0*/
Real_t CE_QGGMRF_Spatial_SecondDerivative(ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  return MRF_Q/(TomoInputsPtr->Sigma_S_Q*ScannedObjectPtr->C_S);
}

/*Computes the qGGMRF spatial prior second derivative at delta = 0*/
Real_t CE_QGGMRF_Temporal_SecondDerivative(ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  return MRF_Q/(TomoInputsPtr->Sigma_T_Q*ScannedObjectPtr->C_T);
}

/*Computes the voxel update and returns it. V is the present value of voxel.
THETA1 and THETA2 are the values used in voxel update. Spatial_Nhood and Time_Nhood gives the 
values of voxels in the neighborhood of V. Time_BDFlag and Spatial_BDFlag are masks which determine 
whether a neighbor should be included in the neighorhood or not.*/
Real_t CE_FunctionalSubstitution(Real_t V, Real_t THETA1, Real_t THETA2, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1])
{
  Real_t u,temp1=0,temp2=0,temp_const,RefValue=0,Delta0;
  Real_t QGGMRF_Params;
  int32_t i,j,k;
  RefValue = V;

  /*Need to Loop this for multiple iterations of substitute function*/
    for (i=0; i < NHOOD_Y_MAXDIM; i++)
    for (j=0; j < NHOOD_X_MAXDIM; j++)
      for (k=0; k < NHOOD_Z_MAXDIM; k++)
      {
        if(Spatial_BDFlag[i][j][k] == true && (i != (NHOOD_Y_MAXDIM-1)/2 || j != (NHOOD_X_MAXDIM-1)/2 || k != (NHOOD_Z_MAXDIM-1)/2))
        {
	  Delta0  = (RefValue - Spatial_Nhood[i][j][k]);
          if(Delta0 != 0)
          QGGMRF_Params = CE_QGGMRF_Spatial_Derivative(Delta0,ScannedObjectPtr,TomoInputsPtr)/(Delta0);
          else {
            QGGMRF_Params = CE_QGGMRF_Spatial_SecondDerivative(ScannedObjectPtr,TomoInputsPtr);
          }

          temp_const = TomoInputsPtr->Spatial_Filter[i][j][k]*QGGMRF_Params;
          temp1 += temp_const*Spatial_Nhood[i][j][k];
          temp2 += temp_const;
        }
      }

	for (i=0; i < NHOOD_TIME_MAXDIM - 1; i++)
	{
		if(Time_BDFlag[i] == true)
        	{
           		Delta0  = (RefValue - Time_Nhood[i]);
           		if(Delta0 != 0)
           			QGGMRF_Params = CE_QGGMRF_Temporal_Derivative(Delta0,ScannedObjectPtr,TomoInputsPtr)/(Delta0);
           		else {
             			QGGMRF_Params = CE_QGGMRF_Temporal_SecondDerivative(ScannedObjectPtr,TomoInputsPtr);
           		}
 
           	temp_const = TomoInputsPtr->Time_Filter[0]*QGGMRF_Params;
           	temp1 += temp_const*Time_Nhood[i];
           	temp2 += temp_const;	
		}
	}
	
      u=(temp1+ (THETA2*V) - THETA1)/(temp2 + THETA2);
  
      RefValue = RefValue + TomoInputsPtr->alpha*(u-RefValue);
#ifdef POSITIVITY_CONSTRAINT
      if (RefValue <= 0)
	RefValue = 0;
#endif
  return RefValue;
}

/*computes the value of cost function. 'ErrorSino' is the error sinogram*/
Real_t computeCost(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino)
{
  Real_t cost=0,temp=0, forward=0, prior=0;
  Real_t delta;
  int32_t i,j,k,p,N_z;
  bool j_minus, k_minus, i_plus, j_plus, k_plus, p_plus;  

 
  #pragma omp parallel for private(j, k, temp) reduction(+:cost)
  for (i = 0; i < SinogramPtr->N_p; i++)
    for (j = 0; j < SinogramPtr->N_r; j++)
    	for (k = 0; k < SinogramPtr->N_t; k++)
	{
        	temp = ErrorSino[i][j][k] * sqrt(TomoInputsPtr->Weight[i][j][k]);
		if (SinogramPtr->ProjSelect[i][j][k] == true)
			temp = temp*temp;
		else
			temp = 2.0*TomoInputsPtr->ErrorSinoDelta*TomoInputsPtr->ErrorSinoThresh*fabs(temp) + TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoThresh*(1.0-2.0*TomoInputsPtr->ErrorSinoDelta);
		cost += temp;
	}

  cost /= 2.0;

/*When computing the cost of the prior term it is important to make sure that you don't include the cost of any pair of neighbors more than once. In this code, a certain sense of causality is used to compute the cost. We also assume that the weghting kernel given by 'Filter' is symmetric. Let i, j and k correspond to the three dimensions. If we go forward to i+1, then all neighbors at j-1, j, j+1, k+1, k, k-1 are to be considered. However, if for the same i, if we go forward to j+1, then all k-1, k, and k+1 should be considered. For same i and j, only the neighbor at k+1 is considred.*/

  temp = 0;
    N_z = ScannedObjectPtr->N_z + 2;
    if (TomoInputsPtr->node_rank == TomoInputsPtr->node_num-1)
	N_z = ScannedObjectPtr->N_z + 1;
    #pragma omp parallel for private(delta, p, j, k, j_minus, k_minus, p_plus, i_plus, j_plus, k_plus) reduction(+:temp)
    for (i = 0; i < ScannedObjectPtr->N_time; i++)
    for (p = 1; p < ScannedObjectPtr->N_z + 1; p++)
    for (j = 0; j < ScannedObjectPtr->N_y; j++)
    {
      for (k = 0; k < ScannedObjectPtr->N_x; k++)
      {
	j_minus = (j - 1 >= 0)? true : false;
	k_minus = (k - 1 >= 0)? true : false;
	
	p_plus = (p + 1 < N_z)? true : false;

	i_plus = (i + 1 < ScannedObjectPtr->N_time)? true : false;
	j_plus = (j + 1 < ScannedObjectPtr->N_y)?	 true : false;
	k_plus = (k + 1 < ScannedObjectPtr->N_x)? 	 true : false;
	
        if(k_plus == true) {
          delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p][j][k + 1]);
          temp += TomoInputsPtr->Spatial_Filter[1][1][2] * CE_QGGMRF_Spatial_Value(delta,ScannedObjectPtr,TomoInputsPtr);
        }

        if(j_plus == true) {
          if(k_minus == true) {
            delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p][j + 1][k - 1]);
            temp += TomoInputsPtr->Spatial_Filter[1][2][0] * CE_QGGMRF_Spatial_Value(delta,ScannedObjectPtr,TomoInputsPtr);
          }

          delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p][j + 1][k]);
          temp += TomoInputsPtr->Spatial_Filter[1][2][1] * CE_QGGMRF_Spatial_Value(delta,ScannedObjectPtr,TomoInputsPtr);

          if(k_plus == true) {
            delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p][j + 1][k + 1]);
            temp += TomoInputsPtr->Spatial_Filter[1][2][2] * CE_QGGMRF_Spatial_Value(delta,ScannedObjectPtr,TomoInputsPtr);
          }
        }

	if (p_plus == true)
	{    
          if(j_minus == true)
          {    
            delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j - 1][k];
            temp += TomoInputsPtr->Spatial_Filter[2][0][1] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
          } 
        
	  delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p+1][j][k];
          temp += TomoInputsPtr->Spatial_Filter[2][1][1] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);

          if(j_plus == true)
          {    
            delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p+1][j + 1][k];
            temp += TomoInputsPtr->Spatial_Filter[2][2][1] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
          }    

          if(j_minus == true)
          {    
            if(k_minus == true)
            {    
              delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j - 1][k - 1];
              temp += TomoInputsPtr->Spatial_Filter[2][0][0] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
            }    

            if(k_plus == true)
            {    
              delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j - 1][k + 1];
              temp += TomoInputsPtr->Spatial_Filter[2][0][2] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
            }    
          }    

          if(k_minus == true)
          {    
            delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j][k - 1];
            temp += TomoInputsPtr->Spatial_Filter[2][1][0] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
          }    

          if(j_plus == true)
          {    
            if(k_minus == true)
            {    
              delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j + 1][k - 1];
              temp += TomoInputsPtr->Spatial_Filter[2][2][0] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
            }    

            if(k_plus == true)
            {    
              delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j + 1][k + 1];
              temp += TomoInputsPtr->Spatial_Filter[2][2][2] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
            }    
          }    

          if(k_plus == true)
          {    
            delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j][k + 1];
            temp += TomoInputsPtr->Spatial_Filter[2][1][2] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
          }    
        }    

        if(i_plus == true) {
            delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i+1][p][j][k]);
            temp += TomoInputsPtr->Time_Filter[0] * CE_QGGMRF_Temporal_Value(delta,ScannedObjectPtr,TomoInputsPtr);
	 }
        }
      }

 /*Use MPI reduction operation to add the forward and prior costs from all nodes*/
 MPI_Reduce(&cost, &forward, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 MPI_Reduce(&temp, &prior, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

 if (TomoInputsPtr->node_rank == 0)
 {
 	fprintf(TomoInputsPtr->debug_file_ptr, "costCompute: forward cost=%f\n",cost);
 	fprintf(TomoInputsPtr->debug_file_ptr, "costCompute: prior cost =%f\n",temp);
 	cost = forward + prior + TomoInputsPtr->node_num*SinogramPtr->N_p*SinogramPtr->N_r*SinogramPtr->N_t*log(TomoInputsPtr->var_est)/2;
 }
 /*Broadcase the value of cost to all nodes*/
 MPI_Bcast(&cost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return cost;
}

/*Updates all voxels from time slices 'time_begin - time_end', x-y slices split along the z-axis from 'slice_begin-slice_end'.
x_rand_select and y_rand_select lists the voxels for each time slice and z-axis slice which needs to be updated.
ErrorSino - Error sinogram
DetectorResponse_XY - detector response in the x-y plane
VoxelLineResponse - Gives the detector response along z-axis for each voxel along z
Iter - The iteration number
zero_count - the number of zero attenuation coefficients
MagUpdateMap - contains the magnitude of each voxel update in the previous iteration which updated that voxel
Mask - All voxels contained in the 'Mask' (contains true or false values for each voxel) are updated*/
Real_t updateVoxels (int32_t time_begin, int32_t time_end, int32_t slice_begin, int32_t slice_end, int32_t xy_begin, int32_t xy_end, int32_t* x_rand_select, int32_t* y_rand_select, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, Real_t** DetectorResponse_XY, AMatrixCol* VoxelLineResponse, int32_t Iter, long int *zero_count, Real_t** MagUpdateMap, uint8_t*** Mask)
{
  Real_t UpdatedVoxelValue, ProjectionEntry;
  int32_t i_r, i_t;
  int32_t p,q,r,slice,i_new,j_new,k_new,idxr,idxq,idxp,index_xy;
  Real_t V,THETA1,THETA2;
  bool ZSFlag;
  int32_t sino_view;
  int32_t z_min, z_max;
  Real_t total_vox_mag = 0.0;

  z_min = 0;
  z_max = ScannedObjectPtr->N_z + 1;
  if (TomoInputsPtr->node_rank == 0)
	z_min = 1;
  if (TomoInputsPtr->node_rank == TomoInputsPtr->node_num - 1)
	z_max = ScannedObjectPtr->N_z;

    Real_t THETASelTemp, Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM]; 
    Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1]; 
    bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM];
    bool Time_BDFlag[NHOOD_TIME_MAXDIM-1];

  int32_t maxview = find_max(ScannedObjectPtr->ProjNum, ScannedObjectPtr->N_time);
  AMatrixCol* AMatrixPtr = (AMatrixCol*)get_spc(maxview, sizeof(AMatrixCol));
  uint8_t AvgNumXElements = (uint8_t)ceil(3*ScannedObjectPtr->delta_xy/SinogramPtr->delta_r);
  
  for (p = 0; p < maxview; p++)
  {
  	AMatrixPtr[p].values = (Real_t*)get_spc(AvgNumXElements,sizeof(Real_t));
  	AMatrixPtr[p].index  = (int32_t*)get_spc(AvgNumXElements,sizeof(int32_t));
  }  

   /*printf ("time_begin = %d, time_end = %d, slice_begin = %d, slice_end = %d\n", time_begin, time_end, slice_begin, slice_end);*/
   for (i_new = time_begin; i_new <= time_end; i_new++) 
   {
      for (index_xy = xy_begin; index_xy <= xy_end; index_xy++) 
      {
        /*printf ("Entering index\n"); */
	k_new = x_rand_select[index_xy];
        j_new = y_rand_select[index_xy];
    	MagUpdateMap[j_new][k_new] = 0;  
           /*	printf ("Entering mask\n");*/ 
	  for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++)
    	  {
		sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
		calcAMatrixColumnforAngle(SinogramPtr, ScannedObjectPtr, DetectorResponse_XY, &(AMatrixPtr[p]), j_new, k_new, sino_view);
    	  }
          for (slice = slice_begin; slice <= slice_end; slice++) {
           /*	printf ("Entering slice\n"); */
            /*For a given (i,j,k) store its 26 point neighborhood*/           
	    if (Mask[j_new][k_new][slice] == 1)
	    {   
	 	if (i_new - 1 >= 0){
			Time_Nhood[0] = ScannedObjectPtr->Object[i_new-1][slice+1][j_new][k_new];
			Time_BDFlag[0] = true;
		}
		else 
		{
			Time_Nhood[0] = 0.0;
			Time_BDFlag[0] = false;
		}

	    	if (i_new + 1 < ScannedObjectPtr->N_time){
			Time_Nhood[1] = ScannedObjectPtr->Object[i_new+1][slice+1][j_new][k_new];
			Time_BDFlag[1] = true;
		}
		else
		{
			Time_Nhood[1] = 0.0;
			Time_BDFlag[1] = false;
		}
	
	
	 for (p = 0; p < NHOOD_Z_MAXDIM; p++)
	 {
		idxp = slice + p;
		if (idxp >= z_min && idxp <= z_max)
		{
	 		for (q = 0; q < NHOOD_Y_MAXDIM; q++)
         		{
	 			idxq = j_new + q - 1;
                		if(idxq >= 0 && idxq < ScannedObjectPtr->N_y)
         			{
					for (r = 0; r < NHOOD_X_MAXDIM; r++)
					{
		    				idxr = k_new + r - 1;
                    				if(idxr >= 0 && idxr < ScannedObjectPtr->N_x){
	                				Spatial_Nhood[p][q][r] = ScannedObjectPtr->Object[i_new][idxp][idxq][idxr];
        	        				Spatial_BDFlag[p][q][r] = true;
                    				}
						else
						{
	                				Spatial_Nhood[p][q][r] = 0.0;
                    					Spatial_BDFlag[p][q][r] = false;
						}
					}
				}
		 		else
				{
         				for (r = 0; r < NHOOD_X_MAXDIM; r++){
	                			Spatial_Nhood[p][q][r] = 0.0;
                    				Spatial_BDFlag[p][q][r] = false;
					}
				}
                	}
		}
		else
        	{ 
			for (q = 0; q < NHOOD_Y_MAXDIM; q++){
				for (r = 0; r < NHOOD_X_MAXDIM; r++){
	              			Spatial_Nhood[p][q][r] = 0.0;
                   			Spatial_BDFlag[p][q][r] = false;
				}
			}
               }
	}

        Spatial_Nhood[(NHOOD_Y_MAXDIM-1)/2][(NHOOD_X_MAXDIM-1)/2][(NHOOD_Z_MAXDIM-1)/2] = 0.0;
        V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]; /*Store the present value of the voxel*/

#ifdef ZERO_SKIPPING
			  /*Zero Skipping Algorithm*/
			 ZSFlag = true;
			 if(V == 0.0 && Iter > 1) /*Iteration starts from 1. Iteration 0 corresponds to initial cost before ICD*/
			  {
					  
					if (Time_Nhood[0] > 0.0 || Time_Nhood[1] > 0.0)
						ZSFlag = false;
			
					for(p = 0; p < NHOOD_Y_MAXDIM; p++)
						for(q = 0; q < NHOOD_X_MAXDIM; q++)
					  		for(r = 0; r < NHOOD_Z_MAXDIM; r++)
							  	if(Spatial_Nhood[p][q][r] > 0.0)
							  	{
									  ZSFlag = false;
								 	  break;
							  	}
			  }
			  else
			  {
				  ZSFlag = false;
			  }
#else
			  ZSFlag = false; /*do ICD on all voxels*/
#endif /*#ifdef ZERO_SKIPPING*/
	if(ZSFlag == false)
	{
		THETA1 = 0.0;
		THETA2 = 0.0;
		for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
		sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
		for (q = 0; q < AMatrixPtr[p].count; q++)
		{
               	    	i_r = (AMatrixPtr[p].index[q]);
        	    	ProjectionEntry = (AMatrixPtr[p].values[q]);
			for (r = 0; r < VoxelLineResponse[slice].count; r++)
			{ 
				i_t = VoxelLineResponse[slice].index[r];
				if (SinogramPtr->ProjSelect[sino_view][i_r][i_t] == true)
				{
	            			THETA2 += (VoxelLineResponse[slice].values[r]*VoxelLineResponse[slice].values[r]*ProjectionEntry*ProjectionEntry*TomoInputsPtr->Weight[sino_view][i_r][i_t]);
               				THETA1 += -(ErrorSino[sino_view][i_r][i_t]*VoxelLineResponse[slice].values[r]*ProjectionEntry*TomoInputsPtr->Weight[sino_view][i_r][i_t]);
            			}
				else
				{
					THETASelTemp = TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoDelta*sqrt(TomoInputsPtr->Weight[sino_view][i_r][i_t])/fabs(ErrorSino[sino_view][i_r][i_t]);
	            			THETA2 += (VoxelLineResponse[slice].values[r]*VoxelLineResponse[slice].values[r]*ProjectionEntry*ProjectionEntry*THETASelTemp);
            				THETA1 += -(ErrorSino[sino_view][i_r][i_t]*VoxelLineResponse[slice].values[r]*ProjectionEntry*THETASelTemp);
				}
            		}
		}
    }


            /*Solve the 1-D optimization problem
            TODO : What if theta1 = 0 ? Then this will give error*/

            UpdatedVoxelValue = CE_FunctionalSubstitution(V, THETA1, THETA2, ScannedObjectPtr, TomoInputsPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag);
              
            ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] = UpdatedVoxelValue;
	
	    for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
		sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
		for (q = 0; q < AMatrixPtr[p].count; q++)
        	{
               	    	i_r = (AMatrixPtr[p].index[q]);
        	    	ProjectionEntry = (AMatrixPtr[p].values[q]);
			for (r = 0; r < VoxelLineResponse[slice].count; r++)
			{ 
				i_t = VoxelLineResponse[slice].index[r];
	        		ErrorSino[sino_view][i_r][i_t] -= (ProjectionEntry*VoxelLineResponse[slice].values[r]*(ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] - V));
	   			if (fabs(ErrorSino[sino_view][i_r][i_t]*sqrt(TomoInputsPtr->Weight[sino_view][i_r][i_t])) < TomoInputsPtr->ErrorSinoThresh)
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = true;
				else
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = false;
	   		}
		}
	    }
	    MagUpdateMap[j_new][k_new] += fabs(ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] - V);
	    total_vox_mag += fabs(ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]);
 	}
		else
		    (*zero_count)++;
       }
       }
       }
}

    
     for (p=0; p<maxview; p++)
     {
     	free(AMatrixPtr[p].values);
     	free(AMatrixPtr[p].index);
     }
     free(AMatrixPtr);
      return (total_vox_mag);
}

/*Upsamples the (N_time x N_z x N_y x N_x) size 'Init' by a factor of 2 along the x-y plane and stores it in 'Object'*/
void upsample_bilinear_2D (Real_t**** Object, Real_t**** Init, int32_t N_time, int32_t N_z, int32_t N_y, int32_t N_x)
{
	int32_t i, j, k, m;
	Real_t **buffer;
	
	#pragma omp parallel for private(buffer, m, j, k)
	for (i=0; i < N_time; i++)
	for (m=0; m < N_z; m++)
	{
        	buffer = (Real_t**)multialloc(sizeof(Real_t), 2, N_y, 2*N_x);
                for (j=0; j < N_y; j++){
			buffer[j][0] = Init[i][m][j][0];
                        buffer[j][1] = (3.0*Init[i][m][j][0] + Init[i][m][j][1])/4.0;
                        buffer[j][2*N_x - 1] = Init[i][m][j][N_x - 1];
                        buffer[j][2*N_x - 2] = (Init[i][m][j][N_x - 2] + 3.0*Init[i][m][j][N_x - 1])/4.0;
                        for (k=1; k < N_x - 1; k++){
				buffer[j][2*k] = (Init[i][m][j][k-1] + 3.0*Init[i][m][j][k])/4.0;
                                buffer[j][2*k + 1] = (3.0*Init[i][m][j][k] + Init[i][m][j][k+1])/4.0;
                        }
                }

                for (k=0; k < 2*N_x; k++){
                        Object[i][m][0][k] = buffer[0][k];
                        Object[i][m][1][k] = (3.0*buffer[0][k] + buffer[1][k])/4.0;
                        Object[i][m][2*N_y-1][k] = buffer[N_y-1][k];
                        Object[i][m][2*N_y-2][k] = (buffer[N_y-2][k] + 3.0*buffer[N_y-1][k])/4.0;
                }

                for (j=1; j<N_y-1; j++){
                        for (k=0; k<2*N_x; k++){
				Object[i][m][2*j][k] = (buffer[j-1][k] + 3.0*buffer[j][k])/4.0;
                                Object[i][m][2*j + 1][k] = (3*buffer[j][k] + buffer[j+1][k])/4.0;
                        }
                }
                multifree(buffer,2);
        }

}

/*Upsamples the (N_z x N_y x N_x) size 'Init' by a factor of 2 along the x-y plane and stores it in 'Object'*/
void upsample_object_bilinear_2D (Real_t*** Object, Real_t*** Init, int32_t N_z, int32_t N_y, int32_t N_x)
{
	int32_t j, k, slice;
	Real_t **buffer;
	
	
  		buffer = (Real_t**)multialloc(sizeof(Real_t), 2, N_y, 2*N_x);
		for (slice=0; slice < N_z; slice++){
			for (j=0; j < N_y; j++){
				buffer[j][0] = Init[slice][j][0];
				buffer[j][1] = (3.0*Init[slice][j][0] + Init[slice][j][1])/4.0;
				buffer[j][2*N_x - 1] = Init[slice][j][N_x - 1];
				buffer[j][2*N_x - 2] = (Init[slice][j][N_x - 2] + 3.0*Init[slice][j][N_x - 1])/4.0;
				for (k=1; k < N_x - 1; k++){
					buffer[j][2*k] = (Init[slice][j][k-1] + 3.0*Init[slice][j][k])/4.0;
					buffer[j][2*k + 1] = (3.0*Init[slice][j][k] + Init[slice][j][k+1])/4.0;
				}
			}

			for (k=0; k < 2*N_x; k++){
				Object[slice+1][0][k] = buffer[0][k];
				Object[slice+1][1][k] = (3.0*buffer[0][k] + buffer[1][k])/4.0;
				Object[slice+1][2*N_y-1][k] = buffer[N_y-1][k];
				Object[slice+1][2*N_y-2][k] = (buffer[N_y-2][k] + 3.0*buffer[N_y-1][k])/4.0;	
			}		

			for (j=1; j<N_y-1; j++){
				for (k=0; k<2*N_x; k++){				
					Object[slice+1][2*j][k] = (buffer[j-1][k] + 3.0*buffer[j][k])/4.0;
					Object[slice+1][2*j + 1][k] = (3*buffer[j][k] + buffer[j+1][k])/4.0;
				}
			}
		}
		multifree(buffer,2);

}


/*Upsamples the (N_time x N_z x N_y x N_x) size 'Init' by a factor of 2 along the in 3D x-y-z coordinates and stores it in 'Object'*/

void upsample_object_bilinear_3D (Real_t*** Object, Real_t*** Init, int32_t N_z, int32_t N_y, int32_t N_x)
{
	int32_t j, k, slice;
	Real_t ***buffer2D, ***buffer3D;
	
  		buffer2D = (Real_t***)multialloc(sizeof(Real_t), 3, N_z, N_y, 2*N_x);
  		buffer3D = (Real_t***)multialloc(sizeof(Real_t), 3, N_z, 2*N_y, 2*N_x);
		for (slice=0; slice < N_z; slice++){
			for (j=0; j < N_y; j++){
				buffer2D[slice][j][0] = Init[slice][j][0];
				buffer2D[slice][j][1] = (3.0*Init[slice][j][0] + Init[slice][j][1])/4.0;
				buffer2D[slice][j][2*N_x - 1] = Init[slice][j][N_x - 1];
				buffer2D[slice][j][2*N_x - 2] = (Init[slice][j][N_x - 2] + 3.0*Init[slice][j][N_x - 1])/4.0;
				for (k=1; k < N_x - 1; k++){
					buffer2D[slice][j][2*k] = (Init[slice][j][k-1] + 3.0*Init[slice][j][k])/4.0;
					buffer2D[slice][j][2*k + 1] = (3.0*Init[slice][j][k] + Init[slice][j][k+1])/4.0;
				}
			}
			for (k=0; k < 2*N_x; k++){
				buffer3D[slice][0][k] = buffer2D[slice][0][k];
				buffer3D[slice][1][k] = (3.0*buffer2D[slice][0][k] + buffer2D[slice][1][k])/4.0;
				buffer3D[slice][2*N_y-1][k] = buffer2D[slice][N_y-1][k];
				buffer3D[slice][2*N_y-2][k] = (buffer2D[slice][N_y-2][k] + 3.0*buffer2D[slice][N_y-1][k])/4.0;	
			}	
			for (j=1; j<N_y-1; j++)
				for (k=0; k<2*N_x; k++){				
					buffer3D[slice][2*j][k] = (buffer2D[slice][j-1][k] + 3.0*buffer2D[slice][j][k])/4.0;
					buffer3D[slice][2*j + 1][k] = (3*buffer2D[slice][j][k] + buffer2D[slice][j+1][k])/4.0;
				}
		}
	
		for (j=0; j<2*N_y; j++)
			for (k=0; k<2*N_x; k++){				
				Object[1][j][k] = buffer3D[0][j][k];
				Object[2][j][k] = (3.0*buffer3D[0][j][k] + buffer3D[1][j][k])/4.0;
				Object[2*N_z][j][k] = buffer3D[N_z-1][j][k];
				Object[2*N_z-1][j][k] = (3.0*buffer3D[N_z-1][j][k] + buffer3D[N_z-2][j][k])/4.0;
			}
	
		for (slice=1; slice < N_z-1; slice++)	
			for (j=0; j<2*N_y; j++)
				for (k=0; k<2*N_x; k++){				
					Object[2*slice+1][j][k] = (buffer3D[slice-1][j][k] + 3.0*buffer3D[slice][j][k])/4.0;
					Object[2*slice+2][j][k] = (3.0*buffer3D[slice][j][k] + buffer3D[slice+1][j][k])/4.0;		
				}
	
		multifree(buffer2D,3);
		multifree(buffer3D,3);

}


/*'InitObject' intializes the Object to be reconstructed to either 0 or an interpolated version of the previous reconstruction. It is used in multi resolution reconstruction in which after every coarse resolution reconstruction the object should be intialized with an interpolated version of the reconstruction following which the object will be reconstructed at a finer resolution. 
--initICD--
If 1, initializes the object to 0
If 2, the code uses bilinear interpolation to initialize the object if the previous reconstruction was at a lower resolution

The function also initializes the magnitude update map 'MagUpdateMap' from the previous coarser resolution 
reconstruction. */

void initObject(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t**** MagUpdateMap)
{
	char object_file[100];
	int dimTiff[4];
	int32_t i, j, k, l;
	Real_t ***Init, ****UpMapInit, object_init_value;
	char initfile[100] = "ObjectInit";
	char MagUpdateMapFile[100] = "mag_update_map";

	sprintf (initfile, "%s_n%d", initfile, TomoInputsPtr->node_rank);
	sprintf (MagUpdateMapFile, "%s_n%d", MagUpdateMapFile, TomoInputsPtr->node_rank);

	object_init_value = convert_HU2um(OBJECT_INIT_VAL);
	for (i = 0; i < ScannedObjectPtr->N_time; i++)
	for (j = 0; j < ScannedObjectPtr->N_z; j++)
	for (k = 0; k < ScannedObjectPtr->N_y; k++)
	for (l = 0; l < ScannedObjectPtr->N_x; l++)
		ScannedObjectPtr->Object[i][j][k][l] = 0;
	
	if (TomoInputsPtr->initICD > 3){
		fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: initObject: initICD value not recognized\n");
		exit(1);
	}
	else if (TomoInputsPtr->initICD == 1 || TomoInputsPtr->initICD == 2 || TomoInputsPtr->initICD == 3)
	{
		if (TomoInputsPtr->initICD == 1)
		{
			for (i = 0; i < ScannedObjectPtr->N_time; i++)
			{
				sprintf(object_file, "%s_n%d_time_%d", OBJECT_FILENAME,TomoInputsPtr->node_rank, i);
				Read4mBin (object_file, 1, ScannedObjectPtr->N_z, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x, &(ScannedObjectPtr->Object[i][1][0][0]), TomoInputsPtr->debug_file_ptr);
			}
			if (TomoInputsPtr->initMagUpMap == 1)
				Read4mBin (MagUpdateMapFile, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x, &(MagUpdateMap[0][0][0][0]), TomoInputsPtr->debug_file_ptr);
		}
		else 
		{
			if (TomoInputsPtr->initICD == 3)
			{
/*				Init = (Real_t****)multialloc(sizeof(Real_t), 4, ScannedObjectPtr->N_time, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2, ScannedObjectPtr->N_z/2);*/
				Init = (Real_t***)multialloc(sizeof(Real_t), 3, ScannedObjectPtr->N_z/2, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
				fprintf(TomoInputsPtr->debug_file_ptr, "initObject: Interpolating object using 3D bilinear interpolation\n");
				for (i = 0; i < ScannedObjectPtr->N_time; i++)
				{
					sprintf(object_file, "%s_n%d_time_%d", OBJECT_FILENAME, TomoInputsPtr->node_rank, i);
					Read4mBin (object_file, 1, ScannedObjectPtr->N_z/2, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2, &(Init[0][0][0]), TomoInputsPtr->debug_file_ptr);
					upsample_object_bilinear_3D (ScannedObjectPtr->Object[i], Init, ScannedObjectPtr->N_z/2, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
				}
				multifree(Init,3);
				/*multifree(Init,4);			*/

				fprintf(TomoInputsPtr->debug_file_ptr, "initObject: Done with interpolating object using 3D bilinear interpolation\n");
				if (TomoInputsPtr->initMagUpMap == 1)
				{
					UpMapInit = (Real_t****)multialloc(sizeof(Real_t), 4, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
					Read4mBin (MagUpdateMapFile, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2, &(UpMapInit[0][0][0][0]), TomoInputsPtr->debug_file_ptr);
					fprintf(TomoInputsPtr->debug_file_ptr, "initObject: Interpolating magnitude update map using 3D bilinear interpolation\n");
					upsample_bilinear_2D (MagUpdateMap, UpMapInit, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
					multifree(UpMapInit,4);	
				}
			}
			else
			{
/*				Init = (Real_t****)multialloc(sizeof(Real_t), 4, ScannedObjectPtr->N_time, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2, ScannedObjectPtr->N_z);*/
				Init = (Real_t***)multialloc(sizeof(Real_t), 3, ScannedObjectPtr->N_z, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
				for (i = 0; i < ScannedObjectPtr->N_time; i++)
				{
					sprintf(object_file, "%s_n%d_time_%d", OBJECT_FILENAME,TomoInputsPtr->node_rank,i);
					Read4mBin (object_file, 1, ScannedObjectPtr->N_z, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2, &(Init[0][0][0]), TomoInputsPtr->debug_file_ptr);
					upsample_object_bilinear_2D (ScannedObjectPtr->Object[i], Init, ScannedObjectPtr->N_z, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
				}
				multifree(Init,3);
/*				multifree(Init,4);*/
				
				if (TomoInputsPtr->initMagUpMap == 1)
				{
					UpMapInit = (Real_t****)multialloc(sizeof(Real_t), 4, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
					Read4mBin (MagUpdateMapFile, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2, &(UpMapInit[0][0][0][0]), TomoInputsPtr->debug_file_ptr);
					upsample_bilinear_2D (MagUpdateMap, UpMapInit, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
					multifree(UpMapInit,4);	
				}
			}
		}
	}

	dimTiff[0] = ScannedObjectPtr->N_time; dimTiff[1] = TomoInputsPtr->num_z_blocks; dimTiff[2] = ScannedObjectPtr->N_y; dimTiff[3] = ScannedObjectPtr->N_x;
	sprintf(object_file, "%s_Init", MagUpdateMapFile);
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (object_file, dimTiff, 0, 1, 2, 3, &(MagUpdateMap[0][0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	
	for (i = 0; i < ScannedObjectPtr->N_time; i++)
	{
		sprintf (object_file, "%s_time_%d", initfile, i);
		dimTiff[0] = 1; dimTiff[1] = ScannedObjectPtr->N_z; dimTiff[2] = ScannedObjectPtr->N_y; dimTiff[3] = ScannedObjectPtr->N_x;
		if (TomoInputsPtr->Write2Tiff == 1)
			WriteMultiDimArray2Tiff (object_file, dimTiff, 0, 1, 2, 3, &(ScannedObjectPtr->Object[i][1][0][0]), 1, TomoInputsPtr->debug_file_ptr);
	}

}	

/*'initErrorSinogram' is used to initialize the error sinogram before start of ICD. It computes e = y - Ax - d. Ax is computed by forward projecting the obkject x.*/
void initErrorSinogam (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t** DetectorResponse, Real_t*** ErrorSino, AMatrixCol* VoxelLineResponse)
{
	Real_t pixel, val, avg=0;
	int32_t dimTiff[4], i, j, k, p, m, n, sino_idx, idx, slice;
	AMatrixCol* AMatrixPtr = (AMatrixCol*)get_spc(ScannedObjectPtr->N_time, sizeof(AMatrixCol));
  	uint8_t AvgNumXElements = (uint8_t)ceil(3*ScannedObjectPtr->delta_xy/SinogramPtr->delta_r);
	char error_file[100] = "error_sinogram";	

	sprintf(error_file, "%s_n%d", error_file, TomoInputsPtr->node_rank);
	for (i = 0; i < ScannedObjectPtr->N_time; i++)
	{
		AMatrixPtr[i].values = (Real_t*)get_spc(AvgNumXElements, sizeof(Real_t));
		AMatrixPtr[i].index = (int32_t*)get_spc(AvgNumXElements, sizeof(int32_t));
	}

	memset(&(ErrorSino[0][0][0]), 0, SinogramPtr->N_p*SinogramPtr->N_t*SinogramPtr->N_r*sizeof(Real_t));	

	#pragma omp parallel for private(j, k, p, sino_idx, slice, pixel, idx, val, m, n)
	for (i=0; i<ScannedObjectPtr->N_time; i++)
	{
		for (j=0; j<ScannedObjectPtr->N_y; j++)
		{
			for (k=0; k<ScannedObjectPtr->N_x; k++){	
	    			for (p=0; p<ScannedObjectPtr->ProjNum[i]; p++){
	    				sino_idx = ScannedObjectPtr->ProjIdxPtr[i][p];
	    				calcAMatrixColumnforAngle(SinogramPtr, ScannedObjectPtr, DetectorResponse, &(AMatrixPtr[i]), j, k, sino_idx); 
	    				for (slice=0; slice<ScannedObjectPtr->N_z; slice++){
/*		printf("count = %d, idx = %d, val = %f\n",  VoxelLineResponse[slice].count, VoxelLineResponse[slice].index[0], VoxelLineResponse[slice].values[0]);*/
	    					pixel = ScannedObjectPtr->Object[i][slice+1][j][k]; /*slice+1 to account for extra z slices required for MPI*/	
	    					for (m=0; m<AMatrixPtr[i].count; m++){
							idx=AMatrixPtr[i].index[m];
							val=AMatrixPtr[i].values[m];
							for (n=0; n<VoxelLineResponse[slice].count; n++){
								ErrorSino[sino_idx][idx][VoxelLineResponse[slice].index[n]] += pixel*val*VoxelLineResponse[slice].values[n];
							}
	        				}
	   				}	
	   			}
			}
		}
	}
	
        #pragma omp parallel for private(j, k) reduction(+:avg)
	for(i=0; i < SinogramPtr->N_p; i++)
    	for (j = 0; j < SinogramPtr->N_r; j++)
    	for (k = 0; k < SinogramPtr->N_t; k++)
	{
			ErrorSino[i][j][k]=SinogramPtr->Projection[i][j][k] - ErrorSino[i][j][k] - SinogramPtr->ProjOffset[j][k];
	   		if (fabs(ErrorSino[i][j][k]*sqrt(TomoInputsPtr->Weight[i][j][k])) < TomoInputsPtr->ErrorSinoThresh)
				SinogramPtr->ProjSelect[i][j][k] = true;
			else
				SinogramPtr->ProjSelect[i][j][k] = false;
/*			if (ErrorSino[i][j][k]*sqrt(TomoInputsPtr->Weight[i][j][k]) < -30)
				TomoInputsPtr->Weight[i][j][k] = 0;*/
			avg+=ErrorSino[i][j][k];
	}
	avg=avg/(SinogramPtr->N_r*SinogramPtr->N_t*SinogramPtr->N_p);
	fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Average of Error Sinogram is %f\n", avg);
	
	dimTiff[0] = 1; dimTiff[1] = SinogramPtr->N_p; dimTiff[2] = SinogramPtr->N_r; dimTiff[3] = SinogramPtr->N_t;
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (error_file, dimTiff, 0, 3, 1, 2, &(ErrorSino[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);

	for (i = 0; i < ScannedObjectPtr->N_time; i++)
	{
		free(AMatrixPtr[i].values);
		free(AMatrixPtr[i].index);
	}
	free (AMatrixPtr);
}

/*randomly select the voxels lines which need to be updated along the x-y plane for each z-block and time slice*/
void randomly_select_x_y (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, uint8_t**** Mask)
{
        int32_t i, j, num,n, Index, col, row, *Counter, ArraySize, block;

        ArraySize = ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
        Counter = (int32_t*)get_spc(ArraySize, sizeof(int32_t));

        for (i=0; i<ScannedObjectPtr->N_time; i++)
        for (block=0; block<TomoInputsPtr->num_z_blocks; block++)
        {   
                ArraySize = ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
                for (Index = 0; Index < ArraySize; Index++)
                        Counter[Index] = Index;
   
		TomoInputsPtr->UpdateSelectNum[i][block] = 0; 
                for (j=0; j<ScannedObjectPtr->N_x*ScannedObjectPtr->N_y; j++){
                        Index = floor(random2() * ArraySize);
                        Index = (Index == ArraySize)?ArraySize-1:Index;
			col = Counter[Index] % ScannedObjectPtr->N_x;
			row = Counter[Index] / ScannedObjectPtr->N_x;
			for (n = block*(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks); n < (block+1)*(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks); n++)
				if (Mask[i][row][col][n] == 1)
				{
					num = TomoInputsPtr->UpdateSelectNum[i][block];
                        		TomoInputsPtr->x_rand_select[i][block][num] = col;
                        		TomoInputsPtr->y_rand_select[i][block][num] = row;
					(TomoInputsPtr->UpdateSelectNum[i][block])++;
					break;
				}
                        Counter[Index] = Counter[ArraySize - 1]; 
                        ArraySize--;
                }   
        }   

        free(Counter);

} 

/*Updates the variance parameter \sigma*/
void update_variance_parameter (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino)
{
	int32_t k, i, j;
	Real_t temp_acc = 0, temp = 0;
	
  	#pragma omp parallel for private(i, j, temp) reduction(+:temp_acc)
	for (k = 0; k < SinogramPtr->N_p; k++)
	for (i = 0; i < SinogramPtr->N_r; i++)
	for (j = 0; j < SinogramPtr->N_t; j++)
	{
		TomoInputsPtr->Weight[k][i][j] = TomoInputsPtr->Weight[k][i][j]*TomoInputsPtr->var_est;
		if (SinogramPtr->ProjSelect[k][i][j] == true)
			temp = ErrorSino[k][i][j]*ErrorSino[k][i][j]*TomoInputsPtr->Weight[k][i][j];
		else
			temp = fabs(ErrorSino[k][i][j])*TomoInputsPtr->ErrorSinoDelta*TomoInputsPtr->ErrorSinoThresh*sqrt(TomoInputsPtr->Weight[k][i][j]*TomoInputsPtr->var_est);
		temp_acc += temp;
	}
 
	MPI_Allreduce(&temp_acc, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	TomoInputsPtr->var_est = temp/(TomoInputsPtr->node_num*SinogramPtr->N_p*SinogramPtr->N_r*SinogramPtr->N_t);

  	#pragma omp parallel for private(i, j)
	for (k = 0; k < SinogramPtr->N_p; k++)
	for (i = 0; i < SinogramPtr->N_r; i++)
	for (j = 0; j < SinogramPtr->N_t; j++)
	{	
		TomoInputsPtr->Weight[k][i][j] /= TomoInputsPtr->var_est;
	   	if (fabs(ErrorSino[k][i][j]*sqrt(TomoInputsPtr->Weight[k][i][j])) < TomoInputsPtr->ErrorSinoThresh)
			SinogramPtr->ProjSelect[k][i][j] = true;
		else
			SinogramPtr->ProjSelect[k][i][j] = false;
	}

} 

/*Updates the projection offset error parameter d_i*/
void update_Sinogram_Offset (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino)
{
/*	int32_t i, j, k;

	fprintf(TomoInputsPtr->debug_file_ptr, "update_Sinogram_Offset: Updating offset error in projection\n");*/
/*	memset(&(num[0][0]), 0, SinogramPtr->N_r*SinogramPtr->N_t*sizeof(Real_t));
	for (k = 0; k < SinogramPtr->N_p; k++)
	for (i = 0; i < SinogramPtr->N_r; i++)
	for (j = 0; j < SinogramPtr->N_t; j++)
	{
		ErrorSino[k][i][j] += SinogramPtr->ProjOffset[i][j];
		num[i][j] += ErrorSino[k][i][j]*TomoInputsPtr->Weight[k][i][j];*/
/*		den[i][j] += TomoInputsPtr->Weight[k][i][j];*/
/*	}
			
	for (i = 0; i < SinogramPtr->N_r; i++)
	for (j = 0; j < SinogramPtr->N_t; j++)
		SinogramPtr->ProjOffset[i][j] = num[i][j]/den[i][j];
		
	for (k = 0; k < SinogramPtr->N_p; k++)
	for (i = 0; i < SinogramPtr->N_r; i++)
	for (j = 0; j < SinogramPtr->N_t; j++)
		ErrorSino[k][i][j] -= SinogramPtr->ProjOffset[i][j]; */

	Real_t numerator, temp, denominator;
	int32_t i, j, k;


	#pragma omp parallel for private(j, k, numerator, denominator, temp)
	for (i = 0; i < SinogramPtr->N_r; i++)
	for (j = 0; j < SinogramPtr->N_t; j++)
	{
		numerator = 0;
		denominator = 0;
		for (k = 0; k < SinogramPtr->N_p; k++)
		{
			temp = TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoDelta*sqrt(TomoInputsPtr->Weight[k][i][j])/fabs(ErrorSino[k][i][j]);
			ErrorSino[k][i][j] += SinogramPtr->ProjOffset[i][j];
			if (SinogramPtr->ProjSelect[k][i][j] == true)
			{
				numerator += ErrorSino[k][i][j]*TomoInputsPtr->Weight[k][i][j];
				denominator += TomoInputsPtr->Weight[k][i][j];
			}
			else
			{
				numerator += temp*ErrorSino[k][i][j];
				denominator += temp;
			} 
		}
			
		/*if (denominator != 0)	*/
		SinogramPtr->ProjOffset[i][j] = numerator/denominator;
		/*else
		{
			int32_t dimTiff[4];
			char projselect_file[100] = "ProjSelectErred";
			sprintf(projselect_file, "%s_n%d", projselect_file, TomoInputsPtr->node_rank);
			dimTiff[0] = 1; dimTiff[1] = SinogramPtr->N_p; dimTiff[2] = SinogramPtr->N_r; dimTiff[3] = SinogramPtr->N_t;
			WriteBoolArray2Tiff (projselect_file, dimTiff, 0, 3, 1, 2, &(SinogramPtr->ProjSelect[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
			printf("ERROR: denominator is zero!!!!\n");
			exit(1);
		}*/
			
		for (k = 0; k < SinogramPtr->N_p; k++)
		{
			ErrorSino[k][i][j] -= SinogramPtr->ProjOffset[i][j]; 
	   		if (fabs(ErrorSino[k][i][j]*sqrt(TomoInputsPtr->Weight[k][i][j])) < TomoInputsPtr->ErrorSinoThresh)
				SinogramPtr->ProjSelect[k][i][j] = true;
			else
				SinogramPtr->ProjSelect[k][i][j] = false;
		}
	}
}

/*Implements mutithreaded shared memory parallelization using OpenMP and splits work among
threads. Each thread gets a certain time slice and z block to update.
Multithreading is done within the z-blocks assigned to each node.
ErrorSino - Error sinogram
Iter - Present iteration number
MagUpdateMap - Magnitude update map containing the magnitude of update of each voxel
Mask - If a certain element is true then the corresponding voxel is updated*/
int updateVoxelsTimeSlices(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr,  Real_t** DetectorResponse, AMatrixCol* VoxelLineResponse, Real_t*** ErrorSino, int32_t Iter, Real_t**** MagUpdateMap, uint8_t**** Mask)
{
	Real_t AverageUpdate = 0, tempUpdate, avg_update_percentage, total_vox_mag = 0.0, vox_mag = 0.0;
	int32_t xy_start, xy_end, i, j, K, block, idx, **z_start, **z_stop;
	Real_t tempTotPix = 0, total_pix = 0;
	long int **zero_count, total_zero_count = 0;	
	int32_t** thread_num = (int32_t**)multialloc(sizeof(int32_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
	MPI_Request *send_reqs, *recv_reqs;

	send_reqs = (MPI_Request*)get_spc(ScannedObjectPtr->N_time, sizeof(MPI_Request)); 
	recv_reqs = (MPI_Request*)get_spc(ScannedObjectPtr->N_time, sizeof(MPI_Request)); 
	z_start = (int32_t**)multialloc(sizeof(int32_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
	z_stop = (int32_t**)multialloc(sizeof(int32_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
    
	randomly_select_x_y (ScannedObjectPtr, TomoInputsPtr, Mask);
    
	zero_count = (long int**)multialloc(sizeof(Real_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
/*	offset_numerator = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
	memset(&(offset_denominator[0][0]), 0, SinogramPtr->N_r*SinogramPtr->N_t*sizeof(Real_t));	
    
	for (k = 0; k < SinogramPtr->N_p; k++)
	for (i = 0; i < SinogramPtr->N_r; i++)
	for (j = 0; j < SinogramPtr->N_t; j++)
		offset_denominator[i][j] += TomoInputsPtr->Weight[k][i][j]; */
 
	memset(&(zero_count[0][0]), 0, ScannedObjectPtr->N_time*TomoInputsPtr->num_z_blocks*sizeof(long int));	
/*	K = ScannedObjectPtr->N_time*ScannedObjectPtr->N_z*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
	K = (K - total_zero_count)/(ScannedObjectPtr->gamma*K);*/
	K = floor(1.0/ScannedObjectPtr->gamma);
	fprintf(TomoInputsPtr->debug_file_ptr, "Number of NHICD iterations is %d\n", K);
	for (j = 0; j < K; j++)
	{
		total_vox_mag = 0.0;
		#pragma omp parallel for collapse(2) private(i, block, idx, xy_start, xy_end) reduction(+:total_vox_mag)
		for (i = 0; i < ScannedObjectPtr->N_time; i++)
		for (block = 0; block < TomoInputsPtr->num_z_blocks; block = block + 2)
		{
			idx = (i % 2 == 0) ? block: block + 1;		
			z_start[i][idx] = idx*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);	
			z_stop[i][idx] = (idx + 1)*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks) - 1;	
			z_stop[i][idx] = (idx >= TomoInputsPtr->num_z_blocks - 1) ? ScannedObjectPtr->N_z - 1: z_stop[i][idx]; 
			xy_start = j*floor(TomoInputsPtr->UpdateSelectNum[i][idx]/K);
			xy_end = (j + 1)*floor(TomoInputsPtr->UpdateSelectNum[i][idx]/K) - 1;
			xy_end = (j == K - 1) ? TomoInputsPtr->UpdateSelectNum[i][idx] - 1: xy_end;
			total_vox_mag += updateVoxels (i, i, z_start[i][idx], z_stop[i][idx], xy_start, xy_end, TomoInputsPtr->x_rand_select[i][idx], TomoInputsPtr->y_rand_select[i][idx], SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse, VoxelLineResponse, Iter, &(zero_count[i][idx]), MagUpdateMap[i][idx], Mask[i]);
			thread_num[i][idx] = omp_get_thread_num();
/*			printf ("i = %d, idx = %d, z_start = %d, z_stop = %d, xy_start = %d, xy_end = %d\n", i, idx, z_start, z_stop, xy_start, xy_end);*/
		}
		
		MPI_Send_Recv_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 0);
		if (TomoInputsPtr->updateProjOffset > 1)
			update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
		MPI_Wait_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 0);

		#pragma omp parallel for collapse(2) private(i, block, idx, xy_start, xy_end) reduction(+:total_vox_mag)
		for (i = 0; i < ScannedObjectPtr->N_time; i++)
		for (block = 0; block < TomoInputsPtr->num_z_blocks; block = block + 2)
		{			
			idx = (i % 2 == 0) ? block + 1: block;		
			z_start[i][idx] = idx*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);	
			z_stop[i][idx] = (idx + 1)*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks) - 1;	
			z_stop[i][idx] = (idx >= TomoInputsPtr->num_z_blocks - 1) ? ScannedObjectPtr->N_z - 1: z_stop[i][idx]; 
			xy_start = j*floor(TomoInputsPtr->UpdateSelectNum[i][idx]/K);
			xy_end = (j + 1)*floor(TomoInputsPtr->UpdateSelectNum[i][idx]/K) - 1;
			xy_end = (j == K - 1) ? TomoInputsPtr->UpdateSelectNum[i][idx] - 1: xy_end;
			total_vox_mag += updateVoxels (i, i, z_start[i][idx], z_stop[i][idx], xy_start, xy_end, TomoInputsPtr->x_rand_select[i][idx], TomoInputsPtr->y_rand_select[i][idx], SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse, VoxelLineResponse, Iter, &(zero_count[i][idx]), MagUpdateMap[i][idx], Mask[i]);
			thread_num[i][idx] = omp_get_thread_num();
			/*printf ("i = %d, idx = %d, z_start = %d, z_stop = %d, xy_start = %d, xy_end = %d\n", i, idx, z_start, z_stop, xy_start, xy_end);*/
		}
        
		MPI_Send_Recv_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 1);
		if (TomoInputsPtr->updateProjOffset > 1)
			update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
		MPI_Wait_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 1);

		VSC_based_Voxel_Line_Select(ScannedObjectPtr, TomoInputsPtr, MagUpdateMap);			
/*		fprintf(TomoInputsPtr->debug_file_ptr, "Number of NHICD voxel lines to be updated in iteration %d is %d\n", j, num_voxel_lines);*/
		if (Iter > 1 && TomoInputsPtr->no_NHICD == 0)
		{
			#pragma omp parallel for collapse(2) private(i, block, idx)
			for (i = 0; i < ScannedObjectPtr->N_time; i++)
			for (block = 0; block < TomoInputsPtr->num_z_blocks; block = block + 2)
			{				
				idx = (i % 2 == 0) ? block: block + 1;		
				z_start[i][idx] = idx*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);	
				z_stop[i][idx] = (idx + 1)*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks) - 1;	
				z_stop[i][idx] = (idx >= TomoInputsPtr->num_z_blocks - 1) ? ScannedObjectPtr->N_z - 1: z_stop[i][idx]; 
				updateVoxels (i, i, z_start[i][idx], z_stop[i][idx], 0, TomoInputsPtr->NHICDSelectNum[i][idx]-1, TomoInputsPtr->x_NHICD_select[i][idx], TomoInputsPtr->y_NHICD_select[i][idx], SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse, VoxelLineResponse, Iter, &(zero_count[i][idx]), MagUpdateMap[i][idx], Mask[i]);
				thread_num[i][idx] = omp_get_thread_num();
				/*printf ("i = %d, idx = %d, z_start = %d, z_stop = %d, xy_start = %d, xy_end = %d\n", i, idx, z_start, z_stop, xy_start, xy_end);*/
			}
			
			MPI_Send_Recv_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 0);
			if (TomoInputsPtr->updateProjOffset > 1)
				update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
			MPI_Wait_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 0);
			
			#pragma omp parallel for collapse(2) private(i, block, idx)
			for (i = 0; i < ScannedObjectPtr->N_time; i++)
			for (block = 0; block < TomoInputsPtr->num_z_blocks; block = block + 2)
			{				
				idx = (i % 2 == 0) ? block + 1: block;		
				z_start[i][idx] = idx*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);	
				z_stop[i][idx] = (idx + 1)*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks) - 1;	
				z_stop[i][idx] = (idx >= TomoInputsPtr->num_z_blocks - 1) ? ScannedObjectPtr->N_z - 1: z_stop[i][idx]; 
				updateVoxels (i, i, z_start[i][idx], z_stop[i][idx], 0, TomoInputsPtr->NHICDSelectNum[i][idx]-1, TomoInputsPtr->x_NHICD_select[i][idx], TomoInputsPtr->y_NHICD_select[i][idx], SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse, VoxelLineResponse, Iter, &(zero_count[i][idx]), MagUpdateMap[i][idx], Mask[i]);
				thread_num[i][idx] = omp_get_thread_num();
				/*printf ("i = %d, idx = %d, z_start = %d, z_stop = %d, xy_start = %d, xy_end = %d\n", i, idx, z_start, z_stop, xy_start, xy_end);*/
			}
	
			MPI_Send_Recv_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 1);
			if (TomoInputsPtr->updateProjOffset > 1)
				update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
			MPI_Wait_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 1);
		}
	}
 
	if (TomoInputsPtr->updateVar == 1)
		update_variance_parameter (SinogramPtr, TomoInputsPtr, ErrorSino);
        
	fprintf(TomoInputsPtr->debug_file_ptr, "updateVoxelsTimeSlices: Time Slice, Z Start, Z End - Thread : ");
  	total_pix = 0;
        for (i=0; i<ScannedObjectPtr->N_time; i++){
        	for (block=0; block<TomoInputsPtr->num_z_blocks; block++){
			total_pix += TomoInputsPtr->UpdateSelectNum[i][block]*(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);
        		for (j=0; j<TomoInputsPtr->UpdateSelectNum[i][block]; j++){
                		AverageUpdate += MagUpdateMap[i][block][TomoInputsPtr->y_rand_select[i][block][j]][TomoInputsPtr->x_rand_select[i][block][j]];
			}
                	total_zero_count += zero_count[i][block];
               	 	fprintf(TomoInputsPtr->debug_file_ptr, "%d,%d,%d-%d; ", i, z_start[i][block], z_stop[i][block], thread_num[i][block]);
		}
	}
  

	MPI_Allreduce(&AverageUpdate, &tempUpdate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&total_pix, &tempTotPix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&total_vox_mag, &vox_mag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        AverageUpdate = tempUpdate/(tempTotPix);
        AverageUpdate = convert2Hounsfield(AverageUpdate);

        fprintf(TomoInputsPtr->debug_file_ptr, "\nupdateVoxelsTimeSlices: Average voxel update over all voxels is %f, total voxels is %f\n", AverageUpdate, tempTotPix);
        fprintf(TomoInputsPtr->debug_file_ptr, "updateVoxelsTimeSlices: Zero count is %ld\n", total_zero_count);
	
	multifree(zero_count,2);
	multifree(thread_num,2);
	multifree(z_start,2);
	multifree(z_stop,2);
	free(send_reqs);
	free(recv_reqs);
/*	multifree(offset_numerator,2);
	multifree(offset_denominator,2);*/
        avg_update_percentage = 100*tempUpdate/vox_mag;
        fprintf(TomoInputsPtr->debug_file_ptr, "updateVoxelsTimeSlices: percentage change in average magnitude is %f\n", avg_update_percentage);
	
	if (avg_update_percentage < TomoInputsPtr->StopThreshold)
        {   
                fprintf(TomoInputsPtr->debug_file_ptr, "updateVoxelsTimeSlices: Reached average magnitude of change in voxels threshold\n");
		return (1);
	}
	return(0);
}



	/*ICD_BackProject calls the ICD optimization function repeatedly till the stopping criteria is met.*/
int ICD_BackProject(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
#ifndef NO_COST_CALCULATE
	Real_t cost, cost_0_iter, cost_last_iter, percentage_change_in_cost = 0;
	char costfile[100]="cost";
#endif
	Real_t ***ErrorSino, x, y, **H_r, *H_t;
	int32_t j, flag, Iter, i, k, l;
	int dimTiff[4];
	char MagUpdateMapFile[100] = "mag_update_map";
	char VarEstFile[100] = "variance_estimate";
	char scaled_error_file[100] = "scaled_errorsino";
	time_t start;
	char detect_file[100]="DetectorResponse";
	char projselect_file[100] = "ProjSelect";
	uint8_t  ****Mask, AvgNumZElements;
	AMatrixCol *VoxelLineResponse;

#ifdef POSITIVITY_CONSTRAINT
	fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Enforcing positivity constraint\n");
#endif

	sprintf(MagUpdateMapFile, "%s_n%d", MagUpdateMapFile, TomoInputsPtr->node_rank);
	sprintf(scaled_error_file, "%s_n%d", scaled_error_file, TomoInputsPtr->node_rank);
	sprintf(detect_file, "%s_n%d", detect_file, TomoInputsPtr->node_rank);
	sprintf(projselect_file, "%s_n%d", projselect_file, TomoInputsPtr->node_rank);

	Real_t**** MagUpdateMap = (Real_t****)multialloc(sizeof(Real_t), 4, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x);
	memset(&(MagUpdateMap[0][0][0][0]), 0, ScannedObjectPtr->N_time*TomoInputsPtr->num_z_blocks*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x*sizeof(Real_t));	
	
	omp_set_num_threads(TomoInputsPtr->num_threads);
	fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Number of CPU cores is %d\n", (int)omp_get_num_procs());
	/*	fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Number of threads is %d\n", TomoInputsPtr->num_threads) ;*/

	Mask = (uint8_t****)multialloc(sizeof(uint8_t), 4, ScannedObjectPtr->N_time, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x, ScannedObjectPtr->N_z);
	for (i = 0; i < ScannedObjectPtr->N_time; i++)
		for (j = 0; j < ScannedObjectPtr->N_y; j++)
			for (k = 0; k < ScannedObjectPtr->N_x; k++){
				x = ScannedObjectPtr->x0 + ((Real_t)k + 0.5)*ScannedObjectPtr->delta_xy;
				y = ScannedObjectPtr->y0 + ((Real_t)j + 0.5)*ScannedObjectPtr->delta_xy;
				if (x*x + y*y < TomoInputsPtr->radius_obj*TomoInputsPtr->radius_obj)
				{
					for (l = 0; l < ScannedObjectPtr->N_z; l++)
						Mask[i][j][k][l] = 1;
				}
				else
				{
					for (l = 0; l < ScannedObjectPtr->N_z; l++)
						Mask[i][j][k][l] = 0;
				}
			}   
		
	H_r = (Real_t **)multialloc(sizeof(Real_t), 2, SinogramPtr->N_p, DETECTOR_RESPONSE_BINS + 1);
	H_t = (Real_t *)get_spc(DETECTOR_RESPONSE_BINS + 1, sizeof(Real_t));
	DetectorResponseProfile (H_r, H_t, SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
		
	VoxelLineResponse = (AMatrixCol*)get_spc(ScannedObjectPtr->N_z, sizeof(AMatrixCol));
	AvgNumZElements = (uint8_t)((ScannedObjectPtr->delta_z/SinogramPtr->delta_t) + 2);
	for (j = 0; j < ScannedObjectPtr->N_z; j++) 
	{
		VoxelLineResponse[j].count = 0; 
		VoxelLineResponse[j].values = (Real_t*)get_spc(AvgNumZElements, sizeof(Real_t));
		VoxelLineResponse[j].index = (int32_t*)get_spc(AvgNumZElements, sizeof(int32_t));
	}

	storeVoxelLineResponse(H_t, VoxelLineResponse, ScannedObjectPtr, SinogramPtr);

	dimTiff[0] = 1; dimTiff[1] = 1; dimTiff[2] = SinogramPtr->N_p; dimTiff[3] = DETECTOR_RESPONSE_BINS+1;
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (detect_file, dimTiff, 0, 1, 2, 3, &(H_r[0][0]), 0, TomoInputsPtr->debug_file_ptr);

	ErrorSino = (Real_t***)multialloc(sizeof(Real_t), 3, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t);
	start=time(NULL);
	initObject(SinogramPtr, ScannedObjectPtr, TomoInputsPtr, MagUpdateMap);
	fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Time taken to read object = %fmins\n", difftime(time(NULL),start)/60.0);

	if (TomoInputsPtr->only_Edge_Updates == 1)
		select_edge_voxels (ScannedObjectPtr, TomoInputsPtr, MagUpdateMap, Mask);

	initErrorSinogam(SinogramPtr, ScannedObjectPtr, TomoInputsPtr, H_r, ErrorSino, VoxelLineResponse);
	multifree(SinogramPtr->Projection,3);
	fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Time taken to initialize object and compute error sinogram = %fmins\n", difftime(time(NULL),start)/60.0);

#ifndef NO_COST_CALCULATE
	sprintf(costfile, "%s_n%d", costfile, TomoInputsPtr->node_rank);
	cost=computeCost(SinogramPtr,ScannedObjectPtr,TomoInputsPtr,ErrorSino);
	cost_0_iter = cost; 
	cost_last_iter = cost;
	fprintf(TomoInputsPtr->debug_file_ptr, "-------------ICD_BackProject: ICD Iter=Before ICD, cost=%f------------\n",cost);
	Write2Bin (costfile, 1, 1, 1, 1, &cost, TomoInputsPtr->debug_file_ptr);
#endif /*Cost calculation endif*/
		
	start=time(NULL);
	for (Iter = 1; Iter <= TomoInputsPtr->NumIter; Iter++)
	{	
		flag = updateVoxelsTimeSlices (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, H_r, VoxelLineResponse, ErrorSino, Iter, MagUpdateMap, Mask);
		if (TomoInputsPtr->WritePerIter == 1)
			write_ObjectProjOff2TiffBinPerIter (SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
#ifndef NO_COST_CALCULATE
		cost = computeCost(SinogramPtr,ScannedObjectPtr,TomoInputsPtr,ErrorSino);
		percentage_change_in_cost = ((cost - cost_last_iter)/(cost - cost_0_iter))*100.0;
		fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Percentage change in cost is %f\n", percentage_change_in_cost);
		fprintf(TomoInputsPtr->debug_file_ptr, "-------------ICD_BackProject: ICD Iter=%d, cost=%f, time since start of ICD = %fmins------------\n",Iter,cost,difftime(time(NULL),start)/60.0);
		Append2Bin (costfile, 1, 1, 1, 1, &cost, TomoInputsPtr->debug_file_ptr);
		if(cost > cost_last_iter){
			printf("ERROR: ICD_BackProject: Cost value increased\n");
			fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: ICD_BackProject: Cost value increased\n");
			exit(1);
		}
		cost_last_iter = cost;
		/*if (percentage_change_in_cost < TomoInputsPtr->cost_thresh && flag != 0 && Iter > 1){*/
		if (flag != 0 && Iter > 1){
			 fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Convergence criteria is met\n");
			 break;
		}
#else
		fprintf(TomoInputsPtr->debug_file_ptr, "-------------ICD_BackProject: ICD Iter=%d, time since start of ICD = %fmins------------\n",Iter,difftime(time(NULL),start)/60.0);
		if (flag != 0 && Iter > 1){
			fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Convergence criteria is met\n");
			break;
		}
#endif
/*		if (TomoInputsPtr->updateProjOffset > 0)
		{	
			update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
#ifndef NO_COST_CALCULATE
			Append2Bin (costfile, 1, 1, 1, 1, &cost, TomoInputsPtr->debug_file_ptr);
			cost = computeCost(SinogramPtr,ScannedObjectPtr,TomoInputsPtr,ErrorSino);
	    		fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Cost after projection offset update is %f\n", cost);
			if(cost > cost_last_iter){
				printf("ERROR: ICD_BackProject: Cost value increased\n");
				fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: ICD_BackProject: Cost value increased\n");
				exit(1);
			}
			cost_last_iter = cost;
#endif
		}*/
		flag = fflush(TomoInputsPtr->debug_file_ptr);
		if (flag != 0 )
			fprintf (TomoInputsPtr->debug_file_ptr, "WARNING: Cannot flush buffer for debug.log\n");
	}

	fprintf (TomoInputsPtr->debug_file_ptr, "ICD_BackProject: The estimated variance is %f\n", TomoInputsPtr->var_est);
	sprintf(VarEstFile, "%s_n%d", VarEstFile, TomoInputsPtr->node_rank);
	Write2Bin (VarEstFile, 1, 1, 1, 1, &(TomoInputsPtr->var_est), TomoInputsPtr->debug_file_ptr);
	Write2Bin (MagUpdateMapFile, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x, &(MagUpdateMap[0][0][0][0]), TomoInputsPtr->debug_file_ptr);
	dimTiff[0] = 1; dimTiff[1] = SinogramPtr->N_p; dimTiff[2] = SinogramPtr->N_r; dimTiff[3] = SinogramPtr->N_t;
	for (i = 0; i < SinogramPtr->N_p; i++)
	for (j = 0; j < SinogramPtr->N_r; j++)
	for (k = 0; k < SinogramPtr->N_t; k++)
		ErrorSino[i][j][k] *= sqrt(TomoInputsPtr->Weight[i][j][k]);

	if (TomoInputsPtr->Write2Tiff == 1)
	{
		WriteMultiDimArray2Tiff (scaled_error_file, dimTiff, 0, 3, 1, 2, &(ErrorSino[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
		WriteBoolArray2Tiff (projselect_file, dimTiff, 0, 3, 1, 2, &(SinogramPtr->ProjSelect[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	}
		
	multifree(ErrorSino,3);
	multifree(H_r,2);
	free(H_t);

	for (j = 0; j < ScannedObjectPtr->N_z; j++) 
	{
		free(VoxelLineResponse[j].values);
		free(VoxelLineResponse[j].index);
	}
	free(VoxelLineResponse);

	fprintf(TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Finished running ICD_BackProject\n");
	multifree(Mask,4);
	multifree(MagUpdateMap, 4);
	return(0);
}

