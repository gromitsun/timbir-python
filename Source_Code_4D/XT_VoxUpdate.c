
#include <stdio.h>
#include "XT_Structures.h"
#include "XT_ICD_update.h"
#include "XT_AMatrix.h"

void compute_voxel_update_AMat1D (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new)
{
  	int32_t p, q, r, sino_view;
	Real_t V,THETA1,THETA2,THETASelTemp;
	Real_t UpdatedVoxelValue, ProjectionEntry;
  	int32_t i_r, i_t;
        V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]; /*Store the present value of the voxel*/
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
}



void compute_voxel_update_AMat2D (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new)
{
  	int32_t p, q, r, sino_view;
	Real_t V,THETA1,THETA2,THETASelTemp,***AMatrix2D;
	Real_t UpdatedVoxelValue, ProjectionEntry;
  	int32_t i_r, i_t;
	int32_t *r_ax_start, *r_ax_num, *t_ax_start, *t_ax_num;
        V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]; /*Store the present value of the voxel*/

	r_ax_start = (int32_t*)get_spc(ScannedObjectPtr->ProjNum[i_new], sizeof(int32_t));
	r_ax_num = (int32_t*)get_spc(ScannedObjectPtr->ProjNum[i_new], sizeof(int32_t));
	t_ax_start = (int32_t*)get_spc(ScannedObjectPtr->ProjNum[i_new], sizeof(int32_t));
	t_ax_num = (int32_t*)get_spc(ScannedObjectPtr->ProjNum[i_new], sizeof(int32_t));
	AMatrix2D = (Real_t***)get_spc(ScannedObjectPtr->ProjNum[i_new], sizeof(Real_t**));

	THETA1 = 0.0;
	THETA2 = 0.0;
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
	sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
	compute_2DAMatrix_4m_1D(&(AMatrix2D[p]), &(AMatrixPtr[p]), &(VoxelLineResponse[slice]), &(r_ax_start[p]), &(r_ax_num[p]), &(t_ax_start[p]), &(t_ax_num[p]));
	for (q = 0; q < r_ax_num[p]; q++)
	{
      	    	i_r = r_ax_start[p] + q;
		for (r = 0; r < t_ax_num[p]; r++)
		{ 
       	    		ProjectionEntry = (AMatrix2D[p][q][r]);
			i_t = t_ax_start[p] + r;
			if (SinogramPtr->ProjSelect[sino_view][i_r][i_t] == true)
			{
	           		THETA2 += (ProjectionEntry*ProjectionEntry*TomoInputsPtr->Weight[sino_view][i_r][i_t]);
               			THETA1 += -(ErrorSino[sino_view][i_r][i_t]*ProjectionEntry*TomoInputsPtr->Weight[sino_view][i_r][i_t]);
            		}
			else
			{
				THETASelTemp = TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoDelta*sqrt(TomoInputsPtr->Weight[sino_view][i_r][i_t])/fabs(ErrorSino[sino_view][i_r][i_t]);
	            		THETA2 += (ProjectionEntry*ProjectionEntry*THETASelTemp);
            			THETA1 += -(ErrorSino[sino_view][i_r][i_t]*ProjectionEntry*THETASelTemp);
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
		for (q = 0; q < r_ax_num[p]; q++)
        	{
               	    	i_r = r_ax_start[p] + q;
			for (r = 0; r < t_ax_num[p]; r++)
			{ 
				i_t = t_ax_start[p] + r;
        	    		ProjectionEntry = (AMatrix2D[p][q][r]);
	        		ErrorSino[sino_view][i_r][i_t] -= (ProjectionEntry*(ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] - V));
	   			if (fabs(ErrorSino[sino_view][i_r][i_t]*sqrt(TomoInputsPtr->Weight[sino_view][i_r][i_t])) < TomoInputsPtr->ErrorSinoThresh)
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = true;
				else
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = false;
	   		}
		}
	}

	free(t_ax_start);
	free(t_ax_num);
	free(r_ax_start);
	free(r_ax_num);

	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++)
		if (r_ax_num[p] != 0 && t_ax_num[p] != 0) 
			multifree(AMatrix2D[p], 2);
	free(AMatrix2D);

}



