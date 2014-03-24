
#include <stdio.h>
#include "XT_Structures.h"
#include "XT_ICD_update.h"
#include "XT_AMatrix.h"
#include <math.h>
#include "allocate.h"



void compute_voxel_update_AMat1D (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new)
{
  	int32_t p, q, r, sino_view, z_overlap_num;
	Real_t V,THETA1,THETA2,THETASelTemp;
	Real_t UpdatedVoxelValue, ProjectionEntry;
  	/*int32_t i_r;*/
        V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]; /*Store the present value of the voxel*/
	z_overlap_num = SinogramPtr->z_overlap_num;
        int32_t i_tBeginning=slice*z_overlap_num;
        float errorSinoThresh=(float)TomoInputsPtr->ErrorSinoThresh;
        float errorSinoDelta=(float)TomoInputsPtr->ErrorSinoDelta;
        int32_t* scannedObjectProjIdxNumber=ScannedObjectPtr->ProjIdxPtr[i_new];

	THETA1 = 0.0;
	THETA2 = 0.0;
        /* 
	fprintf(TomoInputsPtr->debug_file_ptr, "thread num %d sinoview ",omp_get_thread_num());
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
          fprintf(TomoInputsPtr->debug_file_ptr," %d ",ScannedObjectPtr->ProjIdxPtr[i_new][p]);
        }
        fprintf(TomoInputsPtr->debug_file_ptr, "\n i_rBeginning is ");
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){  
          fprintf(TomoInputsPtr->debug_file_ptr," %d ",AMatrixPtr[p].index[0]);       
        }
        fprintf(TomoInputsPtr->debug_file_ptr,"\n ");
        */
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
	sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
        int32_t i_rBeginning=(AMatrixPtr[p].index[0]);

        Real_t* TomoInputsWeightArray=&TomoInputsPtr->Weight[sino_view][i_rBeginning][i_tBeginning];
        Real_t* errorSinoValueArray=&ErrorSino[sino_view][i_rBeginning][i_tBeginning];
        bool* ProjSelectArray=&SinogramPtr->ProjSelect[sino_view][i_rBeginning][i_tBeginning];

	for (q = 0; q < AMatrixPtr[p].count; q++)
	{
       	    	ProjectionEntry = (AMatrixPtr[p].values[q]*SinogramPtr->delta_r);              


		for (r = 0; r < z_overlap_num; r++)
		{ 
			if (*ProjSelectArray == true)
			{

	           		THETA2 += (ProjectionEntry*ProjectionEntry*(*TomoInputsWeightArray));
               			THETA1 += -((*errorSinoValueArray)*ProjectionEntry*(*TomoInputsWeightArray));
            		}
			else
			{
				THETASelTemp = errorSinoThresh*errorSinoDelta*sqrt(*TomoInputsWeightArray)/fabs(*errorSinoValueArray);
	            		THETA2 += (ProjectionEntry*ProjectionEntry*THETASelTemp);
            			THETA1 += -((*errorSinoValueArray)*ProjectionEntry*THETASelTemp);
			}
                        TomoInputsWeightArray++;
                        errorSinoValueArray++;
                        ProjSelectArray++;
            	}
                TomoInputsWeightArray=TomoInputsWeightArray-z_overlap_num+SinogramPtr->N_t;
                errorSinoValueArray=errorSinoValueArray-z_overlap_num+SinogramPtr->N_t;         
                ProjSelectArray=ProjSelectArray-z_overlap_num+SinogramPtr->N_t;    
	}
        }


            /*Solve the 1-D optimization problem
            TODO : What if theta1 = 0 ? Then this will give error*/


        UpdatedVoxelValue = CE_FunctionalSubstitution(V, THETA1, THETA2, ScannedObjectPtr, TomoInputsPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag);
              
        ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] = UpdatedVoxelValue;
	
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
		sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
                int32_t i_rBeginning=(AMatrixPtr[p].index[0]);
                Real_t* TomoInputsWeightArray=&TomoInputsPtr->Weight[sino_view][i_rBeginning][i_tBeginning];
                Real_t* errorSinoValueArray=&ErrorSino[sino_view][i_rBeginning][i_tBeginning];
                bool* ProjSelectArray=&SinogramPtr->ProjSelect[sino_view][i_rBeginning][i_tBeginning];
		for (q = 0; q < AMatrixPtr[p].count; q++)
        	{
               	    	/*i_r = (AMatrixPtr[p].index[q]);*/
        	    	ProjectionEntry = (AMatrixPtr[p].values[q]*SinogramPtr->delta_r);
                        /*Real_t* TomoInputsWeightArray=&TomoInputsPtr->Weight[sino_view][i_r][i_tBeginning];
                        Real_t* errorSinoValueArray=&ErrorSino[sino_view][i_r][i_tBeginning];
                        bool* ProjSelectArray=&SinogramPtr->ProjSelect[sino_view][i_r][i_tBeginning];*/


			for (r = 0; r < z_overlap_num; r++)
			{ 
	        		*errorSinoValueArray -= (ProjectionEntry*(UpdatedVoxelValue - V));
			        *ProjSelectArray = (fabs((*errorSinoValueArray)*sqrt(*TomoInputsWeightArray)) < errorSinoThresh);
                                TomoInputsWeightArray++;
                                errorSinoValueArray++;
                                ProjSelectArray++;
	   		}
                TomoInputsWeightArray=TomoInputsWeightArray-z_overlap_num+SinogramPtr->N_t;
                errorSinoValueArray=errorSinoValueArray-z_overlap_num+SinogramPtr->N_t;         
                ProjSelectArray=ProjSelectArray-z_overlap_num+SinogramPtr->N_t; 
		}
	}
}



void compute_voxel_update_AMat2D (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponseSlice, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new)
{
  	int32_t p, q, r, sino_view;
	Real_t V,THETA1,THETA2,THETASelTemp,***AMatrix2D, **AMatrixTemp;
	Real_t UpdatedVoxelValue, ProjectionEntry;
        float errorSinoThresh=TomoInputsPtr->ErrorSinoThresh;
        float errorSinoDelta=TomoInputsPtr->ErrorSinoDelta;
  	int32_t i_r;
        int32_t i_tBeginning;
        int32_t* scannedObjectProjIdxNumber=ScannedObjectPtr->ProjIdxPtr[i_new];
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
	/*	printf("Start 2D AMatrix computation, p = %d, i_new = %d, j_new = %d, k_new = %d\n", p, i_new, j_new, k_new);*/
		compute_2DAMatrix_4m_1D(&(AMatrixTemp), &(AMatrixPtr[p]), (VoxelLineResponseSlice), &(r_ax_start[p]), &(r_ax_num[p]), &(t_ax_start[p]), &(t_ax_num[p]));
		compute_LapMatrix_4m_AMatrix(SinogramPtr, &(AMatrix2D[p]), &(AMatrixTemp), &(r_ax_start[p]), &(r_ax_num[p]), &(t_ax_start[p]), &(t_ax_num[p]));
                i_tBeginning= t_ax_start[p];

	/*	printf("End 2D AMatrix computation\n");*/
		for (q = 0; q < r_ax_num[p]; q++)
		{
      	    		i_r = r_ax_start[p] + q;
                        Real_t* TomoInputsWeightArray=&TomoInputsPtr->Weight[sino_view][i_r][i_tBeginning];
                        Real_t* errorSinoValueArray=&ErrorSino[sino_view][i_r][i_tBeginning];
                        bool* ProjSelectArray=&SinogramPtr->ProjSelect[sino_view][i_r][i_tBeginning];

			for (r = 0; r < t_ax_num[p]; r++)
			{ 
       	    			ProjectionEntry = (AMatrix2D[p][q][r]);
				if (*ProjSelectArray == true)
				{
	           			THETA2 += (ProjectionEntry*ProjectionEntry*(*TomoInputsWeightArray));
               				THETA1 += -((*errorSinoValueArray)*ProjectionEntry*(*TomoInputsWeightArray));
            			}
				else
				{
					THETASelTemp = errorSinoThresh*errorSinoDelta*sqrt(*TomoInputsWeightArray)/fabs(*errorSinoValueArray);
	            			THETA2 += (ProjectionEntry*ProjectionEntry*THETASelTemp);
            				THETA1 += -((*errorSinoValueArray)*ProjectionEntry*THETASelTemp);
				}
                                TomoInputsWeightArray++;
                                errorSinoValueArray++;
                                ProjSelectArray++;
            		}
		}
        }


            /*Solve the 1-D optimization problem
            TODO : What if theta1 = 0 ? Then this will give error*/


        UpdatedVoxelValue = CE_FunctionalSubstitution(V, THETA1, THETA2, ScannedObjectPtr, TomoInputsPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag);
              
        ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] = UpdatedVoxelValue;
	
/*	printf("updating e vector\n");*/
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
		sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
                i_tBeginning= t_ax_start[p];
		for (q = 0; q < r_ax_num[p]; q++)
        	{
               	    	i_r = r_ax_start[p] + q;
                        Real_t* TomoInputsWeightArray=&TomoInputsPtr->Weight[sino_view][i_r][i_tBeginning];
                        Real_t* errorSinoValueArray=&ErrorSino[sino_view][i_r][i_tBeginning];
                        bool* ProjSelectArray=&SinogramPtr->ProjSelect[sino_view][i_r][i_tBeginning];
			for (r = 0; r < t_ax_num[p]; r++)
			{ 
        	    		ProjectionEntry = (AMatrix2D[p][q][r]);
	        		*errorSinoValueArray -= (ProjectionEntry*(UpdatedVoxelValue - V));
			        *ProjSelectArray = (fabs((*errorSinoValueArray)*sqrt(*TomoInputsWeightArray)) < errorSinoThresh);
                                TomoInputsWeightArray++;
                                errorSinoValueArray++;
                                ProjSelectArray++;
	   		}
		}
	}

	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++)
		if (r_ax_num[p] != 0 && t_ax_num[p] != 0) 
			multifree(AMatrix2D[p], 2);
	free(AMatrix2D);
	free(t_ax_start);
	free(t_ax_num);
	free(r_ax_start);
	free(r_ax_num);


}



void compute_voxel_update (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new)
{
#ifdef PHASE_CONTRAST_TOMOGRAPHY
	compute_voxel_update_AMat2D (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, AMatrixPtr, &(VoxelLineResponse[slice]), Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag, i_new, slice, j_new, k_new);
#else
	compute_voxel_update_AMat1D (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, AMatrixPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag, i_new, slice, j_new, k_new);
#endif
}





