#include <stdio.h>
#include "XT_Structures.h"
#include "XT_AMatrix.h"
#include "allocate.h"



void forward_project_voxel_AMat1D (Sinogram* SinogramPtr, Real_t voxel_val, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, int32_t sino_idx, int32_t slice)
{
	int32_t m, idxBeginning, n, z_overlap_num;
	Real_t val;
	z_overlap_num = SinogramPtr->z_overlap_num;
        int32_t i_tBeginning=slice*z_overlap_num;
        idxBeginning=AMatrixPtr->index[0];
        Real_t* errorSinoValueArray=&ErrorSino[sino_idx][idxBeginning][i_tBeginning];       
	for (m = 0; m < AMatrixPtr->count; m++)
	{
		/*idx = AMatrixPtr->index[m];*/
		val = AMatrixPtr->values[m]*SinogramPtr->delta_r;
		/*val = AMatrixPtr->values[m];
		for (n = 0; n < VoxelLineResponse->count; n++){*/
		for (n = 0; n < z_overlap_num; n++){
			/*ErrorSino[sino_idx][idx][VoxelLineResponse->index[n]] += voxel_val*val*VoxelLineResponse->values[n];*/
			/*ErrorSino[sino_idx][idx][slice*z_overlap_num + n] += voxel_val*val;*/
                        *errorSinoValueArray += voxel_val*val;
                        errorSinoValueArray++;
		}
                errorSinoValueArray=errorSinoValueArray-z_overlap_num+SinogramPtr->N_t;         
	}
}

	
void forward_project_voxel_AMat2D (Sinogram* SinogramPtr, Real_t voxel_val, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, int32_t sino_idx)
{
	int32_t m, n;
	int32_t r_ax_start, r_ax_num, t_ax_start, t_ax_num;
	Real_t **AMatrix2D, **AMatrixTemp;	

	compute_2DAMatrix_4m_1D(&(AMatrixTemp), AMatrixPtr, VoxelLineResponse, &r_ax_start, &r_ax_num, &t_ax_start, &t_ax_num);
	compute_LapMatrix_4m_AMatrix(SinogramPtr, &(AMatrix2D), &(AMatrixTemp), &(r_ax_start), &(r_ax_num), &(t_ax_start), &(t_ax_num));
 
	for (m = 0; m < r_ax_num; m++)
	{
		for (n = 0; n < t_ax_num; n++){
			ErrorSino[sino_idx][r_ax_start+m][t_ax_start+n] += voxel_val*AMatrix2D[m][n];
		}
	}
	
	if (r_ax_num != 0 && t_ax_num != 0)
		multifree(AMatrix2D,2);
}

	
void forward_project_voxel (Sinogram* SinogramPtr, Real_t voxel_val, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, int32_t sino_idx, int32_t slice)
{
#ifdef PHASE_CONTRAST_TOMOGRAPHY
	forward_project_voxel_AMat2D (SinogramPtr, voxel_val, ErrorSino, AMatrixPtr, VoxelLineResponse, sino_idx);
#else
	forward_project_voxel_AMat1D (SinogramPtr, voxel_val, ErrorSino, AMatrixPtr, VoxelLineResponse, sino_idx, slice);
#endif
}

