

#include "XT_Constants.h"
#include <stdio.h>
#include "XT_Structures.h"
#include "XT_ICD_update.h"
#include "XT_AMatrix.h"
#include <math.h>
#include "allocate.h"


/*Updates all voxels from time slices 'time_begin - time_end', x-y slices split along the z-axis from 'slice_begin-slice_end'.
x_rand_select and y_rand_select lists the voxels for each time slice and z-axis slice which needs to be updated.
ErrorSino - Error sinogram
DetectorResponse_XY - detector response in the x-y plane
VoxelLineResponse - Gives the detector response along z-axis for each voxel along z
Iter - The iteration number
zero_count - the number of zero attenuation coefficients
MagUpdateMap - contains the magnitude of each voxel update in the previous iteration which updated that voxel
Mask - All voxels contained in the 'Mask' (contains true or false values for each voxel) are updated*/
Real_t updateVoxels_AttTomo (int32_t time_begin, int32_t time_end, int32_t slice_begin, int32_t slice_end, int32_t xy_begin, int32_t xy_end, int32_t* x_rand_select, int32_t* y_rand_select, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, Real_t** DetectorResponse_XY, AMatrixCol* VoxelLineResponse, int32_t Iter, long int *zero_count, Real_t** MagUpdateMap, uint8_t*** Mask)
{
  int32_t p,q,r,slice,i_new,j_new,k_new,idxr,idxq,idxp,index_xy;
  Real_t V;
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
  Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM];
  Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1];
  bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM];
  bool Time_BDFlag[NHOOD_TIME_MAXDIM-1];
  int32_t maxview = find_max(ScannedObjectPtr->ProjNum, ScannedObjectPtr->N_time);
  /*printf ("maxview = %d, size of AMatrixCol = %d\n",maxview,sizeof(AMatrixCol));*/
  AMatrixCol* AMatrixPtr = (AMatrixCol*)get_spc(maxview, sizeof(AMatrixCol));
  uint8_t AvgNumXElements = (uint8_t)ceil(3*ScannedObjectPtr->delta_xy/SinogramPtr->delta_r);
  
  for (p = 0; p < maxview; p++)
  {
    AMatrixPtr[p].values = (Real_t*)get_spc(AvgNumXElements,sizeof(Real_t));
    AMatrixPtr[p].index = (int32_t*)get_spc(AvgNumXElements,sizeof(int32_t));
  }
  #ifdef DEBUG_HIGH
  fprintf (TomoInputsPtr->debug_file_ptr, "time_begin = %d, time_end = %d, slice_begin = %d, slice_end = %d\n", time_begin, time_end, slice_begin, slice_end);
  #endif
  for (i_new = time_begin; i_new <= time_end; i_new++)
  {
    for (index_xy = xy_begin; index_xy <= xy_end; index_xy++)
    {
      /* printf ("Entering index\n");*/
      k_new = x_rand_select[index_xy];
      j_new = y_rand_select[index_xy];
      MagUpdateMap[j_new][k_new] = 0;
      /* printf ("Entering mask\n");*/
      int32_t sum=0;
      for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++)
      {
        sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
        calcAMatrixColumnforAngle(SinogramPtr, ScannedObjectPtr, DetectorResponse_XY, &(AMatrixPtr[p]), j_new, k_new, sino_view);
        sum+=AMatrixPtr[p].count;
      }
      
      int32_t z_overlap_num = SinogramPtr->z_overlap_num;
      int32_t i_tBeginning=slice_begin*z_overlap_num;
      int32_t sino_viewBegin=ScannedObjectPtr->ProjIdxPtr[i_new][0];
      int32_t NtNrMul=SinogramPtr->N_t*SinogramPtr->N_r;
      int32_t distance= SinogramPtr->N_t;
      Real_t deltar=SinogramPtr->delta_r;
      float errorSinoThresh=(float)TomoInputsPtr->ErrorSinoThresh;
      float errorSinoDelta=(float)TomoInputsPtr->ErrorSinoDelta;
      Real_t projectionValueArray[sum];
      Real_t* projectionValueArrayPointer=&projectionValueArray[0];
      Real_t weightValueArray[sum*z_overlap_num*(slice_end-slice_begin+1)];
      Real_t* weightValueArrayPointer=&weightValueArray[0];
      Real_t* TomoInputsWeightArrayBegin=&TomoInputsPtr->Weight[sino_viewBegin][0][i_tBeginning];
      
      
      bool selectValueArray[sum*z_overlap_num*(slice_end-slice_begin+1)];
      bool* selectValueArrayPointer=&selectValueArray[0];
      bool* ProjSelectArrayBegin=&SinogramPtr->ProjSelect[sino_viewBegin][0][i_tBeginning];
      
      Real_t errorValueArray[sum*z_overlap_num*(slice_end-slice_begin+1)];
      Real_t* errorValueArrayPointer=&errorValueArray[0];
      Real_t* errorSinoValueArrayBegin=&ErrorSino[sino_viewBegin][0][i_tBeginning];             
                      
      for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
        int32_t i_rBeginning=(AMatrixPtr[p].index[0]);
        Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;

        Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;        
        bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;        
        memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
        projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;
        for (q = 0; q < AMatrixPtr[p].count; q++)
        {
          memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num*(slice_end-slice_begin+1));
          weightValueArrayPointer=weightValueArrayPointer+z_overlap_num*(slice_end-slice_begin+1);
          TomoInputsWeightArray=TomoInputsWeightArray+distance;

	  memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num*(slice_end-slice_begin+1));
	  errorValueArrayPointer=errorValueArrayPointer+z_overlap_num*(slice_end-slice_begin+1);
          errorSinoValueArray=errorSinoValueArray+distance;
                    
          memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num*(slice_end-slice_begin+1));	      
          selectValueArrayPointer=selectValueArrayPointer+z_overlap_num*(slice_end-slice_begin+1);
          ProjSelectArray=ProjSelectArray+distance;          
        }
      }
      
      for(p=0;p<sum;p++)
        projectionValueArray[p]=projectionValueArray[p]*deltar;
      projectionValueArrayPointer=&projectionValueArray[0];
      weightValueArrayPointer=&weightValueArray[0];

      errorValueArrayPointer=&errorValueArray[0];      
      selectValueArrayPointer=&selectValueArray[0];      
      
      Real_t UpdatedVoxelValueArray[slice_end-slice_begin+1];
      Real_t ObjectArray[slice_end-slice_begin+1];
      Real_t THETA1[slice_end-slice_begin+1];
      Real_t THETA2[slice_end-slice_begin+1];
      int zeroSkipping[slice_end-slice_begin+1];
      for(p=0;p<(slice_end-slice_begin+1);p++){
        ObjectArray[p]=ScannedObjectPtr->Object[i_new][p+slice_begin+1][j_new][k_new];
      }
      
      for (q=0;q<=(slice_end-slice_begin);q++){
        weightValueArrayPointer=&weightValueArray[0]+q*z_overlap_num;
        errorValueArrayPointer=&errorValueArray[0]+q*z_overlap_num;
        selectValueArrayPointer=&selectValueArray[0]+q*z_overlap_num;
        projectionValueArrayPointer=&projectionValueArray[0];
	Real_t TH2=0.0;
	Real_t TH1=0.0;
        Real_t THETASelTemp;             
        for (p = 0; p < sum; p++)
        {
	  Real_t projectionEntry=*projectionValueArrayPointer;          	  
	  for (r = 0; r < z_overlap_num; r++)
	  {
	     if (*selectValueArrayPointer )
	     {

	         TH2 += (projectionEntry*projectionEntry*(*weightValueArrayPointer));
                 TH1 += -((*errorValueArrayPointer)*projectionEntry*(*weightValueArrayPointer));
             }
	     else
	     {
	         THETASelTemp = errorSinoThresh*errorSinoDelta*sqrt(*weightValueArrayPointer)/fabs(*errorValueArrayPointer);
	         TH2 += (projectionEntry*projectionEntry*THETASelTemp);
                 TH1 += -((*errorValueArrayPointer)*projectionEntry*THETASelTemp);
	     }
	     weightValueArrayPointer++;
	     errorValueArrayPointer++;
	     selectValueArrayPointer++;
          }
        weightValueArrayPointer+=(slice_end-slice_begin)*z_overlap_num;
        errorValueArrayPointer+=(slice_end-slice_begin)*z_overlap_num;
        selectValueArrayPointer+=(slice_end-slice_begin)*z_overlap_num;
        projectionValueArrayPointer++;	                        	          
        }
        THETA2[q]=TH2;
        THETA1[q]=TH1;        
      }      
      
      
      projectionValueArrayPointer=&projectionValueArray[0];
      weightValueArrayPointer=&weightValueArray[0];

      errorValueArrayPointer=&errorValueArray[0];      
      selectValueArrayPointer=&selectValueArray[0]; 
      
            
      for (slice = slice_begin; slice <= slice_end; slice++) {
        /*	printf ("Entering slice\n");*/
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
                  for (r = 0; r <NHOOD_X_MAXDIM; r++){
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
          /*V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new];*/ /*Store the present value of the voxel*/
          #ifdef ZERO_SKIPPING
          /*Zero Skipping Algorithm*/
          ZSFlag = true;
          if(ObjectArray[slice-slice_begin] == 0.0 && Iter > 1) /*Iteration starts from 1. Iteration 0 corresponds to initial cost before ICD*/
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
            zeroSkipping[slice-slice_begin]=1;
        	UpdatedVoxelValueArray[slice-slice_begin] = CE_FunctionalSubstitution(ObjectArray[slice-slice_begin], THETA1[slice-slice_begin], THETA2[slice-slice_begin], ScannedObjectPtr, TomoInputsPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag);
          }
          else
          (*zero_count)++;
        }
      }
      projectionValueArrayPointer=&projectionValueArray[0];
      weightValueArrayPointer=&weightValueArray[0];
      errorValueArrayPointer=&errorValueArray[0];            
      for(p=0;p<(slice_end-slice_begin+1);p++){
        if(zeroSkipping[p]==1){
          ScannedObjectPtr->Object[i_new][slice_begin+p+1][j_new][k_new]=UpdatedVoxelValueArray[p];
          MagUpdateMap[j_new][k_new] += fabs(UpdatedVoxelValueArray[p] - ObjectArray[p]);
          total_vox_mag += fabs(UpdatedVoxelValueArray[p]);
        }
      }
      sino_viewBegin = ScannedObjectPtr->ProjIdxPtr[i_new][0];
      int32_t i_tBegin = slice_begin*z_overlap_num;
      int32_t s=0;            	            
      for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
        int32_t i_rBegin = (AMatrixPtr[p].index[0]);
        Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBegin)*SinogramPtr->N_t;        
        bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBegin)*SinogramPtr->N_t;        
        	  
	for (q = 0; q < AMatrixPtr[p].count; q++)
        {
            Real_t ProjectionEntry = (*projectionValueArrayPointer);
	    for (s=0;s<=(slice_end-slice_begin);s++){
	      if(zeroSkipping[s]==1){ 	            	   	    
	        Real_t temp=(ProjectionEntry*(UpdatedVoxelValueArray[s] - ObjectArray[s]));
	        for (r = 0; r <z_overlap_num; r++)
	        {
	            *errorSinoValueArray =(*errorValueArrayPointer)- temp;
		    *ProjSelectArray = (fabs((*errorValueArrayPointer- temp)*sqrt(*weightValueArrayPointer)) < errorSinoThresh);
		    weightValueArrayPointer++;
		    errorValueArrayPointer++;
		    errorSinoValueArray++;
		    ProjSelectArray++;	
	        }
	      }
	      else{
	       weightValueArrayPointer+=z_overlap_num;
	       errorValueArrayPointer+=z_overlap_num;
	       errorSinoValueArray+=z_overlap_num;
	       ProjSelectArray+=z_overlap_num;
	      }  
	    }
	    projectionValueArrayPointer++;
	    errorSinoValueArray=errorSinoValueArray-z_overlap_num*(slice_end-slice_begin+1)+distance;
	    ProjSelectArray=ProjSelectArray-z_overlap_num*(slice_end-slice_begin+1)+distance;
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



Real_t compute_voxel_update_AMat2D (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new)
{
  	int32_t p, q, r, sino_view;
	Real_t V,THETA1,THETA2,THETASelTemp,***AMatrix2D, *AMatrix2DLine;
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
	/*	printf("Start 2D AMatrix computation, p = %d, i_new = %d, j_new = %d, k_new = %d\n", p, i_new, j_new, k_new);*/
		t_ax_start[p] = slice*SinogramPtr->z_overlap_num;
		t_ax_num[p] = SinogramPtr->z_overlap_num;

		compute_2DAMatrixLine(SinogramPtr, &(AMatrix2DLine), &(AMatrixPtr[p]), &(r_ax_start[p]), &(r_ax_num[p]));
		compute_LapMatrix_4m_AMatrix(SinogramPtr, &(AMatrix2D[p]), &(AMatrix2DLine), &(r_ax_start[p]), &(r_ax_num[p]), &(t_ax_start[p]), &(t_ax_num[p]));
	/*	printf("End 2D AMatrix computation\n");*/
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
	
/*	printf("updating e vector\n");*/
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

	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++)
		if (r_ax_num[p] != 0 && t_ax_num[p] != 0) 
			multifree(AMatrix2D[p], 2);
	free(AMatrix2D);
	free(t_ax_start);
	free(t_ax_num);
	free(r_ax_start);
	free(r_ax_num);
	
        return 0;

}


Real_t updateVoxels_PhCon_Tomo (int32_t time_begin, int32_t time_end, int32_t slice_begin, int32_t slice_end, int32_t xy_begin, int32_t xy_end, int32_t* x_rand_select, int32_t* y_rand_select, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, Real_t** DetectorResponse_XY, AMatrixCol* VoxelLineResponse, int32_t Iter, long int *zero_count, Real_t** MagUpdateMap, uint8_t*** Mask)
{
  int32_t p,q,r,slice,i_new,j_new,k_new,idxr,idxq,idxp,index_xy;
  Real_t V;
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

    Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM]; 
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
		compute_voxel_update_AMat2D (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, AMatrixPtr, VoxelLineResponse, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag, i_new, slice, j_new, k_new);
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

Real_t updateVoxels (int32_t time_begin, int32_t time_end, int32_t slice_begin, int32_t slice_end, int32_t xy_begin, int32_t xy_end, int32_t* x_rand_select, int32_t* y_rand_select, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, Real_t** DetectorResponse_XY, AMatrixCol* VoxelLineResponse, int32_t Iter, long int *zero_count, Real_t** MagUpdateMap, uint8_t*** Mask)
{
	Real_t total_vox_mag;
#ifdef PHASE_CONTRAST_TOMOGRAPHY
	total_vox_mag = updateVoxels_PhCon_Tomo (time_begin, time_end, slice_begin, slice_end, xy_begin, xy_end, x_rand_select, y_rand_select, SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse_XY, VoxelLineResponse, Iter, zero_count, MagUpdateMap, Mask);
#else
	total_vox_mag = updateVoxels_AttTomo (time_begin, time_end, slice_begin, slice_end, xy_begin, xy_end, x_rand_select, y_rand_select, SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse_XY, VoxelLineResponse, Iter, zero_count, MagUpdateMap, Mask);
#endif
	return (total_vox_mag);
}
