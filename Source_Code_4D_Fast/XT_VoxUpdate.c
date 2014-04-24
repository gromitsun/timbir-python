

#include "XT_Constants.h"
#include <stdio.h>
#include "XT_Structures.h"
#include "XT_ICD_update.h"
#include "XT_AMatrix.h"
#include <math.h>
#include "allocate.h"



Real_t compute_voxel_update_AMat1D (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new,Real_t* projectionValueArrayPointer,Real_t* weightValueArrayPointer,int32_t slice_begin,int32_t slice_end,bool* selectValueArrayPointer,Real_t* errorValueArrayPointer,Real_t V,Real_t THETA1,Real_t THETA2)
{
	Real_t UpdatedVoxelValue;
        /*V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new];*/ 

        /*
	for (q = 0; q < sum; q++)
	{
	        Real_t projectionEntry=*projectionValueArrayPointer;
		for (r = 0; r < z_overlap_num; r++)
		{
			if (*selectValueArrayPointer )
			{

	           		THETA2 += (projectionEntry*projectionEntry*(*weightValueArrayPointer));
               			THETA1 += -((*errorValueArrayPointer)*projectionEntry*(*weightValueArrayPointer));
            		}
			else
			{
				THETASelTemp = errorSinoThresh*errorSinoDelta*sqrt(*weightValueArrayPointer)/fabs(*errorValueArrayPointer);
	            		THETA2 += (projectionEntry*projectionEntry*THETASelTemp);
            			THETA1 += -((*errorValueArrayPointer)*projectionEntry*THETASelTemp);
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
         */     
        
        UpdatedVoxelValue = CE_FunctionalSubstitution(V, THETA1, THETA2, ScannedObjectPtr, TomoInputsPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag);
              

	
/*
        weightValueArrayPointer=copyweightValueArrayPointer;
        errorValueArrayPointer=&errorValueArray[0];
        selectValueArrayPointer=&selectValueArray[0];
        projectionValueArrayPointer=copyprojectionValueArrayPointer;	 
	 
	  for(q=0;q<ScannedObjectPtr->ProjNum[i_new];q++){
	    int32_t i_rBegin=i_rArray[2*q];
	    int32_t i_rEnd=i_rArray[2*q+1];
	    for(t=i_rBegin;t<=i_rEnd;t++){
	      Real_t temp=(*projectionValueArrayPointer)*(UpdatedVoxelValue - V);	    
	      for(r=0;r<z_overlap_num;r++){
	        ErrorSino[sino_viewBegin+q][t][i_tBeginning + r]=(*errorValueArrayPointer)-temp;
	        SinogramPtr->ProjSelect[sino_viewBegin+q][t][i_tBeginning +r]=(fabs((*errorValueArrayPointer-temp)*sqrt(*weightValueArrayPointer)) < errorSinoThresh);
                weightValueArrayPointer++;	      	        
	        errorValueArrayPointer++;
	        selectValueArrayPointer++;
	      }
	    weightValueArrayPointer+=(slice_end-slice_begin)*z_overlap_num;  
	    projectionValueArrayPointer++; 	      
	    }
	    
	  }	 	
*/	
	return UpdatedVoxelValue;	





/*
  	int32_t p, q, r, t, sino_viewBegin, z_overlap_num;
	Real_t V,THETA1,THETA2,THETASelTemp;	
	Real_t UpdatedVoxelValue;
        V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]; 
	z_overlap_num = SinogramPtr->z_overlap_num;
        int32_t i_tBeginning=slice*z_overlap_num;
        float errorSinoThresh=(float)TomoInputsPtr->ErrorSinoThresh;
        float errorSinoDelta=(float)TomoInputsPtr->ErrorSinoDelta;
        Real_t deltar=SinogramPtr->delta_r;

	THETA1 = 0.0;
	THETA2 = 0.0;
        sino_viewBegin=ScannedObjectPtr->ProjIdxPtr[i_new][0];
        int32_t sino_viewEnd=ScannedObjectPtr->ProjIdxPtr[i_new][ScannedObjectPtr->ProjNum[i_new]-1];
        
        Real_t* TomoInputsWeightArrayBegin=&TomoInputsPtr->Weight[sino_viewBegin][0][i_tBeginning];
        Real_t* errorSinoValueArrayBegin=&ErrorSino[sino_viewBegin][0][i_tBeginning];
        bool* ProjSelectArrayBegin=&SinogramPtr->ProjSelect[sino_viewBegin][0][i_tBeginning]; 
        int32_t NtNrMul=SinogramPtr->N_t*SinogramPtr->N_r; 
        int32_t distance=  SinogramPtr->N_t;

        int32_t i_rArray[ScannedObjectPtr->ProjNum[i_new]*2];
        int32_t sum=0;
        for(p=0;p<ScannedObjectPtr->ProjNum[i_new];p++){
          i_rArray[p*2]=AMatrixPtr[p].index[0];
          i_rArray[p*2+1]=AMatrixPtr[p].index[0]+AMatrixPtr[p].count-1;
          sum+=AMatrixPtr[p].count;
        }
        int32_t maxI_r=0;
        int32_t pMax=0;
        int32_t pMiddle=0;
        int32_t startI_r=AMatrixPtr[0].index[0];
        int32_t minI_r=startI_r;
        int32_t pMin=0;        
        for(p=0;p<ScannedObjectPtr->ProjNum[i_new];p++){
          if(i_rArray[p*2]>maxI_r){
            maxI_r=i_rArray[p*2];
            pMax=p;
          }  
          if(i_rArray[p*2]<minI_r){
            minI_r=i_rArray[p*2];
            pMin=p;  
          }  
        }
       int32_t pMiddleCount=0;
       p=pMiddle;        
        if(pMax>0&&pMax<(ScannedObjectPtr->ProjNum[i_new]-1)){
         pMiddle=pMax;
         while(AMatrixPtr[p].index[0]==maxI_r){
           pMiddleCount++;
           p++;
         }
       }  
       else{
         pMiddle=pMin;
         while(AMatrixPtr[p].index[0]==minI_r){
           pMiddleCount++;
           p++;
         }
       }  
         

         
                

        fprintf(TomoInputsPtr->debug_file_ptr, "pMiddleCount is %d i_r-2 %d i_r-1 %d i_r %d i_r+pMiddleCount %d i_r+pMiddleCount+1 %d \n",pMiddleCount,AMatrixPtr[pMiddle-2].index[0],AMatrixPtr[pMiddle-1].index[0],AMatrixPtr[pMiddle].index[0],AMatrixPtr[pMiddle+pMiddleCount].index[0],AMatrixPtr[pMiddle+pMiddleCount+1].index[0]);              
        
        Real_t weightValueArray[sum*z_overlap_num];
        Real_t errorValueArray[sum*z_overlap_num];
        bool selectValueArray[sum*z_overlap_num];
        Real_t projectionValueArray[sum];
        Real_t* weightValueArrayPointer=&weightValueArray[0];
        Real_t* errorValueArrayPointer=&errorValueArray[0];
        bool* selectValueArrayPointer=&selectValueArray[0];
        Real_t* projectionValueArrayPointer=&projectionValueArray[0];

        if(pMiddle<(ScannedObjectPtr->ProjNum[i_new]/2)){
	  for (p = 0; p <= (pMiddle-1); p++){
            int32_t i_rBeginning=(AMatrixPtr[p].index[0]);          	                   
            Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          }
          int32_t i_rBeginning=(AMatrixPtr[pMiddle].index[0]);
          Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
          Real_t* errorSinoValueArray=errorSinoValueArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
          bool* ProjSelectArray=ProjSelectArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;          	                             
	  for (p = pMiddle; p < (pMiddle+pMiddleCount); p++){
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          } 
	  for (p = (pMiddle+1); p <= 2*pMiddle; p++){
            int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
            Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          } 
	  for (p = (2*pMiddle+1); p <=( ScannedObjectPtr->ProjNum[i_new]-1); p++){
            int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
            Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          }                                     
        }
        else{
	  for (p = 0; p <= (ScannedObjectPtr->ProjNum[i_new]-1-(2*(ScannedObjectPtr->ProjNum[i_new]-1-pMiddle)+1)); p++){
            int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
            Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          }
	  for (p = (ScannedObjectPtr->ProjNum[i_new]-(2*(ScannedObjectPtr->ProjNum[i_new]-1-pMiddle)+1)); p <=(pMiddle-1); p++){
            int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
            Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          }
	  for (p = pMiddle; p <(pMiddle+1); p++){
            int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
            Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          }
	  for (p = (pMiddle+1); p <=(ScannedObjectPtr->ProjNum[i_new]-1); p++){
            int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
            Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          }                                      
        }
*/

        
        /*
        if(pMiddle<(ScannedObjectPtr->ProjNum[i_new]/2)){
	  for (p = 0; p <= (pMiddle-1); p++){
                int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
            	int32_t size=AMatrixPtr[p].count;	                   
            	Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	Real_t* TomoInputsWeightArray2=TomoInputsWeightArray+2*(pMiddle-p)*NtNrMul;
            	Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	Real_t* errorSinoValueArray2=errorSinoValueArray+2*(pMiddle-p)*NtNrMul;
            	bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	bool* ProjSelectArray2=ProjSelectArray+2*(pMiddle-p)*NtNrMul;
 
 	    	memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*size);
 	    	projectionValueArrayPointer=projectionValueArrayPointer+size;
 	    	memcpy(projectionValueArrayPointer,&(AMatrixPtr[2*pMiddle-p].values[0]),sizeof(Real_t)*size);
 	    	projectionValueArrayPointer=projectionValueArrayPointer+size; 	                      
	    	for (q = 0; q < size; q++)
	    	{
	      		memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	      		weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              		TomoInputsWeightArray=TomoInputsWeightArray+distance;
	      		memcpy(weightValueArrayPointer,TomoInputsWeightArray2,sizeof(Real_t)*z_overlap_num);
	      		weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              		TomoInputsWeightArray2=TomoInputsWeightArray2+distance;  
              
                          	      
	      		memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	      		errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              		errorSinoValueArray=errorSinoValueArray+distance;
	      		memcpy(errorValueArrayPointer,errorSinoValueArray2,sizeof(Real_t)*z_overlap_num);
	      		errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              		errorSinoValueArray2=errorSinoValueArray2+distance;               
              
                      	      
	      		memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
              		selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              		ProjSelectArray=ProjSelectArray+distance;
	      		memcpy(selectValueArrayPointer,ProjSelectArray2,sizeof(bool)*z_overlap_num);	      
              		selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              		ProjSelectArray2=ProjSelectArray2+distance;
              	}		  
          }
          int32_t i_rBeginning=(AMatrixPtr[pMiddle].index[0]);	                   
          Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
          Real_t* errorSinoValueArray=errorSinoValueArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
          bool* ProjSelectArray=ProjSelectArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	  memcpy(projectionValueArrayPointer,&(AMatrixPtr[pMiddle].values[0]),sizeof(Real_t)*AMatrixPtr[pMiddle].count);
 	  projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[pMiddle].count;                  
	  for (q = 0; q < AMatrixPtr[pMiddle].count; q++)
	  {
	      memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	      weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	      memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	      errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              errorSinoValueArray=errorSinoValueArray+distance;         	      
	      memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
              selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              ProjSelectArray=ProjSelectArray+distance;
	  }
	  for (p = (2*pMiddle+1); p <= (ScannedObjectPtr->ProjNum[i_new]-1); p++){
            i_rBeginning=(AMatrixPtr[p].index[0]);	                   
            TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	    for (q = 0; q < AMatrixPtr[p].count; q++)
	    {
	        memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	        weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
                TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	        memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	        errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
                errorSinoValueArray=errorSinoValueArray+distance;         	      
	        memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
                selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
                ProjSelectArray=ProjSelectArray+distance;
	    }
          }
	                    
        }
        else{
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
          int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
          Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
          Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
          bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	  memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	  projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	  for (q = 0; q < AMatrixPtr[p].count; q++)
	  {
	      memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	      weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	      memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	      errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              errorSinoValueArray=errorSinoValueArray+distance;         	      
	      memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
              selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              ProjSelectArray=ProjSelectArray+distance;
	  }
        }        
        
        }                
        */
        
        

        /*
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
          int32_t i_rBeginning=(AMatrixPtr[p].index[0]);	                   
          Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
          Real_t* errorSinoValueArray=errorSinoValueArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
          bool* ProjSelectArray=ProjSelectArrayBegin+p*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	  memcpy(projectionValueArrayPointer,&(AMatrixPtr[p].values[0]),sizeof(Real_t)*AMatrixPtr[p].count);
 	  projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[p].count;                  
	  for (q = 0; q < AMatrixPtr[p].count; q++)
	  {
	      memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	      weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	      memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	      errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              errorSinoValueArray=errorSinoValueArray+distance;         	      
	      memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
              selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              ProjSelectArray=ProjSelectArray+distance;
	  }
        }
        */
        /*
            for(p=0;p<sum;p++)
	      projectionValueArray[p]=projectionValueArray[p]*deltar;        

        
       
        weightValueArrayPointer=&weightValueArray[0];
        errorValueArrayPointer=&errorValueArray[0];
        selectValueArrayPointer=&selectValueArray[0];
        projectionValueArrayPointer=&projectionValueArray[0];

	for (q = 0; q < sum; q++)
	{
	        Real_t projectionEntry=*projectionValueArrayPointer;
		for (r = 0; r < z_overlap_num; r++)
		{
			if (*selectValueArrayPointer )
			{

	           		THETA2 += (projectionEntry*projectionEntry*(*weightValueArrayPointer));
               			THETA1 += -((*errorValueArrayPointer)*projectionEntry*(*weightValueArrayPointer));
            		}
			else
			{
				THETASelTemp = errorSinoThresh*errorSinoDelta*sqrt(*weightValueArrayPointer)/fabs(*errorValueArrayPointer);
	            		THETA2 += (projectionEntry*projectionEntry*THETASelTemp);
            			THETA1 += -((*errorValueArrayPointer)*projectionEntry*THETASelTemp);
			}
			weightValueArrayPointer++;
			errorValueArrayPointer++;
			selectValueArrayPointer++;
            	}
            	projectionValueArrayPointer++;	    
	}


        
        UpdatedVoxelValue = CE_FunctionalSubstitution(V, THETA1, THETA2, ScannedObjectPtr, TomoInputsPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag);
              

	

        weightValueArrayPointer=&weightValueArray[0];
        errorValueArrayPointer=&errorValueArray[0];
        selectValueArrayPointer=&selectValueArray[0];
        projectionValueArrayPointer=&projectionValueArray[0];	 
	 
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
		int32_t sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
		for (q = 0; q < AMatrixPtr[p].count; q++)
        	{
               	    	int32_t i_r = (AMatrixPtr[p].index[q]);
        	    	Real_t ProjectionEntry = (AMatrixPtr[p].values[q]*deltar);
			for (r = 0; r < z_overlap_num; r++)
			{ 
				int32_t i_t = slice*z_overlap_num + r;
	        		ErrorSino[sino_view][i_r][i_t] -= (ProjectionEntry*(UpdatedVoxelValue - V));
	   			if (fabs(ErrorSino[sino_view][i_r][i_t]*sqrt(TomoInputsPtr->Weight[sino_view][i_r][i_t])) < TomoInputsPtr->ErrorSinoThresh)
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = true;
				else
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = false;
	   		}
		}
	}	 	
	
	return UpdatedVoxelValue;	
        */

















/*
  	int32_t p, q, r, t, sino_viewBegin, z_overlap_num;
	Real_t V,THETA1,THETA2,THETASelTemp;	
	Real_t UpdatedVoxelValue;
        V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]; 
	z_overlap_num = SinogramPtr->z_overlap_num;
        int32_t i_tBeginning=slice*z_overlap_num;
        float errorSinoThresh=(float)TomoInputsPtr->ErrorSinoThresh;
        float errorSinoDelta=(float)TomoInputsPtr->ErrorSinoDelta;
        Real_t deltar=SinogramPtr->delta_r;

	THETA1 = 0.0;
	THETA2 = 0.0;
        sino_viewBegin=ScannedObjectPtr->ProjIdxPtr[i_new][0];
        
        Real_t* TomoInputsWeightArrayBegin=&TomoInputsPtr->Weight[sino_viewBegin][0][i_tBeginning];
        Real_t* errorSinoValueArrayBegin=&ErrorSino[sino_viewBegin][0][i_tBeginning];
        bool* ProjSelectArrayBegin=&SinogramPtr->ProjSelect[sino_viewBegin][0][i_tBeginning]; 
        int32_t NtNrMul=SinogramPtr->N_t*SinogramPtr->N_r; 
        int32_t distance=  SinogramPtr->N_t;

        int32_t i_rArray[ScannedObjectPtr->ProjNum[i_new]*2];
        int32_t sum=0;
        for(p=0;p<ScannedObjectPtr->ProjNum[i_new];p++){
          i_rArray[p*2]=AMatrixPtr[p].index[0];
          i_rArray[p*2+1]=AMatrixPtr[p].index[0]+AMatrixPtr[p].count-1;
          sum+=AMatrixPtr[p].count;
        }
        
        int32_t maxI_r=0;
        int32_t pMax=0;
        int32_t pMiddle=0;
        int32_t maxCount=0;
        int32_t startI_r=AMatrixPtr[0].index[0];
        int32_t minI_r=startI_r;
        int32_t pMin=0;        
        for(p=0;p<ScannedObjectPtr->ProjNum[i_new];p++){
          if(i_rArray[p*2]>maxI_r){
            maxI_r=i_rArray[p*2];
            pMax=p;
            maxCount=AMatrixPtr[p].count;
          }  
          if(i_rArray[p*2]<minI_r){
            minI_r=i_rArray[p*2];
            pMin=p;  
          }  
        }
        if(pMax>0&&pMax<ScannedObjectPtr->ProjNum[i_new])
         pMiddle=pMax;
       else
         pMiddle=pMin;
         
       int32_t copyCount=0;                 
       if((2*pMiddle)>=ScannedObjectPtr->ProjNum[i_new])
         copyCount=pMiddle+1;
       else
         copyCount= ScannedObjectPtr->ProjNum[i_new] - pMiddle;        
        
        
        
        
        Real_t weightValueArray[sum*z_overlap_num];
        Real_t errorValueArray[sum*z_overlap_num];
        bool selectValueArray[sum*z_overlap_num];
        Real_t projectionValueArray[sum];
        Real_t* weightValueArrayPointer=&weightValueArray[0];
        Real_t* errorValueArrayPointer=&errorValueArray[0];
        bool* selectValueArrayPointer=&selectValueArray[0];
        Real_t* projectionValueArrayPointer=&projectionValueArray[0];
        

        fprintf(TomoInputsPtr->debug_file_ptr, "reach 1 \n");
        
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
	  if(p==0){
            int32_t i_rBeginning=(AMatrixPtr[pMiddle].index[0]);	                   
            Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            Real_t* errorSinoValueArray=errorSinoValueArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            bool* ProjSelectArray=ProjSelectArrayBegin+pMiddle*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 
 	    memcpy(projectionValueArrayPointer,&(AMatrixPtr[pMiddle].values[0]),sizeof(Real_t)*AMatrixPtr[pMiddle].count);
 	    projectionValueArrayPointer=projectionValueArrayPointer+AMatrixPtr[pMiddle].count;                  
	    for (q = 0; q < AMatrixPtr[pMiddle].count; q++)
	    {
	      memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	      weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              TomoInputsWeightArray=TomoInputsWeightArray+distance;	      
	      memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	      errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              errorSinoValueArray=errorSinoValueArray+distance;         	      
	      memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
              selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              ProjSelectArray=ProjSelectArray+distance;
	    }
	  }
	  else{
	    if((pMiddle-p)>=0 && (pMiddle+p)<ScannedObjectPtr->ProjNum[i_new]){
                int32_t i_rBeginning=(AMatrixPtr[pMiddle-p].index[0]);
            	int32_t size=AMatrixPtr[pMiddle-p].count;	                   
            	Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+(pMiddle-p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	Real_t* TomoInputsWeightArray2=TomoInputsWeightArray+2*p*NtNrMul;
            	Real_t* errorSinoValueArray=errorSinoValueArrayBegin+(pMiddle-p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	Real_t* errorSinoValueArray2=errorSinoValueArray+2*p*NtNrMul;
            	bool* ProjSelectArray=ProjSelectArrayBegin+(pMiddle-p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	bool* ProjSelectArray2=ProjSelectArray+2*p*NtNrMul;
 
 	    	memcpy(projectionValueArrayPointer,&(AMatrixPtr[pMiddle-p].values[0]),sizeof(Real_t)*size);
 	    	projectionValueArrayPointer=projectionValueArrayPointer+size;
 	    	memcpy(projectionValueArrayPointer,&(AMatrixPtr[pMiddle+p].values[0]),sizeof(Real_t)*size);
 	    	projectionValueArrayPointer=projectionValueArrayPointer+size; 	                      
	    	for (q = 0; q < size; q++)
	    	{
	      		memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	      		weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              		TomoInputsWeightArray=TomoInputsWeightArray+distance;
	      		memcpy(weightValueArrayPointer,TomoInputsWeightArray2,sizeof(Real_t)*z_overlap_num);
	      		weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              		TomoInputsWeightArray2=TomoInputsWeightArray2+distance;  
              
                          	      
	      		memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	      		errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              		errorSinoValueArray=errorSinoValueArray+distance;
	      		memcpy(errorValueArrayPointer,errorSinoValueArray2,sizeof(Real_t)*z_overlap_num);
	      		errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              		errorSinoValueArray2=errorSinoValueArray2+distance;               
              
                      	      
	      		memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
              		selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              		ProjSelectArray=ProjSelectArray+distance;
	      		memcpy(selectValueArrayPointer,ProjSelectArray2,sizeof(bool)*z_overlap_num);	      
              		selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              		ProjSelectArray2=ProjSelectArray2+distance;
              	}	              
            }
            else if((pMiddle-p)>=0 && (pMiddle+p)>=ScannedObjectPtr->ProjNum[i_new]){
             	int32_t i_rBeginning=(AMatrixPtr[pMiddle-p].index[0]);
             	int32_t size=AMatrixPtr[pMiddle-p].count;	                   
            	Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+(pMiddle-p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	Real_t* errorSinoValueArray=errorSinoValueArrayBegin+(pMiddle-p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	bool* ProjSelectArray=ProjSelectArrayBegin+(pMiddle-p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 	    	memcpy(projectionValueArrayPointer,&(AMatrixPtr[pMiddle-p].values[0]),sizeof(Real_t)*size);
 	    	projectionValueArrayPointer=projectionValueArrayPointer+size;
	    	for (q = 0; q < size; q++)
	    	{
	      		memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	      		weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              		TomoInputsWeightArray=TomoInputsWeightArray+distance;
	      		memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	      		errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              		errorSinoValueArray=errorSinoValueArray+distance; 
	      		memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
              		selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              		ProjSelectArray=ProjSelectArray+distance;              		             			    	
	    	} 	    	            	            	            	             		                              
            }
            else if ((pMiddle-p)<0 && (pMiddle+p)<ScannedObjectPtr->ProjNum[i_new]){
             	int32_t i_rBeginning=(AMatrixPtr[pMiddle+p].index[0]);
             	int32_t size=AMatrixPtr[pMiddle+p].count;	                   
            	Real_t* TomoInputsWeightArray=TomoInputsWeightArrayBegin+(pMiddle+p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	Real_t* errorSinoValueArray=errorSinoValueArrayBegin+(pMiddle+p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
            	bool* ProjSelectArray=ProjSelectArrayBegin+(pMiddle+p)*NtNrMul+(i_rBeginning)*SinogramPtr->N_t;
 	    	memcpy(projectionValueArrayPointer,&(AMatrixPtr[pMiddle+p].values[0]),sizeof(Real_t)*size);
 	    	projectionValueArrayPointer=projectionValueArrayPointer+size;
	    	for (q = 0; q < size; q++)
	    	{
	      		memcpy(weightValueArrayPointer,TomoInputsWeightArray,sizeof(Real_t)*z_overlap_num);
	      		weightValueArrayPointer=weightValueArrayPointer+z_overlap_num;
              		TomoInputsWeightArray=TomoInputsWeightArray+distance;
	      		memcpy(errorValueArrayPointer,errorSinoValueArray,sizeof(Real_t)*z_overlap_num);
	      		errorValueArrayPointer=errorValueArrayPointer+z_overlap_num;
              		errorSinoValueArray=errorSinoValueArray+distance; 
	      		memcpy(selectValueArrayPointer,ProjSelectArray,sizeof(bool)*z_overlap_num);	      
              		selectValueArrayPointer=selectValueArrayPointer+z_overlap_num;
              		ProjSelectArray=ProjSelectArray+distance;              		             			    	
	    	}             
            }
            else{
              break;
            }              
	  }
	  
        }
            for(p=0;p<sum;p++)
	      projectionValueArray[p]=projectionValueArray[p]*deltar;               
        
        fprintf(TomoInputsPtr->debug_file_ptr, "reach 2 \n");        	        	

        
       
        weightValueArrayPointer=&weightValueArray[0];
        errorValueArrayPointer=&errorValueArray[0];
        selectValueArrayPointer=&selectValueArray[0];
        projectionValueArrayPointer=&projectionValueArray[0];

	for (q = 0; q < sum; q++)
	{
	        Real_t projectionEntry=*projectionValueArrayPointer;
		for (r = 0; r < z_overlap_num; r++)
		{
			if (*selectValueArrayPointer )
			{

	           		THETA2 += (projectionEntry*projectionEntry*(*weightValueArrayPointer));
               			THETA1 += -((*errorValueArrayPointer)*projectionEntry*(*weightValueArrayPointer));
            		}
			else
			{
				THETASelTemp = errorSinoThresh*errorSinoDelta*sqrt(*weightValueArrayPointer)/fabs(*errorValueArrayPointer);
	            		THETA2 += (projectionEntry*projectionEntry*THETASelTemp);
            			THETA1 += -((*errorValueArrayPointer)*projectionEntry*THETASelTemp);
			}
			weightValueArrayPointer++;
			errorValueArrayPointer++;
			selectValueArrayPointer++;
            	}
            	projectionValueArrayPointer++;	    
	}
        fprintf(TomoInputsPtr->debug_file_ptr, "reach 3 \n");               

        
        UpdatedVoxelValue = CE_FunctionalSubstitution(V, THETA1, THETA2, ScannedObjectPtr, TomoInputsPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag);
              

	

	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
		int32_t sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
		for (q = 0; q < AMatrixPtr[p].count; q++)
        	{
               	    	int32_t i_r = (AMatrixPtr[p].index[q]);
        	    	Real_t ProjectionEntry = (AMatrixPtr[p].values[q]*SinogramPtr->delta_r);
			for (r = 0; r < z_overlap_num; r++)
			{ 
				int32_t i_t = slice*z_overlap_num + r;
	        		ErrorSino[sino_view][i_r][i_t] -= (ProjectionEntry*(UpdatedVoxelValue - V));
	   			if (fabs(ErrorSino[sino_view][i_r][i_t]*sqrt(TomoInputsPtr->Weight[sino_view][i_r][i_t])) < TomoInputsPtr->ErrorSinoThresh)
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = true;
				else
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = false;
	   		}
		}
	}	 	
	
	return UpdatedVoxelValue;	

   */
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



Real_t compute_voxel_update (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t*** ErrorSino, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new,Real_t* projectionValueArrayPointer,Real_t* weightValueArrayPointer,int32_t slice_begin,int32_t slice_end,bool* selectValueArrayPointer, Real_t* errorValueArrayPointer,Real_t V,Real_t THETA1,Real_t THETA2)
{
#ifdef PHASE_CONTRAST_TOMOGRAPHY
	return compute_voxel_update_AMat2D (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, AMatrixPtr, VoxelLineResponse, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag, i_new, slice, j_new, k_new);
#else
	return compute_voxel_update_AMat1D (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, AMatrixPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag, i_new, slice, j_new, k_new,projectionValueArrayPointer,weightValueArrayPointer,slice_begin,slice_end,selectValueArrayPointer,errorValueArrayPointer,V,THETA1,THETA2);
#endif
}
