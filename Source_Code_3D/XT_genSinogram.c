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





/*#include <iostream>*/
/*#include "TiffUtilities.h"*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "XT_Structures.h"
#include "XT_Constants.h"
#include "allocate.h"
#include <math.h>
#include "XT_IOMisc.h"
#include "XT_AMatrix.h"
#include "XT_Profile.h"
#include "randlib.h"

/*function read16bitTiff() is defined in TiffUtilities.cpp*/
unsigned short* read16bitTiff(char* file);

/*'gen_projection_filename' generates the name of the projection files given to us complete with the folder heirarchy.*/
void gen_projection_filename (char *file, int time_step, int slice_num)
{
	char time_str1[TIME_DIGITS+1]="\0",time_str2[TIME_DIGITS+1]="\0";
	char slice_str1[ANGLE_DIGITS+1]="\0",slice_str2[ANGLE_DIGITS+1]="\0";
	int buf,num_dig;

	for (buf=time_step,num_dig=0; buf!=0; buf=buf/10) 
		num_dig=num_dig+1;

	sprintf(time_str1,"%d",time_step);	
	strncat(time_str2, REPEAT_DIGITS, TIME_DIGITS-num_dig);

	for (buf=slice_num,num_dig=0; buf!=0; buf=buf/10) 
		num_dig=num_dig+1;
	
	sprintf(slice_str1,"%d",slice_num);	
	strncat(slice_str2, REPEAT_DIGITS, ANGLE_DIGITS-num_dig);

	sprintf(file, "%s%s%s%s_%s%s%s%s_%s%s.tif",FOLDER_PATH,FILE_FIRSTNAME,time_str2,time_str1,PROJECTION_INNER_FOLDER,FILE_FIRSTNAME,time_str2,time_str1,slice_str2,slice_str1);
	/*fprintf(TomoInputsPtr->debug_file_ptr, "\n%s,%s,%d\n",REPEAT_DIGITS,file,TIME_DIGITS-time_step/10-1);		*/
}


void initPhantomStructures(Sinogram* Sino, ScannedObject* ScanObj, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	*Sino = *SinogramPtr;
	*ScanObj = *ScannedObjectPtr;

    	ScanObj->N_x = TomoInputsPtr->phantom_N_xy;
        ScanObj->N_y = TomoInputsPtr->phantom_N_xy;
    	ScanObj->N_z = TomoInputsPtr->phantom_N_z/SinogramPtr->total_t_slices*SinogramPtr->N_t;
    
    	ScanObj->delta_xy = Sino->Length_R/ScanObj->N_x;
    	ScanObj->delta_z = Sino->Length_T/ScanObj->N_z;
 
	Sino->OffsetR = (ScanObj->delta_xy/sqrt(2.0)+Sino->delta_r/2.0)/DETECTOR_RESPONSE_BINS;
	Sino->OffsetT = ((ScanObj->delta_z/2) + Sino->delta_t/2)/DETECTOR_RESPONSE_BINS;
}

void genSinogramFromPhantom (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	FILE *fp;
	Sinogram Phantom_Sino;
	ScannedObject Phantom_ScanObj;
	AMatrixCol *VoxelLineResponse;
	long int stream_offset;
	int32_t i, j, k, m, n, idx, size, t, slice, result; 
        float ***object;
	Real_t pixel, val, **H_r, *H_t;
  	uint8_t AvgNumXElements, AvgNumZElements;
	char phantom_file[1000] = PHANTOM_FILENAME;
	char projection_file[] = PROJECTION_FILENAME;
	char weight_file[] = WEIGHT_MATRIX_FILENAME;
	char detect_file[] = "detector_forwardproj";
	int dimTiff[4];

	if (TomoInputsPtr->phantom_N_xy % SinogramPtr->N_r != 0 || TomoInputsPtr->phantom_N_z % SinogramPtr->total_t_slices != 0){
		printf("ERROR: genSinogramFromPhantom: N_r = %d does not divide phantom_N_xy = %d or N_t = %d does not divide phantom_N_z = %d\n", SinogramPtr->N_r, TomoInputsPtr->phantom_N_xy, SinogramPtr->total_t_slices, TomoInputsPtr->phantom_N_z);
		exit (1);
	}
	
	/*printf("\nsizes = %d, %d, %d, %d\n", sizeof(int), sizeof(int64_t), sizeof(long int), sizeof(long long));*/

	initPhantomStructures(&Phantom_Sino, &Phantom_ScanObj, SinogramPtr, ScannedObjectPtr, TomoInputsPtr);	

	AMatrixCol* AMatrixPtr = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
	/*AvgNumXElements over estimates the total number of entries in a single column of A matrix when indexed by both voxel and angle*/
  	AvgNumXElements = (uint8_t)ceil(3*Phantom_ScanObj.delta_xy/(Phantom_Sino.delta_r));
  	AMatrixPtr->values = (Real_t*)get_spc((int32_t)AvgNumXElements,sizeof(Real_t));
  	AMatrixPtr->index  = (int32_t*)get_spc((int32_t)AvgNumXElements,sizeof(int32_t));

	object = (float***)multialloc(sizeof(float), 3, TomoInputsPtr->phantom_N_xy, TomoInputsPtr->phantom_N_xy, TomoInputsPtr->phantom_N_z);
	memset(&(Phantom_Sino.Projection[0][0][0]), 0, Phantom_Sino.N_p*Phantom_Sino.N_t*Phantom_Sino.N_r*sizeof(Real_t));	

	H_r = (Real_t **)multialloc(sizeof(Real_t), 2, Phantom_Sino.N_p, DETECTOR_RESPONSE_BINS+1);
	H_t = (Real_t *)get_spc(DETECTOR_RESPONSE_BINS + 1, sizeof(Real_t));
	DetectorResponseProfile (H_r, H_t, &Phantom_Sino, &Phantom_ScanObj, TomoInputsPtr);
	
  	AvgNumZElements = (uint8_t)((Phantom_ScanObj.delta_z/Phantom_Sino.delta_t) + 2);
	
	VoxelLineResponse = (AMatrixCol*)get_spc(Phantom_ScanObj.N_z,sizeof(AMatrixCol));
	for (t = 0; t < Phantom_ScanObj.N_z; t++){
    		VoxelLineResponse[t].values = (Real_t*)get_spc(AvgNumZElements, sizeof(Real_t));
    		VoxelLineResponse[t].index = (int32_t*)get_spc(AvgNumZElements, sizeof(int32_t));
	}
	storeVoxelLineResponse(H_t, VoxelLineResponse, &Phantom_ScanObj, &Phantom_Sino);
	
	sprintf(phantom_file, "%s%s.bin", PHANTOM_FOLDER, PHANTOM_FILENAME);
	fp = fopen (phantom_file, "r" );
 	if (fp==NULL) {fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: genSinogramFromPhantom: error in reading file %s\n",phantom_file); exit (1);}		
	size=TomoInputsPtr->phantom_N_xy*TomoInputsPtr->phantom_N_xy*TomoInputsPtr->phantom_N_z;

	fprintf(TomoInputsPtr->debug_file_ptr, "DebugMsg: genSinogramFromPhantom: Forward projecting phantoms - \n");	
	for (i=0; i<Phantom_Sino.N_p; i++){
		stream_offset = (long int)i*(long int)TomoInputsPtr->phantom_N_z*(long int)Phantom_ScanObj.N_y*(long int)Phantom_ScanObj.N_x;
		stream_offset += (long int)(TomoInputsPtr->phantom_N_z/SinogramPtr->total_t_slices)*(long int)SinogramPtr->slice_begin;
		result = fseek (fp, stream_offset*sizeof(float), SEEK_SET);

  		if (result != 0) 
		{fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: Error in seeking file %s, i = %d, stream_offset = %ld\n",phantom_file,i,stream_offset);}

		result = fread (&(object[0][0][0]), sizeof(float), size, fp);
  		if (result != size) 
		{fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: Reading file %s, Number of elements read does not match required, number of elements read=%d\n",phantom_file,result);}

		for (j=0; j<TomoInputsPtr->phantom_N_xy; j++)
		for (k=0; k<TomoInputsPtr->phantom_N_xy; k++){	
	   	    	calcAMatrixColumnforAngle(&Phantom_Sino, &Phantom_ScanObj, H_r, AMatrixPtr, j, k, i); 
                	for (slice=0; slice<TomoInputsPtr->phantom_N_z; slice++){
			    	pixel = (Real_t)object[j][k][slice];
	     	          	for (m=0; m<AMatrixPtr->count; m++){
                            		idx=AMatrixPtr->index[m];
                            		val=AMatrixPtr->values[m];
                            		for (n=0; n<VoxelLineResponse[slice].count; n++)
                                    		Phantom_Sino.Projection[i][idx][VoxelLineResponse[slice].index[n]]+=pixel*val*VoxelLineResponse[slice].values[n];
	     			}
			}
	  	 }
	}

	dimTiff[0] = 1; dimTiff[1] = 1; dimTiff[2] = SinogramPtr->N_p; dimTiff[3] = DETECTOR_RESPONSE_BINS+1;
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (detect_file, dimTiff, 0, 1, 2, 3, &(H_r[0][0]), 0, TomoInputsPtr->debug_file_ptr);
        
	free(AMatrixPtr->values);
	free(AMatrixPtr->index);
        free(VoxelLineResponse->values);
        free(VoxelLineResponse->index);
	multifree(H_r,2);
	free(H_t);
        free(AMatrixPtr);
        free(VoxelLineResponse);
	fclose(fp);	
	multifree(object,3);

	for (i=0; i < SinogramPtr->N_p; i++)
	for (j=0; j < SinogramPtr->N_r; j++){
	for (slice=0; slice < SinogramPtr->N_t; slice++)
		val = SinogramPtr->Projection[i][j][slice];
		val = EXPECTED_COUNTS_FOR_PHANTOM_DATA*exp(-val);
		/*TomoInputsPtr->Weight[i][j] = val + sqrt(val)*random2();*/
		if (TomoInputsPtr->No_Projection_Noise == 1)
			TomoInputsPtr->Weight[i][j][slice] = fabs(val);
		else
			TomoInputsPtr->Weight[i][j][slice] = fabs(val + sqrt(val)*normal());
		
		SinogramPtr->Projection[i][j][slice] = log(EXPECTED_COUNTS_FOR_PHANTOM_DATA/TomoInputsPtr->Weight[i][j][slice]);			
	}

	Write2Bin (projection_file, 1, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t, &(SinogramPtr->Projection[0][0][0]), TomoInputsPtr->debug_file_ptr);
	dimTiff[0] = 1; dimTiff[1] = SinogramPtr->N_p; dimTiff[2] = SinogramPtr->N_r; dimTiff[3] = SinogramPtr->N_t;
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (projection_file, dimTiff, 0, 3, 1, 2, &(SinogramPtr->Projection[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	Write2Bin (weight_file, 1, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t, &(TomoInputsPtr->Weight[0][0][0]), TomoInputsPtr->debug_file_ptr);
}

/*'genSinogram_fromBin' initializes the sinogram and forward model weights from a bin file generated in a previous run. */
int genSinogram_fromBin(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	char sinofile[100]=PROJECTION_FILENAME;
	char weightfile[100]=WEIGHT_MATRIX_FILENAME;
	Real_t *weight, *pixels;
	int32_t i,k,size,slice;
	Real_t sino_avg=0, weight_avg=0;

	sprintf (sinofile, "%s_n%d", sinofile, TomoInputsPtr->node_rank);
	sprintf (weightfile, "%s_n%d", weightfile, TomoInputsPtr->node_rank);

	size=SinogramPtr->N_p*SinogramPtr->N_r*SinogramPtr->N_t;
	pixels = (Real_t*)get_spc(size,sizeof(Real_t));
	weight = (Real_t*)get_spc(size,sizeof(Real_t));

	Read4mBin (sinofile, 1, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t, &(SinogramPtr->Projection[0][0][0]), TomoInputsPtr->debug_file_ptr);
	Read4mBin (weightfile, 1, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t, &(TomoInputsPtr->Weight[0][0][0]), TomoInputsPtr->debug_file_ptr);

	for (i=0; i<SinogramPtr->N_p; i++){
		for (k=0; k<SinogramPtr->N_r; k++){
			for (slice=0; slice<SinogramPtr->N_t; slice++){
				sino_avg += SinogramPtr->Projection[i][k][slice];
				weight_avg += TomoInputsPtr->Weight[i][k][slice];
			}
		}
	}

	sino_avg = sino_avg/(SinogramPtr->N_p*SinogramPtr->N_r*SinogramPtr->N_t);
	weight_avg = weight_avg/(SinogramPtr->N_p*SinogramPtr->N_r*SinogramPtr->N_t);

	fprintf(TomoInputsPtr->debug_file_ptr, "DebugMsg: genSinogram_frombin: sinogram average = %f, Weight average = %f\n", sino_avg, weight_avg);

	free(pixels);
	free(weight);
	return(0);	
}

int genSinogram_subsampleCounts(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	FILE *fp;
	int result;
	char brightfile[100]=BRIGHT_FIELD_FILENAME;
	char weightfile[100]=WEIGHT_MATRIX_FILENAME;
	Real_t temp, ***weight, **bright, weight_temp, bright_temp;
	int32_t p, q, N_t, N_r, mult_xy, mult_z, i,k,size,slice;
	Real_t sino_avg=0, weight_avg=0;
	
	sprintf (brightfile, "%s_n%d", brightfile, TomoInputsPtr->node_rank);
	sprintf (weightfile, "%s_n%d", weightfile, TomoInputsPtr->node_rank);

	N_t = SinogramPtr->N_t;
	N_r = SinogramPtr->N_r;
	mult_xy = ScannedObjectPtr->mult_xy;
	mult_z = ScannedObjectPtr->mult_z;
	
	sprintf(brightfile, "%s.bin", brightfile);
	fp = fopen (brightfile, "r" );
 	if (fp==NULL) {fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: genSinogram_subsampleCounts: error in reading file %s\n",brightfile); exit (1);}		
	size=N_r*N_t*mult_xy*mult_z;
	bright = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r*mult_xy, SinogramPtr->N_t*mult_z);
	result = fread (&(bright[0][0]), sizeof(Real_t), size, fp);
  	if (result != (int)size) 
	{fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: Reading file %s, Number of elements read does not match required, number of elements read=%d\n",brightfile,result);}
	fclose(fp);	

	sprintf(weightfile, "%s.bin", weightfile);
	fp = fopen (weightfile, "r" );
 	if (fp==NULL) {fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: error in reading file %s\n",weightfile); exit (1);}		
	size=SinogramPtr->N_p*N_r*N_t*mult_xy*mult_z;
	weight = (Real_t***)multialloc(sizeof(Real_t), 3, SinogramPtr->N_p, SinogramPtr->N_r*mult_xy, SinogramPtr->N_t*mult_z);
	result = fread (&(weight[0][0][0]), sizeof(Real_t), size, fp);
  	if (result != (int)size) 
	{fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: Reading file %s, Number of elements read does not match required, number of elements read=%d\n",weightfile,result);}
	fclose(fp);

	for (i=0; i<SinogramPtr->N_p; i++)
	{
		for (k=0; k<N_r; k++)
		{
			for (slice=0; slice<N_t; slice++)
			{
				weight_temp = 0; bright_temp = 0;
				for (p = 0; p < mult_xy; p++)
				{
					for (q = 0; q < mult_z; q++)
					{
						weight_temp += weight[i][k*mult_xy+p][slice*mult_z+q];
						bright_temp += bright[k*mult_xy+p][slice*mult_z+q];
					}
				}
				temp = log(bright_temp/weight_temp);
				SinogramPtr->Projection[i][k][slice] = BH_QUAD_COEF*temp*temp + temp;
				
				TomoInputsPtr->Weight[i][k][slice] = weight_temp;
				sino_avg += SinogramPtr->Projection[i][k][slice];
				weight_avg += TomoInputsPtr->Weight[i][k][slice];
			}
		}
	}

	sino_avg = sino_avg/(SinogramPtr->N_p*SinogramPtr->N_r*SinogramPtr->N_t);
	weight_avg = weight_avg/(SinogramPtr->N_p*SinogramPtr->N_r*SinogramPtr->N_t);

	fprintf(TomoInputsPtr->debug_file_ptr, "DebugMsg: genSinogram_fromCounts: sinogram average = %f, Weight average = %f\n", sino_avg, weight_avg);

	multifree(bright,2);
	multifree(weight,3);
	return(0);	
}

