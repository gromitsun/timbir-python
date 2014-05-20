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



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "allocate.h"
#include "XT_Constants.h"
#include "randlib.h"
#include <getopt.h>
#include "XT_Structures.h"
#include <ctype.h>
#include "XT_IOMisc.h"
#include "XT_HDFIO.h"

/*Parses the input text file and extracts words delimited by any non-alphanumeric character.
In this code, it is used to parse view_info.txt and extract the time and angle information.
--Inputs--
fp - Text file pointer
str - pointer to memory where extracted text is stored
numbytes - number of bytes read
--Outputs--
returns the last character which acted as a delimiter for the extracted text*/
char readFileStream(FILE **fp, char* str, int32_t *numbytes)
{
	char lastchar;
	
	lastchar = EOF;
	while (!feof(*fp)){
		*str = fgetc(*fp);
		(*numbytes)++;
		if(isspace(*str) == 0){
			if (isalnum(*str) == 0 && *str != '.')
			{
				lastchar = *str;
				*str = '\0';
				break;
			}
			str++;
		}
	}
	return (lastchar);
}

/*For each time slice in the reconstruction, the function copies the corresponding view indices to a new array (which is then usedafter copying). 
--Inputs--
time - index of time slice
ViewIndex - contains the indices of the views which are assumed to be the forward projections of reconstruction at index 'time'
ViewNum - Number of such views*/
void copyViewIndexMap (ScannedObject* ScannedObjectPtr, int32_t time, int32_t* ViewIndex, int32_t ViewNum)
{
	int32_t i;
        ScannedObjectPtr->ProjIdxPtr[time] = (int32_t*)get_spc(ViewNum, sizeof(int32_t));
	for(i=0; i<ViewNum; i++){
		ScannedObjectPtr->ProjIdxPtr[time][i] = ViewIndex[i];
	}
	ScannedObjectPtr->ProjNum[time] = ViewNum;
#ifdef DEBUG_EN
	fprintf(TomoInputsPtr->debug_file_ptr, "DebugMsg: copyViewIndexMap: Total number of views at time %d relative to object is %d\n", time, ViewNum);
#endif
}

/*Maps the projections to reconstruction time slices.
A projection at time 'projection_time' is assigned to a certain time slice in the reconstruction
 if the start and end time of the time slice includes the 'projection_time'. */
void initSparseAnglesOfObject(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	int32_t  i;
	int32_t *ViewIndex, ViewNum;
	Real_t projection_time, object_time2; 
	int32_t object_time_idx = 0;

	ViewIndex = (int32_t*)get_spc(SinogramPtr->N_p, sizeof(int32_t));
	ViewNum = 0;
	
	object_time2 = ScannedObjectPtr->Rtime0 + ScannedObjectPtr->delta_Rtime;
	for (i=0; i<SinogramPtr->N_p; i++){
			projection_time = SinogramPtr->TimePtr[i];
			while ( projection_time >= object_time2 && object_time_idx + 1 < ScannedObjectPtr->N_time)
			{
				copyViewIndexMap (ScannedObjectPtr, object_time_idx, ViewIndex, ViewNum);
			/*	object_time1 = object_time2;*/
				object_time2 += ScannedObjectPtr->delta_Rtime;
				object_time_idx++;
				ViewNum = 0;	
			}
			ViewIndex[ViewNum] = i;
			ViewNum++;
	}

	copyViewIndexMap (ScannedObjectPtr, object_time_idx, ViewIndex, ViewNum);
	free(ViewIndex);
}

/*Function which parses view_info.txt containing information about the views
and the corresponding times at which they were acquired.
All text before '-' is considered to be the time and all text before a ',' or '\n' is considered to the view angle*/
void initSparseAnglesfrmFile(Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr)
{
        char filename[] = SPARSE_ANGLES_LIST_FILE;
        char str[100], lastchar;
        int32_t numbytes = 0, index = 0;
        Real_t time = 0;
        FILE *fp;

        fp = fopen (filename, "r" );
        if (fp==NULL) {fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: Error in reading file %s\n", filename); exit(1);}

        while ((lastchar = readFileStream(&fp, str, &numbytes)) != EOF){
                if (lastchar == '-')
                        time = (Real_t)atof(str);
                else if(lastchar == ',' || lastchar == '\n'){
                        if(str[0] != '\0'){
                                SinogramPtr->ViewPtr[index] = (Real_t)atof(str)*M_PI/180.0;
                                SinogramPtr->TimePtr[index] = time;
                                index++;
                        }
                }
                else
                {
                        fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: initSparseAnglesfrmFile: File %s is not correctly formatted to be read\n", filename);
                        exit(1);
                }
        }

	if(index != SinogramPtr->N_p){
		fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: Number of projections read is %d and expected is %d\n", index, SinogramPtr->N_p);
		exit(1);
	}
}

/*Populates the Views and times of each projection from the text file view_info.txt into arrays*/
void initRandomAngles (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	int32_t i;
	int32_t sino_view;
        SinogramPtr->ViewPtr = (Real_t*)get_spc(SinogramPtr->N_p, sizeof(Real_t));
        SinogramPtr->TimePtr = (Real_t*)get_spc(SinogramPtr->N_p, sizeof(Real_t));
        ScannedObjectPtr->ProjIdxPtr = (int32_t**)get_spc(ScannedObjectPtr->N_time, sizeof(int32_t*));
	ScannedObjectPtr->ProjNum = (int32_t*)get_spc(ScannedObjectPtr->N_time, sizeof(int32_t));
	memset(&(ScannedObjectPtr->ProjNum[0]), 0, ScannedObjectPtr->N_time*sizeof(int32_t));

        initSparseAnglesfrmFile(SinogramPtr, TomoInputsPtr);
	initSparseAnglesOfObject(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	
	int32_t k;
#ifdef DEBUG_EN
        fprintf(TomoInputsPtr->debug_file_ptr, "initRandomAngles: The initialized angle indices of sinogram are ...\n");
        for (i=0; i<SinogramPtr->N_p; i++){
                    fprintf(TomoInputsPtr->debug_file_ptr, "time %f - ", SinogramPtr->TimePtr[i]);
                    fprintf(TomoInputsPtr->debug_file_ptr, "%f, ", SinogramPtr->ViewPtr[i]);
                }
                fprintf(TomoInputsPtr->debug_file_ptr, "\n");
#endif
	fprintf(TomoInputsPtr->debug_file_ptr, "initRandomAngles: The initialized angle indices of sinogram as corresponding to the scanned object are ...\n");
	for (i=0; i<ScannedObjectPtr->N_time; i++){
		fprintf(TomoInputsPtr->debug_file_ptr, "Object %d : ", i);
		for (k=0; k<ScannedObjectPtr->ProjNum[i]; k++){
			sino_view = ScannedObjectPtr->ProjIdxPtr[i][k];
			fprintf(TomoInputsPtr->debug_file_ptr, "%.1f, ", SinogramPtr->ViewPtr[sino_view]*180/M_PI);
		}
		fprintf(TomoInputsPtr->debug_file_ptr, "\n");
	}	

}

/*Computes the Euclidean distance of a voxel to its neighboring voxels
--Inputs--
i, j, k, l are the indices of neighoring voxels
--Outputs--
Returns the distance */
Real_t distance2node(uint8_t i, uint8_t j, uint8_t k, uint8_t l)
{
	return(sqrt(pow((Real_t)i-1.0, 2.0)+pow((Real_t)j-1.0, 2.0)+pow((Real_t)k-1.0, 2.0)+pow((Real_t)l-1.0, 2.0)));
}

/*Initializes the weights w_{ij} used in the qGGMRF models*/
void initFilter (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	uint8_t i,j,k;
	Real_t temp1,sum=0,prior_const=0;
	prior_const = ScannedObjectPtr->delta_xy*ScannedObjectPtr->delta_xy*ScannedObjectPtr->delta_xy*ScannedObjectPtr->delta_Rtime;
/*Filter coefficients of neighboring pixels are inversely proportional to the distance from the center pixel*/
	TomoInputsPtr->Time_Filter[0] = 1.0/distance2node(0,1,1,1);
	sum += 2.0*TomoInputsPtr->Time_Filter[0];
	
	for (i=0; i<NHOOD_Y_MAXDIM; i++)
	for (j=0; j<NHOOD_X_MAXDIM; j++)
	for (k=0; k<NHOOD_Z_MAXDIM; k++){
	if(i!=(NHOOD_Y_MAXDIM)/2 || j!=(NHOOD_X_MAXDIM-1)/2 || k!=(NHOOD_Z_MAXDIM-1)/2)
	{
		temp1 = 1.0/distance2node(1,i,j,k);
		TomoInputsPtr->Spatial_Filter[i][j][k] = temp1;
		sum=sum+temp1;
	}
	else
		TomoInputsPtr->Spatial_Filter[i][j][k]=0;
	}

	for (i=0; i<NHOOD_Y_MAXDIM; i++)
	for (j=0; j<NHOOD_X_MAXDIM; j++)
	for (k=0; k<NHOOD_Z_MAXDIM; k++){
		TomoInputsPtr->Spatial_Filter[i][j][k] = prior_const*TomoInputsPtr->Spatial_Filter[i][j][k]/sum;
	}

	TomoInputsPtr->Time_Filter[0] = prior_const*TomoInputsPtr->Time_Filter[0]/sum;

#ifdef DEBUG_EN
	sum=0;
	for (i=0; i<NHOOD_Y_MAXDIM; i++)
		for (j=0; j<NHOOD_X_MAXDIM; j++)
			for (k=0; k<NHOOD_Z_MAXDIM; k++)
			{
				sum+=TomoInputsPtr->Spatial_Filter[i][j][k];
				fprintf(TomoInputsPtr->debug_file_ptr, "initFilter: Filter i=%d, j=%d, k=%d, coeff = %f\n", i,j,k,TomoInputsPtr->Spatial_Filter[i][j][k]/prior_const);
			}
			fprintf(TomoInputsPtr->debug_file_ptr, "initFilter: Filter i=0 is %f\n", TomoInputsPtr->Time_Filter[0]/prior_const);
			fprintf(TomoInputsPtr->debug_file_ptr, "initFilter: Sum of filter coefficients is %f\n",(sum+2.0*TomoInputsPtr->Time_Filter[0])/prior_const);	
			fprintf(TomoInputsPtr->debug_file_ptr, "initFilter: delta_xy*delta_xy*delta_z*delta_tau = %f\n",prior_const);	
#endif /*#ifdef DEBUG_EN*/


}

/*Initializes the sines and cosines of angles at which projections are acquired. It is then used when computing the voxel profile*/
void calculateSinCos(Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr)
{
  int32_t i;

  SinogramPtr->cosine=(Real_t*)get_spc(SinogramPtr->N_p, sizeof(Real_t));
  SinogramPtr->sine=(Real_t*)get_spc(SinogramPtr->N_p, sizeof(Real_t));

  for(i=0;i<SinogramPtr->N_p;i++)
  {
    SinogramPtr->cosine[i]=cos(SinogramPtr->ViewPtr[i]);
    SinogramPtr->sine[i]=sin(SinogramPtr->ViewPtr[i]);
  }
  fprintf(TomoInputsPtr->debug_file_ptr, "calculateSinCos: Calculated sines and cosines of angles of rotation\n");
}


/*Initializes the variables in the three major structures used throughout the code -
Sinogram, ScannedObject, TomoInputs. It also allocates memory for several variables.*/
void initStructures (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	/*Initializing Sinogram parameters*/
	int32_t  prev_mult_xy, prev_mult_z, p, q, i, j, k, dim[4];
	Real_t **OffsetTemp;
	char projOffset_file[100] = PROJ_OFFSET_FILENAME;
	char VarEstFile[100] = "variance_estimate";
	Real_t Lap_Kernel[3][3] = {{-0.5,-1,-0.5},{-1,6,-1},{-0.5,-1,-0.5}};

	if (TomoInputsPtr->initICD != 0)
	{
		sprintf(VarEstFile, "%s_n%d", VarEstFile, TomoInputsPtr->node_rank);	
		Read4mBin (VarEstFile, 1, 1, 1, 1, &(TomoInputsPtr->var_est), TomoInputsPtr->debug_file_ptr);
	}

	sprintf(projOffset_file, "%s_n%d", projOffset_file, TomoInputsPtr->node_rank);	
	SinogramPtr->Length_T = SinogramPtr->Length_T/TomoInputsPtr->node_num;
	SinogramPtr->N_t = SinogramPtr->total_t_slices/TomoInputsPtr->node_num;

	if (SinogramPtr->N_t < 3)
	{
		printf ("ERROR: initStructures: Number of spatial slices per node cannot be less than 3\n");
		exit(1);
	}
	
	if (SinogramPtr->N_t % (int32_t)ScannedObjectPtr->delta_z != 0){
		printf("ERROR: initStructures: Cannot do reconstruction since delta_z = %d does not divide %d\n", (int32_t)ScannedObjectPtr->delta_z, SinogramPtr->N_t);
		exit(1);
	}

	if (SinogramPtr->N_r % (int32_t)ScannedObjectPtr->delta_xy != 0){
		printf("ERROR: initStructures: Cannot do reconstruction since delta_xy = %d does not divide %d\n", (int32_t)ScannedObjectPtr->delta_xy, SinogramPtr->N_r);
		exit(1);
	}
	
	SinogramPtr->delta_r = SinogramPtr->Length_R/(SinogramPtr->N_r);
	SinogramPtr->delta_t = SinogramPtr->Length_T/(SinogramPtr->N_t);
	SinogramPtr->R0 = -TomoInputsPtr->RotCenter*SinogramPtr->delta_r;
	SinogramPtr->RMax = (SinogramPtr->N_r-TomoInputsPtr->RotCenter)*SinogramPtr->delta_r;
	SinogramPtr->T0 = -SinogramPtr->Length_T/2.0;
	SinogramPtr->TMax = SinogramPtr->Length_T/2.0;
	
	/*Initializing parameters of the object to be reconstructed*/
	ScannedObjectPtr->Length_X = SinogramPtr->Length_R;
    	ScannedObjectPtr->Length_Y = SinogramPtr->Length_R;
	ScannedObjectPtr->Length_Z = SinogramPtr->Length_T;
    	ScannedObjectPtr->N_x = (int32_t)(SinogramPtr->N_r/ScannedObjectPtr->delta_xy);
	ScannedObjectPtr->N_y = (int32_t)(SinogramPtr->N_r/ScannedObjectPtr->delta_xy);
	ScannedObjectPtr->N_z = (int32_t)(SinogramPtr->N_t/ScannedObjectPtr->delta_z);	
	ScannedObjectPtr->delta_xy = ScannedObjectPtr->delta_xy*SinogramPtr->delta_r;
	ScannedObjectPtr->delta_z = ScannedObjectPtr->delta_z*SinogramPtr->delta_t;

	SinogramPtr->Projection = (Real_t***)multialloc(sizeof(Real_t), 3, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t);
	SinogramPtr->ProjSelect = (bool***)multialloc(sizeof(bool), 3, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t);
	for (i = 0; i < SinogramPtr->N_p; i++)
	for (j = 0; j < SinogramPtr->N_r; j++)
	for (k = 0; k < SinogramPtr->N_t; k++)
		SinogramPtr->ProjSelect[i][j][k] = true;

	SinogramPtr->ProjOffset = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
	TomoInputsPtr->Weight = (Real_t***)multialloc(sizeof(Real_t), 3, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t);
	memset(&(SinogramPtr->ProjOffset[0][0]), 0, SinogramPtr->N_t*SinogramPtr->N_r*sizeof(Real_t));

	prev_mult_xy = 1;
	prev_mult_z = 1;
	if (TomoInputsPtr->initICD == 3)
	{
		prev_mult_xy = 2;
		prev_mult_z = 2;
	}
	else if (TomoInputsPtr->initICD == 2)
		prev_mult_xy = 2;
	
	if (TomoInputsPtr->updateProjOffset == 1 || TomoInputsPtr->updateProjOffset == 3)
	{
		if (TomoInputsPtr->sinobin == 3)
		{
			OffsetTemp = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r/prev_mult_xy, SinogramPtr->N_t/prev_mult_z);
			Read4mBin (projOffset_file, 1, 1, SinogramPtr->N_r/prev_mult_xy, SinogramPtr->N_t/prev_mult_z, &(OffsetTemp[0][0]), TomoInputsPtr->debug_file_ptr);
			for (i = 0; i < SinogramPtr->N_r/prev_mult_xy; i++)
			for (j = 0; j < SinogramPtr->N_t/prev_mult_z; j++)
			for (p = 0; p < prev_mult_xy; p++)
			for (q = 0; q < prev_mult_z; q++)
				SinogramPtr->ProjOffset[i*prev_mult_xy + p][j*prev_mult_z + q] = OffsetTemp[i][j];
			multifree(OffsetTemp,2);	
		}	
		else
		{
			Read4mBin (projOffset_file, 1, 1, SinogramPtr->N_r, SinogramPtr->N_t, &(SinogramPtr->ProjOffset[0][0]), TomoInputsPtr->debug_file_ptr);
		}
		dim[0] = 1; dim[1] = 1; dim[2] = SinogramPtr->N_r; dim[3] = SinogramPtr->N_t;
		if (TomoInputsPtr->Write2Tiff == 1)
			WriteMultiDimArray2Tiff (projOffset_file, dim, 0, 1, 2, 3, &(SinogramPtr->ProjOffset[0][0]), 0, TomoInputsPtr->debug_file_ptr);
	}		

	if (ScannedObjectPtr->delta_xy != ScannedObjectPtr->delta_z)
		printf ("WARNING: initStructures: delta_xy is not equal to delta_z. The spatial invariance of prior does not hold.\n");

	ScannedObjectPtr->x0 = SinogramPtr->R0;
    	ScannedObjectPtr->z0 = SinogramPtr->T0;
    	ScannedObjectPtr->y0 = -ScannedObjectPtr->Length_Y/2.0;
    	ScannedObjectPtr->BeamWidth = SinogramPtr->delta_r; /*Weighting of the projections at different points of the detector*/
/*	ScannedObjectPtr->Object = (Real_t****)multialloc(sizeof(Real_t), 4, ScannedObjectPtr->N_time, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x, ScannedObjectPtr->N_z);*/
	ScannedObjectPtr->Object = (Real_t****)get_spc(ScannedObjectPtr->N_time, sizeof(Real_t***));
	for (i = 0; i < ScannedObjectPtr->N_time; i++)
	{
		ScannedObjectPtr->Object[i] = (Real_t***)multialloc(sizeof(Real_t), 3, ScannedObjectPtr->N_z + 2, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x);
	}
	
	/*OffsetR is stepsize of the distance between center of voxel of the object and the detector pixel, at which projections are computed*/
	SinogramPtr->OffsetR = (ScannedObjectPtr->delta_xy/sqrt(2.0)+SinogramPtr->delta_r/2.0)/DETECTOR_RESPONSE_BINS;
	SinogramPtr->OffsetT = ((ScannedObjectPtr->delta_z/2) + SinogramPtr->delta_t/2)/DETECTOR_RESPONSE_BINS;

	/*TomoInputs holds the input parameters and some miscellaneous variables*/
	TomoInputsPtr->Sigma_S_Q = pow((ScannedObjectPtr->Sigma_S*ScannedObjectPtr->delta_xy),MRF_Q);
	TomoInputsPtr->Sigma_S_Q_P = pow(ScannedObjectPtr->Sigma_S*ScannedObjectPtr->delta_xy,MRF_Q-ScannedObjectPtr->MRF_P);	
	TomoInputsPtr->Sigma_T_Q = pow((ScannedObjectPtr->Sigma_T*ScannedObjectPtr->delta_Rtime),MRF_Q);
	TomoInputsPtr->Sigma_T_Q_P = pow(ScannedObjectPtr->Sigma_T*ScannedObjectPtr->delta_Rtime,MRF_Q-ScannedObjectPtr->MRF_P);	

	TomoInputsPtr->num_z_blocks = 2*floor(TomoInputsPtr->num_threads/ScannedObjectPtr->N_time);
	if (TomoInputsPtr->num_z_blocks < 2)
		TomoInputsPtr->num_z_blocks = 2;
	else if (TomoInputsPtr->num_z_blocks > ScannedObjectPtr->N_z)
		TomoInputsPtr->num_z_blocks = ScannedObjectPtr->N_z;

	if (TomoInputsPtr->num_threads % TomoInputsPtr->num_z_blocks != 0)
		fprintf(TomoInputsPtr->debug_file_ptr, "WARNING: initStructures: Number of threads are not optimized for the number of z slices");


/*	TomoInputsPtr->BoundaryFlag = (uint8_t***)multialloc(sizeof(uint8_t), 3, 3, 3, 3);*/
        TomoInputsPtr->x_rand_select = (int32_t***)multialloc(sizeof(int32_t), 3, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y*ScannedObjectPtr->N_x);
        TomoInputsPtr->y_rand_select = (int32_t***)multialloc(sizeof(int32_t), 3, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y*ScannedObjectPtr->N_x);
        TomoInputsPtr->x_NHICD_select = (int32_t***)multialloc(sizeof(int32_t), 3, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y*ScannedObjectPtr->N_x);
        TomoInputsPtr->y_NHICD_select = (int32_t***)multialloc(sizeof(int32_t), 3, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y*ScannedObjectPtr->N_x);
        TomoInputsPtr->UpdateSelectNum = (int32_t**)multialloc(sizeof(int32_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
        TomoInputsPtr->NHICDSelectNum = (int32_t**)multialloc(sizeof(int32_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);


	ScannedObjectPtr->NHICD_Iterations = 10;
	for (i=0; i<=2; i++)
		for (j=0; j<=2; j++)
			SinogramPtr->Lap_Kernel[i][j] = Lap_Kernel[i][j]/(SinogramPtr->delta_r*SinogramPtr->delta_r);			

	fprintf(TomoInputsPtr->debug_file_ptr, "initStructures: Number of z blocks is %d\n", TomoInputsPtr->num_z_blocks);
	initFilter (ScannedObjectPtr, TomoInputsPtr);
	initRandomAngles (SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	calculateSinCos (SinogramPtr, TomoInputsPtr);
	fprintf(TomoInputsPtr->debug_file_ptr, "initStructures: Initilialized the structures, Sinogram and ScannedObject\n");
}


/*Function which parses the command line input to the C code and initializes several variables.*/
void argsParser (int argc, char **argv, Sinogram* SinogramPtr, ScannedObject *ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	int option_index;
	char c, debug_filename[100];
        Real_t slope_HU, c_HU;
	static struct option long_options[] =
        {
               {"p",       required_argument, 0, 'a'}, /*All 5 options below are prior model parameters*/
               {"sigma_s",  required_argument, 0, 'b'},
               {"c_s",  required_argument, 0, 'c'},
               {"sigma_t",  required_argument, 0, 'd'},
               {"c_t",  required_argument, 0, 'e'}, 
               {"delta_xy",  required_argument, 0, 'f'}, /*normalized value w.r.t. delta_r*/
               {"delta_z",  required_argument, 0, 'z'}, /*normalized value w.r.t. delta_t*/
               {"length_r",  required_argument, 0, 'g'}, /*length of detector across r axis in mm*/
               {"num_threads",  required_argument, 0, 'h'}, /*Number of threads*/
               {"length_t",  required_argument, 0, '3'}, /*length of detector across t axis in mm*/
               {"voxel_thresh",    required_argument, 0, 'j'}, /*Optimization stops once average magnitude of update exceeds this threshold*/
               {"iter",    required_argument, 0, 'k'}, /*Maximum number of iterations of optimization algorithm*/
               {"rotation_center",    required_argument, 0, 'l'}, /*Center of rotation of object*/
               {"alpha",    required_argument, 0, 'm'}, /*Over-relaxation value. Can be between 1 and 2*/
               {"time_reg",    no_argument, 0, 'n'}, /*If true enables time regularization*/
               {"sinobin",    required_argument, 0, 'o'}, /*If 1, reads projection and weight values from bin file. If 3, then reads in bright field counts and weight*/
               {"initICD",    required_argument, 0, 'p'}, /*parameter which specifies the type of interpolation used when initializing object*/
               {"writeTiff",    required_argument, 0, 'q'}, /*Enables writes to Tiff file. If 1 write tiff file. If 2 in addition write object and offset to bin to tiff every iteration*/
               {"NoNoise",    no_argument, 0, 'r'}, /*Does not add noise to projections if true*/
               {"Rtime0",  required_argument, 0, 's'}, /*Time of start of piecewise constant object*/
               {"Rtime_delta",  required_argument, 0, 't'}, /*Length of subinterval of the piecewise constant object along time*/
               {"Rtime_num",  required_argument, 0, 'u'}, /*Number of reconstructions*/
               {"num_projections",  required_argument, 0, 'v'}, /*number of projections used in reconstruction*/
               {"N_r",  required_argument, 0, 'w'}, /*N_r is the number of detector bins along r-axis. parallel to x-axis*/
               {"dont_reconstruct",    no_argument, 0, 'x'}, /*if used, program will not do reconstruction*/
               {"cost_thresh",  required_argument, 0, 'y'}, /*threshold on cost which decides convergence*/
               {"radius_obj",  required_argument, 0, 'i'}, /*Radius of object. Code will update voxels only in the circular region*/
               {"detector_slice_begin",  required_argument, 0, '1'}, /*First slice to be reconstructed*/
               {"detector_slice_num",  required_argument, 0, '2'}, /*Number of slices to be reconstructed*/
               {"phantom_N_xy",  required_argument, 0, '4'}, /*Resolution of phantom along x-y dimension*/
               {"phantom_N_z",  required_argument, 0, '5'}, /*Resolution of phantom along z dimension*/
               {"N_t",  required_argument, 0, '6'},/*Number of detector bins along t-axis (parallel to z axis)*/ 
               {"updateProjOffset",  required_argument, 0, '7'}, 
		/*If 1, just initializes the offset error. If 2, does not initialize but updates error offset. If 3, initializes and updates error offset*/
               {"no_NHICD",  no_argument, 0, '8'},/*if set, the code does not use NHICD*/ 
               {"WritePerIter",  no_argument, 0, '+'}, /*Writes the object and projection offset data to binary file every iteration of ICD*/
               {"only_Edge_Updates",  no_argument, 0, '-'}, /*Only updates the edges of the object*/
               {"zingerT",  required_argument, 0, '*'}, /*Threshold T of generalized Huber function above which measurement is considered anamolous*/
               {"zingerDel",  required_argument, 0, '^'}, /*Threshold \delta of generalized Huber function*/
               {"initMagUpMap",  no_argument, 0, '&'}, /*if set, initializes the magnitude update map from coarse resolution reconstruction*/
               {"readSino4mHDF",  no_argument, 0, '>'}, /*if set, reads the count values from HDF file and computes the projections*/ 
               {"do_VarEst",  no_argument, 0, '%'},/*if set, estimates the variance parameter*/ 
               {"Est_of_Var",  required_argument, 0, '('}, /*contains an initial estimate of variance parameter*/
               {0, 0, 0, 0}
         };

	printf("DebugMsg: Reading command line arguments ....\n");
	slope_HU = (HOUNSFIELD_WATER_MAP-HOUNSFIELD_AIR_MAP)/(WATER_MASS_ATT_COEFF*WATER_DENSITY-AIR_MASS_ATT_COEFF*AIR_DENSITY)/HFIELD_UNIT_CONV_CONST;
	c_HU = -slope_HU*(AIR_MASS_ATT_COEFF*AIR_DENSITY*HFIELD_UNIT_CONV_CONST);

	TomoInputsPtr->sinobin = 0;
	TomoInputsPtr->initICD = 0;
	TomoInputsPtr->Write2Tiff = 0;
	TomoInputsPtr->time_reg = 0;
	TomoInputsPtr->reconstruct = 1;
	TomoInputsPtr->No_Projection_Noise = 0;
	TomoInputsPtr->phantom_N_xy = 0;
	TomoInputsPtr->phantom_N_z = 0;
	TomoInputsPtr->cost_thresh = 100;
	TomoInputsPtr->updateProjOffset = 0;
	TomoInputsPtr->no_NHICD = 0;
	TomoInputsPtr->WritePerIter = 0;
	TomoInputsPtr->only_Edge_Updates = 0;
	TomoInputsPtr->initMagUpMap = 0;
	TomoInputsPtr->updateVar = 0;
	TomoInputsPtr->var_est = 1;
	TomoInputsPtr->readSino4mHDF = 0;
	while(1)
	{		
	   c = getopt_long (argc, argv, "a:b:c:d:e:f:z:g:h:3:j:k:l:m:no:p:q:rs:t:u:v:w:xy:i:1:2:4:5:6:7:8+-*:^:&>%(:", long_options, &option_index);
     
           /* Detect the end of the options. */
           if (c == -1) break;

	  switch (c) { 
		case  0 : printf("ERROR: Argument not recognized\n");		break;
		case 'a': ScannedObjectPtr->MRF_P=(Real_t)atof(optarg);		break;
		case 'b': ScannedObjectPtr->Sigma_S = (Real_t)(atof(optarg) - c_HU)/slope_HU;	break;
		/*case 'b': ScannedObjectPtr->Sigma_S = (Real_t)atof(optarg);	break;*/
		case 'c': ScannedObjectPtr->C_S = (Real_t)atof(optarg);	break;
		case 'd': ScannedObjectPtr->Sigma_T = (Real_t)(atof(optarg) - c_HU)/slope_HU;	break;
		/*case 'd': ScannedObjectPtr->Sigma_T = (Real_t)atof(optarg);	break; */
		case 'e': ScannedObjectPtr->C_T = (Real_t)atof(optarg);	break;
		case 'f': ScannedObjectPtr->mult_xy=(Real_t)atof(optarg);	break;
		case 'z': ScannedObjectPtr->mult_z=(Real_t)atof(optarg);	break;
		case 'g': SinogramPtr->Length_R = (Real_t)atof(optarg);	break;
		case 'h': TomoInputsPtr->num_threads = (int32_t)atoi(optarg);	break;
		case '3': SinogramPtr->Length_T = (Real_t)atof(optarg);	break;
		case 'j': TomoInputsPtr->StopThreshold=(Real_t)atof(optarg);	break;
		case 'k': TomoInputsPtr->NumIter=(int32_t)atoi(optarg);	break;
		case 'l': TomoInputsPtr->RotCenter=(Real_t)atof(optarg);	break;
		case 'm': TomoInputsPtr->alpha=(Real_t)atof(optarg);		break;
		case 'n': TomoInputsPtr->time_reg = 1;				break;
		case 'o': TomoInputsPtr->sinobin = (uint8_t)atoi(optarg);	break;
		case 'p': TomoInputsPtr->initICD = (uint8_t)atoi(optarg);	break;
		case 'q': TomoInputsPtr->Write2Tiff = (uint8_t)atoi(optarg);			break;
		case 'r': TomoInputsPtr->No_Projection_Noise = 1;		break;
		case 's': ScannedObjectPtr->Rtime0 = (Real_t)atof(optarg); break;
		case 't': ScannedObjectPtr->delta_Rtime = (Real_t)atof(optarg); break;	
		case 'u': ScannedObjectPtr->N_time = (int32_t)atoi(optarg); break;	
		case 'v': SinogramPtr->N_p = (int32_t)atoi(optarg); 		break;	
		case 'w': SinogramPtr->N_r = (int32_t)atoi(optarg); 		break;	
		case 'x': TomoInputsPtr->reconstruct = 0; 		break;	
		case 'y': TomoInputsPtr->cost_thresh = (Real_t)atof(optarg); 		break;	
		case 'i': TomoInputsPtr->radius_obj = (Real_t)atof(optarg); 		break;	
		case '1': SinogramPtr->slice_begin = (int32_t)atoi(optarg); 		break;	
		case '2': SinogramPtr->slice_num = (int32_t)atoi(optarg); 		break;	
		case '4': TomoInputsPtr->phantom_N_xy = (int32_t)atoi(optarg); 	break;	
		case '5': TomoInputsPtr->phantom_N_z = (int32_t)atoi(optarg); 		break;	
		case '6': SinogramPtr->total_t_slices = (int32_t)atoi(optarg); 		break;	
		case '7': TomoInputsPtr->updateProjOffset = (int32_t)atoi(optarg); 		break;	
		case '8': TomoInputsPtr->no_NHICD = 1; 		break;	
		case '+': TomoInputsPtr->WritePerIter = 1;	break;
		case '-': TomoInputsPtr->only_Edge_Updates = 1;	break;
		case '*': TomoInputsPtr->ErrorSinoThresh = (Real_t)atof(optarg);	break;
		case '^': TomoInputsPtr->ErrorSinoDelta = (Real_t)atof(optarg);	break;
		case '&': TomoInputsPtr->initMagUpMap = 1;	break;
		case '>': TomoInputsPtr->readSino4mHDF = 1;	break;
		case '%': TomoInputsPtr->updateVar = 1; break;
		case '(': TomoInputsPtr->var_est = (Real_t)atof(optarg); break;
		case '?': printf("ERROR: argsParser: Cannot recognize argument %s\n",optarg); break;
		}
	}

	sprintf(debug_filename ,"DEBUG_n%d_delta_xy_%d_delta_z_%d.log", TomoInputsPtr->node_rank, (int)ScannedObjectPtr->mult_xy, (int)ScannedObjectPtr->mult_z);
	if (TomoInputsPtr->readSino4mHDF == 1)
		sprintf(debug_filename ,"DEBUG_HDF_Read_n%d.log", TomoInputsPtr->node_rank);
	TomoInputsPtr->debug_file_ptr = fopen(debug_filename, "w" );
	printf ("Refer to %s for more information\n", debug_filename);
/*	TomoInputsPtr->debug_file_ptr = stdout;*/
	
	fprintf(TomoInputsPtr->debug_file_ptr, "argsParser: p = %.2f, sigma_s = %f, sigma_t = %f, c_s = %f, c_t = %f, mult_xy = %f, mult_z = %f, Length_R = %.2f, Length_T = %.2f, stop threshold = %.2f, number of iterations = %d, center of rotation = %.2f, alpha = %.2f, time regularization = %d, read sinogram from bin = %d, init ICD = %d, Write Tiff file = %d, Don't add noise = %d, Reconstruction start time = %f, Reconstruction time gap = %f, number of reconstructions = %d, N_p = %d, N_r = %d, reconstruct = %d, cost_thresh = %f, PHANTOM_FILENAME = %s, Slice Begin = %d, Slice Num = %d, Phantom X-Y Resolution = %d, Phantom Z Resolution = %d, N_t = %d, radius of object = %f, Update additive offset error = %d, no_NHICD = %d, Write Tiff and Bin every Iteration = %d, only Edge Updates = %d, Zinger threshold T = %f, Zinger Delta = %f, Read projection from HDF = %d, Update Variance = %d, Variance Estimate = %f\n",ScannedObjectPtr->MRF_P, ScannedObjectPtr->Sigma_S, ScannedObjectPtr->Sigma_T, ScannedObjectPtr->C_S, ScannedObjectPtr->C_T, ScannedObjectPtr->mult_xy, ScannedObjectPtr->mult_z, SinogramPtr->Length_R, SinogramPtr->Length_T, TomoInputsPtr->StopThreshold, TomoInputsPtr->NumIter, TomoInputsPtr->RotCenter, TomoInputsPtr->alpha, TomoInputsPtr->time_reg, TomoInputsPtr->sinobin, TomoInputsPtr->initICD, TomoInputsPtr->Write2Tiff, TomoInputsPtr->No_Projection_Noise, ScannedObjectPtr->Rtime0, ScannedObjectPtr->delta_Rtime, ScannedObjectPtr->N_time, SinogramPtr->N_p, SinogramPtr->N_r, TomoInputsPtr->reconstruct, TomoInputsPtr->cost_thresh, PHANTOM_FILENAME, SinogramPtr->slice_begin, SinogramPtr->slice_num, TomoInputsPtr->phantom_N_xy, TomoInputsPtr->phantom_N_z, SinogramPtr->total_t_slices, TomoInputsPtr->radius_obj, TomoInputsPtr->updateProjOffset, TomoInputsPtr->no_NHICD, TomoInputsPtr->WritePerIter, TomoInputsPtr->only_Edge_Updates, TomoInputsPtr->ErrorSinoThresh, TomoInputsPtr->ErrorSinoDelta, TomoInputsPtr->readSino4mHDF, TomoInputsPtr->updateVar, TomoInputsPtr->var_est);
	
#ifdef READ_PROJECTION_DATA_4M_HDF
	if (TomoInputsPtr->readSino4mHDF == 1)
		gen_projection_4m_HDF (SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
#endif

	if (TomoInputsPtr->sinobin == 3)
	{
		SinogramPtr->total_t_slices = SinogramPtr->total_t_slices/ScannedObjectPtr->mult_z;
		SinogramPtr->N_r = SinogramPtr->N_r/ScannedObjectPtr->mult_xy;
		ScannedObjectPtr->delta_xy = 1;
		ScannedObjectPtr->delta_z = 1;
		TomoInputsPtr->RotCenter /= ScannedObjectPtr->mult_xy;
		SinogramPtr->z_overlap_num = 1;
	}
	else
	{	
		ScannedObjectPtr->delta_xy = ScannedObjectPtr->mult_xy;
		ScannedObjectPtr->delta_z = ScannedObjectPtr->mult_z;
		SinogramPtr->z_overlap_num = ScannedObjectPtr->mult_z;
	}
  	
	if(argc-optind>0){
		fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: argsParser: Argument list has an error\n");
		/*exit(1);*/
	}
	
}
