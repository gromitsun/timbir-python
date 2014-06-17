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
#include <stdlib.h>
#include "XT_Structures.h"
#include <math.h>
#include "allocate.h"
#include "XT_IOMisc.h"

void filter_image (Real_t** in, Real_t** out, Real_t h[5][5], int32_t rows, int32_t cols, int32_t filter_len)
{
	int32_t i, j, k, l;
	
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
		{
			out[i][j] = 0;
			for (k = -(filter_len-1)/2; k <= (filter_len-1)/2; k++)
				for (l = -(filter_len-1)/2; l <= (filter_len-1)/2; l++)
				{
					if (i + k >= 0 && j + l >= 0 && i + k < rows && j + l < cols)
						out[i][j] += h[k+(filter_len-1)/2][l+(filter_len-1)/2]*in[i+k][j+l];
				}
		}
}

void image_dilation (Real_t** in, Real_t** out, int32_t rows, int32_t cols, int32_t length)
{
	int32_t i, j, k, l;

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
		{
			out[i][j] = in[i][j];
			for (k = -(length-1)/2; k <= (length-1)/2; k++)
			{
				for (l = -(length-1)/2; l <= (length-1)/2; l++)
				{
					if (i + k >= 0 && j + l >= 0 && i + k < rows && j + l < cols)					
					{
						if (in[i+k][j+l] > out[i][j])
						{
							out[i][j] = in[i+k][j+l];
						}
					}
				}
			}		
		}		
}

void edge_detect (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, uint8_t**** Mask)
{
	Real_t hx[5][5] = {{3, 0, -3, 0, 0},{10, 0, -10, 0, 0},{3, 0, -3, 0, 0}};
	Real_t hy[5][5] = {{3, 10, 3, 0, 0},{0, 0, 0, 0, 0},{-3, -10, -3, 0, 0}};
	Real_t hblur[5][5], x, y;
	Real_t **temp, **outx, **outy, var = 3.0, sum, thresh = 0.0002;
	int32_t i, j, k, rows, cols, n, m, l;
	char dilated_file[100] = "dilated";
	char edges_file[100] = "edges", filename[100];
	
	rows = ScannedObjectPtr->N_y;
	cols = ScannedObjectPtr->N_x;

	for (i = -2; i <= 2; i++)
		for (j = -2; j <= 2; j++)
		{
			hblur[i+2][j+2] = exp(-(i*i+j*j)/(2*var));
			sum += hblur[i+2][j+2];
		}	
	
	for (i = -2; i <= 2; i++)
		for (j = -2; j <= 2; j++)
			hblur[i+2][j+2] /= sum;	


	temp = (Real_t**)multialloc(sizeof(Real_t), 2, rows, cols);
	outx = (Real_t**)multialloc(sizeof(Real_t), 2, rows, cols);
	outy = (Real_t**)multialloc(sizeof(Real_t), 2, rows, cols);

	for (n = 0; n < ScannedObjectPtr->N_time; n++)
	for (m = 0; m < ScannedObjectPtr->N_z; m++)
	{
		for (k = 0; k < ScannedObjectPtr->N_y; k++)
			for (l = 0; l < ScannedObjectPtr->N_x; l++)
				outx[k][l] = ScannedObjectPtr->Object[n][k][l][m];
		filter_image (outx, temp, hblur, rows, cols, 5); 
		filter_image (temp, outx, hx, rows, cols, 3); 
		filter_image (temp, outy, hy, rows, cols, 3);
		for (i = 0; i < rows; i++)
			for (j = 0; j < cols; j++)
			{
				x = -ScannedObjectPtr->Length_X/2 + ((Real_t)j + 0.5)*ScannedObjectPtr->delta_xy;
				y = -ScannedObjectPtr->Length_Y/2 + ((Real_t)i + 0.5)*ScannedObjectPtr->delta_xy;
				temp[i][j] = 0; outx[i][j] = 0;
				if (x*x + y*y < TomoInputsPtr->radius_obj*TomoInputsPtr->radius_obj)
				{
					temp[i][j] = sqrt(outx[i][j]*outx[i][j]+outy[i][j]*outy[i][j]);
					if (temp[i][j] > thresh)
						outx[i][j] = 1;
				}
			}

		sprintf(filename,"%s_time_%d_z_%d.tif",edges_file, n, m);	
		Write2Tiff(filename, rows, cols, temp, 0, TomoInputsPtr->debug_file_ptr);
		image_dilation (outx, outy, rows, cols, 5);
		sprintf(filename,"%s_time_%d_z_%d.tif", dilated_file, n, m);	
		Write2Tiff(filename, rows, cols, outy, 0, TomoInputsPtr->debug_file_ptr);

		for (k = 0; k < ScannedObjectPtr->N_y; k++)
			for (l = 0; l < ScannedObjectPtr->N_x; l++)
				Mask[n][k][l][m] = (uint8_t)outy[k][l];
	}
	
	multifree(outx,2);
	multifree(outy,2);
	multifree(temp,2);
}

void select_edge_voxels (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t**** MagUpdateMap, uint8_t**** Mask)
{
	int32_t i, j, k, l, n, flag;

	edge_detect (ScannedObjectPtr, TomoInputsPtr, Mask);

	for (i = 0; i < ScannedObjectPtr->N_time; i++)
		for (k = 0; k < ScannedObjectPtr->N_y; k++)
			for (l = 0; l < ScannedObjectPtr->N_x; l++)
				for (j = 0; j < TomoInputsPtr->num_z_blocks; j++)
				{
					flag = 0;
					for (n = j*(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks); n < (j+1)*(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks); n++)
						if (Mask[i][k][l][n] == 1)
							flag = 1;
					if (flag == 0)
						MagUpdateMap[i][j][k][l] = 0;
				}
}
