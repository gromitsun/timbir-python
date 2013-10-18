/* ============================================================================
 * Copyright (c) 2013 K. Aditya Mohan (Purdue University)
 * Copyright (c) 2013 Singanallur Venkatakrishnan (Purdue University)
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
 * Neither the name of K. Aditya Mohan, Singanallur Venkatakrishnan, Purdue
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
#include "XT_AMatrix.h"
#include "XT_Profile.h"
#include "XT_Constants.h"
#include "allocate.h"
#include <math.h>
#include "XT_IOMisc.h"

/*'calcAMatrixColumnforAngle' computes the A matrix column for any voxel and angle of choice. This function is called repeatedy during ICD optimization. The indexing in angle is required since this code implements sparse angle sampling.*/
void calcAMatrixColumnforAngle (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, Real_t** DetectorResponse, AMatrixCol *Ai, int32_t row, int32_t col, int32_t proj_idx)
{
  int32_t j;
  Real_t x,y;
  Real_t r;
  Real_t rmax,rmin;
  Real_t R_Center,delta_r;
  Real_t w1,w2,f1;
  int32_t index_min,index_max,index_delta_r;/*stores the detector index in which the profile lies*/
  int32_t count = 0;

  x = ScannedObjectPtr->x0 + ((Real_t)col+0.5)*ScannedObjectPtr->delta_xy;/*0.5 is for center of voxel. x_0 is the left corner*/
  y = ScannedObjectPtr->y0 + ((Real_t)row+0.5)*ScannedObjectPtr->delta_xy;/*0.5 is for center of voxel. y_0 is the left corner*/
/*'r' is the center of the voxel as projected on the detector. 'rmin' and 'rmax' gives the leftmost and rightmost distance at which the voxel of choice may have non-zero projection. We include some overhead*/
    r = x*SinogramPtr->cosine[proj_idx] - y*SinogramPtr->sine[proj_idx];

    rmin = r - ScannedObjectPtr->delta_xy;
    rmax = r + ScannedObjectPtr->delta_xy;

    if(rmax < SinogramPtr->R0 || rmin > SinogramPtr->RMax){
	Ai->count = 0;
        return;
    }

    index_min = (int32_t)(floor(((rmin - SinogramPtr->R0)/SinogramPtr->delta_r)));
    index_max = (int32_t)(floor((rmax - SinogramPtr->R0)/SinogramPtr->delta_r));
/*index_min and index_max is the quantized version of rmin and rmax*/

    if(index_max >= SinogramPtr->N_r)
      index_max = SinogramPtr->N_r - 1;

    if(index_min < 0)
      index_min = 0;
	
    for(j = index_min;j <= index_max; j++)/*Check*/
    {

      /*Accounting for different sensitivity across the detector*/
      R_Center = (SinogramPtr->R0 + (((Real_t)j) + 0.5) *(SinogramPtr->delta_r));/*the 0.5 is to get to the center of the detector*/

      /*Find the difference between the center of detector and center of projection and compute the Index to look up into*/
      delta_r = fabs(r - R_Center);
      index_delta_r = (int32_t)(floor((delta_r/SinogramPtr->OffsetR)));

      if (index_delta_r >= 0 && index_delta_r < DETECTOR_RESPONSE_BINS)
      {
            w1 = delta_r - index_delta_r*SinogramPtr->OffsetR;
	    w2 = SinogramPtr->OffsetR - w1;

            f1 = (w2*DetectorResponse[proj_idx][index_delta_r] + w1*DetectorResponse[proj_idx][index_delta_r+1])/SinogramPtr->OffsetR;

            if(f1 > 0)
            {
              Ai->values[count] = f1;
              Ai->index[count] = j;
              count++;
            }
       }
    }
    Ai->count = count;
}


