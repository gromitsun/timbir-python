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


/*Function to compute the system matrix on the fly for 
  MBIR Recon*/

#ifndef XT_MATRIX_H
#define XT_MATRIX_H

#include "XT_Structures.h"

void findAMatrix (Sinogram* SinogramPtr, ScannedObject *ScannedObjectPtr, AMatrixCol ***AMatrix);
void AMatrix_free(AMatrixCol*** AMatrix, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr);
void calcAMatrixColumnforAngle (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, Real_t** DetectorResponse, AMatrixCol *Ai, int32_t row, int32_t col, int32_t proj_idx);
void compute_2DAMatrix_4m_1D(Real_t*** AMatrix2D, AMatrixCol* AMatrixPtr, AMatrixCol* VoxelLineResponse, int32_t* r_ax_start, int32_t* r_ax_count, int32_t* t_ax_start, int32_t* t_ax_count);
#endif /*#define XT_MATRIX_H*/
