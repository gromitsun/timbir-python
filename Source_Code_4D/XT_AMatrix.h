/*Function to compute the system matrix on the fly for 
  MBIR Recon*/

#ifndef XT_MATRIX_H
#define XT_MATRIX_H

#include "XT_Structures.h"

void findAMatrix (Sinogram* SinogramPtr, ScannedObject *ScannedObjectPtr, AMatrixCol ***AMatrix);
void AMatrix_free(AMatrixCol*** AMatrix, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr);
void calcAMatrixColumnforAngle (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, Real_t** DetectorResponse, AMatrixCol *Ai, int32_t row, int32_t col, int32_t proj_idx);
#endif /*#define XT_MATRIX_H*/
