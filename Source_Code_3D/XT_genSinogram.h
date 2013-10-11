#ifndef GEN_SINOGRAM_H
#define GEN_SINOGRAM_H
#include "XT_Structures.h"


int genSinogram_fromBin(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);
void genSinogramFromPhantom (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);
int genSinogram_subsampleCounts(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);

int genSinogram(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);
void ForwardProject(AMatrixCol** AMatrix, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, Real_t**** img1);

#endif /*#ifndef GEN_SINOGRAM_H*/
