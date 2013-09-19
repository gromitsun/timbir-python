#ifndef XT_PROFILE_H
#define XT_PROFILE_H
#include "XT_Structures.h"
void calculateVoxelProfile(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t** VoxProfile);
void initializeBeamProfile(ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t *BeamProfile);
void storeVoxelLineResponse(Real_t* H_t,  AMatrixCol* VoxelLineResponse, ScannedObject* ScannedObjectPtr, Sinogram* SinogramPtr);
void DetectorResponseProfile (Real_t** H_r, Real_t* H_t, Sinogram* SinogramPtr, ScannedObject *ScannedObjectPtr, TomoInputs* TomoInputsPtr);
#endif /*#ifndef XT_PROFILE_H*/
