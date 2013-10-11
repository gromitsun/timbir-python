#ifndef XT_INIT_H
#define XT_INIT_H

void initStructures (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);
void argsParser (int argc, char **argv, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);

#endif /*#ifndef XT_INIT_H*/
