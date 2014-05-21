#ifndef XT_HDF_IO_H
#define XT_HDF_IO_H

#ifdef READ_PROJECTION_DATA_4M_HDF
#include "hdf5.h"
void gen_projection_4m_HDF (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);
void read_4m_HDF (TomoInputs* TomoInputsPtr, char filename[], Real_t* object, hsize_t start[], hsize_t num[]);
void write_2_existHDF (TomoInputs* TomoInputsPtr, char filename[], Real_t* object, hsize_t start[], hsize_t num[]);
void write_2_newHDF (TomoInputs* TomoInputsPtr, char filename[], Real_t* object, hsize_t start[], hsize_t num[], int32_t rank, hsize_t dims[]);
#endif

#endif
