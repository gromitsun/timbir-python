#ifndef XT_IOMISC_H
#define XT_IOMISC_H

Real_t convert_HU2um(Real_t val);
void Write2Bin (char *filename, int dim1, int dim2, int dim3, int dim4, Real_t* img, FILE *debug_file_ptr);
void Append2Bin (char *filename, int dim1, int dim2, int dim3, int dim4, Real_t* img, FILE *debug_file_ptr);
void WriteMultiDimArray2Tiff (char *filename, int dim[4], int dim2loop_1, int dim2loop_2, int dim2write_1, int dim2write_2, Real_t* img, int hounsfield_flag, FILE* debug_file_ptr);
void Read4mBin (char *filename, int dim1, int dim2, int dim3, int dim4, Real_t* img, FILE *debug_file_ptr);
void write_ObjectProjOff2TiffBinPerIter (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);
void WriteUint82Tiff(char* filename, int height, int width, uint8_t** imgin, int hounsfield_flag, FILE *debug_file_ptr);
void WriteInt32Tiff(char* filename, int height, int width, int32_t** imgin, int hounsfield_flag, FILE *debug_file_ptr);
void Write2Tiff(char* filename, int height, int width, Real_t** img, int hounsfield_flag, FILE *debug_file_ptr);

#endif
