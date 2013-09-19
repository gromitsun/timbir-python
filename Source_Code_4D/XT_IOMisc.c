#include <stdio.h>
#include "XT_Constants.h"
#include "tiff.h"
#include "allocate.h"
#include "XT_Structures.h"

Real_t convert_HU2um (Real_t val)
{
       Real_t slope_HU=(HOUNSFIELD_WATER_MAP-HOUNSFIELD_AIR_MAP)/(WATER_MASS_ATT_COEFF*WATER_DENSITY-AIR_MASS_ATT_COEFF*AIR_DENSITY)/HFIELD_UNIT_CONV_CONST;
       Real_t c_HU=-slope_HU*(AIR_MASS_ATT_COEFF*AIR_DENSITY*HFIELD_UNIT_CONV_CONST);

       return((val-c_HU)/slope_HU);
}

void Append2Bin (char *filename, int dim1, int dim2, int dim3, int dim4, Real_t* img, FILE *debug_file_ptr)
{
  FILE *fp;
  char file[100];
  sprintf(file,"%s.bin", filename);
  fp=fopen (file,"ab");
  if (fp == NULL){fprintf(debug_file_ptr, "ERROR: Write2Bin: Cannot open file %s\n", file);exit(1);}
  fwrite(img,sizeof(Real_t),dim1*dim2*dim3*dim4,fp);	
  fclose(fp);
}

void Write2Bin (char *filename, int dim1, int dim2, int dim3, int dim4, Real_t* img, FILE *debug_file_ptr)
{
  FILE *fp;
  char file[100];
  sprintf(file,"%s.bin", filename);
  fp=fopen (file,"wb");
  if (fp == NULL){fprintf(debug_file_ptr, "ERROR: Write2Bin: Cannot open file %s\n", file);exit(1);}
  fwrite(img,sizeof(Real_t),dim1*dim2*dim3*dim4,fp);	
  fclose(fp);
}

void Read4mBin (char *filename, int dim1, int dim2, int dim3, int dim4, Real_t* img, FILE *debug_file_ptr)
{
  	char file[100];
  	sprintf(file,"%s.bin", filename);
	FILE* fp;
	int32_t size, result;

	size = dim1*dim2*dim3*dim4;
	fp = fopen (file, "rb" );
	if (fp==NULL) {fprintf(debug_file_ptr, "ERROR: error in reading file %s\n", file); exit (1);}	
	result = fread (img, sizeof(Real_t), size, fp);
  	if (result != size) 
	{fprintf(debug_file_ptr, "ERROR: Number of elements read does not match required, required = %d, read = %d\n", size, result);}
	fclose(fp);
}

void Write2Tiff(char* filename, int height, int width, Real_t** img, int hounsfield_flag, FILE *debug_file_ptr)
{
	struct TIFF_img out_img;
	int i,j;
	Real_t pixel;
	Real_t maxpix,minpix,avgpix,scale;
	FILE *fp;
	char file[100];
	Real_t slope,c;

	if(hounsfield_flag==1){
		slope=(HOUNSFIELD_WATER_MAP-HOUNSFIELD_AIR_MAP)/(WATER_MASS_ATT_COEFF*WATER_DENSITY-AIR_MASS_ATT_COEFF*AIR_DENSITY)/HFIELD_UNIT_CONV_CONST;
		c=-slope*(AIR_MASS_ATT_COEFF*AIR_DENSITY*HFIELD_UNIT_CONV_CONST);
		for ( i = 0; i < height; i++ )
	 	 for ( j = 0; j < width; j++ ) 
			img[i][j]=slope*img[i][j]+c;
#ifdef DEBUG_EN
		fprintf(debug_file_ptr, "Write2Tiff: file=%s, Hounsfeld slope = %f, c = %f\n",filename, slope,c);
#endif
	}	
	get_TIFF ( &out_img, height, width, 'g' );
	avgpix=0;
		
	maxpix=img[0][0];minpix=maxpix;
	 for ( i = 0; i < height; i++ )
 	 for ( j = 0; j < width; j++ ) {
		maxpix=(maxpix<img[i][j])?img[i][j]:maxpix;
		minpix=(minpix>img[i][j])?img[i][j]:minpix;
		avgpix+=img[i][j];	
	}
	avgpix/=(height*width);
	fprintf(debug_file_ptr, "Write2Tiff:file=%s,maxpix=%f,minpix=%f,avgpix=%f,height=%d,width=%d\n",filename,maxpix,minpix,avgpix,height,width);

	if(hounsfield_flag==1){
		maxpix=HOUNSFIELD_MAX;
		minpix=HOUNSFIELD_MIN;
	}
	scale=255/(maxpix-minpix);
	 for ( i = 0; i < height; i++ )
 	 for ( j = 0; j < width; j++ ) {
    		pixel = (int32_t)((img[i][j]-minpix)*scale);
    		if(pixel>255) {
      			out_img.mono[i][j] = 255;
    		}
    		else {
      			if(pixel<0) out_img.mono[i][j] = 0;
      			else out_img.mono[i][j] = pixel;
    		}
  	}
	sprintf(file,"%s.tif", filename);	
	if ( ( fp = fopen ( file, "wb" ) ) == NULL ) {
    		fprintf ( stderr, "cannot open file %s\n",filename);
	        exit ( 1 );
  	}

	/* write image */
	if ( write_TIFF ( fp, &out_img ) ) {
    		fprintf ( stderr, "error writing TIFF file %s\n", filename );
		exit ( 1 );
  	}

  /* close image file */
  fclose ( fp );
}

void WriteUint82Tiff(char* filename, int height, int width, bool** imgin, int hounsfield_flag, FILE *debug_file_ptr)
{
	Real_t** img;
	int32_t i, j;

	img = (Real_t**)multialloc(sizeof(Real_t), 2, height, width);
	for (i = 0; i < height; i++)
	for (j = 0; j < height; j++)
		img[i][j] = (Real_t)imgin[i][j];
	Write2Tiff(filename, height, width, img, hounsfield_flag, debug_file_ptr);
	multifree(img,2);
}

void WriteInt32Tiff(char* filename, int height, int width, int32_t** imgin, int hounsfield_flag, FILE *debug_file_ptr)
{
	Real_t** img;
	int32_t i, j;

	img = (Real_t**)multialloc(sizeof(Real_t), 2, height, width);
	for (i = 0; i < height; i++)
	for (j = 0; j < height; j++)
		img[i][j] = (Real_t)imgin[i][j];
	Write2Tiff(filename, height, width, img, hounsfield_flag, debug_file_ptr);
	multifree(img,2);
}


void WriteMultiDimArray2Tiff (char *filename, int dim[4], int dim2loop_1, int dim2loop_2, int dim2write_1, int dim2write_2, Real_t* img, int hounsfield_flag, FILE* debug_file_ptr)
{
	char file[100];
	Real_t** img_temp;
	int i,j,temp,k,l;

	if (dim2loop_1 > dim2loop_2){
		temp = dim2loop_2;
		dim2loop_2 = dim2loop_1;
		dim2loop_1 = temp;
		printf("WARNING: WriteMultiDimArray2Tiff: Swapping variables dim2write_1 and dim2write_2\n");
	}
 		
	fprintf(debug_file_ptr, "WriteMultiDimArray2Tiff: Writing to tiff file %s*.tif \n", filename);
    img_temp = (Real_t**) multialloc(sizeof(Real_t), 2, dim[dim2write_1], dim[dim2write_2]);
	for (i=0; i<dim[dim2loop_1]; i++){
	for (j=0; j<dim[dim2loop_2]; j++){
	for (k=0; k<dim[dim2write_1]; k++){
	for (l=0; l<dim[dim2write_2]; l++){
		if (dim2write_1 == 2 && dim2write_2 == 3)
			img_temp[k][l] = img[((i*dim[1]+j)*dim[2]+k)*dim[3]+l];
		else if (dim2write_1 == 1 && dim2write_2 == 3)
			img_temp[k][l] = img[((i*dim[1]+k)*dim[2]+j)*dim[3]+l]; 	
		else if (dim2write_1 == 0 && dim2write_2 == 3)
			img_temp[k][l] = img[((k*dim[1]+i)*dim[2]+j)*dim[3]+l]; 
		else if (dim2write_1 == 1 && dim2write_2 == 2)
			img_temp[k][l] = img[((i*dim[1]+k)*dim[2]+l)*dim[3]+j]; 	
		else if (dim2write_1 == 0 && dim2write_2 == 2)
			img_temp[k][l] = img[((k*dim[1]+i)*dim[2]+l)*dim[3]+j]; 	
		else if (dim2write_1 == 0 && dim2write_2 == 1)
			img_temp[k][l] = img[((k*dim[1]+l)*dim[2]+i)*dim[3]+j]; 
		else
			printf("ERROR: WriteMultiDimArray2Tiff: Dimensions not recognized\n");
	}}
		sprintf(file, "%s_d1_%d_d2_%d",filename,i,j);
		Write2Tiff(file, dim[dim2write_1], dim[dim2write_2], img_temp, hounsfield_flag, debug_file_ptr);
	}
	}

	multifree (img_temp,2);
}

void write_ObjectProjOff2TiffBinPerIter (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	int32_t i, dim[4];
	char object_file[100];
	char projOffset_file[100] = PROJ_OFFSET_FILENAME;

	sprintf(projOffset_file, "%s_n%d", projOffset_file, TomoInputsPtr->node_rank);
	for (i = 0; i < ScannedObjectPtr->N_time; i++)
	{
		sprintf (object_file, "%s_n%d_time_%d", OBJECT_FILENAME, TomoInputsPtr->node_rank, i);
		Write2Bin (object_file, 1, ScannedObjectPtr->N_z, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x, &(ScannedObjectPtr->Object[i][1][0][0]), TomoInputsPtr->debug_file_ptr);
		dim[0] = 1; dim[1] = ScannedObjectPtr->N_z; dim[2] = ScannedObjectPtr->N_y; dim[3] = ScannedObjectPtr->N_x;
		if (TomoInputsPtr->Write2Tiff == 1)
			WriteMultiDimArray2Tiff (object_file, dim, 0, 1, 2, 3, &(ScannedObjectPtr->Object[i][1][0][0]), 1, TomoInputsPtr->debug_file_ptr);
	}
	dim[0] = 1; dim[1] = 1; dim[2] = SinogramPtr->N_r; dim[3] = SinogramPtr->N_t;
	Write2Bin (projOffset_file, 1, 1, SinogramPtr->N_r, SinogramPtr->N_t, &(SinogramPtr->ProjOffset[0][0]), TomoInputsPtr->debug_file_ptr);
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (projOffset_file, dim, 0, 1, 2, 3, &(SinogramPtr->ProjOffset[0][0]), 0, TomoInputsPtr->debug_file_ptr);

}

