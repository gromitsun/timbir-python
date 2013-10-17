#include "grid_readwrite.h"
using namespace std;

/* Inputs: */
void GridReadWrite::WriteFloatv2(double* data, string output_name, long image_size, float pixel_size)
{
	float x_resolution = 10.0/pixel_size, y_resolution = 10.0/pixel_size;
	std::cout<<"Entering tiff writer "<<std::endl;

	TIFF* result_tif = TIFFOpen(output_name.c_str(), "w");
	if(!result_tif){
		this->error_flag = 1;
		return;
	}
	TIFFSetField(result_tif, TIFFTAG_IMAGEWIDTH, image_size);  
	TIFFSetField(result_tif, TIFFTAG_IMAGELENGTH, image_size);    
	TIFFSetField(result_tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(result_tif, TIFFTAG_SAMPLEFORMAT, 3); //floating point = 3
	TIFFSetField(result_tif, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(result_tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(result_tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(result_tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	if(pixel_size < 1.0){
		TIFFSetField(result_tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER);
		TIFFSetField(result_tif, TIFFTAG_XRESOLUTION, x_resolution);
		TIFFSetField(result_tif, TIFFTAG_YRESOLUTION, y_resolution);
	}
	float* buf = NULL;
	tsize_t linebytes =  image_size*sizeof(float);
	buf =(float *)_TIFFmalloc(linebytes);
	
	for(uint32 row = 0; row < image_size; row++){
		for(int col = 0; col < image_size; col++){
		  buf[col] = float(data[col+row*image_size])/pixel_size;
		}
		if (TIFFWriteScanline(result_tif, buf, row, 0) < 0){
			cout << "Write Error!\n";
			this->error_flag = 1;
			break;
		}
	}
	TIFFClose(result_tif);
	if(buf){
		_TIFFfree(buf);
	}
	
}



//Constructor

GridReadWrite::GridReadWrite(void)
{
}

//Destructor

GridReadWrite::~GridReadWrite(void)
{
}






