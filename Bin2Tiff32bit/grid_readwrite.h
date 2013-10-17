#ifndef GRID_READWRITE_H
#define GRID_READWRITE_H


#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include "tiff.h"
#include "tiffio.h"

class GridReadWrite
{
private:
	float min, max, slope, offset;
	int num_zeros_each_side;
	
public:
	bool error_flag;
	uint32 width, padded_width, num_projections;

	GridReadWrite(void);
	~GridReadWrite(void);
	
        void WriteFloatv2(double* data, std::string output_name, long image_size, float pixel_size);
	

};
#endif
