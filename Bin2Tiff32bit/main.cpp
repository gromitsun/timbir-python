/*Code to take a bin file containing double precision floats and 
convert to to a collection of 32 bit tiffs 
Based on code from Justin Blair's grid rec wrappers
*/

/*  Compiling on carver
 *  module load tiff
 *  setenv GRIDREC "-L${TIFF_DIR}/lib -I${TIFF_DIR}/include -ltiff"
 *  g++ -o Bin2Tiff main.cpp grid_readwrite.cpp $GRIDREC 
 * Running:
 * ./Bin2Tiff <input file - FULL path> <Output folder name - Full path and the folder SHOULD exist> <Output Base File name> <image size in pixels - assume length and width are same> <first z slice in bin file> <number of slices> <pixel size in micro meter>
 Example: 
./Bin2Tiff.x $GSCRATCH/Results/Edison/MBIR_sigs_666666.666667_sigt_1_r_1_K_1_N_theta_1024_N_p_256_zinger_50.0_0.1/object_n0_time_0.bin $GSCRATCH/TestTiffWriter/ Parikh_mbir_recon_ 2560 1000 32 .65
 */

#define INTERP 1

#include <algorithm>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include "tiff.h"
#include "tiffio.h"
#include "grid_readwrite.h"

using namespace std;

//Concatenates a string and a number 
string GetName(string name_base, int img_num){
	stringstream stream;
	string number;
	stream << name_base << img_num << ".tif";
	return stream.str();
}

int main(int argc, char** argv){
        
		TIFFSetWarningHandler(NULL);
		TIFFSetErrorHandler(NULL);

		int first_img_num;
		int num_slices;
		float pixel_size;
		long image_size;
		float **reconstructed_1_rows=0; 		
                float *test_data;
		GridReadWrite* fileIO;
                string input_file;
		string output_base, output_name;
		string output_path;

		/* Parse the input paramters from the command line. */
               
                input_file.assign(argv[1]);
		output_path.assign(argv[2]);
		output_base.assign(argv[3]);
                image_size = atof(argv[4]);
                first_img_num = atoi(argv[5]);
                num_slices = atoi(argv[6]);
                pixel_size = atof(argv[7]);		

                /*Print parameters*/
		std::cout<<" Input file :"<<input_file<<std::endl;
                std::cout<<" Output path :"<<output_path<<std::endl;
                std::cout<<" Output base name :"<<output_base<<std::endl;
                std::cout<<"First image="<<first_img_num<<std::endl;
                std::cout<<"Num slices="<<num_slices<<std::endl;
		std::cout<<"pix size="<<pixel_size<<std::endl;

		/*Read Input bin file*/
                FILE *fp;

		//Convert the string to character array
		char* inpFile = new char[input_file.size() + 1];
		std::copy(input_file.begin(), input_file.end(), inpFile);
		inpFile[input_file.size()] = '\0'; // don't forget the terminating 0

                fp=fopen(inpFile,"rb");         
                if(fp == NULL){
		  fprintf(stderr, "Input file opening failed. Verify the full path is correct\n");
		  return 1;
		}

		//Generate the full output path to write tiffs 
		if(output_path.substr(output_path.length() - 1) != "/" && output_path.substr(output_path.length() - 1) != "\\")
                {
			output_path.append("/");
		}
		
		//Initialize an object of the readWrite class 
		fileIO = new GridReadWrite;
	
		std::cout<<"Allocating memeory for a individual tiff\n";
   	   
	
		//Allocate temporary buffers to store each image of the reconstruction
		
                double *tmp_buffer = (double*)calloc(image_size*image_size,sizeof(double)); 
                float *buffer = (float*)calloc(image_size*image_size,sizeof(float));
               
		//std::cout<<"Value test float"<<buffer[0]<<std::endl;
                //std::cout<<"Value test double"<<tmp_buffer[0]<<std::endl;

		/* Loop through the slices writing 1 at a time*/	
		for(int img = first_img_num; img < first_img_num+num_slices; img++) {
              	        std::cout<<"Starting to write the image slice index : "<<img<<std::endl;
			output_name = GetName(output_base, img);
  			printf("\nWriting %s." ,output_name.c_str());

			//Read a single image
                         fread(tmp_buffer,sizeof(double),image_size*image_size,fp);
               
		       //Casting from double in the bin file to float
                       for (int i =0; i < image_size*image_size; i++) 
		           buffer[i] = (float)tmp_buffer[i];

                       fileIO->WriteFloatv2(tmp_buffer, output_path+output_name, image_size, pixel_size);
		}
		printf("\nComplete!\n");
		return 0;
}


