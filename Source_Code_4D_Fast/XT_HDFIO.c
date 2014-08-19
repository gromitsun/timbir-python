#include <stdio.h>
#include "XT_IOMisc.h"
#include <math.h>
#include "XT_Constants.h"
#include "XT_Structures.h"
#include <stdlib.h>
#include "allocate.h"


#ifdef READ_PROJECTION_DATA_4M_HDF
#include "hdf5.h"
/*Instead of reading projection and weight data from binary files, use the function below to read those directly from
HDF files.*/
void gen_projection_4m_HDF (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	hid_t data_file_id, white_file_id, dark_file_id, white_dataset, dark_dataset, proj_dataset;
	hid_t white_dataspace, dark_dataspace, proj_dataspace, white_memspace, dark_memspace, proj_memspace;
	hsize_t dims_white[3], dims_dark[3], dims_proj[3], white_offset[3], dark_offset[3], proj_offset[3], white_count[3], dark_count[3], proj_count[3], mem_offset[3]; 
   	herr_t status;
	int32_t i, j, k, m, n, white_rank, dark_rank, proj_rank, extras_r, true_length_r, ratio_r, ratio_t, total_t_slices, dim[4];	
	uint16_t ***white_img, ***dark_img, ***proj_img;
	Real_t temp, ***Projection, ***Weight, **white_2D_img, **dark_2D_img, **wd_dwnsmpl_img;
	char wd_filename[100] = "bright";
	char weight_filename[100] = "weight";
	char proj_filename[100] = "projection";

	char data_hdf_filename[] = DATA_HDF_FILENAME;
	char white_hdf_filename[] = WHITE_HDF_FILENAME;
	char dark_hdf_filename[] = DARK_HDF_FILENAME;

	sprintf(wd_filename, "%s_n%d", wd_filename, TomoInputsPtr->node_rank);
	sprintf(weight_filename, "%s_n%d", weight_filename, TomoInputsPtr->node_rank);
	sprintf(proj_filename, "%s_n%d", proj_filename, TomoInputsPtr->node_rank);

	/*HDF file pointers*/
	data_file_id = H5Fopen(data_hdf_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	white_file_id = H5Fopen(white_hdf_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	dark_file_id = H5Fopen(dark_hdf_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	/*dataset pointers*/
	white_dataset = H5Dopen(white_file_id, "/exchange/data_white", H5P_DEFAULT);
	dark_dataset = H5Dopen(dark_file_id, "/exchange/data_dark", H5P_DEFAULT);
	proj_dataset = H5Dopen(data_file_id, "/exchange/data", H5P_DEFAULT);

	white_dataspace = H5Dget_space (white_dataset);    /* dataspace handle */
	dark_dataspace = H5Dget_space (dark_dataset);    /* dataspace handle */
	proj_dataspace = H5Dget_space (proj_dataset);    /* dataspace handle */
	/*Gives the number of dimensions in a dataset*/
    	white_rank = H5Sget_simple_extent_ndims (white_dataspace);
    	dark_rank = H5Sget_simple_extent_ndims (dark_dataspace);
    	proj_rank = H5Sget_simple_extent_ndims (proj_dataspace);
	
	if (white_rank != 3 || dark_rank != 3 || proj_rank != 3)
	{
		printf("ERROR: gen_projection_4m_HDF: The rank of one of the datasets in one of the HDF files is not 3\n");
		exit(1);
	}

	/*finds the dimension of the dataset and stores them in dims_wd and dims_proj*/	
	status = H5Sget_simple_extent_dims (white_dataspace, dims_white, NULL);
	status = H5Sget_simple_extent_dims (dark_dataspace, dims_dark, NULL);
	status = H5Sget_simple_extent_dims (proj_dataspace, dims_proj, NULL);
	fprintf(TomoInputsPtr->debug_file_ptr, "gen_projection_4m_HDF: size of white (/exchange/data_white) dataset is %zux%zux%zu\n", dims_white[0], dims_white[1], dims_white[2]);
	fprintf(TomoInputsPtr->debug_file_ptr, "gen_projection_4m_HDF: size of dark (/exchange/data_dark) dataset is %zux%zux%zu\n", dims_dark[0], dims_dark[1], dims_dark[2]);
	fprintf(TomoInputsPtr->debug_file_ptr, "gen_projection_4m_HDF: size of count (/exchange/data) dataset is %zux%zux%zu\n", dims_proj[0], dims_proj[1], dims_proj[2]);

	if (dims_white[2] != dims_proj[2] || dims_white[2] < SinogramPtr->N_r)
	{
		printf("ERROR: gen_projection_4m_HDF: dims_white[2] = %zu, dims_proj[2] = %zu, N_r = %d\n", dims_white[2], dims_proj[2], SinogramPtr->N_r);
		exit(1);
	}
	
	if (dims_white[1] != dims_proj[1] || SinogramPtr->slice_num % SinogramPtr->total_t_slices != 0)
	{
		printf("ERROR: gen_projection_4m_HDF: dims_white[1] = %zu, dims_proj[1] = %zu, SinogramPtr->slice_num = %d, N_t = %d\n", dims_white[1], dims_proj[1], SinogramPtr->slice_num, SinogramPtr->total_t_slices);
		exit(1);
	}
	
	if (dims_dark[2] != dims_proj[2] || dims_dark[2] < SinogramPtr->N_r)
	{
		printf("ERROR: gen_projection_4m_HDF: dims_dark[2] = %zu, dims_proj[2] = %zu, N_r = %d\n", dims_dark[2], dims_proj[2], SinogramPtr->N_r);
		exit(1);
	}
	
	if (dims_dark[1] != dims_proj[1] || SinogramPtr->slice_num % SinogramPtr->total_t_slices != 0)
	{
		printf("ERROR: gen_projection_4m_HDF: dims_dark[1] = %zu, dims_proj[1] = %zu, SinogramPtr->slice_num = %d, N_t = %d\n", dims_dark[1], dims_proj[1], SinogramPtr->slice_num, SinogramPtr->total_t_slices);
		exit(1);
	}

        extras_r = dims_proj[2] % SinogramPtr->N_r;
        true_length_r = dims_proj[2] - extras_r;
	SinogramPtr->Length_R = SinogramPtr->Length_R*true_length_r/dims_proj[2];	
/*	TomoInputsPtr->radius_obj = TomoInputsPtr->radius_obj*true_length_r/dims_proj[2];*/	
 
	white_offset[0] = 1;
	white_offset[1] = SinogramPtr->slice_begin + TomoInputsPtr->node_rank*SinogramPtr->slice_num/TomoInputsPtr->node_num;
    	white_offset[2] = extras_r/2;
	
	dark_offset[0] = 1;
	dark_offset[1] = SinogramPtr->slice_begin + TomoInputsPtr->node_rank*SinogramPtr->slice_num/TomoInputsPtr->node_num;
    	dark_offset[2] = extras_r/2;
	
	proj_offset[0] = PROJECTION_HDF_START;
	proj_offset[1] = SinogramPtr->slice_begin + TomoInputsPtr->node_rank*SinogramPtr->slice_num/TomoInputsPtr->node_num;
    	proj_offset[2] = extras_r/2;

	white_count[0] = dims_white[0] - 2;
	white_count[1] = SinogramPtr->slice_num/TomoInputsPtr->node_num;
	white_count[2] = true_length_r;
	
	dark_count[0] = dims_dark[0] - 2;
	dark_count[1] = SinogramPtr->slice_num/TomoInputsPtr->node_num;
	dark_count[2] = true_length_r;

	proj_count[0] = SinogramPtr->N_p;
	proj_count[1] = SinogramPtr->slice_num/TomoInputsPtr->node_num;
	proj_count[2] = true_length_r;

	ratio_r = true_length_r/SinogramPtr->N_r;
	ratio_t = SinogramPtr->slice_num/SinogramPtr->total_t_slices;
	total_t_slices = SinogramPtr->total_t_slices/TomoInputsPtr->node_num;

	fprintf(TomoInputsPtr->debug_file_ptr, "gen_projection_4m_HDF: true_length_r = %d, total_t_slices = %d, extras_r = %d, node_rank = %d\n", true_length_r, total_t_slices, extras_r, TomoInputsPtr->node_rank);
	white_img = (uint16_t***)multialloc(sizeof(uint16_t), 3, white_count[0], white_count[1], white_count[2]);
	dark_img = (uint16_t***)multialloc(sizeof(uint16_t), 3, dark_count[0], dark_count[1], dark_count[2]);
	proj_img = (uint16_t***)multialloc(sizeof(uint16_t), 3, proj_count[0], proj_count[1], proj_count[2]);
	white_2D_img = (Real_t**)multialloc(sizeof(Real_t), 2, white_count[2], white_count[1]);
	dark_2D_img = (Real_t**)multialloc(sizeof(Real_t), 2, dark_count[2], dark_count[1]);
	wd_dwnsmpl_img = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, total_t_slices);
	Projection = (Real_t***)multialloc(sizeof(Real_t), 3, SinogramPtr->N_p, SinogramPtr->N_r, total_t_slices);
	Weight = (Real_t***)multialloc(sizeof(Real_t), 3, SinogramPtr->N_p, SinogramPtr->N_r, total_t_slices);

	/*Selects ROI in the dataset which should be read into arrays*/
    	status = H5Sselect_hyperslab (white_dataspace, H5S_SELECT_SET, white_offset, NULL, white_count, NULL);
    	status = H5Sselect_hyperslab (dark_dataspace, H5S_SELECT_SET, dark_offset, NULL, dark_count, NULL);
    	status = H5Sselect_hyperslab (proj_dataspace, H5S_SELECT_SET, proj_offset, NULL, proj_count, NULL);
    	white_memspace = H5Screate_simple (3, white_count, NULL);  
    	dark_memspace = H5Screate_simple (3, dark_count, NULL);  
    	proj_memspace = H5Screate_simple (3, proj_count, NULL);   
	mem_offset[0] = 0; mem_offset[1] = 0; mem_offset[2] = 0;
    	status = H5Sselect_hyperslab (white_memspace, H5S_SELECT_SET, mem_offset, NULL, white_count, NULL);
    	status = H5Sselect_hyperslab (dark_memspace, H5S_SELECT_SET, mem_offset, NULL, dark_count, NULL);
    	status = H5Sselect_hyperslab (proj_memspace, H5S_SELECT_SET, mem_offset, NULL, proj_count, NULL);

	fprintf(TomoInputsPtr->debug_file_ptr,"gen_projection_4m_HDF: Reading HDF5 dataset ...\n");
    	status = H5Dread(white_dataset, H5T_NATIVE_UINT16, white_memspace, white_dataspace, H5P_DEFAULT, &(white_img[0][0][0]));
    	status = H5Dread(dark_dataset, H5T_NATIVE_UINT16, dark_memspace, dark_dataspace, H5P_DEFAULT, &(dark_img[0][0][0]));
    	status = H5Dread(proj_dataset, H5T_NATIVE_UINT16, proj_memspace, proj_dataspace, H5P_DEFAULT, &(proj_img[0][0][0]));

	fprintf(TomoInputsPtr->debug_file_ptr,"gen_projection_4m_HDF: Read the white, dark and projection datasets\n");
	fprintf(TomoInputsPtr->debug_file_ptr, "gen_projection_4m_HDF: ratio_r = %d, ratio_t = %d\n", ratio_r, ratio_t);
	for (j = 0; j < white_count[2]; j++)
	for (k = 0; k < white_count[1]; k++)
	{
		white_2D_img[j][k] = 0;
		for (i = 0; i < white_count[0]; i++)
		{
			white_2D_img[j][k] += white_img[i][k][j];
		}
		white_2D_img[j][k] /= white_count[0];
	}

	for (j = 0; j < dark_count[2]; j++)
        for (k = 0; k < dark_count[1]; k++)
        {
                dark_2D_img[j][k] = 0;
                for (i = 0; i < dark_count[0]; i++)
                {
                        dark_2D_img[j][k] += dark_img[i][k][j];
                }
                dark_2D_img[j][k] /= dark_count[0];
        }

	fprintf(TomoInputsPtr->debug_file_ptr,"gen_projection_4m_HDF: Generated the downsampled versions of white and dark images\n");
		
	for (i = 0; i < SinogramPtr->N_p; i++)
	for (j = 0; j < SinogramPtr->N_r; j++)
	for (k = 0; k < total_t_slices; k++)
	{
		Projection[i][j][k] = 0;
		Weight[i][j][k] = 0;
		wd_dwnsmpl_img[j][k] = 0;
		for (m = 0; m < ratio_r; m++)
			for (n = 0; n < ratio_t; n++)
			{
				Weight[i][j][k] += fabs(proj_img[i][k*ratio_t + n][j*ratio_r + m] - dark_2D_img[j*ratio_r + m][k*ratio_t + n]);
				wd_dwnsmpl_img[j][k] += fabs(white_2D_img[j*ratio_r + m][k*ratio_t + n] - dark_2D_img[j*ratio_r + m][k*ratio_t + n]);
			}
		temp = log((wd_dwnsmpl_img[j][k])/(Weight[i][j][k]));
		Projection[i][j][k] = BH_QUAD_COEF*temp*temp + temp;
	}

	fprintf(TomoInputsPtr->debug_file_ptr,"gen_projection_4m_HDF: Generated projections and weight data with beamhardening coefficient of %f\n", (float)BH_QUAD_COEF);
	Write2Bin (proj_filename, 1, SinogramPtr->N_p, SinogramPtr->N_r, total_t_slices, &(Projection[0][0][0]), TomoInputsPtr->debug_file_ptr);
	Write2Bin (weight_filename, 1, SinogramPtr->N_p, SinogramPtr->N_r, total_t_slices, &(Weight[0][0][0]), TomoInputsPtr->debug_file_ptr);
	Write2Bin (wd_filename, 1, 1, SinogramPtr->N_r, total_t_slices, &(wd_dwnsmpl_img[0][0]), TomoInputsPtr->debug_file_ptr);
	
	if (TomoInputsPtr->Write2Tiff == 1)
	{
		dim[0] = 1; dim[1] = SinogramPtr->N_p; dim[2] = SinogramPtr->N_r; dim[3] = total_t_slices;
		WriteMultiDimArray2Tiff (proj_filename, dim, 0, 3, 1, 2, &(Projection[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
		WriteMultiDimArray2Tiff (weight_filename, dim, 0, 3, 1, 2, &(Weight[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
		dim[0] = 1; dim[1] = 1; dim[2] = SinogramPtr->N_r; dim[3] = total_t_slices;
		/*WriteMultiDimArray2Tiff (wd_filename, dim, 0, 1, 2, 3, &(wd_dwnsmpl_img[0][0]), 0, TomoInputsPtr->debug_file_ptr);*/
		WriteMultiDimArray2Tiff (wd_filename, dim, 0, 1, 2, 3, &(wd_dwnsmpl_img[0][0]), 0, TomoInputsPtr->debug_file_ptr);
	}
	fprintf(TomoInputsPtr->debug_file_ptr,"gen_projection_4m_HDF: Wrote projections, weight and bright field data to binary files\n");
	
	multifree(white_img, 3);
	multifree(dark_img, 3);
	multifree(proj_img, 3);
	multifree(wd_dwnsmpl_img, 2);
	multifree(white_2D_img, 2);
	multifree(dark_2D_img, 2);
	multifree(Projection, 3);
	multifree(Weight, 3);

	H5Sclose(white_memspace);
	H5Sclose(dark_memspace);
	H5Sclose(proj_memspace);
	H5Sclose(white_dataspace);
	H5Sclose(dark_dataspace);
	H5Sclose(proj_dataspace);
	H5Dclose(white_dataset);
	H5Dclose(dark_dataset);
	H5Dclose(proj_dataset);
	H5Fclose(white_file_id);
	H5Fclose(dark_file_id);
	H5Fclose(data_file_id);
} 


void read_4m_HDF (TomoInputs* TomoInputsPtr, char filename[], Real_t* object, hsize_t start[], hsize_t num[])
{
	hid_t file_id, dataset;
	hid_t dataspace, memspace;
   	herr_t status;
	int32_t rank, i;

	/*HDF file pointers*/
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	/*dataset pointers*/
	dataset = H5Dopen(file_id, HDFPATH, H5P_DEFAULT);
	dataspace = H5Dget_space (dataset);    /* dataspace handle */
	/*Gives the number of dimensions in a dataset*/
    	rank = H5Sget_simple_extent_ndims (dataspace);
	
	hsize_t dims[rank];
	status = H5Sget_simple_extent_dims (dataspace, dims, NULL);
	for (i=0; i<rank; i++)
		if (dims[i] < start[i] + num[i])
			fprintf(TomoInputsPtr->debug_file_ptr, "read_4m_HDF: The size of dimension %d in %s is %zu and is less than the required %zu\n", i, HDFPATH, dims[i], start[i] + num[i]);

    	status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, start, NULL, num, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, dataspace, H5P_DEFAULT, object);

	H5Sclose(dataspace);
	H5Dclose(dataset);
	H5Fclose(file_id);
} 

void write_2_existHDF (TomoInputs* TomoInputsPtr, char filename[], Real_t* object, hsize_t start[], hsize_t num[])
{
	hid_t file_id, dataset;
	hid_t dataspace;
   	herr_t status;
	int32_t rank, i;

	/*HDF file pointers*/
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	/*dataset pointers*/
	dataset = H5Dopen(file_id, HDFPATH, H5P_DEFAULT);
	dataspace = H5Dget_space (dataset);    /* dataspace handle */
	/*Gives the number of dimensions in a dataset*/
    	rank = H5Sget_simple_extent_ndims (dataspace);
	
	hsize_t dims[rank];
	status = H5Sget_simple_extent_dims (dataspace, dims, NULL);
	for (i=0; i<rank; i++)
		if (dims[i] < start[i] + num[i])
			fprintf(TomoInputsPtr->debug_file_ptr, "read_4m_HDF: The size of dimension %d in %s is %zu and is less than the required %zu\n", i, HDFPATH, dims[i], start[i] + num[i]);

    	status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, start, NULL, num, NULL);
    	status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, dataspace, H5P_DEFAULT, object);

	H5Sclose(dataspace);
	H5Dclose(dataset);
	H5Fclose(file_id);
} 



void write_2_newHDF (TomoInputs* TomoInputsPtr, char filename[], Real_t* object, hsize_t start[], hsize_t num[], int32_t rank, hsize_t dims[])
{
	hid_t file_id, dataset;
	hid_t dataspace;
   	herr_t status;

	/*HDF file pointers*/
	/* Create a new file using default properties. */
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	/*dataset pointers*/
	/* Create the data space for the dataset. */
	dataspace = H5Screate_simple(rank, dims, NULL);
	  /* Create the dataset. */
	dataset = H5Dcreate1(file_id, HDFPATH, H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);	

    	status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, start, NULL, num, NULL);
    	status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, dataspace, H5P_DEFAULT, object);

	H5Sclose(dataspace);
	H5Dclose(dataset);
	H5Fclose(file_id);
} 


#endif

