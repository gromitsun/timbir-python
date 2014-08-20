#include <stdio.h>
#include "allocate.h"
#include "hdf5.h"
#include <stdlib.h>

void create_dataset (hid_t FILE_RD, hid_t FILE_WR, char HDF_PATH[], int z_start_extract, int z_num_extract);
int main ()
{
	char Path2Source[] = "/scratch/lustreD/m/mohank/Argonne_Datasets/APS_Beamtime_042014/APS14_60_AlSi.hdf";
	char Path2Dest[] = "/scratch/lustreD/m/mohank/Argonne_Datasets/APS_Beamtime_042014/AlSi_data_c.hdf";
	int z_start_extract = 300;
	int z_num_extract = 8;
	int contain_whites = 1;
	int contain_darks = 1;

	char HDF_WHITES_PATH[] = "/exchange/data_white";
	char HDF_DARKS_PATH[] = "/exchange/data_dark";
	char HDF_DATA_PATH[] = "/exchange/data";
	hid_t FILE_RD, FILE_WR;
	FILE_RD = H5Fopen(Path2Source, H5F_ACC_RDONLY, H5P_DEFAULT);
	FILE_WR = H5Fcreate(Path2Dest, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (contain_whites == 1)
		create_dataset (FILE_RD, FILE_WR, HDF_WHITES_PATH, z_start_extract, z_num_extract);
/*
	if (contain_darks == 1)
		create_dataset (FILE_RD, FILE_WR, HDF_DARKS_PATH, z_start_extract, z_num_extract);
			
	create_dataset (FILE_RD, FILE_WR, HDF_DATA_PATH, z_start_extract, z_num_extract);
*/	
	H5Fclose(FILE_RD);
	H5Fclose(FILE_WR);

	return(0);
}

void create_dataset (hid_t FILE_RD, hid_t FILE_WR, char HDF_PATH[], int z_start_extract, int z_num_extract)
{
	hsize_t dims[3], offset[3], count[3];
	hid_t rd_dataset, wr_dataset, rd_dataspace, wr_dataspace, rd_memspace, datatype;
	herr_t status;
	int rank;
	uint16_t ***data;
	
	rd_dataset = H5Dopen(FILE_RD, HDF_PATH, H5P_DEFAULT);
	rd_dataspace = H5Dget_space (rd_dataset);    /* dataspace handle */
    	rank = H5Sget_simple_extent_ndims (rd_dataspace);
	printf("Rank of %s is %d\n", HDF_PATH, rank);
	status = H5Sget_simple_extent_dims (rd_dataspace, dims, NULL);
	printf("Size of %s is %zux%zux%zu\n", HDF_PATH, dims[0], dims[1], dims[2]);
	offset[0] = 0; offset[1] = z_start_extract; offset[2] = 0;
	count[0] = dims[0]; count[1] = z_num_extract; count[2] = dims[2];
    	status = H5Sselect_hyperslab (rd_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    	rd_memspace = H5Screate_simple (rank, count, NULL);  
	data = (uint16_t***)multialloc(sizeof(uint16_t), 3, count[0], count[1], count[2]);
    	status = H5Dread(rd_dataset, H5T_NATIVE_UINT16, rd_memspace, rd_dataspace, H5P_DEFAULT, &(data[0][0][0]));
	wr_dataspace = H5Screate_simple(rank, count, NULL);
	datatype = H5Tcopy(H5T_NATIVE_UINT16);
	wr_dataset = H5Dcreate2(FILE_WR, HDF_PATH, datatype, wr_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	
	offset[0] = 0; offset[1] = 0; offset[2] = 0;
    	/*status = H5Sselect_hyperslab (wr_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);*/
    	status = H5Dwrite(wr_dataset, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(data[0][0][0]));
	H5Sclose(rd_memspace);
	H5Sclose(rd_dataspace);
	H5Sclose(wr_dataspace);
	H5Dclose(rd_dataset);
	H5Dclose(wr_dataset);
	multifree(data,3);
}
