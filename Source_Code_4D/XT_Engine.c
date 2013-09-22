#include <stdio.h>
#include <math.h>
#include "XT_Structures.h"
#include "XT_AMatrix.h"
#include "XT_genSinogram.h"
#include "allocate.h"
#include "XT_ICD_update.h"
#include "randlib.h"
#include "XT_Init.h"
#include "XT_Constants.h"
#include <time.h>
#include "XT_IOMisc.h"
#include <ctype.h>
#include <mpi.h>

void freeMemory(Sinogram* SinogramPtr, ScannedObject *ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	int32_t i;
	multifree(ScannedObjectPtr->ProjIdxPtr, 2);
	free(ScannedObjectPtr->ProjNum);

	for (i = 0; i < ScannedObjectPtr->N_time; i++)
		multifree(ScannedObjectPtr->Object[i],3);
	free(ScannedObjectPtr->Object);
/*	multifree(ScannedObjectPtr->Object,4);*/

	multifree(TomoInputsPtr->x_rand_select,3);
	multifree(TomoInputsPtr->y_rand_select,3);
	multifree(TomoInputsPtr->x_NHICD_select,3);
	multifree(TomoInputsPtr->y_NHICD_select,3);
	multifree(TomoInputsPtr->UpdateSelectNum,2);
	multifree(TomoInputsPtr->NHICDSelectNum,2);
	multifree(SinogramPtr->ProjOffset,2);
	multifree(SinogramPtr->ProjSelect,3);
	multifree(TomoInputsPtr->Weight,3);	
	free(SinogramPtr);
	free(ScannedObjectPtr);
	free(TomoInputsPtr);
}

/*computes sinogram*/
int computeWriteSinogram(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	char proj_file[100]=PROJECTION_FILENAME;
	char weight_file[100]=WEIGHT_MATRIX_FILENAME;
	int dim[4];

	sprintf(proj_file, "%s_n%d", proj_file, TomoInputsPtr->node_rank);
	sprintf(weight_file, "%s_n%d", weight_file, TomoInputsPtr->node_rank);
	if(TomoInputsPtr->sinobin == 1)
		genSinogram_fromBin(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	else if(TomoInputsPtr->sinobin == 2)
		genSinogramFromPhantom (SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	else if (TomoInputsPtr->sinobin == 3)
		genSinogram_subsampleCounts(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	else
	{
		printf("ERROR: computeWriteSinogram: Cannot recognize sinobin value\n");
		return (1);
	}

	dim[0] = 1;
	dim[1] = SinogramPtr->N_p;
	dim[2] = SinogramPtr->N_r;
	dim[3] = SinogramPtr->N_t;
	if (TomoInputsPtr->Write2Tiff == 1)	
	{
		WriteMultiDimArray2Tiff (proj_file, dim, 0, 3, 1, 2, &(SinogramPtr->Projection[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
		WriteMultiDimArray2Tiff (weight_file, dim, 0, 3, 1, 2, &(TomoInputsPtr->Weight[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	}

	return (0);
}

int main(int argc, char **argv)
{
	Sinogram *SinogramPtr;
	ScannedObject *ScannedObjectPtr;
	TomoInputs* TomoInputsPtr;
	time_t start;	
	int flag;

	start = time(NULL);
	srandom2(761521);
/*	srandom2(761522);*/

/*	fprintf(TomoInputsPtr->debug_file_ptr, "newline=%d, space=%d, '-'=%d, ','=%d, number=%d, alpha=%d\n", isalnum('\n'), isalnum(' '), isalnum('-'), isalnum(','), isalnum('1'), isalnum('a'));*/

	printf ("Running.....\n");
	MPI_Init(&argc, &argv);
	
	SinogramPtr = (Sinogram*)get_spc(1,sizeof(Sinogram));
	ScannedObjectPtr = (ScannedObject*)get_spc(1,sizeof(ScannedObject));
	TomoInputsPtr = (TomoInputs*)get_spc(1,sizeof(TomoInputs));
	MPI_Comm_size(MPI_COMM_WORLD, &(TomoInputsPtr->node_num));
	MPI_Comm_rank(MPI_COMM_WORLD, &(TomoInputsPtr->node_rank));
	
	argsParser (argc, argv, SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	initStructures (SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	fprintf(TomoInputsPtr->debug_file_ptr, "main: Time elapsed is %f minutes\n", difftime(time(NULL), start)/60.0);

	computeWriteSinogram(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	fprintf(TomoInputsPtr->debug_file_ptr, "main: Time elapsed is %f minutes\n", difftime(time(NULL), start)/60.0);

	if (TomoInputsPtr->reconstruct == 1) {
		flag=ICD_BackProject(SinogramPtr, ScannedObjectPtr, TomoInputsPtr); if(flag==1) return(1);
		fprintf(TomoInputsPtr->debug_file_ptr, "main: Time elapsed is %f minutes\n", difftime(time(NULL), start)/60.0);
		if (TomoInputsPtr->WritePerIter == 0)
			write_ObjectProjOff2TiffBinPerIter (SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	}	
	
	fprintf(TomoInputsPtr->debug_file_ptr, "XT_Engine: Will free remaining memory\n");
	fclose(TomoInputsPtr->debug_file_ptr);
	freeMemory(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	printf("Program will exit now\n");

	MPI_Finalize(); 
	return(0);
}
