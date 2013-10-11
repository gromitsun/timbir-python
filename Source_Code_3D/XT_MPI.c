#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "XT_Structures.h"

void MPI_Send_Recv_Z_Slices (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, MPI_Request* send_reqs, MPI_Request* recv_reqs, uint8_t select)
{
	int32_t i, num, N_z, off1, off2;
	N_z = ScannedObjectPtr->N_z;
	num = ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;

	if (select == 0)
	{
		off1 = 0;
		off2 = 1;
	}
	else
	{
		off1 = 1;
		off2 = 0;
	}

	for (i = off1; i < ScannedObjectPtr->N_time; i = i + 2)
	{
		if (TomoInputsPtr->node_rank > 0)
			MPI_Isend(&(ScannedObjectPtr->Object[i][1][0][0]), num, MPI_DOUBLE, TomoInputsPtr->node_rank - 1, i, MPI_COMM_WORLD, &(send_reqs[i]));
		
		if (TomoInputsPtr->node_rank < TomoInputsPtr->node_num - 1)	
			MPI_Irecv(&(ScannedObjectPtr->Object[i][N_z+1][0][0]), num, MPI_DOUBLE, TomoInputsPtr->node_rank + 1, i, MPI_COMM_WORLD, &(recv_reqs[i]));	
	}
		
	for (i = off2; i < ScannedObjectPtr->N_time; i = i + 2)
	{
		if (TomoInputsPtr->node_rank > 0)
			MPI_Irecv(&(ScannedObjectPtr->Object[i][0][0][0]), num, MPI_DOUBLE, TomoInputsPtr->node_rank - 1, i, MPI_COMM_WORLD, &(recv_reqs[i]));	
		
		if (TomoInputsPtr->node_rank < TomoInputsPtr->node_num - 1)	
			MPI_Isend(&(ScannedObjectPtr->Object[i][N_z][0][0]), num, MPI_DOUBLE, TomoInputsPtr->node_rank + 1, i, MPI_COMM_WORLD, &(send_reqs[i]));
	}
}
			

void MPI_Wait_Z_Slices (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, MPI_Request* send_reqs, MPI_Request* recv_reqs, uint8_t select)
{
	int32_t i, num, N_z, off1, off2;
	N_z = ScannedObjectPtr->N_z;
	num = ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;

	if (select == 0)
	{
		off1 = 0;
		off2 = 1;
	}
	else
	{
		off1 = 1;
		off2 = 0;
	}

	for (i = off1; i < ScannedObjectPtr->N_time; i = i + 2)
	{
		if (TomoInputsPtr->node_rank > 0)
			MPI_Wait(&(send_reqs[i]), MPI_STATUSES_IGNORE);
		
		if (TomoInputsPtr->node_rank < TomoInputsPtr->node_num - 1)	
			MPI_Wait(&(recv_reqs[i]), MPI_STATUSES_IGNORE);
	}
		
	for (i = off2; i < ScannedObjectPtr->N_time; i = i + 2)
	{
		if (TomoInputsPtr->node_rank > 0)
			MPI_Wait(&(recv_reqs[i]), MPI_STATUSES_IGNORE);
		
		if (TomoInputsPtr->node_rank < TomoInputsPtr->node_num - 1)	
			MPI_Wait(&(send_reqs[i]), MPI_STATUSES_IGNORE);
	}
}
			


      
