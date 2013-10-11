#include <mpi.h>
#include "XT_Structures.h"
void MPI_Send_Recv_Z_Slices (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, MPI_Request* send_reqs, MPI_Request* recv_reqs, uint8_t select);
void MPI_Wait_Z_Slices (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, MPI_Request* send_reqs, MPI_Request* recv_reqs, uint8_t select);
