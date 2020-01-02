#include <mpi.h>
#include "REWLCommunication.hpp"

void REWLCommunication::initializeCommunication(MPI_Comm mpiCommunicator)
{
  comm = mpiCommunicator;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
}
