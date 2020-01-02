#ifndef REWLCOMMUNICATION_H
#define REWLCOMMUNICATION_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Real.hpp"
//#include "TypeTraits.hpp"

class REWLCommunication {

public :

  MPI_Comm comm;
  int rank;
  int size;

  // member functions
  void initializeCommunication(MPI_Comm mpiCommunicator);

  template<typename T> void swapScalar(T &data, int partner);
  template<typename T> void swapVector(T data[], int nElements, int partner);

  template<typename T> void sendScalar(T &data, int partner);
  template<typename T> void recvScalar(T &data, int partner);

};


// Definitions of template functions for MPI communications

template <typename T>
void REWLCommunication::swapScalar(T &data, int partner)
{

  MPI_Status status;
  MPI_Sendrecv_replace(&data, 1, TypeTraits<T>::mpiType(), partner, 1, partner, 1, comm, &status);

}


template <typename T>
void REWLCommunication::swapVector(T data[], int nElements, int partner)
{

  MPI_Status status;  
  MPI_Sendrecv_replace(&data[0], nElements, TypeTraits<T>::mpiType(), partner, 1, partner, 1, comm, &status);

}


template <typename T>
void REWLCommunication::sendScalar(T &data, int partner)
{
  MPI_Send(&data, 1, TypeTraits<T>::mpiType(), partner, 3, comm);
}


template <typename T>
void REWLCommunication::recvScalar(T &data, int partner)
{

  MPI_Status status;
  MPI_Recv(&data, 1, TypeTraits<T>::mpiType(), partner, 3, comm, &status);

}



#endif
