/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMSCOMMUNICATION_H
#define LSMSCOMMUNICATION_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <vector>

#include "Main/SystemParameters.hpp"

#include "SingleSite/AtomData.hpp"
#include "Potential/PotentialShifter.hpp"

class TmatCommType {
public:
  int remoteNode;
  int numTmats;
  std::vector<int> tmatStoreIdx;
  std::vector<int> globalIdx;
  std::vector<MPI_Request> communicationRequest;
};

class LSMSCommunication {
public:
  int rank;
  int size;
  MPI_Comm comm;

  int numTmatTo, numTmatFrom;
  std::vector<TmatCommType> tmatTo, tmatFrom;
};

// mixing.hpp needs to be included after the definition of LSMSCommunication, since some mix
#include "Main/mixing.hpp"

void initializeCommunication(LSMSCommunication &comm);
void initializeCommunication(LSMSCommunication &comm, MPI_Comm mpiCommunicator);
void finalizeCommunication(void);
void exitLSMS(LSMSCommunication &comm, int errorCode);

void communicateParameters(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                           CrystalParameters &crystal, MixingParameters &mix, AlloyMixingDesc &alloyDesc);

void communicateSingleAtomData(LSMSCommunication &comm, int from, int to,
                               int &local_id, AtomData &atom, int tag=0);

void communicatePotentialShiftParameters(LSMSCommunication &comm, PotentialShifter &ps);

void expectTmatCommunication(LSMSCommunication &comm, LocalTypeInfo &local);
void sendTmats(LSMSCommunication &comm, LocalTypeInfo &local);
void finalizeTmatCommunication(LSMSCommunication &comm);

void printCommunicationInfo(FILE *f, LSMSCommunication &comm);

template<typename T>
void globalMax(LSMSCommunication &comm,T &a)
{
  T r;
  MPI_Allreduce(&a,&r,1,TypeTraits<T>::mpiType(),MPI_MAX,comm.comm);
  a=r;
/*
  fprintf(stderr,"Unsupported type in globalMax\n");
  exit(2);
*/
}

void globalAnd(LSMSCommunication &comm,bool &a);

template<typename T>
void globalSum(LSMSCommunication &comm,T &a)
{
  T r;
  MPI_Allreduce(&a,&r,1,TypeTraits<T>::mpiType(),MPI_SUM,comm.comm);
  a=r;
/*
  fprintf(stderr,"Unsupported type in globalMax\n");
  exit(2);
*/
}

template<typename T>
void globalSum(LSMSCommunication &comm,T *a, int n)
{
  T *r=new T[n];
  MPI_Allreduce(a,r,n,TypeTraits<T>::mpiType(),MPI_SUM,comm.comm);
  for(int i=0; i<n; i++) a[i]=r[i];
  delete[] r;
/*
  fprintf(stderr,"Unsupported type in globalMax\n");
  exit(2);
*/
}

long long calculateFomScale(LSMSCommunication &comm, LocalTypeInfo &local);

#endif
