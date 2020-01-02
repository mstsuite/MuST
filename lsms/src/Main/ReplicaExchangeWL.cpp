#include <mpi.h>
#include <cmath>
#include <cstdio>
#include "Communication/REWLCommunication.hpp"
#include "ReplicaExchangeWL.hpp"

// Constructor
REWL::REWL(MPI_Comm comm, int numberOfWindows, int numberOfWalkersPerWindow, unsigned rngSeed)
{

  REWLcomm.initializeCommunication(comm);
  numWalkers = REWLcomm.size;
  numWindows = numberOfWindows;
  numWalkersPerWindow = numberOfWalkersPerWindow;

  myID = REWLcomm.rank;
  myWindow = int( floor(double(myID) / double(numWalkersPerWindow)) );

  //printf("YingWai's check: Inside REWL constructor. world_rank = %5d, myWindow = %5d, myID = %5d\n", global_rank, myWindow, myID);

  partnerID = -1;
  partnerWindow = -1;

  swapDirection = 0;
  upExchanges = 0;
  downExchanges = 0;

  // Initialize random number generator
  rng.seed(rngSeed);

}


// Destructor
REWL::~REWL()
{

}

//Member functions

void REWL::assignSwapPartner()
{

  partnerID = -1;
  partnerWindow = -1;

  // Find the partner's window
  switch (swapDirection) {
  case 0 :                              // 0-1 & 2-3 & .... & [numprocs-1]-nobody
  {
    if ((myWindow % 2 == 0) && (myWindow < (numWindows - 1))) {
      partnerWindow = myWindow + 1;
      upExchanges++;
    }   
    else if (myWindow % 2 != 0) {
      partnerWindow = myWindow - 1;
      downExchanges++;
    }
    swapDirection = 1;
    //printf("YingWai's check: Walker. %5d, partner is %5d, swapDirection is %5d\n", myWindow, partnerWindow, swapDirection);
    break;
  }
  case 1 :                              // 0-nobody & 1-2 & 3-4 & ....
  {
    if ((myWindow % 2 == 0) && (myWindow != 0)) {
      partnerWindow = myWindow - 1;
      downExchanges++;
    }   
    else if ((myWindow % 2 != 0) && (myWindow < (numWindows - 1))) {
      partnerWindow = myWindow + 1;
      upExchanges++;
    } 
    swapDirection = 0;
    //printf("YingWai's check: Walker. %5d, partner is %5d, swapDirection is %5d\n", myWindow, partnerWindow, swapDirection);
  }
  }

  // Find the partner's ID
  if (numWalkersPerWindow == 1) {
    partnerID = partnerWindow;
  }
  else {
    // to be implemented...
  }

}


void REWL::swapEnergy(double &energyForSwap)
{
 
  if (partnerID != -1)       // Swap energy with partner  
    REWLcomm.swapScalar(energyForSwap, partnerID);

}


bool REWL::determineAcceptance(double myDOSRatio)
{

  double partnerDOSRatio {0.0};
  int change {0};

  if (myWindow % 2 == 1) {          // Receiver and calculator 
    REWLcomm.recvScalar(partnerDOSRatio, partnerID);

    double acceptProb = myDOSRatio * partnerDOSRatio; 

    double rand = rnd(rng);
    if (rand < acceptProb) change = 1;

    REWLcomm.sendScalar(change, partnerID);
  }
  else {                            // Sender and non-calculator
    REWLcomm.sendScalar(myDOSRatio, partnerID);
    REWLcomm.recvScalar(change, partnerID);
  }

  if (change) return true;
  else return false;

}


void REWL::swapConfig(double evecsForSwap[], int numElements)
{
  if (partnerID != -1)
    REWLcomm.swapVector(evecsForSwap, numElements, partnerID);

}


void REWL::swapPotentialShifts(double potentialShiftsForSwap[], int numElements)
{
  if (partnerID != -1)
    REWLcomm.swapVector(potentialShiftsForSwap, numElements, partnerID);

}
