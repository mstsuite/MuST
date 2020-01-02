#include "Main/SystemParameters.hpp"
#include "LSMSCommunication.hpp"

int distributeTypes(CrystalParameters &crystal, LSMSCommunication &comm)
{
  int q=crystal.num_types/comm.size;
  int r=crystal.num_types%comm.size;

  int nodeNumber=0;
  int idx=0;
  for(int i=0; i<r; i++)
  {
    for(int j=0; j<=q; j++)
    {
      crystal.types[idx].node=nodeNumber;
      crystal.types[idx].local_id=j;
      idx++;
    }
    nodeNumber++;
  }

  for(int i=r; i<comm.size; i++)
  {
    for(int j=0; j<q; j++)
    {
      crystal.types[idx].node=nodeNumber;
      crystal.types[idx].local_id=j;
      idx++;
    }
    nodeNumber++;
  }

  if(idx!=crystal.num_types)
  {
    printf("Error distributing atoms to nodes! (This should not happen.)\n");
    exit(1);
  }

  int num_local=0;
  for(int i=0; i<crystal.num_types; i++)
    if(crystal.types[i].node==comm.rank) num_local++;

  return num_local;
}
