#include <vector>
#include <cmath>

#include <stdio.h>

#include "Real.hpp"
#include "SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

class LIZInfoType{
public:
  int idx;
  Real p1,p2,p3;
  Real dSqr;
};

// for buildLIZ see LSMS_1 neighbors_c.f
int buildLIZ(CrystalParameters &crystal, int idx,std::vector<LIZInfoType> &LIZ)
{
  const Real rtol=1.0e-8;
  const int n_max=5;
  int nrsclu=0;
  Real r1,r2,r3,p1,p2,p3,atdistsqr,shift_1,shift_2,shift_3;
  Real rcirclu=crystal.types[crystal.type[idx]].rLIZ;
  Real rcirclusqr=rcirclu*rcirclu;
  int n0,n1,n2,n3;
  int i[(2*n_max+1)*(2*n_max+1)*(2*n_max+1)];
  int j[(2*n_max+1)*(2*n_max+1)*(2*n_max+1)];
  int k[(2*n_max+1)*(2*n_max+1)*(2*n_max+1)];

  r1=std::sqrt(crystal.bravais(0,0)*crystal.bravais(0,0) +
               crystal.bravais(1,0)*crystal.bravais(1,0) +
               crystal.bravais(2,0)*crystal.bravais(2,0));
  r2=std::sqrt(crystal.bravais(0,1)*crystal.bravais(0,1) +
               crystal.bravais(1,1)*crystal.bravais(1,1) +
               crystal.bravais(2,1)*crystal.bravais(2,1));
  r3=std::sqrt(crystal.bravais(0,2)*crystal.bravais(0,2) +
               crystal.bravais(1,2)*crystal.bravais(1,2) +
               crystal.bravais(2,2)*crystal.bravais(2,2));

/*
  n1=std::max(n_max,int(rcirclu/r1+0.9));
  n2=std::max(n_max,int(rcirclu/r2+0.9));
  n3=std::max(n_max,int(rcirclu/r3+0.9));
*/
  n1=std::max(1,int(rcirclu/r1+0.9));
  n2=std::max(1,int(rcirclu/r2+0.9));
  n3=std::max(1,int(rcirclu/r3+0.9));

  if((2*n1+1)*(2*n2+1)*(2*n3+1) > (2*n_max+1)*(2*n_max+1)*(2*n_max+1))
  {
    fprintf(stderr,"FATAL ERROR: buildLIZ: n0 exeeds upper limit for site %d!",idx);
    exit(1);
  }
  n0=0;
  for(int i0=-n1; i0<=n1; i0++)
    for(int j0=-n2; j0<=n2; j0++)
      for(int k0=-n3; k0<=n3; k0++)
      {
        i[n0]=i0; j[n0]=j0; k[n0]=k0; n0++;
      }
  if(std::abs(rcirclu)<rtol) // special case: a one atom LIZ
  {
    LIZ[0].idx=idx; LIZ[0].dSqr=0.0;
    LIZ[0].p1=0.0; LIZ[0].p2=0.0; LIZ[0].p3=0.0;
    nrsclu=1;
    return nrsclu;
  } else {
    nrsclu=0;
    for(int m=0; m<n0; m++)
    {
      shift_1=Real(i[m])*crystal.bravais(0,0) +
        Real(j[m])*crystal.bravais(0,1) +
        Real(k[m])*crystal.bravais(0,2) -
        crystal.position(0,idx);
      shift_2=Real(i[m])*crystal.bravais(1,0) +
        Real(j[m])*crystal.bravais(1,1) +
        Real(k[m])*crystal.bravais(1,2) -
        crystal.position(1,idx);
      shift_3=Real(i[m])*crystal.bravais(2,0) +
        Real(j[m])*crystal.bravais(2,1) +
        Real(k[m])*crystal.bravais(2,2) -
        crystal.position(2,idx);
      for(int n=0; n<crystal.num_atoms; n++)
      {
        p1=shift_1+crystal.position(0,n);
        p2=shift_2+crystal.position(1,n);
        p3=shift_3+crystal.position(2,n);
        atdistsqr=p1*p1+p2*p2+p3*p3;
        if(atdistsqr<=rcirclusqr)
        {
          LIZ[nrsclu].idx=n;
          LIZ[nrsclu].p1=p1; LIZ[nrsclu].p2=p2; LIZ[nrsclu].p3=p3;
          LIZ[nrsclu++].dSqr=atdistsqr;
        }
      }
    }
  }
  return nrsclu;
}

bool nodeIsInList(int val, std::vector<LIZInfoType> &list, int len, CrystalParameters &crystal, int &ret)
{
  for(int i=ret; i<len; i++) if(crystal.types[crystal.type[list[i].idx]].node==val) {ret=i; return true;}
  return false;
}

class NodeIdxInfo {
public:
  int node, localIdx, globalIdx;
};
 
bool nodeLess_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y) {return x.node<y.node;}
bool globalLess_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y) {return x.globalIdx<y.globalIdx;}
bool localLess_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y) {return x.localIdx<y.localIdx;}
bool globalEq_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y) {return x.globalIdx==y.globalIdx;}
bool globalAndNodeEq_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y)
{return (x.globalIdx==y.globalIdx) && (x.node==y.node);}
bool localAndNodeEq_NodeIndexInfo(const NodeIdxInfo &x, const NodeIdxInfo &y)
{return (x.localIdx==y.localIdx) && (x.node==y.node);}
bool dSqrLess_LIZInfoType(const LIZInfoType &x, const LIZInfoType &y) {return x.dSqr<y.dSqr;}

void buildLIZandCommLists(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                          CrystalParameters &crystal, LocalTypeInfo &local)
{
  int numNodes = comm.size;
  std::vector<NodeIdxInfo> toList, fromList;
  std::vector<LIZInfoType> tempLIZ;
  tempLIZ.resize(4096);
  int tempNumLIZ;
  int ret;
  int fromListN,toListN;
  int toCounts[4096];
  int fromCounts[4096];
  int num_store=local.num_local;

// the first num_local entries in tmatStore contain the local tmats
  local.tmatStoreGlobalIdx.resize(4096);

  for(int i=0; i<local.num_local; i++)
  {
    crystal.types[local.global_id[i]].store_id=i;
    local.tmatStoreGlobalIdx[i]=local.global_id[i];
  }

  for(int i=0; i<4096; i++) {toCounts[i]=fromCounts[i]=0;}

  fromListN=toListN=0;
  fromList.resize(4096*local.num_local);
  toList.resize(4096*local.num_local);
// loop over all sites:
  for(int i=0; i<crystal.num_atoms; i++)
  {
    int type_id=crystal.type[i];
    int node=crystal.types[type_id].node;
    
    
    if(node==comm.rank) // local atom type
    {
      tempNumLIZ=buildLIZ(crystal,i,tempLIZ);
// set LIZ
      int local_id=crystal.types[type_id].local_id;
      std::sort(tempLIZ.begin(),tempLIZ.begin()+tempNumLIZ,dSqrLess_LIZInfoType);
      local.atom[local_id].numLIZ=tempNumLIZ;
      local.atom[local_id].LIZGlobalIdx.resize(tempNumLIZ);
      local.atom[local_id].LIZStoreIdx.resize(tempNumLIZ);
      local.atom[local_id].LIZDist.resize(tempNumLIZ);
      local.atom[local_id].LIZlmax.resize(tempNumLIZ);
      local.atom[local_id].LIZPos.resize(3,tempNumLIZ);
      local.atom[local_id].nrmat=0;
      for(int j=0; j<tempNumLIZ; j++)
      {
        local.atom[local_id].LIZGlobalIdx[j]=tempLIZ[j].idx;
        local.atom[local_id].LIZDist[j]=std::sqrt(tempLIZ[j].dSqr);
        local.atom[local_id].LIZPos(0,j)=tempLIZ[j].p1;
        local.atom[local_id].LIZPos(1,j)=tempLIZ[j].p2;
        local.atom[local_id].LIZPos(2,j)=tempLIZ[j].p3;
// calculate the lmax for the various shells
        int lkeep=crystal.types[type_id].lmax;
        local.atom[local_id].lmax = lkeep;
        for(int n1=0; n1<4; n1++)
          if(local.atom[local_id].LIZDist[j]>crystal.types[type_id].rsteps[n1]) lkeep--;
        local.atom[local_id].LIZlmax[j]=lkeep;
        local.atom[local_id].nrmat+=(lkeep+1)*(lkeep+1);
// add to commTmatFrom
        if(crystal.types[crystal.type[tempLIZ[j].idx]].node!=comm.rank) // need this from remote node
        {
          fromList[fromListN].node=crystal.types[crystal.type[tempLIZ[j].idx]].node;
          fromList[fromListN].localIdx=crystal.types[crystal.type[tempLIZ[j].idx]].local_id;
          fromList[fromListN].globalIdx=crystal.type[tempLIZ[j].idx];
          fromListN++;
        }
      }
    } else { // non local atom
// before building the LIZ we should first check if it is actually needed
// i.e find min distance between atom and local atoms and compare with max liz radius
      tempNumLIZ=buildLIZ(crystal,i,tempLIZ);

      ret=0;
      while(nodeIsInList(comm.rank,tempLIZ,tempNumLIZ,crystal,ret)) // the remote node needs the tmat from our node
      {
        toList[toListN].node=crystal.types[crystal.type[i]].node;          // the node that needs a tmat from us
        toList[toListN].localIdx=crystal.types[crystal.type[tempLIZ[ret].idx]].local_id;  // our local index == tmatStore entry
        toList[toListN].globalIdx=crystal.type[i];                         // the type that needs this tmat
        toListN++; ret++;
      }
    }
  }
// sort toList and fromList
  std::sort(fromList.begin(),fromList.begin()+fromListN,globalLess_NodeIndexInfo);
  std::sort(toList.begin(),toList.begin()+toListN,localLess_NodeIndexInfo);
  std::stable_sort(toList.begin(),toList.begin()+toListN,nodeLess_NodeIndexInfo);
// remove duplicates (only works on sorted lists!):
// for FROM remove all duplicate atom types
  std::vector<NodeIdxInfo>::iterator it=std::unique(fromList.begin(),fromList.begin()+fromListN,globalEq_NodeIndexInfo);
  fromList.resize(it-fromList.begin());
// for TO remove all identical atoms to the same node
  it=std::unique(toList.begin(),toList.begin()+toListN,localAndNodeEq_NodeIndexInfo);
  toList.resize(it-toList.begin());

//  fromList.resize(fromListN);
//  toList.resize(toListN);
  std::sort(fromList.begin(),fromList.end(),nodeLess_NodeIndexInfo);
  std::sort(toList.begin(),toList.end(),nodeLess_NodeIndexInfo);

// count the nodes in toList and fromList
  if(lsms.global.iprint>0) printf("toList.size()=%zu\n",toList.size());
  int numToNodes=0;
  if(toList.size()>0)
  {
    numToNodes=1;
    int h=toList[0].node;
    for(int i=0; i<toList.size(); i++)
    {
      if(lsms.global.iprint>0) printf("toList[%d].node=%d .localIdx=%d .globalIdx=%d\n",i,toList[i].node,
                              toList[i].localIdx,toList[i].globalIdx);
      // if(h!=toList[i].node) {h=toList[i].node; numToNodes++;} else toCounts[numToNodes-1]++;
      if(h!=toList[i].node) {h=toList[i].node; numToNodes++; toCounts[numToNodes-1]=1;} else toCounts[numToNodes-1]++;
    }
  }

  if(lsms.global.iprint>0) printf("fromList.size()=%zu\n",fromList.size());
  int numFromNodes=0;
  if(fromList.size()>0)
  {
    numFromNodes=1;
    int h=fromList[0].node;
    for(int i=0; i<fromList.size(); i++)
    {
      if(lsms.global.iprint>0) printf("fromList[%d].node=%d .localIdx=%d .globalIdx=%d\n",i,fromList[i].node,
                              fromList[i].localIdx,fromList[i].globalIdx);
      if(h!=fromList[i].node) {h=fromList[i].node; numFromNodes++; fromCounts[numFromNodes-1]=1;} else fromCounts[numFromNodes-1]++;
    }
  }

  comm.numTmatTo=numToNodes;
  comm.tmatTo.resize(numToNodes);
  comm.numTmatFrom=numFromNodes;
  comm.tmatFrom.resize(numFromNodes);

  int k=0;
  for(int i=0; i<numToNodes; i++)
  {
    comm.tmatTo[i].tmatStoreIdx.resize(toCounts[i]);
    comm.tmatTo[i].globalIdx.resize(toCounts[i]);
    comm.tmatTo[i].communicationRequest.resize(toCounts[i]);
    comm.tmatTo[i].remoteNode=toList[k].node;
    comm.tmatTo[i].numTmats=toCounts[i];
    for(int j=0; j<toCounts[i]; j++)
    {
      comm.tmatTo[i].globalIdx[j]=local.global_id[toList[k].localIdx];
      comm.tmatTo[i].tmatStoreIdx[j]=toList[k++].localIdx;
    }
  }

  k=0;
  if(lsms.global.iprint>0) printf("numFromNodes=%d\n",numFromNodes);
  for(int i=0; i<numFromNodes; i++)
  {
    if(lsms.global.iprint>0) printf("fromCounts[%d]=%d\n",i,fromCounts[i]);
    comm.tmatFrom[i].tmatStoreIdx.resize(fromCounts[i]);
    comm.tmatFrom[i].globalIdx.resize(fromCounts[i]);
    comm.tmatFrom[i].communicationRequest.resize(fromCounts[i]);
    comm.tmatFrom[i].remoteNode=fromList[k].node;
    comm.tmatFrom[i].numTmats=fromCounts[i];
    for(int j=0; j<fromCounts[i]; j++)
    {
      int g=fromList[k++].globalIdx;
      comm.tmatFrom[i].globalIdx[j]=g;
      if(crystal.types[g].store_id<0)
      {
        local.tmatStoreGlobalIdx[num_store]=g;
        crystal.types[g].store_id=num_store++;
      }
//      if(comm.rank==0)
//        printf("  i=%d j=%d : g=%d num_store=%d crystal.types[g].store_id=%d\n",i,j,g,num_store,crystal.types[g].store_id);
      comm.tmatFrom[i].tmatStoreIdx[j]=crystal.types[g].store_id;
    }
  }
  local.tmatStoreGlobalIdx.resize(num_store);
  int kkrsz2=2*(crystal.maxlmax+1)*(crystal.maxlmax+1);
  local.blkSizeTmatStore=kkrsz2*kkrsz2;
  local.lDimTmatStore=local.blkSizeTmatStore*lsms.energyContour.groupSize();
  local.tmatStore.resize(local.lDimTmatStore,num_store);

// set the StorIdx for the local atom LIZs
  for(int i=0; i<local.num_local; i++)
    for(int j=0; j<local.atom[i].numLIZ; j++)
      local.atom[i].LIZStoreIdx[j]=crystal.types[crystal.type[local.atom[i].LIZGlobalIdx[j]]].store_id;
}
