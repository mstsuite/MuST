/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <mpi.h>
#include "LSMSCommunication.hpp"

// #define USE_ISEND

void initializeCommunication(LSMSCommunication &comm)
{
  MPI_Init(NULL,NULL);
  comm.comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm.comm, &comm.rank);
  MPI_Comm_size(comm.comm, &comm.size);
}

void initializeCommunication(LSMSCommunication &comm, MPI_Comm mpiCommunicator)
{
  comm.comm=mpiCommunicator;
  MPI_Comm_rank(comm.comm, &comm.rank);
  MPI_Comm_size(comm.comm, &comm.size);
}

void finalizeCommunication(void)
{
  MPI_Finalize();
}

void exitLSMS(LSMSCommunication &comm, int errorCode)
{
  MPI_Abort(comm.comm, errorCode);
}

void communicateParameters(LSMSCommunication &comm, LSMSSystemParameters &lsms, 
                           CrystalParameters &crystal, MixingParameters &mix,
                           AlloyMixingDesc &alloyDesc)
{
  const int s = sizeof(LSMSSystemParameters) + 9*sizeof(Real) + sizeof(int) 
                + 10 + sizeof(MixingParameters) + 5*sizeof(int)
                + sizeof(int); // <-- +1 for no. alloy classes
  char buf[s];
  int nalloy_classes;

  if(comm.rank==0)
  {
    int pos=0;
    MPI_Pack(lsms.systemid,80,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(lsms.title,80,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.lsmsMode,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(lsms.potential_file_in,128,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(lsms.potential_file_out,128,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.pot_in_type,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.pot_out_type,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.alloy_in_type,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.alloy_out_type,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(lsms.infoEvecFileIn,128,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(lsms.infoEvecFileOut,128,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(lsms.localAtomDataFile,128,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.num_atoms,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.nspin,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.relativity,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.nrelc,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.nrelv,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.n_spin_cant,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.n_spin_pola,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.mtasa,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.xcFunctional[0],numFunctionalIndices,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.fixRMT,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.nscf,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.writeSteps,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.temperature,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.clight,1,MPI_DOUBLE,buf,s,&pos,comm.comm);

    MPI_Pack(&lsms.energyContour.grid,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.energyContour.npts,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.energyContour.ebot,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.energyContour.etop,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.energyContour.eibot,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.energyContour.eitop,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.energyContour.maxGroupSize,1,MPI_INT,buf,s,&pos,comm.comm);

    MPI_Pack(&lsms.adjustContourBottom,1,MPI_DOUBLE,buf,s,&pos,comm.comm);

    MPI_Pack(&lsms.mixing,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.alphaDV,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.rmsTolerance,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.zblockLUSize,1,MPI_INT,buf,s,&pos,comm.comm);

    MPI_Pack(&lsms.global.iprpts,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.global.ipcore,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.global.iprint,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.global.print_node,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.global.default_iprint,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.global.istop,32,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(&lsms.global.GPUThreads,32,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack((int *)&lsms.global.linearSolver,32,MPI_INT,buf,s,&pos,comm.comm);

    MPI_Pack(&crystal.num_types,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&crystal.bravais(0,0),9,MPI_DOUBLE,buf,s,&pos,comm.comm);

    nalloy_classes = alloyDesc.size();
    MPI_Pack(&nalloy_classes,1,MPI_INT,buf,s,&pos,comm.comm);

// MixingParameters
    // MPI_CXX_BOOL is not always available
    // MPI_Pack(&mix.quantity[0],mix.numQuantities,MPI_CXX_BOOL,buf,s,&pos,comm.comm);
    // copy to temporary int array and send this
    int tmpQuantity[mix.numQuantities];
    for(int i=0; i<mix.numQuantities; i++)
      if(mix.quantity[i])
        tmpQuantity[i] = 1;
      else
        tmpQuantity[i] = 0; 
    MPI_Pack(&tmpQuantity[0],mix.numQuantities,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&mix.algorithm[0],mix.numQuantities,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&mix.mixingParameter[0],mix.numQuantities,MPI_DOUBLE,buf,s,&pos,comm.comm);
  }
  MPI_Bcast(buf,s,MPI_PACKED,0,comm.comm);
  if(comm.rank!=0)
  {
    int pos=0;
    MPI_Unpack(buf,s,&pos,lsms.systemid,80,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,lsms.title,80,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.lsmsMode,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,lsms.potential_file_in,128,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,lsms.potential_file_out,128,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.pot_in_type,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.pot_out_type,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.alloy_in_type,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.alloy_out_type,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,lsms.infoEvecFileIn,128,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,lsms.infoEvecFileOut,128,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,lsms.localAtomDataFile,128,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.num_atoms,1,MPI_INT,comm.comm);
    crystal.num_atoms=lsms.num_atoms;
    MPI_Unpack(buf,s,&pos,&lsms.nspin,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.relativity,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.nrelc,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.nrelv,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.n_spin_cant,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.n_spin_pola,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.mtasa,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.xcFunctional[0],numFunctionalIndices,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.fixRMT,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.nscf,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.writeSteps,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.temperature,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.clight,1,MPI_DOUBLE,comm.comm);

    MPI_Unpack(buf,s,&pos,&lsms.energyContour.grid,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.energyContour.npts,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.energyContour.ebot,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.energyContour.etop,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.energyContour.eibot,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.energyContour.eitop,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.energyContour.maxGroupSize,1,MPI_INT,comm.comm);

    MPI_Unpack(buf,s,&pos,&lsms.adjustContourBottom,1,MPI_DOUBLE,comm.comm);

    MPI_Unpack(buf,s,&pos,&lsms.mixing,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.alphaDV,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.rmsTolerance,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.zblockLUSize,1,MPI_INT,comm.comm);

    MPI_Unpack(buf,s,&pos,&lsms.global.iprpts,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.global.ipcore,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.global.iprint,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.global.print_node,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.global.default_iprint,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.global.istop,32,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,&lsms.global.GPUThreads,32,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,(int *)&lsms.global.linearSolver,32,MPI_INT,comm.comm);

    MPI_Unpack(buf,s,&pos,&crystal.num_types,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&crystal.bravais(0,0),9,MPI_DOUBLE,comm.comm);
    crystal.resize(crystal.num_atoms);
    crystal.resizeTypes(crystal.num_types);

    MPI_Unpack(buf,s,&pos,&nalloy_classes,1,MPI_INT,comm.comm);
    alloyDesc.resize(nalloy_classes);

// MixingParameters
    // MPI_CXX_BOOL is not always available
    // MPI_Unpack(buf,s,&pos,&mix.quantity[0],mix.numQuantities,MPI_CXX_BOOL,comm.comm);
    // recieve temporary int array and copy
    int tmpQuantity[mix.numQuantities];
    MPI_Unpack(buf,s,&pos,&tmpQuantity[0],mix.numQuantities,MPI_INT,comm.comm);
    for(int i=0; i<mix.numQuantities; i++)
      if(tmpQuantity[i]==1)
        mix.quantity[i] = true;
      else
        mix.quantity[i] = false; 
    MPI_Unpack(buf,s,&pos,&mix.algorithm[0],mix.numQuantities,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&mix.mixingParameter[0],mix.numQuantities,MPI_DOUBLE,comm.comm);
  }
  lsms.rank = comm.rank;
  MPI_Bcast(&crystal.position(0,0),3*crystal.num_atoms,MPI_DOUBLE,0,comm.comm);
  MPI_Bcast(&crystal.evecs(0,0),3*crystal.num_atoms,MPI_DOUBLE,0,comm.comm);
  MPI_Bcast(&crystal.type[0],crystal.num_atoms,MPI_INT,0,comm.comm);

// This is dangerous and assumes homogeneous nodes:
  MPI_Bcast(&crystal.types[0],crystal.num_types*sizeof(AtomType),MPI_BYTE,0,comm.comm);

  // exchange size of alloy classes and components within each class
  int *ncomps = (int*) malloc(sizeof(int)*nalloy_classes);

  if( comm.rank == 0 )
    for(int i = 0; i < alloyDesc.size(); i++)
      ncomps[i] = alloyDesc[i].size();

  MPI_Bcast(ncomps,nalloy_classes,MPI_INT,0,comm.comm);

  if( comm.rank != 0 )
    for(int i = 0; i < alloyDesc.size(); i++)
      alloyDesc[i].resize(ncomps[i]);

  for(int i = 0; i < alloyDesc.size(); i++)
    MPI_Bcast(alloyDesc[i].data(),ncomps[i]*sizeof(AtomType),MPI_BYTE,0,comm.comm); // <-- see danger above
  free(ncomps);

// get maximum lmax
  crystal.maxlmax=0;
  for(int i=0; i<crystal.num_types; i++)
    if(crystal.types[i].lmax>crystal.maxlmax) crystal.maxlmax=crystal.types[i].lmax; 
  lsms.maxlmax=crystal.maxlmax;

// set lsms.commRank to comm.rank
  lsms.commRank = comm.rank;
}

void communicateSingleAtomData(LSMSCommunication &comm, int from, int to, int &local_id, AtomData &atom, int tag)
{
  const int maxPts = 3051;
  const int maxCore = 30;
  int s = sizeof(AtomData) + sizeof(Real)*(2*3*maxPts+2*maxCore) + sizeof(int)*3*2*maxCore + sizeof(int);
  char buf[s];
  int t;

  if(comm.rank == from)
  {
    int pos=0;
    MPI_Pack(&local_id,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.jmt,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.jws,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.xstart,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.rmt,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.rInscribed,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(atom.header,80,MPI_CHAR,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.alat,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.efermi,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.vdif,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.ztotss,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.zcorss,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.zsemss,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.zvalss,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.qtotws,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.mtotws,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(atom.evec,3,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(atom.evecNew,3,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(atom.evecOut,3,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(atom.xvalws,2,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.localEnergy,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.localMadelungEnergy,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.alloy_class,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.omegaMT,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.omegaWS,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.rws,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.lmax,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.nspin,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.forceZeroMoment,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.numc,1,MPI_INT,buf,s,&pos,comm.comm);
    t=atom.vr.n_row();
    MPI_Pack(&t,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.vr(0,0),t,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.vr(0,1),t,MPI_DOUBLE,buf,s,&pos,comm.comm);

    MPI_Pack(&atom.rhotot(0,0),t,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.rhotot(0,1),t,MPI_DOUBLE,buf,s,&pos,comm.comm);

    MPI_Pack(&atom.corden(0,0),t,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.corden(0,1),t,MPI_DOUBLE,buf,s,&pos,comm.comm);

    MPI_Pack(&atom.b_con[0],3,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.b_basis[0],9,MPI_DOUBLE,buf,s,&pos,comm.comm);
    t=atom.ec.n_row();
    MPI_Pack(&t,1,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.ec(0,0),t,MPI_DOUBLE,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.ec(0,1),t,MPI_DOUBLE,buf,s,&pos,comm.comm);

    MPI_Pack(&atom.nc(0,0),t,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.nc(0,1),t,MPI_INT,buf,s,&pos,comm.comm);

    MPI_Pack(&atom.lc(0,0),t,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.lc(0,1),t,MPI_INT,buf,s,&pos,comm.comm);

    MPI_Pack(&atom.kc(0,0),t,MPI_INT,buf,s,&pos,comm.comm);
    MPI_Pack(&atom.kc(0,1),t,MPI_INT,buf,s,&pos,comm.comm);


    MPI_Send(buf,s,MPI_PACKED,to,tag,comm.comm);
  }
  if(comm.rank==to)
  {
    MPI_Status status;
    MPI_Recv(buf,s,MPI_PACKED,from,tag,comm.comm,&status);

    int pos=0;
    MPI_Unpack(buf,s,&pos,&local_id,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.jmt,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.jws,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.xstart,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.rmt,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.rInscribed,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,atom.header,80,MPI_CHAR,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.alat,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.efermi,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.vdif,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.ztotss,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.zcorss,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.zsemss,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.zvalss,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.qtotws,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.mtotws,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,atom.evec,3,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,atom.evecNew,3,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,atom.evecOut,3,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,atom.xvalws,2,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.localEnergy,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.localMadelungEnergy,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.alloy_class,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.omegaMT,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.omegaWS,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.rws,1,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.lmax,1,MPI_INT,comm.comm);
    atom.kkrsz = (atom.lmax+1)*(atom.lmax+1);
    MPI_Unpack(buf,s,&pos,&atom.nspin,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.forceZeroMoment,1,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.numc,1,MPI_INT,comm.comm);

    MPI_Unpack(buf,s,&pos,&t,1,MPI_INT,comm.comm);
    if(t!=atom.vr.n_row()) atom.resizePotential(t);
    MPI_Unpack(buf,s,&pos,&atom.vr(0,0),t,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.vr(0,1),t,MPI_DOUBLE,comm.comm);

    MPI_Unpack(buf,s,&pos,&atom.rhotot(0,0),t,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.rhotot(0,1),t,MPI_DOUBLE,comm.comm);

    MPI_Unpack(buf,s,&pos,&atom.corden(0,0),t,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.corden(0,1),t,MPI_DOUBLE,comm.comm);

    MPI_Unpack(buf,s,&pos,&atom.b_con[0],3,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.b_basis[0],9,MPI_DOUBLE,comm.comm);

    MPI_Unpack(buf,s,&pos,&t,1,MPI_INT,comm.comm);
    if(t!=atom.nc.n_row()) atom.resizeCore(t);
    MPI_Unpack(buf,s,&pos,&atom.ec(0,0),t,MPI_DOUBLE,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.ec(0,1),t,MPI_DOUBLE,comm.comm);

    MPI_Unpack(buf,s,&pos,&atom.nc(0,0),t,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.nc(0,1),t,MPI_INT,comm.comm);

    MPI_Unpack(buf,s,&pos,&atom.lc(0,0),t,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.lc(0,1),t,MPI_INT,comm.comm);

    MPI_Unpack(buf,s,&pos,&atom.kc(0,0),t,MPI_INT,comm.comm);
    MPI_Unpack(buf,s,&pos,&atom.kc(0,1),t,MPI_INT,comm.comm);
  }
}


void communicatePotentialShiftParameters(LSMSCommunication &comm, PotentialShifter &ps)
{

  //const int s = sizeof(bool) + sizeof(double) * 2;
  const int s = sizeof(int) + sizeof(double) * 2;
  char buf[s];

  if (comm.rank==0)
  {
    int pos = 0;

    // MPI_CXX_BOOL is not always available
    //MPI_Pack(&ps.vSpinShiftFlag, 1, MPI_CXX_BOOL, buf, s, &pos, comm.comm);
    int tmpFlag = 0;
    if (ps.vSpinShiftFlag) tmpFlag = 1;
    
    MPI_Pack(&tmpFlag, 1, MPI_INT, buf, s, &pos, comm.comm);
    MPI_Pack(&ps.minShift, 1, MPI_DOUBLE, buf, s, &pos, comm.comm);
    MPI_Pack(&ps.maxShift, 1, MPI_DOUBLE, buf, s, &pos, comm.comm);
    
  }

  MPI_Bcast(buf, s, MPI_PACKED, 0, comm.comm);

  if (comm.rank!=0)
  {
    int pos = 0;

    // MPI_CXX_BOOL is not always available
    //MPI_Unpack(buf, s, &pos, &ps.vSpinShiftFlag, 1, MPI_CXX_BOOL, comm.comm);
    int tmpFlag = 0;

    MPI_Unpack(buf, s, &pos, &tmpFlag, 1, MPI_INT, comm.comm);
    if (tmpFlag)
      ps.vSpinShiftFlag = true;
    else
      ps.vSpinShiftFlag = false;
     
    MPI_Unpack(buf, s, &pos, &ps.minShift, 1, MPI_DOUBLE, comm.comm);
    MPI_Unpack(buf, s, &pos, &ps.maxShift, 1, MPI_DOUBLE, comm.comm);

  }  


}


void expectTmatCommunication(LSMSCommunication &comm, LocalTypeInfo &local)
{
// prepost all recieves for tmats from remote nodes
  for(int i=0; i<comm.numTmatFrom; i++)
  {
    int from=comm.tmatFrom[i].remoteNode;
    for(int j=0; j<comm.tmatFrom[i].numTmats; j++)
    {
      // printf("Node %d: expect tmat %d from %d\n",comm.rank,comm.tmatFrom[i].globalIdx[j],from);
      MPI_Irecv(&local.tmatStore(0,comm.tmatFrom[i].tmatStoreIdx[j]),2*local.lDimTmatStore,
                MPI_DOUBLE,from,comm.tmatFrom[i].globalIdx[j],comm.comm,
                &comm.tmatFrom[i].communicationRequest[j]);
    }
  }
}

void sendTmats(LSMSCommunication &comm, LocalTypeInfo &local)
{
  for(int i=0; i<comm.numTmatTo; i++)
  {
    int to=comm.tmatTo[i].remoteNode;
    for(int j=0; j<comm.tmatTo[i].numTmats; j++)
    {
      // printf("Node %d: send tmat %d to %d\n",comm.rank,comm.tmatTo[i].globalIdx[j],to);
#ifdef USE_ISEND
      MPI_Isend(&local.tmatStore(0,comm.tmatTo[i].tmatStoreIdx[j]),2*local.lDimTmatStore,
                MPI_DOUBLE,to,comm.tmatTo[i].globalIdx[j],comm.comm,
                &comm.tmatTo[i].communicationRequest[j]);
#else
      MPI_Send(&local.tmatStore(0,comm.tmatTo[i].tmatStoreIdx[j]),2*local.lDimTmatStore,
               MPI_DOUBLE,to,comm.tmatTo[i].globalIdx[j],comm.comm);
#endif      
    }
  }
}
void finalizeTmatCommunication(LSMSCommunication &comm)
{
  MPI_Status status;
  for(int i=0; i<comm.numTmatFrom; i++)
  {
    int from=comm.tmatFrom[i].remoteNode;
    for(int j=0; j<comm.tmatFrom[i].numTmats; j++)
    {
      // printf("Finalize recieve request %d from node %d\n",j,from);
      MPI_Wait(&comm.tmatFrom[i].communicationRequest[j],&status);
    }
  }
#ifdef USE_ISEND
  for(int i=0; i<comm.numTmatTo; i++)
  {
    int to=comm.tmatTo[i].remoteNode;
    for(int j=0; j<comm.tmatTo[i].numTmats; j++)
    {
      // printf("Finalize send request %d to node %d\n",j,to);
      MPI_Wait(&comm.tmatTo[i].communicationRequest[j],&status);
    }
  }
#endif
}

void printCommunicationInfo(FILE *f, LSMSCommunication &comm)
{
  fprintf(f,"Communication: rank no. %d of %d\n",comm.rank,comm.size);
  fprintf(f,"Sending tmats to %d remote nodes:\n",comm.numTmatTo);
  for(int i=0; i<comm.numTmatTo; i++)
  {
    fprintf(f,"Node %d :",comm.tmatTo[i].remoteNode);
    for(int j=0; j<comm.tmatTo[i].numTmats; j++)
      fprintf(f," %d[%d]",comm.tmatTo[i].globalIdx[j],comm.tmatTo[i].tmatStoreIdx[j]);
    fprintf(f,"\n");
  }
  fprintf(f,"Recieving tmats from %d remote nodes:\n",comm.numTmatFrom);
  for(int i=0; i<comm.numTmatFrom; i++)
  {
    fprintf(f,"Node %d :",comm.tmatFrom[i].remoteNode);
    for(int j=0; j<comm.tmatFrom[i].numTmats; j++)
      fprintf(f," %d[%d]",comm.tmatFrom[i].globalIdx[j],comm.tmatFrom[i].tmatStoreIdx[j]);
    fprintf(f,"\n");
  }
}


void globalAnd(LSMSCommunication &comm,bool &a)
{
  int ia=0;
  int r;
  if(a) ia=1;
  MPI_Allreduce(&ia,&r,1,MPI_INT,MPI_PROD,comm.comm);
  a=true;
  if(r==0) a=false;
}


long long calculateFomScale(LSMSCommunication &comm, LocalTypeInfo &local)
{
  // fomScale = \sum_#atoms (LIZ * (lmax+1)^2)^3
  long long fomLocal = 0;
  long long fom  = 0;

  for(int i=0; i<local.num_local; i++)
  {
    // printf("nrmat = %d\n",local.atom[i].nrmat);
    fomLocal += local.atom[i].nrmat*local.atom[i].nrmat*local.atom[i].nrmat;
  }
  // printf("fomLocal = %ld\n",fomLocal);
  MPI_Allreduce(&fomLocal,&fom,1,MPI_LONG_LONG,MPI_SUM,comm.comm);
  // printf("fom = %ld\n",fom);
  return fom;
}
