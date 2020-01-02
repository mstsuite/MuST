#include <stdio.h>
#include <hdf5.h>
#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "SingleSite/AtomData.hpp"
#include "SingleSite/readSingleAtomData.hpp"
#include "SingleSite/writeSingleAtomData.hpp"
#include "PotentialIO.hpp"
#include "HDF5io.hpp"
#include "Main/initializeAtom.hpp"

int loadAlloyBank(LSMSCommunication &comm, LSMSSystemParameters &lsms, AlloyMixingDesc &alloyDesc, AlloyAtomBank &alloyBank) {

  // for Wang-Landau for metallic alloys
  // written as a copy from 'loadPotentials' routine

  // resize alloy bank to match input in alloy descriptor
  // allocate room for potential and energy levels
  alloyBank.resize(alloyDesc.size());
  for(int i = 0; i < alloyBank.size(); i++) {
    alloyBank[i].resize(alloyDesc[i].size());
    for(int j = 0; j < alloyBank[i].size(); j++) {
      alloyBank[i][j].resizePotential(lsms.global.iprpts);
      alloyBank[i][j].resizeCore(lsms.global.ipcore);
    }
  }
  
  if( lsms.alloy_in_type > 1 || lsms.alloy_in_type < -1 ) 
    return 1; // unknown potential type

  if( lsms.alloy_in_type == -1 ) 
    return initializeNewAlloyBank(comm,lsms,alloyDesc,alloyBank);

  int natoms = 0;
  for(int i = 0; i < alloyDesc.size(); i++)
    natoms += alloyDesc[i].size();

  if( comm.rank == 0 ) {

    hid_t fid,fid_1;
    int id,fname_l;
    char fname[256];

    if( lsms.alloy_in_type == 0 ) { // HDF5

      fid=H5Fopen(lsms.alloy_file_in,H5F_ACC_RDONLY,H5P_DEFAULT);
      if( fid < 0 ) {
        printf("loadAlloyBank can't open HDF5 file '%s'\n",lsms.alloy_file_in);
        exit(1);
      }

      int i = -1;
      read_scalar<int>(fid,"LSMS",i);
      printf("Reading LSMS HDF5 input file format %d\n",i);
      if( i != 1 ) {
        printf("Attempting to read alloy bank with potential file version %d\nThis version of LSMS reads version 1 only!\n",i);
        exit(1);
      }

      i = -1;
      read_scalar<int>(fid,"NAtoms",i);
      printf("Reading data for %d atoms.\n",i);
      if( i != natoms ) {
        printf("Attempting to read alloy bank with potentials for %d atoms.\nPotential file contains %d atoms!\n",natoms,i);
        exit(1);
      }
    }

    // loop over all alloy mixing classes
    // loop over all components within a mixing class
    for(int id = 0, i = 0; i < alloyDesc.size(); i++) 
    for(int j = 0; j < alloyDesc[i].size(); j++, id++) {

      printf("reading alloy bank potential (id=%d) for atom type %d.\n",id,i);
      if( lsms.alloy_in_type == 0 ) // LSMS_1 style HDF5
      {
        // Atoms in the LSMS_1 HDF5 file are numbered starting from 000001
        snprintf(fname, 250, "%06d", id+1);
        // printf("Reading data from group '%s'\n",fname);
        fid_1 = H5Gopen2(fid, fname, H5P_DEFAULT);
        if( fid_1 < 0 ) {
          printf("Can't open group '%s'\n",fname);
          exit(1);
        }

        readSingleAtomData_hdf5(fid_1, alloyBank[i][j]);
        H5Gclose(fid_1);

      } else if( lsms.alloy_in_type == 1 ) { // BIGCELL style Text

        snprintf(fname,250,"%s.%d",lsms.alloy_file_in,id);
        fname_l=strlen(fname);
        printf("BIGCELL format file '%s' for alloy bank\n",fname);
        // fflush(stdout);

        readSingleAtomData_bigcell(fname, alloyBank[i][j]);
      }
    }

    // tune alloy potentials
    for(int id = 0, i = 0; i < alloyBank.size(); i++) 
    for(int j = 0; j < alloyBank[i].size(); j++, id++) { 
  
      AtomData &atom = alloyBank[i][j];
      atom.generateRadialMesh();
      atom.ztotss=(Real)alloyDesc[i][j].Z;
      atom.zcorss=(Real)alloyDesc[i][j].Zc;
      atom.zsemss=(Real)alloyDesc[i][j].Zs;
      atom.zvalss=(Real)alloyDesc[i][j].Zv;
      atom.lmax=alloyDesc[i][j].lmax;
      atom.kkrsz=(atom.lmax+1)*(atom.lmax+1);
      atom.alloy_class = i;
   
      // define potential past old grid as constant
      for(int is=0; is < lsms.n_spin_pola; is++) {
  
        Real vc = atom.vr(atom.jmt-1,is)/atom.r_mesh[atom.jmt];
        for(int ir = atom.jmt; ir < lsms.global.iprpts; ir++)
          atom.vr(ir,is) = 0;
      }
    }
  }

  // distribute alloy potentials over all nodes
  int dummy;
  if( comm.rank == 0 ) {

    for(int node = 1; node < comm.size; node++)
    for(int id = 0, i = 0; i < alloyBank.size(); i++) 
    for(int j = 0; j < alloyBank[i].size(); j++, id++) 
      communicateSingleAtomData(comm, 0, node, dummy, alloyBank[i][j], id);
  }
  else {

    for(int id = 0, i = 0; i < alloyBank.size(); i++) 
    for(int j = 0; j < alloyBank[i].size(); j++, id++) {
      communicateSingleAtomData(comm, 0, comm.rank, dummy, alloyBank[i][j], id);
      alloyBank[i][j].generateRadialMesh();
    }
  }

  return 0;
}


