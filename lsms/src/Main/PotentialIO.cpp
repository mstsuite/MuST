/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
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

int loadPotentials(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{
  AtomData pot_data;
  if(lsms.pot_in_type>1 || lsms.pot_in_type<-1) return 1; // unknown potential type

  if(lsms.pot_in_type==-1)   initializeNewPotentials(comm,lsms,crystal,local); // start without potential files
//  if(lsms.pot_in_type==1) // we can't resize the data structures inside the fortran read routine
  {
    pot_data.resizePotential(lsms.global.iprpts);
    pot_data.resizeCore(lsms.global.ipcore);
  }

  // safe the input rmt, jmt, jws:
  Real rfixgrid[local.num_local];
  int jmtIn[local.num_local], jwsIn[local.num_local];
  for(int i=0; i<local.num_local; i++)
  {
    rfixgrid[i] = local.atom[i].rmt;
    jmtIn[i] = local.atom[i].jmt;
    jwsIn[i] = local.atom[i].jws;
  }

  if(comm.rank==0)
  {
    hid_t fid,fid_1;
    int id,fname_l;
    char fname[256];

    if(lsms.pot_in_type==0) // HDF5
    {
      fid=H5Fopen(lsms.potential_file_in,H5F_ACC_RDONLY,H5P_DEFAULT);
      if(fid<0)
      {
        printf("loadPotentials can't open HDF5 file '%s'\n",lsms.potential_file_in);
        exit(1);
      }
      int i=-1;
      read_scalar<int>(fid,"LSMS",i);
      printf("Reading LSMS HDF5 input file format %d\n",i);
      if(i!=1)
      {
        printf("Attempting to read potential file version %d\nThis version of LSMS reads version 1 only!\n",i);
        exit(1);
      }
      i=-1;
      read_scalar<int>(fid,"NAtoms",i);
      printf("Reading data for %d atoms.\n",i);
      if(i!=crystal.num_types)
      {
        printf("Attempting to read potentials for %d atoms.\nPotential file contains %d atoms!\n",crystal.num_types,i);
        exit(1);
      }
    }

// loop over all atom types:
    for(int i=0; i<crystal.num_types; i++)
    {
      id=crystal.types[i].pot_in_idx;
      if(id<0) id=i;
      // printf("reading potential (id=%d) for atom type %d.\n",id,i);
      if(lsms.pot_in_type==0) // LSMS_1 style HDF5
      {
        // Atoms in the LSMS_1 HDF5 file are numbered starting from 000001
        snprintf(fname,250,"%06d",i+1);
        fid_1=H5Gopen2(fid,fname,H5P_DEFAULT);
        // printf("Reading data from group '%s'\n",fname);
        if(fid_1<0)
        {
          printf("Can't open group '%s'\n",fname);
          exit(1);
        }
        if(crystal.types[i].node==comm.rank)
        {
          readSingleAtomData_hdf5(fid_1,local.atom[crystal.types[i].local_id]);
          local.atom[crystal.types[i].local_id].evec[0]=crystal.evecs(0,i);
          local.atom[crystal.types[i].local_id].evec[1]=crystal.evecs(1,i);
          local.atom[crystal.types[i].local_id].evec[2]=crystal.evecs(2,i);
        } else {
          readSingleAtomData_hdf5(fid_1,pot_data);
          pot_data.evec[0]=crystal.evecs(0,i);
          pot_data.evec[1]=crystal.evecs(1,i);
          pot_data.evec[2]=crystal.evecs(2,i);
          communicateSingleAtomData(comm, comm.rank, crystal.types[i].node, crystal.types[i].local_id, pot_data);
        }
        H5Gclose(fid_1);
      } else if(lsms.pot_in_type==1) { // BIGCELL style Text
        snprintf(fname,250,"%s.%d",lsms.potential_file_in,id);
        fname_l=strlen(fname);
        // printf("BIGCELL format file '%s'\n",fname);
        if(crystal.types[i].node==comm.rank)
        {
          readSingleAtomData_bigcell(fname,local.atom[crystal.types[i].local_id]);
          local.atom[crystal.types[i].local_id].evec[0]=crystal.evecs(0,i);
          local.atom[crystal.types[i].local_id].evec[1]=crystal.evecs(1,i);
          local.atom[crystal.types[i].local_id].evec[2]=crystal.evecs(2,i);
        } else {
          readSingleAtomData_bigcell(fname,pot_data);
          pot_data.evec[0]=crystal.evecs(0,i);
          pot_data.evec[1]=crystal.evecs(1,i);
          pot_data.evec[2]=crystal.evecs(2,i);
          communicateSingleAtomData(comm, comm.rank, crystal.types[i].node, crystal.types[i].local_id, pot_data);
        }
      }
    }
    // close the file
    if(lsms.pot_in_type==0)
    {
      H5Fclose(fid);
    }
  } else { // comm.rank!=0
    int local_id;
    for(int i=0; i<local.num_local; i++)
    {
      communicateSingleAtomData(comm, 0, comm.rank, local_id, pot_data);
      local.atom[local_id]=pot_data;
    }
  }

// adjust the local potentials:
  for(int i=0; i<local.num_local; i++)
  {
    if(local.atom[i].nspin != lsms.n_spin_pola)
    {
      printf("Input potential: %d spins != lsms setting: %d spins\n",
             local.atom[i].nspin, lsms.n_spin_pola);
      local.atom[i].changeNspin(lsms.n_spin_pola);
    }
    local.atom[i].forceZeroMoment = crystal.types[local.global_id[i]].forceZeroMoment;
    // if(local.atom[i].forceZeroMoment)
    // {
    //   local.atom[i].averageSpins();
    // }

    local.atom[i].generateRadialMesh();
    local.atom[i].ztotss=(Real)crystal.types[local.global_id[i]].Z;
    local.atom[i].zcorss=(Real)crystal.types[local.global_id[i]].Zc;
    local.atom[i].zsemss=(Real)crystal.types[local.global_id[i]].Zs;
    local.atom[i].zvalss=(Real)crystal.types[local.global_id[i]].Zv;
    local.atom[i].alloy_class=crystal.types[local.global_id[i]].alloy_class;
    local.atom[i].lmax=crystal.types[local.global_id[i]].lmax;
    local.atom[i].kkrsz=(local.atom[i].lmax+1)*
                        (local.atom[i].lmax+1);
    if(lsms.fixRMT==0)
    {
      local.atom[i].rmt=local.atom[i].rInscribed;
    }
// define potential past old grid as constant
    for(int is=0; is<lsms.n_spin_pola; is++)
    {
      Real vc=local.atom[i].vr(local.atom[i].jmt-1,is)/local.atom[i].r_mesh[local.atom[i].jmt];
      for(int ir=local.atom[i].jmt; ir<lsms.global.iprpts; ir++)
        local.atom[i].vr(ir,is)=0;
        //local.atom[i].vr(ir,is)=vc*local.atom[i].r_mesh[ir];
    }
  }

// calculate global properties: chempot, zvaltss, ...
  lsms.zvaltss=0.0;
  lsms.chempot=0.0;
  for(int i=0; i<local.num_local; i++)
  {
    lsms.zvaltss+=local.atom[i].zvalss*Real(local.n_per_type[i]);
    lsms.chempot+=local.atom[i].efermi*Real(local.n_per_type[i]);
  }
  Real fspace[2];
  fspace[0]=lsms.zvaltss;
  fspace[1]=lsms.chempot;
  
  globalSum(comm,fspace,2);
  
  lsms.zvaltss=fspace[0]; // /Real(lsms.num_atoms);
  lsms.chempot=fspace[1]/Real(lsms.num_atoms);
  
  return 0;
}

void initialAtomSetup(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{
  for(int i=0; i<local.num_local; i++)
  {
    // local.atom[i].generateRadialMesh();
    local.atom[i].ztotss=(Real)crystal.types[local.global_id[i]].Z;
    local.atom[i].zcorss=(Real)crystal.types[local.global_id[i]].Zc;
    local.atom[i].zsemss=(Real)crystal.types[local.global_id[i]].Zs;
    local.atom[i].zvalss=(Real)crystal.types[local.global_id[i]].Zv;
    local.atom[i].lmax=crystal.types[local.global_id[i]].lmax;
    local.atom[i].kkrsz=(local.atom[i].lmax+1)*
                        (local.atom[i].lmax+1);
  }

// calculate global properties: chempot, zvaltss, ...
  lsms.zvaltss=0.0;
  lsms.chempot=0.0;
  for(int i=0; i<local.num_local; i++)
  {
    lsms.zvaltss+=local.atom[i].zvalss*Real(local.n_per_type[i]);
    lsms.chempot+=local.atom[i].efermi*Real(local.n_per_type[i]);
  }
  Real fspace[2];
  fspace[0]=lsms.zvaltss;
  fspace[1]=lsms.chempot;
  
  globalSum(comm,fspace,2);
  
  lsms.zvaltss=fspace[0]; // /Real(lsms.num_atoms);
  lsms.chempot=fspace[1]/Real(lsms.num_atoms);
}

int writePotentials(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{
  AtomData pot_data;
  if(lsms.pot_out_type<0) return 0; // don't write potential
  if(lsms.pot_out_type>1) return 1; // unknown potential type

// update the Fermi energies
  for(int i=0; i<local.num_local; i++) local.atom[i].efermi=lsms.chempot;

  if(comm.rank==0)
  {
    hid_t fid,fid_1;
    int id,fname_l;
    char fname[256];

    if(lsms.pot_out_type==0) // HDF5
    {
      fid=H5Fcreate(lsms.potential_file_out,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      if(fid<0)
      {
        printf("writePotentials can't create HDF5 file '%s'\n",
               lsms.potential_file_out);
        exit(1);
      }
      // create LSMS and NAtoms tags for hdf5 file
      write_scalar<int>(fid,"LSMS",1);
      write_scalar<int>(fid,"NAtoms",crystal.num_types);

    }

// loop over all atom types:
    for(int i=0; i<crystal.num_types; i++)
    {
      id=i;
      // id=crystal.types[i].pot_out_idx; // this is more complicated!! need to average over all atoms with same id!
      // if(id<0) id=i;
      // printf("write potential (id=%d) for atom type %d.\n",id,i);
      if(lsms.pot_out_type==0) // LSMS_1 style HDF5
      {
        // Atoms in the LSMS_1 HDF5 file are numbered starting from 000001
        snprintf(fname,250,"%06d",i+1);
        fid_1=H5Gcreate2(fid,fname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        if(crystal.types[i].node==comm.rank)
        {
          writeSingleAtomData_hdf5(fid_1,local.atom[crystal.types[i].local_id],i+1);
        } else {
          communicateSingleAtomData(comm, crystal.types[i].node, comm.rank, crystal.types[i].local_id, pot_data,i);
          writeSingleAtomData_hdf5(fid_1,pot_data,i+1);
        }
        H5Gclose(fid_1);
      } else if(lsms.pot_out_type==1) { // BIGCELL style Text
        snprintf(fname,250,"%s.%d",lsms.potential_file_out,id);
        fname_l=strlen(fname);
        // printf("BIGCELL format file '%s'\n",fname);
        if(crystal.types[i].node==comm.rank)
        {
          writeSingleAtomData_bigcell(fname,local.atom[crystal.types[i].local_id]);
        } else {
          int local_id;
          communicateSingleAtomData(comm, crystal.types[i].node, comm.rank, local_id, pot_data,i);
          if(local_id!=crystal.types[i].local_id) printf("WARNING: local_id doesn't match in writePotentials!\n");
          writeSingleAtomData_bigcell(fname,pot_data);
        }
      }
    }
    // close the file
    if(lsms.pot_out_type==0)
    {
      H5Fclose(fid);
    }
    
  } else { // comm.rank!=0
    for(int i=0; i<local.num_local; i++)
    {
      communicateSingleAtomData(comm, comm.rank, 0, i, local.atom[i],local.global_id[i]);
    }
  }

  return 0;
}

