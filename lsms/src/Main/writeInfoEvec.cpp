/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

// write out the magnetism and constraint info for each site
void writeSingleEvec(FILE *f,int z, int i, Real posX, Real posY, Real posZ, AtomData &atom)
{
// Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  e_mix  B_x B_y B_z  vSpinShift
  fprintf(f,"%3d %8d  %21.15lf  %21.15lf  %21.15lf  %12.6lf %12.6lf  %21.15lf  %21.15lf  %21.15lf  %6.2lf  %21.15lf  %21.15lf  %21.15lf  %8.4lf  %21.15lf  %21.15lf  %21.15lf\n",
          z,i, posX, posY, posZ,
          atom.qtotws, atom.mtotws,
          atom.evec[0], atom.evec[1], atom.evec[2],
          -1.0,
          atom.b_con[0], atom.b_con[1], atom.b_con[2],
          atom.vSpinShift,
          atom.evecNew[0], atom.evecNew[1], atom.evecNew[2]);
}

// write out the magnetism and constraint info for each site
void writeSingleLocalAtomData(FILE *f,int z, int i, Real posX, Real posY, Real posZ, AtomData &atom)
{
// Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  B_x B_y B_z  vSpinShift localVolume localEnergy
  fprintf(f,"%3d %8d  %21.15lf  %21.15lf  %21.15lf  %12.6lf %12.6lf  %21.15lf  %21.15lf  %21.15lf  %21.15lf  %21.15lf  %21.15lf  %8.4lf  %.15lf  %.15lf \n",
          z,i, posX, posY, posZ,
          atom.qtotws, atom.mtotws,
          atom.evec[0], atom.evec[1], atom.evec[2],
          atom.b_con[0], atom.b_con[1], atom.b_con[2],
          atom.vSpinShift,
          atom.omegaWS, atom.localEnergy+atom.localMadelungEnergy);
}

void readSingleEvec(FILE *f,int &z, int &i, Real &posX, Real &posY, Real &posZ, AtomData &atom)
{
  Real tmp1, tmp2;
// Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  e_mix  B_x B_y B_z  vSpinShift
  fscanf(f,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n"
         ,&z,&i, &posX, &posY, &posZ, &tmp1, &atom.mtotws,
         &atom.evec[0], &atom.evec[1], &atom.evec[2],
         &tmp2,
         &atom.b_con[0], &atom.b_con[1], &atom.b_con[2],
         &atom.vSpinShift,
         &atom.evecOut[0], &atom.evecOut[1], &atom.evecOut[2]);
}

int writeInfoEvec(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local, Real eband, const char *name)
{
  AtomData pot_data;

  if(comm.rank==0)
  {
    FILE *outf=fopen(name,"w");

    fprintf(outf,"%.15lf %.15lf %.15lf\n", lsms.totalEnergy, eband, lsms.chempot);
    
// loop over all atom types:
    for(int i=0; i<crystal.num_types; i++)
    {
      if(crystal.types[i].node==comm.rank)
      {
        writeSingleEvec(outf,crystal.types[i].Z,i,
                        crystal.position(0,crystal.types[i].first_instance), // posX
                        crystal.position(1,crystal.types[i].first_instance), // posY
                        crystal.position(2,crystal.types[i].first_instance), // posZ
                        local.atom[crystal.types[i].local_id]);
      } else {
        int local_id;
        communicateSingleAtomData(comm, crystal.types[i].node, comm.rank, local_id, pot_data,i);
        if(local_id!=crystal.types[i].local_id) printf("WARNING: local_id doesn't match in writePotentials!\n");
        writeSingleEvec(outf,crystal.types[i].Z,i,
                        crystal.position(0,crystal.types[i].first_instance), // posX
                        crystal.position(1,crystal.types[i].first_instance), // posY
                        crystal.position(2,crystal.types[i].first_instance), // posZ
                        pot_data);
      }
    }
    fclose(outf);
  } else { // comm.rank!=0
    for(int i=0; i<local.num_local; i++)
    {
      communicateSingleAtomData(comm, comm.rank, 0, i, local.atom[i],local.global_id[i]);
    }
  }

  return 0;
}

int readInfoEvec(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local, const char *name)
{
  AtomData pot_data;

  int Z,ii;
  Real posX,posY,posZ;
  Real ef, etot, eband;

  if(comm.rank==0)
  {
    FILE *inf=fopen(name,"r");

    fscanf(inf,"%lf %lf %lf\n", &etot, &eband, &ef);
    
// loop over all atom types:
    for(int i=0; i<crystal.num_types; i++)
    {
      if(crystal.types[i].node==comm.rank)
      {
        readSingleEvec(inf,Z,ii,
                       posX,
                       posY,
                       posZ,
                       local.atom[crystal.types[i].local_id]);
// still need to set evec in crystal!
      } else {
        int local_id;
        if(local_id!=crystal.types[i].local_id) printf("WARNING: local_id doesn't match in writePotentials!\n");
        readSingleEvec(inf,Z,ii,
                       posX,
                       posY,
                       posZ,
                       pot_data);
        communicateSingleAtomData(comm, comm.rank, crystal.types[i].node, local_id, pot_data,ii);
      }
    }
    fclose(inf);
  } else { // comm.rank!=0
    for(int i=0; i<local.num_local; i++)
    {
      communicateSingleAtomData(comm, 0, comm.rank, i, local.atom[i],local.global_id[i]);
    }
  }

  return 0;
}

 int writeLocalAtomData(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local, Real eband, const char *name)
{
  AtomData pot_data;

  if(comm.rank==0)
  {
    FILE *outf=fopen(name,"w");

    fprintf(outf,"# totalEnergy  bandEnergy  fermiEnergy electrostaticEnergy\n");
    fprintf(outf,"# Z global_id x y z  qtotws  mtotws  evec_x evec_y evec_z  e_mix  B_x B_y B_z  vSpinShift localVolume localEnergy\n");
    
    fprintf(outf,"%.15lf %.15lf %.15lf %.15lf\n", lsms.totalEnergy, eband, lsms.chempot, lsms.u0);
    
// loop over all atom types:
    for(int i=0; i<crystal.num_types; i++)
    {
      if(crystal.types[i].node==comm.rank)
      {
        writeSingleLocalAtomData(outf,crystal.types[i].Z,i,
                        crystal.position(0,crystal.types[i].first_instance), // posX
                        crystal.position(1,crystal.types[i].first_instance), // posY
                        crystal.position(2,crystal.types[i].first_instance), // posZ
                        local.atom[crystal.types[i].local_id]);
      } else {
        int local_id;
        communicateSingleAtomData(comm, crystal.types[i].node, comm.rank, local_id, pot_data,i);
        if(local_id!=crystal.types[i].local_id) printf("WARNING: local_id doesn't match in writePotentials!\n");
        writeSingleLocalAtomData(outf,crystal.types[i].Z,i,
                        crystal.position(0,crystal.types[i].first_instance), // posX
                        crystal.position(1,crystal.types[i].first_instance), // posY
                        crystal.position(2,crystal.types[i].first_instance), // posZ
                        pot_data);
      }
    }
    fclose(outf);
  } else { // comm.rank!=0
    for(int i=0; i<local.num_local; i++)
    {
      communicateSingleAtomData(comm, comm.rank, 0, i, local.atom[i],local.global_id[i]);
    }
  }

  return 0;
}

