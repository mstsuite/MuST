/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <stdio.h>
#include <string>
#include "SystemParameters.hpp"
#include "Potential/getXCName.hpp"

const char *potentialTypeName[]=
{
  "HDF5 (LSMS_1 format)",
  "Text (BIGCELL format)"
};

void printLSMSGlobals(FILE *f,LSMSSystemParameters &lsms)
{
  fprintf(f,"LSMS Globals:\n");
  fprintf(f,"  iprpts=%d\n",lsms.global.iprpts);
  fprintf(f,"  ipcore=%d\n",lsms.global.ipcore);
  fprintf(f,"  iprint=%d\n",lsms.global.iprint);
  fprintf(f,"  print_node=%d\n",lsms.global.print_node);
  fprintf(f,"  default_iprint=%d\n",lsms.global.default_iprint);
  fprintf(f,"  istop=%32s\n",lsms.global.istop);
  if(lsms.zblockLUSize>0) fprintf(f,"  zblockLUSize=%d\n",lsms.zblockLUSize);
  fprintf(f,"  linearSolver=%d\n",lsms.global.linearSolver);
}

void printLSMSSystemParameters(FILE *f,LSMSSystemParameters &lsms)
{
  fprintf(f,"Title    : '%s'\n",lsms.title);
  fprintf(f,"Id       : '%s'\n",lsms.systemid);
  fprintf(f,"No. Atoms: %d\n",lsms.num_atoms);
  fprintf(f,"Mtasa    : %d\n",lsms.mtasa);
  if(lsms.pot_in_type>=0)
    fprintf(f,"Pot. In  : %s [%s]\n",lsms.potential_file_in,potentialTypeName[lsms.pot_in_type]);
  else
    fprintf(f,"Pot. In  : generated.\n");
  if(lsms.pot_out_type>=0)
    fprintf(f,"Pot. Out : %s [%s]\n",lsms.potential_file_out,potentialTypeName[lsms.pot_out_type]);
  else
    fprintf(f,"Pot. Out : not written.\n");
  fprintf(f,"spin_cant: %d\n",lsms.n_spin_cant);
  fprintf(f,"Relativity: ");
  switch(lsms.relativity)
  {
  case none: fprintf(f,"non relativistic\n"); break;
  case scalar: fprintf(f,"scalar relativistic\n"); break;
  case full: fprintf(f,"fully relativistic\n"); break;
  default: fprintf(f,"!!!!UNKNOWN!!!!\n");
  }
  fprintf(f,"xcFunctional:"); for(int i=0; i<numFunctionalIndices; i++) fprintf(f," %d",lsms.xcFunctional[i]);
  std::string name;
  getXCName(lsms,name);
  fprintf(f," [%s]",name.c_str());
  fprintf(f,"\n");
  fprintf(f,"Electron Temperature: %lgK\n",lsms.temperature);
  fprintf(f,"RMS Tolerance: %lg\n",lsms.rmsTolerance);
}

void printCrystalParameters(FILE *f, CrystalParameters &crystal)
{
  fprintf(f,"Bravais lattice : %lf %lf %lf\n",crystal.bravais(0,0),crystal.bravais(1,0),crystal.bravais(2,0));
  fprintf(f,"                : %lf %lf %lf\n",crystal.bravais(0,1),crystal.bravais(1,1),crystal.bravais(2,1));
  fprintf(f,"                : %lf %lf %lf\n",crystal.bravais(0,2),crystal.bravais(1,2),crystal.bravais(2,2));

  for(int i=0; i<crystal.num_atoms; i++)
  {
    fprintf(f,"Atom %7d : %lf %lf %lf\n",i+1,crystal.position(0,i),crystal.position(1,i),crystal.position(2,i));
    fprintf(f,"             : evec = %lf %lf %lf\n",crystal.evecs(0,i),crystal.evecs(1,i), crystal.evecs(2,i));
    fprintf(f,"             : type = %d\n",crystal.type[i]+1);
  }
  for(int i=0; i<crystal.num_types; i++)
  {
    fprintf(f,"Type %7d : %s [%d.%d]\n",i+1,crystal.types[i].name,crystal.types[i].node,crystal.types[i].local_id);
    fprintf(f,"             : pot_in_idx = %d (Potential file index for input)\n",crystal.types[i].pot_in_idx);
    fprintf(f,"             : lmax = %d, rLIZ = %lf\n",crystal.types[i].lmax,crystal.types[i].rLIZ);
    fprintf(f,"             : rsteps = %lf %lf %lf %lf\n",crystal.types[i].rsteps[0],crystal.types[i].rsteps[1],
                                                          crystal.types[i].rsteps[2],crystal.types[i].rsteps[3]);
  }
}

void printLocalTypeInfo(FILE *f, LocalTypeInfo &local)
{
  fprintf(f,"No. Local Atom Types : %d\n",local.num_local);
  for(int i=0; i<local.num_local; i++)
  {
    fprintf(f,"Local Atom %2d : %.80s\n",i,local.atom[i].header);
    fprintf(f,"              : %lf %lf %lf %lf\n",local.atom[i].vr(0,0),local.atom[i].vr(0,1),
            local.atom[i].vr(50,0), local.atom[i].vr(50,1));
    fprintf(f,"              : %lf %lf\n",local.atom[i].ec(0,0),local.atom[i].ec(0,1));
    fprintf(f,"              : n_per_type = %d\n",local.n_per_type[i]);
  }
}

void printAlloyParameters(FILE *f, AlloyMixingDesc &alloy)
{
  fprintf(f,"No. of Alloy Classes : %lu\n",alloy.size());
  for(int i = 0; i < alloy.size(); i++) {
    fprintf(f,"Alloy Class #%1d contains %lu components\n",i+1,alloy[i].size());
    for(int j = 0; j < alloy[i].size(); j++) {
      fprintf(f, "  atom : %s\n", alloy[i][j].name);
      fprintf(f, "    conc : %.8f\n", alloy[i][j].conc);
      fprintf(f, "    Zval : Z=%3d Zc=%3d Zs=%3d, Zv=%3d\n", 
        alloy[i][j].Z, alloy[i][j].Zc, alloy[i][j].Zs, alloy[i][j].Zv);
    }
  }
}

void printLIZInfo(FILE * f, AtomData &atom)
{
  fprintf(f," No.        dist       Global Idx      Store Idx\n");
  for(int i=0; i<atom.numLIZ; i++)
    fprintf(f,"%4d  %12.7lf %7d   %4d\n",i,atom.LIZDist[i],atom.LIZGlobalIdx[i],atom.LIZStoreIdx[i]);

} 

