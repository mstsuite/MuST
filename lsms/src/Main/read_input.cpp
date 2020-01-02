/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <stdlib.h>

#include "lua.hpp"
//#include "lua.h"
//#include "lauxlib.h"
//#include "lualib.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "PhysicalConstants.hpp"

#include "SystemParameters.hpp"
#include "mixing.hpp"
#include "LSMSMode.hpp"
#include "LuaInterface/LuaSupport.hpp"
#include "../Potential/PotentialShifter.hpp"

#include <iostream>

void repeatBasisCell(LSMSSystemParameters &lsms, CrystalParameters &crystal, int nx, int ny, int nz, int unique)
{
  int numBasis=crystal.num_atoms;
  int numSites=numBasis*nx*ny*nz;

  Matrix<Real> basis,basis_evecs;
  std::vector<int> basis_type;
  basis=crystal.position;
  basis_evecs=crystal.evecs;
  basis_type=crystal.type;

  crystal.position.resize(3,numSites);
  crystal.evecs.resize(3,numSites);
  crystal.type.resize(numSites);

  if(unique)
  {
    crystal.types.resize(numSites);
    crystal.num_types=numSites;
    for(int i=0; i<numBasis; i++) if(crystal.types[i].pot_in_idx<0) crystal.types[i].pot_in_idx=i;
  }

  int i=0;
  for(int ix=0; ix<nx; ix++)
    for(int iy=0; iy<ny; iy++)
      for(int iz=0; iz<nz; iz++)
      {
        for(int ib=0; ib<numBasis; ib++)
        {
          crystal.position(0,i)=basis(0,ib)+ix*crystal.bravais(0,0)+iy*crystal.bravais(0,1)+iz*crystal.bravais(0,2);
          crystal.position(1,i)=basis(1,ib)+ix*crystal.bravais(1,0)+iy*crystal.bravais(1,1)+iz*crystal.bravais(1,2);
          crystal.position(2,i)=basis(2,ib)+ix*crystal.bravais(2,0)+iy*crystal.bravais(2,1)+iz*crystal.bravais(2,2);
          if(unique)
          {
            crystal.type[i]=i;
            if(i>=numBasis)
            {
              crystal.types[i]=crystal.types[ib];
              crystal.types[i].number_of_instances=1;
              crystal.types[i].first_instance=i;
            }
          } else {
            crystal.type[i]=basis_type[ib];
          }
          crystal.evecs(0,i)=basis_evecs(0,ib);
          crystal.evecs(1,i)=basis_evecs(1,ib);
          crystal.evecs(2,i)=basis_evecs(2,ib);
          i++;
        }
      }
  crystal.num_atoms=numSites;
  crystal.bravais(0,0)*=nx;
  crystal.bravais(1,0)*=nx;
  crystal.bravais(2,0)*=nx;
  crystal.bravais(0,1)*=ny;
  crystal.bravais(1,1)*=ny;
  crystal.bravais(2,1)*=ny;
  crystal.bravais(0,2)*=nz;
  crystal.bravais(1,2)*=nz;
  crystal.bravais(2,2)*=nz;

  lsms.num_atoms=numSites;
}

int readInput(lua_State *L, LSMSSystemParameters &lsms, CrystalParameters &crystal, MixingParameters &mix, PotentialShifter &potentialShifter,
    AlloyMixingDesc& alloyDesc)
{

// c     read in the structure  identity.................................
//       read(10,'(a)') systemid
//                      CHARACTER*50

  luaGetStrN(L,"systemid",lsms.systemid,50);
  snprintf(lsms.potential_file_in,128,"v_%s",lsms.systemid);
  snprintf(lsms.potential_file_out,128,"w_%s",lsms.systemid);
  luaGetStrN(L,"potential_file_in",lsms.potential_file_in,128);
  luaGetStrN(L,"potential_file_out",lsms.potential_file_out,128);
  lsms.pot_in_type=0; // HDF5
  luaGetInteger(L,"pot_in_type",&lsms.pot_in_type);
  lsms.pot_out_type=-1; // don't write potential
  luaGetInteger(L,"pot_out_type",&lsms.pot_out_type);

// c     read in the standard output target switch.......................
//       read(10,'(a)') output_to_screen  

// c     ================================================================
// c     read in subroutine stop level
//       read(10,'(a)') istop
// c     write(6,'(a)') istop

  char ctmp[32]; strncpy(ctmp,"main",32);
  luaGetStrN(L,"istop",ctmp,32); lsms.global.setIstop(ctmp);

  lsms.global.linearSolver=defaultSolver;
  luaGetInteger(L,"linearSolver",(int *)&lsms.global.linearSolver);

// c     ================================================================
// c     read in print level for a particular node and rest of the nodes.
//       read(10,*    ) node_print,print_instr,nprint
// c     write(6,'(3i5)') node_print,print_instr,nprint
  luaGetInteger(L,"iprpts",&lsms.global.iprpts);
  luaGetInteger(L,"ipcore",&lsms.global.ipcore);
  luaGetInteger(L,"print_node",&lsms.global.print_node);
  luaGetInteger(L,"default_iprint",&lsms.global.default_iprint);
  luaGetInteger(L,"iprint",&lsms.global.iprint);
#ifdef _OPENMP
  lsms.global.GPUThreads=std::min(16,omp_get_max_threads());
#else
  lsms.global.GPUThreads=1;
#endif
  luaGetInteger(L,"gpu_threads",&lsms.global.GPUThreads);
// c     ================================================================
// c     read in the number of atoms in the system.......................
//       read(10,*    ) num_atoms
  lsms.num_atoms=0;
  luaGetInteger(L,"num_atoms", &lsms.num_atoms);
//       if(num_atoms.lt.1 .or. num_atoms.gt.max_atoms) then
//          write(6,'(/,'' RDIN_ATOM_FT::'',
//      >               '' num_atoms exceeds the upper limit'')')
//          write(6,'(  ''                num_atoms:'',i5)')num_atoms
//          write(6,'(  ''                max_atoms:'',i5)')max_atoms
//          call fstop(sname)
//       endif

  lsms.relativity=scalar;
  lsms.nrelc=lsms.nrelv=0;
  char rel_str[80];
  rel_str[0]='s'; rel_str[1]=0; // defualt is scalar relativistic
  luaGetStrN(L,"relativity",rel_str,40);
  switch(rel_str[0])
  {
    case 'n': case 'N': lsms.relativity=none; lsms.nrelv=10; lsms.nrelc=10; break;
    case 's': case 'S': lsms.relativity=scalar; lsms.nrelv=0; lsms.nrelc=0; break;
    case 'f': case 'F': lsms.relativity=full; lsms.nrelv=0; lsms.nrelc=0; break;
  }
  rel_str[0]='d'; rel_str[1]=0; // default 
  luaGetStrN(L,"core_relativity",rel_str,40);
  switch(rel_str[0])
  {
    case 'n': case 'N': lsms.nrelc=10; break;
    case 'f': case 'F': lsms.nrelc=0; break;
  }
  lsms.clight=cphot*std::pow(10.0,lsms.nrelv);

  lsms.mtasa = 0;
  luaGetInteger(L,"mtasa",&lsms.mtasa);
  lsms.fixRMT = 0;
  luaGetInteger(L,"fixRMT",&lsms.fixRMT);
//       else if( mtasa .lt. -3 ) then
//            write(6,'(/,'' RDIN_ATOM_FT:: mtasa andor rmt0'',i5,d13.4)')
//      >     mtasa,rmt0
// 	   call fstop(sname)
//       endif

// c     read spin polarization index....................................
//       read(10,*    ) nspin,i_vdif,iexch
 
/* nspin = 1 : non spin polarized
           2 : spin polarized, collinear
           3 : spin polarized, non-collinear
           4 : fully relativistic
*/
  lsms.nspin=3;
  luaGetInteger(L,"nspin",&lsms.nspin);
  if(lsms.relativity==full)
  {
    lsms.nspin=3;
  }

  if(lsms.nspin>1) lsms.n_spin_pola=2; else lsms.n_spin_pola=1;
  if(lsms.nspin>2) lsms.n_spin_cant=2; else lsms.n_spin_cant=1;

// read exchange correlation functional specification:
  for(int i=0; i<numFunctionalIndices; i++) lsms.xcFunctional[i]=-1;
  lsms.xcFunctional[0]=0;
  lsms.xcFunctional[1]=1; // default functional is von Barth-Hedin
  for(int i=0; i<numFunctionalIndices; i++)
  {
    luaGetIntegerPositionInTable(L,"xcFunctional",i+1,&lsms.xcFunctional[i]);
  }


//       if (nspin.lt.0 .or.(nspin.le.3.and. nspin.gt.ipspin+1)) then
//          write(6,'('' RDIN_ATOM_FT:: Wrong input for nspin'')')
//          call fstop(sname)
//       endif

// Read Bravais Vectors:
  for(int i=0; i<3; i++)
  {
    luaGetPositionInTable(L,"bravais",i+1);
    for(int j=0; j<3; j++) luaGetRealPositionFromStack(L,j+1,&crystal.bravais(j,i));
    lua_pop(L,2);
  }

  /* printf("after reading bravais desc\n");
  luaStackDump(L); */


// Read atom positions and evecs
  crystal.resize(lsms.num_atoms);
  crystal.resizeTypes(lsms.num_atoms);
  crystal.num_atoms=lsms.num_atoms;
  crystal.num_types=0;
  for(int i=0; i<crystal.num_atoms; i++)
  {
    int t; // the type to be assigned

    luaGetPositionInTable(L,"site",i+1);
    luaGetFieldFromStack(L,"pos");
    for(int j=0; j<3; j++) luaGetRealPositionFromStack(L,j+1,&crystal.position(j,i));
    lua_pop(L,1);

    if(!luaGetIntegerFieldFromStack(L,"type",&t) || (t-1)==i)
    {
      luaGetStrNFromStack(L,"atom",crystal.types[crystal.num_types].name,4);
      luaGetIntegerFieldFromStack(L,"pot_in_idx",&crystal.types[crystal.num_types].pot_in_idx);
      luaGetIntegerFieldFromStack(L,"lmax",&crystal.types[crystal.num_types].lmax);
      luaGetIntegerFieldFromStack(L,"Z",&crystal.types[crystal.num_types].Z);
      luaGetIntegerFieldFromStack(L,"Zc",&crystal.types[crystal.num_types].Zc);
      luaGetIntegerFieldFromStack(L,"Zs",&crystal.types[crystal.num_types].Zs);
      luaGetIntegerFieldFromStack(L,"Zv",&crystal.types[crystal.num_types].Zv);
      luaGetIntegerFieldFromStack(L,"forceZeroMoment",&crystal.types[crystal.num_types].forceZeroMoment);
      luaGetIntegerFieldFromStack(L,"alloy_class",&crystal.types[crystal.num_types].alloy_class);
      crystal.types[crystal.num_types].alloy_class--; // <-- zero-based indexing
      luaGetRealFieldFromStack(L,"rLIZ",&crystal.types[crystal.num_types].rLIZ);
      luaGetFieldFromStack(L,"rsteps");
      for(int j=0; j<4; j++) luaGetRealPositionFromStack(L,j+1,&crystal.types[crystal.num_types].rsteps[j]);
      lua_pop(L,1);
      crystal.types[crystal.num_types].rad = 2.0;
      luaGetRealFieldFromStack(L,"rad",&crystal.types[crystal.num_types].rad);
      crystal.types[crystal.num_types].first_instance=i;
      crystal.types[crystal.num_types].number_of_instances=1;
      crystal.type[i]=crystal.num_types;
      crystal.num_types++;
    } else if(t<i && t>=0) {
      crystal.type[i]=crystal.type[t-1];
      crystal.types[crystal.type[i]].number_of_instances++;
    } else {
      fprintf(stderr,"Illegal type reference for atom %d : %d!\n",i+1,t);
      exit(1);
    }

    luaGetFieldFromStack(L,"evec");
    for(int j=0; j<3; j++) luaGetRealPositionFromStack(L,j+1,&crystal.evecs(j,i));
    lua_pop(L,1);
    lua_pop(L,2);
  }

  int xRepeat=1;
  luaGetInteger(L,"xRepeat",&xRepeat);
  int yRepeat=1;
  luaGetInteger(L,"yRepeat",&yRepeat);
  int zRepeat=1;
  luaGetInteger(L,"zRepeat",&zRepeat);
  int makeTypesUnique=1;
  luaGetInteger(L,"makeTypesUnique",&makeTypesUnique);

  repeatBasisCell(lsms, crystal, xRepeat, yRepeat, zRepeat,makeTypesUnique);

  /* printf("after reading atomic site desc\n");
  luaStackDump(L); */

  // print site values as check
  /* std::cout << "num_atoms = " << crystal.num_atoms << std::endl;
  for(int i = 0; i < crystal.num_atoms; i++) {
    std::cout << "atom[" << i <<  "] = " << crystal.types[i].name << std::endl;
    std::cout << "Z  = " << crystal.types[i].Z << std::endl;
  } */
 
  // for Wang-Landau for alloys

  // read in alloy file information
  snprintf(lsms.alloy_file_in,128,"bank_v_%s",lsms.systemid);
  snprintf(lsms.alloy_file_out,128,"bank_w_%s",lsms.systemid);
  luaGetStrN(L,"alloy_file_in",lsms.alloy_file_in,128);
  luaGetStrN(L,"alloy_file_out",lsms.alloy_file_out,128);
  lsms.alloy_in_type=0; // HDF5
  luaGetInteger(L,"alloy_in_type",&lsms.alloy_in_type);
  lsms.alloy_out_type=-1; // don't write potential
  luaGetInteger(L,"alloy_out_type",&lsms.alloy_out_type);

  // read in number of alloy classes
  int nalloy_class = 0;
  luaGetInteger(L,"nalloy_class",&nalloy_class);
  alloyDesc.resize(nalloy_class);

  // read in atomic species for each class
  for(int i = 0; i < alloyDesc.size(); i++) {
    luaGetPositionInTable(L,"alloy_class",i+1);

    int ncomp;
    luaGetIntegerFieldFromStack(L,"ncomp",&ncomp);
    alloyDesc[i].resize(ncomp);

    luaGetFieldFromStack(L,"comp");
    for(int j = 0; j < alloyDesc[i].size(); j++) {
      lua_pushinteger(L,j+1);
      lua_gettable(L,-2);

      luaGetStrNFromStack(L,"atom",alloyDesc[i][j].name,4);
      luaGetRealFieldFromStack(L,"conc",&alloyDesc[i][j].conc);
      luaGetIntegerFieldFromStack(L,"Z",&alloyDesc[i][j].Z);
      luaGetIntegerFieldFromStack(L,"Zc",&alloyDesc[i][j].Zc);
      luaGetIntegerFieldFromStack(L,"Zs",&alloyDesc[i][j].Zs);
      luaGetIntegerFieldFromStack(L,"Zv",&alloyDesc[i][j].Zv);
      lua_pop(L,1);
    } 
    lua_pop(L,3);
  } 

  // print values and quit to double-check
  /* std::cout << "nalloy_class = " << alloyDesc.size() << std::endl;
  for(int i = 0; i < alloyDesc.size(); i++) {
    std::cout << "class[" << i << "].ncomp = " << alloyDesc[i].size() << std::endl;
    for(int j = 0; j < alloyDesc[i].size(); j++) {
      std::cout << "atom = " << alloyDesc[i][j].name << std::endl;
      std::cout << "Z  = " << alloyDesc[i][j].Z << std::endl;
      std::cout << "Zc = " << alloyDesc[i][j].Zc << std::endl;
      std::cout << "Zs = " << alloyDesc[i][j].Zs << std::endl;
      std::cout << "Zv = " << alloyDesc[i][j].Zv << std::endl;
    }
  } 
  
  printf("after reading alloy desc\n");
  luaStackDump(L); */

// c     ================================================================
// c     read in a title to identify the system .........................
//       read(10,'(a)') system_title
  luaGetStrN(L,"system_title",lsms.title,80);
// c     ================================================================
// c     Read number of Gaussian points for r and theta integrations.....
//       read(10,*    ) ngaussr,ngaussq
// c     ================================================================
// c     Read in name of the atom........................................
// c     Read in the atom position vector................................
// c     Read in cut off radius for the LIZ of the atom..................
// c     Read in radius steps for lmax, lmax-1, lmax-2, lmax-3 ..........
// c     ================================================================
// c
// c     ================================================================
// c     read in names of info_table & info_evec files:..................
// c     ================================================================
//       read(10,'(a)') text
//       read(10,'(2a30)')info_table,info_evec

// c     ================================================================
// c     read in parameters that control energy inregration:.............
// c     ================================================================
// c     igrid   : specifies energy contour  1=slow......................
// c     igrid   : specifies energy contour  2=gaussian..................
// c     igrid   : specifies energy contour  3=Don's Fermi function poles
// c
// c     for igrid =1...[Zero Temperature Calculations Only].............
// c     ebot    : bottom of contour: real axis [may be mod. by semcor]..
// c     etop    : top of contour: real axis [usually reset to chempot]..
// c     eitop   : top    of contour on imaginary axis...................
// c     eibot   : bottom of contour on imaginary axis...................
// c     npts    : number of energy points per 0.1 ry....................
// c     kelvin  : not used..............................................
// c
// c     for igrid =2...[Zero Temperature Calculations Only].............
// c     ebot    : bottom of contour: real axis [may be mod. by semcor]..
// c     etop    : top of contour: real axis [usually reset to chempot]..
// c     eitop   : not used..............................................
// c     eibot   : not used..............................................
// c     npts    : number of Gaussian distributed energy points..........
// c     kelvin  : not used..............................................
// c
// c     for igrid =3...[Finite Temperature Calculations Only]...........
// c     ebot    : bottom of contour: real axis [may be mod. by semcor]..
// c                                            [then reset in congauss].
// c     etop    : top of contour on real axis [usually reset to chempot]
// c     eitop   : not used..............................................
// c     eibot   : not used..............................................
// c     npts    : not used..............................................
// c     kelvin  : temperature in kelvin..................................
// c     nument  : # of gaussian points on elliptical contour for Entropy
// c     ================================================================

  luaGetIntegerFieldInTable(L,"energyContour","grid",&lsms.energyContour.grid);
  luaGetRealFieldInTable(L,"energyContour","ebot",&lsms.energyContour.ebot);
  luaGetRealFieldInTable(L,"energyContour","etop",&lsms.energyContour.etop);
  luaGetRealFieldInTable(L,"energyContour","eibot",&lsms.energyContour.eibot);
  luaGetRealFieldInTable(L,"energyContour","eitop",&lsms.energyContour.eitop);
  luaGetIntegerFieldInTable(L,"energyContour","npts",&lsms.energyContour.npts);
  lsms.energyContour.maxGroupSize=50;
  luaGetIntegerFieldInTable(L,"energyContour","maxGroupSize",&lsms.energyContour.maxGroupSize);

  lsms.adjustContourBottom = -1.0;
  luaGetReal(L,"adjustContourBottom",&lsms.adjustContourBottom);

  lsms.temperature = 0.0;
  luaGetReal(L,"temperature",&lsms.temperature);

  // read in calculation mode. Default: main
  lsms.lsmsMode = LSMSMode::main;
  {
    char h[80];
    if(luaGetStrN(L,"lsmsMode",h,80))
    {
      if(strncmp(h,"main",80)==0) lsms.lsmsMode=LSMSMode::main;
      else if(strncmp(h,"liz0",80)==0) lsms.lsmsMode=LSMSMode::liz0;
      else
      {
        printf("Unknown lsmsMode: '%s'\n Defaulting to 'main'.\n",h);
        lsms.lsmsMode=LSMSMode::main;
      }
    }
  }

// c
// c     ================================================================
// c     read in controls for performing SCF calculation:................
// c     ================================================================
// c     nscf        : maximum number of scf iterations requested........
  lsms.nscf = 1;
  luaGetInteger(L, "nscf", &lsms.nscf);

// read the frequency of writing the potential during an scf calculation (writeSteps)
  lsms.writeSteps=30000;
  luaGetInteger(L, "writeSteps", &lsms.writeSteps);

// c     alpdv       : mixing parameter for chg. den. or potential.......
// c     alpma       : mixing parameter for moment density...............
// c     alpev       : mixing parameter for moment orientation...........
// c     mix_quant   : mixing charge density[potential] => 0[1]..........
// c     mix_algor   : mixing simple[DGAnderson] => 0[1].................
// lsms.mixing = 4*mix_algor+mix_quant; -1: no mixing
  lsms.mixing = -1;

  char quantity[80], algorithm[80];
  int numberOfMixQuantities = 0;

  // nullify all quantities; frozen potential will be set by default
  for (int i = 0; i < mix.numQuantities; i++)
  {
    mix.quantity[i] = false;
    mix.algorithm[i] = MixingParameters::noAlgorithm;
    mix.mixingParameter[i] = 0.0;
  }

  luaGetInteger(L, "numberOfMixQuantities", &numberOfMixQuantities);

  for (int i = 0; i < numberOfMixQuantities; i++)
  {
    luaGetPositionInTable(L, "mixing", i+1);
    luaGetStrNFromStack(L, "quantity", quantity, 50);

    int quantityIdx = -1;
    if (strcmp("no_mixing", quantity) == 0)
      quantityIdx = MixingParameters::no_mixing;
    else if (strcmp("charge", quantity) == 0) 
      quantityIdx = MixingParameters::charge;
    else if (strcmp("potential", quantity) == 0) 
      quantityIdx = MixingParameters::potential;
    else if (strcmp("moment_magnitude", quantity) == 0) 
      quantityIdx = MixingParameters::moment_magnitude;
    else if (strcmp("moment_direction", quantity) == 0)
      quantityIdx = MixingParameters::moment_direction;

    if (quantityIdx >= 0) mix.quantity[quantityIdx] = true;

    luaGetStrNFromStack(L, "algorithm", algorithm, 50);

    if (strcmp("simple", algorithm) == 0)
      mix.algorithm[quantityIdx] = MixingParameters::simple;
    else if (strcmp("broyden", algorithm) == 0)
      mix.algorithm[quantityIdx] = MixingParameters::broyden;

    luaGetRealFieldFromStack(L, "mixing_parameter", &mix.mixingParameter[quantityIdx]);

    lua_pop(L,2);

  }

  lsms.rmsTolerance = 1.0e-8;
  luaGetReal(L,"rmsTolerance",&lsms.rmsTolerance);

  int potentialShiftSwitch = 0;
  potentialShifter.vSpinShiftFlag = false;
// check if potentialShift has been assigned
  lua_getglobal(L,"potentialShift"); 
  if(lua_istable(L,-1))
  {
    lua_pop(L,1);
    luaGetIntegerFieldInTable(L, "potentialShift","switch", &potentialShiftSwitch);
    if (potentialShiftSwitch)
    {
      potentialShifter.vSpinShiftFlag = true;
      luaGetRealFieldInTable(L, "potentialShift","shift_min", &potentialShifter.minShift);
      luaGetRealFieldInTable(L, "potentialShift","shift_max", &potentialShifter.maxShift);
    }
  } else {
    lua_pop(L,1);
  }

  // check to read evec and constraints from file:
  lsms.infoEvecFileIn[0]=0;
  luaGetStrN(L,"infoEvecFileIn",lsms.infoEvecFileIn,120);
  snprintf(lsms.infoEvecFileOut,120,"info_evec_out");
  luaGetStrN(L,"infoEvecFileOut",lsms.infoEvecFileOut,120);

  lsms.localAtomDataFile[0]=0;
  luaGetStrN(L,"localAtomDataFile",lsms.localAtomDataFile,120);
  
  // read default block size for zblock_lu
  lsms.zblockLUSize=0;
  luaGetInteger(L,"zblockLUSize",&lsms.zblockLUSize);
// c     iharris = 0 : do not calculate harris energy....................
// c     iharris = 1 : calculate harris energy using updated chem. potl..
// c     iharris >=2 : calculate harris energy at fixed chem. potl.......
// c     i_potwrite  : the number of iterations between potential writes.
// c     movie   = 0 : no movie data will be written.....................
// c             = 1 : movie data will be written........................
// c     ctq         : coefficient of torque ............................
// c     ================================================================
//       read(10,'(a)') text
//       read(10, *   ) nscf,alpdv,alpma,alpev,mix_quant,mix_algor,
//      >               iharris,i_potwrite,movie
// c     ================================================================
// c     check consistencey of these parameters..........................
// c     ================================================================
//       if(i_potwrite.gt.nscf) then
//          i_potwrite=nscf
//       endif
//       if(i_potwrite.lt.-nscf) then
//          i_potwrite=-nscf
//       end if
//       if(mix_quant.ne.0 .and. mix_quant.ne.1) then
//          write(6,'('' RDIN_ATOM_FT::'',
//      >             '' Incorrect input data for mix_quant ='',
//      >             1i3)')mix_quant
//          call fstop(sname)
//       else if(mix_algor.gt.2) then 
//          write(6,'('' RDIN_ATOM_FT::'',
//      >             '' Incorrect input data for mix_algor ='',
//      >             1i3)')mix_algor
//          call fstop(sname)
//       endif
// c     ================================================================
// c     if calculating the harris energy make sure that mix_quant switch
// c     is set mix the charge density...................................
// c     ================================================================
//       if(iharris.ne.0) then
//          mix_quant=0
//       endif
// c
// c     ================================================================
// c     read in quantities that control Spin Dynamics :.................
// c     ================================================================
// c     nstep        : number of time steps [for SCF only : nstep=1 ]...
// c     tstep        : time step........................................
// c     etol         : over-ride scf convergence tolerence : energy ....
// c     ptol         : over-ride scf convergence tolerence : pressure ..
// c     eftol        : over-ride scf convergence tolerence : Fermi engy.
// c     rmstol       : over-ride scf convergence tolerence : rmstol ....
// c     ----------------------------------------------------------------
// c     etol,ptol,eftol, & rmstol can be relaxed for SD.................
// c     ================================================================
//       read(10,'(a)') text
// c     write(6,'(a70)')text
//       read(10,*) ntstep,tstep,etol,ptol,eftol,rmstol
// c     write(6,'('' RDIN_ATOM_FT: ntstep,tstep,etol,ptol,eftol,rmstol'',
// c    >i4,f10.4,3x,4d10.5)') ntstep,tstep,etol,ptol,eftol,rmstol
// c
// c     ================================================================
// c     read controls for calculation of torque & exchange interactions.
// c     ================================================================
// c     ctq        : 
// c     j_ij       : 
// c     ================================================================
//       read(10,'(a)') text
//       read(10, *   ) ctq,j_ij
// c     ================================================================
// c     check consistencey of these parameters..........................
// c     ================================================================
//       if(ctq.lt.0.0d0) then
//          write(6,'('' RDIN_ATOM_FT::'',
//      >             '' Incorrect input data for ctq < 0.0'',
//      >             1f10.5)')ctq
//          call fstop(sname)
//       else if(j_ij.ne.0 .and. j_ij.ne.1) then
//          write(6,'('' RDIN_ATOM_FT::'',
//      >             '' Incorrect input data for j_ij <> 0, and <> 1'',
//      >             1i5)')j_ij
//          call fstop(sname)
//       endif
 
  /* printf("End of read input\n"); 
  luaStackDump(L); */
  return 0;
}
