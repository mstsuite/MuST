/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Misc/Coeficients.hpp"
#include "VORPOL.hpp"

void setupVorpol(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local,SphericalHarmonicsCoeficients &shc)
{
  int ipvp=50; // parameter from LSMS_1.9 vorpol.h
  int iprcrit=260;
  int ipnode=ipvp*(ipvp-1);
  int ipcorn=(ipvp*(ipvp-1)*(ipvp-2))/6;
  int ipedge=(ipvp*ipvp-1)/2;
  const Real sphereVolumeFactor=4.0*M_PI/3.0;

  std::vector<Real> atom_position_1(crystal.num_atoms);
  std::vector<Real> atom_position_2(crystal.num_atoms);
  std::vector<Real> atom_position_3(crystal.num_atoms);
  std::vector<Real> rad(crystal.num_atoms);
  for(int i=0; i<crystal.num_atoms; i++)
  {
    atom_position_1[i]=crystal.position(0,i);
    atom_position_2[i]=crystal.position(1,i);
    atom_position_3[i]=crystal.position(2,i);
    rad[i]=2.0; // ignore the scaling of atomic cell volumes for now
  }

  crystal.omega=
    (crystal.bravais(1,0)*crystal.bravais(2,1)
     -crystal.bravais(2,0)*crystal.bravais(1,1))*crystal.bravais(0,2)+
    (crystal.bravais(2,0)*crystal.bravais(0,1)
     -crystal.bravais(0,0)*crystal.bravais(2,1))*crystal.bravais(1,2)+
    (crystal.bravais(0,0)*crystal.bravais(1,1)
     -crystal.bravais(1,0)*crystal.bravais(0,1))*crystal.bravais(2,2);
  crystal.omega=std::abs(crystal.omega);

  for(int i=0; i<local.num_local; i++)
  {
    int lmax=2*local.atom[i].lmax;
    local.atom[i].rInscribed=-1.0;
    local.atom[i].voronoi.rInscribedSphere=-1.0;
    local.atom[i].voronoi.wylm.resize((2*lmax+1)*(lmax+1),lsms.ngaussr,iprcrit-1);
    local.atom[i].voronoi.gwwylm.resize(lsms.ngaussr,iprcrit-1);
    local.atom[i].voronoi.grwylm.resize(lsms.ngaussr,iprcrit-1);
    int my_atom=local.global_id[i]+1;
    int num_atoms=crystal.num_atoms;
    setup_vorpol_(&my_atom,&num_atoms,
                  &atom_position_1[0],
                  &atom_position_2[0], &atom_position_3[0],
                  &crystal.bravais(0,0),
                  &lmax,&shc.clm[0],&lsms.ngaussq,&lsms.ngaussr,
                  &local.atom[i].voronoi.rInscribedSphere,&local.atom[i].voronoi.omegaInt,
                  local.atom[i].voronoi.dipint,&rad[0],
                  &ipvp,&ipnode,&ipcorn,&ipedge,&iprcrit,
                  &local.atom[i].voronoi.gwwylm(0,0),&local.atom[i].voronoi.grwylm(0,0),
                  &local.atom[i].voronoi.ncrit,&local.atom[i].voronoi.wylm(0,0,0),
                  &lsms.global.iprint,lsms.global.istop,32);
    local.atom[i].rInscribed=local.atom[i].voronoi.rInscribedSphere;
// set rmt according to value of fixRMT
    if(lsms.fixRMT==0)
    {
      local.atom[i].rmt=local.atom[i].rInscribed;
      local.atom[i].generateNewMesh = true;
    }
    // do we need to change rmt for mtasa==1?

    local.atom[i].omegaMT=sphereVolumeFactor*std::pow(local.atom[i].rmt,3);
    local.atom[i].omegaWS=local.atom[i].voronoi.omegaInt+local.atom[i].omegaMT;
    local.atom[i].rws=std::pow(local.atom[i].omegaWS/sphereVolumeFactor,1.0/3.0);

    switch(lsms.mtasa)
    {
    case 1:
      local.atom[i].rmt=local.atom[i].rws;
      local.atom[i].omegaMT=local.atom[i].omegaWS;
      local.atom[i].rInscribed=local.atom[i].rws;
      break;
    case 2:
      local.atom[i].rmt=local.atom[i].rws;
      local.atom[i].omegaMT=local.atom[i].omegaWS;
      local.atom[i].rInscribed=local.atom[i].voronoi.rInscribedSphere;
      break;
    default: // MT
      local.atom[i].rInscribed=local.atom[i].voronoi.rInscribedSphere;
    }
 }
}


void calculateVolumes(LSMSCommunication &comm, LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{
  Real volumeTotal = 0.0;
  Real volumeMT = 0.0;
  const Real sphereVolumeFactor = 4.0 * M_PI / 3.0;
  
  for (int i=0; i<local.num_local; i++)
  {
    volumeTotal += local.atom[i].omegaWS * local.n_per_type[i];
    switch (lsms.mtasa)
    {
      case 1:
        volumeMT += local.atom[i].omegaMT;
        break;
      default:
        volumeMT += sphereVolumeFactor * std::pow(local.atom[i].rInscribed, 3);
    }
  }

  globalSum(comm, volumeTotal);
  globalSum(comm, volumeMT);

  lsms.volumeTotal = volumeTotal;
  lsms.volumeNorm = crystal.omega / volumeTotal;
  lsms.volumeInterstitial = volumeTotal - volumeMT;

  if(lsms.global.iprint >= 0)
  {
    printf("\n");
    printf("Total cell volume     = %20.13f\n", lsms.volumeTotal);
    printf("Volume renorm. factor = %20.13f\n", lsms.volumeNorm);
    printf("WS cell volume        = %20.13f\n", local.atom[0].omegaWS);
    printf("Interstitial volume   = %20.13f\n", lsms.volumeInterstitial);
    printf("WS sphere radius      = %20.13f\n", local.atom[0].rws);
  }
}
