#ifndef LSMS_ATOM_H
#define LSMS_ATOM_H

#include "RadialGrid.hpp"
#include "Real.hpp"
#include "Matrix.hpp"

typedef class Atom : public RadialPotential {
public:
  inline Atom(int iprpts, int ipcore)
  {
    g->generate(0.001,0.001,iprpts,1,1); sync();
    rhotot.resize(iprpts,2);
    nc.resize(ipcore,2); lc.resize(ipcore,2); kc.resize(ipcore,2);
    ec.resize(ipcore,2);
  }
  int id;
  int lmax,kkrsz;
  Real alat,efermi;
  Real ztotss, zcorss;
  Real xvalws[2];
  int nspin;
  Real evec[3];
  // RadialGrid *g;
  // RadialPotential *vr;
  Matrix<Real> rhotot;
  Real vdif;
  char header[80];
// Core state information
  int numc;
  Matrix<int> nc,lc,kc;
  Matrix<Real> ec;
};

//     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//       subroutine single_pot_read(present_atom,
//      >                   nspin,alat,efermi,evec,
//      >                   jmt,rmt,jws,xstart,
//      >                   vr,vdif,rhotot,xvalws,
//      >                   ztotss,zcorss,
//      >                   numc,nc,lc,kc,ec,
//      >                   header,loc_id,iprint,istop)
// c     ================================================================

extern "C" {
  void single_pot_read(int *present_atom,
                       int *nspin, Real *alat, Real *efermi, Real *evec,
                       int *jmt, Real *rmt, int *jws, Real *xstart,
                       Real *vr, Real *vdif,Real *rhotot,Real *xvalws,
                       Real *ztotss, Real *zcorss,
                       int *numc, int *nc, int *lc,int *kc, Real *ec,
                       char *header,int *loc_id,int *iprint,char *istop);
}

void singlePotentialRead(int loc_id, Atom &a)
{
  int jmt, jws;
  Real rmt, xtart;
  int iprint=0;
  char istop[32]='Not Here';

  single_pot_read(&a.id,
                  &a.nspin, &a.alat, &a.efermi, a.evec,
                  &jmt, &rmt, &jws, &xstart,
                  &a.vr(0,0), &a.vdif, &a.rhotot(0,0), a.xvalws,
                  &a.ztotss, &a.zcorss,
                  &a.numc, &a.nc(0,0), &a.lc(0,0), &a.kc(0,0), &a.ec(0,0),
                  a.header, &loc_id, &iprint, istop);
}

#endif
