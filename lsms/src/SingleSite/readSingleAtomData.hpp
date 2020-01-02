#ifndef READSINGLEATOMDATA_H
#define READSINGLEATOMDATA_H

#include <hdf5.h>
#include "AtomData.hpp"

int readSingleAtomData_hdf5(hid_t loc_id, AtomData &atom);

int readSingleAtomData_bigcell(const char *fname, AtomData &atom);

/*
      subroutine readSingleAtomData_bigcell(fname_c,fname_l,
     >     header,jmt,jws,xstart,
     >     rmt,alat,efermi,vdif,ztotss,zcorss,
     >     nspin,numc,xvalws,
     >     vr,rhotot,v_dim,
     >     ec,nc,lc,kc,c_dim)
      implicit none

      byte fname_c(128)
      integer fname_l,i
      integer jmt,jws,nspin,numc,v_dim,c_dim
      integer nc(c_dim,*),lc(c_dim,*),kc(c_dim,*)
      real*8 xvalws(2)
      real*8 xstart,rmt,alat,efermi,vdif,ztotss,zcorss
      real*8 vr(v_dim,*),rhotot(v_dim,*)
      real*8 ec(c_dim,*)
      character*80 header,jtitle
      character*127 fname
      integer funit
*/

extern "C"
{
void f_readsingleatomdata_bigcell_(const char *fname_c,int *fname_l,
                                   char *header,int *jmt,int *jws, double *xstart,
                                   double *rmt,double *alat,double *efermi,
                                   double *vdif,double *ztotss,double *zcorss,
                                   int *nspin,int *numc,double *xvalws,
                                   double *vr,double *rhotot,double *corden,int *v_dim,
                                   double *ec,int *nc,int *lc,int *kc,int *c_dim,int);
};
#endif
