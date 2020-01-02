#ifndef WRITESINGLEATOMDATA_H
#define WRITESINGLEATOMDATA_H

#include <hdf5.h>
#include "AtomData.hpp"

int writeSingleAtomData_hdf5(hid_t loc_id, AtomData &atom, int present_atom);

int writeSingleAtomData_bigcell(const char *fname, AtomData &atom);

extern "C"
{
void f_writesingleatomdata_bigcell_(const char *fname_c,int *fname_l,
                                    char *header_c, int *header_l,int *jmt,int *jws, double *xstart,
                                   double *rmt,double *alat,double *efermi,
                                   double *vdif,double *ztotss,double *zcorss,
                                   int *nspin,int *numc,double *xvalws,
                                   double *vr,double *rhotot,double *corden,int *v_dim,
                                   char *atname, double *zsemss, double *zvalss,
                                   double *ec,int *nc,int *lc,int *kc,int *c_dim,int);
};
#endif
