#include <string.h>
#include "AtomData.hpp"
#include "readSingleAtomData.hpp"

int readSingleAtomData_bigcell(const char *fname, AtomData &atom)
{
  int fname_l,v_dim,c_dim;

  fname_l=strlen(fname);
  v_dim=atom.vr.l_dim();
  c_dim=atom.ec.l_dim();

  f_readsingleatomdata_bigcell_(fname,&fname_l,
                               atom.header,&atom.jmt,&atom.jws,&atom.xstart,
                               &atom.rmt,&atom.alat,&atom.efermi,&atom.vdif,&atom.ztotss,&atom.zcorss,
                               &atom.nspin,&atom.numc,atom.xvalws,
                                &atom.vr(0,0),&atom.rhotot(0,0),&atom.corden(0,0), &v_dim,
                                &atom.ec(0,0),&atom.nc(0,0),&atom.lc(0,0),&atom.kc(0,0),&c_dim,80);
  atom.evec[0]=atom.evec[1]=0.0; atom.evec[2]=1.0;
  return -1;
}
