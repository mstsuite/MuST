#include <string.h>
#include "AtomData.hpp"
#include "writeSingleAtomData.hpp"

int writeSingleAtomData_bigcell(const char *fname, AtomData &atom)
{
  int fname_l,v_dim,c_dim,header_l;
  char atname[4]="xx";

  fname_l=strlen(fname);
  header_l=80;
  v_dim=atom.vr.l_dim();
  c_dim=atom.ec.l_dim();

  f_writesingleatomdata_bigcell_(fname,&fname_l,
                                 atom.header,&header_l,&atom.jmt,&atom.jws,&atom.xstart,
                                 &atom.rmt,&atom.alat,&atom.efermi,&atom.vdif,&atom.ztotss,&atom.zcorss,
                                 &atom.nspin,&atom.numc,atom.xvalws,
                                 &atom.vr(0,0),&atom.rhotot(0,0),&atom.corden(0,0), &v_dim,
                                 atname,&atom.zsemss,&atom.zvalss,
                                 &atom.ec(0,0),&atom.nc(0,0),&atom.lc(0,0),&atom.kc(0,0),&c_dim,2);
  return -1;
}
