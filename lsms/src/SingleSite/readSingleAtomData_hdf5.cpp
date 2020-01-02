#define H5_USE_16_API
#include <hdf5.h>
#include "AtomData.hpp"
#include "Main/HDF5io.hpp"

//--------------------------------------------------------------------------
int readSingleAtomData_hdf5(hid_t loc_id, AtomData &atom)
{

//      include 'atom_param.h'

  hid_t dti,dtf,dtc;

  int present_atom;
  int num_read;
  char gname[40];

  dti = H5Tcopy(H5T_NATIVE_INT);
  dtf = H5Tcopy(H5T_NATIVE_DOUBLE);
  dtc = H5Tcopy(H5T_NATIVE_CHAR);

  read_scalar<int>(loc_id,"PresentAtom",present_atom,dti);
  read_scalar<int>(loc_id,"jmt",atom.jmt,dti);
  read_scalar<int>(loc_id,"jws",atom.jws,dti);
  read_scalar<int>(loc_id,"Nspin",atom.nspin,dti);
  read_scalar<int>(loc_id,"NumC",atom.numc,dti);
  read_scalar<double>(loc_id,"alat",atom.alat,dtf);
  read_scalar<double>(loc_id,"Efermi",atom.efermi,dtf);
  read_scalar<double>(loc_id,"Vdif",atom.vdif,dtf);
  read_scalar<double>(loc_id,"Ztot",atom.ztotss,dtf);
  read_scalar<double>(loc_id,"Zcore",atom.zcorss,dtf);
  read_scalar<double>(loc_id,"Xstart",atom.xstart,dtf);
  read_scalar<double>(loc_id,"rmt",atom.rmt,dtf);

  int npts=std::max(atom.jmt,atom.jws);
  if(npts>atom.vr.n_row())
  {
    printf("No. of radial grid point in potential file (jmt, jws) is larger then lsms.global.iprpts!\n");
    exit(1);
    //  atom.resizePotential(npts);
  }
  // if(atom.numc<atom.nc.n_row())
  atom.resizeCore(atom.numc);

  num_read=read_vector<double>(loc_id,"xvalws",atom.xvalws,atom.nspin,dtf);
  if(atom.nspin!=num_read) printf("WARNING in single_pot_read, xvalws\n");

  for(int ns=0; ns<atom.nspin; ns++)
  {
    snprintf(gname,40,"V%1.1d", ns+1);
    num_read=read_vector<double>(loc_id,gname,&atom.vr(0,ns),atom.jmt,dtf);
    if(atom.jmt!=num_read) printf("WARNING in single_pot_read, vr\n");
    snprintf(gname,40,"rhotot%1.1d", ns+1);
    num_read=read_vector<double>(loc_id,gname,&atom.rhotot(0,ns),atom.jws,dtf);
    if(atom.jws!=num_read) printf("WARNING in single_pot_read, rhotot\n");

//  read core states if numc>0
    if(atom.numc>0)
    {
      snprintf(gname,40,"ec%1.1d", ns+1);
      num_read=read_vector<double>(loc_id,gname,&atom.ec(0,ns),atom.numc,dtf);
      if(atom.numc!=num_read) printf("WARNING in single_pot_read, ec\n");
      snprintf(gname,40,"nc%1.1d", ns+1);
      num_read=read_vector<int>(loc_id,gname,&atom.nc(0,ns),atom.numc,dti);
      if(atom.numc!=num_read) printf("WARNING in single_pot_read, nc\n");
      snprintf(gname,40,"lc%1.1d", ns+1);
      num_read=read_vector<int>(loc_id,gname,&atom.lc(0,ns),atom.numc,dti);
      if(atom.numc!=num_read) printf("WARNING in single_pot_read, lc\n");
      snprintf(gname,40,"kc%1.1d", ns+1);
      num_read=read_vector<int>(loc_id,gname,&atom.kc(0,ns),atom.numc,dti);
      if(atom.numc!=num_read) printf("WARNING in single_pot_read, kc\n");
    }
  }
  num_read=read_vector<double>(loc_id,"evec",atom.evec,3,dtf);
  if(3!=num_read) printf("WARNING in single_pot_read, evec\n");
  num_read=read_vector<char>(loc_id,"Header",atom.header,80,dtc);
  if(80!=num_read) printf("WARNING in single_pot_read, header\n");

  H5Tclose(dti);
  H5Tclose(dtf);
  H5Tclose(dtc);

  return present_atom;
}


