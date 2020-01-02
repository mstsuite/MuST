#include <string.h>
#include "AtomData.hpp"
#include <hdf5.h>
#include "writeSingleAtomData.hpp"
#include "Main/HDF5io.hpp"

int writeSingleAtomData_hdf5(hid_t loc_id, AtomData &atom, int present_atom)
{
  char gname[128];
  
  // printf("WARNING: Writing of hdf5 potentials not implemented yet!\n");
  write_scalar<int>(loc_id,"PresentAtom",present_atom);
  write_scalar<int>(loc_id,"jmt",atom.jmt);
  write_scalar<int>(loc_id,"jws",atom.jws);
  write_scalar<int>(loc_id,"Nspin",atom.nspin);
  write_scalar<int>(loc_id,"NumC",atom.numc);
  write_scalar<double>(loc_id,"alat",atom.alat);
  write_scalar<double>(loc_id,"Efermi",atom.efermi);
  write_scalar<double>(loc_id,"Vdif",atom.vdif);
  write_scalar<double>(loc_id,"Ztot",atom.ztotss);
  write_scalar<double>(loc_id,"Zcore",atom.zcorss);
  write_scalar<double>(loc_id,"Xstart",atom.xstart);
  write_scalar<double>(loc_id,"rmt",atom.rmt);
  // printf("Wrote Scalars\n");
  write_vector<double>(loc_id,"xvalws",&atom.xvalws[0],atom.nspin);
  // printf("Wrote xvalws (%lf %lf)\n",atom.xvalws[0],atom.xvalws[1]);
  
  for(int ns=0; ns<atom.nspin; ns++)
  {
    snprintf(gname,100,"V%1d",ns+1);
    //printf("vr(0,0)=%lf  vr(jmt-1,0)=%lf\n",atom.vr(0,ns),atom.vr(atom.jmt-1,ns));
    write_vector<double>(loc_id,gname,&atom.vr(0,ns),atom.jmt);
    // printf("Wrote %s\n",gname);
    //return 0;
    snprintf(gname,100,"rhotot%1d",ns+1);
    write_vector<double>(loc_id,gname,&atom.rhotot(0,ns),atom.jws);
    // printf("Wrote %s\n",gname);
    if(atom.numc>0) // write core states if they exist
    {
      snprintf(gname,100,"ec%1d",ns+1);
      write_vector<double>(loc_id,gname,&atom.ec(0,ns),atom.numc);
      //printf("Wrote %s\n",gname);
      snprintf(gname,100,"nc%1d",ns+1);
      write_vector<int>(loc_id,gname,&atom.nc(0,ns),atom.numc);
      //printf("Wrote %s\n",gname);
      snprintf(gname,100,"lc%1d",ns+1);
      write_vector<int>(loc_id,gname,&atom.lc(0,ns),atom.numc);
      //printf("Wrote %s\n",gname);
      snprintf(gname,100,"kc%1d",ns+1);
      write_vector<int>(loc_id,gname,&atom.kc(0,ns),atom.numc);
      //printf("Wrote %s\n",gname);
    }
  }
    write_vector<char>(loc_id,"Header",&atom.header[0],80);
    // printf("Wrote Header\n");
    for(int j=0; j<100; j++) gname[j]=' ';
    write_vector<char>(loc_id,"JTitle",&gname[0],80);
    //printf("Wrote JTitle\n",gname);
    write_vector<double>(loc_id,"evec",&atom.evec[0],3);
    //printf("Wrote evec\n",gname);

    // printf("Wrote present_atom=%d\n",present_atom);
  
  return 0;
}
