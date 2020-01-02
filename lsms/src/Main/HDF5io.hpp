/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_HDF5_IO
#define LSMS_HDF5_IO

// #define H5_USE_16_API
#include <hdf5.h>
#include <hdf5_hl.h>

//------------------------------------------------------------------
template<typename T>
int read_scalar(hid_t loc_id, const char *name, T &value, hid_t dti)
{
  hid_t space_id,dset_id;
  int hdferr,num_read;

  // printf("Attempting to read scalar '%s'\n",name);

  num_read=-1;
  dset_id=H5Dopen2(loc_id,name,H5P_DEFAULT);

  hdferr= H5Dread(dset_id,dti,H5S_ALL,H5S_ALL,H5P_DEFAULT,&value);
  if(hdferr<0) return hdferr;

  num_read=1;

  H5Dclose(dset_id);

  return num_read;
}

template<typename T>
int read_scalar(hid_t loc_id, const char *name, T &value)
{
  hid_t space_id,dset_id;
  int hdferr,num_read;

  // printf("Attempting to read scalar '%s'\n",name);

  num_read=-1;
  dset_id=H5Dopen2(loc_id,name,H5P_DEFAULT);

  hdferr= H5Dread(dset_id,TypeTraits<T>::hdf5Type(),H5S_ALL,H5S_ALL,H5P_DEFAULT,&value);
  if(hdferr<0) return hdferr;

  num_read=1;

  H5Dclose(dset_id);

  return num_read;
}


//------------------------------------------------------------------
template<typename T>
int read_vector(hid_t loc_id,const char *name,T *value,int len,hid_t dti)
{
  hid_t space_id,dset_id,dti_file;
  int hdferr,ndims;
  hsize_t dims[2], maxdims[2];

  // printf("Attempting to read vector '%s'\n",name);

  dset_id= H5Dopen2(loc_id,name,H5P_DEFAULT);

  space_id= H5Dget_space(dset_id);
  ndims= H5Sget_simple_extent_ndims(space_id);
  if(ndims!=1) return -1;
  H5Sget_simple_extent_dims(space_id,dims,maxdims);

  if(dims[0]>len) return -1;

  hdferr= H5Dread(dset_id,dti,H5S_ALL,H5S_ALL,H5P_DEFAULT,value);
  if(hdferr<0) return -1;

  H5Sclose(space_id);
  H5Dclose(dset_id);

  return dims[0];
}

template<typename T>
int read_vector(hid_t loc_id,const char *name,T *value,int len)
{
  hid_t space_id,dset_id,dti_file;
  int hdferr,ndims;
  hsize_t dims[2], maxdims[2];

  // printf("Attempting to read vector '%s'\n",name);

  dset_id= H5Dopen2(loc_id,name,H5P_DEFAULT);

  space_id= H5Dget_space(dset_id);
  ndims= H5Sget_simple_extent_ndims(space_id);
  if(ndims!=1) return -1;
  H5Sget_simple_extent_dims(space_id,dims,maxdims);

  if(dims[0]>len) return -1;

  hdferr= H5Dread(dset_id,TypeTraits<T>::hdf5Type(),H5S_ALL,H5S_ALL,H5P_DEFAULT,value);
  if(hdferr<0) return -1;

  H5Sclose(space_id);
  H5Dclose(dset_id);

  return dims[0];
}

template<typename T>
int write_scalar(hid_t loc_id, const char *name, T value)
{
  hid_t space_id, dset_id;
  hsize_t dims[]={1};

  // printf("Attempting to write scalar '%s'\n",name);

  space_id = H5Screate(H5S_SCALAR);
  dset_id = H5Dcreate2(loc_id, name, TypeTraits<T>::hdf5Type(), space_id,
                       H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  
  H5Dwrite(dset_id, TypeTraits<T>::hdf5Type(),H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);

  H5Dclose(dset_id);
  H5Sclose(space_id);
  return 1;
}

template<typename T>
int write_vector(hid_t loc_id, const char *name, T *value, int len)
{
  hid_t space_id, dset_id;
  hsize_t dims[2];

  dims[0]=len;
  space_id = H5Screate_simple(1, dims, NULL);
  dset_id = H5Dcreate2(loc_id, name, TypeTraits<T>::hdf5Type(), space_id,
                       H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Dwrite(dset_id,TypeTraits<T>::hdf5Type(),H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
  H5Dclose(dset_id);
  H5Sclose(space_id);

  return len;
}

#endif
