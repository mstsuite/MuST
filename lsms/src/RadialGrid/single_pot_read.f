c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine single_pot_read(present_atom,
     >                   nspin,alat,efermi,evec,
     >                   jmt,rmt,jws,xstart,
     >                   vr,vdif,rhotot,xvalws,
     >                   ztotss,zcorss,
     >                   numc,nc,lc,kc,ec,
     >                   header,loc_id,iprint,istop)
c     ================================================================


      use hdf5

      implicit none

      include 'atom_param.h'

      integer(hid_t) loc_id,dti,dtf,dtc
      integer num_read,ns,hdferr
      integer present_atom,nspin,jws,jmt,numc
      real*8 alat,efermi,evec(3),vdif,ztotss,zcorss
      real*8 xstart,rmt
      real*8 vr(iprpts,2), rhotot(iprpts,2), xvalws(2)
      real*8 ec(ipcore,2)
      integer nc(ipcore,2),lc(ipcore,2),kc(ipcore,2)
      character gname*40
      character header*80
      character*32 istop
      character*32 sname
      integer iprint

      parameter (sname='single_pot_read')

      call h5tcopy_f(H5T_NATIVE_INTEGER,dti,hdferr)
      call h5tcopy_f(H5T_NATIVE_DOUBLE,dtf,hdferr)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,dtc,hdferr)

      call read_scalar_int(loc_id,'PresentAtom',
     >                     present_atom,dti,num_read)
      call read_scalar_int(loc_id,'jmt',jmt,dti,num_read)
      call read_scalar_int(loc_id,'jws',jws,dti,num_read)
      call read_scalar_int(loc_id,'Nspin',nspin,dti,num_read)
      call read_scalar_int(loc_id,'NumC',numc,dti,num_read)
      call read_scalar_float(loc_id,'alat',alat,dtf,num_read)
      call read_scalar_float(loc_id,'Efermi',efermi,dtf,num_read)
      call read_scalar_float(loc_id,'Vdif',vdif,dtf,num_read)
      call read_scalar_float(loc_id,'Ztot',ztotss,dtf,num_read)
      call read_scalar_float(loc_id,'Zcore',zcorss,dtf,num_read)
      call read_scalar_float(loc_id,'Xstart',xstart,dtf,num_read)
      call read_scalar_float(loc_id,'rmt',rmt,dtf,num_read)


      call read_vector_float(loc_id,'xvalws',
     >                       xvalws,nspin,dtf,num_read)
          if(nspin.ne.num_read) write(*,*)
     >      'WARNING in single_pot_read, xvalws'

      do ns=1,nspin
          write(gname,'("V",i1.1)') ns
          call read_vector_float(loc_id,gname,vr(1,ns),
     >                           jmt,dtf,num_read)
          if(jmt.ne.num_read) write(*,*) 
     >      'WARNING in single_pot_read, vr'
          write(gname,'("rhotot",i1.1)') ns
          call read_vector_float(loc_id,gname,rhotot(1,ns),
     >                           jws,dtf,num_read)
          if(jws.ne.num_read) 
     >      write(*,*) 'WARNING in single_pot_read, rhotot'
          write(gname,'("ec",i1.1)') ns
!  read core states if numc>0
          if(numc.gt.0) then
            call read_vector_float(loc_id,gname,ec(1,ns),
     >                           numc,dtf,num_read)
            if(numc.ne.num_read)
     >         write(*,*) 'WARNING in single_pot_read, ec'
             write(gname,'("nc",i1.1)') ns
             call read_vector_int(loc_id,gname,nc(1,ns),
     >                            numc,dti,num_read)
            if(numc.ne.num_read) write(*,*)
     >        'WARNING in single_pot_read, nc'
            write(gname,'("lc",i1.1)') ns
            call read_vector_int(loc_id,gname,lc(1,ns),numc,
     >                           dti,num_read)
            if(numc.ne.num_read) 
     >        write(*,*) 'WARNING in single_pot_read, lc'
            write(gname,'("kc",i1.1)') ns
            call read_vector_int(loc_id,gname,kc(1,ns),
     >                           numc,dti,num_read)
            if(numc.ne.num_read)
     >         write(*,*) 'WARNING in single_pot_read, kc'
          end if
      end do
      call read_vector_float(loc_id,'evec',evec,3,dtf,num_read)
          if(3.ne.num_read) 
     >    write(*,*) 'WARNING in single_pot_read, evec'
      call read_vector_char(loc_id,'Header',header,80,dtc,num_read)
          if(80.ne.num_read) 
     >    write(*,*) 'WARNING in single_pot_read, header'

 10   continue

      call h5tclose_f(dti,hdferr)
      call h5tclose_f(dtf,hdferr)
      call h5tclose_f(dtc,hdferr)

      end
!------------------------------------------------------------------
      subroutine read_scalar_int(loc_id,name,value,dti,num_read)

      use hdf5
      implicit none

      integer(hid_t) loc_id,space_id,dset_id,dti
      integer hdferr,value,num_read
      integer(hsize_t) dims(2)
      character(*) name

      dims(1)=1
      dims(2)=0
      num_read=-1
      call h5dopen_f(loc_id,name,dset_id,hdferr)
      if(hdferr.lt.0) return

      call h5dread_f(dset_id,dti,value,dims,hdferr)
      if(hdferr.lt.0) return

      num_read=1

      call h5dclose_f(dset_id,hdferr)

      end
!------------------------------------------------------------------
      subroutine read_scalar_float(loc_id,name,value,dtf,num_read)

      use hdf5
      implicit none

      integer(hid_t) loc_id,space_id,dset_id,dtf
      integer hdferr,num_read
      real*8 value
      integer(hsize_t) dims(2)
      character(*) name

      num_read=-1
      call h5dopen_f(loc_id,name,dset_id,hdferr)
      if(hdferr.lt.0) return

      call h5dread_f(dset_id,dtf,value,dims,hdferr)
      if(hdferr.lt.0) return

      num_read=1

      call h5dclose_f(dset_id,hdferr)

      end
!------------------------------------------------------------------
      subroutine read_vector_int(loc_id,name,value,len,dti,num_read)

      use hdf5
      implicit none

      integer(hid_t) loc_id,space_id,dset_id,dti,dti_file
      integer hdferr,value(len),num_read,len,ndims
      integer(hsize_t) dims(2),maxdims(2)
      character(*) name

      num_read=-1
      call h5dopen_f(loc_id,name,dset_id,hdferr)
      if(hdferr.lt.0) return
      call h5dget_space_f(dset_id,space_id,hdferr)
      call h5sget_simple_extent_ndims_f(space_id,ndims,hdferr)
      if(ndims.ne.1) return
      call h5sget_simple_extent_dims_f(space_id,dims,maxdims,hdferr)

      if(dims(1).gt.len) return

      call h5dread_f(dset_id,dti,value,dims,hdferr)
      if(hdferr.lt.0) return

      num_read=dims(1)

      call h5sclose_f(space_id,hdferr)
      call h5dclose_f(dset_id,hdferr)

      end
!------------------------------------------------------------------
      subroutine read_vector_float(loc_id,name,value,len,dtf,num_read)

      use hdf5
      implicit none

      integer(hid_t) loc_id,space_id,dset_id,dtf,dtf_file
      integer hdferr,num_read,len,ndims
      real*8 value(len)
      integer(hsize_t) dims(2),maxdims(2)
      character(*) name

      num_read=-1
      call h5dopen_f(loc_id,name,dset_id,hdferr)
      if(hdferr.lt.0) return
      call h5dget_space_f(dset_id,space_id,hdferr)
      call h5sget_simple_extent_ndims_f(space_id,ndims,hdferr)
      if(ndims.ne.1) return
      call h5sget_simple_extent_dims_f(space_id,dims,maxdims,hdferr)

      if(dims(1).gt.len) return

      call h5dread_f(dset_id,dtf,value,dims,hdferr)
      if(hdferr.lt.0) return

      num_read=dims(1)

      call h5sclose_f(space_id,hdferr)
      call h5dclose_f(dset_id,hdferr)

      end
!------------------------------------------------------------------
      subroutine read_vector_char(loc_id,name,value,len,dtc,num_read)

      use hdf5
      implicit none

      integer(hid_t) loc_id,space_id,dset_id,dtc
      integer hdferr,num_read,len,ndims
      character(*) value
      integer(hsize_t) dims(2),maxdims(2)
      character(*) name

      num_read=-1
      call h5dopen_f(loc_id,name,dset_id,hdferr)
      if(hdferr.lt.0) return
      call h5dget_space_f(dset_id,space_id,hdferr)
      call h5sget_simple_extent_ndims_f(space_id,ndims,hdferr)
      if(ndims.ne.1) return
      call h5sget_simple_extent_dims_f(space_id,dims,maxdims,hdferr)

      if(dims(1).gt.len) return

      call h5dread_f(dset_id,dtc,value,dims,hdferr)
      if(hdferr.lt.0) return

      num_read=dims(1)

      call h5sclose_f(space_id,hdferr)
      call h5dclose_f(dset_id,hdferr)

      end

