module HDF5Module
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
   use hdf5
!
public :: initHDF5,     &
          endHDF5,      &
          hdf5_fcreate, &
          hdf5_fopen,   &
          hdf5_fclose,  &
          hdf5_gcreate, &
          hdf5_gopen,   &
          hdf5_gclose,  & 
          hdf5_read,    &
          hdf5_write
!
   interface hdf5_write
      module procedure hdf5_write_s
      module procedure hdf5_write_i0, hdf5_write_i1 ! , hdf5_write_i2
      module procedure hdf5_write_r0, hdf5_write_r1 ! , hdf5_write_r2
   end interface hdf5_write
!
   interface hdf5_read
      module procedure hdf5_read_s
      module procedure hdf5_read_i0, hdf5_read_i1 ! , hdf5_read_i2
      module procedure hdf5_read_r0, hdf5_read_r1 ! , hdf5_read_i2
   end interface hdf5_read
!
private
!
   logical :: Initialized
!
   integer(kind=IntKind), parameter :: str_len = 150
!
   type path_struct
      character(len=str_len) :: path
      integer(kind=IntKind) :: key
      type(path_struct), pointer :: next
   end type path_struct
!
   type(path_struct), pointer :: dir_list
!
   type hdf5_id_list
      integer (kind=IntKind)       :: key
      integer (kind=hid_t)         :: h5d_id
      type (hdf5_id_list), pointer :: next
   end type hdf5_id_list
!
   type(hdf5_id_list), pointer :: hdf5_list
!
   logical :: fopen_flag
   logical :: gopen_flag
!
   integer (kind=hid_t) :: dti, dtr, dts
   integer (kind=hid_t) :: dti_file, dtr_file, dts_file
   integer (kind=hid_t) :: current_file_id, current_group_id
!
   integer (kind=IntKind) :: NumActivePathKey
   integer (kind=IntKind) :: NumActiveIDs
   integer (kind=IntKind) :: NumGenID
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initHDF5()  
!  ===================================================================
   implicit none
!
   integer(kind=IntKind) :: hdferr
!
   call h5open_f(hdferr)
!
   call h5tcopy_f( H5T_NATIVE_INTEGER, dti, hdferr )
   call h5tcopy_f( H5T_NATIVE_DOUBLE, dtr, hdferr )
   call h5tcopy_f( H5T_NATIVE_CHARACTER, dts, hdferr )
!
!  The following lines should be a separate routine which specify 
!  the encription type of the data in the file, or 
!  some variables should be passed in the initHDF5 to fix them. 
!
   call h5tcopy_f( H5T_IEEE_F64BE, dtr_file, hdferr )
   call h5tcopy_f( H5T_STD_I32BE, dti_file, hdferr )
   call h5tcopy_f( H5T_NATIVE_CHARACTER, dts_file, hdferr )
!
   nullify( hdf5_list )
   NumActivePathKey = 0
   NumActiveIDs = 0
   NumGenID = 0
   fopen_flag = .false.
   gopen_flag = .false.
!
   Initialized = .true.
!
   end subroutine initHDF5
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endHDF5
!  ===================================================================
   implicit none
!
   integer(kind=IntKind) :: hdferr
!   
   call h5tclose_f(dti,hdferr)
   call h5tclose_f(dtr,hdferr)
   call h5tclose_f(dts,hdferr)
!
   call h5tclose_f(dti_file,hdferr)
   call h5tclose_f(dtr_file,hdferr)
   call h5tclose_f(dts_file,hdferr)
!
   call delHDF5()
   NumGenID = 0
!
   call h5close_f(hdferr)
!
   Initialized = .false.
!
   end subroutine endHDF5
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delHDF5(key)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), optional :: key
!
   type (hdf5_id_list), pointer :: phdf5
   type (hdf5_id_list), pointer :: phdf5_id
!
   integer(kind=IntKind) :: i   
!
   phdf5_id => hdf5_list
   if ( present(key) ) then
      phdf5 => hdf5_list
      do i = 1, NumActiveIDs
         if ( key /= phdf5_id%key ) then
            phdf5 => phdf5_id
            phdf5_id => phdf5_id%next
         else
            if ( associated(phdf5_id%next) ) then
               phdf5%next => phdf5_id%next
               deallocate(phdf5_id)
            else
               nullify(phdf5%next)
               deallocate(phdf5_id)
            endif
            NumActiveIDs = NumActiveIDs -1
            return
         endif
      enddo
      call ErrorHandler("endHDF5",'open key not found', key)
   else
      do i = 1, NumActiveIDs
         phdf5 => phdf5_id%next
         deallocate(phdf5_id)
         phdf5_id => phdf5
      enddo
      NumActiveIDs = 0
   endif
!
   end subroutine delHDF5
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getHDF5_id(key)                            result(h5d_id)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind) :: i
   integer(kind=hid_t)   :: h5d_id
!
   type (hdf5_id_list), pointer :: phdf5_l
!   
   phdf5_l => hdf5_list
   do i = 1, NumActiveIDs
      if ( key == phdf5_l%key ) then
         h5d_id = phdf5_l%h5d_id
         exit
      endif
      phdf5_l => phdf5_l%next
   enddo
!
   end function getHDF5_id
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine addHDF5_id(key,h5d_id) 
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(out) :: key
   integer(kind=hid_t), intent(in) :: h5d_id
!
   integer(kind=IntKind) :: i
!
   type(hdf5_id_list), pointer :: phdf5_id
!   
   key = NumGenID+1
   if ( NumActiveIDs == 0 ) then
      allocate( hdf5_list )
      hdf5_list%key = key
      hdf5_list%h5d_id = h5d_id
      nullify( hdf5_list%next )
      NumActiveIDs = NumActiveIDs+1
   else
!
      phdf5_id => hdf5_list
      do i = 1,NumActiveIDs-1
         phdf5_id => phdf5_id%next
      enddo
      allocate( phdf5_id%next )
      phdf5_id => phdf5_id%next
      phdf5_id%key = key
      phdf5_id%h5d_id = h5d_id
      nullify(phdf5_id%next)
      NumActiveIDs = NumActiveIDs+1
   endif
!
   end subroutine addHDF5_id
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine add_path( dir_name, key, parent_key )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: dir_name
!
   integer(kind=IntKind), optional :: parent_key
!
   integer(kind=IntKind), intent(out) :: key
!
   character(len=str_len) :: p_path   
!
   integer(kind=IntKind) :: i, str_size
!
   type(path_struct), pointer :: p_dir
!
   if ( NumActivePathKey==0 ) then
      allocate(dir_list)
      nullify(dir_list%next)
      call make_dir(trim(dir_name))
       
      dir_list%path = trim(dir_name)//'/'
      dir_list%key = NumActivePathKey+1     
   else
      p_dir => dir_list
      do i = 1,NumActivePathKey-1
         p_dir => p_dir%next
      enddo
      allocate(p_dir%next)
      p_dir => p_dir%next
      nullify(p_dir%next)
      if ( present(parent_key) ) then
         p_path = get_path(parent_key)
         dir_list%path = trim(p_path)//trim(dir_name)
         call make_dir(p_path//trim(dir_name))
!         call make_dir(dir_list%path)
      else
         dir_list%path = trim(dir_name)
         call make_dir(trim(dir_name))
!         call make_dir(dir_list%path)
      endif
      dir_list%path = trim(dir_list%path)//'/'
   endif
!
   NumActivePathKey = NumActivePathKey+1 
!
   end subroutine add_path
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine del_path_list(key)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), optional :: key
!   
   integer(kind=IntKind) :: i
!
   type(path_struct), pointer :: p_dir, p_dir_tmp
!
   p_dir => dir_list
   if ( .not.present(key) ) then
      do i = 1,NumActivePathKey
         p_dir_tmp => p_dir%next
         deallocate(p_dir)
         p_dir => p_dir_tmp
      enddo
      NumActivePathKey = 0
   else
      if ( NumActivePathKey==1 ) then
         deallocate(p_dir)
         NumActivePathKey = 0
         return
      endif
      do i = 1,NumActivePathKey-1
         p_dir_tmp => p_dir
         p_dir => p_dir%next
         if ( key == p_dir%key ) then
            p_dir_tmp%next => p_dir%next
            deallocate(p_dir)
         endif
      enddo
      NumActivePathKey = NumActivePathKey-1
   endif
!
   end subroutine del_path_list
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function get_path(key)                               result(p_path)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: key
!   
   character(len=str_len) :: p_path
!
   integer(kind=IntKind) :: i
!
   type(path_struct), pointer :: p_dir
!
   p_dir => dir_list
   do i = 1,NumActivePathKey
      if ( p_dir%key == key ) then
         p_path = p_dir%path
         exit
      endif
   enddo
!
   end function get_path
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function get_pathkey(name)                           result(key)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind) :: key
!   
   character(len=*), intent(in) :: name
!
   character(len=str_len) :: path
   type(path_struct), pointer :: p_dir
!
   integer(kind=IntKind) :: i, str_size
!
   str_size = len( trim(name) )+1
   p_dir => dir_list
   do i = 1,NumActivePathKey
      path = p_dir%path
      if ( trim(path) == trim(name)//'/' ) then
         key = p_dir%key
         exit
      endif
   enddo
!
   end function get_pathkey
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine make_dir(name)                              
!  ===================================================================
   implicit none
!
   character(len=*), target :: name
   integer(kind=IntKind) :: access, status
!
   status=access('name',' ')
   if ( status==0 ) then
      status = access('name','rx')
      if (status/=0) then 
         call system('mkdir '//name)
      endif 
   else
      call system('mkdir '//name)
   endif
!
   end subroutine make_dir
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine parent_create(path_name,key)
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: path_name
!
   integer(kind=IntKind), intent(out) :: key
!   
   call add_path( path_name, key )
!
   end subroutine parent_create
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine child_create(parent_key,child_name,key)
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: child_name
!
   integer(kind=IntKind), intent(in)  :: parent_key
   integer(kind=IntKind), intent(out) :: key
!  
   call add_path( trim(child_name), key, parent_key )
!
   end subroutine child_create
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_fcreate(path_key,file_name,key)
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: file_name
!
   integer(kind=IntKind), optional    :: path_key
   integer(kind=IntKind), intent(out) :: key
!
   character(len=str_len) :: path
!
   integer(kind=hid_t) :: file_id
!
   integer(kind=IntKind) :: hdferr
!   
   if ( present(path_key) ) then
      path = get_path(path_key)
      call h5fopen_f( trim(path)//trim(file_name), H5F_ACC_RDONLY_F,  & 
                   file_id, hdferr )
   else
      call h5fopen_f( trim(file_name), H5F_ACC_RDONLY_F,              & 
                      file_id, hdferr )
   endif
   call addHDF5_id(key,file_id)
!
   end subroutine hdf5_fcreate
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_fopen(path_key,file_name,key)
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: file_name
!
   integer(kind=IntKind), optional    :: path_key
   integer(kind=IntKind), intent(out) :: key
!
   character(len=str_len) :: path
!
   integer(kind=hid_t) :: file_id
!
   integer(kind=IntKind) :: hdferr
!   
   if ( present(path_key) ) then
      path = get_path(path_key)
      call h5fopen_f( trim(path)//trim(file_name), H5F_ACC_RDONLY_F,  & 
                   file_id, hdferr )
   else
      call h5fopen_f( trim(file_name), H5F_ACC_RDONLY_F,              & 
                      file_id, hdferr )
   endif
   call addHDF5_id(key,file_id)
!
   end subroutine hdf5_fopen
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_fclose(key)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=hid_t)   :: file_id
!
   integer(kind=IntKind) :: hdferr
!   
   file_id = getHDF5_id(key)
   call h5fclose_f( file_id, hdferr )
   call delHDF5(key)
!
   end subroutine hdf5_fclose
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_gcreate(key,cname,child_key)
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: cname
!
   integer(kind=IntKind), intent(in) :: key
   integer(kind=IntKind), intent(out) :: child_key
!
   integer(kind=hid_t) :: h5d_id, h5d_cid
!
   integer(kind=IntKind) :: hdferr
!   
   h5d_id = getHDF5_id(key)
   call h5gcreate_f( h5d_id, cname, h5d_cid, hdferr )
   call addHDF5_id(child_key,h5d_cid)
!
   end subroutine hdf5_gcreate
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_gopen(key,cname,child_key)
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: cname
!
   integer(kind=IntKind), intent(in)  :: key
   integer(kind=IntKind), intent(out) :: child_key 
!
   integer(kind=hid_t) :: h5d_id, h5d_cid
!
   integer(kind=IntKind) :: hdferr
! 
   h5d_id = getHDF5_id(key)
   call h5gopen_f( h5d_id, cname, h5d_cid, hdferr )
   call addHDF5_id( child_key, h5d_cid )
!
   end subroutine hdf5_gopen
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_gclose(key)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=hid_t) :: h5d_id
!
   integer(kind=IntKind) :: hdferr
!
   h5d_id = getHDF5_id(key)
   call h5gclose_f( h5d_id, hdferr )
   call delHDF5(key)
!
   end subroutine hdf5_gclose
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_read_s( key, name, value, vlen, num_read )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in)  :: name
   character(len=*), target :: value
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), intent(in) :: vlen
   integer(kind=IntKind), intent(out) :: num_read
!
   integer(kind=hid_t)   :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2), maxdims(2)
   integer(kind=IntKind) :: hdferr, ndims
!
   num_read = -1
   h5d_id = getHDF5_id(key)
   call h5dopen_f( h5d_id, name, dset_id, hdferr )
   if ( hdferr < 0 ) then
      return
   endif
!
   call h5dget_space_f( dset_id, space_id, hdferr )
   call h5sget_simple_extent_ndims_f( space_id, ndims, hdferr )
   if (  ndims /= 1 ) then
      return
   endif
!
   call h5sget_simple_extent_dims_f( space_id, dims, maxdims, hdferr )
   if ( dims(1) > vlen ) then
      return
   endif
!
   call h5dread_f( dset_id, dts, value, dims, hdferr )
   if ( hdferr<0 ) then
      return
   endif
!
   num_read=dims(1)
!
   call h5sclose_f( space_id, hdferr )
   call h5dclose_f( dset_id, hdferr )
!
   end subroutine hdf5_read_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_read_i0( key, name, value, num_read )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: name
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), intent(out) :: value
   integer(kind=IntKind), intent(inout) :: num_read
!
   integer(kind=hid_t)        :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2)
   integer(kind=IntKind) :: hdferr
!
   dims(1)=1
   dims(2)=0
   num_read=-1
   h5d_id = getHDF5_id(key)
   call h5dopen_f( h5d_id, name, dset_id, hdferr )
   if ( hdferr < 0 ) then
      return
   endif
!
   call h5dread_f( dset_id, dti, value, dims, hdferr )
   if ( hdferr<0 ) then
      return
   endif
!
   num_read = 1
!
   call h5dclose_f( dset_id, hdferr )
!
   end subroutine hdf5_read_i0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_read_i1( key, name, value, vlen, num_read )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: name
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), target :: value(:)
   integer(kind=IntKind), intent(in) :: vlen
   integer(kind=IntKind), intent(inout) :: num_read
!
   integer(kind=hid_t)   :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2), maxdims(2)
   integer(kind=IntKind) :: hdferr, ndims
!
   num_read=-1
   h5d_id = getHDF5_id(key)
   call h5dopen_f( h5d_id, name, dset_id, hdferr )
   if ( hdferr < 0 ) then
      return
   endif
!
   call h5dget_space_f(dset_id,space_id,hdferr)
   call h5sget_simple_extent_ndims_f( space_id, ndims, hdferr )
   if (  ndims /= 1 ) then
      return
   endif
!
   call h5sget_simple_extent_dims_f(space_id,dims,maxdims,hdferr)
   if ( dims(1) > vlen ) then
      return
   endif
!
   call h5dread_f( dset_id, dti, value, dims, hdferr )
   if ( hdferr < 0 ) then
      return
   endif
!
   num_read = dims(1)
!
   call h5sclose_f( space_id, hdferr)
   call h5dclose_f( dset_id, hdferr )
!
   end subroutine hdf5_read_i1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_read_r0( key, name, value, num_read )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: name
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), intent(inout) :: num_read
!
   real(kind=RealKind), intent(out) :: value
!
   integer(kind=hid_t)        :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2)
   integer(kind=IntKind) :: hdferr
!
   dims(1)=1
   dims(2)=0
   num_read=-1
   h5d_id = getHDF5_id(key)
   call h5dopen_f( h5d_id, name, dset_id, hdferr )
   if ( hdferr < 0 ) then
      return
   endif
!
   call h5dread_f( dset_id, dtr, value, dims, hdferr )
   if ( hdferr < 0 ) then
      return
   endif
!
   num_read = 1
!
   call h5dclose_f( dset_id, hdferr )
!
   end subroutine hdf5_read_r0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_read_r1( key, name, value, vlen, num_read )
!  ===================================================================
   implicit none
!
!
   character(len=*), intent(in) :: name
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), intent(in) :: vlen
   integer(kind=IntKind), intent(inout) :: num_read
!
   real(kind=RealKind), target :: value(:)
!
   integer(kind=hid_t)        :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2), maxdims(2)
   integer(kind=IntKind) :: hdferr, ndims
!
   num_read=-1
   h5d_id = getHDF5_id(key) 
   call h5dopen_f( h5d_id, name, dset_id, hdferr )
   if ( hdferr < 0 ) then
      return
   endif
!
   call h5dget_space_f( dset_id, space_id, hdferr)
   call h5sget_simple_extent_ndims_f( space_id, ndims, hdferr )
   if (  ndims /= 1 ) then
      return
   endif
!
   call h5sget_simple_extent_dims_f(space_id,dims,maxdims,hdferr)
   if ( dims(1) > vlen ) then
      return
   endif
!
   call h5dread_f( dset_id, dtr, value, dims, hdferr )
   if ( hdferr < 0 ) then
      return
   endif
!
   num_read=dims(1)
!
   call h5sclose_f( space_id, hdferr )
   call h5dclose_f( dset_id, hdferr )
!
   end subroutine hdf5_read_r1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_write_s( key, name, value, vlen )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: name
   character(len=*), intent(in) :: value
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), intent(in) :: vlen
!
   integer(kind=hid_t)        :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2), maxdims(2)
   integer(kind=IntKind) :: hdferr, ndims
!
   dims(1) = vlen
   call h5screate_simple_f( 1, dims, space_id, hdferr )
!
   h5d_id = getHDF5_id(key)
   call h5dcreate_f( h5d_id, name, dts_file, space_id, dset_id, hdferr )
!
   call h5dwrite_f( dset_id, dts, value, dims, hdferr )
!
   call h5dclose_f( dset_id, hdferr )
   call h5sclose_f( space_id, hdferr )
!
   end subroutine hdf5_write_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_write_i0( key, name, value )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: name
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), intent(in) :: value
!
   integer(kind=hid_t)        :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2)
   integer(kind=IntKind) :: hdferr
!
   call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
!   
   h5d_id = getHDF5_id(key)
   call h5dcreate_f( h5d_id, name, dtr_file, space_id, dset_id, hdferr )
!
   call h5dwrite_f(dset_id,dti,value,dims,hdferr)
!   
   call h5dclose_f(dset_id,hdferr)
   call h5sclose_f(space_id,hdferr)
!
   end subroutine hdf5_write_i0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_write_i1( key, name, value, vlen )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: name
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), intent(in) :: vlen
   integer(kind=IntKind), target :: value(:)
!
   integer(kind=hid_t)        :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2)
   integer(kind=IntKind) :: hdferr
!
   dims(1)=vlen
   call h5screate_simple_f( 1, dims, space_id, hdferr )
!   
   h5d_id = getHDF5_id(key)
   call h5dcreate_f( h5d_id, name, dtr_file, space_id, dset_id, hdferr )
!
   call h5dwrite_f( dset_id, dti, value, dims, hdferr )
!   
   call h5dclose_f( dset_id, hdferr )
   call h5sclose_f( space_id, hdferr )
!
   end subroutine hdf5_write_i1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_write_r0( key, name, value )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: name
!
   integer(kind=IntKind), intent(in) :: key
!
   real(kind=RealKind), intent(in) :: value
!
   integer(kind=hid_t)        :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2)
   integer(kind=IntKind) :: hdferr
!
   call h5screate_f( H5S_SCALAR_F, space_id, hdferr)
!   
   h5d_id = getHDF5_id(key)
   call h5dcreate_f( h5d_id, name, dtr_file, space_id, dset_id, hdferr )
!
   call h5dwrite_f( dset_id, dtr, value, dims, hdferr )
!
   call h5dclose_f( dset_id, hdferr)
   call h5sclose_f( space_id, hdferr)
!
   end subroutine hdf5_write_r0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hdf5_write_r1( key, name, value, vlen )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: name
!
   integer(kind=IntKind), intent(in) :: key
!
   integer(kind=IntKind), intent(in) :: vlen
   real(kind=RealKind), target :: value(:)
!
   integer(kind=hid_t)        :: h5d_id, space_id, dset_id
   integer(hsize_t)      :: dims(2)
   integer(kind=IntKind) :: hdferr
!
   dims(1) = vlen
   call h5screate_simple_f( 1, dims, space_id, hdferr )
!   
   h5d_id = getHDF5_id(key)
   call h5dcreate_f( h5d_id, name, dtr_file, space_id, dset_id, hdferr )
!
   call h5dwrite_f( dset_id, dtr, value, dims, hdferr )
!   
   call h5dclose_f( dset_id, hdferr )
   call h5sclose_f( space_id, hdferr )
!
   end subroutine hdf5_write_r1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine copyStrToArray( string, length, sarray, asize, ind )
!  ===================================================================
   implicit none
!
   character(len=*), intent(in) :: string
   character(len=1), target     :: sarray(:)
!
   integer(kind=IntKind), intent(in) :: length, asize, ind
!
   integer(kind=IntKind) :: i
!
   do i = 1,length
      sarray(ind+i) = string(i:i)
   enddo
!
   end subroutine copyStrToArray
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine copyArrayToStr( sarray, asize, string, length, ind )
!  ===================================================================
   implicit none
!
   character(len=*), intent(out) :: string
   character(len=1), target     :: sarray(:)
!
   integer(kind=IntKind), intent(in) :: length, asize, ind
!
   integer(kind=IntKind) :: i
!
   do i = 1,asize
      string(ind+i:ind+i) = sarray(i)
   enddo
!
   end subroutine copyArrayToStr
!  ==================================================================

end module HDF5Module
