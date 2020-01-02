!  *******************************************************************
!  *                                                                 *
!  *                                                                 *
!  *                                                                 *
!  *                                                                 *
!  *******************************************************************
module SpinRotationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, HALF, ONE, CZERO, CONE, SQRTm1,  &
                               TEN2m8, TEN2m6
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initSpinRotation,  &
          endSpinRotation,   &
          printSpinRotation, &
          resetSpinRotation, &
          rotateLtoG,        &
          rotateGtoL,        &
          transformDensityMatrix
!
   interface transformDensityMatrix
      module procedure transformDM1, transformDM2, transformDM3
   end interface
!
private
   integer (kind=IntKind) :: NumLocalSpinDirs
!
   logical :: Initialized = .false.
!
   real (kind=RealKind), parameter :: tol = TEN2m6
   real (kind=RealKind), allocatable :: evec(:,:)
!
   complex (kind=CmplxKind), allocatable :: u(:,:)
   complex (kind=CmplxKind), allocatable :: ud(:,:)
   complex (kind=CmplxKind), allocatable :: wx(:,:)
   complex (kind=CmplxKind), allocatable :: wy(:,:)
   complex (kind=CmplxKind), allocatable :: wz(:,:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSpinRotation(na,e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind) :: i
!
   real (kind=RealKind) :: e(3,na)
!
   if (na < 1) then
      call ErrorHandler('initSpinRotation',                            & 
                        'No. Local Spin Directions < 1',na)
   endif
!
   NumLocalSpinDirs=na
!
   allocate( u(4,na), ud(4,na), evec(3,na))
   allocate( wx(4,na), wy(4,na), wz(4,na))
!
   do i=1,na
      evec(1:3,i) = e(1:3,i)
      u(1,i) = sqrt(HALF*(ONE+evec(3,i)))
      if( abs(u(1,i)) < tol ) then
         u(1,i) = CZERO
         u(2,i) = CONE
         u(3,i) = CONE
      else
         u(2,i) =-HALF*(evec(1,i)+SQRTm1*evec(2,i))/u(1,i)
         u(3,i) = HALF*(evec(1,i)-SQRTm1*evec(2,i))/u(1,i)
      endif
      u(4,i)=u(1,i)
      ud(1,i)=conjg(u(1,i))
      ud(2,i)=conjg(u(3,i))
      ud(3,i)=conjg(u(2,i))
      ud(4,i)=conjg(u(4,i))
!     ----------------------------------------------------------------
      call u_sigma_u(u(1:4,i),ud(1:4,i),wx(1:4,i),wy(1:4,i),wz(1:4,i))
!     ----------------------------------------------------------------
   enddo
!
   Initialized = .true.
!
   end subroutine initSpinRotation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSpinRotation()
!  ===================================================================
   implicit none
!
   deallocate( u, ud, evec, wx, wy, wz)
   Initialized = .false.
   NumLocalSpinDirs = 0
!
   end subroutine endSpinRotation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine resetSpinRotation(id,e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: e(3)
!
   if (id < 1 .or. id > NumLocalSpinDirs ) then
      call ErrorHandler( 'resetSpinRotation',                          & 
                       'Wrong Local Spin Directions',NumLocalSpinDirs)
   endif
!
   if ( .not.Initialized ) then
      call ErrorHandler("resetSpinRotation",                           &
                   'initSpinRotation have to be initialized first')
   endif
!
   evec(1:3,id) = e(1:3)
   u(1,id) = sqrt(HALF*(ONE+evec(3,id)))
   if ( abs(u(1,id)) < tol ) then
      u(1,id) = CZERO
      u(2,id) = CONE
      u(3,id) = CONE
   else
      u(2,id) =-HALF*(evec(1,id)+SQRTm1*evec(2,id))/u(1,id)
      u(3,id) = HALF*(evec(1,id)-SQRTm1*evec(2,id))/u(1,id)
   endif
   u(4,id)  = u(1,id)
   ud(1,id) = conjg(u(1,id))
   ud(2,id) = conjg(u(3,id))
   ud(3,id) = conjg(u(2,id))
   ud(4,id) = conjg(u(4,id))
!  -------------------------------------------------------------------
   call u_sigma_u(u(1:4,id),ud(1:4,id),wx(1:4,id),wy(1:4,id),wz(1:4,id))
!  -------------------------------------------------------------------
!
   end subroutine resetSpinRotation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine rotateLtoG(id,nr,nc,ml1,ml2,mg)
!  ===================================================================
!
!  *******************************************************************
!  *  input::                                                        *
!  *       id:  the local spin index, same as the local atom index   *
!  *       nr:  number of rows of matrix in L-space                  *
!  *       nc:  number of columns of matrix in L-space               *
!  *       ml1: (1,1) element of matrix in the local spin space      *
!  *            dimension(nr,nc)                                     *
!  *       ml2: (2,2) element of matrix in the local spin space      *
!  *            dimension(nr,nc)                                     *
!  *                                                                 *
!  *  output::                                                       *
!  *       mg:  matrix in the global spin space                      *
!  *            dimension(2*nr,2*nc)                                 *
!  *                                                                 *
!  *  mg = ud * ml * u, where                                        *
!  *       u:   transformation matrix for spin space                 *
!  *       ud:  complex conjugate and transpose of u matrix          *
!  *******************************************************************
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: nr, nc
   integer (kind=IntKind) :: i, j
!
   complex (kind=CmplxKind), intent(in) :: ml1(:,:), ml2(:,:)
   complex (kind=CmplxKind), intent(out) :: mg(:,:)
!
   complex (kind=CmplxKind) :: a
   complex (kind=CmplxKind) :: b
!
   if (.not.Initialized) then
      call ErrorHandler('rotateLtoG','Need to initialize the module first')
   else if (id < 1 .or. id > NumLocalSpinDirs) then
      call ErrorHandler('rotateLtoG','invalid id',id)
   endif
!
!  ===================================================================
!  zero out mg
!  ===================================================================
   do j=1,2*nc
      do i=1,2*nr
         mg(i,j) = CZERO
      enddo
   enddo
!
   a=ud(1,id)*u(1,id)
   b=ud(3,id)*u(2,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,ml1(1,i),1,mg(1,i),1)
      call zaxpy(nr,b,ml2(1,i),1,mg(1,i),1)
!     ----------------------------------------------------------------
   enddo
   a=ud(2,id)*u(1,id)
   b=ud(4,id)*u(2,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,ml1(1,i),1,mg(nr+1,i),1)
      call zaxpy(nr,b,ml2(1,i),1,mg(nr+1,i),1)
!     ----------------------------------------------------------------
   enddo
   a=ud(1,id)*u(3,id)
   b=ud(3,id)*u(4,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,ml1(1,i),1,mg(1,nc+i),1)
      call zaxpy(nr,b,ml2(1,i),1,mg(1,nc+i),1)
!     ----------------------------------------------------------------
   enddo
   a=ud(2,id)*u(3,id)
   b=ud(4,id)*u(4,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,ml1(1,i),1,mg(nr+1,nc+i),1)
      call zaxpy(nr,b,ml2(1,i),1,mg(nr+1,nc+i),1)
!     ----------------------------------------------------------------
   enddo
!
   end subroutine rotateLtoG
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine rotateGtoL(id,nr,nc,mg,ml)
!  ===================================================================
!
!  *******************************************************************
!  *  input::                                                        *
!  *       id:  the local spin index, same as the local atom index   *
!  *       nr:  number of rows of matrix in G-space                  *
!  *       nc:  number of columns of matrix in G-space               *
!  *       mg:  matrix in the global spin space                      *
!  *            dimension(2*nr,2*nc)                                 *
!  *                                                                 *
!  *  output::                                                       *
!  *       ml:  matrix in the local spin space                       *
!  *            dimension(nr,nc,4)                                   *
!  *                                                                 *
!  *  ml = u * mg * ud, where                                        *
!  *       u:   transformation matrix for spin space                 *
!  *       ud:  complex conjugate and transpose of u matrix          *
!  *******************************************************************
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: nr, nc
   integer (kind=IntKind) :: i, j, k
!
   complex (kind=CmplxKind), intent(in) :: mg(:,:)
   complex (kind=CmplxKind), target, intent(out) :: ml(nr,nc,4)
!
   complex (kind=CmplxKind) :: a
   complex (kind=CmplxKind) :: b
!
   if (.not.Initialized) then
      call ErrorHandler('rotateGtoL','Need to initialize the module first')
   else if (id < 1 .or. id > NumLocalSpinDirs) then
      call ErrorHandler('rotateGtoL','invalid id',id)
   endif
!
!  ===================================================================
!  zero out ml.
!  ===================================================================
   do k=1,4
      do j=1,nc
         do i=1,nr
            ml(i,j,k) = CZERO
         enddo
      enddo
   enddo
!
   a=u(1,id)*ud(1,id)
   b=u(3,id)*ud(1,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,mg(1,i),   1,ml(1,i,1),1)
      call zaxpy(nr,b,mg(nr+1,i),1,ml(1,i,1),1)
!     ----------------------------------------------------------------
   enddo
   a=u(1,id)*ud(2,id)
   b=u(3,id)*ud(2,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,mg(1,nc+i),   1,ml(1,i,1),1)
      call zaxpy(nr,b,mg(nr+1,nc+i),1,ml(1,i,1),1)
!     ----------------------------------------------------------------
   enddo
!
   a=u(2,id)*ud(1,id)
   b=u(4,id)*ud(1,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,mg(1,i),   1,ml(1,i,2),1)
      call zaxpy(nr,b,mg(nr+1,i),1,ml(1,i,2),1)
!     ----------------------------------------------------------------
   enddo
   a=u(2,id)*ud(2,id)
   b=u(4,id)*ud(2,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,mg(1,nc+i),   1,ml(1,i,2),1)
      call zaxpy(nr,b,mg(nr+1,nc+i),1,ml(1,i,2),1)
!     ----------------------------------------------------------------
   enddo
!
   a=u(1,id)*ud(3,id)
   b=u(3,id)*ud(3,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,mg(1,i),   1,ml(1,i,3),1)
      call zaxpy(nr,b,mg(nr+1,i),1,ml(1,i,3),1)
!     ----------------------------------------------------------------
   enddo
   a=u(1,id)*ud(4,id)
   b=u(3,id)*ud(4,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,mg(1,nc+i),   1,ml(1,i,3),1)
      call zaxpy(nr,b,mg(nr+1,nc+i),1,ml(1,i,3),1)
!     ----------------------------------------------------------------
   enddo
!
   a=u(2,id)*ud(3,id)
   b=u(4,id)*ud(3,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,mg(1,i),   1,ml(1,i,4),1)
      call zaxpy(nr,b,mg(nr+1,i),1,ml(1,i,4),1)
!     ----------------------------------------------------------------
   enddo
   a=u(2,id)*ud(4,id)
   b=u(4,id)*ud(4,id)
   do i=1,nc
!     ----------------------------------------------------------------
      call zaxpy(nr,a,mg(1,nc+i),   1,ml(1,i,4),1)
      call zaxpy(nr,b,mg(nr+1,nc+i),1,ml(1,i,4),1)
!     ----------------------------------------------------------------
   enddo
!
   end subroutine rotateGtoL
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine u_sigma_u(cu,cud,usu_x,usu_y,usu_z)
!  ===================================================================
   implicit   none
!
   complex (kind=CmplxKind), intent(in) :: cu(4)
   complex (kind=CmplxKind), intent(in) :: cud(4)
   complex (kind=CmplxKind), intent(out) :: usu_x(4)
   complex (kind=CmplxKind), intent(out) :: usu_y(4)
   complex (kind=CmplxKind), intent(out) :: usu_z(4)
!
!  *******************************************************************
!  calculate:
!                  +                   +                   +
!     u * sigma * u ,     u * sigma * u ,     u * sigma * u
!              x                   y                   z
!  *******************************************************************
!
   usu_x(1)=cu(3)*cud(1)+cu(1)*cud(2)
   usu_x(2)=cu(4)*cud(1)+cu(2)*cud(2)
   usu_x(3)=cu(3)*cud(3)+cu(1)*cud(4)
   usu_x(4)=cu(4)*cud(3)+cu(2)*cud(4)
!
   usu_y(1)=(cu(3)*cud(1)-cu(1)*cud(2))*SQRTm1
   usu_y(2)=(cu(4)*cud(1)-cu(2)*cud(2))*SQRTm1
   usu_y(3)=(cu(3)*cud(3)-cu(1)*cud(4))*SQRTm1
   usu_y(4)=(cu(4)*cud(3)-cu(2)*cud(4))*SQRTm1
!
   usu_z(1)=cu(1)*cud(1)-cu(3)*cud(2)
   usu_z(2)=cu(2)*cud(1)-cu(4)*cud(2)
   usu_z(3)=cu(1)*cud(3)-cu(3)*cud(4)
   usu_z(4)=cu(2)*cud(3)-cu(4)*cud(4)
!
   end subroutine u_sigma_u
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine transformDM1(id,gm)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: id
!
   complex (kind=CmplxKind), intent(inout):: gm(1:4)
   complex (kind=CmplxKind) :: t1, t2, t3, t4
!
   if (.not.Initialized) then
      call ErrorHandler('transformDensityMatrix',                       &
                        'Need to initialize the module first')
   else if (id < 1 .or. id > NumLocalSpinDirs) then
      call ErrorHandler('transformDensityMatrix','invalid id',id)
   endif
!
   t1=gm(1)+gm(4)
   t2=(gm(1)*wx(1,id)+gm(4)*wx(4,id))+(gm(2)*wx(3,id)+gm(3)*wx(2,id))
   t3=(gm(1)*wy(1,id)+gm(4)*wy(4,id))+(gm(2)*wy(3,id)+gm(3)*wy(2,id))
   t4=(gm(1)*wz(1,id)+gm(4)*wz(4,id))+(gm(2)*wz(3,id)+gm(3)*wz(2,id))
   gm(1)=t1
   gm(2)=t2
   gm(3)=t3
   gm(4)=t4
!
   end subroutine transformDM1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine transformDM2(id,n1,gm)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: id, n1
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), intent(inout) :: gm(1:n1,1:4)
   complex (kind=CmplxKind) :: t1, t2, t3, t4
!
   if (.not.Initialized) then
      call ErrorHandler('transformDensityMatrix',                       &
                        'Need to initialize the module first')
   else if (id < 1 .or. id > NumLocalSpinDirs) then
      call ErrorHandler('transformDensityMatrix','invalid id',id)
   else if (n1 < 1) then
      call ErrorHandler('transformDensityMatrix','invalid first dim',n1)
   endif
!
   do i=1,n1
      t1=gm(i,1)+gm(i,4)
      t2=(gm(i,1)*wx(1,id)+gm(i,4)*wx(4,id))+(gm(i,2)*wx(3,id)+gm(i,3)*wx(2,id))
      t3=(gm(i,1)*wy(1,id)+gm(i,4)*wy(4,id))+(gm(i,2)*wy(3,id)+gm(i,3)*wy(2,id))
      t4=(gm(i,1)*wz(1,id)+gm(i,4)*wz(4,id))+(gm(i,2)*wz(3,id)+gm(i,3)*wz(2,id))
      gm(i,1)=t1
      gm(i,2)=t2
      gm(i,3)=t3
      gm(i,4)=t4
   enddo
!
   end subroutine transformDM2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine transformDM3(id,n1,n2,gm)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: id, n1, n2
   integer (kind=IntKind) :: i, j
!
   complex (kind=CmplxKind), target, intent(inout) :: gm(1:n1,1:n2,1:4)
   complex (kind=CmplxKind) :: t1, t2, t3, t4
!
   if (.not.Initialized) then
      call ErrorHandler('transformDensityMatrix',                     &
                        'Need to initialize the module first')
   else if (id < 1 .or. id > NumLocalSpinDirs) then
      call ErrorHandler('transformDensityMatrix','invalid id',id)
   else if (n1 < 1) then
      call ErrorHandler('transformDensityMatrix','invalid first dim',n1)
   else if (n2 < 1) then
      call ErrorHandler('transformDensityMatrix','invalid second dim',n2)
   endif
!
   do j=1,n2
      do i=1,n1
         t1=gm(i,j,1)+gm(i,j,4)
         t2=(gm(i,j,1)*wx(1,id)+gm(i,j,4)*wx(4,id))+                  &
            (gm(i,j,2)*wx(3,id)+gm(i,j,3)*wx(2,id))
         t3=(gm(i,j,1)*wy(1,id)+gm(i,j,4)*wy(4,id))+                  &
            (gm(i,j,2)*wy(3,id)+gm(i,j,3)*wy(2,id))
         t4=(gm(i,j,1)*wz(1,id)+gm(i,j,4)*wz(4,id))+                  &
            (gm(i,j,2)*wz(3,id)+gm(i,j,3)*wz(2,id))
         gm(i,j,1)=t1
         gm(i,j,2)=t2
         gm(i,j,3)=t3
         gm(i,j,4)=t4
      enddo
   enddo
!
   end subroutine transformDM3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printSpinRotation(id)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in), optional :: id
   integer (kind=IntKind) :: id_fst, id_lst, i
!
   write(6,'(/,80(''-''))')
   write(6,'(/,23x,a)')'**********************************'
   write(6,'( 23x,a )')'* Output from printSpinRotation *'
   write(6,'(23x,a,/)')'**********************************'
!
   if (present(id)) then
      id_fst=id
      id_lst=id
   else
      id_fst=1
      id_lst=NumLocalSpinDirs
   endif
!
   do i = id_fst,id_lst
      write(6,'(1x,a,3(1x,f12.8))')  "Local Evec:: ", evec(1:3,i)
      write(6,'(1x,a,8(1x,f12.8),/)') "U - elem  :: ", u(1:4,i)
      write(6,'(2x,2(3x,a,1x,2f12.8))')   "wx(1) =", wx(1,i),"wx(3) =",wx(3,i)
      write(6,'(2x,2(3x,a,1x,2f12.8),/)') "wx(2) =", wx(2,i),"wx(4) =",wx(4,i)
      write(6,'(2x,2(3x,a,1x,2f12.8))')   "wy(1) =", wy(1,i),"wy(3) =",wy(3,i)
      write(6,'(2x,2(3x,a,1x,2f12.8),/)') "wy(2) =", wy(2,i),"wy(4) =",wy(4,i)
      write(6,'(2x,2(3x,a,1x,2f12.8))')  "wz(1) =", wz(1,i),"wz(3) =",wz(3,i)
      write(6,'(2x,2(3x,a,1x,2f12.8),/)') "wz(2) =", wz(2,i),"wz(4) =",wz(4,i)
   enddo
!
   write(6,'(/,80(''-''))')
   end subroutine printSpinRotation
!  ===================================================================
end module SpinRotationModule
