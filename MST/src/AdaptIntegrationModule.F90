module AdaptIntegrationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO
!
public :: initAdaptIntegration,    &
          endAdaptIntegration,     &
          setupAdaptMesh,          &
          getAdaptMesh,            &
          getAdaptMeshWeight,      &
          getAdaptIntegration,     &
          getUniMeshIntegration,   &
          getUniformIntegration,   &
          getAuxDataAdaptIntegration, &
          getPeakPos,              &
          getWeightedIntegration !added by xianglin
!
private
   integer (kind=IntKind) :: AuxArraySize
   integer (kind=IntKind) :: AdaptMeshSize
   integer (kind=IntKind) :: GroupID, NumPEsInGroup, MyPEinGroup
!
   real (kind=RealKind), allocatable, target :: AdaptMesh(:)
   real (kind=RealKind), allocatable, target :: AdaptMeshWeight(:)
!
   real (kind=RealKind) :: LowEnd, HighEnd, peak_pos
!
   complex (kind=CmplxKind), allocatable :: AuxArray(:)
   complex (kind=CmplxKind), allocatable, target :: AuxArrayIntegral(:)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initAdaptIntegration(n,gid,rel_B)
!  ===================================================================
   use GroupCommModule, only : getNumPEsInGroup, getMyPEinGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in), optional :: gid
   logical, intent(in), optional :: rel_B
!
   AuxArraySize = n
   if (present(rel_B)) then
      if (rel_B) then
            AuxArraySize = 4*n
      endif
   endif
   allocate( AuxArray(AuxArraySize), AuxArrayIntegral(AuxArraySize) )
!
   if (present(gid)) then
      GroupID = gid
      NumPEsInGroup = getNumPEsInGroup(GroupID)
      MyPEinGroup = getMyPEinGroup(GroupID)
   else
      GroupID = -1
      NumPEsInGroup = 1
      MyPEinGroup = 0
   endif
!
   end subroutine initAdaptIntegration
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endAdaptIntegration()
!  ===================================================================
   implicit none
!
   deallocate( AuxArray, AuxArrayIntegral )
!
   if (allocated( AdaptMesh )) then
      deallocate( AdaptMesh )
   endif
   if (allocated( AdaptMeshWeight )) then
      deallocate( AdaptMeshWeight )
   endif
!
   end subroutine endAdaptIntegration
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAdaptMesh(n) result(xi)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out) :: n
!
   real (kind=RealKind), pointer :: xi(:)
!
   xi => AdaptMesh(:)
   n = AdaptMeshSize
!
   end function getAdaptMesh
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAdaptMeshWeight(n) result(wi)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out) :: n
!
   real (kind=RealKind), pointer :: wi(:)
!
   wi => AdaptMeshWeight(:)
   n = AdaptMeshSize
!
   end function getAdaptMeshWeight
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPeakPos() result(peakpos)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: peakpos
!
   peakpos = peak_pos
!
   end function getPeakPos
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function qFit(x1,x2,x3,f1,f2,f3)
!  ===================================================================
     implicit none
     real (kind=RealKind), dimension(3) :: qFit
     real (kind=RealKind), intent(in) :: x1,x2,x3
     real (kind=RealKind), intent(in) :: f1,f2,f3
!  ===================================================================
     real (kind=RealKind) :: dx21,dx31,dx32,df21,df31

     dx32 = x3 - x2
     dx31 = x3 - x1
     dx21 = x2 - x1
     df21 = f2 - f1
     df31 = f3 - f1
     qFit(3) = (df31/dx31 - df21/dx21)/dx32
     qFit(2) = df21/dx21 - qFit(3)*(x1 + x2)
     qFit(1) = f1 - x1*(qFit(2)+ qFit(3)*x1)

   end function qFit
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function qInt(iengy, irho)
     implicit none
     real (kind=RealKind) :: qInt
     real (kind=RealKind), intent(in), dimension(:),target :: iengy
     real (kind=RealKind), intent(in), dimension(:),target :: irho
!  ===================================================================
     integer (kind=IntKind) :: i
     real (kind=RealKind) :: l1,l2,dint,sum
     real (kind=RealKind), dimension(3) :: coef
     real (kind=RealKind), dimension(:),pointer :: x,f

     x => iengy
     f => irho
     l1 = x(1)
     sum = ZERO
     do i=2,size(x)-1
       l2 = (x(i) + x(i+1))/2
       coef = qFit(x(i-1),x(i),x(i+1),f(i-1),f(i),f(i+1))
       dint = coef(1)*(l2-l1)+coef(2)*(l2**2-l1**2)/2
       sum = sum + dint+coef(3)*(l2**3-l1**3)/3
       l1 = l2
     enddo
     l2 = x(size(x))
     dint = coef(1)*(l2-l1)+coef(2)*(l2**2-l1**2)/2
     qInt = sum + dint+coef(3)*(l2**3-l1**3)/3

   end function qInt
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupAdaptMesh(x0,x1,nmesh,info,func)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   real (kind=RealKind), intent(in) :: x0, x1
   integer (kind=IntKind), intent(in) :: nmesh
   integer (kind=IntKind), intent(in) :: info(*)

   integer (kind=IntKind) i, imax
   real (kind=RealKind) dx, phasemax
   real (kind=RealKind), dimension(nmesh) :: x,xnu,y,dy
   real (kind=RealKind), dimension(:), allocatable :: adx
   real (kind=RealKind), dimension(3) :: coef
!
   interface
      function func(info,x) result(y)
         use KindParamModule, only : IntKind, RealKind
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: x
         real (kind=RealKind) :: y
      end function func
   end interface
!
   if (allocated(AdaptMesh) ) then
      deallocate(AdaptMesh, AdaptMeshWeight)
   endif
   allocate( AdaptMesh(nmesh), AdaptMeshWeight(nmesh) )
   AdaptMeshSize = nmesh
!
   if (nmesh < 3) then
      call ErrorHandler('setupAdaptMesh','nmesh < 3',nmesh,.true.)
   endif
!
   if (x0 < x1) then
      LowEnd = x0
      HighEnd = x1
   else if (x0 > x1) then
      LowEnd = x1
      HighEnd = x0
   else
      call ErrorHandler('setupAdaptMesh','The variable range is ill defined',x0,x1,.true.)
   endif
!! 
!! The range is divided by 30, creating 31 uniformly distributed energy mesh points
!! and phase shifts are calculated on this uniform mesh for further refinement
   dx = (x1 - x0)/(nmesh-1)
   x(1) = x0
   y(1) = func(info,x(1))
   do i=2,nmesh
     x(i)=x(i-1) + dx
     y(i) = func(info,x(i))
   end do

!! First derivatives of the phase shift is calculated with help of cubic interpolation
   coef = qFit(x(1),x(2),x(3),y(1),y(2),y(3))
   dy(1) = abs(coef(2) + 2.0d0*coef(3)*x(1))
   do i=2,nmesh-1
     coef = qFit(x(i-1),x(i),x(i+1),y(i-1),y(i),y(i+1))
     dy(i) = abs(coef(2) + 2.0d0*coef(3)*x(i))
   enddo
   dy(nmesh) = abs(coef(2) + 2.0d0*coef(3)*x(nmesh))

!! Location of maximum of first derivative of the phase shift is found here
!! This is actually peak for the density of states
   imax = maxloc(dy,1)
!! Now non-uniform mesh is created using the steepnes or the first 
!! derivative of the phase shift
   if (imax == nmesh) then
      coef = qFit(x(imax-2),x(imax-1),x(imax),dy(imax-2),dy(imax-1),dy(imax))
   else if (imax == 1) then
      coef = qFit(x(imax),x(imax+1),x(imax+2),dy(imax),dy(imax+1),dy(imax+2))
   else
      coef = qFit(x(imax-1),x(imax),x(imax+1),dy(imax-1),dy(imax),dy(imax+1))
   endif
   x(imax) = -coef(2)/coef(3)/2.0d0
   phasemax = coef(1)+coef(2)*x(imax)+coef(3)*x(imax)**2
   phasemax = phasemax+1.d-08

   allocate(adx(imax))
   do i=1,imax
!     adx(i) = dlog(ONE + phasemax - dabs(dy(i)))
     adx(i) = dabs(phasemax - dabs(dy(i)))
   enddo
   adx = adx*(x(imax) - x(1))/sum(adx)
   xnu(1) = x(1)
   do i=2,imax-1
     xnu(i) = xnu(i-1) + adx(i-1)
   enddo
   deallocate(adx)

   xnu(imax) = x(imax)
   peak_pos = xnu(imax)
   if (nmesh-imax > 0) then   ! Added if-condition by Yang 6/8/2017
      allocate(adx(nmesh-imax))
      do i=1,nmesh-imax
!        adx(i) = dlog(ONE + phasemax - dabs(dy(i+imax)))
         adx(i) = dabs(phasemax - dabs(dy(i+imax)))
      enddo
      adx = adx*(x(nmesh) - x(imax))/sum(adx)
      do i=1,nmesh-imax
        xnu(i+imax) = xnu(i+imax-1) + adx(i)
      enddo
      deallocate(adx)
   endif
!
   AdaptMesh = xnu
!  
   end subroutine setupAdaptMesh
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAdaptIntegration(info,func,ngrid) result(fint)
!  ===================================================================
   use MathParamModule, only : ZERO, CZERO
   use ErrorHandlerModule, only : ErrorHandler
   use IntegrationModule, only : calIntegration
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind), intent(out) :: ngrid
   integer (kind=IntKind) :: i,j,imax, n_loc
   real (kind=RealKind) :: fint, rho0
   real (kind=RealKind), dimension(AdaptMeshSize) :: rho
   real (kind=RealKind), allocatable :: ehalf(:),rhohalf(:)
   complex (kind=CmplxKind), allocatable :: aa_mesh(:,:)
   complex (kind=CmplxKind) :: aa_int(AdaptMeshSize)
!
   interface
      function func(info,x,aux,xi,redundant) result(y)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: x
         real (kind=RealKind), intent(in), optional :: xi
         real (kind=RealKind) :: y
         complex (kind=CmplxKind), intent(out), target :: aux(:)
         logical, intent(in), optional :: redundant
      end function func
   end interface
!
   allocate( aa_mesh(AdaptMeshSize,AuxArraySize) )
   aa_mesh = CZERO; rho = ZERO
!
   ngrid = AdaptMeshSize
!
   if (NumPEsInGroup <= ngrid) then
      n_loc = ngrid/NumPEsInGroup
   else
      n_loc = 1
   endif
   do i = MyPEinGroup+1, n_loc*NumPEsInGroup, NumPEsInGroup
     if (i <= ngrid) then
        rho(i) = func(info,AdaptMesh(i),AuxArray)
        do j = 1, AuxArraySize
          aa_mesh(i,j) = AuxArray(j)
        enddo
     else ! run redundant calculations. The results are not stored
        rho0 = func(info,AdaptMesh(ngrid),AuxArray)
     endif
   enddo
   if (NumPEsInGroup > 1) then
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,rho,ngrid)
      call GlobalSumInGroup(GroupID,aa_mesh,ngrid,AuxArraySize)
!     ----------------------------------------------------------------
    endif
   do i = n_loc*NumPEsInGroup+1, ngrid
     rho(i) = func(info,AdaptMesh(i),AuxArray)
     do j = 1, AuxArraySize
       aa_mesh(i,j) = AuxArray(j)
     enddo
   enddo
!
   imax = maxloc(rho,1)
   allocate(ehalf(imax),rhohalf(imax))
   do i=1,imax
     ehalf(i) = AdaptMesh(i)
     rhohalf(i) = rho(i)
   enddo
   fint = qInt(ehalf,rhohalf)
   deallocate(ehalf, rhohalf)

   allocate(ehalf(ngrid-imax+1),rhohalf(ngrid-imax+1))
   do i=1,ngrid-imax+1
     ehalf(i) = AdaptMesh(i+imax-1)
     rhohalf(i) = rho(i+imax-1)
   enddo
   fint = fint+qInt(ehalf,rhohalf)
   deallocate(ehalf, rhohalf)

   do j = 1, AuxArraySize
!     ----------------------------------------------------------------
      call calIntegration(0,ngrid,AdaptMesh,aa_mesh(1:ngrid,j),aa_int)
!     ----------------------------------------------------------------
      AuxArrayIntegral(j) = aa_int(ngrid)
   enddo
!
   deallocate( aa_mesh )
!
   end function getAdaptIntegration
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getUniformIntegration(n,x0,x,info,func,nm) result(fint)
!  ===================================================================
!  Note: For using different number of processors, the integration
!        mesh can potentially be different, thus the result could be 
!        slightly different. Increasing the number of uniform energy
!        grid points will in general help reducing the difference.
!  ===================================================================
   use MathParamModule, only : CZERO, ZERO, TEN2m8
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use IntegrationModule, only : calIntegration
!
   use MPPModule, only : setCommunicator, resetCommunicator
   use MPPModule, only : nbrecvMessage, nbsendMessage, WaitMessage, bcastMessage
!
   use GroupCommModule, only : GlobalSumInGroup, GlobalMaxInGroup
   use GroupCommModule, only : getGroupCommunicator
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind), intent(out) :: nm
   integer (kind=IntKind) :: i, j, n_loc, m_loc, imax, ns, ms
   integer (kind=IntKind) :: msgid0, msgid1, msgid2, msgid3, comm
!
   real (kind=RealKind), intent(in) :: x0, x
   real (kind=RealKind) :: fint, del, f0, peak_val
   real (kind=RealKind) :: x_mesh(n)
   real (kind=RealKind), allocatable :: f_mesh(:), y_mesh(:)
!
   complex (kind=CmplxKind), allocatable :: aa_mesh(:,:), aa_int(:)
   complex (kind=CmplxKind), allocatable :: a_buf(:), b_buf(:)
!
   interface
      function func(info,x,aux,xi,redundant) result(y)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: x
         real (kind=RealKind), intent(in), optional :: xi
         real (kind=RealKind) :: y
         complex (kind=CmplxKind), intent(out), target :: aux(:)
         logical, intent(in), optional :: redundant
      end function func
   end interface
!
   if (n < 2) then
!     ----------------------------------------------------------------
      call ErrorHandler('getUniformIntegration','Number of mesh < 2',n)
!     ----------------------------------------------------------------
   else if (x0 > x) then
!     ----------------------------------------------------------------
      call ErrorHandler('getUniformIntegration','x0 > x1',x0,x)
!     ----------------------------------------------------------------
   endif
!
   nm = n
   del = (x-x0)/real(n-1,kind=RealKind)
   do i = 1, n
      x_mesh(i) = x0 + (i-1)*del
   enddo
!
   if (n >= NumPEsInGroup) then
      n_loc = n/NumPEsInGroup ! number of energy points mapped onto the local process
      ns = MyPEinGroup*n_loc
   else ! This may happen if number of CPU Cores > n. No parallelization is done
      n_loc = n
      ns = 1
   endif
!
   if (MyPEinGroup == 0) then
      ms = 0
   else
      ms = 1
   endif
!
   allocate( f_mesh(n_loc+ms), y_mesh(n_loc+ms) )
   allocate( aa_mesh(n_loc+ms,AuxArraySize), aa_int(n_loc+ms) )
   allocate(a_buf(AuxArraySize), b_buf(AuxArraySize))
   f_mesh = ZERO; y_mesh = ZERO
   aa_mesh = CZERO; aa_int = CZERO
!
   do i = 1, n_loc
      f_mesh(ms+i) = func(info,x_mesh(ns+i),AuxArray)
      do j = 1, AuxArraySize
         aa_mesh(ms+i,j) = AuxArray(j)
      enddo
   enddo
   imax = maxloc(f_mesh,1)
   peak_pos = x_mesh(ns+imax-ms)
   peak_val = f_mesh(imax)
!
   if (n_loc < n) then
      comm = getGroupCommunicator(GroupID)
!     ----------------------------------------------------------------
      call setCommunicator(comm,MyPEinGroup,NumPEsInGroup)
!     ----------------------------------------------------------------
      if (MyPEinGroup > 0) then
!        -------------------------------------------------------------
         msgid0 = nbrecvMessage(f_mesh(1),10001,MyPEinGroup-1)
         msgid1 = nbrecvMessage(a_buf,AuxArraySize,10002,MyPEinGroup-1)
!        -------------------------------------------------------------
      endif
      if (MyPEinGroup < NumPesInGroup-1) then
!        -------------------------------------------------------------
         msgid2 = nbsendMessage(f_mesh(n_loc+ms),10001,MyPEinGroup+1)
!        -------------------------------------------------------------
         do j = 1, AuxArraySize
            b_buf(j) = aa_mesh(n_loc+ms,j)
         enddo
!        -------------------------------------------------------------
         msgid3 = nbsendMessage(b_buf,AuxArraySize,10002,MyPEinGroup+1)
!        -------------------------------------------------------------
      endif
      if (MyPEinGroup > 0) then
!        -------------------------------------------------------------
         call waitMessage(msgid0)
         call waitMessage(msgid1)
!        -------------------------------------------------------------
         do j = 1, AuxArraySize
            aa_mesh(1,j) = a_buf(j)
         enddo
      endif
      if (MyPEinGroup < NumPEsInGroup-1) then
!        -------------------------------------------------------------
         call waitMessage(msgid2)
         call waitMessage(msgid3)
!        -------------------------------------------------------------
      endif
!     ----------------------------------------------------------------
      call resetCommunicator()
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   call calIntegration(0,n_loc+ms,x_mesh(ns+1-ms:),f_mesh,y_mesh)
!  -------------------------------------------------------------------
   fint = y_mesh(n_loc+ms)
!
   do j = 1, AuxArraySize
!     ----------------------------------------------------------------
      call calIntegration(0,n_loc+ms,x_mesh(ns+1-ms:),aa_mesh(:,j),aa_int)
!     ----------------------------------------------------------------
      AuxArrayIntegral(j) = aa_int(n_loc+ms)
   enddo
!
   if (n_loc < n) then ! In case parallelization is performed
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,fint)
!     ----------------------------------------------------------------
      f0 = peak_val
!     ----------------------------------------------------------------
      call GlobalMaxInGroup(GroupID,f0)
!     ----------------------------------------------------------------
      if (abs(f0-peak_val) > TEN2m8) then
         peak_pos = ZERO
      endif
      peak_val = f0
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,peak_pos)
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,AuxArrayIntegral,AuxArraySize)
!     ----------------------------------------------------------------
   endif
!
   ns = n_loc*NumPEsInGroup
   if (ns < n) then
!     ----------------------------------------------------------------
      call setCommunicator(comm,MyPEinGroup,NumPEsInGroup)
!     ----------------------------------------------------------------
      if (MyPEinGroup == NumPEsInGroup-1) then
         f_mesh(1) = f_mesh(n_loc+ms)
         do j = 1, AuxArraySize
            b_buf(j) = aa_mesh(n_loc+ms,j)
         enddo
      endif
      call bcastMessage(f_mesh(1),NumPEsInGroup-1)
      call bcastMessage(b_buf,AuxArraySize,NumPEsInGroup-1)
!     ----------------------------------------------------------------
      call resetCommunicator()
!     ----------------------------------------------------------------
      do j = 1, AuxArraySize
          aa_mesh(1,j) = b_buf(j)
      enddo
      do i = 1, n-ns
         f_mesh(i+1) = func(info,x_mesh(ns+i),AuxArray,redundant=.true.)
         do j = 1, AuxArraySize
            aa_mesh(i+1,j) = AuxArray(j)
         enddo
      enddo
      imax = maxloc(f_mesh(2:),1)
      if (peak_val < f_mesh(imax)) then
         peak_val = f_mesh(imax)
         peak_pos = x_mesh(ns+imax-1)
      endif
!     ----------------------------------------------------------------
      call calIntegration(0,n-ns+1,x_mesh(ns:),f_mesh,y_mesh)
!     ----------------------------------------------------------------
      fint = fint + y_mesh(n-ns+1)
!
      do j = 1, AuxArraySize
!        -------------------------------------------------------------
         call calIntegration(0,n-ns+1,x_mesh(ns:),aa_mesh(:,j),aa_int)
!        -------------------------------------------------------------
         AuxArrayIntegral(j) = AuxArrayIntegral(j) + aa_int(n-ns+1)
      enddo
   endif
!
   deallocate( aa_mesh, aa_int, f_mesh, y_mesh )
   deallocate(a_buf, b_buf)
!
   end function getUniformIntegration
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAuxDataAdaptIntegration() result(fint)
!  ===================================================================
   complex (kind=CmplxKind), pointer :: fint(:)
!
   fint => AuxArrayIntegral(1:AuxArraySize)
!
   end function getAuxDataAdaptIntegration
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getUniMeshIntegration(n,x0,x,info,func,nm) result(fint)
!  ===================================================================
!  Use Boole's rule, also known as Bode's rule, to perform the integration
!  Note: For using different number of processors, the integration
!        mesh could be different, thus the result could be different
!  ===================================================================
   use MathParamModule, only : CZERO, ZERO, TWO, SEVEN
   use ErrorHandlerModule, only : ErrorHandler
   use GroupCommModule, only : GlobalSumInGroup, GlobalCollectInGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind), intent(out) :: nm
   integer (kind=IntKind) :: i, j, imax, m, m4
!
   real (kind=RealKind), intent(in) :: x0, x
   real (kind=RealKind) :: fint, h, f_mesh, x1, x2, x3, x4, x5, p, fac
   real (kind=RealKind), allocatable :: peak_dos(:,:)
!
   interface
      function func(info,x,aux,xi,redundant) result(y)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: x
         real (kind=RealKind), intent(in), optional :: xi
         real (kind=RealKind) :: y
         complex (kind=CmplxKind), intent(out), target :: aux(:)
         logical, intent(in), optional :: redundant
      end function func
   end interface
!
   if (n < 4) then
!     ----------------------------------------------------------------
      call ErrorHandler('getUniMeshIntegration','Number of mesh < 4',n)
!     ----------------------------------------------------------------
   else if (x0 > x) then
!     ----------------------------------------------------------------
      call ErrorHandler('getUniMeshIntegration','x0 > x1',x0,x)
!     ----------------------------------------------------------------
   endif
!
   m = n/4
   j = NumPEsInGroup
   do while (m > j)
      j = j + NumPEsInGroup
   enddo
   m = j    ! m is a multiple of NumPEsInGroup, thus it depends on the 
            ! total number of processors
   m4 = 4*m ! Divide [x0,x] into 4*m intervals. Every 4 intervals or 5 points
            ! forms a segment. Thus, there are m joining segments in total.
            ! The integration over each segment is performed using Bode's rule.
   h = (x-x0)/real(m4,kind=RealKind)
   nm = m4+1  ! nm = the total number of mesh points
!
   allocate(peak_dos(2,NumPEsInGroup))
   peak_dos = ZERO
!
   f_mesh = func(info,x0,AuxArray,redundant=.true.)
   peak_dos(1,MyPEinGroup+1) = f_mesh
   peak_dos(2,MyPEinGroup+1) = x0
   fac = SEVEN/real(NumPEsInGroup,RealKind)
   fint = fac*f_mesh
   AuxArrayIntegral = fac*AuxArray
   do i = MyPEinGroup+1, m, NumPEsInGroup
      j = (i-1)*4; x1 = x0 + j*h
      x2 = x0 + (j+1)*h; x3 = x0 + (j+2)*h; x4 = x0 + (j+3)*h; x5 = x0 + (j+4)*h
!
      f_mesh = func(info,x2,AuxArray)
      if (peak_dos(1,MyPEinGroup+1) < f_mesh) then
         peak_dos(1,MyPEinGroup+1) = f_mesh
         peak_dos(2,MyPEinGroup+1) = x2
      endif
      fint = fint + 32.0d0*f_mesh
      AuxArrayIntegral = AuxArrayIntegral + 32.0d0*AuxArray
!
      f_mesh = func(info,x3,AuxArray)
      if (peak_dos(1,MyPEinGroup+1) < f_mesh) then
         peak_dos(1,MyPEinGroup+1) = f_mesh
         peak_dos(2,MyPEinGroup+1) = x3
      endif
      fint = fint + 12.0d0*f_mesh
      AuxArrayIntegral = AuxArrayIntegral + 12.0d0*AuxArray
!
      f_mesh = func(info,x4,AuxArray)
      if (peak_dos(1,MyPEinGroup+1) < f_mesh) then
         peak_dos(1,MyPEinGroup+1) = f_mesh
         peak_dos(2,MyPEinGroup+1) = x4
      endif
      fint = fint + 32.0d0*f_mesh
      AuxArrayIntegral = AuxArrayIntegral + 32.0d0*AuxArray
!
      f_mesh = func(info,x5,AuxArray)
      if (peak_dos(1,MyPEinGroup+1) < f_mesh) then
         peak_dos(1,MyPEinGroup+1) = f_mesh
         peak_dos(2,MyPEinGroup+1) = x5
      endif
      if (i == m) then
         fint = fint + SEVEN*f_mesh
         AuxArrayIntegral = AuxArrayIntegral + SEVEN*AuxArray
      else
         fint = fint + 14.0d0*f_mesh
         AuxArrayIntegral = AuxArrayIntegral + 14.0d0*AuxArray
      endif
   enddo

   if (NumPEsInGroup > 1) then
!     ----------------------------------------------------------------
      call GlobalCollectInGroup(GroupID,peak_dos,2)
      call GlobalSumInGroup(GroupID,fint)
      call GlobalSumInGroup(GroupID,AuxArrayIntegral,AuxArraySize)
!     ----------------------------------------------------------------
   endif
   p = ZERO
   peak_pos = ZERO
   do i = 1, NumPEsInGroup
      if (p < peak_dos(1,i)) then
         p = peak_dos(1,i)
         peak_pos = peak_dos(2,i)
      endif
   enddo
   fac = h/22.5d0
   fint = fac*fint
   AuxArrayIntegral = fac*AuxArrayIntegral
!
   deallocate( peak_dos )
!
   end function getUniMeshIntegration
!  ===================================================================

!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWeightedIntegration(n,x0,x,info,func,n_pole,jost_pole) result(fint)
!  ===================================================================
   use MathParamModule, only : CZERO, ZERO
   use ErrorHandlerModule, only : ErrorHandler
   use IntegrationModule, only : calIntegration
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind), intent(in) :: n_pole
   complex (kind=CmplxKind) :: jost_pole(*)
!   integer (kind=IntKind), intent(out) :: nm
   integer (kind=IntKind) :: i, j, n_loc, imax
!
   real (kind=RealKind), intent(in) :: x0, x
   real (kind=RealKind) :: fint, del, f0
   real (kind=RealKind) :: x_mesh(n), f_mesh(n), y_mesh(n)
!
!  complex (kind=CmplxKind) :: aa_mesh(n,AuxArraySize)
   complex (kind=CmplxKind), allocatable :: aa_mesh(:,:)
   complex (kind=CmplxKind) :: aa_int(n)
!
   interface
      function func(info,x,aux,xi,redundant) result(y)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: x
         real (kind=RealKind), intent(in), optional :: xi
         real (kind=RealKind) :: y
         complex (kind=CmplxKind), intent(out), target :: aux(:)
         logical, intent(in), optional :: redundant
      end function func
   end interface
!
   if (n < 2) then
!     ----------------------------------------------------------------
      call ErrorHandler('getWeightedIntegration','Number of mesh < 2',n)
!     ----------------------------------------------------------------
   else if (x0 > x) then
!     ----------------------------------------------------------------
      call ErrorHandler('getWeightedIntegration','x0 > x1',x0,x)
!     ----------------------------------------------------------------
   endif
!
   allocate( aa_mesh(n,AuxArraySize) )
   aa_mesh = CZERO; f_mesh = ZERO
!
   call root_find( n,x0,x,x_mesh,n_pole,jost_pole )
!   open(202,file='E_mesh',action='write')
!   write(202,*) x_mesh
!   write(202,*) 'AuxArraySize=', AuxArraySize
!   close(202)

!   open(201,file='Poles',action='write')
!   write(201,*) "n_pole=",n_pole
!   do i=1,n_pole
!      write(201,*) jost_pole(i)
!   enddo
!   close(201)

   if (n >= NumPEsInGroup) then
      n_loc = n/NumPEsInGroup
   else ! This may happen if number of CPU Cores > n
      n_loc = 1
   endif
   do i = MyPEinGroup+1, n_loc*NumPEsInGroup, NumPEsInGroup
      if (i <= n) then
         f_mesh(i) = func(info,x_mesh(i),AuxArray)
         do j = 1, AuxArraySize
            aa_mesh(i,j) = AuxArray(j)
         enddo
      else ! run redundant calculations. The results are not stored
         f0 = func(info,x_mesh(n),AuxArray)
      endif
   enddo

 
   imax = maxloc(f_mesh,1)
   peak_pos = x_mesh(imax)
   if (NumPEsInGroup > 1) then
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,f_mesh,n)
      call GlobalSumInGroup(GroupID,aa_mesh,n,AuxArraySize)
!     ----------------------------------------------------------------
   endif
   do i = n_loc*NumPEsInGroup+1, n
      f_mesh(i) = func(info,x_mesh(i),AuxArray,redundant=.true.)
      do j = 1, AuxArraySize
         aa_mesh(i,j) = AuxArray(j)
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call calIntegration(0,n,x_mesh,f_mesh,y_mesh)
!  -------------------------------------------------------------------
   fint = y_mesh(n)
!
   do j = 1, AuxArraySize
!     ----------------------------------------------------------------
      call calIntegration(0,n,x_mesh,aa_mesh(1:n,j),aa_int)
!     ----------------------------------------------------------------
      AuxArrayIntegral(j) = aa_int(n)
   enddo
!
!   open(203,file='AuxArrayIntegral',action='write')
!   write(203,*) AuxArrayIntegral(1:AuxArraySize/4)
!   close(203)
   open(204,file='fmesh',action='write')
   write(204,'(a20)') 'x_mesh (energy)'
   write(204,*) x_mesh
   write(204,'(a20)') 'f_mesh (DOS)'
   write(204,*) f_mesh
   write(204,'(a20)') 'y_mesh (IDOS)'
   write(204,*) y_mesh
   close(204)
   deallocate( aa_mesh )
!
   end function getWeightedIntegration

!  ===================================================================
   subroutine root_find( N_energy,E_min,E_max,E_mesh,n_pole,jost_pole )
!  ===================================================================
   use KindParamModule, only : RealKind, IntKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
   implicit none
!
   integer (kind=IntKind), intent(in) :: N_energy
   real (RealKind), intent (in) :: E_min,E_max
   integer (kind=IntKind), intent(in) :: n_pole
   complex (kind=CmplxKind) :: jost_pole(*)
   real (RealKind), intent (out) :: E_mesh(N_energy)
   real (kind=RealKind) :: x, x1, x_min, x_max
   real (kind=RealKind) :: y_min,y_max, y_mesh(N_energy),dy, fun_min,fun_max,fun_mid
   real (kind=RealKind), parameter :: epsilon=1.0d-8
   integer (kind=IntKind) :: i,j,m,n
   ! find y_mesh(N_energy)
   if (E_min > E_max) then
      call ErrorHandler('root_newton', 'E_min > E_max')
   endif
!
   y_min=weight(E_min,0.d0,n_pole,jost_pole)
   y_max=weight(E_max,0.d0,n_pole,jost_pole)
   if (N_energy>1) then
      dy=(y_max-y_min)/(N_energy-1)
      do j=1,N_energy
         y_mesh(j)=y_min+(j-1)*dy
      enddo
   else if (N_energy==1) then
      y_mesh(1)=y_min
   else
      call ErrorHandler('root_newton', 'N_energy<0')
   endif
!
   ! find E_mesh
   E_mesh(1)=E_min
   E_mesh(N_energy)=E_max
   !bisection method to find roots
   do j=2,N_energy-1
      x_min=E_mesh(j-1)
      x_max=E_max
      fun_min=y_mesh(j-1)-y_mesh(j)
      fun_max=y_mesh(N_energy)-y_mesh(j)
      do i=1,100
         if ( ABS(fun_max-fun_min)<epsilon ) then
            exit
         else
            x=(x_min+x_max)/2.d0
            fun_mid=weight(x,y_mesh(j),n_pole,jost_pole)
            if (fun_min*fun_mid<0.d0) then
               x_max=x
               fun_max=fun_mid
            else
               x_min=x
               fun_min=fun_mid
            endif
            !print*,'j=',j,' i=',i, ' x=',x
         endif
      enddo
      E_mesh(j)=(x_min+x_max)/2.d0
   enddo
!
   end subroutine root_find

!  ===================================================================
   function weight(x,y0,n_pole,jost_pole) result(y)
!  ===================================================================
   use MathParamModule, only : PI
   use KindParamModule, only : RealKind, IntKind, CmplxKind
   implicit none
   real (kind=RealKind), intent(in) :: x,y0
   integer (kind=IntKind), intent(in) :: n_pole
   complex (kind=CmplxKind), intent(in) :: jost_pole(*)
   real (kind=RealKind) :: y
   integer (kind=IntKind) :: i,j
   real (kind=RealKind), parameter :: back_ground=2.d0
   real (kind=RealKind) :: x0,lambda
   y=back_ground*x-y0
   do i=1,N_pole
      x0=real( jost_pole(i) )
      lambda=ABS(imag( jost_pole(i) ))
      y=y+1.d0/PI*ATAN( (x-x0)/lambda )+.5d0
   enddo
   end function weight
!  ===================================================================
end module AdaptIntegrationModule
