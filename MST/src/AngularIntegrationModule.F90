!*********************************************************************
!* The purpose of this module is to establish a angular grid on      *
!* spherical surface and allow for an angular integration of a real  *
!* function multiplied by spherical harmonics over the surface.      *
!*                                                                   *
!* The public functions include                                      *
!*     getSphericalGridData:  returns data space for the underlying  *
!*                            spherical grid. The data is real type. *
!*     getTheta:              given an index of SphericalGridData,   *
!*                            returns the corresponding theta angle  *
!*     getPhi:                given an index of SphericalGridData,   *
!*                            returns the corresponding phi angle    *
!*     calAngularIntegration: with the data on the underlying sphere *
!*                            surface for a radius being established,*
!*                            it calculates the integration of the   *
!*                            data, multiplied by spherical harmonics*
!*                            over the spherical surface, i.e., the  *
!*                            solid angle.                           *
!*********************************************************************
module AngularIntegrationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
   use MathParamModule, only : CZERO, ZERO, TEN2m6, TEN2m8, Y0, HALF, ONE, TWO, PI, PI2
   use IntegerFactorsModule, only : kofj
!
   implicit none
!
public :: initAngularIntegration,    &
          endAngularIntegration,     &
          getNumSphericalGridPoints, &
          getSphericalGridData,      &
          getTheta,                  &
          getPhi,                    &
          getUnitVec,                &
          calAngularIntegration,     &
          retrieveSphHarmExpanData
!
   interface getSphericalGridData
      module procedure getSphericalGridData0, getSphericalGridData1
   end interface getSphericalGridData
!
   interface calAngularIntegration
      module procedure calAngularIntegration_s, calAngularIntegration_p
   end interface calAngularIntegration
!
   interface retrieveSphHarmExpanData
      module procedure retrieveSphHarmExpanData_s, retrieveSphHarmExpanData_p
   end interface retrieveSphHarmExpanData
!
private
   integer (kind=IntKind) :: lmax_max = -1
   integer (kind=IntKind) :: jmax_max = 0
   integer (kind=IntKind) :: kmax_max = 0
   integer (kind=IntKind) :: n_theta = 0
   integer (kind=IntKind) :: n_phi = 0
   integer (kind=IntKind) :: nr_max = 0
   integer (kind=IntKind) :: n_func = 0
   integer (kind=IntKind) :: num_angles = 0
   integer (kind=IntKind) :: d1_size = 0
   integer (kind=IntKind), allocatable :: ing2ip(:), ing2it(:)
!
   real (kind=RealKind), allocatable :: theta(:), phi(:), wght_theta(:), wght_phi(:)
   real (kind=RealKind), allocatable :: upos(:,:), fact(:)
   real (kind=RealKind), allocatable, target :: angular_data(:,:)
!
   complex (kind=CmplxKind), allocatable :: ylm(:,:)
   complex (kind=CmplxKind), allocatable :: expansion_data(:,:)
!
   logical :: Initialized = .false.
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initAngularIntegration(num_theta,num_phi,lmax,num_r,num_func)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: num_theta,num_phi,lmax
   integer (kind=IntKind), intent(in) :: num_r,num_func
!
   n_theta = num_theta
   n_phi = num_phi
   num_angles = num_theta*num_phi
   nr_max = num_r
   n_func = num_func
   lmax_max = lmax
   jmax_max = (lmax+1)*(lmax+2); jmax_max = jmax_max/2
   kmax_max = (lmax+1)**2
   d1_size = nr_max*n_func
!
   allocate(theta(n_theta), phi(n_phi), wght_theta(n_theta), wght_phi(n_phi))
   allocate(ylm(num_angles,kmax_max),upos(3,num_angles),fact(n_theta))
   allocate(angular_data(d1_size,num_angles))
   allocate(expansion_data(d1_size,jmax_max))
   allocate(ing2ip(num_angles), ing2it(num_angles))
!
!  -------------------------------------------------------------------
   call setAngularData()
!  -------------------------------------------------------------------
   angular_data = ZERO
   expansion_data = CZERO
!
   Initialized = .true.
!
   end subroutine initAngularIntegration
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endAngularIntegration()
!  ===================================================================
   implicit none
!
   deallocate(theta, phi, wght_theta, wght_phi, ylm, angular_data)
   deallocate(upos, fact, ing2ip, ing2it, expansion_data)
!
   lmax_max = -1
   jmax_max = 0
   kmax_max = 0
   n_theta = 0
   n_phi = 0
   nr_max = 0
   n_func = 0
   num_angles = 0
   d1_size = 0
!
   Initialized = .false.
!
   end subroutine endAngularIntegration
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setAngularData()
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer (kind=IntKind) :: it, ip, ing, kl
!
   real (kind=RealKind) :: cost, sint, cosp, sinp
!
   complex (kind=CmplxKind) :: ylm_ngl(kmax_max)
!
!  ===================================================================
!  Setting up a uniform grid with weight along phi
!  ===================================================================
   do ip = 1, n_phi
      wght_phi(ip) = PI2/n_phi
      phi(ip) = (ip-1)*wght_phi(ip)+HALF*wght_phi(ip)
   enddo
!
!  ===================================================================
!  Setting up a Gaussian quadrature grid with weight along theta
!  ===================================================================
   call gauleg(ZERO, PI, theta, wght_theta, n_theta)
!  -------------------------------------------------------------------
!
   ing = 0
   do it = 1, n_theta
      cost  = cos(theta(it))
      sint  = sin(theta(it))
      do ip = 1, n_phi
         ing = ing + 1
         ing2ip(ing) = ip
         ing2it(ing) = it
         sinp = sin(phi(ip))
         cosp = cos(phi(ip))
         upos(1,ing) = cosp*sint
         upos(2,ing) = sinp*sint
         upos(3,ing) = cost
!        -------------------------------------------------------
         call calYlm(upos(:,ing),lmax_max,ylm_ngl)
!        -------------------------------------------------------
         do kl = 1, kmax_max
            ylm(ing,kl) = ylm_ngl(kl)
         enddo
      enddo
      fact(it)  = sint*wght_theta(it)
   enddo
!
   end subroutine setAngularData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumSphericalGridPoints() result(n)
!  ===================================================================
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumSphericalGridPoints',                  &
                        'AngularIntegrationModule is not initialized')
   endif
!
   n = num_angles
!
   end function  getNumSphericalGridPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphericalGridData0() result(f)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: f(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getSphericalGridData',                       &
                        'AngularIntegrationModule is not initialized')
   endif
!
   f => angular_data(:,:)
!
   end function getSphericalGridData0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphericalGridData1(ing) result(f)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ing
!
   real (kind=RealKind), pointer :: f(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getSphericalGridData',                       &
                        'AngularIntegrationModule is not initialized')
   else if (ing < 1 .or. ing > num_angles) then
      call ErrorHandler('getSphericalGridData',                       &
                        'Angular data index is out of range',ing)
   endif
!
   f => angular_data(:,ing)
!
   end function getSphericalGridData1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getUnitVec(ing,theta_w,phi_w) result(pos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ing
   integer (kind=IntKind) :: ip, it
!   
   real (kind=RealKind), intent(out), optional :: theta_w, phi_w
   real (kind=RealKind) :: pos(3)
!
   if (.not.Initialized) then
      call ErrorHandler('getUnitVec',                                 &
                        'AngularIntegrationModule is not initialized')
   else if (ing < 1 .or. ing > num_angles) then
      call ErrorHandler('getUnitVec',                                 &
                        'Angular data index is out of range',ing)
   endif
!
   pos(1:3) = upos(1:3,ing)
!
   if (present(theta_w)) then
      it = ing2it(ing)
      theta_w = wght_theta(it)
   endif
!
   if (present(phi_w)) then
      ip = ing2ip(ing)
      phi_w = wght_phi(ip)
   endif
!
   end function getUnitVec
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTheta(ing,w) result(t)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ing
   integer (kind=IntKind) :: it
!   
   real (kind=RealKind), intent(out), optional :: w
   real (kind=RealKind) :: t
!
   if (.not.Initialized) then
      call ErrorHandler('getTheta',                                   &
                        'AngularIntegrationModule is not initialized')
   else if (ing < 1 .or. ing > num_angles) then
      call ErrorHandler('getTheta',                                   &
                        'Angular data index is out of range',ing)
   endif
!
   it = ing2it(ing)
   t = theta(it)
!
   if (present(w)) then
      w = wght_theta(it)
   endif
!
   end function getTheta
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPhi(ing,w) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ing
   integer (kind=IntKind) :: ip
!   
   real (kind=RealKind), intent(out), optional :: w
   real (kind=RealKind) :: p
!
   if (.not.Initialized) then
      call ErrorHandler('getPhi',                                     &
                        'AngularIntegrationModule is not initialized')
   else if (ing < 1 .or. ing > num_angles) then
      call ErrorHandler('getPhi',                                     &
                        'Angular data index is out of range',ing)
   endif
!
   ip = ing2ip(ing)
   p = phi(ip)
!
   if (present(w)) then
      w = wght_phi(ip)
   endif
!
   end function getPhi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calAngularIntegration_s()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: jl, kl, nf, it, ip, n, ing
!
   complex (kind=CmplxKind) :: pXC_rphi(d1_size)
   complex (kind=CmplxKind) :: cfac
!
   if (.not.Initialized) then
      call ErrorHandler('calAngularIntegration',                      &
                        'AngularIntegrationModule is not initialized')
   endif
!
   expansion_data = CZERO
   do jl = 1, jmax_max
      ing = 0
      kl = kofj(jl)
      do it = 1, n_theta
         pXC_rphi = CZERO
         do ip = 1, n_phi
            ing = ing + 1
            cfac = conjg(ylm(ing,kl))*wght_phi(ip)
            do n = 1, d1_size
               pXC_rphi(n) = pXC_rphi(n) + cfac*angular_data(n,ing)
            enddo
         enddo
         do n = 1, d1_size
            expansion_data(n,jl) = expansion_data(n,jl) + fact(it)*pXC_rphi(n)
         enddo
!        cfac = fact(it)
!        call zaxpy(d1_size,cfac,pXC_rphi,1,expansion_data(1,jl),1)
      enddo
   enddo
!
   end subroutine calAngularIntegration_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calAngularIntegration_p(MyPEinThisGroup,NumPEsInThisGroup,ThisGID)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: MyPEinThisGroup, NumPEsInThisGroup, ThisGID
!
   integer (kind=IntKind) :: jl, kl, nf, it, ip, n, ing
!
   complex (kind=CmplxKind) :: pXC_rphi(d1_size)
   complex (kind=CmplxKind) :: cfac
!
   if (.not.Initialized) then
      call ErrorHandler('calAngularIntegration',                      &
                        'AngularIntegrationModule is not initialized')
   endif
!
   expansion_data =CZERO
   do jl = MyPEinThisGroup+1, jmax_max, NumPEsInThisGroup
      ing = 0
      kl = kofj(jl)
      do it = 1, n_theta
         pXC_rphi = CZERO
         do ip = 1, n_phi
            ing = ing + 1
            cfac = conjg(ylm(ing,kl))*wght_phi(ip)
            do n = 1, d1_size
               pXC_rphi(n) = pXC_rphi(n) + cfac*angular_data(n,ing)
            enddo
         enddo
         do n = 1, d1_size
            expansion_data(n,jl) = expansion_data(n,jl) + fact(it)*pXC_rphi(n)
         enddo
!        cfac = fact(it)
!        call zaxpy(d1_size,cfac,pXC_rphi,1,expansion_data(1,jl),1)
      enddo
   enddo
   if (NumPEsInThisGroup > 1) then
!     ----------------------------------------------------------------
      call GlobalSumInGroup(ThisGID,expansion_data,d1_size,jmax_max)
!     ----------------------------------------------------------------
   endif
!
   end subroutine calAngularIntegration_p
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine retrieveSphHarmExpanData_s(numr,jmax,nf,func_L)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: numr, jmax, nf
   integer (kind=IntKind) :: jl, ir, n, size1, size2
!
   complex (kind=CmplxKind), intent(out) :: func_L(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('retrieveSphHarmExpanData',                   &
                        'AngularIntegrationModule is not initialized')
   endif
!
   size1 = size(func_L,1); size2 = size(func_L,2)
!
   if (size1 < numr) then
      call ErrorHandler('retrieveSphHarmExpanData','size1 < numr',size1,numr)
   else if (nr_max < numr) then
      call ErrorHandler('retrieveSphHarmExpanData','nr_max < numr',nr_max,numr)
   else if (size2 < jmax) then
      call ErrorHandler('retrieveSphHarmExpanData','size2 < jmax',size2,jmax)
   else if (kmax_max < kofj(jmax)) then
      call ErrorHandler('retrieveSphHarmExpanData','kmax_max < kofj(jmax)', &
                        kmax_max,kofj(jmax))
   else if (n_func < nf) then
      call ErrorHandler('retrieveSphHarmExpanData','n_func < nf',n_func,nf)
   else if (nf < 1) then
      call ErrorHandler('retrieveSphHarmExpanData','nf < 1',nf)
   endif
!
   func_L = CZERO
   do jl = 1, jmax
      do ir = 1, numr
         n = n_func*(ir-1) + nf
         func_L(ir,jl) = expansion_data(n,jl)
      enddo
   enddo
!  call zcopy(numr*jmax,expansion_data(nf,1),n_func,func_L,1)
!
   end subroutine retrieveSphHarmExpanData_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine retrieveSphHarmExpanData_p(numr,jmax,nf,func_L,         &
                                         MyPEinThisGroup,NumPEsInThisGroup,ThisGID)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: numr, jmax, nf
   integer (kind=IntKind), intent(in) :: MyPEinThisGroup, NumPEsInThisGroup, ThisGID
   integer (kind=IntKind) :: jl, ir, n, size1, size2, m
!
   complex (kind=CmplxKind), intent(out) :: func_L(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('retrieveSphHarmExpanData',                   &
                        'AngularIntegrationModule is not initialized')
   endif
!
   size1 = size(func_L,1); size2 = size(func_L,2)
!
   if (size1 < numr) then
      call ErrorHandler('retrieveSphHarmExpanData','size1 < numr',size1,numr)
   else if (nr_max < numr) then
      call ErrorHandler('retrieveSphHarmExpanData','nr_max < numr',nr_max,numr)
   else if (size2 < jmax) then
      call ErrorHandler('retrieveSphHarmExpanData','size2 < jmax',size2,jmax)
   else if (kmax_max < kofj(jmax)) then
      call ErrorHandler('retrieveSphHarmExpanData','kmax_max < kofj(jmax)', &
                        kmax_max,kofj(jmax))
   else if (n_func < nf) then
      call ErrorHandler('retrieveSphHarmExpanData','n_func < nf',n_func,nf)
   else if (nf < 1) then
      call ErrorHandler('retrieveSphHarmExpanData','nf < 1',nf)
   endif
!
   func_L = CZERO
!
   do jl = MyPEinThisGroup+1, jmax, NumPEsInThisGroup
      do ir = 1, numr
         n = n_func*(ir-1) + nf
         func_L(ir,jl) = expansion_data(n,jl)
      enddo
   enddo
!  ===================================================================
!  Using zcopy instead ...
!  ===================================================================
!! m = mod(jmax,NumPEsInThisGroup)
!! n = (jmax-m)/NumPEsInThisGroup
!! if (MyPEinThisGroup == NumPEsInThisGroup-1) then
!!    call zcopy(numr*(n+m),expansion_data(nf,MyPEinThisGroup*n+1),n_func, &
!!               func_L(1,MyPEinThisGroup*n+1),1)
!! else
!!    call zcopy(numr*n,expansion_data(nf,MyPEinThisGroup*n+1),n_func,&
!!               func_L(1,MyPEinThisGroup*n+1),1)
!! endif
!
   if (NumPEsInThisGroup > 1) then
!     ----------------------------------------------------------------
      call GlobalSumInGroup(ThisGID,func_L,size1,jmax)
!     ----------------------------------------------------------------
   endif
!
   end subroutine retrieveSphHarmExpanData_p
!  ===================================================================
end module AngularIntegrationModule
