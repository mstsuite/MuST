program testParallelFFT
!  ********************************************************************
!  test MPI FFTW parallel fast Fourier transform routines using Fortran 2003
!
!  The function to be Fourier transformed is:
!
!  f(x,y,z) = exp(alpha*x+beta*y+gamma*z), 0<x<a, 0<y<b, 0<z<c
!
!  F(kx,ky,kz) = (exp(alpha*a)-1)/(alpha+i*kx)*
!                (exp(beta*b)-1)/(beta+i*ky)*(exp(gamma*c)-1)/(gamma+i*kz)
!
!  where kx = 2*pi*jx/a, ky = 2*pi*jy/b, kz = 2*pi*jz/c, with
!        jx, jy, jz = ..., -2, -1, 0, +1, +2, ...
!  ********************************************************************
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use TimerModule, only : initTimer, getTime
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, CONE, PI2, SQRTm1
!
   use PrimeFactorsModule, only : getSubFactors
!
   use MPPModule, only : initMPP, endMPP, getCommunicator, bcastMessage
   use MPPModule, only : MyPE, NumPEs, syncAllPEs, GlobalSum
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use Uniform3DGridModule, only : initUniform3DGrid, endUniform3DGrid
   use Uniform3DGridModule, only : createUniform3DGrid, printUniform3DGrid
   use Uniform3DGridModule, only : createProcessorMesh
   use Uniform3DGridModule, only : distributeUniformGrid
!
   use ParallelFFTModule, only : initParallelFFT, endParallelFFT,         &
                                 getLocalArraySize, getNumLocalGridPoints,&
                                 allocateFunctionSpace,                   &
                                 aliasFunctionSpace,                      &
                                 getGridPointCoord, performTransformR2C,  &
                                 getProcessorMesh
   use ParallelFFTModule, only : getGlobalGridIndexRange
!
   implicit   none
!
   logical :: found
!
   integer (kind=IntKind), parameter :: man_limit = 64*64*64
!
   integer (kind=IntKind) :: local_n, local_start, nump_local
   integer (kind=IntKind) :: nx, ny, nz, nt, nt_fftw
!
   integer (kind=IntKind) :: i, j, k, idk, idr, idyz, idm, ik
   integer (kind=IntKind) :: ifc, jfc, kfc
   integer (kind=IntKind) :: fs, fe, numkx, numky, numkz
   integer (kind=IntKind) :: jconj, kconj
   integer (kind=IntKind) :: comm, m
   integer (kind=IntKind) :: ibuf(3)
!
!  ===================================================================
!  p3dfft related variables
!  ===================================================================
   integer (kind=IntKind), pointer :: factors(:,:)
   integer (kind=IntKind) :: ig, jg, kg
   integer (kind=IntKind) :: memsize, m1, m2, m3, n1, n2, n3
   integer (kind=IntKind) :: grid_start(3), grid_end(3), gir(3,3)
   real (kind=RealKind), pointer :: rfunc_fftwr(:)
   real (kind=RealKind), pointer :: p_rfunc_fftwr(:,:,:)
   real (kind=RealKind) :: p(3), kvec_fftw(3)
!  ===================================================================
!
   real (kind=RealKind) :: t0, vfac, vol
   real (kind=RealKind) :: a, b, c
   real (kind=RealKind) :: x, y, z, dx, dy, dz
   real (kind=RealKind) :: kx, ky, kz
   real (kind=RealKind) :: ka, kb, kc
   real (kind=RealKind) :: alpha, beta, gamma, fbuf(6), cell(3,3)
   real (kind=RealKind), allocatable :: func(:), kvec(:,:), rvec(:,:)
   real (kind=RealKind), allocatable :: kvec_nr(:,:)
!
   complex (kind=CmplxKind) :: ftx, fty, ftz, expfac
   complex (kind=CmplxKind), allocatable :: ft(:)
   complex (kind=CmplxKind), allocatable :: ftm(:)
   complex (kind=CmplxKind), allocatable :: ft_a(:), ft_nr(:)
   complex (kind=CmplxKind), pointer :: ft_fftw_ip(:)
!
!  -------------------------------------------------------------------
   call initMPP()
!  -------------------------------------------------------------------
!
   if (MyPE == 0) then
      fbuf = ZERO; ibuf = 0
!
      do while (fbuf(1) <= TEN2m6 .or. fbuf(2) <= TEN2m6 .or. fbuf(3) <= TEN2m6)
         write(6,'(1x,a,$)')'Box size in x, y, z direction: '
         read(5,*)fbuf(1:3)
      enddo
      write(6,'(f10.5,'','',f10.5,'','',f10.5)')fbuf(1:3)
!
      do while (ibuf(1) < 1 .or. ibuf(2) < 1 .or. ibuf(3) < 1)
         write(6,'(1x,a,$)')'Number of mesh points in x, y, z direction: '
         read(5,*)ibuf(1:3)
      enddo
      write(6,'(i5,'','',i5,'','',i5)')ibuf(1:3)
!
      write(6,'(/,1x,a)')    &
       'This Code Will Perform Fast Fourier Transform of an Exponential Function'
      write(6,'(//,1x,a)')'                 exp(alpha*x + beta*y + gamma*z)'
!
      write(6,'(//,1x,a,$)')'Enter alpha, beta, gamma: '
      read(5,*)fbuf(4:6)
      write(6,'(f10.5,'','',f10.5,'','',f10.5)')fbuf(4:6)
   endif
!  -------------------------------------------------------------------
   call bcastMessage(fbuf,6,0)
   call bcastMessage(ibuf,3,0)
!  -------------------------------------------------------------------
   a = fbuf(1); b = fbuf(2); c = fbuf(3)
   alpha = fbuf(4); beta = fbuf(5); gamma = fbuf(6)
   nx = ibuf(1); ny = ibuf(2); nz = ibuf(3)
!
   nt = nx*ny*nz
   nt_fftw = nx*ny*2*(nz/2+1)
   allocate(func(nt), ft(ny*nz), ftm(nt))
   allocate(rvec(1:3,1:nt), kvec(1:3,1:nt))
   vol = a*b*c
   vfac = nt
!
   ifc = nx/2+1
   jfc = ny/2+1
   kfc = nz/2+1
!
   idr = 0
   dx = a/real(nx,kind=RealKind)
   dy = b/real(ny,kind=RealKind)
   dz = c/real(nz,kind=RealKind)
   do k = 1, nz
      z = (k-1)*dz
      do j = 1, ny
         y = (j-1)*dy
         do i = 1, nx
            x = (i-1)*dx
            idr = idr + 1
            rvec(1,idr) = x
            rvec(2,idr) = y
            rvec(3,idr) = z
            func(idr) = exp(alpha*x+beta*y+gamma*z)
         enddo
      enddo
   enddo
!
   call initTimer()
   t0 = getTime()
!
!  ====================================================================
!  Hand coded version
!  ====================================================================
   if (nt < man_limit) then
      idk = 0
      ka = PI2/a
      kb = PI2/b
      kc = PI2/c
      ftm(1:nt) = CZERO
      do k = 1, nz
         if (nz == 1) then
            kz = ZERO
         else
            kz = (k-kfc)*kc
         endif
         do j = 1, ny
            if (ny == 1) then
               ky = ZERO
            else
               ky = (j-jfc)*kb
            endif
            do i = 1, nx
               if (nx == 1) then
                  kx = ZERO
               else
                  kx = (i-ifc)*ka
               endif
               idk = idk + 1
               kvec(1,idk) = kx
               kvec(2,idk) = ky
               kvec(3,idk) = kz
               do idr = 1, nt
                  expfac = SQRTm1*(kx*rvec(1,idr)+ky*rvec(2,idr)+kz*rvec(3,idr))
                  ftm(idk) = ftm(idk) + func(idr)*exp(expfac)
               enddo
            enddo
         enddo
      enddo
      if (MyPE == 0) then
         write(6,'('' time for hand coded ft ='',1f10.5)')getTime()-t0
      endif
      ftm(1:nt) = ftm(1:nt)/vfac
   endif
!
   t0 = getTime()
!  ====================================================================
!  Numerical Recipes version 2.0 code
!  --------------------------------------------------------------------
   call rlft3(func,ft,nx,ny,nz,1)
!  --------------------------------------------------------------------
   func(1:nt) = func(1:nt)/vfac
   ft(1:ny*nz) = ft(1:ny*nz)/vfac
!
   allocate(ft_a(nt), ft_nr(nt), kvec_nr(1:3,1:nt))
!
!  ===================================================================
!  determine the k-vector for the Numerical Recipes rlft3 routine
!  ===================================================================
   idk = 0; idyz = 0
   ka = PI2/a
   kb = PI2/b
   kc = PI2/c
   do k = 1, nz
      if (k < kfc .or. nz == 1) then
         kz = (k-1)*kc
      else
         kz = (k-1-nz)*kc
      endif
!
      if (k == 1) then
         kconj = 1
      else
         kconj = nz + 2 - k
      endif
!
      if (abs(gamma) < TEN2m6) then
         if (k == 1) then
            ftz = c
         else
            ftz = CZERO
         endif
      else
         ftz = (exp(gamma*c)-CONE)/cmplx(gamma,kz,kind=CmplxKind)
      endif
!
      do j = 1, ny
         if (j < jfc .or. ny == 1) then
            ky = (j-1)*kb
         else
            ky = (j-1-ny)*kb
         endif
!
         if (j == 1) then
            jconj = 1
         else
            jconj = ny + 2 - j
         endif
!
         if (abs(beta) < TEN2m6) then
            if (j == 1) then
               fty = b
            else
               fty = CZERO
            endif
         else
            fty = (exp(beta*b)-CONE)/cmplx(beta,ky,kind=CmplxKind)
         endif
!
         idyz = idyz + 1
         do i = 1, nx
            if (i < ifc .or. nx == 1) then
               kx = (i-1)*ka
            else
               kx = (i-1-nx)*ka
            endif
            if (abs(alpha) < TEN2m6) then
               if (i == 1) then
                  ftx = a
               else
                  ftx = CZERO
               endif
            else
               ftx = (exp(alpha*a)-CONE)/cmplx(alpha,kx,kind=CmplxKind)
            endif
            idk = idk + 1
            kvec_nr(1,idk) = kx
            kvec_nr(2,idk) = ky
            kvec_nr(3,idk) = kz
            if (i == ifc) then
               ft_nr(idk) = ft(idyz)
            else if (i < ifc) then
               idm = 2*i + (j-1)*nx + (k-1)*ny*nx
               ft_nr(idk) = cmplx(func(idm-1),func(idm),kind=CmplxKind)
            else
               idm = 2*(nx+2-i) + (jconj-1)*nx + (kconj-1)*ny*nx
               ft_nr(idk) = cmplx(func(idm-1),-func(idm),kind=CmplxKind)
            endif
            ft_a(idk) = ftx*fty*ftz/vol  ! ft_a contains the analytical results
         enddo
      enddo
   enddo
!
   if (MyPE == 0) then
      write(6,'('' time for calling rlft3 ='',1f10.5)')getTime()-t0
   endif
!
   if (nt < man_limit .and. MyPE == 0) then
      do idm = 1, nt
         found = .false.
         LOOP_idk: do idk = 1, nt
            if ( abs(kvec_nr(1,idk)-kvec(1,idm)) < TEN2m6 .and.          &
                 abs(kvec_nr(2,idk)-kvec(2,idm)) < TEN2m6 .and.          &
                 abs(kvec_nr(3,idk)-kvec(3,idm)) < TEN2m6 ) then
               if (abs(ftm(idm)-ft_nr(idk)) > TEN2m6) then
                  write(6,'(a)')'WARNING: abs(ftm-ft_nr) >> 0'
                  write(6,'(3f10.5,3(2x,2d15.8))')kvec_nr(1:3,idk),      &
                                                  ftm(idm), ft_nr(idk), ft_a(idk)
                  stop 'Error'
!              else
!                 write(6,'(3f10.5,3(2x,2d15.8))')kvec_nr(1:3,idk),      &
!                                                 ftm(idm), ft_nr(idk), ft_a(idk)
               endif
               found = .true.
               exit LOOP_idk
            endif
         enddo LOOP_idk
         if (.not. found) then
            write(6,'(a)')'k-vectors can not be matched!'
            stop 'ERROR'
         endif
      enddo
      write(6,'(a)')'Testing rlft3 is successful!'
   endif
!  -------------------------------------------------------------------
   call syncAllPEs()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Start testing ParallelFFT
!  ===================================================================
   if (mod(ny*nz,NumPEs) /= 0) then
      call ErrorHandler('ft_p3dfft','ny*nz is not divisable by NumPEs',ny*nz,NumPEs)
   endif
!
!  ===================================================================
!  Set up uniform 3D grid
!  ===================================================================
   cell = ZERO
   cell(1,1) = a
   cell(2,2) = b
   cell(3,3) = c
!  -------------------------------------------------------------------
   t0 = getTime()
   call initUniform3DGrid('main',0)
   if (MyPE == 0) then
      write(6,'('' time for calling initUniform3DGrid ='',1f10.5)')getTime()-t0
   endif
   t0 = getTime()
   call createUniform3DGrid('FFT', nx, ny, nz, cell)
   if (MyPE == 0) then
      write(6,'('' time for calling createUniform3DGrid ='',1f10.5)')getTime()-t0
   endif
!  ===================================================================
!  Initialize ParallelFFT module
!  -------------------------------------------------------------------
   call initParallelFFT()
!  -------------------------------------------------------------------
   t0 = getTime()
   call createProcessorMesh('FFT', getProcessorMesh())
   if (MyPE == 0) then
      write(6,'('' time for calling createProcessorMesh ='',1f10.5)')getTime()-t0
   endif
!  -------------------------------------------------------------------
   if (MyPE == 0) then
!     ----------------------------------------------------------------
      call printUniform3DGrid('FFT')
!     ----------------------------------------------------------------
   endif
!
   t0 = getTime()
   if (MyPE == 0) then
      write(6,'('' time for calling initParallelFFT ='',1f10.5)')getTime()-t0
   endif
!  -------------------------------------------------------------------
   call allocateFunctionSpace(rfunc_fftwr,memsize)
   p_rfunc_fftwr => aliasFunctionSpace(rfunc_fftwr,m1,m2,m3)
   n1 = getNumLocalGridPoints('R',1)
   n2 = getNumLocalGridPoints('R',2)
   n3 = getNumLocalGridPoints('R',3)
   gir = getGlobalGridIndexRange('R')
!  -------------------------------------------------------------------
   grid_start(1:3) = gir(1:3,1)
   grid_end(1:3) = gir(1:3,2)
!  -------------------------------------------------------------------
   call distributeUniformGrid('FFT',grid_start,grid_end)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize rfunc_fftwr, the array to be transformed
!  ===================================================================
   t0 = getTime()
   rfunc_fftwr = ZERO
   idk = 0
   do k = 1, n3
      do j = 1, n2
         do i = 1, n1
            idk = idk + 1
            p = getGridPointCoord('R',i,j,k)
!p = getGridPointCoord('R',idk)
!write(6,'(a,3f12.5,4x,3f12.5)')'r_fft,r = ',p(1:3),rvec(1:3,idk)
            p_rfunc_fftwr(i,j,k) = exp(alpha*p(1)+beta*p(2)+gamma*p(3))
         enddo
      enddo
   enddo
   if (MyPE == 0) then
      write(6,'('' time for initialize the real function ='',1f10.5)')getTime()-t0
   endif
!
!  -------------------------------------------------------------------
   call allocateFunctionSpace(ft_fftw_ip,nump_local)
!  -------------------------------------------------------------------
!
   t0 = getTime()
!  ===================================================================
!  transform from physical space to wavenumber space
!  Its result, ft_fftw_ip, has already been multiplied by 1/ng
!  -------------------------------------------------------------------
   call performTransformR2C(rfunc_fftwr,ft_fftw_ip)
!  -------------------------------------------------------------------
   if (MyPE == 0) then
      write(6,'('' time for in-place fftw real ='',1f10.5)')getTime()-t0
!
!     do idk = 1, nump_local
!        kvec_fftw = getGridPointCoord('K',idk)
!        write(6,'(a,3f12.5,4x,3f12.5)')'kvec_fft,kvec = ',kvec_fftw(1:3),kvec(1:3,idk)
!     enddo
   endif
!
   call syncAllPEs()
! 
   if (nt < man_limit) then
      do idk = 1, nump_local
         found = .false.
         kvec_fftw = getGridPointCoord('K',idk)
         LOOP_idm_para: do idm = 1, nt
            if ( abs(kvec_fftw(1)-kvec_nr(1,idm)) < TEN2m6 .and.          &
                 abs(kvec_fftw(2)-kvec_nr(2,idm)) < TEN2m6 .and.          &
                 abs(kvec_fftw(3)-kvec_nr(3,idm)) < TEN2m6 ) then
               if (abs(ft_nr(idm)-ft_fftw_ip(idk)) > TEN2m6) then
                  write(6,'(a)')'ERROR: bad transformation'
                  write(6,'(2i4,2x,3f10.5,2(2x,2d15.8))')MyPE,idk,kvec_nr(1:3,idm),   &
                                                         ft_nr(idm), ft_fftw_ip(idk)
!                 stop 'Error'
!              else
!                 write(6,'(i4,2x,3f10.5,2(2x,2d15.8))')MyPE,kvec_nr(1:3,idm),      &
!                                                       ft_nr(idm), ft_fftw_ip(idk)
               endif
               found = .true.
               exit LOOP_idm_para
            endif
         enddo LOOP_idm_para
         if (.not. found) then
            write(6,'(a,2i4,2x,3f10.5)')'The parallel k-vector can not be matched:',MyPE,idk,kvec_fftw(1:3)
            stop 'Error!'
         endif
      enddo
      do idm = 1, nt
         found = .false.
         LOOP_idk_serial: do idk = 1, nump_local
            kvec_fftw = getGridPointCoord('K',idk)
            if ( abs(kvec_fftw(1)-kvec_nr(1,idm)) < TEN2m6 .and.          &
                 abs(kvec_fftw(2)-kvec_nr(2,idm)) < TEN2m6 .and.          &
                 abs(kvec_fftw(3)-kvec_nr(3,idm)) < TEN2m6 ) then
               if (abs(ft_nr(idm)-ft_fftw_ip(idk)) > TEN2m6) then
                  write(6,'(a)')'ERROR: bad transformation'
                  write(6,'(2i4,2x,3f10.5,2(2x,2d15.8))')MyPE,idk,kvec_nr(1:3,idm),   &
                                                         ft_nr(idm), ft_fftw_ip(idk)
!                 stop 'Error'
!              else
!                 write(6,'(i4,2x,3f10.5,2(2x,2d15.8))')MyPE,kvec_nr(1:3,idm),      &
!                                                       ft_nr(idm), ft_fftw_ip(idk)
               endif
               found = .true.
               exit LOOP_idk_serial
            endif
         enddo LOOP_idk_serial
         if (.not. found) then
            k = 0
         else 
            k = 1
         endif
         call GlobalSum(k)
         if (k == 0) then
            write(6,'(a,2i4,2x,3f10.5)')'The serial k-vector can not be matched:',MyPE,idm,kvec_nr(1:3,idm)
            stop 'Error!'
         endif
      enddo
   endif
!
   deallocate( ft_fftw_ip, rfunc_fftwr )
!
!  ===================================================================
!  Free work space
!  -------------------------------------------------------------------
   call endParallelFFT()
   call endUniform3DGrid()
!  -------------------------------------------------------------------
!
   deallocate(func, ft, ftm, rvec, kvec)
   deallocate(kvec_nr, ft_a, ft_nr)
!
   call endMPP()
!
   stop 'Ok'
!
end program testParallelFFT
