program ft_mpi_fftw
!  ********************************************************************
!  test MPI FFTW parallel fast Fourier transform routines using Fortran 2003
!
!       It appears that MPI FFTW is problematic.
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
   use iso_c_binding
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use TimerModule, only : initTimer, getTime
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, CONE, PI2, SQRTm1
!
   use MPPModule, only : initMPP, endMPP, getCommunicator, bcastMessage
   use MPPModule, only : MyPE, NumPEs, syncAllPEs
!
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit   none
!
!  Note: fftw3.f03 and fftw3-mpi.f03 can not be used together.
include 'mpif.h'
#ifdef FFTW
   include 'fftw3-mpi.f03'
#endif
!
   logical :: found
!
   type (C_PTR) :: plan
!
   integer (C_INTPTR_T) :: c_nx, c_ny, c_nz, c_nt
   integer (C_INTPTR_T) :: c_local_n, c_local_start, c_alloc_local
!
   type (C_PTR) :: c_func_fftwr, c_cfunc_fftwr
   real (kind=C_DOUBLE), pointer :: rfunc_fftwr(:,:,:)
   complex (kind=C_DOUBLE_COMPLEX), pointer :: cfunc_fftwr(:,:,:)
!
   integer (kind=IntKind), parameter :: man_limit = 64*64*64
!
   integer (kind=IntKind) :: local_n, local_start, nump_local
   integer (kind=IntKind) :: nx, ny, nz, nt, nt_fftw
!
   integer (kind=IntKind) :: i, j, k, idk, idr, idyz, idm
   integer (kind=IntKind) :: ifc, jfc, kfc
   integer (kind=IntKind) :: jconj, kconj
   integer (kind=IntKind) :: comm
   integer (kind=IntKind) :: ibuf(3)
!
   real (kind=RealKind) :: t0, vfac, vol
   real (kind=RealKind) :: a, b, c
   real (kind=RealKind) :: x, y, z, dx, dy, dz
   real (kind=RealKind) :: kx, ky, kz
   real (kind=RealKind) :: ka, kb, kc
   real (kind=RealKind) :: alpha, beta, gamma, fbuf(6)
   real (kind=RealKind), allocatable :: func(:), kvec(:,:), rvec(:,:)
   real (kind=RealKind), allocatable :: kvec_nr(:,:), kvec_fftw(:,:)
!
   complex (kind=CmplxKind) :: ftx, fty, ftz, expfac
   complex (kind=CmplxKind), allocatable :: ft(:)
   complex (kind=CmplxKind), allocatable :: ftm(:)
   complex (kind=CmplxKind), allocatable :: ft_a(:), ft_nr(:)
   complex (kind=CmplxKind), allocatable :: ft_fftw_ip(:)
!
!  -------------------------------------------------------------------
   call initMPP()
!  -------------------------------------------------------------------
!
   if (MyPE == 0) then
      write(6,'(a)')'Data type checking...'
      write(6,'(4x,2(a,i2))')'C_INT            = ',C_INT,           ';   IntKind = ',IntKind
      write(6,'(4x,2(a,i2))')'C_INTPTR_T       = ',C_INTPTR_T,      ';   IntKind = ',IntKind
      write(6,'(4x,2(a,i2))')'C_DOUBLE         = ',C_DOUBLE  ,      ';  RealKind = ',RealKind
      write(6,'(4x,2(a,i2))')'C_DOUBLE_COMPLEX = ',C_DOUBLE_COMPLEX,'; CmplxKind = ',CmplxKind
      write(6,'(/)')
!
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
   if (MyPE == 0) then
      write(6,'('' time for calling rlft3 ='',1f10.5)')getTime()-t0
   endif
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
#ifdef FFTW
!  ===================================================================
!  Testing MPI-FFTW.....
!  -------------------------------------------------------------------
   call fftw_mpi_init()
!  -------------------------------------------------------------------
!
   comm = getCommunicator()
   c_nx = nx; c_ny = ny; c_nz = nz; c_nt = nt
!
!  ===================================================================
!  get local data size and allocate
!  -------------------------------------------------------------------
!  c_alloc_local = fftw_mpi_local_size_3d(c_nz, c_ny, c_nx/2+1, MPI_COMM_WORLD, c_local_n, c_local_start)
   c_alloc_local = fftw_mpi_local_size_3d(c_nz, c_ny, c_nx, MPI_COMM_WORLD, c_local_n, c_local_start)
   local_n = c_local_n; local_start = c_local_start
!
   c_func_fftwr = fftw_alloc_real(2*c_alloc_local)
   c_cfunc_fftwr = fftw_alloc_complex(c_alloc_local)
   call c_f_pointer(c_cfunc_fftwr, cfunc_fftwr, [c_nx,c_ny,c_local_n])
!  call c_f_pointer(c_func_fftwr, rfunc_fftwr, [(c_nx+2),c_ny,c_local_n])
!  call c_f_pointer(c_cfunc_fftwr, cfunc_fftwr, [(c_nx/2+1),c_ny,c_local_n])
!  call c_f_pointer(c_func_fftwr, cfunc_fftwr, [c_nx/2+1,c_ny,c_local_n])
!  -------------------------------------------------------------------
   
!  ===================================================================
!  create plan for in-place r2c DFT
!  -------------------------------------------------------------------
!  plan = fftw_mpi_plan_dft_r2c_3d(c_nz, c_ny, c_nx, rfunc_fftwr, cfunc_fftwr, MPI_COMM_WORLD, FFTW_MEASURE)
   plan = fftw_mpi_plan_dft_3d(c_nz, c_ny, c_nx, cfunc_fftwr, cfunc_fftwr, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize rfunc_fftwr to some function
!  ===================================================================
   do k = 1, local_n
      z = (local_start+k-1)*dz
      do j = 1, ny
         y = (j-1)*dy
         do i = 1, nx
            x = (i-1)*dx
            idm = i + (j-1)*(nx+2) + (k-1)*ny*(nx+2)         ! beware of the padding space
!           rfunc_fftwr(idm) = exp(alpha*x+beta*y+gamma*z) ! used for in-place transformation
!           rfunc_fftwr(i,j,k) = exp(alpha*x+beta*y+gamma*z)
            cfunc_fftwr(i,j,k) = exp(alpha*x+beta*y+gamma*z)
         enddo
      enddo
   enddo
print *,'MyPE, local_start, local_n = ',MyPE, local_start, local_n
call syncAllPEs()
!  ===================================================================
!
!  ===================================================================
!  compute transforms as many times as desired
!  -------------------------------------------------------------------
   t0 = getTime()
   call fftw_mpi_execute_dft(plan,cfunc_fftwr,cfunc_fftwr)
!  call fftw_mpi_execute_dft_r2c(plan,rfunc_fftwr,cfunc_fftwr)
   call fftw_destroy_plan(plan)
   if (MyPE == 0) then
      write(6,'('' time for in-place fftw real ='',1f10.5)')getTime()-t0
   endif
!  -------------------------------------------------------------------
   rfunc_fftwr = rfunc_fftwr/vfac
!
!  -------------------------------------------------------------------
   call syncAllPEs()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  determine the k-vector for the MPI-FFTW r2c routine and check results
!  ===================================================================
   nump_local = nx*ny*local_n
   allocate( kvec_fftw(3,nump_local) , ft_fftw_ip(nump_local) )
!
   ka = PI2/a
   kb = PI2/b
   kc = PI2/c
   idk = 0
   do k = 1, local_n
      if (k < kfc .or. nz == 1) then
         kz = (k+local_start-1)*kc
      else
         kz = (k+local_start-1-nz)*kc
      endif
!
      if (k == 1) then
         kconj = 1
      else
         kconj = local_n + 2 - k
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
         do i = 1, nx
            if (i < ifc .or. nx == 1) then
               kx = (i-1)*ka
            else
               kx = (i-1-nx)*ka
            endif
            idk = idk + 1
            kvec_fftw(1,idk) = kx
            kvec_fftw(2,idk) = ky
            kvec_fftw(3,idk) = kz
if (.false.) then
!           ==========================================================
!           for the in-place case
!           ==========================================================
            if (i <= ifc) then
               idm = 2*i + (j-1)*(nx+2) + (k-1)*ny*(nx+2)
!              ft_fftw_ip(idk) = cmplx(rfunc_fftwr(idm-1), -rfunc_fftwr(idm), &
               ft_fftw_ip(idk) = cmplx(rfunc_fftwr(2*i-1,j,k), -rfunc_fftwr(2*i,j,k), &
                                       kind=CmplxKind)
            else
               idm = 2*(nx+2-i) + (jconj-1)*(nx+2) + (kconj-1)*ny*(nx+2)
!              ft_fftw_ip(idk) = cmplx(rfunc_fftwr(idm-1), rfunc_fftwr(idm),  &
               ft_fftw_ip(idk) = cmplx(rfunc_fftwr(2*(nx+2-i)-1,j,k), rfunc_fftwr(2*(nx+2-i),j,k),  &
                                       kind=CmplxKind)
            endif
endif
         enddo
      enddo
   enddo
!
!  do idk = 1, nump_local
   idk = 0
   do k = 1, local_n
   do j = 1, ny
   do i = 1, nx
      idk = idk + 1
      found = .false.
      LOOP_idm_para: do idm = 1, nt
         if ( abs(kvec_fftw(1,idk)-kvec_nr(1,idm)) < TEN2m6 .and.          &
              abs(kvec_fftw(2,idk)-kvec_nr(2,idm)) < TEN2m6 .and.          &
              abs(kvec_fftw(3,idk)-kvec_nr(3,idm)) < TEN2m6 ) then
!           if (abs(ftm(idm)-ft_fftw_ip(idk)) > TEN2m6) then
        !   if (abs(ftm(idm)-cfunc_fftwr(i,j,k)) > TEN2m6) then
        !      write(6,'(a)')'ERROR: bad transformation'
        !      write(6,'(3f10.5,2(2x,2d15.8))')kvec_nr(1:3,idm),      &
        !                                      ftm(idm), cfunc_fftwr(i,j,k)
!       !                                      ftm(idm), ft_fftw_ip(idk)
        !      stop 'Error'
        !   else
               write(6,'(3f10.5,2(2x,2d15.8))')kvec_nr(1:3,idm),      &
                                               ftm(idm), cfunc_fftwr(i,j,k)
!                                              ftm(idm), ft_fftw_ip(idk)
        !   endif
            found = .true.
            exit LOOP_idm_para
         endif
      enddo LOOP_idm_para
      if (.not. found) then
         write(6,'(a)')'k-vectors can not be matched in MPI FFTW case!'
         stop 'ERROR'
      endif
   enddo
   enddo
   enddo
!
   deallocate( kvec_fftw, ft_fftw_ip )
   nullify( rfunc_fftwr, cfunc_fftwr)
!  ===================================================================
!
!  ===================================================================
!  Clean up the MPI-FFTW space
!  -------------------------------------------------------------------
   call fftw_mpi_cleanup()
   call fftw_free(c_func_fftwr)
   call fftw_free(c_cfunc_fftwr)
!  -------------------------------------------------------------------
#endif
!
   deallocate(func, ft, ftm, rvec, kvec)
   deallocate(kvec_nr, ft_a, ft_nr)
!
   call endMPP()
!
   stop 'Ok'
!
end program ft_mpi_fftw
