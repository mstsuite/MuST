program fourier_transform
!  ********************************************************************
!  test Fourier transform routines
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
   implicit   none
!
#ifdef FFTW
   include 'fftw3.f'
#endif
!
   logical :: found
!
   integer (kind=8) :: plan
   integer (kind=IntKind) :: nx, ny, nz, nt, nt_fftw
!
   integer (kind=IntKind) :: i, j, k, idk, idr, idyz, idm
   integer (kind=IntKind) :: ifc, jfc, kfc
   integer (kind=IntKind) :: jconj, kconj
!
   real (kind=RealKind) :: t0, vfac, vol
   real (kind=RealKind) :: a, b, c
   real (kind=RealKind) :: x, y, z, dx, dy, dz
   real (kind=RealKind) :: kx, ky, kz
   real (kind=RealKind) :: ka, kb, kc
   real (kind=RealKind) :: alpha, beta, gamma
   real (kind=RealKind), allocatable :: func(:), kvec(:,:), rvec(:,:)
   real (kind=RealKind), allocatable :: kvec_nr(:,:), kvec_fftw(:,:)
   real (kind=RealKind), allocatable :: func_fftwr(:)
!
   complex (kind=CmplxKind) :: ftx, fty, ftz, expfac
   complex (kind=CmplxKind), allocatable :: ft(:)
   complex (kind=CmplxKind), allocatable :: ftm(:)
   complex (kind=CmplxKind), allocatable :: ft_fftwc(:), ft_fftwr(:)
   complex (kind=CmplxKind), allocatable :: ft_a(:), ft_nr(:), ft_fftw(:)
   complex (kind=CmplxKind), allocatable :: ft_fftw_ip(:)
!
   a = ZERO; b = ZERO; z = ZERO
   do while (a <= TEN2m6 .or. b <= TEN2m6 .or. c <= TEN2m6)
      write(6,'(1x,a,$)')'Box size in x, y, z direction: '
      read(5,*)a, b, c
   enddo
   write(6,'(f10.5,'','',f10.5,'','',f10.5)')a, b, c
!
   nx = 0; ny = 0; nz = 0
   do while (nx < 1 .or. ny < 1 .or. nz < 1)
      write(6,'(1x,a,$)')'Number of mesh points in x, y, z direction: '
      read(5,*)nx, ny, nz
   enddo
   write(6,'(i5,'','',i5,'','',i5)')nx, ny, nz
!
   write(6,'(/,1x,a)')    &
    'This Code Will Perform Fast Fourier Transform of an Exponential Function'
   write(6,'(//,1x,a)')'                 exp(alpha*x + beta*y + gamma*z)'
!
   alpha = ZERO; beta = ZERO; gamma = ZERO
   write(6,'(//,1x,a,$)')'Enter alpha, beta, gamma: '
   read(5,*)alpha, beta, gamma
   write(6,'(f10.5,'','',f10.5,'','',f10.5)')alpha, beta, gamma
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
   write(6,'('' time for hand coded ft ='',1f10.5)')getTime()-t0
   ftm(1:nt) = ftm(1:nt)/vfac
!
!  do idk = 1, nt
!     write(6,'(3f10.5,4x,2d15.8)')kvec(1:3,idk),ftm(idk)
!  enddo
!
   t0 = getTime()
!  ====================================================================
!  Numerical Recipes version 2.0 code
!  --------------------------------------------------------------------
   call rlft3(func,ft,nx,ny,nz,1)
!  --------------------------------------------------------------------
   write(6,'('' time for calling rlft3 ='',1f10.5)')getTime()-t0
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
            ft_a(idk) = ftx*fty*ftz/vol
         enddo
      enddo
   enddo
!
#ifdef FFTW
!  ===================================================================
!  FFTW version 3.x code
!  ===================================================================
!
   allocate( kvec_fftw(3,nt), ft_fftw(nt), ft_fftwc(nt) )
   allocate( func_fftwr((nx+2)*ny*nz), ft_fftw_ip(nt) )
   idr = 0
   do k = 1, nz
      z = (k-1)*dz
      do j = 1, ny
         y = (j-1)*dy
         do i = 1, nx
            x = (i-1)*dx
            idr = idr + 1
            func(idr) = exp(alpha*x+beta*y+gamma*z)
            ft_fftwc(idr) = cmplx(func(idr),ZERO,kind=CmplxKind)
!
            idm = i + (j-1)*(nx+2) + (k-1)*ny*(nx+2)      ! beware of the padding space
            func_fftwr(idm) = exp(alpha*x+beta*y+gamma*z) ! used for in-place transformation
         enddo
      enddo
   enddo
!
   t0 = getTime()
!  -------------------------------------------------------------------
   call dfftw_plan_dft_3d( plan, nx, ny, nz, ft_fftwc, ft_fftwc,      &
                           FFTW_BACKWARD, FFTW_ESTIMATE )
   call dfftw_execute( plan )
   call dfftw_destroy_plan( plan )
   ft_fftwc(1:nt) = ft_fftwc(1:nt)/vfac
!  -------------------------------------------------------------------
   write(6,'('' time for fftw cmplx ='',1f10.5)')getTime()-t0
!
   allocate( ft_fftwr((nx/2+1)*ny*nz) )
!
   t0 = getTime()
!  -------------------------------------------------------------------
   call dfftw_plan_dft_r2c_3d( plan, nx, ny, nz, func, ft_fftwr,      &
                                 FFTW_ESTIMATE )
   call dfftw_execute( plan )
   call dfftw_destroy_plan( plan )
   ft_fftwr(:) = ft_fftwr(:)/vfac
!  -------------------------------------------------------------------
   write(6,'('' time for fftw real ='',1f10.5)')getTime()-t0
!
   t0 = getTime()
!  -------------------------------------------------------------------
   call dfftw_plan_dft_r2c_3d( plan, nx, ny, nz, func_fftwr, func_fftwr, &
                               FFTW_ESTIMATE )
   call dfftw_execute( plan )
   call dfftw_destroy_plan( plan )
   func_fftwr(:) = func_fftwr(:)/vfac
!  -------------------------------------------------------------------
   write(6,'('' time for in-place fftw real ='',1f10.5)')getTime()-t0
!
!  ===================================================================
!
!  ===================================================================
!  determine the k-vector for the fftw r2c routine
!  ===================================================================
   idk = 0
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
            if (i <= ifc) then
               idm = i + (j-1)*(nx/2+1) + (k-1)*ny*(nx/2+1)
               ft_fftw(idk) = conjg(ft_fftwr(idm))
!              =======================================================
!                         for the in-place case
!              =======================================================
               idm = 2*i + (j-1)*(nx+2) + (k-1)*ny*(nx+2)
               ft_fftw_ip(idk) = cmplx(func_fftwr(idm-1), -func_fftwr(idm), &
                                       kind=CmplxKind)
            else
               idm = (nx+2-i) + (jconj-1)*(nx/2+1) + (kconj-1)*ny*(nx/2+1)
               ft_fftw(idk) = ft_fftwr(idm)
!              =======================================================
!                         for the in-place case
!              =======================================================
               idm = 2*(nx+2-i) + (jconj-1)*(nx+2) + (kconj-1)*ny*(nx+2)
               ft_fftw_ip(idk) = cmplx(func_fftwr(idm-1), func_fftwr(idm),  &
                                       kind=CmplxKind)
            endif
         enddo
      enddo
   enddo
!
   LOOP_idm_w: do idm = 1, nt
      found = .false.
      LOOP_idk_w: do idk = 1, nt
         if ( abs(kvec(1,idk)-kvec_nr(1,idm)) < TEN2m6 .and.          &
              abs(kvec(2,idk)-kvec_nr(2,idm)) < TEN2m6 .and.          &
              abs(kvec(3,idk)-kvec_nr(3,idm)) < TEN2m6 ) then
            if (abs(ftm(idk)-ft_nr(idm)) > TEN2m6  .or.               &
                abs(ftm(idk)-ft_fftwc(idm)) > TEN2m6 .or.             &
                abs(ftm(idk)-ft_fftw(idm)) > TEN2m6 .or.              &
                abs(ftm(idk)-ft_fftw_ip(idm)) > TEN2m6) then
               write(6,'(a)')'ERROR: bad transformation'
               write(6,'(3f10.5,3(2x,2d15.8))')kvec_nr(1:3,idm),      &
                                               ftm(idk), ft_nr(idm),  &
                                               ft_fftw_ip(idm)
               stop 'Error'
            else
               write(6,'(3f10.5,3(2x,2d15.8))')kvec_nr(1:3,idm),      &
                                               ftm(idk), ft_nr(idm),  &
                                               ft_fftw_ip(idm)
            endif
            found = .true.
            exit LOOP_idk_w
         endif
      enddo LOOP_idk_w
      if (.not. found) then
         write(6,'(a)')'k-vectors can not be matched!'
         stop 'ERROR'
      endif
   enddo LOOP_idm_w
!
   deallocate( kvec_fftw, ft_fftw, ft_fftwr, ft_fftwc, func_fftwr, ft_fftw_ip )
#endif
!
   write(6,'(/)')
   LOOP_idm: do idm = 1, nt
      found = .false.
      LOOP_idk: do idk = 1, nt
         if ( abs(kvec(1,idk)-kvec_nr(1,idm)) < TEN2m6 .and.          &
              abs(kvec(2,idk)-kvec_nr(2,idm)) < TEN2m6 .and.          &
              abs(kvec(3,idk)-kvec_nr(3,idm)) < TEN2m6 ) then
            if (abs(ftm(idk)-ft_nr(idm)) > TEN2m6) then
               write(6,'(a)')'WARNING: abs(ftm-ft_nr) >> 0'
               write(6,'(3f10.5,3(2x,2d15.8))')kvec_nr(1:3,idm),      &
                                               ftm(idk), ft_nr(idm), ft_a(idm)
               stop 'Error'
            else
               write(6,'(3f10.5,3(2x,2d15.8))')kvec_nr(1:3,idm),      &
                                               ftm(idk), ft_nr(idm), ft_a(idm)
            endif
            found = .true.
            exit LOOP_idk
         endif
      enddo LOOP_idk
      if (.not. found) then
         write(6,'(a)')'k-vectors can not be matched!'
         stop 'ERROR'
      endif
   enddo LOOP_idm
!
   deallocate(func, ft, ftm, rvec, kvec)
   deallocate(kvec_nr, ft_a, ft_nr)
!
   stop 'Ok'
!
end program fourier_transform
