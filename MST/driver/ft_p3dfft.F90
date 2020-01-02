program ft_p3dfft
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
#ifdef P3DFFT
   use p3dfft
#endif
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
   integer (kind=IntKind) :: dims(2)
   integer (kind=IntKind) :: iproc, jproc, ig, jg, kg
   integer (kind=IntKind) :: istart(3),iend(3),isize(3)
   integer (kind=IntKind) :: fstart(3),fend(3),fsize(3)
   integer (kind=IntKind) :: tstart(3),tend(3),memsize(3)
   real (kind=RealKind), allocatable :: rfunc_fftwr(:)
   real (kind=RealKind), pointer :: p_rfunc_fftwr(:,:,:)
!  ===================================================================
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
#ifdef P3DFFT
!  ===================================================================
!  Testing P3DFFT
!  ===================================================================
   if (mod(ny*nz,NumPEs) /= 0) then
      call ErrorHandler('ft_p3dfft','ny*nz is not divisable by NumPEs',ny*nz,NumPEs)
   endif
!
!  ===================================================================
!  NumPEs is divided into a iproc x jproc stencle
!  ===================================================================
   factors => getSubFactors(NumPEs,2,m) ! break the total number of processes into two factors
   k = ny + nz + 1
   iproc = factors(1,1)
   jproc = factors(2,1)
   found = .false.
   do i = 1, m
      if (mod(ny,factors(1,i)) == 0 .and. mod(nz,factors(2,i)) == 0 .and.  &
          mod(ny,factors(2,i)) == 0 .and. mod(nx/2,factors(1,i)) == 0) then
         j = (nx/2+ny)/factors(1,i) + (ny+nz)/factors(2,i)
         if (j < k .and. factors(1,i) <= factors(2,i)) then
            iproc = factors(1,i)
            jproc = factors(2,i)
            found = .true.
            k = j
         endif
      endif
   enddo
!
   if (.not.found) then
      write(6,'(a,2i5)')'A parallel processor grid is not established for the given real space mesh',iproc,jproc
   endif
!
   if (MyPE == 0) then
      write(6,'(a,2i5)')'iproc, jproc = ',iproc, jproc
   endif
!
   dims(1) = iproc
   dims(2) = jproc
!
!  ===================================================================
!  Set up work structures for P3DFFT
!  -------------------------------------------------------------------
   comm = getCommunicator()
   call p3dfft_setup(dims,nx,ny,nz,comm)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Get dimensions for the original array of real numbers, X-pencils
!  -------------------------------------------------------------------
   call p3dfft_get_dims(istart,iend,isize,1)
!  -------------------------------------------------------------------
   write(6,'(a,10i5)')'MyPE, istart, iend, isize   = ',MyPE,istart(1),iend(1),isize(1), &
                                                            istart(2),iend(2),isize(2), &
                                                            istart(3),iend(3),isize(3)
!
!  ===================================================================
!  Get dimensions for the R2C-forward-transformed array of complex numbers
!  Z-pencils (depending on how the library was compiled, the first
!  dimension could be either X or Z)
!  -------------------------------------------------------------------
   call p3dfft_get_dims(fstart,fend,fsize,2)
!  -------------------------------------------------------------------
   write(6,'(a,10i5)')'MyPE, fstart, fend, fsize   = ',MyPE,fstart(1),fend(1),fsize(1), &
                                                            fstart(2),fend(2),fsize(2), &
                                                            fstart(3),fend(3),fsize(3)
   if (fend(1) < fstart(1) .or. fend(2) < fstart(2) .or. fend(3) < fstart(3)) then
      call ErrorHandler('main','fend < fstart')
   endif
!
!  ===================================================================
!  Since we are allocating the same array for input and output,
!  we need to make sure it has enough space. This is achieved
!  by calling get_dims with option 3 to get memsize and use memsize
!  to allocate space.
!  -------------------------------------------------------------------
   call p3dfft_get_dims(tstart,tend,memsize,3)
!  -------------------------------------------------------------------
   write(6,'(a,10i5)')'MyPE, tstart, tend, memsize = ',MyPE,tstart(1),tend(1),memsize(1), &
                                                            tstart(2),tend(2),memsize(2), &
                                                            tstart(3),tend(3),memsize(3)
!
   allocate(rfunc_fftwr(memsize(1)*memsize(2)*memsize(3)))
   p_rfunc_fftwr => aliasArray3_r(rfunc_fftwr,memsize(1),memsize(2),memsize(3))
!
!  ===================================================================
!  initialize rfunc_fftwr, the array to be transformed
!  ===================================================================
   rfunc_fftwr = ZERO
   do k = 1, isize(3)
      z = (k+istart(3)-2)*dz
      do j = 1, isize(2)
         y = (j+istart(2)-2)*dy
         do i = 1, isize(1)
            x = (i+istart(1)-2)*dx
!           idm = i + (j-1)*(nx+2) + (k-1)*ny*(nx+2)         ! beware of the padding space
!           rfunc_fftwr(idk) = exp(alpha*x+beta*y+gamma*z) ! used for in-place transformation
            p_rfunc_fftwr(i,j,k) = exp(alpha*x+beta*y+gamma*z)
         enddo
      enddo
   enddo
!
   t0 = getTime()
!  ===================================================================
!  transform from physical space to wavenumber space
!  -------------------------------------------------------------------
   call ftran_r2c(p_rfunc_fftwr,p_rfunc_fftwr,'fft')
!  -------------------------------------------------------------------
   if (MyPE == 0) then
      write(6,'('' time for in-place fftw real ='',1f10.5)')getTime()-t0
   endif
!
   rfunc_fftwr = rfunc_fftwr/vfac
!
!  ===================================================================
!  determine the k-vector for the MPI-FFTW r2c routine and check results
!  translate the in-place FFT results and place them into ft_fftw_ip
!  ===================================================================
   if (fend(1) < ifc) then
      fe = fend(1)
   else
      fe = fend(1) - 1
   endif
   if (fstart(1) == 1) then
      fs = 2
   else
      fs = fstart(1)
   endif
   numkx = fe-fs+1+fsize(1)
!
   nump_local = numkx*fsize(2)*fsize(3)
   write(6,'(/,a,8i8,/)')'MyPE, fs,fe,xstart,xend,xsize,numx,numlocal: ',&
                          MyPE,fs,fe,fstart(1),fend(1),fsize(1),numkx,nump_local
   write(6,'(/,a,3i8,/)')'MyPE, isize, numlocal: ',&
                          MyPE,isize(1)*isize(2)*isize(3),nump_local
!
   allocate( kvec_fftw(3,nump_local) , ft_fftw_ip(nump_local) )
   kvec_fftw = ZERO
   ft_fftw_ip = CZERO
!
   ka = PI2/a
   kb = PI2/b
   kc = PI2/c
   idk = 0
   LOOP_kg: do kg = fstart(3), fend(3)
      k = kg-fstart(3)+1
!     if (kg < kfc .or. nz == 1) then
!        kz = (kg-1)*kc
!     else
!        kz = (kg-1-nz)*kc  ! negative kz
!     endif
      LOOP_jg: do jg = fstart(2), fend(2)
         j = jg-fstart(2)+1
!        if (jg < jfc .or. ny == 1) then
!           ky = (jg-1)*kb
!        else
!           ky = (jg-1-ny)*kb  ! negative ky
!        endif
!
!        =============================================================
!        Loop over 1, 2, ..., Nx/2+1 values of kx
!        since kx's are distributed, each processor loops over the kx's 
!        mapped on it.
!        =============================================================
         LOOP_ig: do ig = fstart(1), fend(1) ! numkx+fstart(1)-1
            i = ig-fstart(1)+1
            if (ig < ifc .or. nx == 1) then
               kx = (ig-1)*ka
            else
               kx = (ig-1-nx)*ka  ! negative kx
            endif
!           ==========================================================
!           determine ky, kz that are mapped on the current processor
!           ==========================================================
            if (jg < jfc .or. ny == 1) then
               ky = (jg-1)*kb
            else
               ky = (jg-1-ny)*kb  ! negative ky
            endif
            if (kg < kfc .or. nz == 1) then
               kz = (kg-1)*kc
            else
               kz = (kg-1-nz)*kc  ! negative kz
            endif
!
            idk = idk + 1
            kvec_fftw(1,idk) = kx
            kvec_fftw(2,idk) = ky
            kvec_fftw(3,idk) = kz
!
!           ==========================================================
!           for the in-place case
!           ==========================================================
            idm = 2*i+(j-1)*2*fsize(1)+(k-1)*fsize(2)*2*fsize(1)
!           write(6,'(2i4,2x,3f10.5,2(2x,d15.8))')MyPE,idm,kvec_fftw(1:3,idk),rfunc_fftwr(idm-1),-rfunc_fftwr(idm)
            ft_fftw_ip(idk) = cmplx(rfunc_fftwr(idm-1), -rfunc_fftwr(idm), kind=CmplxKind)
         enddo LOOP_ig
!
!        =============================================================
!        Loop over Nx/2+2, Nx/2+3, ..., Nx values of kx, which are negative
!        since kx's are distributed, each processor loops over the kx's 
!        that those corresponding -kx (>0) are mapped on it.
!        =============================================================
         do ig = fe, fs, -1
            i = ig-fstart(1)+1
!
            kx = -(ig-1)*ka      ! negative kx
!           ==========================================================
!           determine ky, kz that are mapped on the current processor
!           ==========================================================
            if (jg <= jfc) then
               ky = -(jg-1)*kb   ! negative ky
            else
               ky = (ny-jg+1)*kb
            endif
            if (kg <= kfc) then
               kz = -(kg-1)*kc   ! negative kz
            else
               kz = (nz-kg+1)*kc
            endif
!
            idk = idk + 1
            kvec_fftw(1,idk) = kx
            kvec_fftw(2,idk) = ky
            kvec_fftw(3,idk) = kz
!
!           ==========================================================
!           for the in-place case
!           ==========================================================
            idm = 2*i+(j-1)*2*fsize(1)+(k-1)*fsize(2)*2*fsize(1)
!           write(6,'(2i4,2x,3f10.5,2(2x,d15.8),a)')MyPE,idm,kvec_fftw(1:3,idk),rfunc_fftwr(idm-1),rfunc_fftwr(idm),' *'
            ft_fftw_ip(idk) = cmplx(rfunc_fftwr(idm-1), rfunc_fftwr(idm), kind=CmplxKind)
         enddo
      enddo LOOP_jg
   enddo LOOP_kg
!
   if (MyPE == 0) then
      do idk = 1, nump_local
         write(6,'(a,3f12.5,4x,3f12.5)')'kvec_fft,kvec = ',kvec_fftw(1:3,idk),kvec(1:3,idk)
      enddo
   endif
!
   call syncAllPEs()
! 
   if (nt < 8*man_limit) then
      do idk = 1, nump_local
         found = .false.
         LOOP_idm_para: do idm = 1, nt
            if ( abs(kvec_fftw(1,idk)-kvec_nr(1,idm)) < TEN2m6 .and.          &
                 abs(kvec_fftw(2,idk)-kvec_nr(2,idm)) < TEN2m6 .and.          &
                 abs(kvec_fftw(3,idk)-kvec_nr(3,idm)) < TEN2m6 ) then
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
            write(6,'(a,2i4,2x,3f10.5)')'The parallel k-vector can not be matched:',MyPE,idk,kvec_fftw(1:3,idk)
         endif
      enddo
      do idm = 1, nt
         found = .false.
         LOOP_idk_serial: do idk = 1, nump_local
            if ( abs(kvec_fftw(1,idk)-kvec_nr(1,idm)) < TEN2m6 .and.          &
                 abs(kvec_fftw(2,idk)-kvec_nr(2,idm)) < TEN2m6 .and.          &
                 abs(kvec_fftw(3,idk)-kvec_nr(3,idm)) < TEN2m6 ) then
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
         endif
      enddo
   endif
!
   nullify(p_rfunc_fftwr)
   deallocate( kvec_fftw, ft_fftw_ip, rfunc_fftwr )
!
!  ===================================================================
!  Free work space
!  -------------------------------------------------------------------
   call p3dfft_clean()
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
contains
   include '../lib/arrayTools.F90'
end program ft_p3dfft
