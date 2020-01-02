program  madpot
!  ********************************************************************
!  calculate the Madelung potential created by point charges with uniform
!  background
!
!  Here we assume that the electroic multipole moments are ZERO
!  The electrostatic potential is essentially V_{inter} + V_{Mad}
!  ********************************************************************
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use TimerModule, only : initTimer, getTime
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, TWO, CZERO, SQRTm1, PI4
   use MathParamModule, only : PI2, THIRD, SIX
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : calYlm
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use DataServiceCenterModule, only : getDataStorage,                &
                                       RealType, RealMark,            &
                                       IntegerType, IntegerMark,      &
                                       CharacterType, CharacterMark
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use MadelungModule, only : initMadelung, endMadelung
   use MadelungModule, only : getMadelungMatrix
   use MadelungModule, only : getDLMatrix
   use MadelungModule, only : getDLFactor
   use MadelungModule, only : getJelliumPotential
   use MadelungModule, only : getRSpaceLattice
   use MadelungModule, only : getNumRSVectors
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra, &
                               getVolume, genPolyhedron, &
                               getInscrSphRadius, getOutscrSphRadius, &
                               getWignerSeitzRadius, printPolyhedron
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid, & 
                                getGrid, getMaxNumRmesh
   use RadialGridModule, only : genRadialGrid, printRadialGrid
!
   use StepFunctionModule, only : initStepFunction,     &
                                  endStepFunction,      &
                                  getVolumeIntegration
!
   use InterpolationModule, only : FitInterp
!
   use IntegrationModule, only : calIntegration
!
   use MPPModule, only : initMPP, endMPP
!
   implicit   none
!
   logical :: testMadGen
!
   character (len=60) :: fname
!
   integer (kind=IntKind) :: num_atoms
   integer (kind=IntKind) :: numlocal_atoms
   integer (kind=IntKind) :: lmax_rho
   integer (kind=IntKind) :: lmax_pot
   integer (kind=IntKind) :: kmax_pot
   integer (kind=IntKind) :: jmax_pot
   integer (kind=IntKind) :: lmax_mad
   integer (kind=IntKind) :: jmax_mad
   integer (kind=IntKind) :: funit
   integer (kind=IntKind) :: jend, jmt, nr_max
!
   integer (kind=IntKind), parameter :: ndivin = 1001
   integer (kind=IntKind), parameter :: ndivout = 0
   integer (kind=IntKind), parameter :: nmult = 1
!
   integer (kind=IntKind) :: i, j, k, jl, kl, n, l, m, status, nr, ir
   integer (kind=IntKind) :: jlp, klp, lp, mp, jlt, nrs, m1m
   integer (kind=IntKind), allocatable :: gindex(:), lmax_step(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind) :: l_pot, m_pot, jl_pot, kl_pot, kl_rho, l_rho, m_rho, jl_rho, mdif, lsum
!
   real (kind=RealKind) :: t0
   real (kind=RealKind) :: a0, dr, rt, rdist, v, r, pot, pot2, rfac, rb, dps
   real (kind=RealKind) :: r2l, signm, Y0inv, vec(3), vecb(3)
   real (kind=RealKind) :: ztot, rho0, vol, vm, a0fac, Z0, vtmp, V0_mad1, V0_mad2, v0_mt
   real (kind=RealKind), allocatable :: lfac(:), vmad0_r(:), rho_0(:), rho_0r(:), rr(:), sqrt_r(:), V1_r(:)
   real (kind=RealKind), allocatable :: xpos(:), ypos(:), zpos(:), zj(:)
   real (kind=RealKind), pointer :: Bravais(:,:)
   real (kind=RealKind), pointer :: atom_position(:,:)
   real (kind=RealKind), pointer :: dlf(:)
   real (kind=RealKind), pointer :: rslat(:,:)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), pointer :: madmat(:)
!
   complex (kind=CmplxKind) :: vl, sumat, sumjl
   complex (kind=CmplxKind), pointer :: dlm(:,:)
   complex (kind=CmplxKind), allocatable :: Ylm(:)
   complex (kind=CmplxKind), allocatable :: Ylmb(:)
   complex (kind=CmplxKind), allocatable :: vmad0(:,:), v_tilt(:,:)
!
   type (GridStruct), pointer :: Grid
!
   interface
      subroutine printDensity_r(denType, id, nrs, r_mesh, den_r)
      use KindParamModule, only : IntKind, RealKind
      implicit none
      character (len=*), intent(in) :: denType
      integer (kind=IntKind), intent(in) :: id, nrs
      real (kind=RealKind), intent(in), target :: r_mesh(1:nrs)
      real (kind=RealKind), intent(in), target :: den_r(1:nrs)
      end subroutine printDensity_r
   end interface
!
   interface
      subroutine readPositionData(fname,NumAtomsIn,NumAtomsOut)
         use KindParamModule, only : IntKind
         character (len=*), intent(in) :: fname
         integer (kind=IntKind), intent(in), optional :: NumAtomsIn
         integer (kind=IntKind), intent(out), optional :: NumAtomsOut
      end subroutine readPositionData
   end interface
!
!  -------------------------------------------------------------------
   call initMPP()
!  -------------------------------------------------------------------
!
   num_atoms = 0
   do while (num_atoms < 1)
      write(6,'(1x,a,$)')'Number of atoms: '
      read(5,*)num_atoms
   enddo
   write(6,'(i4)')num_atoms
!
!  lmax_rho = -1
!  do while (lmax_rho < 0)
!     write(6,'(/,1x,a,$)')'Lmax for electron density: '
!     read(5,*)lmax_rho
!  enddo
!  write(6,'(i4)')lmax_rho
   lmax_rho = 0
!
   lmax_pot = -1
   do while (lmax_pot < 0)
      write(6,'(/,1x,a,$)')'Lmax for potential expansion: '
      read(5,*)lmax_pot
   enddo
   write(6,'(i4)')lmax_pot
!
   status = -1
   do while (status < 0)
      write(6,'(/,1x,a,$)')'Position data file name: '
      read(5,'(a)',iostat=status)fname
   enddo
   write(6,'(a)')trim(fname)
!
!  -------------------------------------------------------------------
   call readPositionData(fname,NumAtomsIn=num_atoms)
!  -------------------------------------------------------------------
!
   write(6,'(/,1x,a)')'This Code Will Calculate Electrostatic Potential'
!
   nr = -1
   do while (nr < 1)
      write(6,'(/,1x,a,$)')'On how many mesh points: '
      read(5,*)nr
   enddo
!
   allocate(xpos(nr), ypos(nr), zpos(nr))
!
   if (nr == 1) then
      status = -1
      do while (status < 0)
         write(6,'(/,1x,a,$)')'On which point (x, y, z): '
         read(5,*,iostat=status)xpos(1), ypos(1), zpos(1)
      enddo
   else
      i = 0; j = 0; k = 0
      do while (i == 0 .and. j == 0 .and. k == 0)
         write(6,'(/,1x,a)')'Along which direction to calculate potential: '
         write(6,'(1x,a,$)')'Type in three integers (e.g., 1, 0, 0): '
         read(5,*,iostat=status)i, j, k
      enddo
!
      rdist = -ONE
      do while (rdist <= TEN2m6)
         write(6,'(/,1x,a,$)')'Up to distance: '
         read(5,*)rdist
      enddo
!
      rt = i*i+j*j+k*k
      dr = rdist/sqrt(rt)/real(nr-1,RealKind)
      xpos(1) = i*TEN2m6
      ypos(1) = j*TEN2m6
      zpos(1) = k*TEN2m6
      do n = 2, nr
         xpos(n) = i*(n-1)*dr
         ypos(n) = j*(n-1)*dr
         zpos(n) = k*(n-1)*dr
      enddo
      open(unit=12,file='madpot.dat',status='unknown',form='formatted')
      write(12,'(a,3(i2,a))')'# Direction: [',i,',',j,',',k,']'
      write(12,'(a,i5)')'# Number of Points: ',nr
      write(12,'(a,i5)')'# Lmax for Potential: ',lmax_pot
   endif
!
   Bravais => getDataStorage('Bravais Vector',3,3,RealMark)
   atom_position => getDataStorage('Atomic Position',3,num_atoms,RealMark)
!
   write(6,'(/,a)')'Bravais Lattice:'
   write(6,'(3f15.8)')(Bravais(1:3,j),j=1,3)
!
   write(6,'(/,a)')'Atomic Position:'
   do i=1,num_atoms
      write(6,'(3f15.8)')atom_position(1:3,i)
   enddo
   write(6,'(/)')
!
!  ===================================================================
!  initilize kmax:
!  ===================================================================
   jmax_pot = (lmax_pot+1)*(lmax_pot+2)/2
   lmax_mad = max(lmax_rho,lmax_pot)
   jmax_mad = (lmax_mad+1)*(lmax_mad+2)/2
   kmax_pot = (lmax_pot+1)*(lmax_pot+1)
!
   allocate(Ylm(1:kmax_pot), lfac(0:lmax_mad))
   allocate(Ylmb(1:kmax_pot))
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(4*lmax_mad)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(2*lmax_mad,'none',-1)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize AtomInfo parameters
!  ===================================================================
!
   allocate( gindex(1:num_atoms), zj(1:num_atoms) )
   allocate( ngr(1:num_atoms), ngt(1:num_atoms), lmax_step(1:num_atoms) )
!
   do i=1,num_atoms
      ngr(i) = 80
      ngt(i) = 60
      lmax_step(i) = lmax_pot
   enddo
!
   do i=1,num_atoms
      gindex(i) = i
      zj(i) = ONE     ! force each atom having unit charge
   enddo
!
   vol = (Bravais(2,1)*Bravais(3,2)-Bravais(3,1)*Bravais(2,2))*Bravais(1,3)+ &
         (Bravais(3,1)*Bravais(1,2)-Bravais(1,1)*Bravais(3,2))*Bravais(2,3)+ &
         (Bravais(1,1)*Bravais(2,2)-Bravais(2,1)*Bravais(1,2))*Bravais(3,3)
   vol = abs(vol)
!
   ztot = ZERO
   do i=1,num_atoms
      ztot = ztot + zj(i)
   enddo
!
   lfac(0) = ONE
   do l = 1, lmax_mad
      lfac(l) = (2*l+1)*lfac(l-1)
   enddo
!
   testMadGen = .false.
!   testMadGen = .true.
!  -------------------------------------------------------------------
   call initRadialGrid(num_atoms, 'none', 0)
   if ( .not.testMadGen ) then
!  --------------------------------------------------------------------
   call initPolyhedra(num_atoms,Bravais,'none',0)
!  --------------------------------------------------------------------
   do i=1,num_atoms
!     ----------------------------------------------------------------
      call genPolyhedron(i,i,num_atoms,atom_position)
!     ----------------------------------------------------------------
!      if ( mod(i-1,10000)==0 ) then
         call printPolyhedron(i)
!      endif
!     ----------------------------------------------------------------
      call genRadialGrid(i,getInscrSphRadius(i), getInscrSphRadius(i),&
                         getWignerSeitzRadius(i), getOutscrSphRadius(i), ndivin,ndivout,nmult)
!     ----------------------------------------------------------------
   enddo
!  --------------------------------------------------------------------
   call initStepFunction(num_atoms,maxval(lmax_step),lmax_step,ngr,ngt,'none',0)
!  --------------------------------------------------------------------
   endif
   nr_max = getMaxNumRmesh()
   call initTimer()
   t0 = getTime()
!   numlocal_atoms = 2
   call initMadelung(num_atoms,num_atoms,gindex,lmax_rho,lmax_pot,     &
                     Bravais,atom_position,0)
 !  write(6,*) "End init Madelung"
   allocate( vmad0(nr_max,jmax_pot), v_tilt(nr_max,jmax_pot), vmad0_r(0:nr_max) )
   allocate( V1_r(0:nr_max), sqrt_r(0:nr_max) )
   allocate( rho_0(0:nr_max), rho_0r(0:nr_max), rr(0:nr_max) )
!
   Y0inv = sqrt(PI4)
   rho0 = ztot/vol
   Z0 = ztot/num_atoms
   rho_0(0:nr_max) = rho0*Y0inv
   write(6,'(a,f16.8)')'rho0 = ',rho0
   write(6,'(a,f16.8)')'Z0   = ',Z0
   write(6,'(a,f16.8)')'vol   = ', vol
   write(6,'(a,f16.8)')'vol   = ', getVolume(1)
!
   write(6,'('' time for set up MadelungMatrix ='',1f10.5)')getTime()-t0
!  --------------------------------------------------------------------
   rslat => getRSpaceLattice(a0fac)
   nrs = getNumRSVectors()
!
   do n=1,num_atoms
!      if ( mod(i-1,10000)==0 ) then
         write(6,'(a,i5,a,3f15.8)')'Atom Index:',n,',  Position:',        &
                                   atom_position(1:3,n)
!      endif
!     ----------------------------------------------------------------
      madmat => getMadelungMatrix(n)
      dlm => getDLMatrix(n,a0)
!     ----------------------------------------------------------------
!      if ( mod(i-1,10000)==0 ) then
         write(6,'(a,f12.8)')'a0 = ',a0
         write(6,'(a)')'Madelung Matrix:'
         write(6,'(i5,2x,d15.8)')(i,madmat(i), i=1,num_atoms)
!      endif
!
      vm = ZERO
      do j = 1, num_atoms
         vm = vm + madmat(j)
      enddo
      V0_mad1 = ZERO
      V0_mad2 = ZERO
!
!      if ( mod(i-1,10000)==0 ) then
         write(6,'(a)')'DL Matrix:'
!      endif
   jl = 0
   do l = 0, lmax_pot
      do m = 0, l
         jl = (l+1)*(l+2)/2 - l + m
         kl = (l+1)*(l+1)- l + m
         do lp = 0, lmax_rho
            do mp = 0, lp
               jlp = (lp+1)*(lp+2)/2 - lp + mp
               klp = (lp+1)*(lp+1)-lp+mp
               jlt = (l+lp+1)*(l+lp+1) - (l+lp) + mp - m
               do i = 1, num_atoms
                  if (abs(dlm(i,kl)/a0**l) > TEN2m6 ) then 
!                     if ( mod(i-1,10000)==0 ) then
                     write(6,'(5i5,2x,2d15.8)')i,l,m,lp,mp,dlm(i,jlt)/a0**(l+lp)
!                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
!
!   if ( mod(i-1,10000)==0 ) then
      write(6,'(//,a)')'M_{LL} Matrix:'
!   endif
   jl = 0
   do l = 0, lmax_pot
      do m = 0, l
         jl = (l+1)*(l+2)/2 - l + m
         kl = (l+1)*(l+1)-l + m
         dlf => getDLFactor(jl)
         jlp = 0
         do lp = 0, lmax_rho
            do mp = 0, lp
               jlp = jlp + 1
               klp = (lp+1)*(lp+1)-lp+mp
               jlt = (l+lp+1)*(l+lp+1) - (l+lp) + mp - m
               do i = 1, num_atoms
                  if (abs(dlf(klp)*dlm(i,jlt)/a0**(l+lp)) > TEN2m6 ) then
!                     if ( mod(i-1,10000)==0 ) then
                     write(6,'(5i5,2x,2d15.8)')i,l,m,lp,mp, &
                           dlf(klp)*dlm(i,jlt)/a0**(l+lp)
!                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
      Grid => getGrid(n)
      r_mesh => Grid%r_mesh
      rr(1:nr_max) = r_mesh
      rr(0) = ZERO
      jend = Grid%jend
      jmt = Grid%jmt
!
      sqrt_r(0) = ZERO
      do ir = 1, jend
         sqrt_r(ir) = sqrt(r_mesh(ir))
      enddo
!
      do ir = 1,jend
         rho_0r(ir) = rho0*Y0inv/r_mesh(ir)
      enddo
      vtmp =  getVolumeIntegration(n, jend, r_mesh(1:jend), -1,   &
                                 rho_0(1:jend), truncated=.false.)
      do ir = 1,jmt
         vmad0(ir,1) = TWO*( PI4*THIRD*rho0*Y0inv*(r_mesh(ir)**2)   +   &
                         ( vtmp - PI2*rho0*Y0inv*(r_mesh(ir)**2) )  - Z0*Y0inv/r_mesh(ir) )
         vmad0_r(ir) = real(vmad0(ir,1))*r_mesh(ir)/Y0inv
      enddo
      do ir = jmt+1,jend
         vmad0(ir,1) = TWO*( getVolumeIntegration(n, ir, rr(1:ir), 0,   &
                                 rho_0(1:ir), truncated=.false.)/r_mesh(ir) + &  !+   &
                         ( vtmp  &
                          - getVolumeIntegration(n, ir, rr(1:ir), -1, &
                                 rho_0(1:ir), truncated=.false.) ) - Z0*Y0inv/r_mesh(ir) )
         vmad0_r(ir) = real(vmad0(ir,1))*r_mesh(ir)/Y0inv
      enddo
!      if ( mod(i-1,10000)==0 ) then
         call printDensity_r("Phi1",n,jend,r_mesh(1:jend),vmad0_r(1:jend))
!      endif
!     =========================
!     Using sphere integration
!     =========================
      call FitInterp( 4, sqrt_r(1:4), vmad0_r(1:4), ZERO,             &
                     vmad0_r(0), dps )
!     ----------------------------------------------------------------
      call calIntegration( jend+1, sqrt_r(0:jend), vmad0_r(0:jend),   &
                           V1_r(0:jend), 3)
      V0_mad1 = SIX*V1_r(jend)/(r_mesh(jend)**3)
!      if ( mod(i-1,10000)==0 ) then
         write(6,*) " Potential Zero :: ", V0_mad1/getVolume(n)
!      endif
!     ===================================
!     Using Voronoi Polyhedra integration
!     ===================================
!      V0_mad1 = V0_mad1 + getVolumeIntegration( n, jend,&
!                                        r_mesh, 1, 1,   &
!                                        0, vmad0, truncated=.false. )
!      write(6,*) " Potential Zero :: ", V0_mad1/getVolume(n)
!
      m1m = 1
      LoopJl: do l_pot = 0,lmax_pot
         do m_pot = 0,l_pot
            jl_pot = (l_pot+1)*(l_pot+2)/2 - l_pot + m_pot
            kl_pot = (l_pot+1)*(l_pot+1) - l_pot + m_pot
!
            dlf => getDLFactor(jl_pot)
!
            sumjl = czero
            do kl_rho = 1, 1, -1
               m_rho = 0
               l_rho = 0
               jl_rho = 1
!
               mdif = m_rho-m_pot
!            mdif = m_pot-m_rho
               lsum = l_pot+l_rho
               kl = (lsum+1)*(lsum+1)-lsum+mdif
!
               sumat = czero
               do i = 1,num_atoms
                  sumat = sumat - dlm(i,kl)*Z0*Y0inv
               enddo
               sumjl = sumjl + dlf(kl_rho)*sumat/a0**lsum
            enddo
!
            if ( l_pot==0 ) then
               do i = 1,jend
                  v_tilt(i,jl_pot) = v_tilt(i,jl_pot)  +                            &
                        cmplx( real(TWO*sumjl - PI4*THIRD*rho0*Y0inv*r_mesh(i)**2 - &
                                    TWO*Z0*Y0inv/r_mesh(i), kind=RealKind),      &
                                    ZERO,kind=CmplxKind)
               enddo
            else if ( m_pot == 0 ) then
               do i = 1,jend
                  v_tilt(i,jl_pot) = v_tilt(i,jl_pot) +                  &
                        cmplx( real(TWO*m1m*sumjl*(r_mesh(i)**l_pot),     &
                                    kind=RealKind), ZERO, kind=CmplxKind )
               enddo
            else
               do i = 1,jend
                   v_tilt(i,jl_pot) = v_tilt(i,jl_pot) +                &
                                      TWO*m1m*sumjl*(r_mesh(i)**l_pot)
               enddo
            endif
            jl_pot = jl_pot + 1
            kl_pot = kl_pot + 1
         enddo
         m1m = (-1)*m1m
!
      enddo LoopJl
!
      do i = 1,jend
         vmad0_r(i) = real(v_tilt(i,1),kind=RealKind)/Y0inv
      enddo
!      if ( mod(i-1,10000)==0 ) then
         call printDensity_r("Phi2",n,jend,r_mesh,vmad0_r(1:jend))
!      endif
      V0_mad2 = getVolumeIntegration( n, jend,&
                                     r_mesh, kmax_pot, jmax_pot,   &
                                     0, v_tilt, v0_mt, truncated=.false. )
!      if (nr > 1 .and. mod(i-1,10000)==0 )  then
         write(6,*) " R_MT        ::",  r_mesh(jmt)
         write(6,*) " R_WS        ::",  r_mesh(jend)
         write(6,*) " Phi1        :: ", V0_mad1/getVolume(n)
         write(6,*) " Phi1/Z0     :: ", Y0inv*V0_mad1/(rho0*getVolume(n))
         write(6,*) " Phi2        :: ", V0_mad2/getVolume(n)
         write(6,*) " Phi2/Z0     :: ", Y0inv*V0_mad2/(rho0*getVolume(n))
!      endif
!
      if ( nr > 1 ) then
         if ( mod(i-1,10000)==0 )  then
         write(12,'(a)')'#'
         write(12,'(a,i5)')'# Atom Index: ',n
         write(12,'(a)')'# ============================================='
         write(12,'(a)')'#      r      v_Madelung     v_Slater  DeltaV  '
         write(12,'(a)')'# ---------------------------------------------'
         funit = 12
         endif
      else
         funit = 6
      endif
      do i = 1, nr  ! loop over the mesh points along a line
         vec(1) = xpos(i); vec(2) = ypos(i); vec(3) = zpos(i)
!        -------------------------------------------------------------
         call calYlm(vec,lmax_pot,Ylm)
!        -------------------------------------------------------------
         r = sqrt(xpos(i)*xpos(i)+ypos(i)*ypos(i)+zpos(i)*zpos(i))
         v = ZERO
!         m1m = 1
!         r2l = (ONE-TWO*mod(lmax_pot,2))*(r/a0)**lmax_pot
         do l = 0,lmax_pot
            r2l = (-r/a0)**l
            do m = 0,l
               kl = (l+1)*(l+1) - l - m
               jl = (l+1)*(l+2)/2 - l + m
               dlf => getDLFactor(jl)
               vl = CZERO
!              =======================================================
!              sum over the number of atoms
!              =======================================================
               do j = 1, num_atoms
                  vl = vl + ((-ONE)**m)*dlm(j,kl)*zj(j)
               enddo
!              =======================================================
!              sum over {l,m}: V_{l,m} * Y_{l,m}
!
!              we only sum over non-negative m. So for m > 0, we need
!              to include a factor 2.0.
!              =======================================================
               if (m == 0) then
                  v = v - Y0inv*TWO*r2l*dlf(1)*real(vl*Ylm(kl),RealKind)
               else
                  v = v - Y0inv*TWO*r2l*dlf(1)*real(vl*Ylm(kl),RealKind)
               endif
            enddo
!            m1m = -m1m
         enddo
         v = v - TWO*zj(n)/r - TWO*PI2*THIRD*r*r*rho0
         pot = -zj(n)*getJelliumPotential(n,vec)
!         if ( mod(i-1,10000)==0 ) then
            write(funit,'(f12.5,3(6x,f21.8))')r,v/TWO,pot, pot-v/TWO
!         endif
!        write(funit,'(f12.5,6x,f21.8,6x,f21.8,6x,f21.8)')r,v,v/TWO,pot
!
!        =============================================================
!        use a different formula
!        =============================================================
         pot2 = ZERO
         do j = 1, num_atoms
            do ir = nrs, 1, -1
               vecb(1) = a0fac*rslat(ir,1) + atom_position(1,j)
               vecb(2) = a0fac*rslat(ir,2) + atom_position(2,j)
               vecb(3) = a0fac*rslat(ir,3) + atom_position(3,j)
               rb = sqrt(vecb(1)*vecb(1)+vecb(2)*vecb(2)+vecb(3)*vecb(3))
               if (rb > TEN2m6) then
                  call calYlm(vecb,lmax_pot,Ylmb)
                  do l = lmax_pot, 1, -1
                     kl = (l+1)*(l+1)
                     jl = (l+1)*(l+2)/2
                     rfac = PI4*zj(j)*(r/rb)**l/rb/(2*l+ONE)
                     do m = l, 0, -1
                        if (m == 0) then
                           pot2 = pot2 -                              &
                                  rfac*real(Ylm(kl)*conjg(Ylmb(kl)),RealKind)
                        else
                           pot2 = pot2 -                              &
                              TWO*rfac*real(Ylm(kl)*conjg(Ylmb(kl)),RealKind)
                        endif
                        kl = kl - 1
                        jl = jl - 1
                     enddo
                  enddo
               endif
            enddo
         enddo
         pot2 = pot2 - zj(n)/r - PI2*THIRD*r*r*rho0 - vm
!         if ( mod(i-1,10000)==0 ) then
            write(6,'(f12.5,6x,f21.8,6x,f21.8)')r,pot2,pot
!         endif
      enddo
   enddo
!
   if (nr > 1) then
      close(unit=12)
   endif
!
!  ====================================================================
!  check the accuracy of multipole expansion formula
!  ====================================================================
!
!  vecb(1:3) = Bravais(1:3,2) + Bravais(1:3,3) - Bravais(1:3,1)
   vecb(1:3) = Bravais(1:3,1)
   rb = sqrt(vecb(1)*vecb(1)+vecb(2)*vecb(2)+vecb(3)*vecb(3))
   call calYlm(vecb,lmax_pot,Ylmb)
!
   write(6,'(/,a)')'Expansion accuracy'
   write(6,'(a)')'=================================================='
   do i = 1, nr  ! loop over the mesh points along a line
      vec(1) = xpos(i); vec(2) = ypos(i); vec(3) = zpos(i)
!     ----------------------------------------------------------------
      call calYlm(vec,lmax_pot,Ylm)
!     ----------------------------------------------------------------
      r = sqrt(xpos(i)*xpos(i)+ypos(i)*ypos(i)+zpos(i)*zpos(i))
      v = ONE/sqrt((vec(1)-vecb(1))**2+(vec(2)-vecb(2))**2+(vec(3)-vecb(3))**2)
      pot = ZERO
      do l = lmax_pot, 0, -1
         r2l = (r/rb)**l/rb
         rfac = r2l/(TWO*l+ONE)
!        write(6,'(a,i5,f20.8)')'l, rfac = ',l,rfac
         do m = l, 0, -1
            kl = (l+1)*(l+1)-l+m
            jl = (l+1)*(l+2)/2-l+m
            if (m == 0) then
               pot = pot + rfac*real(Ylm(kl)*conjg(Ylmb(kl)),RealKind)
            else
               pot = pot + TWO*rfac*real(Ylm(kl)*conjg(Ylmb(kl)),RealKind)
            endif
         enddo
      enddo
      pot = pot*PI4
      write(6,'(a,f10.5,2f20.8)')'r,v,pot = ',r,v,pot
   enddo
!
   nullify(Bravais,dlf)
   nullify(dlm,madmat)
   nullify(atom_position)
!
   deallocate(sqrt_r, V1_r )
   deallocate(rho_0, rho_0r )
   deallocate( vmad0, v_tilt, vmad0_r )
   deallocate(gindex)
   deallocate(Ylm, Ylmb, lfac)
   deallocate(xpos, ypos, zpos, zj)
!
!  --------------------------------------------------------------------
   call endPolyhedra()
!  --------------------------------------------------------------------
   call endStepFunction()
!  --------------------------------------------------------------------
   call endMadelung()
!  --------------------------------------------------------------------
   call endGauntFactors()
   call endSphericalHarmonics()
   call endMPP()
!  --------------------------------------------------------------------
!
end program madpot
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDensity_r(denType, id, nrs, r_mesh, den_r)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only: Ten2m8
!
   implicit none
!
   character (len=*), intent(in) :: denType
!
   integer (kind=IntKind), intent(in) :: id, nrs
!
   real (kind=RealKind), intent(in), target :: r_mesh(1:nrs)
!
   real (kind=RealKind), intent(in), target :: den_r(1:nrs)
!
   character(len=20) :: file_den
   integer (kind=IntKind) :: lenDenName, MyPE
   integer (kind=IntKind) :: ir, funit
   integer (kind=IntKind) :: offset = 100000
!
   MyPE = 0
   if (MyPE /= 0) then
     return
   endif
   lenDenName = len(denType)+1
   file_den(1:20) = "                    "
   file_den(1:lenDenName-1) = trim(denType)
   file_den(lenDenName:lenDenName) = "_"
!
   if ( id == 1 ) then
      write(file_den(lenDenName+1:lenDenName+6),'(i6)') offset+MyPE+id
      file_den(lenDenName+1:lenDenName+1) = 'n'
      funit = 55+MyPE+id
      open(unit=funit,file=trim(file_den),status='unknown')
      write(funit,'(a)') "#Ind     r_mesh   (lm):"
      do ir = 1, nrs
         write(funit,'(i4,2(1x,d16.8))') id, r_mesh(ir), den_r(ir)
      enddo
      write(funit,'(a)') " "
      close(funit)
   endif
!
   end subroutine printDensity_r
