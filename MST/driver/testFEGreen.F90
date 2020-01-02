program testFEGreen
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use SystemModule, only : getNumAtoms, getAtomPosition
!
   use ScfDataModule, only : ErBottom, ErTop, EiBottom, EiTop, Temperature
!
   use AtomModule, only : getPhiLmax, getRhoLmax, getStepFuncLmax
!
   use ContourModule, only : getNumEs, getEPoint, setupContour
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, SQRTm1, TEN2m8, PI4, PI
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
   use SphericalHarmonicsModule, only : calYlm, calYlmConjg
!
   use IntegerFactorsModule, only : lofk, mofk, jofk, lofj, mofj, kofj, m1m
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
   use BesselModule, only : SphericalBessel, SphericalHankel
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use RadialGridModule, only : getGrid
!
   use PolyhedraModule, only : getVolume, getInscrSphVolume
!
   use StepFunctionModule, only : getVolumeIntegration
!
   implicit   none
!
   character (len=4) :: istop = 'none'
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: lmax_phi_max
   integer (kind=IntKind) :: kmax_phi_max
   integer (kind=IntKind) :: lmax_max
   integer (kind=IntKind) :: ne
!
   integer (kind=IntKind) :: i, j, kl, ie, l, m, ir, i3, jl3, kl3
   integer (kind=IntKind) :: fstatus
   integer (kind=IntKind), allocatable :: lmax_phi(:), jmax_phi(:), kmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_green(:), jmax_green(:)
!
   integer (kind=IntKind), pointer :: nj3(:,:)
   integer (kind=IntKind), pointer :: kj3(:,:,:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
!
   real (kind=RealKind), parameter :: tol = TEN2m8
   real (kind=RealKind) :: r1, r2, d, dos, dos_mt
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
!
   complex (kind=CmplxKind) :: energy, kappa, x1, x2, green_expand, green, yyp, x
   complex (kind=CmplxKind), pointer :: e_mesh(:)
   complex (kind=CmplxKind), allocatable :: ylm1(:)
   complex (kind=CmplxKind), allocatable :: ylm2(:)
   complex (kind=CmplxKind), allocatable :: bjl(:)
   complex (kind=CmplxKind), allocatable :: bhl(:)
   complex (kind=CmplxKind), allocatable :: dos_r_jl(:,:)
!
   type (GridStruct), pointer :: Grid
!
!  -------------------------------------------------------------------
   call startProcess()
   call setupContour( ErBottom, ErTop, EiBottom, EiTop )
   ne = getNumEs()
   e_mesh => getEPoint()
   NumAtoms = getNumAtoms()
!  -------------------------------------------------------------------
!
   allocate(AtomPosition(1:3,1:NumAtoms), lmax_phi(NumAtoms), jmax_phi(NumAtoms), kmax_phi(NumAtoms))
   allocate(lmax_green(NumAtoms), jmax_green(NumAtoms))
   lmax_phi_max = 0; lmax_max = 0
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
      lmax_phi(i) = getPhiLmax(i)
      lmax_green(i) = 2*lmax_phi(i)
      jmax_green(i) = (lmax_green(i)+1)*(lmax_green(i)+2)/2
      jmax_phi(i) = (lmax_phi(i)+1)*(lmax_phi(i)+2)/2
      kmax_phi(i) = (lmax_phi(i)+1)**2
      lmax_phi_max = max(lmax_phi_max,lmax_phi(i))
      lmax_max = max(lmax_max,2*lmax_phi_max,getStepFuncLmax(i))
   enddo
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(lmax_max,istop,iprint)
!  -------------------------------------------------------------------
   call setupRadGridAndCell(NumAtoms,lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize kmax:
!  ===================================================================
   kmax_phi_max=(lmax_phi_max+1)**2
   allocate(ylm1(kmax_phi_max),ylm2(kmax_phi_max),bjl(0:lmax_phi_max),bhl(0:lmax_phi_max))
!
   write(6,'(/,a)')'Free-electron Green function test....'
   write(6,'(a,2i5,/)')'lmax, kmax = ',lmax_phi,kmax_phi
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   write(6,'(a)')'  Re[Energy]  Im[Energy]     DOS_in_Cell      DOS_in_MT      Vcell*sqrt(|E|)/PI2**2    Vmt*sqrt(|E|)/PI2**2'
   do i=1, NumAtoms
      Grid => getGrid(i)
      allocate(dos_r_jl(Grid%jend,jmax_green(i)))
      do ie = 1, ne
!        =============================================================
!        set energy &  electron momentum: ............................
!        =============================================================
         energy=e_mesh(ie)
         kappa=sqrt(energy)
         dos_r_jl = CZERO
         do ir = 1, Grid%jend
            x = kappa*Grid%r_mesh(ir)
!           ----------------------------------------------------------
            call SphericalBessel(lmax_phi(i),x,bjl)
            call SphericalHankel(lmax_phi(i),x,bhl)
!           ----------------------------------------------------------
            do kl = 1, kmax_phi(i)
               l = lofk(kl)
               do i3 = 1, nj3(kl,kl)
                  kl3 = kj3(i3,kl,kl)
                  if (mofk(kl3) >= 0) then
                     jl3 = jofk(kl3)
                     dos_r_jl(ir,jl3) = dos_r_jl(ir,jl3) + bjl(l)*bhl(l)*cgnt(i3,kl,kl)
                  endif
               enddo
            enddo
         enddo
         dos_r_jl = kappa*dos_r_jl/PI   ! DOS = real(sum{dos_r_jl*Ylm})
!        -------------------------------------------------------------
         dos = getVolumeIntegration( i, Grid%jend, Grid%r_mesh,       &
                                     jmax_green(i), 0, dos_r_jl, dos_mt)
!        -------------------------------------------------------------
         write(6,'(2f12.8,2x,d15.8,2x,d15.8,6x,d15.8,10x,d15.8)')energy,dos,dos_mt, &
                   getVolume(i)*sqrt(abs(e_mesh(ie)))/(4.0d0*PI**2),    &
                   getInscrSphVolume(i)*sqrt(abs(e_mesh(ie)))/(4.0d0*PI**2)
      enddo
      deallocate(dos_r_jl)
   enddo
!
   do i=1,NumAtoms-1
!     ----------------------------------------------------------------
      call calYlm(AtomPosition(1:3,i),lmax_phi_max,ylm1(1:kmax_phi_max))
!     ----------------------------------------------------------------
      r1 =  sqrt(AtomPosition(1,i)**2+AtomPosition(2,i)**2+AtomPosition(3,i)**2)
      do j = i+1, NumAtoms
!        -------------------------------------------------------------
         call calYlmConjg(AtomPosition(1:3,j),lmax_phi_max,ylm2(1:kmax_phi_max))
!        -------------------------------------------------------------
         r2 =  sqrt(AtomPosition(1,j)**2+AtomPosition(2,j)**2+AtomPosition(3,j)**2)
         d = sqrt((AtomPosition(1,i)-AtomPosition(1,j))**2+           &
                  (AtomPosition(2,i)-AtomPosition(2,j))**2+           &
                  (AtomPosition(3,i)-AtomPosition(3,j))**2)
         write(6,'(/,'' i, j, d = '',2i5,2x,d15.8)')i,j,d
         do ie = 1, ne
!           ==========================================================
!           set energy &  electron momentum: ............................
!           ==========================================================
            energy=e_mesh(ie)
            kappa=sqrt(energy)
            green = -exp(SQRTm1*kappa*d)/(PI4*d)
!           ==========================================================
            if (r1 <= r2) then
               x1 = r1*kappa
               x2 = r2*kappa
            else
               x1 = r2*kappa
               x2 = r1*kappa
            endif
!           ----------------------------------------------------------
            call SphericalBessel(lmax_phi_max,x1,bjl)
            call SphericalHankel(lmax_phi_max,x2,bhl)
!           ----------------------------------------------------------
            green_expand = CZERO
            kl = kmax_phi_max+1
            do l = lmax_phi_max, 0, -1
               yyp = CZERO
               do m = l, -l, -1
                  kl = kl - 1
                  yyp = yyp + ylm1(kl)*ylm2(kl)
               enddo
               green_expand = green_expand + bjl(l)*bhl(l)*yyp
            enddo
            green_expand = -SQRTm1*kappa*green_expand
            write(6,'(2d15.8,2x,2d15.8,2x,2d15.8,6x,d15.8)')energy,green,green_expand, &
                                                            abs(green-green_expand)
         enddo
      enddo
   enddo
!
   deallocate(ylm1,ylm2,AtomPosition,bjl,bhl)
!
!  -------------------------------------------------------------------
   call endGauntFactors()
   call endSphericalHarmonics()
   call finishProcess()
!  -------------------------------------------------------------------
!
end program testFEGreen
