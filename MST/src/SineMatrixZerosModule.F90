module SineMatrixZerosModule
!  *******************************************************************
!  * Purpose: determine the zeros of single scattering sine matrix,  *      
!  *          and the residual matrix of sine matrix inverse         *
!  *******************************************************************
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                               Ten2m8, FOURTH, CZERO, TEN2m7
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : syncAllPEs, MyPE
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEInGroup
   use GroupCommModule, only : GlobalSumInGroup, bcastMessageInGroup
!
public :: initSineMatrixZeros,     &
          endSineMatrixZeros,      &
          findSineMatrixZeros,     &
          computeSineZeroDensity,  &
          getNumSineMatrixZeros,   &
          getSineMatrixZero,       &
          getSineZeroDensity,      &
          isSineZeroInEnergyRange, &
          printSineMatrixZerosInfo,&
          isSineMatrixZerosInitialized
!
private
   logical :: isInitialized = .false.
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind) :: eGID, NumPEsInEGroup, MyPEInEGroup
   integer (kind=IntKind) :: lmax_kkr_max, kmax_kkr_max, kmax_kkr_save
   integer (kind=IntKind) :: lmax_rho_max, kmax_rho_max, jmax_rho_max
   integer (kind=IntKind) :: MaxNumRs
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
   real (kind=RealKind), allocatable :: gaunt(:,:,:)
!
   real (kind=RealKind), parameter :: degen_tol = TEN2m6
!
   complex (kind=CmplxKind), allocatable, target :: wspace0(:), wspace1(:)
   complex (kind=CmplxKind), allocatable, target :: wspace2(:), wspace3(:)
   complex (kind=CmplxKind), allocatable, target :: wspace4(:)
!
   type ZeroDensityStruct
      integer (kind=IntKind) :: NumDegens
!
      real (kind=RealKind) :: ZeroE
!
      complex (kind=CmplxKind), allocatable :: ResidualMat(:)
      complex (kind=CmplxKind), allocatable :: Density(:)
      complex (kind=CmplxKind), allocatable :: Deriv_Density(:)
   end type ZeroDensityStruct
!
   type ZeroStruct
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind) :: jmax_rho
      integer (kind=IntKind) :: kmax_kkr
      integer (kind=IntKind) :: NumRs
      integer (kind=IntKind), allocatable :: NumZeros(:)
      integer (kind=IntKind), allocatable :: ZeroIndex(:,:)
!
      real (kind=RealKind), allocatable :: ebot(:)
      real (kind=RealKind), allocatable :: etop(:)
!
      type (ZeroDensityStruct), allocatable :: ZeroState(:,:)
   end type ZeroStruct
!
   type (ZeroStruct), allocatable :: SineZero(:,:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSineMatrixZeros(nla,npola,nspecies,lmax_kkr,lmax_rho,iprint)
!  ===================================================================
   use GauntFactorsModule, only : getK3, getNumK3, getGauntFactor
!
   use RadialGridModule, only : getMaxNumRmesh
!
   use QuadraticMatrixModule, only : initQuadraticMatrix,  &
                                     endQuadraticMatrix,   &
                                     isQuadraticMatrixInitialized
!
   use IntegerFactorsModule, only : initIntegerFactors,               &
                                    isIntegerFactorsInitialized
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nla, npola, iprint
   integer (kind=IntKind), intent(in) :: nspecies(nla)
   integer (kind=IntKind), intent(in) :: lmax_kkr(nla)
   integer (kind=IntKind), intent(in) :: lmax_rho(nla)
   integer (kind=IntKind) :: id, is, n, m, i2, kl1, kl2, kl3, lmax_max
!
   logical, parameter :: isGeneral = .false.
   logical :: yes
!
   LocalNumAtoms = nla
   n_spin_pola = npola
   print_level = iprint
!
   allocate(SineZero(LocalNumAtoms,n_spin_pola))
!
   eGID = getGroupID('Energy Mesh')
   NumPEsInEGroup = getNumPEsInGroup(eGID)
   MyPEinEGroup = getMyPEinGroup(eGID)
   if (iprint >= 0) then
      write(6,'(a,3i5)')'MyPE, MyPEinEGroup, NumPEsInEGroup = ',MyPE,MyPEinEGroup,NumPEsInEGroup
   endif
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   lmax_kkr_max = 0
   lmax_rho_max = 0
   do id = 1,LocalNumAtoms
      lmax_kkr_max = max(lmax_kkr_max,lmax_kkr(id))
      lmax_rho_max = max(lmax_rho_max,lmax_rho(id))
   enddo
   do is = 1, n_spin_pola
      do id = 1,LocalNumAtoms
         SineZero(id,is)%NumSpecies = nspecies(id)
         SineZero(id,is)%jmax_rho = (lmax_rho(id)+1)*(lmax_rho(id)+2)/2
         SineZero(id,is)%kmax_kkr = (lmax_kkr(id)+1)**2
         SineZero(id,is)%NumRs = 0
      enddo
   enddo
!
   kmax_kkr_max = (lmax_kkr_max+1)**2
   kmax_kkr_save = kmax_kkr_max
   kmax_rho_max = (lmax_rho_max+1)**2
   jmax_rho_max = (lmax_rho_max+1)*(lmax_rho_max+2)/2
   lmax_max = max(lmax_kkr_max, lmax_rho_max)
!
   allocate( gaunt(kmax_kkr_max,kmax_kkr_max,kmax_rho_max) )
   gaunt = ZERO
   do kl3 = 1, kmax_rho_max
      do kl1 = 1, kmax_kkr_max
         do i2 = 1, nj3(kl1,kl3)
            kl2 = kj3(i2,kl1,kl3)
            if (kl2 <= kmax_kkr_max) then
               gaunt(kl2,kl1,kl3) = cgnt(i2,kl1,kl3)
            endif
         enddo
      enddo
   enddo
!
!  ===================================================================
!  calculate the charge density associated with each bound state
!  ===================================================================
   MaxNumRs = getMaxNumRmesh()
   allocate( wspace0(kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace1(kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace2(MaxNumRs*kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace3(MaxNumRs*kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace4(MaxNumRs*kmax_kkr_max*kmax_kkr_max) )
!
   do is = 1, n_spin_pola
      do id = 1, LocalNumAtoms
         n = SineZero(id,is)%NumSpecies
         m = SineZero(id,is)%kmax_kkr
         allocate( SineZero(id,is)%NumZeros(n) );    SineZero(id,is)%NumZeros = 0
         allocate( SineZero(id,is)%ZeroIndex(2*m,n) ); SineZero(id,is)%ZeroIndex = 0
         allocate( SineZero(id,is)%ebot(n) ); SineZero(id,is)%ebot = ZERO
         allocate( SineZero(id,is)%etop(n) ); SineZero(id,is)%etop = ZERO
         allocate( SineZero(id,is)%ZeroState(2*m,n) )
      enddo
   enddo
!
   yes = isQuadraticMatrixInitialized(n)
   if (yes .and. n < kmax_kkr_max) then
      call endQuadraticMatrix()
      yes = .false.
   endif
!
   if (.not.yes) then
!     ----------------------------------------------------------------
      call initQuadraticMatrix(kmax_kkr_max,isGeneral)
!     ----------------------------------------------------------------
   endif
   if (.not.isIntegerFactorsInitialized()) then
!     ----------------------------------------------------------------
      call initIntegerFactors(lmax_max)
!     ----------------------------------------------------------------
   endif
!
   isInitialized = .true.
!
   end subroutine initSineMatrixZeros
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSineMatrixZeros()
!  ===================================================================
   use QuadraticMatrixModule, only : endQuadraticMatrix, &
                                     isQuadraticMatrixInitialized
!
   implicit none
!
   integer (kind=IntKind) :: id, ib, ia, is
!
   deallocate( gaunt )
   deallocate( wspace0 )
   deallocate( wspace1 )
   deallocate( wspace2 )
   deallocate( wspace3 )
   deallocate( wspace4 )
!
   do is = 1, n_spin_pola
      do id = 1, LocalNumAtoms
         do ia = 1, SineZero(id,is)%NumSpecies
            do ib = 1, SineZero(id,is)%NumZeros(ia)
               if ( allocated(SineZero(id,is)%ZeroState(ib,ia)%ResidualMat) ) then
                  deallocate( SineZero(id,is)%ZeroState(ib,ia)%ResidualMat )
               endif
               if ( allocated(SineZero(id,is)%ZeroState(ib,ia)%Density) ) then
                  deallocate( SineZero(id,is)%ZeroState(ib,ia)%Density )
               endif
               if ( allocated(SineZero(id,is)%ZeroState(ib,ia)%Deriv_Density) ) then
                  deallocate( SineZero(id,is)%ZeroState(ib,ia)%Deriv_Density )
               endif
            enddo
         enddo
         deallocate(SineZero(id,is)%ZeroIndex)
         deallocate(SineZero(id,is)%ebot)
         deallocate(SineZero(id,is)%etop)
         deallocate(SineZero(id,is)%NumZeros)
         deallocate(SineZero(id,is)%ZeroState)
      enddo
   enddo
   deallocate(SineZero)
!
   if (isQuadraticMatrixInitialized()) then
!     ----------------------------------------------------------------
      call endQuadraticMatrix()
!     ----------------------------------------------------------------
   endif
!
   isInitialized = .false.
!
   end subroutine endSineMatrixZeros
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSineMatrixZerosInitialized() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   y = isInitialized
!
   end function isSineMatrixZerosInitialized
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumSineMatrixZeros(id,ia,is) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
   integer (kind=IntKind) :: n
!
   n = SineZero(id,is)%NumZeros(ia)
!
   end function getNumSineMatrixZeros
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSineMatrixZero(id,ia,is,ips,nd,sorted) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia, ips
   integer (kind=IntKind), intent(out) :: nd
   logical, optional, intent(in) :: sorted
   integer (kind=IntKind) :: ip
!
   real (kind=RealKind) :: p
!
   if (ips < 1 .or. ips > SineZero(id,is)%NumZeros(ia)) then
      call ErrorHandler('getSineMatrixPole','Zero index is out of range',ip)
   endif
!
   ip = ips
   if (present(sorted)) then
      if (sorted) then
         ip = SineZero(id,is)%ZeroIndex(ips,ia)
      endif
   endif
   p = SineZero(id,is)%ZeroState(ip,ia)%ZeroE
   nd = SineZero(id,is)%ZeroState(ip,ia)%NumDegens
!
   end function getSineMatrixZero
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSineZeroDensity(id,ia,is,ibs,NumRs,jmax_rho,derivative,sorted) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is, ibs
   integer (kind=IntKind), intent(out), optional :: NumRs, jmax_rho
   logical, optional, intent(in) :: sorted
   integer (kind=IntKind) :: NumBPs
   integer (kind=IntKind) :: ib
!
   complex (kind=CmplxKind), intent(out), optional, pointer :: derivative(:,:)
   complex (kind=CmplxKind), pointer :: p(:,:)
!
   if (present(NumRs)) then
       NumRs = SineZero(id,is)%NumRs
   endif
!
   if (present(jmax_rho)) then
      jmax_rho = SineZero(id,is)%jmax_rho
   endif
!
   NumBPs = SineZero(id,is)%NumZeros(ia)
   if (NumBPs < 1) then
      nullify(p)
      return
   else if (ibs < 1 .or. ibs > NumBPs) then
      call ErrorHandler('getSineZeroDensity','Invalid sine zero index',ibs)
   endif
!
   ib = ibs
   if (present(sorted)) then
      if (sorted) then
         ib = SineZero(id,is)%ZeroIndex(ibs,ia)
      endif
   endif
!
   if (present(derivative)) then
      derivative => aliasArray2_c(SineZero(id,is)%ZeroState(ib,ia)%Deriv_Density, &
                                  SineZero(id,is)%NumRs,SineZero(id,is)%jmax_rho)
   endif
!
   p => aliasArray2_c(SineZero(id,is)%ZeroState(ib,ia)%Density,                    &
                      SineZero(id,is)%NumRs,SineZero(id,is)%jmax_rho)
!
   end function getSineZeroDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printSineMatrixZerosInfo(id,ia,is)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
   integer (kind=IntKind) :: ip
!
   write(6,'(/,3(a,i3))')'Spin index: ',is,',  Site index: ',id,',  Species index: ',ia
   write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of sine   matrix zeros found within (', &
                                    SineZero(id,is)%ebot(ia),', ',    &
                                    SineZero(id,is)%etop(ia),'): ',   &
                                    SineZero(id,is)%NumZeros(ia)
   do ip = 1, SineZero(id,is)%NumZeros(ia)
      write(6,'(a,i2,5x,a,f20.12)')'Degeneracy = ',SineZero(id,is)%ZeroState(ip,ia)%NumDegens, &
              ', sine   matrix zero energy = ',SineZero(id,is)%ZeroState(ip,ia)%ZeroE
   enddo
!
   end subroutine printSineMatrixZerosInfo
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine findSineMatrixZeros(id,ia,is,eb,et,Delta,EiBound,       &
                                  CheckZeros,PanelOnZero,AccumulationCounts)
!  ===================================================================
   use SSSolverModule, only : solveSingleScattering, getSineMatrix
   use SSSolverModule, only : getSolutionRmeshSize
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix
   use QuadraticMatrixModule, only : solveQuadraticEquation, getEigenValue,   &
                                     solveLinearEquation, getEigenVector,     &
                                     getResidualMatrix
!
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
   integer (kind=IntKind), optional, intent(inout) :: AccumulationCounts
!
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind), intent(in), optional :: Delta
   real (kind=RealKind), intent(in), optional :: EiBound
!
   logical, optional, intent(in) :: CheckZeros
   logical, optional, intent(in) :: PanelOnZero
!
   integer (kind=IntKind) :: ie, iw, kmax_kkr, jmax_rho, NumWindows, info, kl, klp
   integer (kind=IntKind) :: i, j, n, nv, nb, je, ip, nb0, ib, nz0, kmax_sine
   integer (kind=IntKind) :: MyNumWindows
   integer (kind=IntKind), allocatable :: bpdeg(:), degens(:)
   integer (kind=IntKind), allocatable :: nbr(:)
!
   logical :: isZeroInterval = .false.
   logical :: chkzero = .false.
   logical :: found
!
   real (kind=RealKind) :: WindowWidth
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0, err, w, ei
!
   real (kind=RealKind), allocatable :: bpe(:), ebr(:)
   real (kind=RealKind) :: bpe_prev
!
   complex (kind=CmplxKind) :: e, cde, ce0, det
   complex (kind=CmplxKind), pointer :: sine_mat(:,:)
   complex (kind=CmplxKind), pointer :: s0(:,:), s1(:,:), s2(:,:)
   complex (kind=CmplxKind), pointer :: sm(:,:)
   complex (kind=CmplxKind), pointer :: pv(:), evr(:), evl(:), em(:,:)
   complex (kind=CmplxKind), allocatable :: vt(:), diag(:)
   complex (kind=CmplxKind), allocatable :: bmat(:,:)
   complex (kind=CmplxKind), allocatable :: brmat(:,:)
!
   if (eb > et) then
      call ErrorHandler('findSineMatrixPoles','eb > et',eb,et)
   endif
!
   if (present(Delta)) then
      WindowWidth = 2.0d0*Delta
   else
      WindowWidth = 0.01d0
   endif
!
   if (present(EiBound)) then
      ei = EiBound
   else
      ei = 0.002d0
   endif
!
   kmax_kkr = SineZero(id,is)%kmax_kkr
   jmax_rho = SineZero(id,is)%jmax_rho
!
   allocate( bpdeg(2*kmax_kkr), degens(2*kmax_kkr) )
   allocate( bpe(2*kmax_kkr), ebr(2*kmax_kkr) )
   allocate( vt(kmax_kkr), diag(kmax_kkr) )
   allocate( bmat(kmax_kkr*kmax_kkr,2*kmax_kkr) )
   allocate( brmat(kmax_kkr*kmax_kkr,2*kmax_kkr) )
!
   if (.not.present(CheckZeros)) then
      chkzero = .false.
   else
      chkzero = CheckZeros
   endif
!
   if (.not.present(PanelOnZero)) then
      isZeroInterval = .false.
   else
      isZeroInterval = PanelOnZero
   endif
!
   SineZero(id,is)%ebot(ia) = eb
   SineZero(id,is)%etop(ia) = et
!
   if (isZeroInterval) then
      NumWindows = 1
      MyNumWindows = 1
!     de = FOURTH*(et-eb)
      de = HALF*(et-eb)
   else
      NumWindows = ceiling((et-eb)/WindowWidth)
!     NumWindows = NumWindows - mod(NumWindows,4)
      if (NumWindows < NumPEsInEGroup) then
         NumWindows = NumPEsInEGroup
         MyNumWindows = 1
      else
         NumWindows = ceiling(NumWindows/real(NumPEsInEGroup))*NumPEsInEGroup
         MyNumWindows = NumWindows/NumPEsInEGroup
      endif
!     WindowWidth = (et-eb)/real(NumWindows,kind=RealKind)
      if (present(Delta)) then
         de = HALF*Delta
      else
         de = WindowWidth/4.0d0
      endif
   endif
!
   de2 = de*TWO; dede2 = de*de*TWO
   cde = de
!  write(6,'(a,2d15.8,i5)')'de,win,nw = ',de,WindowWidth,NumWindows
!
   s0 => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
   s1 => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
   s2 => aliasArray2_c(wspace2,kmax_kkr,kmax_kkr)
   sm => aliasArray2_c(wspace3,kmax_kkr,kmax_kkr)
!
!  ===================================================================
!  Hopefully, QuadraticMatrixModule could be updated in the future
!  so that kmax_kkr_max is made to be the leading matrix dimension and
!  the following 4 lines of the code become unnecessary
!  ===================================================================
   if (kmax_kkr /= kmax_kkr_save) then
      call endQuadraticMatrix()
      call initQuadraticMatrix(kmax_kkr)
      kmax_kkr_save = kmax_kkr
   endif
!
   nb = 0
   bpe = ZERO
   bpe_prev = ZERO
   bpdeg = 0
   bmat = CZERO
   do iw = 1, MyNumWindows
      w0 = eb + (iw+MyPEInEGroup*MyNumWindows-1)*WindowWidth
      e0 = w0 + (HALF)*WindowWidth
!     write(6,'(a,i3,a,f6.3,a,f6.3,a)')'Window:',iw,'  (',w0,',',w0+WindowWidth,')'
      if (isZeroInterval) then
         e0 = ZERO
      else if ((abs(e0) < Ten2m6 .or. abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6)) then
         if (e0 < ZERO) then
            e0 = e0 - HALF*de
         else
            e0 = e0 + HALF*de
         endif
      endif
      ce0 = e0
!
      if (isZeroInterval) then
         s0(1:kmax_kkr,1:kmax_kkr) = CZERO
      else
         e = cmplx(e0,ZERO,kind=CmplxKind)
!        -------------------------------------------------------------
         call solveSingleScattering(is, id, ce0, CZERO, atom=ia)
!        -------------------------------------------------------------
         sine_mat => getSineMatrix(spin=is,site=id,atom=ia)
!        -------------------------------------------------------------
         kmax_sine = size(sine_mat,dim=1)
         if (kmax_kkr == kmax_sine) then
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sine_mat,1,s0,1)
!           ----------------------------------------------------------
         else if (kmax_kkr < kmax_sine) then
            do kl = 1, kmax_kkr
               do klp = 1, kmax_kkr
                  s0(klp,kl) = sine_mat(klp,kl)
               enddo
            enddo
         else
!           ----------------------------------------------------------
            call ErrorHandler('findSineMatrixZeros','kmax_kkr > kmax_sine',kmax_kkr,kmax_sine)
!           ----------------------------------------------------------
         endif
      endif
!
      e = cmplx(e0+de,ZERO,kind=CmplxKind)
!     ----------------------------------------------------------------
      call solveSingleScattering(is, id, e, CZERO, atom=ia)
!     call solveSingleScattering(is, id, ce0, -cde, atom=ia)
!     ----------------------------------------------------------------
      sine_mat => getSineMatrix(spin=is,site=id,atom=ia)
!     ----------------------------------------------------------------
      kmax_sine = size(sine_mat,dim=1)
      if (kmax_kkr == kmax_sine) then
!        -------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sine_mat,1,s2,1)
!        -------------------------------------------------------------
      else if (kmax_kkr < kmax_sine) then
         do kl = 1, kmax_kkr
            do klp = 1, kmax_kkr
               s2(klp,kl) = sine_mat(klp,kl)
            enddo
         enddo
      else
!        -------------------------------------------------------------
         call ErrorHandler('findSineMatrixZeros','kmax_kkr > kmax_sine',kmax_kkr,kmax_sine)
!        -------------------------------------------------------------
      endif
!
      e = cmplx(e0-de,ZERO,kind=CmplxKind)
!     ----------------------------------------------------------------
      call solveSingleScattering(is, id, e, CZERO, atom=ia)
!     call solveSingleScattering(is, id, ce0, cde, atom=ia)
!     ----------------------------------------------------------------
      sine_mat => getSineMatrix(spin=is,site=id,atom=ia)
      kmax_sine = size(sine_mat,dim=1)
      if (kmax_kkr == kmax_sine) then
!        -------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sine_mat,1,sm,1)
!        -------------------------------------------------------------
      else if (kmax_kkr < kmax_sine) then
         do kl = 1, kmax_kkr
            do klp = 1, kmax_kkr
               sm(klp,kl) = sine_mat(klp,kl)
            enddo
         enddo
      else
!        -------------------------------------------------------------
         call ErrorHandler('findSineMatrixZeros','kmax_kkr > kmax_sine',kmax_kkr,kmax_sine)
!        -------------------------------------------------------------
      endif
!
      s1 = (s2 - sm)/de2
      s2 = (s2 + sm - TWO*s0)/dede2
!
      if (isZeroInterval) then
!        -------------------------------------------------------------
         call solveLinearEquation(s1,s2,info)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call solveQuadraticEquation(s0,s1,s2,info)
!        -------------------------------------------------------------
      endif
!
      if (info /= 0) then
         stop 'Error in s0, s1, s2'
      endif
!
!     ----------------------------------------------------------------
      pv => getEigenValue(nv)
!     ----------------------------------------------------------------
      do ie = 1, nv
         if (.false.) then ! change to .false. to turn off self-checking
!           ==========================================================
!           Check eigenvalues and eigenvectors
!           ==========================================================
            write(6,'(a,i5,2d15.8)')'Index, Eigenvalue = ',ie,pv(ie)
            sm = s0 + s1*pv(ie) + s2*(pv(ie)*pv(ie))
            evr => getEigenVector('R',ie)
            vt = CZERO
            do j = 1, kmax_kkr
               do i = 1, kmax_kkr
                  vt(i) = vt(i) + sm(i,j)*evr(j)
               enddo
            enddo
            do i = 1, kmax_kkr
               err = abs(vt(i))
               if (err > ten2m7) then
                  call ErrorHandler('findSineMatrixPoles','Right-side eigenvector error > 10^7',err)
               endif
            enddo
            write(6,'(a)')'Right-side eigenvector passed!'
!
            evl => getEigenVector('L',ie)
            vt = CZERO
            do j = 1, kmax_kkr
               do i = 1, kmax_kkr
                  vt(j) = vt(j) + evl(i)*sm(i,j)
               enddo
            enddo
            do i = 1, kmax_kkr
               err = abs(vt(i))
               if (err > ten2m7) then
                  call WarningHandler('findSineMatrixZeros','Left-side eigenvector error > 10^7',err)
               endif
            enddo
            write(6,'(a)')'Left-side eigenvector passed!'
         endif
!        =============================================================
         if (abs(aimag(pv(ie))) < Ten2m8) then   ! Zeros on the real energy axis
            pe = real(pv(ie),kind=RealKind) + e0
            if (pe >= w0 .and. pe <= w0+WindowWidth) then
!              -------------------------------------------------------
               em => getResidualMatrix(ie) ! em is the residule matrix of
                                           ! integrating sm^{-1} around its eigenvalue
!              -------------------------------------------------------
               if (size(em,1) /= kmax_kkr) then
                  call ErrorHandler('findSineMatrixPoles','inconsistent matrix size',size(em,1),kmax_kkr)
               endif
               found = .false.
               do i = 1, nb
                  if (abs(pe-bpe(i)) < degen_tol) then
                     bpdeg(i) = bpdeg(i) + 1
!                    -------------------------------------------------
                     call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,bmat(1,i),1)
!                    -------------------------------------------------
                     found = .true.
                     exit
                  endif
               enddo
               if (.not.found) then
                  nb = nb + 1
                  bpe(nb) = pe
                  bpdeg(nb) = 1
!                 ----------------------------------------------------
                  call zcopy(kmax_kkr*kmax_kkr,em,1,bmat(1,nb),1)
!                 ----------------------------------------------------
               endif
            endif
         else if (aimag(pv(ie)) > ZERO .and. aimag(pv(ie)) < ei) then
            pe = real(pv(ie),kind=RealKind) + e0
            if (pe >= w0 .and. pe <= w0+WindowWidth) then
!              call WarningHandler('findSineMatrixPoles','For ei = ',ei)
!              call WarningHandler('findSineMatrixPoles',                &
!                                  'A sine-matrix zero is found at e with 0 < Im[e] < ei', &
!                                  pv(ie)+e0)
            endif
         endif
      enddo
   enddo
!
   if (chkzero) then
      do je = 1, nb
         e0 = bpe(je)
         do ie = -10, 10
            e = e0 + ie*0.001d0
            call solveSingleScattering(is, id, e, CZERO, atom=ia)
            sine_mat => getSineMatrix()
            kmax_sine = size(sine_mat,dim=1)
            call calcDet(sine_mat,kmax_sine,det,diag)
            write(6,'(a,f12.5,2x,2d16.8)')'e,det = ',real(e),det
         enddo
         write(6,'(/)')
      enddo
   endif
!
   allocate(nbr(NumPEsInEGroup))
   nbr = 0
   nbr(MyPEinEGroup+1) = nb
!  -------------------------------------------------------------------
   call GlobalSumInGroup(eGID,nbr,NumPEsInEGroup)
!  -------------------------------------------------------------------
   if (present(AccumulationCounts)) then
      SineZero(id,is)%NumZeros(ia) = AccumulationCounts
   else
      SineZero(id,is)%NumZeros(ia) = 0
   endif
   nz0 = SineZero(id,is)%NumZeros(ia)
   do ip = 1, NumPEsInEGroup
      SineZero(id,is)%NumZeros(ia) = SineZero(id,is)%NumZeros(ia) + nbr(ip)
   enddo
   if (SineZero(id,is)%NumZeros(ia) > 2*SineZero(id,is)%kmax_kkr) then
      call ErrorHandler('findSineMatrixPoles','NumZeros > 2*kmax_kkr', &
                        SineZero(id,is)%NumZeros(ia), 2*SineZero(id,is)%kmax_kkr)
   endif
!
   do n = nz0 + 1, SineZero(id,is)%NumZeros(ia)
      if ( .not.allocated(SineZero(id,is)%ZeroState(n,ia)%ResidualMat) ) then
         allocate( SineZero(id,is)%ZeroState(n,ia)%ResidualMat(kmax_kkr*kmax_kkr) )
      endif
      SineZero(id,is)%ZeroState(n,ia)%ResidualMat = CZERO
!
      if ( .not.allocated(SineZero(id,is)%ZeroState(n,ia)%Density) ) then
         allocate( SineZero(id,is)%ZeroState(n,ia)%Density(MaxNumRs*jmax_rho) )
      endif
      SineZero(id,is)%ZeroState(n,ia)%Density = CZERO
!
      if ( .not.allocated(SineZero(id,is)%ZeroState(n,ia)%Deriv_Density) ) then
         allocate( SineZero(id,is)%ZeroState(n,ia)%Deriv_Density(MaxNumRs*jmax_rho) )
      endif
      SineZero(id,is)%ZeroState(n,ia)%Deriv_Density = CZERO
   enddo
!
   ebr = ZERO
   nb0 = nz0
   do ip = 1, NumPEsInEGroup
      if (nbr(ip) > 0) then
         if (MyPEinEGroup == ip-1) then
            ebr(1:nbr(ip)) = bpe(1:nbr(ip))
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr*nbr(ip),bmat,1,brmat,1)
!           ----------------------------------------------------------
            degens(1:nbr(ip)) = bpdeg(1:nbr(ip))
         endif
!        -------------------------------------------------------------
         call bcastMessageInGroup(eGID,ebr,nbr(ip),ip-1)
         call bcastMessageInGroup(eGID,degens,nbr(ip),ip-1)
         call bcastMessageInGroup(eGID,brmat,kmax_kkr*kmax_kkr,nbr(ip),ip-1)
!        -------------------------------------------------------------
         do ib = 1, nbr(ip)
            SineZero(id,is)%ZeroState(nb0+ib,ia)%ZeroE = ebr(ib)
            SineZero(id,is)%ZeroState(nb0+ib,ia)%NumDegens = degens(ib)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,brmat(1,ib),1,               &
                       SineZero(id,is)%ZeroState(nb0+ib,ia)%ResidualMat,1)
!           ----------------------------------------------------------
         enddo
         nb0 = nb0 + nbr(ip)
      endif
   enddo
!
   SineZero(id,is)%NumRs = getSolutionRmeshSize()
!
!  -------------------------------------------------------------------
   call sortZeros(SineZero(id,is)%NumZeros(ia),                       &
                  SineZero(id,is)%ZeroState(:,ia),                    &
                  SineZero(id,is)%ZeroIndex(:,ia) )
!  -------------------------------------------------------------------
!
   deallocate(vt, diag)
   deallocate(bmat)
   deallocate(brmat)
   deallocate(bpdeg, degens)
   deallocate(bpe, ebr)
   deallocate(nbr)
   nullify(pv, evl, evr, em)
   nullify(s0, s1, s2, sm)
!
   end subroutine findSineMatrixZeros
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSineZeroInEnergyRange(id,ia,is,e1,e2) result(y)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ia
   integer (kind=IntKind) :: ip
!
   real (kind=RealKind), intent(in) :: e1, e2
   real (kind=RealKind) ::  e_low, e_high
!
   logical :: y
!
   if (e1 < e2) then
      e_low = e1
      e_high = e2
   else
      e_low = e2
      e_high = e1
   endif
!
   y = .false.
   do ip = 1, SineZero(id,is)%NumZeros(ia)
      if (SineZero(id,is)%ZeroState(ip,ia)%ZeroE <= e_high .and.      &
          e_low <= SineZero(id,is)%ZeroState(ip,ia)%ZeroE) then
         y = .true.
         exit
      endif
   enddo
!
   end function isSineZeroInEnergyRange
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeSineZeroDensity(id,ia,is)
!  ===================================================================
!  Note: The density associated with the bound state pole includes the
!        degeneracy of the pole, since the residual matrix has already
!        been multiplied by the number of degeneracies.
!  *******************************************************************
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : getRegSolution, getRegSolutionDerivative
   use SSSolverModule, only : getJostInvMatrix
!  use SSSolverModule, only : getOmegaHatMatrix
!  use SSSolverModule, only : computeDOS, getDOS
!
   use RadialGridModule, only : getGrid
!
   use MatrixModule, only : computeAStarT
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk, bofk, m1m
!
   use StepFunctionModule, only : getVolumeIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ia
!
   integer (kind=IntKind) :: ie, ib, info, ip
   integer (kind=IntKind) :: kmax_kkr, jmax_rho, kmax_rho, NumRs, NumBPs, kmax_phi
   integer (kind=IntKind) :: kl, klp, klp_bar, kl1, kl2, kl3, kl3_bar, m3, mp, ir, jl3
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: Deriv_Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: smr(:,:), BSinv(:,:)
   complex (kind=CmplxKind), pointer :: PhiLr(:,:,:), DerPhiLr(:,:,:)
   complex (kind=CmplxKind), pointer :: BPhiLr(:,:,:), DerBPhiLr(:,:,:), PPr(:,:,:)
   complex (kind=CmplxKind), pointer :: jost_inv(:,:)
!  complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind) :: e, cfac, cfac0, cfac1, cfac2, kappa
!
   type (GridStruct), pointer :: Grid
!
   logical :: non_zero
!
   kmax_kkr = SineZero(id,is)%kmax_kkr
   jmax_rho = SineZero(id,is)%jmax_rho
   kmax_rho = kofj(jmax_rho)
   NumRs = SineZero(id,is)%NumRs
   NumBPs = SineZero(id,is)%NumZeros(ia)
!
   if (NumBPs < 1) then
      return
   endif
!
!  ===================================================================
!  Note: Bdensity stores the density multiplied by the number of degeneracies
!        of the bound state.
!  ===================================================================
   smr => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
   BSinv => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
!
   Grid => getGrid(id)
   r_mesh => Grid%r_mesh
   BPhiLr => aliasArray3_c(wspace2,NumRs,kmax_kkr,kmax_kkr)
   PPr => aliasArray3_c(wspace3,NumRs,kmax_kkr,kmax_kkr)
   DerBPhiLr => aliasArray3_c(wspace4,NumRs,kmax_kkr,kmax_kkr)
!
   do ib = 1, NumBPs
      SineZero(id,is)%ZeroState(ib,ia)%Density = CZERO
      SineZero(id,is)%ZeroState(ib,ia)%Deriv_Density = CZERO
   enddo
!
   do ib = MyPEinEGroup+1, NumBPs, NumPEsInEGroup
      Bdensity => aliasArray2_c(SineZero(id,is)%ZeroState(ib,ia)%Density,NumRs,jmax_rho)
      Deriv_Bdensity => aliasArray2_c(SineZero(id,is)%ZeroState(ib,ia)%Deriv_Density,NumRs,jmax_rho)
      e = SineZero(id,is)%ZeroState(ib,ia)%ZeroE
      kappa = sqrt(e)
      cfac0 = HALF*kappa
!     ----------------------------------------------------------------
      call solveSingleScattering(is, id, e, CZERO, atom=ia)
!     ----------------------------------------------------------------
      jost_inv => getJostInvMatrix()
!     ================================================================
!     calculate (residual of sine_mat)^(T*) and store the result in smr
!     ----------------------------------------------------------------
      call computeAStarT(SineZero(id,is)%ZeroState(ib,ia)%ResidualMat,kmax_kkr,kmax_kkr,smr)
!     ----------------------------------------------------------------
!
!     ================================================================
!     calculate jost_inv*ResidualMat^{T*} and store the result in BSinv
!     ----------------------------------------------------------------
      call zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
                 jost_inv,kmax_kkr,smr,kmax_kkr,CZERO,BSinv,kmax_kkr)
!     ----------------------------------------------------------------
!call !zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,jost_inv,kmax_kkr,smat_inv,kmax_kkr,CZERO,BSinv,kmax_kkr)
!
      PhiLr => getRegSolution()
      kmax_phi = size(PhiLr,dim=2)
!     ----------------------------------------------------------------
      call zgemm('n','n',NumRs*kmax_kkr,kmax_kkr,kmax_kkr,CONE,          &
                 PhiLr,NumRs*kmax_phi,BSinv,kmax_kkr,                    &
                 CZERO,BPhiLr,NumRs*kmax_kkr)
!     ----------------------------------------------------------------
      PPr = CZERO
      do klp = 1, kmax_kkr
         mp = mofk(klp)
         klp_bar = bofk(klp)
         cfac = m1m(mp)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               do ir = 1, NumRs
                  PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*BPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar)
               enddo
            enddo
         enddo
      enddo
      do jl3 = 1, jmax_rho
         m3 = mofj(jl3)
         kl3 = kofj(jl3)
         kl3_bar = bofk(kl3)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               cfac1 = cfac0*gaunt(kl1,kl2,kl3)
               cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
               do ir = 1, NumRs
                  Bdensity(ir,jl3) = Bdensity(ir,jl3)           &
                                      + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
!+ cfac1*PPr(ir,kl1,kl2) - conjg(cfac2*PPr(ir,kl1,kl2))
               enddo
            enddo
         enddo
      enddo
!
      DerPhiLr => getRegSolutionDerivative()
!     ----------------------------------------------------------------
      call zgemm('n','n',NumRs*kmax_kkr,kmax_kkr,kmax_kkr,CONE,       &
                 DerPhiLr,NumRs*kmax_phi,BSinv,kmax_kkr,              &
                 CZERO,DerBPhiLr,NumRs*kmax_kkr)
!     ----------------------------------------------------------------
      PPr = CZERO
      do klp = 1, kmax_kkr
         mp = mofk(klp)
         klp_bar = bofk(klp)
         cfac = m1m(mp)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               do ir = 1, NumRs
                  PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*DerBPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar) &
                                                    + cfac*BPhiLr(ir,kl1,klp)*DerPhiLr(ir,kl2,klp_bar)
               enddo
            enddo
         enddo
      enddo
      do jl3 = 1, jmax_rho
         m3 = mofj(jl3)
         kl3 = kofj(jl3)
         kl3_bar = bofk(kl3)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               cfac1 = cfac0*gaunt(kl1,kl2,kl3)
               cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
               do ir = 1, NumRs
                  Deriv_Bdensity(ir,jl3) = Deriv_Bdensity(ir,jl3)  &
                                            + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
!+ cfac1*PPr(ir,kl1,kl2) - conjg(cfac2*PPr(ir,kl1,kl2))
               enddo
            enddo
         enddo
      enddo
!
!     ================================================================
!     Get rid of r^2 from Bdensity and Deriv_Bdensity
!     ================================================================
      do jl3 = 1, jmax_rho
         do ir = 1, NumRs
            Bdensity(ir,jl3) = Bdensity(ir,jl3)/r_mesh(ir)**2
         enddo
         do ir = 1, NumRs
            Deriv_Bdensity(ir,jl3) = Deriv_Bdensity(ir,jl3)/r_mesh(ir)**2
         enddo
      enddo
   enddo
!
   do ib = 1, NumBPs
      ip = mod(ib-1,NumPEsInEGroup)
!     ---------------------------------------------------------------
      call bcastMessageInGroup(eGID,SineZero(id,is)%ZeroState(ib,ia)%Density,NumRs*jmax_rho,ip)
      call bcastMessageInGroup(eGID,SineZero(id,is)%ZeroState(ib,ia)%Deriv_Density,NumRs*jmax_rho,ip)
!     ---------------------------------------------------------------
      Bdensity => aliasArray2_c(SineZero(id,is)%ZeroState(ib,ia)%Density,NumRs,jmax_rho)
   enddo
!
   nullify(jost_inv, BSinv, smr, Grid, r_mesh)
   nullify(BPhiLr, PhiLr, DerBPhiLr, DerPhiLr, PPr)
   nullify(Bdensity, Deriv_Bdensity)
!
   end subroutine computeSineZeroDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calcDet(mat,kmax,det,diag)
!  ===================================================================
   use KindParamModule, only : RealKind, CmplxKind, IntKind
!
   use MathParamModule, only : CONE
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kmax
   integer (kind=IntKind) :: kl
!
   complex (kind=CmplxKind), intent(out) :: det
   complex (kind=CmplxKind), intent(out) :: diag(kmax)
   complex (kind=CmplxKind), intent(in) :: mat(kmax,kmax)
   complex (kind=CmplxKind) :: matU(kmax,kmax)
!
!  ----------------------------------------------------------
   call zcopy(kmax*kmax,mat,1,matU,1)
!  ----------------------------------------------------------
   call GaussianElim(matU,kmax)
!  ----------------------------------------------------------
   det = CONE
   do kl = 1,kmax
      det = det*matU(kl,kl)
      diag(kl) = matU(kl,kl)
   enddo
!
   end subroutine calcDet
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sortZeros(nz,zs,idx)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nz
   integer (kind=IntKind), intent(out) :: idx(nz)
   integer (kind=IntKind) :: i, j, ki, kj
!
   type (ZeroDensityStruct), intent(in) :: zs(nz)
!
   do i = 1, nz
      idx(i) = i
   enddo
!
   do j = 1, nz-1
      kj = idx(j)
      do i = j+1, nz
         ki = idx(i)
         if (zs(ki)%ZeroE < zs(kj)%ZeroE) then
            idx(i) = kj
            kj = ki
         endif
      enddo
      idx(j) = kj
   enddo
!
   end subroutine sortZeros
!  ===================================================================
end module SineMatrixZerosModule
