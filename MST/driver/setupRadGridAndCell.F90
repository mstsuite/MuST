subroutine setupRadGridAndCell(LocalNumAtoms,lmax_max)
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8, TEN2m6, ONE
!
   use ErrorHandlerModule, only : WarningHandler, ErrorHandler
!
   use ScfDataModule, only : getSingleSiteSolverType
   use ScfDataModule, only : ngaussr, ngaussq
!
   use RadialGridModule, only : initRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
   use RadialGridModule, only : getRadialGridRadius
!
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
   use PolyhedraModule, only : printPolyhedron
   use PolyhedraModule, only : getWignerSeitzRadius, getNeighborDistance
!
   use PotentialTypeModule, only : isASAPotential, isMuffinTinPotential, &
                                   isMuffinTinASAPotential, isFullPotential, &
                                   isMuffinTinTestPotential
!
   use StepFunctionModule, only : initStepFunction
   use StepFunctionModule, only : printStepFunction, testStepFunction
!
   use AtomModule, only : getRadialGridData, getMuffinTinRadius
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax, getStepFuncLmax
   use AtomModule, only : getAtomCoreRad, setAtomCoreRad, setMuffinTinRadius
!
   use OutputModule, only : getStandardOutputLevel
!
   use SystemVolumeModule, only : setSystemVolumeMT, updateSystemVolume, &
                                  printSystemVolume
!
   implicit none
!
   logical :: isRmtUpdated
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms, lmax_max
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:), lmax_step(:)
!
   real (kind=RealKind) :: rmt, rend, rws, rinsc, hin, rmt_grid, rc
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!
!  ===================================================================
!  initialize radial grid
!  -------------------------------------------------------------------
   call initRadialGrid(LocalNumAtoms, 'main', 0)
!  -------------------------------------------------------------------
!
   isRmtUpdated = .false.
   do i=1,LocalNumAtoms
!     ----------------------------------------------------------------
      call getRadialGridData(i,ndivin,ndivout,nmult,hin)
!     ----------------------------------------------------------------
      if (getStandardOutputLevel(i) >= 0) then
!        -------------------------------------------------------------
         call printPolyhedron(i)
!        -------------------------------------------------------------
      endif
      rend = max(getOutscrSphRadius(i),getAtomCoreRad(i))
      rmt = getMuffinTinRadius(i)
      if ( rmt < 0.010d0 ) then
         call ErrorHandler('main','rmt < 0.01',rmt)
      endif
      if (isMuffinTinPotential() .or. isMuffinTinTestPotential()) then
         rinsc = getInscrSphRadius(i)
         if (rmt - rinsc > TEN2m6) then
            call WarningHandler('main','rinsc < rmt',rinsc,rmt)
         endif
         rws = getWignerSeitzRadius(i)
!         volume = PI4*THIRD*(rws**3)
!         call setAtomVolWS(i,volume)
!         volume = PI4*THIRD*(rmt**3)
!         call setAtomVolMT(i,volume)
!         volume = PI4*THIRD*(rinsc**3)
!         call setAtomVolINSC(i,volume)
!         volume = getAtomicVPVolume(i)
!         call setAtomVolVP(i,volume)
!         call genRadialGrid( i, rinsc, rinsc, rws, ndivin, ndivout, nmult)
!
!
!        ========================================================
!???     Temporary fix due to inconsistent definitions of pzz/pzj
!???     between the relativistic and non-relativistic solvers
!        ========================================================
         if (getSingleSiteSolverType()==1) then
            rend=rws
         endif
!        ========================================================
!        -------------------------------------------------------------
!011820  call genRadialGrid( i, xstart, rmt, rinsc, rws, rend, ndivin)
!        call genRadialGrid( i, rmt, rinsc, rend, ndivin)
!081220  call genRadialGrid( i, rmt, rinsc, rend, ndivin = ndivin)
!        -------------------------------------------------------------
      else if (isASAPotential() ) then
!        rend =  getWignerSeitzRadius(i)
         rinsc = getWignerSeitzRadius(i)
         rmt = rinsc
!        if ( rmt < 0.010d0 ) then
!           rmt = rinsc
!        endif
         ndivin = ndivin+11
!        -------------------------------------------------------------
!        call genRadialGrid( i, xstart, rmt, rinsc, rinsc, rend, ndivin )
!081220  call genRadialGrid( i, rmt, rinsc, rend, ndivin = ndivin )
!        -------------------------------------------------------------
      else if ( isMuffinTinASAPotential() ) then
         rend =  getWignerSeitzRadius(i)
         rinsc = getWignerSeitzRadius(i)
!         volume = PI4*THIRD*(rend**3)
!         call setAtomVolWS(i,volume)
!         volume = PI4*THIRD*(rmt**3)
!         call setAtomVolMT(i,volume)
!         volume = PI4*THIRD*(rinsc**3)
!         call setAtomVolINSC(i,volume)
!         volume = getAtomicVPVolume(i)
!         call setAtomVolVP(i,volume)
!         call genRadialGrid( i, rmt, rinsc, rinsc, ndivin, ndivout, nmult)
!        -------------------------------------------------------------
!        call genRadialGrid( i, xstart, rmt, rinsc, rend, rend, ndivin )
!081220  call genRadialGrid( i, rmt, rinsc, rend, ndivin = ndivin )
!        -------------------------------------------------------------
      else
         if (getNeighborDistance(i,1)-getOutscrSphRadius(i) < TEN2m8) then
!           ----------------------------------------------------------
            call WarningHandler('setupRadGridAndCell',                &
                                'Ill condition found: Neighbor distance <= Rcs', &
                                getNeighborDistance(i,1),getOutscrSphRadius(i))
!           ----------------------------------------------------------
         endif
         rinsc = getInscrSphRadius(i)
         rws = getWignerSeitzRadius(i)
!        call genRadialGrid( i, rmt, rinsc, rend, ndivin)
!081220  call genRadialGrid( i, rmt, rinsc, rend, ndivin = ndivin)
!        if ( rmt < 0.010d0 ) then
!           rmt = getInscrSphRadius(i)
!        endif
!        rws = getWignerSeitzRadius(i)
!        -------------------------------------------------------------
!        call genRadialGrid( i, rmt, rinsc, rws, rend, ndivin, ndivout, nmult)
!        -------------------------------------------------------------
      endif
      rc = getAtomCoreRad(i)
      if (hin > TEN2m6 .or. rc < ONE) then
!        -------------------------------------------------------------
         call genRadialGrid(id=i, rmt=rmt, rinsc=rinsc, rend=rend,    &
                            ndivin=ndivin, xstep=hin)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call genRadialGrid(id=i, rmt=rmt, rinsc=rinsc, rend=rend,    &
                            ndivin=ndivin, rfix=rc)
!        -------------------------------------------------------------
      endif
      if (getStandardOutputLevel(i) >= 0) then
!        -------------------------------------------------------------
         call printRadialGrid(i)
!        -------------------------------------------------------------
      endif
      rmt_grid = getRadialGridRadius(i,MT=.true.)
      if (abs(rmt_grid-rmt) > TEN2m6) then
!        -------------------------------------------------------------
         call setMuffinTinRadius(i,rmt_grid)
!        -------------------------------------------------------------
         if (.not.isFullPotential()) then
!           ==========================================================
!           For Muffin-tin, ASA, or Muffin-tin-ASA calculations, since
!           potential outside rmt is 0, core radius is set to rmt
!           ----------------------------------------------------------
            call setAtomCoreRad(i,rmt_grid)
!           ----------------------------------------------------------
         endif
!        -------------------------------------------------------------
         call setSystemVolumeMT(i,rmt_grid)
!        -------------------------------------------------------------
         isRmtUpdated = .true.
      endif
   enddo
   if ( isRmtUpdated ) then
!     ----------------------------------------------------------------
      call updateSystemVolume()
!     ----------------------------------------------------------------
      call printSystemVolume()
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  initialize step function module
!  ===================================================================
   allocate( ngr(LocalNumAtoms), ngt(LocalNumAtoms), lmax_step(LocalNumAtoms) )
   do i=1,LocalNumAtoms
      ngr(i) = ngaussr
      ngt(i) = ngaussq
   enddo
!
   do i=1,LocalNumAtoms
      lmax_step(i)  = getStepFuncLmax(i)
   enddo
!
!  -------------------------------------------------------------------
   call initStepFunction(LocalNumAtoms, lmax_max, lmax_step, ngr, ngt, &
                         'main', 0)
!  -------------------------------------------------------------------
   deallocate( ngr, ngt )
!
   do i=1,LocalNumAtoms
      if (getStandardOutputLevel(i) >= 0) then
!        -------------------------------------------------------------
         call printStepFunction(i)
!        -------------------------------------------------------------
         call testStepFunction(i)
!        -------------------------------------------------------------
      endif
   enddo
!
end subroutine setupRadGridAndCell
