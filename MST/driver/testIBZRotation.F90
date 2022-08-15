program testIBZRotation
!  ********************************************************************
!  main to test the Rotation Matrix code.
!      lmax_kkr = 5 fails tol = TEN2m11 test for the structure constant rotation
!      lmax_kkr = 6 fails tol = TEN2m10 test for the structure constant rotation
!  ********************************************************************
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use ChemElementModule, only : MaxLenOfAtomName
!
   use SystemModule, only : getNumAtoms, getBravaisLattice, getAtomPosition, &
                            getAtomName, getAtomicNumber
!
   use ScfDataModule, only : NumKMeshs, Symmetrize, isLSMS
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop
!
   use TimerModule, only : initTimer, getTime
!
   use MathParamModule, only : ZERO, HALF, ONE, CZERO, CONE, SQRTm1, TEN2m6, &
                               TEN2m7, TEN2m8, TEN2m9, TEN2m10, TEN2m11, TEN2m12
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
   use SphericalHarmonicsModule, only : calylm
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
   use WriteMatrixModule,  only : writeMatrix
!
   use BZoneModule, only : getNumKs, getAllKPoints, printBZone, getWeight, getWeightSum
!
   use IBZRotationModule, only : initIBZRotation, endIBZRotation,     &
                                 computeRotationMatrix, printIBZRotationMatrix, &
                                 getNumIBZRotations, isProperRotation,&
                                 getIBZRotationMatrix, getIBZRotationMatrix3D
   use IBZRotationModule, only : checkCrystalSymmetry
   use IBZRotationModule, only : getBasisRotationTable
!
   use AtomModule, only : getPhiLmax
   use AtomModule, only : getLocalNumSpecies, getLocalAtomicNumber
   use AtomModule, only : getStepFuncLmax, getTruncPotLmax
   use AtomModule, only : getPotLmax, getKKRLmax, getRhoLmax
!
   use LatticeModule, only : initLattice, endLattice, getLatticeType
! 
   use IntegerFactorsModule, only : lofk, mofk, lofj, mofj, kofj, jofk, m1m
!
   use StrConstModule, only : initStrConst, endStrConst
   use StrConstModule, only : getStrConstMatrix
!
   use MatrixModule, only : computeUAUt, computeUAUtc, computeUAU, computeUAUts
!
   use Atom2ProcModule, only : getGlobalIndex, getLocalNumAtoms
!
   use OutputModule, only : getStandardOutputLevel
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,   &
                                    calSymmetryFlags
!
   use PotentialTypeModule, only : isFullPotential
!
   use MadelungModule, only : initMadelung, endMadelung
!
   use PotentialModule, only : initPotential, endPotential, getPotential
   use PotentialModule, only : readPotential, getTruncatedPotential
!
   use PolyhedraModule, only : getInscrSphRadius
   use PolyhedraModule, only : getOutscrSphRadius
!
   use StepFunctionModule, only : getStepFunction
!
   use SSSolverModule, only : initSSSolver, endSSSolver, solveSingleScattering
   use SSSolverModule, only : getSineMatrix, getCosineMatrix, getTMatrix, getJostMatrix
   use SSSolverModule, only : getRegSolution
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use RadialGridModule, only : getGrid
!
   implicit   none
!
   character (len=MaxLenOfAtomName), allocatable :: AtomName(:)
!
   logical :: redundant, failed, sin_failed, cos_failed, jost_failed
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: NumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: lmax_phi, kmax_phi, lmax_max, kmax_kkr
   integer (kind=IntKind) :: lmax_pot_max, lmax_rho_max, lmax_kkr_max, lmax_trunc_pot_max
   integer (kind=IntKind) :: lmax_step_max, jmax_step, jmax_kkr, jmax_pot
   integer (kind=IntKind) :: m1, i1, kl1, jl1
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_green(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: lmax_trunc_pot(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)

!
   integer (kind=IntKind) :: i, j, l, m, n, ij, is, id, ig, kl, klp, mp
   integer (kind=IntKind) :: k, nk, nr, ir, jr, ia, ja, nw, nd, jl, jlp
   integer (kind=IntKind), parameter :: MaxRotations = 48
   integer (kind=IntKind), pointer :: rotation_table(:,:)
   integer (kind=IntKind), allocatable :: nshift(:,:,:)
!
   real (kind=RealKind) :: t0, kfac, kr, tw, sumw, r
   real (kind=RealKind) :: rot3d(3,3), bravais(3,3), Ri(3), Rj(3)
   real (kind=RealKind) :: kin(3), kin_new(3), aij(3), krot(3,MaxRotations)
   real (kind=RealKind), pointer :: kvec(:,:), cgnt(:,:,:)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), parameter :: tol = TEN2m7
!
   complex (kind=CmplxKind) :: energy, kappa, cfac
   complex (kind=CmplxKind), pointer :: rotmat(:,:), rotmatc(:,:)
   complex (kind=CmplxKind), pointer :: strcon(:,:)
   complex (kind=CmplxKind), pointer :: jost_mat(:,:), sin_mat(:,:), cos_mat(:,:)
   complex (kind=CmplxKind), pointer :: pot_jl(:,:), PhiLr(:,:,:)
   complex (kind=CmplxKind), allocatable :: pot_ori(:), pot_rot(:)
   complex (kind=CmplxKind), allocatable :: trunc_pot_ori(:), trunc_pot_rot(:)
   complex (kind=CmplxKind), allocatable :: stepf(:), stepf_rot(:)
   complex (kind=CmplxKind), allocatable :: stepfLL(:,:), stepfLL_rot(:,:)
   complex (kind=CmplxKind), allocatable :: strcon_new(:,:)
   complex (kind=CmplxKind), allocatable :: strcon_rot(:,:), strcon_ori(:,:)
   complex (kind=CmplxKind), allocatable :: phi_rot(:,:), phi_ori(:,:)
   complex (kind=CmplxKind), allocatable :: sin_rot(:,:), sin_ori(:,:)
   complex (kind=CmplxKind), allocatable :: cos_rot(:,:), cos_ori(:,:)
   complex (kind=CmplxKind), allocatable :: jost_rot(:,:), jost_ori(:,:)
   complex (kind=CmplxKind), allocatable :: ylm(:), ylm_new(:), ylm_rot(:)
   complex (kind=CmplxKind), allocatable :: emat(:,:)
   complex (kind=CmplxKind), allocatable :: WORK(:), us(:)
!
   type (GridStruct), pointer :: Grid
!
!  *******************************************************************
!
!  -------------------------------------------------------------------
   call startProcess()
   if (isLSMS()) then
      call ErrorHandler('testIBZRotation','Needs to set method to be KKR in the input file')
   endif
   call printBZone()
   bravais = getBravaisLattice()
   call initLattice(bravais)
!  -------------------------------------------------------------------
!
!  ===================================================================
!
   NumAtoms = getNumAtoms()
   allocate(AtomPosition(1:3,1:NumAtoms), AtomName(NumAtoms))
!
   lmax_phi = 0
   do ig=1,NumAtoms
      AtomName(ig) = getAtomName(ig)
      AtomPosition(1:3,ig)=getAtomPosition(ig)
      lmax_phi = max(lmax_phi,getPhiLmax(ig))
   enddo
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_phi)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(lmax_phi,istop,iprint)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize kmax:
!  ===================================================================
   kmax_phi=(lmax_phi+1)**2
   allocate(ylm(kmax_phi), ylm_new(kmax_phi), ylm_rot(kmax_phi))
   allocate(emat(kmax_phi,kmax_phi), strcon_rot(kmax_phi,kmax_phi))
   allocate(WORK(2*kmax_phi*kmax_phi), strcon_ori(kmax_phi,kmax_phi))
   allocate(us(kmax_phi), strcon_new(kmax_phi,kmax_phi))
!
!  --------------------------------------------------------------------
   call initTimer()
!  call initIBZRotation(.false.,getLatticeType(),lmax_phi,Symmetrize)
   call initIBZRotation(.false.,getLatticeType(),lmax_phi,1)
   call computeRotationMatrix(bravais,NumAtoms,AtomPosition,aname=AtomName)
   call printIBZRotationMatrix(Rot3D_Only=.true.)
   if( checkCrystalSymmetry(bravais,NumAtoms,AtomPosition,aname=AtomName) ) then
      write(6,'(/,a,/)')'The crystal system does have the point group symmetry!'
   else
      call ErrorHandler('testIBZRotation','The crystal system does not have the point group symmetry')
   endif
   nr = getNumIBZRotations()
!
!  ===================================================================
!  Test: rotmat is unitary matrix
!  ===================================================================
   do ir = 1, nr
      rot3d = getIBZRotationMatrix3D(ir)
      rotmat => getIBZRotationMatrix('n',ir)
      rotmatc => getIBZRotationMatrix('c',ir)
      do kl = 1, kmax_phi
         do klp = 1, kmax_phi
            if (lofk(kl) /= lofk(klp) .and. abs(rotmat(klp,kl)) > tol) then
               write(6,'(a,2i4,2x,2i4,2x,2d15.8)')'Coupling between l <> lp element:', &
               lofk(klp),mofk(klp),lofk(kl),mofk(kl),rotmat(klp,kl)
               call writeMatrix('Rotation matrix',rotmat,kmax_phi,kmax_phi,tol)
               call ErrorHandler('testIBZRotation','Rotation matrix test failed!')
            endif
         enddo
      enddo
!     ----------------------------------------------------------------
      call zgemm('n','t',kmax_phi,kmax_phi,kmax_phi,CONE,rotmat,   &
                 kmax_phi,rotmatc,kmax_phi,CZERO,emat,kmax_phi)
!     ----------------------------------------------------------------
      do j = 1, kmax_phi
         do i = 1, kmax_phi
            if ((i == j .and. abs(emat(i,j)-CONE) > tol) .or.      &
                (i /= j .and. abs(emat(i,j)) > tol)) then
               write(6,'(/,a,i5)')'Rotmat is not unitary: rotation index = ',ir
               write(6,'(a)')'3D rotation matrix:'
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
               call writeMatrix('Rotation matrix',rotmat,kmax_phi,kmax_phi,tol)
               call ErrorHandler('testIBZRotation','Unitary matrix test 1 failed!')
            endif
         enddo
      enddo
!     ----------------------------------------------------------------
      call zgemm('t','n',kmax_phi,kmax_phi,kmax_phi,CONE,rotmatc,   &
                 kmax_phi,rotmat,kmax_phi,CZERO,emat,kmax_phi)
!     ----------------------------------------------------------------
      do j = 1, kmax_phi
         do i = 1, kmax_phi
            if ((i == j .and. abs(emat(i,j)-CONE) > tol) .or.      &
                (i /= j .and. abs(emat(i,j)) > tol)) then
               write(6,'(/,a,i5)')'Rotmat is not unitary: rotation index = ',ir
               write(6,'(a)')'3D rotation matrix:'
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
               call writeMatrix('Rotation matrix',rotmat,kmax_phi,kmax_phi,tol)
               call ErrorHandler('testIBZRotation','Unitary matrix test 2 failed!')
            endif
         enddo
      enddo
   enddo
   write(6,'(/,1x,a)')'Passed rotation matrix test...'
!
   nk = getNumKs()
   kvec => getAllKPoints(kfac)
!
   allocate(nshift(3,NumAtoms,nr))
   rotation_table => getBasisRotationTable(ntab=nshift)
!
   if (nr > MaxRotations) then
      call ErrorHandler('testBZone','Number of crystal rotations exceeds its physical limit',nr,MaxRotations)
   endif
!
   tw = ZERO ! normalization factor for the weight
   sumw = ZERO
   do k = 1, nk
      nw = 0
      do ir = 1, nr
         rot3d = getIBZRotationMatrix3D(ir)
         krot(1,ir) = rot3d(1,1)*kvec(1,k)+rot3d(1,2)*kvec(2,k)+rot3d(1,3)*kvec(3,k)
         krot(2,ir) = rot3d(2,1)*kvec(1,k)+rot3d(2,2)*kvec(2,k)+rot3d(2,3)*kvec(3,k)
         krot(3,ir) = rot3d(3,1)*kvec(1,k)+rot3d(3,2)*kvec(2,k)+rot3d(3,3)*kvec(3,k)
         if (abs(krot(1,ir)-kvec(1,k))+abs(krot(2,ir)-kvec(2,k))+abs(krot(3,ir)-kvec(3,k)) < TEN2m7) then
            nw = nw + 1
         endif
      enddo
      tw = tw + ONE/real(nw,kind=RealKind)
      sumw = sumw + getWeight(k)
   enddo
!
!  ===================================================================
!  Check the redundancy of each k-point due to rotation operations.
!  -------------------------------------------------------------------
   write(6,'(/,a,i5)')'Number of rotations: ',nr
   write(6,'(a)')'Check the redundancy of each k-point due to rotations...'
   write(6,'(a)')'Note: the redundancy divided by the number of rotations equals the weight.'
   write(6,'(a)')'For symmetrized calculation with IBZ rotations, the k-points generated are in IBZ.'
   write(6,'(a)')'In such case, in the following table, the three numbers under weight should be equal,'
   write(6,'(a,/)')'and the two numbers under normalized weight should be equal.'
   write(6,'(93(''=''))')
   write(6,'(a)')'    k-vec(1)     k-vec(2)     k-vec(3)   redundancy         weight          normalized weight'
   write(6,'(93(''-''))')
   do k = 1, nk
      nw = 0
      do ir = 1, nr
         rot3d = getIBZRotationMatrix3D(ir)
         krot(1,ir) = rot3d(1,1)*kvec(1,k)+rot3d(1,2)*kvec(2,k)+rot3d(1,3)*kvec(3,k)
         krot(2,ir) = rot3d(2,1)*kvec(1,k)+rot3d(2,2)*kvec(2,k)+rot3d(2,3)*kvec(3,k)
         krot(3,ir) = rot3d(3,1)*kvec(1,k)+rot3d(3,2)*kvec(2,k)+rot3d(3,3)*kvec(3,k)
         if (abs(krot(1,ir)-kvec(1,k))+abs(krot(2,ir)-kvec(2,k))+abs(krot(3,ir)-kvec(3,k)) < TEN2m7) then
            nw = nw + 1
         endif
      enddo
!
      nd = nr   ! Number of distinguishable k-vectors among the k-vectors generated by rotating kvec
      do ir = 1, nr-1
         redundant = .false.
         LOOP_jr: do jr = ir+1, nr
            if (abs(krot(1,ir)-krot(1,jr))+abs(krot(2,ir)-krot(2,jr))+abs(krot(3,ir)-krot(3,jr)) < TEN2m7) then
               redundant = .true.
               exit LOOP_jr
            endif
         enddo LOOP_jr
         if (redundant) then
            nd = nd - 1
         endif
      enddo
!     For symmetrized calculation with IBZ rotations, the k-points generated are
!     in IBZ, and we should get nd = nr/nw = getWeight(k), as well as ONE/(nw*tw) = getWeight(k)/sumw
      write(6,'(3(f12.5,1x),2x,i5,5x,2i5,f10.5,2x,2f10.5)')kvec(1:3,k),nw,nr/nw,nd,getWeight(k),ONE/(nw*tw),getWeight(k)/sumw
   enddo
   write(6,'(93(''=''))')
!
!  -------------------------------------------------------------------
   call initStrConst(lmax_phi,NumAtoms,AtomPosition,bravais,istop,iprint)
!  -------------------------------------------------------------------
   energy = (0.6d0,0.05d0)
   kappa = sqrt(energy)
   do k = 1, nk
      kin(1)=kvec(1,k)  !+0.1
      kin(2)=kvec(2,k)  !+0.2
      kin(3)=kvec(3,k)  !+0.3
      write(6,'(/,'' K-point: kx, ky, kz ='',3f15.8)') kin(1),kin(2),kin(3)
!
!     ================================================================
!     Test 1: spherical harmonics rotation.
!     ----------------------------------------------------------------
      call calYlm(kin,lmax_phi,ylm)
!     ----------------------------------------------------------------
      do ir = 1, nr
         rot3d = getIBZRotationMatrix3D(ir)
         rotmat => getIBZRotationMatrix('n',ir)
         rotmatc => getIBZRotationMatrix('c',ir)
         kin_new = ZERO
         do j = 1, 3
            do i = 1, 3
               kin_new(i) = kin_new(i) + rot3d(i,j)*kin(j)
            enddo
         enddo
!        -------------------------------------------------------------
         call calYlm(kin_new,lmax_phi,ylm_new)
!        -------------------------------------------------------------
         ylm_rot = CZERO
         do j = 1, kmax_phi
            do i = 1, kmax_phi
               ylm_rot(i) = ylm_rot(i) + rotmat(i,j)*ylm(j)
            enddo
         enddo
         do i = 1, kmax_phi
            if (abs(ylm_rot(i)-ylm_new(i)) > tol) then
               write(6,'(/,a,i5)')'ylm_rot <> ylm_new: rotation index = ',ir
               write(6,'(a)')'3D rotation matrix:'
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
               l = lofk(i)
               write(6,'(a)')'l = ',l
               do m = -l, l
                  write(6,'(a,2i5,t30,a,2d15.8)')'l, m =',l,m,       &
                          'ylm_rot = ',ylm_rot((l+1)*(l+1)-l+m)
                  write(6,'(a,2i5,t30,a,2d15.8)')'l, m =',l,m,       &
                          'ylm_new = ',ylm_new((l+1)*(l+1)-l+m)
               enddo
               call ErrorHandler('testIBZRotation','Spherical harmonics rotation test failed!')
            endif
         enddo
      enddo
      write(6,'(/,1x,a)')'Passed spherical harmonics rotation test...'
!
!     ================================================================
!     Test 3: rotate the structure constant matrix
!     ================================================================
!     do n = 1, NumAtoms*NumAtoms
!        ia = mod(n-1,NumAtoms)+1
!        ja = (n-1)/NumAtoms+1
!        strcon => getStrConstMatrix(kin,kappa,ia,ja,lmax_phi,lmax_phi)
!        write(6,'(a,i5,a,i5)')'ia = ',ia,', ja = ',ja
!        call writeMatrix('strcon_rot',strcon,kmax_phi,kmax_phi,tol)
!     enddo
!
      LOOP_ir: do ir = 1, nr
         rot3d = getIBZRotationMatrix3D(ir)
!        write(6,'(/,a,i5)')'3D rotation matrix of ir = ',ir
!        write(6,'(3f15.8)')rot3d(1:3,1)
!        write(6,'(3f15.8)')rot3d(1:3,2)
!        write(6,'(3f15.8)')rot3d(1:3,3)
!        write(6,'(a,16i5)')'rotation_table : ',(rotation_table(ia,ir),ia = 1, NumAtoms)
!        -------------------------------------------------------------
!        call writeMatrix('rotmat',rotmat,kmax_phi,kmax_phi,tol)
!        -------------------------------------------------------------
         kin_new = ZERO
         do j = 1, 3
            do i = 1, 3
               kin_new(i) = kin_new(i) + rot3d(i,j)*kin(j)
            enddo
         enddo
         rotmatc => getIBZRotationMatrix('c',ir)
         rotmat => getIBZRotationMatrix('n',ir)
         do n = 1, NumAtoms*NumAtoms
!           ==========================================================
!           Store the structure constant matrix calculated at kin_new 
!           in strcon_new
!           ==========================================================
            ia = mod(n-1,NumAtoms)+1
            ja = (n-1)/NumAtoms+1
            strcon => getStrConstMatrix(kin_new,kappa,ia,ja,lmax_phi,lmax_phi,aij)
            strcon_new = strcon
!
!           ==========================================================
!           test structure constant matrix transformation using rotation_table and nshift
!           ==========================================================
            strcon => getStrConstMatrix(kin,kappa,rotation_table(ia,ir),rotation_table(ja,ir), &
                                        lmax_phi,lmax_phi)
            strcon_ori = strcon
            strcon_rot = CZERO
!           ----------------------------------------------------------
            call computeUAUtc(rotmatc,kmax_phi,kmax_phi,rotmatc,kmax_phi,CONE, &
                              strcon_ori,kmax_phi,CZERO,strcon_rot,kmax_phi,WORK)
!           ----------------------------------------------------------
            if (ia /= ja) then
               Ri(:) = nshift(1,ia,ir)*bravais(:,1)+nshift(2,ia,ir)*bravais(:,2)+nshift(3,ia,ir)*bravais(:,3)
               Rj(:) = nshift(1,ja,ir)*bravais(:,1)+nshift(2,ja,ir)*bravais(:,2)+nshift(3,ja,ir)*bravais(:,3)
               kr = kin(1)*(Ri(1)-Rj(1))+kin(2)*(Ri(2)-Rj(2))+kin(3)*(Ri(3)-Rj(3))
               cfac = exp(SQRTm1*kr)
               strcon_rot = cfac*strcon_rot
            endif
            failed = .false.
            LOOP_j1: do j = 1, kmax_phi
               do i = 1, kmax_phi
                  if (abs(strcon_rot(i,j)-strcon_new(i,j)) > tol) then
                     if (abs(strcon_rot(i,j)) > ONE) then
                        if (abs((strcon_rot(i,j)-strcon_new(i,j))/strcon_rot(i,j)) > tol) then
                           failed = .true.
                        endif
                     else
                        failed = .true.
                     endif
                  endif
                  if (failed) then
                     write(6,'(/,a,2i5,2x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                               i,j,strcon_rot(i,j),strcon_rot(i,j)-strcon_new(i,j)
                     exit LOOP_j1
                  endif
               enddo
            enddo LOOP_j1
            if (failed) then
               write(6,'(/,a,i5,a,i5,a,3f10.5)')'ia = ',ia,', ja = ',ja,', aij = ',aij(1:3)
               write(6,'(a,3f10.5)')'kvec after rotation = ',kin_new(1:3)
               write(6,'(a,i5)')'strcon_rot <> strcon_new: rotation index = ',ir
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
!              -------------------------------------------------------
               call writeMatrix('rotmat',rotmat,kmax_phi,kmax_phi,tol)
!              -------------------------------------------------------
               write(6,'(a,2i5)')'ia -> rotation_table : ',ia,rotation_table(ia,ir)
               write(6,'(a,2i5)')'ja -> rotation_table : ',ja,rotation_table(ja,ir)
               call writeMatrix('strcon_rot',strcon_rot,kmax_phi,kmax_phi,tol)
               call writeMatrix('strcon_new',strcon_new,kmax_phi,kmax_phi,tol)
               call ErrorHandler('testIBZRotation','Structure constant rotation test 2 failed!')
            endif
         enddo
      enddo LOOP_ir
   enddo
   write(6,'(/,1x,a)')'Passed the structure constant rotation test...'
!
!  ===================================================================
!  Test the rotation of the step function and scattering matrices
!  ===================================================================
!
   LocalNumAtoms=getLocalNumAtoms()
!
   allocate(atom_print_level(1:LocalNumAtoms), GlobalIndex(LocalNumAtoms))
   allocate(AtomicNumber(1:NumAtoms))
   allocate(lmax_pot(LocalNumAtoms), lmax_rho(LocalNumAtoms))
   allocate(lmax_kkr(LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms), lmax_green(LocalNumAtoms))
   allocate(lmax_trunc_pot(LocalNumAtoms))
!
   do ig=1,NumAtoms
      AtomicNumber(ig)=getAtomicNumber(ig)
   enddo
!
   do id=1,LocalNumAtoms
      atom_print_level(id) = getStandardOutputLevel(id)
      GlobalIndex(id)=getGlobalIndex(id)
   enddo
   lmax_max = 0
   lmax_kkr_max = 0
   lmax_rho_max = 0
   lmax_pot_max = 0
   lmax_trunc_pot_max = 0
   lmax_step_max = 0
   do id=1,LocalNumAtoms
      lmax_kkr(id) = getKKRLmax(id)
      lmax_rho(id) = getRhoLmax(id)
      lmax_pot(id) = getPotLmax(id)
      lmax_trunc_pot(id) = getTruncPotLmax(id)
      if (isFullPotential()) then
         lmax_step(id)  = max(getStepFuncLmax(id),lmax_kkr(id)*2)
         lmax_max = max( lmax_max, 2*lmax_step(id), 2*lmax_pot(id),     &
                         2*lmax_rho(id), lmax_kkr(id), lmax_trunc_pot(id) )
      else
         lmax_step(id)  = getStepFuncLmax(id)
         lmax_max = max( lmax_max, lmax_step(id), lmax_rho(id), lmax_kkr(id) )
      endif
      lmax_green(id) = lmax_rho(id)
      lmax_kkr_max = max(lmax_kkr_max,lmax_kkr(id))
      lmax_rho_max = max(lmax_rho_max,lmax_rho(id))
      lmax_pot_max = max(lmax_pot_max,lmax_pot(id))
      lmax_step_max = max(lmax_step_max,lmax_step(id))
      lmax_trunc_pot_max = max(lmax_trunc_pot_max,lmax_trunc_pot(id))
   enddo
!
!  -------------------------------------------------------------------
   call endSphericalHarmonics()
   call initSphericalHarmonics(2*lmax_max)
   call endGauntFactors()
   call initGauntFactors(lmax_max,istop,0)
!  -------------------------------------------------------------------
   call setupRadGridAndCell(LocalNumAtoms,lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Test the rotation of the step function
!  ===================================================================
   allocate(stepf((lmax_step_max+1)*(lmax_step_max+2)/2))
   allocate(stepf_rot((lmax_step_max+1)*(lmax_step_max+2)/2))
   allocate(stepfLL((lmax_kkr_max+1)**2,(lmax_kkr_max+1)**2))
   allocate(stepfLL_rot((lmax_kkr_max+1)**2,(lmax_kkr_max+1)**2))
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   do id = 1, LocalNumAtoms
      r=HALF*(getOutscrSphRadius(id)+getInscrSphRadius(id))
      jmax_step = (lmax_step(id)+1)*(lmax_step(id)+2)/2
      kmax_kkr = (lmax_kkr(id)+1)**2
      jmax_kkr = (lmax_kkr(id)+1)*(lmax_kkr(id)+2)/2
      write(6,'(/,a)')'   l   m              step function'
      do jl = 1,jmax_step
!        -------------------------------------------------------------
         stepf(jl) = getStepFunction(id,jl,r)
!        -------------------------------------------------------------
         if (abs(stepf(jl)) > tol) then
            write(6,'(2i4,5x,2d15.8)')lofj(jl),mofj(jl),stepf(jl)
         endif
      enddo
      do kl = 1, kmax_kkr
         do klp = 1, kmax_kkr
            stepfLL(klp,kl) = CZERO
            do i1=1,nj3(klp,kl)
               kl1 = kj3(i1,klp,kl)
               m1 = mofk(kl1)
               jl1 = jofk(kl1)
               if (m1 >= 0) then
                  stepfLL(klp,kl) = stepfLL(klp,kl) + cgnt(i1,klp,kl)*stepf(jl1)
               else
                  stepfLL(klp,kl) = stepfLL(klp,kl) + cgnt(i1,klp,kl)*m1m(m1)*conjg(stepf(jl1))
               endif
            enddo
         enddo
      enddo
!
      do ir = 1, nr
         rot3d = getIBZRotationMatrix3D(ir)
         rotmat => getIBZRotationMatrix('n',ir)
         stepf_rot = CZERO
         do klp = 1, kmax_kkr
            jlp = jofk(klp)
            mp = mofk(klp)
            if (mp < 0) then
               do jl = 1, jmax_kkr
                  kl = kofj(jl)
                  stepf_rot(jl) = stepf_rot(jl) + rotmat(kl,klp)*m1m(mp)*conjg(stepf(jlp))
               enddo
            else
               do jl = 1, jmax_kkr
                  kl = kofj(jl)
                  stepf_rot(jl) = stepf_rot(jl) + rotmat(kl,klp)*stepf(jlp)
               enddo
            endif
         enddo
!
         failed = .false.
!
         do jl = 1, jmax_kkr
            if (abs(stepf_rot(jl)-stepf(jl)) > tol) then
               failed = .true.
               write(6,'(/,a,i5,2x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                     jl,stepf_rot(jl),stepf_rot(jl)-stepf(jl)
               exit
            endif
         enddo
         if (failed) then
            write(6,'(a,i5)')'stepf_rot <> stepf: rotation index = ',ir
            write(6,'(3f15.8)')rot3d(1:3,1)
            write(6,'(3f15.8)')rot3d(1:3,2)
            write(6,'(3f15.8)')rot3d(1:3,3)
            write(6,'(/,a)')'   l   m           step function               step function rot'
            do jl = 1,jmax_kkr
               if (abs(stepf(jl)) > tol .or. abs(stepf_rot(jl)) > tol) then
                  write(6,'(2i4,2x,2d15.8,2x,2d15.8)')lofj(jl),mofj(jl),stepf(jl),stepf_rot(jl)
               endif
            enddo
            call ErrorHandler('testIBZRotation','Step Function rotation test failed!')
         endif
!
         rotmatc => getIBZRotationMatrix('c',ir)
         stepfLL_rot = CZERO
!        -------------------------------------------------------------
         call computeUAUtc(rotmatc,kmax_kkr,kmax_kkr,rotmatc,kmax_kkr,CONE, &
                           stepfLL,kmax_kkr,CZERO,stepfLL_rot,kmax_kkr,WORK)
!        -------------------------------------------------------------
         LOOP_kl1: do kl = 1, kmax_kkr
            do klp = 1, kmax_kkr
               if (abs(stepfLL_rot(klp,kl)-stepfLL(klp,kl)) > tol) then
                  if (abs(stepfLL_rot(klp,kl)) > ONE) then
                     if (abs((stepfLL_rot(klp,kl)-stepfLL(klp,kl))/stepfLL_rot(klp,kl)) > tol) then
                        failed = .true.
                     endif
                  else
                     failed = .true.
                  endif
               endif
               if (failed) then
                  write(6,'(/,a,2i5,2x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                            klp,kl,stepfLL_rot(klp,kl),stepfLL_rot(klp,kl)-stepfLL(klp,kl)
                  exit LOOP_kl1
               endif
            enddo
         enddo LOOP_kl1
         if (failed) then
            write(6,'(a,i5)')'stepfLL_rot <> stepfLL: rotation index = ',ir
            write(6,'(3f15.8)')rot3d(1:3,1)
            write(6,'(3f15.8)')rot3d(1:3,2)
            write(6,'(3f15.8)')rot3d(1:3,3)
!           ----------------------------------------------------------
            call writeMatrix('rotmat',rotmat,kmax_kkr,kmax_kkr,tol)
!           ----------------------------------------------------------
            call writeMatrix('stepfLL_rot',stepfLL_rot,kmax_kkr,kmax_kkr,tol)
            call writeMatrix('stepfLL',stepfLL,kmax_kkr,kmax_kkr,tol)
            call ErrorHandler('testIBZRotation','Step Function Matrix rotation test failed!')
         endif
      enddo
   enddo
   write(6,'(/,a)')'Passed the step function rotation test......'
   deallocate(stepf, stepf_rot, stepfLL, stepfLL_rot)
!
!  -------------------------------------------------------------------
   call initMadelung(LocalNumAtoms,NumAtoms,GlobalIndex,              &
                     lmax_rho_max,lmax_pot_max,bravais,AtomPosition,0)
!  -------------------------------------------------------------------
   call initSystemSymmetry( NumAtoms, LocalNumAtoms, lmax_pot, lmax_step, atom_print_level )
!  -------------------------------------------------------------------
   call calSymmetryFlags()
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot,lmax_step,               &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
   call readPotential()
!
!  ===================================================================
!  Test the rotation of the potential function
!  ===================================================================
   allocate(pot_ori((lmax_pot_max+1)*(lmax_pot_max+2)/2))
   allocate(pot_rot((lmax_pot_max+1)*(lmax_pot_max+2)/2))
   if (isFullPotential()) then
      write(6,'(/,a,i5)')'lmax_trunc_pot_max = ',lmax_trunc_pot_max
      allocate(trunc_pot_ori((lmax_trunc_pot_max+1)*(lmax_trunc_pot_max+2)/2))
      allocate(trunc_pot_rot((lmax_trunc_pot_max+1)*(lmax_trunc_pot_max+2)/2))
   endif
   do is = 1, n_spin_pola
      do id = 1, LocalNumAtoms
         Grid => getGrid(id)
         i = int((Grid%jend + Grid%jinsc)*HALF)
!
         jmax_pot = (lmax_pot(id)+1)*(lmax_pot(id)+2)/2
         pot_jl => getPotential(id,1,is)
         do jl = 1, jmax_pot
            pot_ori(jl) = pot_jl(i,jl)
         enddo
!
         do ir = 1, nr
            rot3d = getIBZRotationMatrix3D(ir)
            rotmat => getIBZRotationMatrix('n',ir)
            pot_rot = CZERO
            do klp = 1, kmax_kkr
               jlp = jofk(klp)
               mp = mofk(klp)
               if (mp < 0) then
                  do jl = 1, jmax_kkr
                     kl = kofj(jl)
                     pot_rot(jl) = pot_rot(jl) + rotmat(kl,klp)*m1m(mp)*conjg(pot_ori(jlp))
                  enddo
               else
                  do jl = 1, jmax_kkr
                     kl = kofj(jl)
                     pot_rot(jl) = pot_rot(jl) + rotmat(kl,klp)*pot_ori(jlp)
                  enddo
               endif
            enddo
            failed = .false.
            do jl = 1, jmax_kkr
               if (abs(pot_rot(jl)-pot_ori(jl)) > tol) then
                  failed = .true.
                  write(6,'(/,a,i5,2x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                        jl,pot_rot(jl),pot_rot(jl)-pot_ori(jl)
                  exit
               endif
            enddo
            if (failed) then
               write(6,'(a,i5)')'pot_rot <> pot_ori: rotation index = ',ir
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
               write(6,'(/,a)')'   l   m       potential function           potential function rot'
               do jl = 1,jmax_kkr
                  if (abs(pot_ori(jl)) > tol .or. abs(pot_rot(jl)) > tol) then
                     write(6,'(2i4,2x,2d15.8,2x,2d15.8)')lofj(jl),mofj(jl),pot_ori(jl),pot_rot(jl)
                  endif
               enddo
               call ErrorHandler('testIBZRotation','Potential function rotation test failed!')
            endif
         enddo
!
         if (isFullPotential()) then
            i = int((Grid%jend - Grid%jinsc)*HALF)
            jmax_pot = (lmax_trunc_pot(id)+1)*(lmax_trunc_pot(id)+2)/2
            write(6,'(/,a,i4)')'Testing the truncation potential rotation, lmax_pot_trunc = ', &
                               lmax_trunc_pot(id)
            pot_jl => getTruncatedPotential(id,1,is)
            do jl = 1, jmax_pot
               trunc_pot_ori(jl) = pot_jl(i,jl)
            enddo
!
            do ir = 1, nr
               rot3d = getIBZRotationMatrix3D(ir)
               rotmat => getIBZRotationMatrix('n',ir)
               trunc_pot_rot = CZERO
               do klp = 1, kmax_kkr
                  jlp = jofk(klp)
                  mp = mofk(klp)
                  if (mp < 0) then
                     do jl = 1, jmax_kkr
                        kl = kofj(jl)
                        trunc_pot_rot(jl) = trunc_pot_rot(jl) + rotmat(kl,klp)*m1m(mp)*conjg(trunc_pot_ori(jlp))
                     enddo
                  else
                     do jl = 1, jmax_kkr
                        kl = kofj(jl)
                        trunc_pot_rot(jl) = trunc_pot_rot(jl) + rotmat(kl,klp)*trunc_pot_ori(jlp)
                     enddo
                  endif
               enddo
               failed = .false.
               do jl = 1, jmax_kkr
                  if (abs(trunc_pot_rot(jl)-trunc_pot_ori(jl)) > tol) then
                     failed = .true.
                     write(6,'(/,a,i5,2x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                           jl,trunc_pot_rot(jl),trunc_pot_rot(jl)-trunc_pot_ori(jl)
                     exit
                  endif
               enddo
               if (failed) then
                  write(6,'(a,i5)')'truncated pot_rot <> truncated pot_ori: rotation index = ',ir
                  write(6,'(3f15.8)')rot3d(1:3,1)
                  write(6,'(3f15.8)')rot3d(1:3,2)
                  write(6,'(3f15.8)')rot3d(1:3,3)
                  write(6,'(/,a)')'   l   m       trunc pot function           trunc pot function rot'
                  do jl = 1,jmax_kkr
                     if (abs(trunc_pot_ori(jl)) > tol .or. abs(trunc_pot_rot(jl)) > tol) then
                        write(6,'(2i4,2x,2d15.8,2x,2d15.8)')lofj(jl),mofj(jl),trunc_pot_ori(jl),trunc_pot_rot(jl)
                     endif
                  enddo
                  call ErrorHandler('testIBZRotation','Truncated potential function rotation test failed!')
               endif
            enddo
         endif
      enddo
   enddo
   write(6,'(/,a)')'Passed the potential function rotation test......'
   deallocate(pot_ori, pot_rot)
   if (isFullPotential()) then
      deallocate(trunc_pot_ori, trunc_pot_rot)
   endif
!
!  ===================================================================
!  Test the rotation of the regular solutions and scattering matrices
!  ===================================================================
!  -------------------------------------------------------------------
   call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, &
                     lmax_kkr, lmax_kkr, lmax_pot, lmax_step, lmax_green,  &
                     n_spin_pola, n_spin_cant, 0, istop, atom_print_level)
!  -------------------------------------------------------------------
   do is = 1, n_spin_pola
      do id = 1, LocalNumAtoms
!        -------------------------------------------------------------
         call solveSingleScattering(is, id, energy, CZERO)
!        -------------------------------------------------------------
         kmax_kkr = (lmax_kkr(id)+1)**2
!
!        =============================================================
         Grid => getGrid(id)
         i = int((Grid%jend + Grid%jinsc)*HALF)
         PhiLr => getRegSolution(is,site=id,atom=1)
!
         allocate(phi_ori(kmax_kkr,kmax_kkr), phi_rot(kmax_kkr,kmax_kkr))
!
         do kl = 1, kmax_kkr
            do klp = 1, kmax_kkr
               phi_ori(klp,kl) = PhiLr(i,klp,kl)
            enddo
         enddo
         do ir = 1, nr
            rot3d = getIBZRotationMatrix3D(ir)
            rotmatc => getIBZRotationMatrix('c',ir)
            rotmat => getIBZRotationMatrix('n',ir)
            phi_rot = CZERO
!           ----------------------------------------------------------
            call computeUAUtc(rotmatc,kmax_kkr,kmax_kkr,rotmatc,kmax_kkr,CONE, &
                              phi_ori,kmax_kkr,CZERO,phi_rot,kmax_kkr,WORK)
!           ----------------------------------------------------------
            failed = .false.
            LOOP_kl2: do kl = 1, kmax_kkr
               do klp = 1, kmax_kkr
                  if (abs(phi_rot(klp,kl)-phi_ori(klp,kl)) > tol) then
                     if (abs(phi_rot(klp,kl)) > ONE) then
                        if (abs((phi_rot(klp,kl)-phi_ori(klp,kl))/phi_rot(klp,kl)) > tol) then
                           failed = .true.
                        endif
                     else
                        failed = .true.
                     endif
                  endif
                  if (failed) then
                     write(6,'(/,a,2i4,2x,2i4,4x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                               lofk(klp),mofk(klp),lofk(kl),mofk(kl),phi_rot(klp,kl),phi_ori(klp,kl)
                     exit LOOP_kl2
                  endif
               enddo
            enddo LOOP_kl2
            if (failed) then
               write(6,'(a,i5)')'phi_rot <> phi_ori: rotation index = ',ir
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
               write(6,'(a)')'phi function ::'
               write(6,'(/,a)')'   l   m    lp  mp            Regular Sol                      Reggular Sol rot'
               do kl = 1, kmax_kkr
                  do klp = 1, kmax_kkr
                     if (abs(phi_rot(klp,kl)) > tol .or. abs(phi_ori(klp,kl)) > tol) then
                        write(6,'(2i4,2x,2i4,2x,2d15.8,2x,2d15.8)')lofk(kl),mofk(kl),lofk(klp),mofk(klp), &
                                                                   phi_ori(klp,kl),phi_rot(klp,kl)
                     endif
                  enddo
               enddo
               call ErrorHandler('testIBZRotation','Regular solution rotation test failed!')
            endif
         enddo
         write(6,'(/,a)')'Passed the regular solution rotation test......'
!
         deallocate(phi_ori, phi_rot)
!        =============================================================
!
!        -------------------------------------------------------------
         allocate(sin_ori(kmax_kkr,kmax_kkr), sin_rot(kmax_kkr,kmax_kkr))
         allocate(cos_ori(kmax_kkr,kmax_kkr), cos_rot(kmax_kkr,kmax_kkr))
         allocate(jost_ori(kmax_kkr,kmax_kkr), jost_rot(kmax_kkr,kmax_kkr))
!        -------------------------------------------------------------
         jost_mat => getJostMatrix()
         sin_mat => getSineMatrix()
         cos_mat => getCosineMatrix()
!        -------------------------------------------------------------
!        call writeMatrix('jost_mat',jost_mat,kmax_kkr,kmax_kkr,TEN2m8)
         do ir = 1, nr
            rot3d = getIBZRotationMatrix3D(ir)
            rotmatc => getIBZRotationMatrix('c',ir)
            rotmat => getIBZRotationMatrix('n',ir)
            sin_ori = sin_mat
            cos_ori = cos_mat
            jost_ori = jost_mat
            sin_rot = CZERO
            cos_rot = CZERO
            jost_rot = CZERO
!           ----------------------------------------------------------
!           call computeUAUt(rotmatc,kmax_kkr,kmax_kkr,rotmat,kmax_kkr,CONE, &
!                            sin_ori,kmax_kkr,CZERO,sin_rot,kmax_kkr,WORK)
            do kl = 1, kmax_kkr
               WORK = CZERO
               do klp = 1, kmax_kkr
                  do kl1 = 1, kmax_kkr
                     WORK(kl1) = WORK(kl1) + sin_mat(kl1,klp)*rotmat(kl,klp)
                  enddo
               enddo
               do kl1 = 1, kmax_kkr
                  do klp = 1, kmax_kkr
                     sin_rot(klp,kl) = sin_rot(klp,kl) + conjg(rotmat(klp,kl1))*WORK(kl1)
                  enddo
               enddo
            enddo
!           ----------------------------------------------------------
            call computeUAUtc(rotmatc,kmax_kkr,kmax_kkr,rotmatc,kmax_kkr,CONE, &
                              cos_ori,kmax_kkr,CZERO,cos_rot,kmax_kkr,WORK)
!           ----------------------------------------------------------
            call computeUAUtc(rotmatc,kmax_kkr,kmax_kkr,rotmatc,kmax_kkr,CONE, &
                              jost_ori,kmax_kkr,CZERO,jost_rot,kmax_kkr,WORK)
!           ----------------------------------------------------------
            sin_failed = .false.
            cos_failed = .false.
            jost_failed = .false.
            LOOP_kl: do kl = 1, kmax_kkr
               do klp = 1, kmax_kkr
                  if (abs(sin_rot(klp,kl)-sin_mat(klp,kl)) > tol) then
                     if (abs(sin_rot(klp,kl)) > ONE) then
                        if (abs((sin_rot(klp,kl)-sin_mat(klp,kl))/sin_rot(klp,kl)) > tol) then
                           sin_failed = .true.
                        endif
                     else
                        sin_failed = .true.
                     endif
                  endif
                  if (sin_failed) then
                     write(6,'(/,a,2i5,2x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                               klp,kl,sin_rot(klp,kl),sin_rot(klp,kl)-sin_mat(klp,kl)
                     exit LOOP_kl
                  endif
                  if (abs(cos_rot(klp,kl)-cos_ori(klp,kl)) > tol) then
                     if (abs(cos_rot(klp,kl)) > ONE) then
                        if (abs((cos_rot(klp,kl)-cos_ori(klp,kl))/cos_rot(klp,kl)) > tol) then
                           cos_failed = .true.
                        endif
                     else
                        cos_failed = .true.
                     endif
                  endif
                  if (cos_failed) then
                     write(6,'(/,a,2i5,2x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                               klp,kl,cos_rot(klp,kl),cos_rot(klp,kl)-cos_ori(klp,kl)
                     exit LOOP_kl
                  endif
                  if (abs(jost_rot(klp,kl)-jost_ori(klp,kl)) > tol) then
                     if (abs(jost_rot(klp,kl)) > ONE) then
                        if (abs((jost_rot(klp,kl)-jost_ori(klp,kl))/jost_rot(klp,kl)) > tol) then
                           jost_failed = .true.
                        endif
                     else
                        jost_failed = .true.
                     endif
                  endif
                  if (jost_failed) then
                     write(6,'(/,a,2i5,2x,2d15.8,2x,2d15.8)')'Failed element: ',  &
                               klp,kl,jost_rot(klp,kl),jost_rot(klp,kl)-jost_ori(klp,kl)
                     exit LOOP_kl
                  endif
               enddo
            enddo LOOP_kl
!           if (sin_failed .or. ir == 5) then
            if (sin_failed) then
               write(6,'(a,i5)')'sin_rot <> sin_ori: rotation index = ',ir
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
!              write(6,'(a)')'rotmat ::'
               write(6,'(a)')'sin matrix ::'
               write(6,'(/,a)')'   l   m    lp  mp            sine matrix                  sine matrix rot'
               do kl = 1, kmax_kkr
                  do klp = 1, kmax_kkr
!                    if (abs(rotmat(kl,klp)) > tol) then
!                       write(6,'(2i4,4x,2i4,2x,2d15.8)')lofk(kl),mofk(kl),lofk(klp),mofk(klp), &
!                                                        rotmat(kl,klp)
!                    endif
                     if (abs(sin_rot(klp,kl)) > tol .or. abs(sin_ori(klp,kl)) > tol) then
                        write(6,'(2i4,2x,2i4,2x,2d15.8,2x,2d15.8)')lofk(kl),mofk(kl),lofk(klp),mofk(klp), &
                                                                   sin_ori(klp,kl),sin_rot(klp,kl)
                     endif
                  enddo
               enddo
!              -------------------------------------------------------
!              call writeMatrix('rotmat',rotmat,kmax_kkr,kmax_kkr,tol)
!              -------------------------------------------------------
!              call writeMatrix('sin_rot',sin_rot,kmax_kkr,kmax_kkr,tol)
!              call writeMatrix('sin_mat',sin_mat,kmax_kkr,kmax_kkr,tol)
               call ErrorHandler('testIBZRotation','Sine matrix rotation test failed!')
            else if (cos_failed) then
               write(6,'(a,i5)')'cos_rot <> cos_ori: rotation index = ',ir
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
!              -------------------------------------------------------
               call writeMatrix('rotmat',rotmat,kmax_kkr,kmax_kkr,tol)
!              -------------------------------------------------------
               call writeMatrix('cos_rot',cos_rot,kmax_kkr,kmax_kkr,tol)
               call writeMatrix('cos_ori',cos_ori,kmax_kkr,kmax_kkr,tol)
               call ErrorHandler('testIBZRotation','Cosine matrix rotation test failed!')
            else if (jost_failed) then
               write(6,'(a,i5)')'jost_rot <> jost_ori: rotation index = ',ir
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
!              -------------------------------------------------------
               call writeMatrix('rotmat',rotmat,kmax_kkr,kmax_kkr,tol)
!              -------------------------------------------------------
               call writeMatrix('jost_rot',jost_rot,kmax_kkr,kmax_kkr,tol)
               call writeMatrix('jost_ori',jost_ori,kmax_kkr,kmax_kkr,tol)
               call ErrorHandler('testIBZRotation','Jost matrix rotation test failed!')
            endif
         enddo
         deallocate(sin_ori, sin_rot, cos_ori, cos_rot, jost_ori, jost_rot)
      enddo
   enddo
   write(6,'(/,a)')'Passed the scattering matrix rotation test......'
!
!  ===================================================================
   deallocate( atom_print_level, GlobalIndex, AtomicNumber, lmax_pot, lmax_rho, &
               lmax_kkr, lmax_step, lmax_green )
   deallocate(AtomPosition, AtomName, nshift)
   deallocate(strcon_new)
   deallocate(ylm, ylm_new, ylm_rot, emat, strcon_rot, WORK, strcon_ori, us)
!
   call endSSSolver()
   call endPotential()
   call endMadelung()
   call endSystemSymmetry()
   call endStrConst()
   call endIBZRotation()
   call endGauntFactors()
   call endSphericalHarmonics()
   call endLattice()
   call finishProcess()
!
end program testIBZRotation
