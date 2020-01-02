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
                            getAtomName
!
   use ScfDataModule, only : NumKMeshs, Symmetrize, isLSMS
!
   use TimerModule, only : initTimer, getTime
!
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1,        &
                               TEN2m7, TEN2m8, TEN2m9, TEN2m10, TEN2m11, TEN2m12
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
   use SphericalHarmonicsModule, only : calylm
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
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
!
   use LatticeModule, only : initLattice, endLattice, getLatticeType
! 
   use IntegerFactorsModule, only : lofk
!
   use StrConstModule, only : initStrConst, endStrConst
   use StrConstModule, only : getStrConstMatrix
!
   use MatrixModule, only : computeUAUt, computeUAUtc, computeUAU, computeUAUts
!
   implicit   none
!
   character (len=4) :: istop = 'none'
   character (len=MaxLenOfAtomName), allocatable :: AtomName(:)
!
   logical :: redundant, failed
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: lmax_phi
   integer (kind=IntKind) :: kmax_phi
!
   integer (kind=IntKind) :: i, j, l, m, n, ij
   integer (kind=IntKind) :: k, nk, nr, ir, jr, ia, ja, nw, nd
   integer (kind=IntKind), parameter :: MaxRotations = 48
   integer (kind=IntKind), pointer :: rotation_table(:,:)
   integer (kind=IntKind), allocatable :: nshift(:,:,:)
!
   real (kind=RealKind) :: t0, kfac, kr, tw, sumw
   real (kind=RealKind) :: rot3d(3,3), bravais(3,3), Ri(3), Rj(3)
   real (kind=RealKind) :: kin(3), kin_new(3), aij(3), krot(3,MaxRotations)
   real (kind=RealKind), pointer :: kvec(:,:)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), parameter :: tol = TEN2m10
!
   complex (kind=CmplxKind) :: energy, kappa, cfac
   complex (kind=CmplxKind), pointer :: rotmat(:,:), rotmatc(:,:)
   complex (kind=CmplxKind), pointer :: strcon(:,:)
   complex (kind=CmplxKind), allocatable :: strcon_new(:,:)
   complex (kind=CmplxKind), allocatable :: strcon_rot(:,:), strcon_ori(:,:)
   complex (kind=CmplxKind), allocatable :: ylm(:), ylm_new(:), ylm_rot(:)
   complex (kind=CmplxKind), allocatable :: emat(:,:)
   complex (kind=CmplxKind), allocatable :: WORK(:), us(:)
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
   do i=1,NumAtoms
      AtomName(i) = getAtomName(i)
      AtomPosition(1:3,i)=getAtomPosition(i)
      lmax_phi = max(lmax_phi,getPhiLmax(i))
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
               write(6,'(a,i5)')'ylm_rot <> ylm_new: rotation index = ',ir
               write(6,'(/,a)')'3D rotation matrix:'
               write(6,'(3f15.8)')rot3d(1:3,1)
               write(6,'(3f15.8)')rot3d(1:3,2)
               write(6,'(3f15.8)')rot3d(1:3,3)
               l = lofk(i)
               write(6,'(/,a)')'l = ',l
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
      write(6,'(/,a)')'Succeeded in spherical harmonics rotation test...'
!
!     ================================================================
!     Test 2: rotmat is unitary matrix
!     ================================================================
      do ir = 1, nr
         rot3d = getIBZRotationMatrix3D(ir)
         rotmat => getIBZRotationMatrix('n',ir)
         rotmatc => getIBZRotationMatrix('c',ir)
!        -------------------------------------------------------------
         call zgemm('n','t',kmax_phi,kmax_phi,kmax_phi,CONE,rotmat,   &
                    kmax_phi,rotmatc,kmax_phi,CZERO,emat,kmax_phi)
!        -------------------------------------------------------------
         do j = 1, kmax_phi
            do i = 1, kmax_phi
               if ((i == j .and. abs(emat(i,j)-CONE) > tol) .or.      &
                   (i /= j .and. abs(emat(i,j)) > tol)) then
                  write(6,'(a,i5)')'Rotmat is not unitary: rotation index = ',ir
                  write(6,'(/,a)')'3D rotation matrix:'
                  write(6,'(3f15.8)')rot3d(1:3,1)
                  write(6,'(3f15.8)')rot3d(1:3,2)
                  write(6,'(3f15.8)')rot3d(1:3,3)
                  call writeMatrix('Rotation matrix',rotmat,kmax_phi,kmax_phi,tol)
                  call ErrorHandler('testIBZRotation','Unitary matrix test failed!')
               endif
            enddo
         enddo
      enddo
      write(6,'(/,a)')'Succeeded in unitary matrix test...'
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
                     failed = .true.
                     exit LOOP_j1
                  endif
               enddo
            enddo LOOP_j1
            if (failed) then
               write(6,'(a,i5,a,i5,a,3f10.5)')'ia = ',ia,', ja = ',ja,', aij = ',aij(1:3)
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
      write(6,'(/,a)')'Succeeded in structure constant rotation test...'
   enddo
!
   deallocate(AtomPosition, AtomName, nshift)
   deallocate(strcon_new)
   deallocate(ylm, ylm_new, ylm_rot, emat, strcon_rot, WORK, strcon_ori, us)
!
   call endStrConst()
   call endIBZRotation()
   call endGauntFactors()
   call endSphericalHarmonics()
   call endLattice()
   call finishProcess()
!
end program testIBZRotation
