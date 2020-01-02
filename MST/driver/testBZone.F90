program testBZone
!  ********************************************************************
!  main to test structure constant code.
!  ********************************************************************
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, &
                                       isDataStorageExisting, &
                                       getDataStorage, RealMark
   use InputModule, only : initInput, endInput, getKeyValue
!
   use SystemModule, only : initSystem, endSystem
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber
!
   use ScfDataModule, only : initScfData, endScfData
   use ScfDataModule, only : getPotentialTypeParam
   use ScfDataModule, only : isReadKmesh, getKmeshFileName
   use ScfDataModule, only : NumKMeshs, kGenScheme, Kdiv, Symmetrize
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType, &
                                   isFullPotential
!
   use BZoneModule, only : initBZone, printBZone, endBZone, getNumKs,   &
                           getAllKPoints, getBZoneLattice, getRotMatrix3D, getNumRotations
   use BZoneModule, only : getAllWeights, getWeightSum
!
   use TimerModule, only : initTimer, getTime
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, SQRTm1
!
   use MPPModule, only : initMPP, endMPP
!
   use Matrix3dModule, only : invm3
!
   implicit   none
!
   character (len=4) :: istop = 'none'
!
   logical :: redundant
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: def_id, info_id, ierr
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: nk = 1
!
   integer (kind=IntKind) :: i, j, n, ir, jr, nr, nw, nd, k, mk
   integer (kind=IntKind), parameter :: MaxRotations = 48
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
!
   real (kind=RealKind), parameter :: tol = TEN2m6
   real (kind=RealKind) :: t0, kfac, tw
   real (kind=RealKind), pointer :: Bravais(:,:)
   real (kind=RealKind), pointer :: kvec(:,:), bz_latt(:,:)
   real (kind=RealKind) :: binv(3,3), rot(3,3), k1, k2, k3, krot(3,MaxRotations)
   real (kind=RealKind), allocatable :: kvec_BZ(:,:)
   real (kind=RealKind), pointer :: wght(:)
   real (kind=RealKind) :: wsum
!
!  *******************************************************************
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
!
!  -------------------------------------------------------------------
   call initMPP()
   call initDataServiceCenter()
   call initInput()
   call readInputs(def_id,info_id)
   call initScfData(def_id)
   call initPotentialType(getPotentialTypeParam())
   call initSystem(def_id)
!  -------------------------------------------------------------------
!  
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      Bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('testProcMapping','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!
   NumAtoms = getNumAtoms()
!
   allocate(AtomPosition(1:3,1:NumAtoms), AtomicNumber(NumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
      AtomicNumber(i)=getAtomicNumber(i)
   enddo
!
!  ===================================================================
!  Initialize the Brillouin zone mesh for k-space integration
!  ===================================================================
   if (isReadKmesh()) then
!     ----------------------------------------------------------------
      call initBZone(getKmeshFileName(),'none',-1)
!     ----------------------------------------------------------------
   else if (NumKMeshs > 0) then
!     ----------------------------------------------------------------
      call initBZone(NumKMeshs,kGenScheme,Kdiv,Symmetrize,Bravais, &
                     NumAtoms,AtomPosition,AtomicNumber,'none',-1)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('main','No K mesh is initialized')
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call printBZone()
   nk = getNumKs()
   kvec => getAllKPoints(kfac)
!  -------------------------------------------------------------------
!
   if (.not.isReadKmesh()) then
!     ================================================================
!     Print (k1,k2,k3), which are the k-points in BZ lattice coordinates:
!           kvec = k1*BZ(:,1) + k2*BZ(:,2) + k3*BZ(:,3)
!     ================================================================
      bz_latt => getBZoneLattice()
      write(6,'(/,a)')'Brillioun zone lattice:'
      write(6,'(3f12.5)')bz_latt(1:3,1)
      write(6,'(3f12.5)')bz_latt(1:3,2)
      write(6,'(3f12.5)')bz_latt(1:3,3)
      nr = getNumRotations()
      if (nr > MaxRotations) then
         call ErrorHandler('testBZone','Number of crystal rotations exceeds its physical limit',nr,MaxRotations)
      endif
!
      call invm3(bz_latt,binv,ierr)
      if (ierr > 0) then
         call ErrorHandler('testBZone','Error in inversing the bravais lattice matrix')
      endif
!
!     ================================================================
!     If the k-points are generated only within IBZ, produce
!     k-points in the entire BZ by rotation
!     ================================================================
      if (Symmetrize > 0) then
         allocate(kvec_BZ(3,nr*nk))
         write(6,'(/,a)')'====== k-points in entire BZ ==============='
         write(6,'(a)')  '     Kx          Ky          Kz'
         write(6,'(a)')  '--------------------------------------------'
         mk = 0
         do k = 1, nk
            mk = mk + 1
            kvec_BZ(1,mk) = kvec(1,k)
            kvec_BZ(2,mk) = kvec(2,k)
            kvec_BZ(3,mk) = kvec(3,k)
            do ir = 1, nr
               rot = getRotMatrix3D(ir)
               k1 = rot(1,1)*kvec(1,k)+rot(1,2)*kvec(2,k)+rot(1,3)*kvec(3,k)
               k2 = rot(2,1)*kvec(1,k)+rot(2,2)*kvec(2,k)+rot(2,3)*kvec(3,k)
               k3 = rot(3,1)*kvec(1,k)+rot(3,2)*kvec(2,k)+rot(3,3)*kvec(3,k)
               redundant = .false.
               LOOP_j: do j = 1, mk
                  if (abs(k1-kvec_BZ(1,j))+abs(k2-kvec_BZ(2,j))+abs(k3-kvec_BZ(3,j)) < tol) then
                     redundant = .true.
                     exit LOOP_j
                  endif
               enddo LOOP_j
               if (.not.redundant) then
                  mk = mk + 1
                  kvec_BZ(1,mk) = k1
                  kvec_BZ(2,mk) = k2
                  kvec_BZ(3,mk) = k3
               endif
            enddo
         enddo
         do j = 1, mk
            write(6,'(3f12.5)')kvec_BZ(1:3,j)
         enddo
         write(6,'(a)')'=============================================='
         deallocate(kvec_BZ)
      endif
!
      tw = ZERO ! normalization factor for the weight
      do k = 1, nk
         nw = 0
         do ir = 1, nr
            rot = getRotMatrix3D(ir)
            krot(1,ir) = rot(1,1)*kvec(1,k)+rot(1,2)*kvec(2,k)+rot(1,3)*kvec(3,k)
            krot(2,ir) = rot(2,1)*kvec(1,k)+rot(2,2)*kvec(2,k)+rot(2,3)*kvec(3,k)
            krot(3,ir) = rot(3,1)*kvec(1,k)+rot(3,2)*kvec(2,k)+rot(3,3)*kvec(3,k)
            if (abs(krot(1,ir)-kvec(1,k))+abs(krot(2,ir)-kvec(2,k))+abs(krot(3,ir)-kvec(3,k)) < TEN2m6) then
               nw = nw + 1
            endif
         enddo
         tw = tw + ONE/real(nw,kind=RealKind)
      enddo
!
      write(6,'(/,a)')'Transforming the k-points coordinates...'
      write(6,'(a)')'The k-points in BZ coordinates v.s. k-points in Cartition coordinates, redundancy, weight, and fraction'
!
      wght => getAllWeights()
      wsum = getWeightSum()
!
      do k = 1, nk
         k1 = binv(1,1)*kvec(1,k)+binv(1,2)*kvec(2,k)+binv(1,3)*kvec(3,k)
         k2 = binv(2,1)*kvec(1,k)+binv(2,2)*kvec(2,k)+binv(2,3)*kvec(3,k)
         k3 = binv(3,1)*kvec(1,k)+binv(3,2)*kvec(2,k)+binv(3,3)*kvec(3,k)
!
         nw = 0   ! Number of rotations that rotate kvec back to itself...
         do ir = 1, nr
            rot = getRotMatrix3D(ir)
            krot(1,ir) = rot(1,1)*kvec(1,k)+rot(1,2)*kvec(2,k)+rot(1,3)*kvec(3,k)
            krot(2,ir) = rot(2,1)*kvec(1,k)+rot(2,2)*kvec(2,k)+rot(2,3)*kvec(3,k)
            krot(3,ir) = rot(3,1)*kvec(1,k)+rot(3,2)*kvec(2,k)+rot(3,3)*kvec(3,k)
            if (abs(krot(1,ir)-kvec(1,k))+abs(krot(2,ir)-kvec(2,k))+abs(krot(3,ir)-kvec(3,k)) < TEN2m6) then
               nw = nw + 1
            endif
         enddo
!
!        nw = the redundancy
!        weight = nr/nw, which is one way to calculate the weight.
!
         nd = nr   ! Number of distinguishable k-vectors among the k-vectors generated by rotating kvec 
         do ir = 1, nr-1
            redundant = .false.
            LOOP_jr: do jr = ir+1, nr
               if (abs(krot(1,ir)-krot(1,jr))+abs(krot(2,ir)-krot(2,jr))+abs(krot(3,ir)-krot(3,jr)) < TEN2m6) then
                  redundant = .true.
                  exit LOOP_jr
               endif
            enddo LOOP_jr
            if (redundant) then
               nd = nd - 1
            endif
         enddo
!        We should get nd = nr/nw...........................
         write(6,'(3f12.5,a,3f12.5,2x,i5,2x,2i5,2x,3f10.5)')k1,k2,k3,'   ==',kvec(1:3,k),nw,nr/nw,nd,ONE/(nw*tw), &
                                                                             wght(k),wght(k)/wsum
      enddo
   endif
!
   call endBZone()
   call endSystem()
   call endPotentialType()
   call endScfData()
   call endInput()
   call endDataServiceCenter()
   call endMPP()
!
end program testBZone
