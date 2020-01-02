module BZoneModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
   use MathParamModule, only : ZERO, ONE, TEN2m6, HALF, THREE, CZERO, PI2
   use MPPModule, only : MyPE
!
public :: initBZone,         &
          endBZone,          &
!         getNumKMeshs,      &
          getNumKs,          &
          getKPoint,         &
          getNumRotations,   &
          getRotMatrix3D,    &
          getWeight,         &
          getWeightSum,      &
          getAllKPoints,     &
          getAllWeights,     &
          getBZoneVolume,    &
          getBZoneLattice,   &
          integBZone,        &
          printBZone
!
   interface initBZone
      module procedure initBZone0, initBZone1
   end interface
!
   interface getKPoint
      module procedure getKPoint0, getKPoint1
   end interface
!
   interface getAllKPoints
      module procedure getAllKPoints0, getAllKPoints1
   end interface
!
   interface getNumKs
      module procedure getNumKs0, getNumKs1
   end interface
!
   interface getWeight
      module procedure getWeight0, getWeight1
   end interface
!
   interface getAllWeights
      module procedure getAllWeights0, getAllWeights1
   end interface
!
   interface getWeightSum
      module procedure getWeightSum0, getWeightSum1
   end interface
!
   interface integBZone
      module procedure integBZoneR0, integBZoneR1, integBZoneC0, integBZoneC1
   end interface
!
private
!
   logical :: Initialized = .false.
   logical :: lattice_defined = .false.
!
   integer (kind=IntKind), allocatable :: NumKs(:)
   integer (kind=IntKind) :: NumMeshs
   integer (kind=IntKind) :: kGenScheme
   integer (kind=IntKind) :: NumRotations
!
   integer (kind=IntKind), parameter :: Symmetrize = 1
   integer (kind=IntKind), parameter :: SpecialKpoints = 0
   integer (kind=IntKind), parameter :: Tetrahedron = 1
   integer (kind=IntKind), parameter :: Direction = 2
   integer (kind=IntKind), parameter :: GammaOnly = -1
!
   real (kind=RealKind), allocatable, target :: KPoints(:,:,:)
   real (kind=RealKind), pointer :: KxPoint(:,:)
   real (kind=RealKind), pointer :: KyPoint(:,:)
   real (kind=RealKind), pointer :: KzPoint(:,:)
   real (kind=RealKind), allocatable, target :: Kweight(:,:)
   real (kind=RealKind), allocatable :: KweightSum(:)
!
   integer (kind=IntKind) :: MessageID = 20000
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind) :: initlz = 0
!
   real (kind=RealKind) :: vol_bz, vol_ws
   real (kind=RealKind) :: Rotation(49,3,3)
   real (kind=RealKind) :: kfac = ZERO
   real (kind=RealKind), target :: BZ_latt(3,3)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initBZone0(kmesh_file,istop,iprint)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: kmesh_file
   character (len=*), intent(in) :: istop
   character (len=1) :: dummy
!
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: ios, im, ik, maxk
!
   stop_routine = istop
   print_level = iprint
!
   initlz = 0
!
   open(unit=12,file=kmesh_file,form='formatted', status='old',iostat=ios)
   if (ios /= 0) then
      call ErrorHandler('initBZone0','invalid file name',kmesh_file)
   endif
   read(12,*)NumMeshs
   allocate(NumKs(NumMeshs))
   maxk=0
   do im=1,NumMeshs
      read(12,*)NumKs(im)
      maxk=max(maxk,NumKs(im))
      do ik = 1, NumKs(im)
         read(12,'(a)')dummy
      enddo
   enddo
   if (maxk == 0) then
      call ErrorHandler('initBZone0','maxk = 0')
   endif
   maxk = maxk*48
   rewind(12)
!
   allocate( KPoints(1:3,1:maxk,1:NumMeshs) )
   KxPoint => KPoints(1,1:maxk,1:NumMeshs)
   KyPoint => KPoints(2,1:maxk,1:NumMeshs)
   KzPoint => KPoints(3,1:maxk,1:NumMeshs)
   allocate(Kweight(maxk,NumMeshs), KweightSum(NumMeshs))
!
   read(12,*)NumMeshs, kfac
   do im=1,NumMeshs
      read(12,*)NumKs(im)
      KweightSum(im) = ZERO
      do ik = 1, NumKs(im)
         read(12,*)KxPoint(ik,im),KyPoint(ik,im),KzPoint(ik,im),      &
                   Kweight(ik,im)
         KweightSum(im) = KweightSum(im) + Kweight(ik,im)
      enddo
   enddo
!
   NumRotations=1
   Rotation(1,1:3,1:3)=ZERO
   Rotation(1,1,1)=ONE
   Rotation(1,2,2)=ONE
   Rotation(1,3,3)=ONE
!
   close(12)
   Initialized = .true.
   lattice_defined = .false.
!
   end subroutine initBZone0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initBZone1(num_meshs,kgen,Kdiv,isym,bravais,           &
                          num_atoms,avec,AtomType,istop,iprint)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: istop
   character (len=11), parameter :: sname = 'initKspace1'
!
   integer (kind=IntKind), intent(in) :: num_meshs
   integer (kind=IntKind), intent(in) :: kgen
   integer (kind=IntKind), intent(in) :: isym
   integer (kind=IntKind), intent(in) :: Kdiv(3,num_meshs)
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: AtomType(num_atoms)
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: im, ik, i, j, ikpsym
   integer (kind=IntKind) :: maxk
   integer (kind=IntKind), allocatable :: itype(:)
   integer (kind=IntKind), allocatable :: wght(:,:)
!
   real (kind=RealKind), intent(in) :: bravais(3,3)
   real (kind=RealKind), intent(in) :: avec(3,num_atoms)
   real (kind=RealKind) :: vws, vbz
   real (kind=RealKind), allocatable :: ax(:), ay(:), az(:), qmesh(:,:,:)
!
   if (num_meshs < 1) then
      call ErrorHandler(sname,'NumMeshs < 1',num_meshs)
   else if (num_atoms < 1) then
      call ErrorHandler(sname,'NumAtoms < 1',num_atoms)
   endif
!
   initlz = 0
   NumMeshs = num_meshs
   allocate(NumKs(NumMeshs))
!maxk
   kGenScheme = kgen
   stop_routine = istop
   print_level = iprint
!
   if (isym >= 0 .and. isym/=1) then
!     ikpsym=-2    ! This change is made on July 23, 2016
      ikpsym=-1    ! by Yang Wang
   else
      ikpsym=isym
   endif
!
   do im=1,NumMeshs
      if (Kdiv(1,im)==0) then
         call ErrorHandler(sname,'Kdiv(1,im) = 0')
      else if (Kdiv(2,im)==0) then
         call ErrorHandler(sname,'Kdiv(2,im) = 0')
      else if (Kdiv(3,im)==0) then
         call ErrorHandler(sname,'Kdiv(3,im) = 0')
      endif
   enddo
!
   allocate(itype(num_atoms), ax(num_atoms), ay(num_atoms), az(num_atoms))
   do i=1,num_atoms
      ax(i)=avec(1,i)
      ay(i)=avec(2,i)
      az(i)=avec(3,i)
      itype(i)=AtomType(i)
   enddo
!
   if (kGenScheme == GammaOnly) then
      NumKs = 1
      allocate( KPoints(3,1,1:NumMeshs) )
      KxPoint => KPoints(1,1:1,1:NumMeshs)
      KyPoint => KPoints(2,1:1,1:NumMeshs)
      KzPoint => KPoints(3,1:1,1:NumMeshs)
      allocate(Kweight(1,NumMeshs), KweightSum(NumMeshs))
      do im=1,NumMeshs
         KweightSum(im) = ZERO
         do ik=1,NumKs(im)
            KxPoint(ik,im)=ZERO
            KyPoint(ik,im)=ZERO
            KzPoint(ik,im)=ZERO
            Kweight(ik,im)=ONE
            KweightSum(im) = KweightSum(im) + Kweight(ik,im)
         enddo
      enddo
      NumRotations=1
      Rotation(1,1:3,1:3)=ZERO
      Rotation(1,1,1)=ONE
      Rotation(1,2,2)=ONE
      Rotation(1,3,3)=ONE
   else if (kGenScheme == Tetrahedron) then
      call ErrorHandler(sname,'Not implemented scheme',kGenScheme)
   else if (kGenScheme == Direction) then
      call ErrorHandler(sname,'Not implemented scheme',kGenScheme)
   else if (kGenScheme == SpecialKpoints) then
!     ================================================================
!     Calculate volume of Unit cell volume............................
!     ================================================================
      vws=(bravais(2,1)*bravais(3,2)-bravais(3,1)*bravais(2,2))*bravais(1,3) &
         +(bravais(3,1)*bravais(1,2)-bravais(1,1)*bravais(3,2))*bravais(2,3) &
         +(bravais(1,1)*bravais(2,2)-bravais(2,1)*bravais(1,2))*bravais(3,3)
!     ----------------------------------------------------------------
      call vectprd(bravais(1,2),bravais(1,3),BZ_latt(1,1))
!     ----------------------------------------------------------------
      call vectprd(bravais(1,3),bravais(1,1),BZ_latt(1,2))
!     ----------------------------------------------------------------
      call vectprd(bravais(1,1),bravais(1,2),BZ_latt(1,3))
!     ----------------------------------------------------------------
      do j = 1,3
         do i = 1,3
            BZ_latt(i,j) = BZ_latt(i,j)/vws
         enddo
      enddo
      kfac = PI2
!     ================================================================
!     Calculate volume of Brillouin zone..............................
!     ================================================================
      vbz = (BZ_latt(2,1)*BZ_latt(3,2)-BZ_latt(3,1)*BZ_latt(2,2))*BZ_latt(1,3) &
           +(BZ_latt(3,1)*BZ_latt(1,2)-BZ_latt(1,1)*BZ_latt(3,2))*BZ_latt(2,3) &
           +(BZ_latt(1,1)*BZ_latt(2,2)-BZ_latt(2,1)*BZ_latt(1,2))*BZ_latt(3,3)
      vol_ws=abs(vws)
      vol_bz=abs(vbz)*kfac**3
!
!     ================================================================
!     allocate temporary space for k-meshs and their weight...........
!     ================================================================
      maxk=0
      do im=1,NumMeshs
         maxk=max(maxk,Kdiv(1,im)*Kdiv(2,im)*Kdiv(3,im))
      enddo
      if (maxk == 0) then
         call ErrorHandler(sname,'maxk = 0')
      endif
      maxk = maxk*48
      allocate(qmesh(3,maxk,NumMeshs),wght(maxk,NumMeshs))
!
!     ================================================================
!     Generate k-point meshs in Brillouin zone........................
!     ----------------------------------------------------------------
      call genkpt1(bravais, BZ_latt,                                   &
                   ax, ay, az,                                         &
                   itype, num_atoms,                                   &
                   ikpsym,                                             &
                   maxk,                                               &
                   Kdiv,                                               &
                   qmesh, wght)
!     ----------------------------------------------------------------
      maxk=0
      do im=1,NumMeshs
         maxk=max(maxk,NumKs(im))
      enddo
      allocate( KPoints(1:3,1:maxk,1:NumMeshs) )
      KxPoint => KPoints(1,1:maxk,1:NumMeshs)
      KyPoint => KPoints(2,1:maxk,1:NumMeshs)
      KzPoint => KPoints(3,1:maxk,1:NumMeshs)
      allocate(Kweight(maxk,NumMeshs), KweightSum(NumMeshs))
      do im=1,NumMeshs
         KweightSum(im) = ZERO
         do ik=1,NumKs(im)
            KxPoint(ik,im)=qmesh(1,ik,im)
            KyPoint(ik,im)=qmesh(2,ik,im)
            KzPoint(ik,im)=qmesh(3,ik,im)
            Kweight(ik,im)=wght(ik,im)
            KweightSum(im) = KweightSum(im) + wght(ik,im)
         enddo
      enddo
      deallocate(qmesh, wght)
   endif
   deallocate(itype,ax, ay, az)
!
   Initialized = .true.
   lattice_defined = .true.
!
   end subroutine initBZone1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endBZone()
!  ===================================================================
!
   deallocate(NumKs)
   deallocate(Kweight, KweightSum, KPoints)
   nullify( KxPoint, KyPoint, KzPoint )
!
   Initialized = .false.
!
   end subroutine endBZone
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKPoint0(ik,kf) result(kp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ik
   real (kind=RealKind), intent(out) :: kf
   real (kind=RealKind) :: kp(3)
!
   if (.not.Initialized) then
      call ErrorHandler('getKPoint','BZoneModule is not initialized')
   else if (ik>NumKs(1) .or. ik < 1) then
      call ErrorHandler('getKPoint','k point index out of range',ik)
   endif
   kp(1) = KxPoint(ik,1)
   kp(2) = KyPoint(ik,1)
   kp(3) = KzPoint(ik,1)
   kf = kfac
!
   end function getKPoint0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKPoint1(ik,im,kf) result(kp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ik, im
   real (kind=RealKind), intent(out) :: kf
   real (kind=RealKind) :: kp(3)
!
   if (.not.Initialized) then
      call ErrorHandler('getKPoint','BZoneModule is not initialized')
   else if (im>NumMeshs .or. im < 1) then
      call ErrorHandler('getKPoint','mesh index out of range',im)
   else if (ik>NumKs(im) .or. ik < 1) then
      call ErrorHandler('getKPoint','k point index out of range',ik)
   endif
   kp(1) = KxPoint(ik,im)
   kp(2) = KyPoint(ik,im)
   kp(3) = KzPoint(ik,im)
   kf = kfac
!
   end function getKPoint1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAllKPoints0(kf) result(kp)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: kp(:,:)
   real (kind=RealKind), intent(out) :: kf
!
   if (.not.Initialized) then
      call ErrorHandler('getAllKPoints','BZoneModule is not initialized')
   endif
   kp => KPoints(1:3,1:NumKs(1),1)
   kf = kfac
!
   end function getAllKPoints0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAllKPoints1(im,kf) result(kp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: im
   real (kind=RealKind), pointer :: kp(:,:)
   real (kind=RealKind), intent(out) :: kf
!
   if (.not.Initialized) then
      call ErrorHandler('getAllKPoints','BZoneModule is not initialized')
   else if (im>NumMeshs .or. im < 1) then
      call ErrorHandler('getAllKPoints','mesh index out of range',im)
   endif
   kp => KPoints(1:3,1:NumKs(im),im)
   kf = kfac
!
   end function getAllKPoints1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumKs0() result(nk)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: nk
!
   if (.not.Initialized) then
      call ErrorHandler('getNumKs','BZoneModule is not initialized')
   endif
   nk = NumKs(1)
!
   end function getNumKs0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumKs1(im) result(nk)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: im
!
   integer (kind=IntKind) :: nk
!
   if (.not.Initialized) then
      call ErrorHandler('getNumKs','BZoneModule is not initialized')
   else if (im>NumMeshs .or. im < 1) then
      call ErrorHandler('getKPoint','mesh index out of range',im)
   endif
   nk = NumKs(im)
!
   end function getNumKs1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumRotations() result(nr)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: nr
!
   if (.not.Initialized) then
      call ErrorHandler('getNumRotations','BZoneModule is not initialized')
   endif
   nr = NumRotations
!
   end function getNumRotations
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRotMatrix3D(ir) result(rm)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: ir
   real (kind=RealKind) :: rm(3,3)
!
   if (.not.Initialized) then
      call ErrorHandler('getNumRotMatrix','BZoneModule is not initialized')
   else if (ir < 1 .or. ir > NumRotations) then
      call ErrorHandler('getRotMatrix','rotation matrix index out of range',ir)
   endif
   rm(1:3,1:3)=Rotation(ir,1:3,1:3)
!
   end function getRotMatrix3D
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBZoneLattice() result(latt)
!  ===================================================================
   implicit none
   real (kind=RealKind), pointer :: latt(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getBZoneLattice','BZoneModule is not initialized')
   else if (.not.lattice_defined) then
      call ErrorHandler('getBZoneLattice','BZ lattice is not defined in initialization')
   endif
   latt => BZ_latt
!
   end function getBZoneLattice
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine vectprd(a,b,c)
!  ====================================================================
!         c = [a,b]
!  ====================================================================
   implicit none
!  ====================================================================
   real (kind=RealKind), intent(in) :: a(3)
   real (kind=RealKind), intent(in) :: b(3)
   real (kind=RealKind), intent(out) :: c(3)
!  ====================================================================
   c(1) = a(2)*b(3) - b(2)*a(3)
   c(2) = a(3)*b(1) - b(3)*a(1)
   c(3) = a(1)*b(2) - b(1)*a(2)
!  ====================================================================
   end subroutine vectprd
!  ====================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine genkpt1(bravais, k_latt,                               &
                       ax, ay, az,                                    &
                       itype, num_atoms,                              &
                       ikpsym,                                        &
                       maxk,                                          &
                       Kdiv,                                          &
                       qmesh, wght)
!  ====================================================================
!
   implicit none
!
   character (len=7), parameter :: sname='genkpt1'
!
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(inout) :: itype(num_atoms)
   integer (kind=IntKind), intent(in) :: ikpsym
   integer (kind=IntKind), intent(in) :: maxk
   integer (kind=IntKind), intent(in) :: Kdiv(3,NumMeshs)
   integer (kind=IntKind), intent(out) :: wght(maxk,NumMeshs)
!
   integer (kind=IntKind), allocatable :: if0(:,:),if00(:,:),if0temp(:)
   integer (kind=IntKind), allocatable :: kcheck(:),jf0(:,:)
!
   integer (kind=IntKind) :: ib(48),ib0(48)
   integer (kind=IntKind), parameter :: ilog=6
   integer (kind=IntKind) :: i, j, n, invadd, imesh, ntot
   integer (kind=IntKind) :: iq1, iq2, iq3
   integer (kind=IntKind) :: ihg,ihc,nrot,li,isy
   integer (kind=IntKind) :: ihg0,ihc0,nc0,li0,isy0,l
!
   real (kind=RealKind), intent(in) :: bravais(3,3), k_latt(3,3)
   real (kind=RealKind), intent(in) :: ax(num_atoms)
   real (kind=RealKind), intent(in) :: ay(num_atoms)
   real (kind=RealKind), intent(in) :: az(num_atoms)
   real (kind=RealKind), intent(out) :: qmesh(3,maxk,NumMeshs)
   real (kind=RealKind) :: kslatt(3,3)
!
   real (kind=RealKind), allocatable :: rx(:,:)
!
   real (kind=RealKind) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
   real (kind=RealKind) :: v(3,48),v0(3,48)
   real (kind=RealKind) :: r0(49,3,3), r(49,3,3)
!
   real (kind=RealKind) :: wvk0(3)
   real (kind=RealKind) :: x0(3)
   real (kind=RealKind) :: proj1, proj2, proj3
!
   real (kind=RealKind), parameter :: eps = TEN2m6
!
   allocate(if0(48,num_atoms),if00(48,num_atoms),if0temp(num_atoms), &
            kcheck(maxk),jf0(num_atoms,num_atoms), rx(3,num_atoms))
!
   wvk0(1:3)=ZERO
   x0(1:3)=ZERO
!
   do i=1,3
      a1(i) =  bravais(i,1)
      a2(i) =  bravais(i,2)
      a3(i) =  bravais(i,3)
      b1(i) =  k_latt(i,1)
      b2(i) =  k_latt(i,2)
      b3(i) =  k_latt(i,3)
   enddo
   do i=1,3
      kslatt(1,i) =  b1(i)
      kslatt(2,i) =  b2(i)
      kslatt(3,i) =  b3(i)
   enddo
!
!  ===================================================================
!  determine group symmetry of the lattice
!
!      ihg .... point group of the primitive lattice, holohedral
!               group number:
!               ihg=1 stands for triclinic system
!               ihg=2 stands for monoclinic system
!               ihg=3 stands for orthorhombic system
!               ihg=4 stands for tetragonal system
!               ihg=5 stands for cubic system
!               ihg=6 stands for trigonal system
!               ihg=7 stands for hexagonal system
!      ihc .... code distinguishing between hexagonal and cubic
!               groups
!               ihc=0 stands for hexagonal groups
!               ihc=1 stands for cubic groups
!      isy .... code indicating whether the space group is
!               symmorphic or nonsymmorphic
!               isy=0 means nonsymmorphic group
!               isy=1 means symmorphic group
!               the group is considered symmorphic if for each
!               operation of the point group the sum of the 3
!               components of abs(v(n)) ( nonprimitive translation,
!               see below) is lt. 0.0005
!      li ..... code indicating whether the point group
!               of the crystal contains inversion or not
!               (operations 13 or 25 in respectively hexagonal
!               or cubic groups).
!               li=0 means: does not contain inversion
!               li.gt.0 means: there is inversion in the point
!                              group of the crystal
!      nrot.... total number of elements in the point group of the
!               crystal
!      ib ..... list of the rotations constituting the point group
!               of the crystal. the numbering is that defined in
!               worlton and warren, i.e. the one materialized in the
!               array r (see below)
!               only the first nrot elements of the aray ib are
!               meaningful
!      v  ..... nonprimitive translations (in the case of nonsymmor-
!               phic groups). v(i,n) is the i-th component
!               of the translation connected with the n-th element
!               of the point group (i.e. with the rotation
!               number ib(n) ).
!               attention: v(i) are not cartesian components,
!               they refer to the system a1,a2,a3.
!      if0 .... the function defined in maradudin, ipatova by
!               eq. (3.2.12): atom transformation table.
!               the element if0(n,kapa) means that the n-th
!               operation of the space group (i.e. operation number
!               ib(n), together with an eventual nonprimitive
!               translation  v(n)) transfers the atom kapa into the
!               atom if0(n,kapa).
!      r ...... list of the 3 x 3 rotation matrices
!               (xyz representation of the o(h) or d(6)h groups)
!               all 48 or 24 matrices are listed.
!  ===================================================================
!
!  ===================================================================
!  Determine symmetry of Bravais Lattice :.........................
!  -------------------------------------------------------------------
   call grpsymm(ilog,bravais,kslatt,num_atoms,num_atoms,ax,ay,az,itype,r, &
                ib,if0,if0temp,v,rx,ihg,ihc,nrot,li,isy)
!  -------------------------------------------------------------------
!
   if(print_level.ge.1) then
      write(ilog,*) 'GENKPT1::  IHG=',ihg,' IHC=',ihc
      write(ilog,*) 'GENKPT1::  ISY=',isy,'  LI=',li
      write(ilog,*) 'GENKPT1::  NROT = ',nrot
   endif
!
!  ===================================================================
!  Check that the symmetry is correct :...............................
!  -------------------------------------------------------------------
   call checksymm(ilog,bravais,kslatt,if0,jf0,ib,r,v,num_atoms,num_atoms,nrot)
!  -------------------------------------------------------------------
   NumRotations = nrot
   Rotation(1:nrot,1:3,1:3) = r(1:nrot,1:3,1:3)
!
   do n=nrot+1,48
      ib(n) = 0
   enddo
!
   invadd=0
!
!  ===================================================================
!  generate the bravais lattice
!  ===================================================================
   if(print_level.ge.1) then
      write(ilog,'('' GENKPT1:: The (unstrained) bravais lattice'',   &
  &                '' - used for generating the largest possible mesh'')')
   endif
!
!  ===================================================================
!  Determine symmetry of Lattice (Bravais + Basis):...................
!  -------------------------------------------------------------------
   call grpsymm(ilog,bravais,kslatt,num_atoms,1,ax,ay,az,itype,r0,    &
                ib0,if00,if0temp,v0,rx,ihg0,ihc0,nc0,li0,isy0)
!  -------------------------------------------------------------------
!
   if(print_level.ge.1) then
      write(ilog,*) 'GENKPT1::  IHG0=',ihg0,' IHC0=',ihc0
      write(ilog,*) 'GENKPT1::  ISY0=',isy0,'  LI0=',li0
      write(ilog,*) 'GENKPT1::  NC0 = ',nc0
   endif
!
!  ===================================================================
!  it is assumed that the same 'type' of symmetry operations
!  (cubic/hexagonal) will apply to the crystal as well as the bravais
!  lattice.
!  ===================================================================
!  if( ikpsym.eq.-2 ) then     ! This change is made on July 23, 2016
   if( ikpsym < 0 ) then       ! by Yang Wang
      NumRotations = 1         ! full zone
   endif
!  ===================================================================
!  / '  Generation of special points',
!  / '  (iq1,iq2,iq3 are the (generalized) monkhorst-pack parameters'
!  (they are not multiplied by 2*pi because b1,2,3  were not,either)
!  ===================================================================
   do imesh = 1,NumMeshs
      iq1 = Kdiv(1,imesh)
      iq2 = Kdiv(2,imesh)
      iq3 = Kdiv(3,imesh)
!
!     ================================================================
!     Determine Special k-points:.....................................
!     ----------------------------------------------------------------
      call spkpt(ilog,iq1,iq2,iq3,wvk0,maxk,a1,a2,a3,b1,b2,b3,        &
                 invadd,nrot,ib,r,ntot,qmesh(1,1,imesh),wght(1,imesh),&
                 kcheck,nc0,ib0,ikpsym)
!     ----------------------------------------------------------------
      NumKs(imesh)=ntot
!
      if(print_level.ge.1) then
          write(ilog,'('' GENKPT1:: IMESH='',i3,'' IQ1='',i3,         &
     &                 '' IQ2='',i3,'' IQ3='',i3,                     &
     &                 '' ==> number of special k-pts ='',i5)')       &
                imesh,iq1,iq2,iq3,NumKs(imesh)
      endif

      if(NumKs(imesh) .le. 0) then
         write(ilog,'(a,2i5)')                                        &
          ' GENKPT1:: dimensions spec. k-points NumKs,maxk ',NumKs(imesh),maxk
         call ErrorHandler('GENKPT','NK <= 0',NumKs(imesh))
      endif
!
!     ================================================================
!     set near-zeroes equal to zero:...................................
!     ================================================================
      do l=1,NumKs(imesh)
         do i = 1,3
            if (abs(qmesh(i,l,imesh)) .lt. eps) then
               qmesh(i,l,imesh) = ZERO
            endif
         enddo
!        =============================================================
!        express special points in basis:.............................
!        =============================================================
         proj1 = ZERO
         proj2 = ZERO
         proj3 = ZERO
         do i = 1,3
            proj1 = proj1 + qmesh(i,l,imesh)*a1(i)
            proj2 = proj2 + qmesh(i,l,imesh)*a2(i)
            proj3 = proj3 + qmesh(i,l,imesh)*a3(i)
         enddo
         do i = 1,3
            qmesh(i,l,imesh) = proj1*b1(i) + proj2*b2(i) + proj3*b3(i)
         enddo
      enddo
   enddo
!
   do i=1,3
      kslatt(i,1) =  b1(i)
      kslatt(i,2) =  b2(i)
      kslatt(i,3) =  b3(i)
   enddo
!
   do n=1,NumRotations
      do i=1,3
         do j=1,3
            if ( ib(n).lt.n ) then
               call ErrorHandler("GENKPT1","IB(N) < N ???",ib(n),n)
            endif
            Rotation(n,j,i) = Rotation(ib(n),j,i)
         enddo
      enddo
   enddo
!
   deallocate(if0,if00,if0temp,kcheck,jf0,rx)
!
   if(stop_routine.eq.sname) then
      call StopHandler(sname)
   endif
!
   end subroutine genkpt1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine grpsymm(ilog,bravais,kslatt,num_atoms,nat,ax,ay,az,itype,r, &
                      ib,if0,if0temp,v,rx,ihg,ihc,nc,li,isy)
!  ===================================================================
   implicit none
!
   character (len=7), parameter :: sname='grpsymm'
!
   integer (kind=IntKind), intent(in) :: ilog
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: nat
   integer (kind=IntKind), intent(inout) :: itype(num_atoms)
   integer (kind=IntKind), intent(out) :: ib(48)
   integer (kind=IntKind), intent(out) :: if0(48,num_atoms)
   integer (kind=IntKind), intent(out) :: if0temp(num_atoms)
   integer (kind=IntKind), intent(out) :: ihg,ihc,nc,li,isy
!
   real (kind=RealKind), intent(in) :: bravais(3,3),kslatt(3,3)
   real (kind=RealKind), intent(in) :: ax(num_atoms)
   real (kind=RealKind), intent(in) :: ay(num_atoms)
   real (kind=RealKind), intent(in) :: az(num_atoms)
   real (kind=RealKind), intent(out) :: r(49,3,3)
   real (kind=RealKind), intent(out) :: v(3,48)
   real (kind=RealKind), intent(out) :: rx(3,num_atoms)
!
!  ===================================================================
!  determine point group of the bravais lattice
!  -------------------------------------------------------------------
   call pgl1(bravais,kslatt,ihc,nc,ib,ihg,r)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  determine space group of the crystal
!  -------------------------------------------------------------------
   call spacegrp( ilog,kslatt,num_atoms,nat,ax,ay,az,                 &
                  itype,r,ib,if0,if0temp,v,rx,ihg,nc,li,isy)
!  -------------------------------------------------------------------

!  ===================================================================
   if(stop_routine.eq.sname) then
      call StopHandler(sname)
   endif
!
   end subroutine grpsymm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine pgl1(a,ai,ihc,nc,ib,ihg,r)
!  ===================================================================
!
!  *******************************************************************
!  written on september 11th, 1979 - from acmi complex
!  auxiliary subroutine to group1
!     subroutine pgl determines the point group of the lattice and the
!     crystal system.
!  subroutines needed: rot1, rlv3
!  a ..... direct lattice vectors
!  ai .... reciprocal lattice vectors
!  *******************************************************************
!
   implicit none
!
   integer (kind=IntKind), intent(out) :: ihc, ihg
   integer (kind=IntKind), intent(out) :: nc
   integer (kind=IntKind), intent(out) :: ib(48)
   integer (kind=IntKind) :: i, j, k, n, nr, lx
!
   real (kind=RealKind), intent(in) :: a(3,3), ai(3,3)
   real (kind=RealKind), intent(out) :: r(49,3,3)
   real (kind=RealKind) :: vr(3), xa(3)
   real (kind=RealKind) :: tr
!
!  ===================================================================
!  ihc is 0 for hexagonal groups and 1 for cubic groups.
!  ===================================================================
   ihc = 0
   nr = 24
   do
      nc = 0
!     ----------------------------------------------------------------
      call rot1 (ihc,r)
!     ----------------------------------------------------------------
      LOOP_n: do n = 1,nr
         ib(n) = 0
         tr = ZERO
!        =============================================================
!        rotate the a1, a2, a3 vectors by rotation no. n
!        =============================================================
         do k = 1,3
            do i = 1,3
               xa(i) = ZERO
               do j = 1,3
                  xa(i) = xa(i) + r(n,i,j)*a(j,k)
               enddo
            enddo
!           ----------------------------------------------------------
            call rlv3(ai,xa,vr,lx)
!           ----------------------------------------------------------
            do i = 1,3
               tr = tr + abs(vr(i))
            enddo
!           ==========================================================
!           if vr.ne.0, then xa cannot be a multiple of a lattice vector
!           ==========================================================
            if (tr .gt. 0.001d0) then
               cycle LOOP_n
            endif
         enddo
!
         nc = nc + 1
         ib(nc) = n
      enddo LOOP_n
!
!     ================================================================
!     ihg stands for holohedral group number.
!     ================================================================
      if (ihc /= 0) then
!        =============================================================
!        cubic group:
!        =============================================================
         if (nc  .lt. 4)  ihg = 1
         if (nc  .eq. 4)  ihg = 2
         if (nc  .gt. 4)  ihg = 3
         if (nc  .eq. 16) ihg = 4
         if (nc  .gt. 16) ihg = 5
         return
      else
!        =============================================================
!        hexagonal group:
!        =============================================================
         if (nc  .eq. 12) ihg = 6
         if (nc  .gt. 12) ihg = 7
         if (nc  .ge. 12) return
!        =============================================================
!        too few operations, try cubic group:
!        =============================================================
         nr  = 48
         ihc = 1
      endif
   enddo         
   end subroutine pgl1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine rot1 (ihc,r)
!  ===================================================================
!
!  *******************************************************************
!                                     written on february 17th, 1976
! generation of the x,y,z-transformation matrices 3x3 for hexagonal and
! cubic groups
! subroutines needed -- none
! this is identical with the subroutine rot of worlton-warren (in the ac
! -complex), only the way of transferring the data was changed
! input data
!      ihc...switch determining if we desire the hexagonal group (ihc=0)
!            or the cubic group (ihc=1)
! output data
!      r...the 3x3 matrices of the desired coordinate representation
!          their numbering corresponds to the symmetry elements as liste
!          in worlton-warren
!          for ihc=0 the first 24 matrices of the array r represent
!                                          the full hexagonal group d(6h
!          for ihc=1 the first 48 matrices of the array r represent
!                                          the full cubic group o(h)
!
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in) :: ihc
   integer (kind=IntKind) :: i, j, k, n, nv
!
   real (kind=RealKind), intent(out) :: r(49,3,3)
   real (kind=RealKind) :: f
!
   do j=1,3
      do i=1,3
         do n=1,49
            r(n,i,j)=ZERO
         enddo
      enddo
   enddo
   if (ihc <= 0) then
!     ================================================================
!     define the generators for the rotation matrices--hexagonal group
!     ================================================================
      f = HALF*sqrt(THREE)
      r(2,1,1) = HALF
      r(2,1,2) = - f
      r(2,2,1) = f
      r(2,2,2) = HALF
      r(7,1,1) = -HALF
      r(7,1,2) = -f
      r(7,2,1) = -f
      r(7,2,2) = HALF
      do n = 1,6
         r(n,3,3)    = ONE
         r(n+18,3,3) = ONE
         r(n+6,3,3)  = - ONE
         r(n+12,3,3) = - ONE
      enddo
!     ================================================================
!     generate the rest of the rotation matrices
!     ================================================================
      do i = 1,2
         r(1,i,i) = ONE
         do j = 1,2
            r(6,i,j) = r(2,j,i)
            do k = 1,2
               r(3,i,j)  = r(3,i,j) +  r(2,i,k)*r(2,k,j)
               r(8,i,j)  = r(8,i,j) +  r(2,i,k)*r(7,k,j)
               r(12,i,j) = r(12,i,j) + r(7,i,k)*r(2,k,j)
            enddo
         enddo
      enddo
      do i = 1,2
         do j = 1,2
            r(5,i,j) = r(3,j,i)
            do k = 1,2
               r(4,i,j)  = r(4,i,j)  + r(2,i,k)*r(3,k,j)
               r(9,i,j)  = r(9,i,j)  + r(2,i,k)*r(8,k,j)
               r(10,i,j) = r(10,i,j) + r(12,i,k)*r(3,k,j)
               r(11,i,j) = r(11,i,j) + r(12,i,k)*r(2,k,j)
            enddo
         enddo
      enddo
      do n = 1,12
         nv = n + 12
         do i = 1,2
            do j = 1,2
               r(nv,i,j) = - r(n,i,j)
            enddo
         enddo
      enddo
   else
!     ================================================================
!     define the generators for the rotation matrices--cubic group
!     ================================================================
      r(9,1,3)  = ONE
      r(9,2,1)  = ONE
      r(9,3,2)  = ONE
      r(19,1,1) = ONE
      r(19,2,3) = - ONE
      r(19,3,2) = ONE
      do i = 1,3
         r(1,i,i) = ONE
         do j = 1,3
            r(20,i,j) = r(19,j,i)
            r(5,i,j)  = r(9,j,i)
            do k  = 1,3
               r(2,i,j)  = r(2,i,j)  + r(19,i,k)*r(19,k,j)
               r(16,i,j) = r(16,i,j) + r(9,i,k)*r(19,k,j)
               r(23,i,j) = r(23,i,j) + r(19,i,k)*r(9,k,j)
            enddo
         enddo
      enddo
      do i = 1,3
         do j = 1,3
            do k = 1,3
               r(6,i,j)  = r(6,i,j)  + r(2,i,k)*r(5,k,j)
               r(7,i,j)  = r(7,i,j)  + r(16,i,k)*r(23,k,j)
               r(8,i,j)  = r(8,i,j)  + r(5,i,k)*r(2,k,j)
               r(10,i,j) = r(10,i,j) + r(2,i,k)*r(9,k,j)
               r(11,i,j) = r(11,i,j) + r(9,i,k)*r(2,k,j)
               r(12,i,j) = r(12,i,j) + r(23,i,k)*r(16,k,j)
               r(14,i,j) = r(14,i,j) + r(16,i,k)*r(2,k,j)
               r(15,i,j) = r(15,i,j) + r(2,i,k)*r(16,k,j)
               r(22,i,j) = r(22,i,j) + r(23,i,k)*r(2,k,j)
               r(24,i,j) = r(24,i,j) + r(2,i,k)*r(23,k,j)
            enddo
         enddo
      enddo
      do i=1,3
         do j=1,3
            do k=1,3
               r(3,i,j)  = r(3,i,j)  + r(5,i,k)*r(12,k,j)
               r(4,i,j)  = r(4,i,j)  + r(5,i,k)*r(10,k,j)
               r(13,i,j) = r(13,i,j) + r(23,i,k)*r(11,k,j)
               r(17,i,j) = r(17,i,j) + r(16,i,k)*r(12,k,j)
               r(18,i,j) = r(18,i,j) + r(16,i,k)*r(10,k,j)
               r(21,i,j) = r(21,i,j) + r(12,i,k)*r(15,k,j)
            enddo
         enddo
      enddo
      do n = 1,24
         nv = n + 24
         do i = 1,3
            do j = 1,3
               r(nv,i,j) = - r(n,i,j)
            enddo
         enddo
      enddo
   endif
!
   end subroutine rot1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine rlv3(ai,xb,vr,il)
!  ===================================================================
   implicit none
!
!  *******************************************************************
!            written on september 11th, 1979 - from acmi complex
!     auxiliary subroutine to group1
!     subroutine rlv removes a direct lattice vector from xb leaving the
!     remainder in vr.  if a nonzero lattice vector was removed, il is
!     made nonzero.  vr stands for v-reference.
!     ai(i,j) are the reciprocal lattice vectors, b(i) = ai(i,j),j=1,2,3
!     vr is not given in cartesian coordinates but
!     in the system a1,a2,a3.     k.k., 23.10.1979
!  *******************************************************************
   integer (kind=IntKind), intent(out) :: il
   integer (kind=IntKind) :: i, j
!
   real (kind=RealKind), intent(in) :: ai(3,3),xb(3)
   real (kind=RealKind), intent(out) :: vr(3)
   real (kind=RealKind) :: del(3)
   real (kind=RealKind) :: ts
   real (kind=RealKind), parameter :: eps = TEN2m6
!
   il = 0
   ts = ZERO
   do i = 1,3
      del(i) = ZERO 
      vr(i) = ZERO
      ts = ts + abs(xb(i))
   enddo
   if (ts .le. eps) return
   do i = 1,3
      do j = 1,3
         vr(i) = vr(i) + ai(i,j)*xb(j)
      enddo
      if (vr(i) .gt.  0.9d0) del(i) = eps
      if (vr(i) .lt. -0.9d0) del(i) = - eps
!     del is added in to eliminate roundoff errors in the function amod.
      vr(i) = vr(i) + del(i)
      il    = il + abs(vr(i))
      vr(i) = - mod(vr(i),ONE)
      vr(i) = vr(i) + del(i)
   enddo
   end subroutine rlv3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checksymm(ilog,a,ai,if0,jf0,ib,r,v,natx,nat,nc)
!  ===================================================================
   implicit none
!
   character (len=9), parameter :: sname='checksymm'
!
   logical :: ltestf,ltest,lcheck
!
   integer (kind=IntKind), intent(in) :: ilog, natx, nat, nc
   integer (kind=IntKind), intent(in) :: if0(48,natx),ib(48)
   integer (kind=IntKind), intent(out) :: jf0(natx,natx)
   integer (kind=IntKind) :: multab(48,48)
   integer (kind=IntKind) :: i, j, k, iopt, jopt, kopt, io, in, jn, is
   integer (kind=IntKind) :: ki, kj, kn, mtot, ms, itest, ii, jj
   integer (kind=IntKind) :: irotn
!
   real (kind=RealKind), intent(in) :: a(3,3),ai(3,3)
   real (kind=RealKind), intent(in) :: v(3,48),r(49,3,3)
   real (kind=RealKind) :: r0(49,3,3)
   real (kind=RealKind) :: tmp(3,3)
   real (kind=RealKind) :: s, s0
   real (kind=RealKind), parameter :: eps = TEN2m6
!
!  ===================================================================
!  test the group theoretical information (m.w. 8-94)
!  ===================================================================
   ltest=.true.
   ltestf=.true.
!
   if(print_level.ge.0) then
      write(ilog,'(/,a)')'********************************'
      write(ilog,'( a )')'* Output from CHECKSYMM        *'
      write(ilog,'(a,/)')'********************************'
   endif
!
!  test if0
   do j=1,nat
      do i=1,nat
         jf0(i,j)=0
      enddo
      do iopt=1,nc
         io=iopt
         in=if0(io,j)
         if(in.le.0.or.in.gt.natx) then
            write(*,*) 'io =',io,', j =',j,', in =',in
            write(*,*) 'natx =',natx
            call ErrorHandler(sname,'Wrong index IN',in)
         endif
         jf0(in,j)=jf0(in,j)+1
      enddo
   enddo
!
   do j=1,nat
      in=jf0(j,j)
!     test that nc/(equivalent atoms) is an integer
      if(mod(nc,in).ne.0) then
         ltestf=.false.
      endif
!
      is=0
      do i=1,nat
         is=is+jf0(i,j)
      enddo
!     test that one atom is generated for each operation
      if(is.ne.nc) then
         ltestf=.false.
      endif
!
      do i=1,nat
!        test that all equivalent atoms are generated equally
         if(jf0(i,j).ne.0) then
            if(jf0(i,j).ne.in) then
               ltestf=.false.
            endif
!           test that all equivalent atoms are generated equivalently
            if(jf0(i,j).ne.jf0(j,i))then
               ltestf=.false.
            endif
         endif
      enddo
   enddo
!
   if(.not.ltestf) then
      write(ilog,'(1x//,1x,70(''-'')/)')
      write(ilog,'(a)') ' ERROR: if0 is incorrect'
      write(ilog,'(a)') ' generated equivalent number of atoms (jf0):'
      do j=1,nat
         write(ilog,'(1x/,'' Atom'',i4,'':'',12i4/,(10x,12i4))')  &
               j,(jf0(i,j),i=1,nat)
      enddo
      write(ilog,'(1x/,1x,70(''-'')/)')
   endif
!  test the closure of the group
!  first put rotation matrices in lattice basis
   do iopt=1,nc
      in=ib(iopt)
      do j=1,3
         do i=1,3
            s=ZERO
            do k=1,3
               s=s+r(in,i,k)*a(k,j)
            enddo
            tmp(i,j)=s
         enddo
      enddo
      do j=1,3
         do i=1,3
            s=ZERO
            do k=1,3
               s=s+ai(i,k)*tmp(k,j)
            enddo
            r0(in,i,j)=s
         enddo
      enddo
   enddo
!  set-up multiplication table
   do iopt=1,nc
      in=ib(iopt)
      do jopt=1,nc
         jn=ib(jopt)
         multab(jopt,iopt)=0
!        get r0(in)r0(jn)
         do j=1,3
            do i=1,3
               s=ZERO
               do k=1,3
                  s=s+r0(in,i,k)*r0(jn,k,j)
               enddo
               tmp(i,j)=s
            enddo
         enddo
!        determine which element of the group this is
         do kopt=1,nc
            kn=ib(kopt)
            lcheck=.true.
            do kj=1,3
               do ki=1,3
                  if(abs(r0(kn,ki,kj)-tmp(ki,kj)).gt.eps) lcheck=.false.
               enddo
            enddo
            if(lcheck) then
               multab(jopt,iopt)=kopt
               irotn=kopt
               go to 901
            endif
         enddo
         ltest=.false.
  901    continue
!        check nonsymorphic part
         LOOP_i: do i=1,3
            s=ZERO
            do j=1,3
               s=s+r0(in,i,j)*v(j,jopt)
            enddo
            s0=s+v(i,iopt)
            s=abs(s0-v(i,irotn))
            do itest=0,2
               if(abs(s).lt.eps) then
                  cycle LOOP_i
               endif
               s=s-ONE
            enddo
            ltest=.false.
            multab(jopt,iopt)=-multab(jopt,iopt)
            exit LOOP_i
         enddo LOOP_i
      enddo
   enddo
!  test the closure property
   mtot=(nc*(nc+1))/2
   do i=1,nc
      ms=0
      do j=1,nc
         ms=ms+abs(multab(j,i))
      enddo
      if(ms.ne.mtot) then
         ltest=.false.
      endif
   enddo
   do j=1,nc
      ms=0
      do i=1,nc
         ms=ms+abs(multab(j,i))
      enddo
      if(ms.ne.mtot) then
         ltest=.false.
      endif
   enddo

   if(print_level >= 0 .or. MyPE == 0) then
      if(ltest.and.ltestf) then
         write(ilog,'(''Symmetry is correct'')')
      else
         write(ilog,'(''***** symmetry checks failed *********'')')
         write(ilog,'(1x/)')
         if(.not.ltest) then
            write(ilog,'(a)') 'Multiplication table (- denotes'//    &
                              ' error in non-symmorphic part)'
            write(ilog,'(1x/)')
            do i=1,nc
               write(ilog,'(1x/,'' i='',i3,'':'',12i4/,(7x,12i4))')   &
                     i,(multab(j,i),j=1,nc)
            enddo
         endif
      endif
!     write out the operations in lattice form:
      if((.not.ltest).or.(.not.ltestf)) then
         write(ilog, '(1x/,70(''-'')/,'' Operations given in lattice form'')')
         do iopt=1,nc
            in=ib(iopt)
            write(ilog,'(1x/,'' operation number'',i3,'' ('',i3,'')'')') &
                  iopt,in
            do ii=1,3
               write(ilog,'(5x,''( '',3f8.3,'' )         ( '',f8.3,'' )'')')&
                    (r0(in,ii,jj),jj=1,3),v(ii,iopt)
            enddo
         enddo
         write(ilog,'(1x/)')
         write(ilog,'(''***** symmetry checks failed *********'')')
      endif
      write(ilog,'(80(''-''))')
   endif
   end subroutine checksymm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine spacegrp(ilog,ai,natx,nat,x_1,x_2,x_3,                &
                       ity,r,ib,if0,if0temp,v,rx,ihg,nc,li,isy)
!  ===================================================================
   implicit none
!
   character (len=8), parameter :: sname='spacegrp'
   character (len=12) :: icst(7)
!
   logical :: lcell
!
   integer (kind=IntKind), intent(in) :: ilog
   integer (kind=IntKind), intent(in) :: natx, nat, ihg
   integer (kind=IntKind), intent(inout) :: nc
   integer (kind=IntKind), intent(inout) :: ity(natx)
   integer (kind=IntKind), intent(inout) :: ib(48)
   integer (kind=IntKind), intent(out) :: if0(48,natx),if0temp(natx)
   integer (kind=IntKind), intent(out) :: li, isy
   integer (kind=IntKind) :: index(48)
   integer (kind=IntKind) :: n, i, k, icell, l, k1, k2, k3, k4, ni, il
!
   real (kind=RealKind), intent(in) :: ai(3,3)
   real (kind=RealKind), intent(in) :: x_1(natx)
   real (kind=RealKind), intent(in) :: x_2(natx)
   real (kind=RealKind), intent(in) :: x_3(natx)
   real (kind=RealKind), intent(in) :: r(49,3,3)
   real (kind=RealKind), intent(out) :: v(3,48)
   real (kind=RealKind), intent(out) :: rx(3,natx)
   real (kind=RealKind) :: vr(3),vt(3),xb(3)
   real (kind=RealKind) :: da, dif, vs, epsil
   real (kind=RealKind), parameter :: eps = TEN2m6
!
   data icst /'triclinic','monoclinic','orthorhombic','tetragonal',  &
              'cubic','trigonal','hexagonal'/
!
   if(print_level.ge.0) then
      write(ilog,'(/,a)')'********************************'
      write(ilog,'( a )')'* Output from SPACEGRP         *'
      write(ilog,'(a,/)')'********************************'
   endif
!
   lcell = .false.
!  zero index array
   do n=1,48
      index(n) = 0
   enddo
!  loop over elements of the group of the lattice
   do n=1,nc
      icell = 0
!     rotate the atoms x to position rx
      l=ib(n)
      do k=1,nat
         do i=1,3
            rx(i,k)=ZERO
            rx(i,k)=rx(i,k)+r(l,i,1)*x_1(k)
            rx(i,k)=rx(i,k)+r(l,i,2)*x_2(k)
            rx(i,k)=rx(i,k)+r(l,i,3)*x_3(k)
         enddo
      enddo
!     consider the first rotated atom k1=1
      LOOP_do: do
         k1=1
!        loop over original atoms in the cell
         LOOP_k2: do k2 = 1,nat
            if (ity(k1) .ne. ity(k2)) then
               cycle LOOP_k2
            endif
!           ==========================================================
!           check to see whether k2 rotates into k1 under n with nonsymmorphic 
!           translation vr
!           ==========================================================
            xb(1)=rx(1,k1)-x_1(k2)
            xb(2)=rx(2,k1)-x_2(k2)
            xb(3)=rx(3,k1)-x_3(k2)
!           ==========================================================
!           subroutine rlv removes a direct lattice vector from xb leaving the
!           remainder in vr.  if a nonzero lattice vector was removed, il is
!           made nonzero.
!           vr stands for v-reference
!           ----------------------------------------------------------
            call rlv3 (ai,xb,vr,il)
!           ----------------------------------------------------------
!           loop over rotated atoms
            LOOP_k3: do k3 = 1,nat
!              loop over original atoms
               LOOP_k4: do k4 = 1,nat
                  if (ity(k3) .ne. ity(k4)) then
                     cycle LOOP_k4
                  endif
!                 ====================================================
!                 check to see whether each atom in the original system 
!                 is found in the rotated system with corresponding 
!                 nonprimitive translation vr. if a single match cannot 
!                 be made, then k2 does not rotate into k1 under n. if
!                 all the atoms match up, then if0(n,k) is the atom 
!                 transformation table with nonsymmorphic translation 
!                 vr. we do not stop here, however, but check to see 
!                 whether there exits more than one nonsymmorphic 
!                 translation vr corresponding to the group element n 
!                 which could be the case in a supercell geometry. in 
!                 this case, we artificially change the internal
!                 symmetry of the cell by 'painting atom number 1 a 
!                 different color'
!                 ====================================================
                  xb(1)=rx(1,k3)-x_1(k4)
                  xb(2)=rx(2,k3)-x_2(k4)
                  xb(3)=rx(3,k3)-x_3(k4)
!                 ====================================================
!                 vt stands for v-test
!                 ----------------------------------------------------
                  call rlv3 (ai,xb,vt,il)
!                 ----------------------------------------------------
                  dif = ZERO
                  do i = 1,3
                     da = abs(vr(i) - vt(i)) + eps
                     dif = dif + mod(da,1.d0)
                  enddo
                  if (dif .gt. 0.001d0) then
                     cycle LOOP_k4
                  endif
!                 ====================================================
!                 atom k4 rotates into k3 with nonsymmorphic translation vr
!                 ====================================================
                  if0temp(k3)=k4
!                 ====================================================
!                 check for the next rotated atom k3
!                 ====================================================
                  cycle LOOP_k3
               enddo LOOP_k4
!              =======================================================
!              cannot find all the atoms in the rotated system for this
!              choice of vr
!              =======================================================
               cycle LOOP_k2
            enddo LOOP_k3
!           ==========================================================
!           an exact match is made: n is an element of the group with
!           nonsymmorphic translation vr which is not given in cartesian 
!           coordinates but in the system a1,a2,a3.
!           ==========================================================
            icell = icell + 1
            index(n) = 1
            do i=1,3
               v(i,n) = vr(i)
            enddo
            do k=1,nat
               if0(n,k) = if0temp(k)
            enddo
         enddo LOOP_k2
         if (icell .gt. 1) then
            if (lcell) then
               call ErrorHandler(sname,'cannot reduce symmetries of the cell')
            endif
            lcell = .true.
            if (print_level >= 1 .or. MyPE == 0) then
               write(ilog,'(5x,''A nonprimitive supercell was chosen with'', &
     &               i5,'' times'',/,5x,                                     &
     &               ''as many more space group operations than that for'',  &
     &               /,5x,''the corresponding primitive unit cell'',/)')icell
            endif
            icell = 0
!           change the internal symmetry by changing the atomic number of atom 1
            ity(1) = -ity(1)
            cycle LOOP_do
         endif
         exit LOOP_do
      enddo LOOP_do
   enddo
!  change back the atomic number of atom 1
   if (lcell) then
      ity(1) = -ity(1)
   endif
   i  = 0
!  inverse of hexagonal group 
   ni = 13
!  inverse of cubic group
   if (ihg .lt. 6) then
      ni = 25
   endif
   li = 0
!  reshuffle
   LOOP_n: do n=1,nc
      if (index(n) .eq. 0) then
         cycle LOOP_n
      endif
      i = i + 1
      ib(i) = ib(n)
      if (ib(i) .eq. ni) li = i
      do k = 1,nat
         if0(i,k) = if0(n,k)
      enddo
      v(1,i)   = v(1,n)
      v(2,i)   = v(2,n)
      v(3,i)   = v(3,n)
   enddo LOOP_n
   nc = i
   if (print_level >= 1 .or. MyPE == 0) then
      if (.not.(ihg .eq. 7 .and. nc .eq. 24 .or.                       &
          ihg .eq. 5 .and. nc .eq. 48)) then
         write(ilog,60) icst(ihg),(ib(i),i=1,nc)
60       format (' The crystal system is ',a,' with operations:'/2(3x,24i3/))
         write (ilog,'(a,a,a,a)') "The point group of the crystal",    &
                                  " is the full ",icst(ihg)," group"
      endif
      write (ilog,'(a,a,a,a)') "The point group of the crystal",       &
                               " is the full ",icst(ihg)," group"
   endif
!  check symmorphic nature of group
   vs = ZERO
   do n = 1,nc
      do i = 1,3
         vs = vs + abs(v(i,n))
      enddo
   enddo
!  epsil = 0.0005d0
   epsil = eps
   if (vs <= epsil) then
      if(print_level.ge.0) then
         write (ilog,'(''Space group is symmorphic'')')
      endif
      isy = 1
   else
      if(print_level.ge.0) then
         write(ilog,'(''Space group is non-symmorphic'',    &
      &        '' OR a non standard origin of coordinates was used.'')')
      endif 
      isy = 0
   endif
   if(print_level.ge.0) then
      write(ilog,'(80(''-''))')
   endif
!
   end subroutine spacegrp
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine spkpt(ilog,iq1,iq2,iq3,wvk0,maxk,a1,a2,a3,b1,b2,b3,     &
                    inv,nc,ib,r,ntot,wvkl,lwght,                      &
                    kcheck,ncbrav,ibrav,istriz)
!  ===================================================================
!
!  *******************************************************************
!  to find the shortest wave-vector.
!  the rotations of the bravais lattice are applied to the monkhorst/pack 
!  mesh in order to find all k-points that are related by symmetry
!
!  input data:
!      iq1,iq2,iq3 .. parameter q of monkhorst and pack,
!               generalized and different for the 3 directions b1,
!               b2 and b3
!      wvk0 ... the 'arbitrary' shift of the whole mesh, denoted k0
!               in macdonald. wvk0 = 0 corresponds to the original
!               scheme of monkhorst and pack.
!               units: 2pi/(units of length  used in a1, a2, a3),
!               i.e. the same  units as the generated special points.
!      maxk ... variable dimension of the (output) arrays wvkl,
!               lwght,lrot, i.e. space reserved for the special
!               points and accessories.
!               it has to be maxk .ge. ntot (total number of special
!               points), but the subroutine does not check on this.
!      istriz . indicates whether additional mesh points should be
!               generated by applying group operations to the mesh
!               to the mesh.
!               abs(istriz) .eq. 1 means symmetrize
!               abs(istriz) .ne. 1 means don't symmetrize
!                   istriz  .lt. 0 means don't exclude equivalent points
!      b1,b2,b3 .. reciprocal lattice vectors, not multiplied by
!               any 2pi (in units reciprocal to those of a1,a2,a3)
!      inv .... code indicating whether we wish  to add the inversion
!               to the point group of the crystal or not (in the
!               case that the point group does not contain any).
!               inv=0 means: do not add inversion
!               inv.ne.0 means: add the inversion
!               inv.ne.0 should be the standard choice when spkpt
!               is used in reciprocal space - in order to make
!               use of the hermiticity of hamiltonian.
!               when used in direct space, the right choice of inv
!               will depend on the nature of the physical problem.
!               in the cases where the inversion is added by the
!               switch inv, the list ib will not be modified but in
!               the output list lrot some of the operations will
!               appear with negative sign; this means that they have
!               to be applied multiplied by inversion.
!      nc ..... total number of elements in the point group of the
!               crystal
!      ib ..... list of the rotations constituting the point group
!               of the crystal. the numbering is that defined in
!               worlton and warren, i.e. the one materialized in the
!               array r (see below)
!               only the first nc elements of the array ib are
!               meaningful
!      r ...... list of the 3 x 3 rotation matrices
!               (xyz representation of the o(h) or d(6)h groups)
!               all 48 or 24 matrices are listed.
!      ncbrav . total number of elements in rbrav
!      ibrav .. list of ncbrav operations of the bravais lattice
!  output data:
!      ntot ... total number of special points
!               if ntot appears negative, this is an error signal
!               which means that the dimension maxk was chosen
!               too small so that the arrays wvkl etc. cannot
!               accomodate all the generated special points.
!               in this case the arrays will be filled up to maxk
!               and further generation of new points will be
!               interrupted.
!      wvkl ... list of special points .
!               cartesian coordinates and not multiplied by 2*pi.
!               only the first ntot vectors are meaningful
!               although no 2 points from the list are equivalent
!               by symmetry, this subroutine still has a kind of
!               'beauty defect': the points finally
!               selected are not necessarily situated in a
!               'compact' irreducible brill.zone; they might lie in
!               different irreducible parts of the b.z. - but they
!               do represent an irreducible set for integration
!               over the entire b.z.
!     lwght ... the list of weights of the corresponding points.
!               these weights are not normalized (just integers)
!      lrot ... for each special point the 'unfolding rotations'
!               are listed. if e.g. the weight of the i-th special
!               point is lwght(i), then the rotations with numbers
!               lrot(j,i), j=1,2,...,lwght(i) will 'spread' this
!               single point from the irreducible part of b.z. into
!               several points in an elementary unit cell
!               (parallelopiped) of the reciprocal space.
!               some operation numbers in the list lrot may appear
!               negative, this means that the corresponding rotation
!               has to be applied with inversion (the latter having
!               been artificially added as symmetry operation in
!               case inv.ne.0).no other effort was taken,to renumber
!               the rotations with minus sign or to extend the
!               list of the point-group operations in the list nb.
!  *******************************************************************
!
   implicit none
!
   character (len=5), parameter :: sname='spkpt'
!
   integer (kind=IntKind), intent(in) :: ilog
   integer (kind=IntKind), intent(in) :: iq1, iq2, iq3
   integer (kind=IntKind), intent(in) :: maxk
   integer (kind=IntKind), intent(in) :: inv, nc, ncbrav, istriz
   integer (kind=IntKind), intent(in) :: ib(48),ibrav(48)
   integer (kind=IntKind), intent(out) :: ntot
   integer (kind=IntKind), intent(out) :: lwght(maxk),kcheck(maxk)
   integer (kind=IntKind) :: iout
   integer (kind=IntKind) :: iqp1, iqp2, iqp3
   integer (kind=IntKind) :: i1, i2, i3
   integer (kind=IntKind) :: kcount, nplane
   integer (kind=IntKind) :: i, j, k, iop, iremov, n, ibsign, m
   integer (kind=IntKind), parameter ::  nrsdir=100
!
   real (kind=RealKind), intent(in) :: wvk0(3)
   real (kind=RealKind), intent(in) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
   real (kind=RealKind), intent(in) :: r(49,3,3)
   real (kind=RealKind), intent(out) :: wvkl(3,maxk)
   real (kind=RealKind) :: ur1, ur2, ur3
   real (kind=RealKind) :: sum, diff
   real (kind=RealKind) :: wvk(3),wva(3),rsdir(4,nrsdir),proja(3),projb(3)
   real (kind=RealKind), parameter :: eps = TEN2m6
!
   iout = ilog
!  ===================================================================
!  define the 1st brillouin zone
!  -------------------------------------------------------------------
   call bzdefi(b1,b2,b3,rsdir,nrsdir,nplane,iout)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  generation of the mesh (they are not multiplied by 2*pi)
!  by the monkhorst/pack algorithm, supplemented by all rotations
!  ===================================================================
!
!  ===================================================================
!  zero index array
!  ===================================================================
   do k=1,maxk
      kcheck(k) = 0
   enddo

!  ===================================================================
!  odd iq's should not end up on origin
!  ===================================================================
   iqp1=iq1
   iqp2=iq2
   iqp3=iq3
   if( mod(iq1,2).eq.1 ) iqp1=iq1-1
   if( mod(iq2,2).eq.1 ) iqp2=iq2-1
   if( mod(iq3,2).eq.1 ) iqp3=iq3-1
!
   kcount = 0
   do i1=1,iq1
      do i2=1,iq2
         do i3=1,iq3
            ur1=dble(1 + iqp1 - 2*i1)/dble(2*iq1)
            ur2=dble(1 + iqp2 - 2*i2)/dble(2*iq2)
            ur3=dble(1 + iqp3 - 2*i3)/dble(2*iq3)
            do i=1,3
               wvk(i) = ur1*b1(i) + ur2*b2(i) + ur3*b3(i) + wvk0(i)
            enddo
!           reduce wvk to the 1st brillouin zone
!           ----------------------------------------------------------
            call bzrduc(wvk,a1,a2,a3,b1,b2,b3,rsdir,nrsdir,nplane,iout)
!           ----------------------------------------------------------
            if (abs(istriz) .eq. 1) then
!              apply all the bravais lattice operations to wvk
               LOOP_iop: do iop = 1,ncbrav
                  do i=1,3
                     wva(i) = ZERO
                     do j = 1,3
                        wva(i) = wva(i) + r(ibrav(iop),i,j)*wvk(j)
                     enddo
                  enddo
!                 ====================================================
!                 check that wva is inside the 1st bz.
!                 ====================================================
                  if (inbz(wva,rsdir,nrsdir,nplane) .eq. 0) then
                     call ErrorHandler(sname,'rotated k point outside 1st bz')
                  endif
!                 ====================================================
!                 include the very first point
!                 ====================================================
                  if (i1.eq.1.and.i2.eq.1.and.i3.eq.1.and.iop.eq.1) then
                     kcount = kcount + 1
                     wvkl(1,kcount) = wva(1)
                     wvkl(2,kcount) = wva(2)
                     wvkl(3,kcount) = wva(3)
                     kcheck(kcount) = 1
                  else
!                    =================================================
!                    has this point been encountered before?
!                    =================================================
                     do j=1,kcount
                        sum = abs(wva(1) - wvkl(1,j)) +               &
                              abs(wva(2) - wvkl(2,j)) +               &
                              abs(wva(3) - wvkl(3,j))
                        if (sum .lt. eps) then
                           cycle LOOP_iop
                        endif
                     enddo
!                    =================================================
!                    this is a new point
!                    =================================================
                     kcount = kcount + 1
                     if (kcount .gt. maxk) then
                         call ErrorHandler(sname,'kcount > maxk',kcount,maxk)
                     endif
                     wvkl(1,kcount) = wva(1)
                     wvkl(2,kcount) = wva(2)
                     wvkl(3,kcount) = wva(3)
                     kcheck(kcount) = 1
                  endif            
               enddo LOOP_iop
            else
               kcount = kcount + 1
               if ( kcount .gt. maxk ) then
                  write(6,'('' 2nd kcount,maxk '',2i10)') kcount,maxk
                  call ErrorHandler(sname,'kcount > maxk',kcount,maxk)
               endif
               wvkl(1,kcount) = wvk(1)
               wvkl(2,kcount) = wvk(2)
               wvkl(3,kcount) = wvk(3)
               kcheck(kcount) = 1
            endif
         enddo
      enddo
   enddo
!
   if(print_level.ge.0) then
      write(ilog,'(/,a)')'********************************'
      write(ilog,'( a )')'* Output from SPKPT            *'
      write(ilog,'(a,/)')'********************************'
      write(ilog,'(''Wavevector mesh contains '',i7,'' points'')') kcount
   endif
!
   if (abs(istriz) .eq. 1) then
!     ================================================================
!     figure out if any special point difference (k - k') is an
!     integral multiple of a reciprocal-space vector
!     ================================================================
      iremov = 0
      do i=1,(kcount-1)
         do k=1,3
            wva(k) = wvkl(k,i)
         enddo
!        =============================================================
!        project wva onto b1,2,3:
!        =============================================================
         proja(1) = ZERO
         proja(2) = ZERO
         proja(3) = ZERO
         do k=1,3
            proja(1) = proja(1) + wva(k)*a1(k)
            proja(2) = proja(2) + wva(k)*a2(k)
            proja(3) = proja(3) + wva(k)*a3(k)
         enddo
!        =============================================================
!        loop over the rest of the mesh points
!        =============================================================
         LOOP_j: do j=(i+1),kcount
            if (kcheck(j) .eq. 0) then
               cycle LOOP_j
            endif
            do k=1,3
               wvk(k) = wvkl(k,j)
            enddo
!
!           ==========================================================
!           project wvk onto b1,2,3:
!           ==========================================================
            projb(1) = ZERO
            projb(2) = ZERO
            projb(3) = ZERO
            do k=1,3
               projb(1) = projb(1) + wvk(k)*a1(k)
               projb(2) = projb(2) + wvk(k)*a2(k)
               projb(3) = projb(3) + wvk(k)*a3(k)
            enddo
!
!           ==========================================================
!           check whether (proja-projb) is integral?
!           ================================================================
            do k=1,3
               diff = proja(k) - projb(k)
               if (abs(nint(diff)-diff) .gt. eps) then
                  cycle LOOP_j
               endif
            enddo
!           ==========================================================
!           diff is integral: remove wvk from mesh
!           ==========================================================
            kcheck(j) = 0
            iremov = iremov + 1
         enddo LOOP_j
      enddo
!
      if (iremov .gt. 0) then
         write(iout,'(2x,''some of these mesh points are related by '',&
              &          ''lattice translation vectors''/1x,i6,        &
              &          '' of the mesh points removed.''/)')iremov
      endif
!     reshuffle
      j=0
      LOOP_i: do i=1,kcount
         if (kcheck(i) .eq. 0) then
            cycle LOOP_i
         endif
         j=j+1
         do k=1,3
            wvkl(k,j) = wvkl(k,i)
         enddo
         kcheck(j) = kcheck(i)
      enddo LOOP_i
      kcount = j
   endif
!
!  ===================================================================
!  in the mesh of wavevectors, now search for equivalent points:
!  the inversion (time reversal !) may be used.
!  ===================================================================
   ntot = 0
   LOOP_k: do k=1,kcount
      if (kcheck(k) .eq. 0) then
         cycle LOOP_k
      endif
      ntot = ntot + 1
      lwght(ntot) = 1
      if (istriz.lt.0 .or. k.eq.kcount) then
         cycle LOOP_k
      endif
!     find all the equivalent points
      do n=1,nc
!        rotate:
         do i=1,3
            wva(i) = ZERO
            do j=1,3
               wva(i) = wva(i) + r(ib(n),i,j)*wvkl(j,k)
            enddo
         enddo
         ibsign = 1
!        =============================================================
!        find out which point this is (could have been discarded!)
!        =============================================================
         LOOP_do2: do
            LOOP_m: do m=k+1,kcount
               if (kcheck(m) .eq. 0) then
                  cycle LOOP_m
               endif
               sum = abs(wva(1)-wvkl(1,m))+abs(wva(2)-wvkl(2,m)) +      &
                     abs(wva(3)-wvkl(3,m))
               if (sum .lt. eps) then
                  lwght(ntot) = lwght(ntot) + 1
                  kcheck(m) = 0
                  exit LOOP_m
               endif
            enddo LOOP_m
!
!           ==========================================================
!           use time-reversal symmetry if the inversion does not exit
!           ==========================================================
            if (inv .eq. 1 .and. ibsign .eq. 1) then
               ibsign = -1
               do i=1,3
                  wva(i) = -wva(i)
               enddo
            else
               exit LOOP_do2
            endif
         enddo LOOP_do2
      enddo
   enddo LOOP_k
!  reshuffle
   j=0
   LOOP_i2: do i=1,kcount
      if (kcheck(i) .eq. 0) then
         cycle LOOP_i2
      endif
      j=j+1
      do k=1,3
         wvkl(k,j) = wvkl(k,i)
      enddo
      kcheck(j) = kcheck(i)
   enddo LOOP_i2
   kcount = j
!
   if(print_level.ge.0) then
      write(ilog,'(80(''-''))')
   endif
!
   if (kcount .ne. ntot) then
       write(6,'('' kcount,ntot '',2i15)') kcount,ntot
       call ErrorHandler(sname,'invalid number of k points')
   endif
!
!  ===================================================================
!  total number of special points: ntot
!  before using the list wvkl as wave vectors, they have to be
!  multiplied by 2*pi
!  the list of weights lwght is not normalized
!  ===================================================================
   if(stop_routine.eq.sname) then
      call StopHandler(sname)
   endif
!
   end subroutine spkpt
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bzdefi(b1,b2,b3,rsdir,nrsdir,nplane,iout)
!  ===================================================================
!
!  *******************************************************************
!  find the vectors whose halves define the 1st brillouin zone
!
!  on output, nplane tells how many elements of rsdir contain
!  normal vectors defining the planes
!  method: starting with the parallelopiped spanned by b1,2,3
!  around the origin, vectors inside a sufficiently large sphere
!  are tested to see whether the planes at 1/2*b will
!  further confine the 1bz.
!  the resulting vectors are not cleaned to avoid redundant planes.
!  *******************************************************************
!
   implicit none
!
   character (len=6), parameter :: sname='bzdefi'
!
   integer (kind=IntKind), intent(in) :: nrsdir
   integer (kind=IntKind), intent(out) :: nplane
   integer (kind=IntKind), intent(in) :: iout
   integer (kind=IntKind) :: i1, i2, i3
   integer (kind=IntKind) :: n1, n2, n3
   integer (kind=IntKind) :: nb1, nb2, nb3
   integer (kind=IntKind) :: nnb1, nnb2, nnb3
   integer (kind=IntKind) :: i, j, n
!
   real (kind=RealKind), intent(in) :: b1(3),b2(3),b3(3)
   real (kind=RealKind), intent(out) :: rsdir(4,nrsdir)
   real (kind=RealKind) :: bvec(3)
   real (kind=RealKind) :: b1len, b2len, b3len, bmax
   real (kind=RealKind) :: project
   real (kind=RealKind), parameter :: eps = TEN2m6
!
!  once initialized, we do not repeat the calculation
   if (initlz .ne. 0) then
      return
   endif
!
   initlz = 1
   b1len = b1(1)**2 + b1(2)**2 + b1(3)**2
   b2len = b2(1)**2 + b2(2)**2 + b2(3)**2
   b3len = b3(1)**2 + b3(2)**2 + b3(3)**2
!
!  ===================================================================
!  lattice containing entirely the brillouin zone
!  ===================================================================
   bmax = b1len + b2len + b3len
   nb1 = int( sqrt(bmax/b1len) + 1.0d-6)
   nb2 = int( sqrt(bmax/b2len) + 1.0d-6)
   nb3 = int( sqrt(bmax/b3len) + 1.0d-6)
   do i = 1,nrsdir
      do j = 1,4
         rsdir(j,i) = ZERO
      enddo
   enddo
!
!  ===================================================================
!  1bz is certainly confined inside the 1/2(b1,b2,b3) parallelopiped
!  ===================================================================
   do i = 1,3
      rsdir(i,1) = b1(i)
      rsdir(i,2) = b2(i)
      rsdir(i,3) = b3(i)
   enddo
   rsdir(4,1) = b1len
   rsdir(4,2) = b2len
   rsdir(4,3) = b3len
!
!  ===================================================================
!  starting confinement: 3 planes
!  ===================================================================
   nplane = 3
!
   nnb1 = 2*nb1 + 1
   nnb2 = 2*nb2 + 1
   nnb3 = 2*nb3 + 1
   do n1 = 1,nnb1
      i1 = nb1 + 1 - n1
      do n2 = 1,nnb2
         i2 = nb2 + 1 - n2
         LOOP_n3: do n3 = 1,nnb3
            i3 = nb3 + 1 - n3
            if (i1.eq.0 .and. i2.eq.0 .and. i3.eq.0) then
               cycle LOOP_n3
            endif
            do i = 1,3
               bvec(i) = dble(i1)*b1(i) + dble(i2)*b2(i) + dble(i3)*b3(i)
            enddo
!           ==========================================================
!           does the plane of 1/2*bvec narrow down the 1bz ?
!           ==========================================================
            do n = 1,nplane
               project = HALF*(rsdir(1,n)*bvec(1) + rsdir(2,n)*bvec(2) + &
                               rsdir(3,n)*bvec(3) ) / rsdir(4,n)
!              =======================================================
!              1/2*bvec is outside the bz - skip this direction
!              the 1.0e-6 takes care of single points touching the bz, and
!              of the -(plane)
!              =======================================================
               if (abs(project) .gt. HALF-eps) then
                  cycle LOOP_n3
               endif
            enddo
!           ==========================================================
!           1/2*bvec further confines the 1bz - include into rsdir
!           ==========================================================
            nplane = nplane + 1
            if (nplane .gt. nrsdir) then
               write (iout,'('' BZDEFI::  *** fatal error ***'')')
               write (iout,'('' BZDEFI::  many planes, nrsdir ='',i5)') nrsdir
               call ErrorHandler(sname,'too many planes',nplane)
            endif
            do i = 1,3
               rsdir(i,nplane) = bvec(i)
            enddo
!           ==========================================================
!           length squared
!           ==========================================================
            rsdir(4,nplane) = bvec(1)**2 + bvec(2)**2 + bvec(3)**2
         enddo LOOP_n3
      enddo
   enddo              
!
!  ===================================================================
!  print information
!  ===================================================================
   if(print_level.ge.1) then
      write(iout,'('' BZDEFI::  1st brillouin zone is confined '',         &
          &        '' by (at most)'',i3,'' planes as defined by the +/- '',&
          &        '' halves of the vectors:'',/,(1x,3f10.4))')            &
          &        nplane,((rsdir(i,n),i=1,3),n=1,nplane)
   endif
!
   if(stop_routine.eq.sname) then
      call StopHandler(sname)
   endif
   end subroutine bzdefi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bzrduc(wvk,a1,a2,a3,b1,b2,b3,rsdir,nrsdir,nplane,iout)
!  ===================================================================
!
!  *******************************************************************
!  reduce wvk to lie entirely within the 1st brillouin zone
!  by adding b-vectors
!  look around +/- "nzones" to locate vector
!  *******************************************************************
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nrsdir, nplane, iout
   integer (kind=IntKind) :: nn1, nn2, nn3
   integer (kind=IntKind) :: n1, n2, n3
   integer (kind=IntKind) :: i1, i2, i3, i
   integer (kind=IntKind), parameter :: yes = 1
   integer (kind=IntKind), parameter :: no = 0
   integer (kind=IntKind), parameter :: nzones=2
   integer (kind=IntKind), parameter :: nnn=2*nzones+1
   integer (kind=IntKind), parameter :: nn=nzones+1
!
   real (kind=RealKind), intent(in) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
   real (kind=RealKind), intent(in) :: rsdir(4,nrsdir)
   real (kind=RealKind), intent(inout) :: wvk(3)
   real (kind=RealKind) :: wva(3),wb(3)
!
!  ===================================================================
!  check if wvk already inside 1bz
!  ===================================================================
   if (inbz(wvk,rsdir,nrsdir,nplane) .eq. yes) then
      return
   endif
!
!  ===================================================================
!  express wvk in the basis of b1,2,3.
!  this permits an estimate of how far wvk is from the 1bz.
!  ===================================================================
   wb(1) = wvk(1)*a1(1) + wvk(2)*a1(2) + wvk(3)*a1(3)
   wb(2) = wvk(1)*a2(1) + wvk(2)*a2(2) + wvk(3)*a2(3)
   wb(3) = wvk(1)*a3(1) + wvk(2)*a3(2) + wvk(3)*a3(3)
   nn1 = nint(wb(1))
   nn2 = nint(wb(2))
   nn3 = nint(wb(3))
!
!  ===================================================================
!  look around the estimated vector for the one truly inside the 1bz
!  ===================================================================
   do n1 = 1,nnn
      i1 = nn - n1 - nn1
      do n2 = 1,nnn
         i2 = nn - n2 - nn2
         do n3 = 1,nnn
            i3 = nn - n3 - nn3
            do i = 1,3
               wva(i) = wvk(i) + dble(i1)*b1(i) + dble(i2)*b2(i) +    &
                        dble(i3)*b3(i)
            enddo
            if (inbz(wva,rsdir,nrsdir,nplane) .eq. yes) then
!              =======================================================
!              the reduced vector
!              =======================================================
               do i = 1,3
                  wvk(i) = wva(i)
               enddo
               return
            endif
         enddo
      enddo
   enddo
!  ===================================================================
!  fatal error
!  ===================================================================
   write(iout,'(''Subroutine bzrduc *** fatal error ***  wavevector '',&
  &             3f10.4,'' could not be reduced to the 1bz'')')wvk
!  ---------------------------------------------------------------------
   call StopHandler('bzrduc')
!  ---------------------------------------------------------------------
   end subroutine bzrduc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function inbz(wvk,rsdir,nrsdir,nplane) result(nbz)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nrsdir, nplane
   integer (kind=IntKind) :: nbz
   integer (kind=IntKind) :: n
   integer (kind=IntKind), parameter :: yes = 1
   integer (kind=IntKind), parameter :: no = 0
!
   real (kind=RealKind), intent(in) :: wvk(3),rsdir(4,nrsdir)
   real (kind=RealKind) :: project
   real (kind=RealKind), parameter :: eps = TEN2m6
!
!  ===================================================================
!  is wvk in the 1st brillouin zone ?
!  check whether wvk lies inside all the planes that define the 1bz.
!  ===================================================================
   nbz = no
   do n = 1,nplane
      project=(rsdir(1,n)*wvk(1)+rsdir(2,n)*wvk(2)+rsdir(3,n)*wvk(3))/rsdir(4,n)
!     ================================================================
!     check wvk is outside the bz
!     ================================================================
      if (abs(project) > HALF + eps) then
         return
      endif
   enddo
   nbz = yes
   end function inbz
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printBZone(Rot3D)
!  ===================================================================
   implicit none
!
   logical, intent(in), optional :: Rot3D
   logical :: print_Rot3D
!
   integer (kind=IntKind) :: i, k
!
   if (.not.Initialized) then
      call WarningHandler('printBZone','BZoneModule is not initialized')
      return
   endif
!
   if (present(Rot3D)) then
      print_Rot3D = Rot3D
   else
      print_Rot3D = .false.
   endif
!
   write(6,'(/,a)')'*******************************'
   write(6,'( a )')'* Output from printBZone      *'
   write(6,'(a,/)')'*******************************'
   write(6,'(''Number of K Meshes = '',i5)')NumMeshs
   do i=1,NumMeshs
      write(6,'(/,''================================================='')')
      write(6,'(''Mesh index         = '',i5)')i
      write(6,'(''Number of K points = '',i5)')NumKs(i)
      write(6,'(''================================================='')')
      write(6,'(''     kx           ky           kz           wght '')')
      write(6,'(''-------------------------------------------------'')')
      do k=1,NumKs(i)
         write(6,'(4(f10.5,3x))')KxPoint(k,i),KyPoint(k,i),KzPoint(k,i), &
                                 Kweight(k,i)
      enddo
   enddo
   write(6,'(''================================================='')')
!
   if (print_Rot3D) then
      write(6,'(/,''Number of Rotations = '',i5)')NumRotations
      write(6,'(''Rotation Matrix = '')')
      write(6,'(80(''=''))')
      do i=1,NumRotations,3
         k = i
         if (i+1 >= NumRotations) then
            exit
         endif
         write(6,'(/,''      ['',f5.2,x,f5.2,x,f5.2,'']        ['',         &
            &      f5.2,x,f5.2,x,f5.2,'']        ['',f5.2,x,f5.2,x,f5.2,'']'')') &
            Rotation(i,1,1),Rotation(i,1,2),Rotation(i,1,3),              &
            Rotation(i+1,1,1),Rotation(i+1,1,2),Rotation(i+1,1,3),        &
            Rotation(i+2,1,1),Rotation(i+2,1,2),Rotation(i+2,1,3)
         write(6,'(''r('',i2,'')=['',f5.2,x,f5.2,x,f5.2,''], r('',i2,'')=['',    &
            &      f5.2,x,f5.2,x,f5.2,''], r('',i2,'')=['',f5.2,x,f5.2,x,f5.2,   &
            &      '']'')') &
            i,Rotation(i,2,1),Rotation(i,2,2),Rotation(i,2,3),            &
            i+1,Rotation(i+1,2,1),Rotation(i+1,2,2),Rotation(i+1,2,3),    &
            i+2,Rotation(i+2,2,1),Rotation(i+2,2,2),Rotation(i+2,2,3)
         write(6,'(''      ['',f5.2,x,f5.2,x,f5.2,'']        ['',         &
            &      f5.2,x,f5.2,x,f5.2,'']        ['',f5.2,x,f5.2,x,f5.2,'']'')') &
            Rotation(i,3,1),Rotation(i,3,2),Rotation(i,3,3),              &
            Rotation(i+1,3,1),Rotation(i+1,3,2),Rotation(i+1,3,3),        &
            Rotation(i+2,3,1),Rotation(i+2,3,2),Rotation(i+2,3,3)
      enddo
      if (NumRotations-k == 0) then
         write(6,'(/,''      ['',f5.2,x,f5.2,x,f5.2,'']'')')              &
            Rotation(k,1,1),Rotation(k,1,2),Rotation(k,1,3)
         write(6,'(''r('',i2,'')=['',f5.2,x,f5.2,x,f5.2,'']'')')          &
            k,Rotation(k,2,1),Rotation(k,2,2),Rotation(k,2,3)
         write(6,'(''      ['',f5.2,x,f5.2,x,f5.2,'']'')')                &
            Rotation(k,3,1),Rotation(k,3,2),Rotation(k,3,3)
      else if (NumRotations-k == 1) then
         write(6,'(/,''      ['',f5.2,x,f5.2,x,f5.2,'']        ['',       &
            &      f5.2,x,f5.2,x,f5.2,'']'')')                            &
            Rotation(k,1,1),Rotation(k,1,2),Rotation(k,1,3),              &
            Rotation(k+1,1,1),Rotation(k+1,1,2),Rotation(k+1,1,3)
         write(6,'(''r('',i2,'')=['',f5.2,x,f5.2,x,f5.2,''], r('',i2,'')=['',  &
            &      f5.2,x,f5.2,x,f5.2,'']'')')                            &
            k,Rotation(k,2,1),Rotation(k,2,2),Rotation(k,2,3),            &
            k+1,Rotation(k+1,2,1),Rotation(k+1,2,2),Rotation(k+1,2,3)
         write(6,'(''      ['',f5.2,x,f5.2,x,f5.2,'']        ['',         &
            &      f5.2,x,f5.2,x,f5.2,'']'')')                            &
            Rotation(k,3,1),Rotation(k,3,2),Rotation(k,3,3),              &
            Rotation(k+1,3,1),Rotation(k+1,3,2),Rotation(k+1,3,3)
      endif
      write(6,'(80(''-''))')
   endif
!
   end subroutine printBZone
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumKMeshs() result(m)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: m
!
   if (.not.Initialized) then
      call ErrorHandler('getNumKMeshs','BZoneModule is not initialized')
   endif
   m = NumMeshs
   end function getNumKMeshs
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWeight0(k) result(w)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: k
   real (kind=RealKind) :: w
!
   if (.not.Initialized) then
      call ErrorHandler('getWeight','BZoneModule is not initialized')
   else if (k < 1 .or. k > NumKs(1)) then
      call ErrorHandler('getWeight','invalid k-point index',k)
   endif
   w = Kweight(k,1)
   end function getWeight0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWeight1(k,m) result(w)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: m
   integer (kind=IntKind), intent(in) :: k
   real (kind=RealKind) :: w
!
   if (.not.Initialized) then
      call ErrorHandler('getWeight','BZoneModule is not initialized')
   else if (m < 1 .or. m > NumMeshs) then
      call ErrorHandler('getWeight','invalid mesh index',m)
   else if (k < 1 .or. k > NumKs(m)) then
      call ErrorHandler('getWeight','invalid k-point index',k)
   endif
   w = Kweight(k,m)
   end function getWeight1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAllWeights0() result(w)
!  ===================================================================
   implicit none
   real (kind=RealKind), pointer :: w(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getAllWeights','BZoneModule is not initialized')
   endif
   w => Kweight(1:NumKs(1),1)
   end function getAllWeights0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAllWeights1(m) result(w)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: m
   real (kind=RealKind), pointer :: w(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getAllWeights','BZoneModule is not initialized')
   else if (m < 1 .or. m > NumMeshs) then
      call ErrorHandler('getAllWeights','invalid mesh index',m)
   endif
   w => Kweight(1:NumKs(m),m)
   end function getAllWeights1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWeightSum0() result(ws)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: ws
!
   if (.not.Initialized) then
      call ErrorHandler('getWeightSum','BZoneModule is not initialized')
   endif
   ws = KweightSum(1)
   end function getWeightSum0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWeightSum1(m) result(ws)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: m
   real (kind=RealKind) :: ws
!
   if (.not.Initialized) then
      call ErrorHandler('getWeightSum','BZoneModule is not initialized')
   else if (m < 1 .or. m > NumMeshs) then
      call ErrorHandler('getWeightSum','invalid mesh index',m)
   endif
   ws = KweightSum(m)
   end function getWeightSum1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBZoneVolume() result(v)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: v
!
   v = vol_bz
!
   end function getBZoneVolume
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine integBZoneR0(f,n,r)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: ik, i
!
   real (kind=RealKind), intent(in) :: f(:,:)
   real (kind=RealKind), intent(out) :: r(n)
!
   if (.not.Initialized) then
      call ErrorHandler('integBZone','BZoneModule is not initialized')
   else if (size(f,1) < n) then
      call ErrorHandler('integBZone','1st Dim of f is out of bound',size(f,1))
   else if (size(f,2) < NumKs(1)) then
      call ErrorHandler('integBZone','2nd Dim of f is out of bound',size(f,2))
   endif
!
   r(1:n)=ZERO
   do ik=1,NumKs(1)
      do i=1,n
         r(i) = r(i)+f(i,ik)*Kweight(ik,1)
      enddo
   enddo
!
   end subroutine integBZoneR0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine integBZoneR1(f,n,im,r)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n, im
   integer (kind=IntKind) :: ik, i
!
   real (kind=RealKind), intent(in) :: f(:,:)
   real (kind=RealKind), intent(out) :: r(n)
!
   if (.not.Initialized) then
      call ErrorHandler('integBZone','BZoneModule is not initialized')
   else if (im < 1 .or. im > NumMeshs) then
      call ErrorHandler('integBZone','invalid mesh index',im)
   else if (size(f,1) < n) then
      call ErrorHandler('integBZone','1st Dim of f is out of bound',size(f,1))
   else if (size(f,2) < NumKs(im)) then
      call ErrorHandler('integBZone','2nd Dim of f is out of bound',size(f,2))
   endif
!
   r(1:n)=ZERO
   do ik=1,NumKs(im)
      do i=1,n
         r(i) = r(i)+f(i,ik)*Kweight(ik,im)
      enddo
   enddo
!
   end subroutine integBZoneR1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine integBZoneC0(f,n,c)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: ik, i
!
   complex (kind=CmplxKind), intent(in) :: f(:,:)
   complex (kind=CmplxKind), intent(out) :: c(n)
!
   if (.not.Initialized) then
      call ErrorHandler('integBZone','BZoneModule is not initialized')
   else if (size(f,1) < n) then
      call ErrorHandler('integBZone','1st Dim of f is out of bound',size(f,1))
   else if (size(f,2) < NumKs(1)) then
      call ErrorHandler('integBZone','2nd Dim of f is out of bound',size(f,2))
   endif
!
   c(1:n)=CZERO
   do ik=1,NumKs(1)
      do i=1,n
         c(i) = c(i)+f(i,ik)*Kweight(ik,1)
      enddo
   enddo
!
   end subroutine integBZoneC0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine integBZoneC1(f,n,im,c)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n, im
   integer (kind=IntKind) :: ik, i
!
   complex (kind=CmplxKind), intent(in) :: f(:,:)
   complex (kind=CmplxKind), intent(out) :: c(n)
!
   if (.not.Initialized) then
      call ErrorHandler('integBZone','BZoneModule is not initialized')
   else if (im < 1 .or. im > NumMeshs) then
      call ErrorHandler('integBZone','invalid mesh index',im)
   else if (size(f,1) < n) then
      call ErrorHandler('integBZone','1st Dim of f is out of bound',size(f,1))
   else if (size(f,2) < NumKs(im)) then
      call ErrorHandler('integBZone','2nd Dim of f is out of bound',size(f,2))
   endif
!
   c(1:n)=CZERO
   do ik=1,NumKs(im)
      do i=1,n
         c(i) = c(i)+f(i,ik)*Kweight(ik,im)
      enddo
   enddo
!
   end subroutine integBZoneC1
!  ===================================================================
end module BZoneModule
