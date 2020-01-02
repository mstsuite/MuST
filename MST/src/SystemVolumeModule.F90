module SystemVolumeModule
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ZERO
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initSystemVolume,           &
          endSystemVolume,            &
          printSystemVolume,          &
          setSystemVolumeMT,          &
          updateSystemVolume,         &
          getSystemVolume,            &
          getTotalVPVolume,           &
          getTotalInscribedVolume,    &  ! Note: for ASA case, this still returns Inscribed volume
          getTotalMTVolume,           &  ! Note: for ASA case, this still returns muffin-tin volume
          getAtomicVPVolume,          &
          getAtomicInscribedVolume,   &  ! Note: for ASA case, this still returns Inscribed volume
          getAtomicMTVolume,          &  ! Note: for ASA case, this still returns muffin-tin volume
          getTotalInterstitialVolume, &
          getTotalInterstitialMTVolume
!
   interface getAtomicVPVolume
      module procedure getAtomicVPVolume_one, getAtomicVPVolume_all
   end interface getAtomicVPVolume
!
   interface getAtomicInscribedVolume
      module procedure getAtomicInscribedVolume_one, getAtomicInscribedVolume_all
   end interface getAtomicInscribedVolume
!
private
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: LocalNumAtoms
!
   real (kind=RealKind) :: SystemVolume
   real (kind=RealKind) :: TotalMuffintinVolume
   real (kind=RealKind) :: TotalInscribedVolume
   real (kind=RealKind) :: TotalVoronoiPolyhedraVolume
   real (kind=RealKind) :: TotalInterstitialVolume
   real (kind=RealKind) :: TotalInterstitialMTVolume
   real (kind=RealKind) :: Bravais(3,3)
   real (kind=RealKind), allocatable, target :: AtomicVoronoiPolyhedraVolume(:)
   real (kind=RealKind), allocatable, target :: AtomicInscribedVolume(:)
   real (kind=RealKind), allocatable :: AtomicMuffintinVolume(:)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSystemVolume()
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup, &
                               GlobalSumInGroup
!
   use MathParamModule, only : THIRD, PI4
!
   use SystemModule, only : getNumAtoms, getBravaisLattice, getAtomPosition
   use SystemModule, only : getRadicalPlaneRatio
!
   use Atom2ProcModule, only : getGlobalIndex
   use Atom2ProcModule, only : getLocalIndex, getAtom2ProcInGroup
   use Atom2ProcModule, only : printAtom2ProcTable, getLocalNumAtoms
   use Atom2ProcModule, only : getMaxLocalNumAtoms
!  
   use PolyhedraModule, only : initPolyhedra
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : getInscrSphRadius, getVolume, getBoxVolume
!
   implicit none
!
   integer (kind=IntKind) :: i, ig, maxnum, id, MyPEinGroup, NumPEsInGroup
!
   real (kind=RealKind), allocatable :: memtemp(:,:)
   real (kind=RealKind), allocatable :: atompos(:,:)
   real (kind=RealKind), allocatable :: radplane(:)
   real (kind=RealKind) :: msgbuf(3)
!
   id = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(id)
   MyPEinGroup = getMyPEinGroup(id)
!
   GlobalNumAtoms = getNumAtoms()
   Bravais(1:3,1:3) = getBravaisLattice()
!
!  -------------------------------------------------------------------
   LocalNumAtoms=getLocalNumAtoms(MyPEinGroup)
!  -------------------------------------------------------------------
   call initPolyhedra(LocalNumAtoms,Bravais,'None',0)
!  -------------------------------------------------------------------
!
   allocate(atompos(3,GlobalNumAtoms))
   allocate(radplane(GlobalNumAtoms))
   do i=1,GlobalNumAtoms
      atompos(1:3,i)=getAtomPosition(i)
      radplane(i)=getRadicalPlaneRatio(i)
   enddo
!
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i,MyPEinGroup)
!     ----------------------------------------------------------------
      call genPolyhedron(i,ig,GlobalNumAtoms,atompos,radplane)
!     ----------------------------------------------------------------
   enddo
   deallocate(atompos,radplane)
!
!  SystemVolume = ( Bravais(2,1)*Bravais(3,2)-Bravais(3,1)*Bravais(2,2)) &
!                  *Bravais(1,3)+                                        &
!                 ( Bravais(3,1)*Bravais(1,2)-Bravais(1,1)*Bravais(3,2)) &
!                  *Bravais(2,3)+                                        &
!                 ( Bravais(1,1)*Bravais(2,2)-Bravais(2,1)*Bravais(1,2)) &
!                  *Bravais(3,3)
!  SystemVolume=abs(SystemVolume) 
   SystemVolume = getBoxVolume()
!
   msgbuf(1:3) = ZERO
   do i=1,LocalNumAtoms
       msgbuf(1) = msgbuf(1) + PI4*THIRD*(getInscrSphRadius(i))**3
       msgbuf(2) = msgbuf(2) + getVolume(i)
       msgbuf(3) = msgbuf(3) + getVolume(i)                           &
                  - PI4*THIRD*(getInscrSphRadius(i))**3
   enddo
!
!  -------------------------------------------------------------------
   call GlobalSumInGroup(id,msgbuf,3)
!  -------------------------------------------------------------------
   TotalInscribedVolume = msgbuf(1)
   TotalMuffintinVolume = TotalInscribedVolume
   TotalVoronoiPolyhedraVolume = msgbuf(2)
   TotalInterstitialVolume = msgbuf(3)
   TotalInterstitialMTVolume = TotalInterstitialVolume
!
   allocate(AtomicVoronoiPolyhedraVolume(GlobalNumAtoms))
   allocate(AtomicMuffintinVolume(GlobalNumAtoms))
   allocate(AtomicInscribedVolume(GlobalNumAtoms))
!
   allocate(memtemp(2,GlobalNumAtoms))
   memtemp(1:2,1:GlobalNumAtoms) = ZERO
!
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i,MyPEinGroup)
      memtemp(1,ig) = getVolume(i)
      memtemp(2,ig) = PI4*THIRD*(getInscrSphRadius(i))**3
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(id,memtemp,2,GlobalNumAtoms)
!  -------------------------------------------------------------------
   do ig=1,GlobalNumAtoms
      AtomicVoronoiPolyhedraVolume(ig) = memtemp(1,ig)
      AtomicInscribedVolume(ig) = memtemp(2,ig)
      AtomicMuffintinVolume(ig) = AtomicInscribedVolume(ig)
   enddo
!
   deallocate(memtemp)
!
   end subroutine initSystemVolume
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSystemVolume()
!  ===================================================================
   use Atom2ProcModule, only : endAtom2Proc
   use PolyhedraModule, only : endPolyhedra
!
   implicit none
!
!  -------------------------------------------------------------------
   call endPolyhedra()
!  -------------------------------------------------------------------
!
   deallocate(AtomicVoronoiPolyhedraVolume)
   deallocate(AtomicInscribedVolume)
   deallocate(AtomicMuffintinVolume)
!
   end subroutine endSystemVolume
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSystemVolume() result(v)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: v
!  
   v = SystemVolume
!  
   end function getSystemVolume
!  ===================================================================
!  
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTotalVPVolume() result(v)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: v
!  
   v = TotalVoronoiPolyhedraVolume
!  
   end function getTotalVPVolume        
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTotalInscribedVolume() result(v)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: v
!  
   v = TotalInscribedVolume
!  
   end function getTotalInscribedVolume   
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTotalMTVolume() result(v)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: v
!  
   v = TotalMuffintinVolume
!  
   end function getTotalMTVolume 
!  ===================================================================
!  
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTotalInterstitialVolume() result(v)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: v
!
   v = TotalInterstitialVolume
!
   end function getTotalInterstitialVolume
!  ===================================================================
!  
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTotalInterstitialMTVolume() result(v)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: v
!
   v = TotalInterstitialMTVolume
!
   end function getTotalInterstitialMTVolume
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicVPVolume_one(ig) result(v)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ig
   real (kind=RealKind) :: v
!
   if (ig<1 .or. ig>GlobalNumAtoms) then
      call ErrorHandler('getAtomicVPVolume','Invalid atom index',ig)
   endif
   v = AtomicVoronoiPolyhedraVolume(ig)
!
   end function getAtomicVPVolume_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicVPVolume_all() result(v)
!  ===================================================================
   implicit none
   real (kind=RealKind), pointer :: v(:)
!
   v => AtomicVoronoiPolyhedraVolume
!
   end function getAtomicVPVolume_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicInscribedVolume_one(ig) result(v)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ig
   real (kind=RealKind) :: v
!
   if (ig<1 .or. ig>GlobalNumAtoms) then
      call ErrorHandler('getAtomicInscribedVolume','Invalid atom index',ig)
   endif
   v = AtomicInscribedVolume(ig)
!
   end function getAtomicInscribedVolume_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicInscribedVolume_all() result(v)
!  ===================================================================
   implicit none
   real (kind=RealKind), pointer :: v(:)
!
   v => AtomicInscribedVolume
!
   end function getAtomicInscribedVolume_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicMTVolume(ig) result(v)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ig
   real (kind=RealKind) :: v
!
   if (ig<1 .or. ig>GlobalNumAtoms) then
      call ErrorHandler('getAtomicMTVolume','Invalid atom index',ig)
   endif
   v = AtomicMuffintinVolume(ig)
!
   end function getAtomicMTVolume
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSystemVolumeMT(id,rmt)
!  ===================================================================
   use Atom2ProcModule, only : getGlobalIndex
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   use MathParamModule, only : THIRD, PI4
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   real (kind=RealKind), intent(in) :: rmt
!
   integer (kind=IntKind) :: i, ig, MyPEinGroup
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('setSystemVolumeMT','Invalid atom index',ig)
   endif
!
   i = getGroupID('Unit Cell')
   MyPEinGroup = getMyPEinGroup(i)
   ig=getGlobalIndex(id,MyPEinGroup)
!
   AtomicMuffintinVolume(ig)=PI4*THIRD*(rmt)**3
!
   end subroutine setSystemVolumeMT
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateSystemVolume()
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup,          &
                                getMyPEinGroup, GlobalSumInGroup
!
   use Atom2ProcModule, only : getLocalIndex, getAtom2ProcInGroup
   use Atom2ProcModule, only : getMaxLocalNumAtoms
   use Atom2ProcModule, only : getGlobalIndex
   use GroupCommModule, only : getGroupID, getMyPEinGroup
!
   implicit none
!
   integer (kind=IntKind) :: i, ig, id
   integer (kind=IntKind) :: maxnum, MyPEinGroup, NumPEsInGroup
!
   real (kind=RealKind), allocatable :: memtemp(:)
!
   id = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(id)
   MyPEinGroup = getMyPEinGroup(id)
   maxnum = getMaxLocalNumAtoms()
!
   allocate(memtemp(GlobalNumAtoms))
   memtemp(1:GlobalNumAtoms) = ZERO
!
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i,MyPEinGroup)
      memtemp(ig) = AtomicMuffintinVolume(ig)
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(id,memtemp,GlobalNumAtoms)
!  -------------------------------------------------------------------
   do ig=1,GlobalNumAtoms
      AtomicMuffintinVolume(ig) = memtemp(ig)
   enddo
!
   TotalInterstitialMTVolume = ZERO
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i,MyPEinGroup)
      TotalInterstitialMTVolume = AtomicVoronoiPolyhedraVolume(ig) - &
                                  AtomicMuffintinVolume(ig)
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(id,TotalInterstitialMTVolume)
!  -------------------------------------------------------------------
   TotalMuffintinVolume = TotalVoronoiPolyhedraVolume - &
                          TotalInterstitialMTVolume
!  -------------------------------------------------------------------
   deallocate(memtemp)
!
   end subroutine updateSystemVolume
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printSystemVolume()
!  ===================================================================
   implicit none
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')'*************************************'
   write(6,'( 24x,a )')'*   Output from printSystemVolume   *'
   write(6,'(24x,a,/)')'*************************************'
   write(6,'(/,80(''=''))')
   write(6,'(''Total System Volume:            '',f19.8)')SystemVolume
   write(6,'(''Total Voronoi Polyhedra Volume: '',f19.8)')              &
                                           TotalVoronoiPolyhedraVolume
   write(6,'(''Total Inscribed Sphere Volume : '',f19.8)')TotalInscribedVolume
   write(6,'(''Total Interstitial Volume     : '',f19.8)')TotalInterstitialVolume
   write(6,'(''Total Muffintin Sphere Volume : '',f19.8)')TotalMuffintinVolume
   write(6,'(''Total Interstitial MT Volume  : '',f19.8)')TotalInterstitialMTVolume
   write(6,'(80(''=''))')
!
   end subroutine printSystemVolume
!  ===================================================================
end module SystemVolumeModule
