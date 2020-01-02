module MediumHostModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MathParamModule, only : TEN2m8, ONE, ZERO
!
public :: initMediumHost,    &
          endMediumHost,     &
          printMediumHost,   &
          getNumSites,       &
          getSitePosition,   &
          getLmaxKKR,        &
          getLmaxPhi,        &
          getMediumLattice,  &
          getGlobalSiteIndex,&
          getLocalNumSites,  &
          getNumSpecies,     &
          getSpeciesContent, &
          getSpeciesAtomicNumber
!
private
   integer (kind=IntKind) :: NumSites       ! number of atomic sites in the medium unit cell
   integer (kind=IntKind) :: MyPEinGroup, NumPEsInGroup, GroupID
   integer (kind=IntKind) :: lmax_kkr_max, lmax_phi_max
   integer (kind=IntKind) :: max_species, max_local_sites
!
   integer (kind=IntKind), allocatable :: NumLocalSites(:) ! number of medium sites on a local process
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: global_index(:,:)
   integer (kind=IntKind), allocatable :: num_species(:)
   integer (kind=IntKind), allocatable :: species_Z(:,:)
!
   real (kind=RealKind) :: bravais_lattice(3,3)
   real (kind=RealKind), allocatable :: site_position(:,:)
   real (kind=RealKind), allocatable :: species_content(:,:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMediumHost(tbl_id)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getMyPEinGroup, getNumPEsInGroup
!
   use ChemElementModule, only :  MaxLenOfAtomName, getZtot, getZval, getName
!
   use InputModule, only : getKeyValue, getTableIndex, getKeyIndexValue
   use InputModule, only : isKeyExisting
!
   use ScfDataModule, only : isKKRCPA, isEmbeddedCluster
!
   use SystemModule, only : getNumAtoms, getLmaxKKR, getLmaxPhi,      &
                            getBravaisLattice,                        &
                            getAtomPosition, getNumAlloyElements,     &
                            getAlloyElementContent, getAlloyElementName
!
   use Atom2ProcModule, only : getLocalNumAtoms, getMaxLocalNumAtoms, &
                               getGlobalIndex
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: tbl_id
   integer (kind=IntKind) :: i, ia, pe, ig
!
   character (len=MaxLenOfAtomName) :: aname
!
   real (kind=RealKind) :: c
!
   if (isKKRCPA()) then
      NumSites = getNumAtoms()
      GroupID = getGroupID('Unit Cell')
      MyPEinGroup = getMyPEinGroup(GroupID)
      NumPEsInGroup = getNumPEsInGroup(GroupID)
      max_local_sites = getMaxLocalNumAtoms()
      lmax_kkr_max = getLmaxKKR()
      lmax_phi_max = getLmaxPhi()
      bravais_lattice = getBravaisLattice()
   else if (isEmbeddedCluster()) then
      call ErrorHandler('initMediumHost','Not yet implemented')
   else
      call ErrorHandler('initMediumHost','Invalid medium type')
   endif
!
   allocate(NumLocalSites(NumPEsInGroup))
   allocate(global_index(max_local_sites,NumPEsInGroup))
   allocate(num_species(NumSites), lmax_kkr(NumSites), lmax_phi(NumSites))
   allocate(site_position(3,NumSites))
!
   if (isKKRCPA()) then
      global_index = -1
      do pe = 0, NumPEsInGroup-1
        NumLocalSites(pe+1) = getLocalNumAtoms(pe)
        do i = 1, NumLocalSites(pe+1)
           global_index(i,pe+1) = getGlobalIndex(i,pe)
        enddo
      enddo
      do ig = 1, NumSites
         num_species(ig) = getNumAlloyElements(ig)
         lmax_kkr(ig) = getLmaxKKR(ig)
         lmax_phi(ig) = getLmaxPhi(ig)
         site_position(1:3,ig) = getAtomPosition(ig)
      enddo
   else if (isEmbeddedCluster()) then
      call ErrorHandler('initMediumHost','Not yet implemented')
   else
      call ErrorHandler('initMediumHost','Invalid medium type')
   endif
!
   max_species = 1
   do ig = 1, NumSites
      max_species = max(max_species,num_species(ig)) 
   enddo
!
   allocate(species_content(max_species,NumSites))
   allocate(species_Z(max_species,NumSites))
!
   species_content = ZERO
   if (isKKRCPA()) then
      do ig = 1, NumSites
         c = ZERO
         do ia = 1, num_species(ig)
            species_content(ia,ig) = getAlloyElementContent(ig,ia)
            aname = getAlloyElementName(ig,ia)
            species_Z(ia,ig) = getZtot(aname)
            c = c + species_content(ia,ig)
         enddo
         if (abs(c - ONE) > TEN2m8) then
            call WarningHandler('initMediumHost','Total species content <> 1',c)
            call ErrorHandler('initMediumHost','Site index',ig)
         endif
      enddo
   else if (isEmbeddedCluster()) then
      call ErrorHandler('initMediumHost','Not yet implemented')
   else
      call ErrorHandler('initMediumHost','Invalid medium type')
   endif
!
   end subroutine initMediumHost
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endMediumHost()
!  ===================================================================
   implicit none
!
   deallocate(NumLocalSites, lmax_kkr, lmax_phi, site_position, global_index)
   deallocate(num_species, species_content, species_Z)
!
   end subroutine endMediumHost
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumSites() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumSites
!
   end function getNumSites
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSitePosition(i) result(posi)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
!
   real (kind=RealKind) :: posi(3)
!
   if (i < 1 .or. i > NumSites) then
      call ErrorHandler('getSitePosition','Medium atom index is out of range', &
                        i,NumSites)
   endif
!
   posi(1:3) = site_position(1:3,i)
!
   end function getSitePosition
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxKKR(i) result(lmax)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: i
   integer (kind=IntKind) :: lmax
!
   if (present(i)) then
      if (i < 1 .or. i > NumSites) then
         call ErrorHandler('getLmaxKKR','Medium atom index is out of range', &
                           i,NumSites)
      endif
      lmax = lmax_kkr(i)
   else
      lmax = lmax_kkr_max
   endif
!
   end function getLmaxKKR
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxPhi(i) result(lmax)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: i
   integer (kind=IntKind) :: lmax
!
   if (present(i)) then
      if (i < 1 .or. i > NumSites) then
         call ErrorHandler('getLmaxPhi','Medium atom index is out of range', &
                           i,NumSites)
      endif
      lmax = lmax_phi(i)
   else
      lmax = lmax_phi_max
   endif
!
   end function getLmaxPhi
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMediumLattice() result(bra)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: bra(3,3)
!
   bra = bravais_lattice
!
   end function getMediumLattice
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalSiteIndex(i,pe) result(j)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind), intent(in), optional :: pe
   integer (kind=IntKind) :: n, j
!
   if (present(pe)) then
      if (pe < 0 .or. pe > NumPEsInGroup-1) then
         call ErrorHandler('getGlobalSiteIndex','The PE index is out of range', &
                           pe,NumPEsInGroup-1)
      endif
      n = pe + 1
   else 
      n = MyPEinGroup + 1
   endif
!
   if (i < 1 .or. i > NumLocalSites(n)) then
      call ErrorHandler('getGlobalSiteIndex','Medium atom local index is out of range', &
                        i,NumLocalSites(n))
   endif
!
   j = global_index(i,n)
!
   end function getGlobalSiteIndex
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalNumSites(pe) result(lna)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: pe
   integer (kind=IntKind) :: n, lna
!
   if (present(pe)) then
      if (pe < 0 .or. pe > NumPEsInGroup-1) then
         call ErrorHandler('getLocalNumSites','The PE index is out of range', &
                           pe,NumPEsInGroup-1)
      endif
      n = pe + 1
   else 
      n = MyPEinGroup + 1
   endif
!
   lna = NumLocalSites(n)
!
   end function getLocalNumSites
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumSpecies(ig) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ig
   integer (kind=IntKind) :: n
!
   if (ig < 1 .or. ig > NumSites) then
      call ErrorHandler('getNumSpecies','Medium atom index is out of range',ig)
   endif
!
   n = num_species(ig)
!
   end function getNumSpecies
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSpeciesContent(ic,ig) result(c)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ic, ig
!
   real (kind=RealKind) :: c
!
   if (ig < 1 .or. ig > NumSites) then
      call ErrorHandler('getSpeciesContent','Medium atom index is out of range',ig)
   else if (ic < 1 .or. ic > num_species(ig)) then
      call ErrorHandler('getSpeciesContent','Species index is out of range',ic)
   endif
!
   c = species_content(ic,ig)
!
   end function getSpeciesContent
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSpeciesAtomicNumber(ic,ig) result(z)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ic, ig
   integer (kind=IntKind) :: z
!
   if (ig < 1 .or. ig > NumSites) then
      call ErrorHandler('getSpeciesAtomicNumber','Medium atom index is out of range',ig)
   else if (ic < 1 .or. ic > num_species(ig)) then
      call ErrorHandler('getSpeciesAtomicNumber','Species index is out of range',ic)
   endif
!
   z = species_Z(ic,ig)
!
   end function getSpeciesAtomicNumber
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printMediumHost(pe)
!  ===================================================================
   use MPPModule, only : MyPE, startRoundTurn, finishMyTurn
   use MPPModule, only : syncAllPEs
   use ChemElementModule, only :  getName
   implicit none
!
   integer (kind=intKind), intent(in), optional :: pe
   integer (kind=intKind) :: print_pe
!
   if (present(pe)) then
      print_pe = pe
   else
      print_pe = 0
   endif
!
   if (print_pe < 0) then
      call startRoundTurn()
      call printData()
      call finishMyTurn()
   else 
      if (print_pe == MyPE) then
         call printData()
      endif
      call syncAllPEs()
   endif
!
   end subroutine printMediumHost
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printData()
!  ===================================================================
   use MPPModule, only : MyPE
   use ChemElementModule, only :  getName
   implicit none
!
   integer (kind=intKind) :: ig, ia, z
!
   write(6,'(/,a)')'======= Print Medium Host Data ======='
   write(6,'(a,i5)')'Total number of atomic sites: ',NumSites
   write(6,'(a,2i5)')'My process ID in the group and whole world: ',   &
                     MyPEinGroup,MyPE
   do ig = 1, NumSites
      write(6,'(a,i5,4x,3f10.5)')'Atomic site index and position: ',   &
                                 ig,site_position(1:3,ig)
      write(6,'(a,i5)')'Number of species on the site: ',num_species(ig)
      do ia = 1, num_species(ig)
         write(6,'(4x,a,a,4x,a,f10.5)')'Species name   : ',            &
                                       getName(species_Z(ia,ig)),      &
                                       'Species content: ',            &
                                       species_content(ia,ig)
      enddo
   enddo
   write(6,'(a,/)')'========================================'
!
   end subroutine printData
!  ===================================================================
end module MediumHostModule
