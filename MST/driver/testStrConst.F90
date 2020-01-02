program testStrConst
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
   use ScfDataModule, only : isReadEmesh, getEmeshFileName
   use ScfDataModule, only : isReadKmesh, getKmeshFileName
   use ScfDataModule, only : NumEs, ContourType, eGridType
   use ScfDataModule, only : NumKMeshs, kGenScheme, Kdiv, Symmetrize
   use ScfDataModule, only : ErBottom, ErTop, EiBottom, EiTop, Temperature
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType, &
                                   isFullPotential
!
   use ContourModule, only : initContour, endContour, getNumEs, getEPoint, &
                             setupContour
!
   use BZoneModule, only : initBZone, printBZone, endBZone, getNumKs, getAllKPoints
!
   use TimerModule, only : initTimer, getTime
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, SQRTm1, TEN2m8
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use MatrixModule, only : computeAStarT
!
   use StrConstModule, only : initStrConst, endStrConst
   use StrConstModule, only : getStrConstMatrix
!
   use MPPModule, only : initMPP, endMPP
!
   use WriteMatrixModule,  only : writeMatrix
!
   implicit   none
!
   character (len=4) :: istop = 'none'
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: lmax_kkr
   integer (kind=IntKind) :: lmax_phi
   integer (kind=IntKind) :: kmax_kkr
   integer (kind=IntKind) :: kmax_phi
   integer (kind=IntKind) :: ne = 1
   integer (kind=IntKind) :: nk = 1
!
   integer (kind=IntKind) :: i, j, n, kl, klp
   integer (kind=IntKind) :: ie
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: fstatus
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
!
   real (kind=RealKind), parameter :: tol = TEN2m8
   real (kind=RealKind) :: t0, kfac
   real (kind=RealKind), pointer :: Bravais(:,:)
!
   real (kind=RealKind), pointer :: kvec(:,:)
   real (kind=RealKind) :: kin(3)
!
! *********************************************************************
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
!
   integer (kind=IntKind), parameter :: max_epts = 10
!
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind), pointer :: e_mesh(:)
   complex (kind=CmplxKind), pointer :: strcon_matrix(:,:)
   complex (kind=CmplxKind), allocatable :: scm(:,:)
   complex (kind=CmplxKind), allocatable :: scmST(:,:)
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
!  ===================================================================
!  Initialize the contour in energy complex plane to find the total
!  number of energy mesh needed
!  ===================================================================
   if (isReadEmesh()) then
!     ----------------------------------------------------------------
      call initContour(getEmeshFileName(), 'none', -1)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call initContour( ContourType, eGridType, NumEs, Temperature, 'none', -1)
!     ----------------------------------------------------------------
      if (ErTop <= ErBottom) then
         ErTop = ErBottom + ONE
      endif
!     ----------------------------------------------------------------
      call setupContour( ErBottom, ErTop, EiBottom, EiTop )
!     ----------------------------------------------------------------
   endif
!  
!  -------------------------------------------------------------------
   ne = getNumEs()
   e_mesh => getEPoint()
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
!  -------------------------------------------------------------------
   fstatus = getKeyValue(info_id,'Default Lmax-T matrix',lmax_kkr)
   fstatus = getKeyValue(info_id,'Default Lmax-Wave Func',lmax_phi)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(lmax_kkr+lmax_phi)
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
   kmax_kkr=(lmax_kkr+1)**2
   kmax_phi=(lmax_phi+1)**2
   allocate(scm(kmax_phi,kmax_phi),scmST(kmax_phi,kmax_phi))
!
!  --------------------------------------------------------------------
   call initTimer()
   call initStrConst(lmax_phi,NumAtoms,AtomPosition,Bravais,istop,iprint)
!  --------------------------------------------------------------------
   do ie=1,ne
      do k=1,nk
!        =============================================================
!        set energy &  electron momentum: ............................
!        =============================================================
         energy=e_mesh(ie)
         kin(1)=kvec(1,k)
         kin(2)=kvec(2,k)
         kin(3)=kvec(3,k)
         write(6,'(/,'' K-point: kx, ky, kz ='',3f15.8)')          &
               kin(1),kin(2),kin(3)
         write(6,'(  '' Energy ='',2d15.8)')energy
!        =============================================================
         kappa=sqrt(energy)
         do n=1,NumAtoms*NumAtoms
            i = mod(n-1,NumAtoms)+1
            j = (n-1)/NumAtoms+1
            t0=getTime()
!           ----------------------------------------------------------
            strcon_matrix => getStrConstMatrix(kin,kappa,i,j,lmax_phi,lmax_phi)
            scm = strcon_matrix
            strcon_matrix => getStrConstMatrix(-kin,kappa,j,i,lmax_phi,lmax_phi)
            call computeAStarT(strcon_matrix,kmax_phi,kmax_phi,scmST)
!           ----------------------------------------------------------
            write(6,'('' time for getStrConstMatrix ='',1f10.5)')     &
                  getTime()-t0
            do kl = 1, kmax_phi
               do klp = 1, kmax_phi
                  if (abs(scmST(kl,klp)-scm(kl,klp)) > tol) then
                     write(6,'(a,4i5)')'i, j, kl, klp = ',i,j,kl,klp
                     call ErrorHandler('testStrConst','scmST <> scm',scmST(kl,klp),scm(kl,klp))
                  endif
               enddo
            enddo
            write(6,'('' i, j, n = '',3i5)')i,j,n
!           write(6,'(2(2i3,2x,2d16.8))')                             &
!              ((i,j,strcon_matrix(klp,kl), klp=1,kmax_kkr),kl=1,kmax_phi)
!           ----------------------------------------------------------
            call writeMatrix('Structure Constants',strcon_matrix,kmax_phi,kmax_phi,tol)
!           ----------------------------------------------------------
         enddo
      enddo
   enddo
!
   deallocate(scm,scmST)
   call endStrConst()
   call endContour()
   call endGauntFactors()
   call endSphericalHarmonics()
   call endSystem()
   call endPotentialType()
   call endScfData()
   call endInput()
   call endDataServiceCenter()
   call endMPP()
!
end program testStrConst
