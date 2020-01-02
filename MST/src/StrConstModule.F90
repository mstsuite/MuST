!  *******************************************************************
!  * Module for calculating KKR structure constant matrix            *
!  * Public functions:                                               *
!  *                                                                 *
!  *    call initStrConst(lmax,Bravais,istop,iprint)                 *
!  *    Purpose: initialize the module for calculating k-space       *
!  *             KKR structure constant matrix for 1 atom per cell   *
!  *    Note:    GauntFactorsModule needs to be initialized first    *
!  *    Input:   lmax = largest possible l-cut off for gij matrix    *
!  *                    (integer)                                    *
!  *             Bravais  = Bravais lattice vector (real array(3,3)) *
!  *             istop = routine name to stop (character string)     *
!  *             iprint= print level (integer)                       *
!  *                                                                 *
!  *    call initStrConst(lmax,num_atoms,atom_posi,Bravais,istop,    *
!  *                      iprint)                                    *
!  *    Purpose: initialize the module for calculating k-space       *
!  *             KKR structure constant matrix                       *
!  *    Note:    GauntFactorsModule needs to be initialized first    *
!  *    Input:   lmax = largest possible l-cut off for gij matrix    *
!  *                    (integer)                                    *
!  *             num_atoms= number of atoms per unit cell (integer)  *
!  *             atom_posi= the position of the atoms in the unit    *
!  *                        cell (real array(1:3,1:num_atoms))       *
!  *             Bravais  = Bravais lattice vector (real array(3,3)) *
!  *             istop = routine name to stop (character string)     *
!  *             iprint= print level (integer)                       *
!  *                                                                 *
!  *    getStrConstMatrix(kvec,kappa,i,j,lmaxi,lmaxj)                *
!  *    Purpose: calculate k-space KKR structure constant matrix     *
!  *    Note:    need to call initStrConst first                     *
!  *    Input:   kvec   = k-vector (real array(3))                   *
!  *             kappa  = sqrt of energy (complex)                   *
!  *             i      = ith atom (integer), 0 < i <= num_atoms     *
!  *             j      = jth atom (integer), 0 < j <= num_atoms     *
!  *             lmaxi  = l-cut off for row index of gij (integer)   *
!  *             lmaxj  = l-cut off for column index of gij (integer)*
!  *             Both lmaxi and lmaxj can not be larger than lmax,   *
!  !             which is previously defined when call initStrConst. *
!  *    Result:  srecon_const = k-space KKR structure constant       *
!  *                    {i, j} matrix (complex array(:,:))  block    *
!  *                    where,  1st dim = kmaxi=(lmaxi+1)**2         *
!  *                            2nd dim = kmaxj=(lmaxj+1)**2         *
!  *                                                                 *
!  *    getStrConstMatrix(kvec,kappa,ni,nj,id,jd,lmaxi,lmaxj)        *
!  *    Purpose: calculate k-space KKR structure constant matrix     *
!  *    Note:    need to call initStrConst first                     *
!  *    Input:   kvec   = k-vector (real array(3))                   *
!  *             kappa  = sqrt of energy (complex)                   *
!  *             ni     = num. of i-array elements, ni < num_atoms   *
!  *             nj     = num. of j-array elements, nj < num_atoms   *
!  *             id     = array of ith atoms, 0 < id <= num_atoms    *
!  *             jd     = array of jth atoms, 0 < jd <= num_atoms    *
!  *             lmaxi  = array of l-cut off for row index of gij    *
!  *             lmaxj  = array of l-cut off for column index of gij *
!  *             Both lmaxi and lmaxj can not be larger than lmax,   *
!  !             which is previously defined when call initStrConst. *
!  *    Result:  srecon_const = k-space KKR structure constant       *
!  *                    {i, j} matrix (complex array(:,:))  block    *
!  *                    where,  1st dim = kmaxi=(lmaxi+1)**2         *
!  *                            2nd dim = kmaxj=(lmaxj+1)**2         *
!  *                                                                 *
!  *    delStrConstMatrix()                                          *
!  *    Purpose: deallocate the structure constanct matrix.          *
!  *                                                                 *
!  *    call endStrConst()                                           *
!  *    Purpose: deallocate the internal allocated arrays and clear  *
!  *             the storage                                         *
!  *    Usage:   should be used whenever one or more atoms change    *
!  *             their position or before the end of the program     *
!  *                                                                 *
!  *                                                                 *
!  *    getFreeElectronPoleSum(kvec,energy)                          *
!  *    Purpose: return -1/pi*sum_G log[(k+G)^2 - E]                 *
!  *    Input:   kvec   = k-vector (real array(3))                   *
!  *             energy = energy in complex type                     *
!  *    Result:  -1/pi*sum_G log[(kvec+G)^2 - energy]                *
!  *                                                                 *
!  *******************************************************************
module StrConstModule
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use KindParamModule, only : CmplxKind
!
   use MathParamModule, only : zero
   use MathParamModule, only : half
   use MathParamModule, only : one
   use MathParamModule, only : two
   use MathParamModule, only : ten2m6
   use MathParamModule, only : ten2m8
   use MathParamModule, only : pi
   use MathParamModule, only : pi2
   use MathParamModule, only : pi4
   use MathParamModule, only : czero
   use MathParamModule, only : cone
   use MathParamModule, only : sqrtm1
!
   use ErrorHandlerModule, only : StopHandler
   use ErrorHandlerModule, only : ErrorHandler
   use ErrorHandlerModule, only : WarningHandler
!
   use SphericalHarmonicsModule, only : calYlmConjg
!
   use PublicTypeDefinitionsModule, only : ScmBlockStruct
!
public :: initStrConst,      &
          endStrConst,       &
          getStrConstMatrix, &
          getFreeElectronPoleSum
!
   interface initStrConst
      module procedure initStrConst_1, initStrConst_m
   end interface
!
   interface getStrConstMatrix
      module procedure getStrConstMatrix0, getStrConstMatrix1
   end interface
!
private
!
   logical :: Initialized = .false.
   logical, parameter :: check_realspace = .false.
!
   character (len=40) :: stop_routine
!
   logical :: real_space_scheme
!
   integer (kind=IntKind) :: iprslat = 0
   integer (kind=IntKind) :: ipknlat = 0
!
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: nrslat
   integer (kind=IntKind) :: nknlat
   integer (kind=IntKind) :: lmax_max
   integer (kind=IntKind) :: lmax_max2, kmax_max2
   integer (kind=IntKind) :: lmax_dlm
   integer (kind=IntKind) :: kmax_dlm
   integer (kind=IntKind), allocatable :: isub(:)
   integer (kind=IntKind), allocatable :: jsub(:)
   integer (kind=IntKind) :: npij
   integer (kind=IntKind) :: numk_inside
   integer (kind=IntKind) :: numk_outside
   integer (kind=IntKind) :: num_emesh
   integer (kind=IntKind) :: print_level
!
   real (kind=RealKind) :: eta
   real (kind=RealKind) :: omegaws
   real (kind=RealKind), allocatable :: rslat_x(:)
   real (kind=RealKind), allocatable :: rslat_y(:)
   real (kind=RealKind), allocatable :: rslat_z(:)
   real (kind=RealKind), allocatable :: rslatsq(:)
   real (kind=RealKind), allocatable :: knlat_x(:)
   real (kind=RealKind), allocatable :: knlat_y(:)
   real (kind=RealKind), allocatable :: knlat_z(:)
   real (kind=RealKind), allocatable :: knlatsq(:)
   real (kind=RealKind), allocatable :: aij(:,:)
!
   real (kind=RealKind) :: Bravais(3,3)
   real (kind=RealKind) :: ScalingFactor
   real (kind=RealKind) :: k_sav(3)
!
   complex (kind=CmplxKind), allocatable :: fac(:)
   complex (kind=CmplxKind), allocatable :: cofk(:)
   complex (kind=CmplxKind), allocatable, target :: cofr(:)
   complex (kind=CmplxKind), allocatable :: kfac_sav(:,:,:)
   complex (kind=CmplxKind), allocatable :: rylm_sav(:,:,:)
   complex (kind=CmplxKind), allocatable :: rfac_sav(:,:,:)
   complex (kind=CmplxKind), allocatable :: dlke(:)
   complex (kind=CmplxKind), allocatable, target :: dterm(:)
   complex (kind=CmplxKind) :: d3term
   complex (kind=CmplxKind) :: kappa_sav
   complex (kind=CmplxKind) :: energy_new
   complex (kind=CmplxKind) :: kappa_new
!
   integer (kind=IntKind), allocatable :: lofk(:)
   integer (kind=IntKind), allocatable :: mofk(:)
!
   complex (kind=CmplxKind), allocatable :: illp(:,:)
   complex (kind=CmplxKind), allocatable :: ilp1(:)
   complex (kind=CmplxKind), allocatable :: facr(:)
   complex (kind=CmplxKind), allocatable, target :: strcon_matrix(:)
!
   type (ScmBlockStruct), allocatable, target :: scm_blocks(:,:)
!
   integer (kind=IntKind) :: ni_sav = 0
   integer (kind=IntKind) :: nj_sav = 0
!
   logical, allocatable :: est_done(:), kst_done(:)
   logical :: kfac_allocated_full = .false.
   logical :: rylm_allocated_full = .false.
   logical :: rfac_allocated_full = .false.
!
   integer (kind=IntKind) :: nknlat_fep = 0
   integer (kind=IntKind) :: ipknlat_fep = 0
   real (kind=RealKind), allocatable :: knlat_x_fep(:)
   real (kind=RealKind), allocatable :: knlat_y_fep(:)
   real (kind=RealKind), allocatable :: knlat_z_fep(:)
   real (kind=RealKind), allocatable :: knlatsq_fep(:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genFactors(lmax)
!  ===================================================================
   implicit   none
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l, m, kl, lp
!
   kl=0
   do l=0,lmax
      do m=-l,l
         kl=kl+1
         lofk(kl)=l
         mofk(kl)=m
      enddo
   enddo
!
   ilp1(0)=sqrtm1
   do l=1,lmax
      ilp1(l)=ilp1(l-1)*sqrtm1
   enddo
   do lp=0,lmax
      do l=0,lmax
         illp(l,lp)=ilp1(l)/ilp1(lp) ! illp = i**(l-lp)
!        illp(l,lp)=-ilp1(l)*ilp1(lp)  ! !!!This needs to be figured out "why?"
      enddo
   enddo
!
   end subroutine genFactors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initStrConst_1(lmax,brav,istop,iprint)
!  ===================================================================
   implicit   none
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(in) :: iprint
!
   integer (kind=IntKind) :: jrs
!
   real (kind=RealKind), intent(in) :: brav(3,3)
!
   real (kind=RealKind) :: bvec(3,3)
   real (kind=RealKind) :: rsij(3)
   real (kind=RealKind) :: sfac, efac
!
   if (Initialized) then
      call WarningHandler('initStrConst','StrConstModule has been initialized')
      return
   else
      Initialized=.true.
   endif
!
   NumAtoms = 1
!
   Bravais(1:3,1:3) = brav(1:3,1:3)
   stop_routine=istop
   print_level=iprint
!
   lmax_max=lmax
   lmax_max2=lmax_max*2
   kmax_max2=(lmax_max2+1)**2
   lmax_dlm=0
!
   allocate( lofk(kmax_max2), mofk(kmax_max2) )
   allocate( illp(0:lmax_max2,0:lmax_max2), ilp1(0:lmax_max2) )
   allocate( facr(0:lmax_max2) )
!  -------------------------------------------------------------------
   call genFactors(lmax_max2)
!  -------------------------------------------------------------------
!
   k_sav(1)=-10000.0d0
   k_sav(2)=-20000.0d0
   k_sav(3)=-30000.0d0
   kappa_sav=(0.0d0,0.0d0)
!
   if (print_level >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,25x,a  )')'********************************'
      write(6,'(  25x,a  )')'*   Output from initStrConst   *'
      write(6,'(  25x,a,/)')'********************************'
   endif
!  ===================================================================
!  calScalingFactor determines the scaling factor, Ewald parameter,
!  iprslat, and ipknlat.
!  -------------------------------------------------------------------
   call calScalingFactor(Bravais,sfac,efac)
!  -------------------------------------------------------------------
   ScalingFactor = sfac
   eta = efac
   if (print_level >= 0) then
      write(6,'('' initStrConst:: Ewald Parameter = '',f10.5)')eta
      write(6,'(''                Scaling Factor  = '',f10.5)')ScalingFactor
   endif
   if (print_level >= 0) then
      write(6,'('' initStrConst:: Ewald Parameter = '',f10.5)')eta
      write(6,'(''                Scaling Factor  = '',f10.5)')ScalingFactor
   endif
!
!  ===================================================================
!  change units so that bravais_lattice and energy are 
!  in the units of ScalingFactor
!  ===================================================================
   bvec(1:3,1:3)=Bravais(1:3,1:3)/ScalingFactor
!  
!  ===================================================================
   allocate(rslat_x(1:iprslat), rslat_y(1:iprslat), rslat_z(1:iprslat), &
            rslatsq(1:iprslat))
   rslat_x(1:iprslat) = ZERO
   rslat_y(1:iprslat) = ZERO
   rslat_z(1:iprslat) = ZERO
   rslatsq(1:iprslat) = ZERO
!  
   allocate(knlat_x(1:ipknlat), knlat_y(1:ipknlat), knlat_z(1:ipknlat),&
            knlatsq(1:ipknlat))
   knlat_x(1:ipknlat) = ZERO
   knlat_y(1:ipknlat) = ZERO
   knlat_z(1:ipknlat) = ZERO
   knlatsq(1:ipknlat) = ZERO
!  ===================================================================
!
!  ===================================================================
!  Allocate memories...............................................
!  -------------------------------------------------------------------
   allocate( isub(1), jsub(1),                                        &
             aij(3,1),fac(-lmax_max2:lmax_max2),                      &
             cofk(0:lmax_max2), cofr(0:lmax_max2),                    &
             kst_done(1), est_done(1),                                &
             dlke(kmax_max2), dterm(kmax_max2) )
!  -------------------------------------------------------------------
   allocate( kfac_sav(1:kmax_max2,ipknlat,1),                         &
             rylm_sav(1:kmax_max2,iprslat,1),                         &
             rfac_sav(0:lmax_max2,iprslat,1) )
!  -------------------------------------------------------------------
   kfac_allocated_full = .false.
   rylm_allocated_full = .false.
   rfac_allocated_full = .false.
!
!  ===================================================================
!  obtain the lattice vectors according to the given Bravais lattice
!  rslat_x, rslat_y, rslat_z, rslatsq, knlat_x, knlat_y, knlat_z, and
!  knlatsq are in the units of ScalingFactor
!  -------------------------------------------------------------------
   call getstruc(bvec,lmax_max2)
!  -------------------------------------------------------------------
!
   npij=1
   aij(1,1)=ZERO
   aij(2,1)=ZERO
   aij(3,1)=ZERO
   isub(1)=1
   jsub(1)=1
   do jrs=1,nrslat
      rsij(1)=rslat_x(jrs)
      rsij(2)=rslat_y(jrs)
      rsij(3)=rslat_z(jrs)
!     ----------------------------------------------------------------
      call calYlmConjg(rsij,lmax_max2,rylm_sav(1:kmax_max2,jrs,1))
!     ----------------------------------------------------------------
   enddo
!
   real_space_scheme = .false.
   kst_done(1)=.false.
   est_done(1)=.false.
!
   ni_sav = 0
   nj_sav = 0
!
   if (print_level >= 0) then
      write(6,'(80(''-''),/)')
   endif
!
   end subroutine initStrConst_1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initStrConst_m(lmax,num_atoms,atom_posi,brav,istop,iprint)
!  ===================================================================
   implicit   none
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: iprint
!
   integer (kind=IntKind) :: i,j,jrs
!
   real (kind=RealKind), intent(in) :: atom_posi(3,num_atoms)
   real (kind=RealKind), intent(in) :: brav(3,3)
!
   real (kind=RealKind), allocatable :: atom_posi_x(:)
   real (kind=RealKind), allocatable :: atom_posi_y(:)
   real (kind=RealKind), allocatable :: atom_posi_z(:)
   real (kind=RealKind) :: bvec(3,3)
   real (kind=RealKind) :: rsij(3)
   real (kind=RealKind) :: sfac, efac, msize
!
   if (Initialized) then
      call WarningHandler('initStrConst','StrConstModule has been initialized')
      return
   else
      Initialized=.true.
   endif
!
   NumAtoms = num_atoms
!
   allocate(atom_posi_x(NumAtoms), atom_posi_y(NumAtoms),atom_posi_z(NumAtoms))
!
   Bravais(1:3,1:3) = brav(1:3,1:3)
   stop_routine=istop
   print_level=iprint
!
   lmax_max=lmax
   lmax_max2=lmax_max*2
   kmax_max2=(lmax_max2+1)**2
   lmax_dlm=0
!
   allocate( strcon_matrix((lmax+1)**2*(lmax+1)**2) )
   allocate( lofk(kmax_max2), mofk(kmax_max2) )
   allocate( illp(0:lmax_max2,0:lmax_max2), ilp1(0:lmax_max2) )
   allocate( facr(0:lmax_max2) )
!  -------------------------------------------------------------------
   call genFactors(lmax_max2)
!  -------------------------------------------------------------------
!
   k_sav(1)=-10000.0d0
   k_sav(2)=-20000.0d0
   k_sav(3)=-30000.0d0
   kappa_sav=(0.0d0,0.0d0)
!
   do i=1,NumAtoms
      atom_posi_x(i)=atom_posi(1,i)
      atom_posi_y(i)=atom_posi(2,i)
      atom_posi_z(i)=atom_posi(3,i)
   enddo
!
   if (print_level >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,25x,a  )')'********************************'
      write(6,'(  25x,a  )')'*   Output from initStrConst   *'
      write(6,'(  25x,a,/)')'********************************'
   endif
!
!  ===================================================================
!  calScalingFactor determines the scaling factor, Ewald parameter,
!  iprslat, and ipknlat.
!  -------------------------------------------------------------------
   call calScalingFactor(Bravais,sfac,efac)
!  -------------------------------------------------------------------
   ScalingFactor = sfac
   eta = efac
   if (print_level >= 0) then
      write(6,'('' initStrConst:: Ewald Parameter = '',f10.5)')eta
      write(6,'(''                Scaling Factor  = '',f10.5)')ScalingFactor
   endif
!
!  ===================================================================
!  change units so that bravais_lattice, atom_posi_* and energy are 
!  in the units of ScalingFactor
!  ===================================================================
   bvec(1:3,1:3)=Bravais(1:3,1:3)/ScalingFactor
   atom_posi_x(1:NumAtoms)=atom_posi_x(1:NumAtoms)/ScalingFactor
   atom_posi_y(1:NumAtoms)=atom_posi_y(1:NumAtoms)/ScalingFactor
   atom_posi_z(1:NumAtoms)=atom_posi_z(1:NumAtoms)/ScalingFactor
!  
!  ===================================================================
   allocate(rslat_x(1:iprslat), rslat_y(1:iprslat), rslat_z(1:iprslat), &
            rslatsq(1:iprslat))
   rslat_x(1:iprslat) = ZERO
   rslat_y(1:iprslat) = ZERO
   rslat_z(1:iprslat) = ZERO
   rslatsq(1:iprslat) = ZERO
!  
   allocate(knlat_x(1:ipknlat), knlat_y(1:ipknlat), knlat_z(1:ipknlat),&
            knlatsq(1:ipknlat))
   knlat_x(1:ipknlat) = ZERO
   knlat_y(1:ipknlat) = ZERO
   knlat_z(1:ipknlat) = ZERO
   knlatsq(1:ipknlat) = ZERO
!  ===================================================================
!
!  ===================================================================
!  Allocate memories...............................................
!  -------------------------------------------------------------------
   allocate( isub(NumAtoms*(NumAtoms-1)+1),                           &
             jsub(NumAtoms*(NumAtoms-1)+1),                           &
             aij(3,NumAtoms*(NumAtoms-1)+1),fac(-lmax_max2:lmax_max2),&
             cofk(0:lmax_max2), cofr(0:lmax_max2),                    &
             kst_done(NumAtoms*(NumAtoms-1)+1),                       &
             est_done(NumAtoms*(NumAtoms-1)+1),                       &
             dlke(kmax_max2), dterm(kmax_max2) )
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Beware!!!
!      The following expression can result negative msize if NumAtoms
!      is large enough, even though msize is declared to be of RealKind.
!              msize = kmax_max2*ipknlat*(NumAtoms*(NumAtoms-1)+1)
!      We need to calculate msize in a different way!!!!!!!!!
!  ===================================================================
   msize=kmax_max2*ipknlat; msize=msize*(NumAtoms*(NumAtoms-1)+1)
   if (msize < 1.0d+06) then
      allocate(kfac_sav(1:kmax_max2,ipknlat,NumAtoms*(NumAtoms-1)+1))
      kfac_allocated_full = .true.
   else
      allocate(kfac_sav(1:kmax_max2,1:ipknlat,2))
      kfac_allocated_full = .false.
   endif
   msize=kmax_max2*iprslat; msize=msize*(NumAtoms*(NumAtoms-1)+1)
   if (msize < 1.0d+06) then
      allocate(rylm_sav(1:kmax_max2,iprslat,NumAtoms*(NumAtoms-1)+1))
      rylm_allocated_full = .true.
   else
      allocate(rylm_sav(1:kmax_max2,1:iprslat,2))
      rylm_allocated_full = .false.
   endif
   msize=(lmax_max2+1)*iprslat; msize=msize*(NumAtoms*(NumAtoms-1)+1)
   if (msize < 1.0d+06) then
      allocate(rfac_sav(0:lmax_max2,iprslat,NumAtoms*(NumAtoms-1)+1))
      rfac_allocated_full = .true.
   else
      allocate(rfac_sav(0:lmax_max2,1:iprslat,2))
      rfac_allocated_full = .false.
   endif
!
!  ===================================================================
!  obtain the lattice vectors according to the given Bravais lattice
!  rslat_x, rslat_y, rslat_z, rslatsq, knlat_x, knlat_y, knlat_z, and
!  knlatsq are in the units of ScalingFactor
!  -------------------------------------------------------------------
   call getstruc(bvec,lmax_max2)
!  -------------------------------------------------------------------
   do jrs=1,nrslat
      rsij(1)=rslat_x(jrs)
      rsij(2)=rslat_y(jrs)
      rsij(3)=rslat_z(jrs)
!     ----------------------------------------------------------------
      call calYlmConjg(rsij,lmax_max2,rylm_sav(1:kmax_max2,jrs,1))
!     ----------------------------------------------------------------
   enddo
!
   npij=1
   aij(1,npij)=zero
   aij(2,npij)=zero
   aij(3,npij)=zero
   isub(npij)=1
   jsub(npij)=1
   do j=1,NumAtoms
      do i=1,NumAtoms
         if(i.ne.j) then
            npij=npij+1
            isub(npij)=i
            jsub(npij)=j
            aij(1,npij)=atom_posi_x(j)-atom_posi_x(i)
            aij(2,npij)=atom_posi_y(j)-atom_posi_y(i)
            aij(3,npij)=atom_posi_z(j)-atom_posi_z(i)
!aij(:,npij) = -aij(:,npij)
!           ==========================================================
!           rylm_sav doesn't depend on energy and k and it can be 
!           calculated once for all
!                           m * -> ->
!               rylm_sav = Y   (Rs+aij)
!                           l
!           ==========================================================
            if (rylm_allocated_full) then
               do jrs=1,nrslat
                  rsij(1)=rslat_x(jrs)+aij(1,npij)
                  rsij(2)=rslat_y(jrs)+aij(2,npij)
                  rsij(3)=rslat_z(jrs)+aij(3,npij)
!                 ----------------------------------------------------
                  call calYlmConjg(rsij,lmax_max2,                    &
                                   rylm_sav(1:kmax_max2,jrs,npij))
!                 ----------------------------------------------------
               enddo
            endif
         endif
      end do
   end do
!
   deallocate(atom_posi_x, atom_posi_y, atom_posi_z)
   real_space_scheme = .false.
!
   do i=1,npij
      est_done(i) = .false.
      kst_done(i) = .false.
   enddo
!
   ni_sav = 0
   nj_sav = 0
!
   if (print_level >= 0) then
      write(6,'(80(''-''),/)')
   endif
!
   end subroutine initStrConst_m
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endStrConst()
!  ===================================================================
   implicit   none
!
   if (Initialized) then
      Initialized=.false.
!     ----------------------------------------------------------------
      deallocate( isub, jsub, aij, fac, cofk, cofr, kfac_sav,         &
                  rylm_sav, rfac_sav, dlke, dterm, kst_done, est_done )
      deallocate( rslat_x, rslat_y, rslat_z, knlat_x, knlat_y, knlat_z,&
                  rslatsq, knlatsq )
!     ----------------------------------------------------------------
   else
      call WarningHandler('endStrConst','StrConstModule is not initialized')
   endif
!
   if (allocated(strcon_matrix)) then 
      deallocate(strcon_matrix)
   endif
!
   if (nknlat_fep > 0) then
      deallocate(knlat_x_fep,knlat_y_fep,knlat_z_fep,knlatsq_fep)
      nknlat_fep = 0
      ipknlat_fep = 0
   endif
!
   if (ni_sav > 0 .or. nj_sav > 0) then
      call deleteScmBlocks()
   endif
!
   deallocate( illp, ilp1, lofk, mofk )
   deallocate( facr )
   iprslat = 0; ipknlat = 0
   ni_sav = 0; nj_sav = 0
!
   kfac_allocated_full = .false.
   rylm_allocated_full = .false.
   rfac_allocated_full = .false.
!
   end subroutine endStrConst
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calScalingFactor(Bravais,sfac,efac)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: Bravais(3,3)
   real (kind=RealKind), intent(out) :: sfac, efac
!
   integer (kind=IntKind) :: nm1, nm2, nm3, iter
   integer (kind=IntKind), parameter :: max_iter = 1000
!
   real (kind=RealKind) :: a0, a1, a2, a3
   real (kind=RealKind) :: vbrar(3,3), vbrak(3,3)
   real (kind=RealKind) :: volr, vfac
   real (kind=RealKind) :: rscut, kncut
!
   real (kind=RealKind), parameter :: tfac = 1.8d0
   real (kind=RealKind), parameter :: fstep = 0.002d0
!
   logical :: done = .false.
!
   a1=sqrt(Bravais(1,1)*Bravais(1,1)+Bravais(2,1)*Bravais(2,1)+ &
           Bravais(3,1)*Bravais(3,1) )
   a2=sqrt(Bravais(1,2)*Bravais(1,2)+Bravais(2,2)*Bravais(2,2)+ &
           Bravais(3,2)*Bravais(3,2) )
   a3=sqrt(Bravais(1,3)*Bravais(1,3)+Bravais(2,3)*Bravais(2,3)+ &
           Bravais(3,3)*Bravais(3,3) )
!
   a0=min(a1,a2,a3)
   efac=half+0.05*max(a1,a2,a3)/a0
!
   sfac=a0/PI2
!
   done = .false.
   iter = 0
   do while (.not. done)
      iter = iter + 1
!     ================================================================
!     scale the Bravais lattice and the reciprical lattice.                         
!     ================================================================
      vbrar(1:3,1:3) = Bravais(1:3,1:3)/sfac
      volr=(vbrar(2,1)*vbrar(3,2)-vbrar(3,1)*vbrar(2,2))*vbrar(1,3)+  &
           (vbrar(3,1)*vbrar(1,2)-vbrar(1,1)*vbrar(3,2))*vbrar(2,3)+  &
           (vbrar(1,1)*vbrar(2,2)-vbrar(2,1)*vbrar(1,2))*vbrar(3,3)
      vfac=PI2/volr
      vbrak(1,1)=vfac*(vbrar(2,2)*vbrar(3,3)-vbrar(3,2)*vbrar(2,3))
      vbrak(2,1)=vfac*(vbrar(3,2)*vbrar(1,3)-vbrar(1,2)*vbrar(3,3))
      vbrak(3,1)=vfac*(vbrar(1,2)*vbrar(2,3)-vbrar(2,2)*vbrar(1,3))
      vbrak(1,2)=vfac*(vbrar(2,3)*vbrar(3,1)-vbrar(3,3)*vbrar(2,1))
      vbrak(2,2)=vfac*(vbrar(3,3)*vbrar(1,1)-vbrar(1,3)*vbrar(3,1))
      vbrak(3,2)=vfac*(vbrar(1,3)*vbrar(2,1)-vbrar(2,3)*vbrar(1,1))
      vbrak(1,3)=vfac*(vbrar(2,1)*vbrar(3,2)-vbrar(3,1)*vbrar(2,2))
      vbrak(2,3)=vfac*(vbrar(3,1)*vbrar(1,2)-vbrar(1,1)*vbrar(3,2))
      vbrak(3,3)=vfac*(vbrar(1,1)*vbrar(2,2)-vbrar(2,1)*vbrar(1,2))
!
!     ================================================================
!     calculate rscut, the radius of real space truncation sphere.....
!     ----------------------------------------------------------------
      call getscut(rsfunc,efac,lmax_max2,                             &
                   vbrar(1:3,1),vbrar(1:3,2),vbrar(1:3,3),rscut,nm1,nm2,nm3)
      call numlat(vbrar,rscut,nm1,nm2,nm3,nrslat)
!     ----------------------------------------------------------------
!
!     ================================================================
!     calculate rscut, the radius of real space truncation sphere.
!     ----------------------------------------------------------------
      call getscut(knfunc,efac,lmax_max2,                             &
                   vbrak(1:3,1),vbrak(1:3,2),vbrak(1:3,3),kncut,nm1,nm2,nm3)
      call numlat(vbrak,kncut,nm1,nm2,nm3,nknlat)
!     ----------------------------------------------------------------
!     write(6,'(a,3i8)')'iter, nrslat, nknlat = ',iter,nrslat,nknlat
      if (iter > max_iter .or. sfac <= 0.1d0) then
!        =============================================================
!        If this message shows up, reduce fstep value. 
!        -------------------------------------------------------------
         call WarningHandler('calScalingFactor',                        &
                             'The scaling factor may not be optimal',sfac)
!        -------------------------------------------------------------
         done = .true.
      else if (nknlat < nrslat/2) then
!        sfac = sfac/tfac
         sfac = sfac-fstep
      else if (nrslat < nknlat/2) then
         sfac = sfac+fstep
      else
         done = .true.
      endif
   enddo
!  write(6,'(a,2i8)')' nrslat, nknlat = ',nrslat,nknlat
!
   iprslat = nrslat
   ipknlat = nknlat
!
   end subroutine calScalingFactor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStrConstMatrix0(k_in,kappa_in,ia,ja,lmaxi,lmaxj,aij_opt) result(scm)
!  ===================================================================
   use MathParamModule, only : ten2m10
!
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
!  *******************************************************************
!  *
!  *  calculate KKR structure constant matrix: strcon_matrix
!  *
!  *  Note:
!  *        * kx, ky, kz, kappa_in, RSpaceBravais, atom_posi must 
!  *          be in the same units (a.u. or d.u.)
!  *
!  *        * the returned value strcon_matrix is in the same units as
!  *          the input quantities.   
!  *        * ia, ja: atom index in the unit cell
!  *        * lmaxi: the lmax for L
!  *        * lmaxj: the lmax for L'
!  *
!  *******************************************************************
!
   implicit   none
!
   character (len=17), parameter :: sname='getStrConstMatrix'
!
   logical :: esame,ksame
!
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: ja
   integer (kind=IntKind), intent(in) :: lmaxi
   integer (kind=IntKind), intent(in) :: lmaxj
!
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: lp
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: kl
   integer (kind=IntKind) :: klp
   integer (kind=IntKind) :: kmaxi, kmaxj
   integer (kind=IntKind) :: nj3
   integer (kind=IntKind), pointer :: kj3(:)
!
   real (kind=RealKind), intent(in) :: k_in(3)
   real (kind=RealKind), intent(out), optional :: aij_opt(3)
!
   real (kind=RealKind) :: kvec(3)
   real (kind=RealKind) :: rs_test
   real (kind=RealKind), pointer :: cgnt(:)
!
   complex (kind=CmplxKind), intent(in) :: kappa_in
!
   complex (kind=CmplxKind) :: pinv
   complex (kind=CmplxKind) :: c
   complex (kind=CmplxKind) :: z
   complex (kind=CmplxKind) :: zinv
!
   complex (kind=CmplxKind), pointer :: hl(:)
!  complex (kind=CmplxKind), pointer :: scm(:,:,:)
   complex (kind=CmplxKind), pointer :: scm(:,:)
!
!
   if(.not.Initialized) then
      call ErrorHandler(sname,'need to call initStrConst first')
   endif
!
   if(abs(kappa_in).lt.ten2m8) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'energy = 0.0')
!     ----------------------------------------------------------------
   else if(ia > NumAtoms .or. ia < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'ia < 1 or ia > num_atoms',ia,NumAtoms)
!     ----------------------------------------------------------------
   else if(ja > NumAtoms .or. ja < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'ja < 1 or ja > num_atoms',ja,NumAtoms)
!     ----------------------------------------------------------------
   else if(lmaxi > lmax_max) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'lmaxi > lmax_max',lmaxi,lmax_max)
!     ----------------------------------------------------------------
   else if(lmaxj > lmax_max) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'lmaxj > lmax_max',lmaxj,lmax_max)
!     ----------------------------------------------------------------
   endif
!
   kmaxi=(lmaxi+1)**2
   kmaxj=(lmaxj+1)**2
!
   scm => aliasArray2_c(strcon_matrix,kmaxi,kmaxj)
!
   if(abs(k_in(1)-k_sav(1))+abs(k_in(2)-k_sav(2))+abs(k_in(3)-k_sav(3)) &
      < ten2m6) then
      ksame=.true.
      if (lmax_dlm < lmaxi+lmaxj) then
         write(6,'(a,3f12.5)')'kvec input: ',k_in(1:3)
         write(6,'(a,3f12.5)')'kvec saved: ',k_sav(1:3)
         call ErrorHandler(sname,'lmaxi or lmaxj has increased',lmaxi,lmaxj)
      endif
   else
      k_sav(1)=k_in(1)
      k_sav(2)=k_in(2)
      k_sav(3)=k_in(3)
      ksame=.false.
      do n=1,npij
         kst_done(n)=.false.
      enddo
   endif
!
   if (abs(kappa_in-kappa_sav) < ten2m6) then
      esame=.true.
      if (lmax_dlm < lmaxi+lmaxj) then
         write(6,'(a,2f12.5)')'Kappa input: ',kappa_in
         write(6,'(a,2f12.5)')'Kappa saved: ',kappa_sav
         call ErrorHandler(sname,'lmaxi or lmaxj has increased',lmaxi,lmaxj)
      endif
   else
      kappa_sav=kappa_in
      esame=.false.
      do n=1,npij
         est_done(n)=.false.
      enddo
   endif
!
!  if (esame .and. ksame) then
!     ----------------------------------------------------------------
!     call WarningHandler(sname,                                      &
!             'matrix has been calculated or invalid kappa or k vector')
!     ----------------------------------------------------------------
!     return
!  endif
!
   lmax_dlm=lmaxi+lmaxj
   kmax_dlm=(lmax_dlm+1)**2
   if (.not.esame) then
      kappa_new=kappa_in*ScalingFactor
      energy_new=kappa_new*kappa_new
!
!     ================================================================
!     determine if real space method can be used:
!               real_space_scheme = .false., use Ewald's method    
!                                   .true., use real space method
!     ================================================================
      real_space_scheme=.false.
      if (check_realspace) then
         hl=>cofr
         z=kappa_new*sqrt(rslatsq(nrslat))
         zinv=cone/z
         hl(0)=-sqrtm1*exp(sqrtm1*z)*zinv
         if(lmax_dlm.gt.0) then
            hl(1)=-sqrtm1*hl(0)*(cone+sqrtm1*zinv)
            rs_test=max(abs(hl(0)),abs(hl(1)))
         else
            rs_test=abs(hl(0))
         endif
         do i=2,lmax_dlm
            hl(i)=(2*i-1)*hl(i-1)*zinv-hl(i-2)
            rs_test=max(abs(hl(i)),rs_test)
         enddo
         nullify( hl )
!        =============================================================
!        Allows for real space approach, if rs_test < 10**(-10)
!        =============================================================
         if(rs_test.lt.ten2m10 .and. check_realspace) then
            real_space_scheme=.true.
            if (print_level >= 0) then
               write(6,'(/,'' Real space method, energy ='',2d15.8,/)')kappa_in*kappa_in
            endif
         endif
      endif
!     real_space_scheme=.true.
!
      if (.not. real_space_scheme) then
!        =============================================================
!        Call cald3 to calculate the third term:
!                                  inf
!        d3term = -sqrt(eta)/pi2 * sum (energy/eta)**n/[(2n-1)*n!]
!                                  n=0
!        -------------------------------------------------------------
         call cald3(energy_new,eta)
!        -------------------------------------------------------------
!
!        =============================================================
!        cofk and cofr are factors of the D1 and D2 terms, repectively.
!
!                                                -l
!        cofk = -4pi/omegaws * exp[energy/eta] * p
!
!                   l+1                  l
!        cofr = (-2)   /sqrt(pi) * (-i/p) ,     where p = sqrt(energy)
!        =============================================================
         pinv=cone/kappa_new
         c=-two*sqrtm1*pinv ! this is how it was programmed originally
!        c=two*sqrtm1*pinv  ! See PRB 55, 12946 (1997)
         cofk(0)=cone
         cofr(0)=cone
         do l=1,lmax_dlm
            cofk(l)=pinv**l
            cofr(l)=c**l
         end do
         c=-pi4/omegaws*exp(energy_new/eta)
         cofk(0:lmax_dlm)=c*cofk(0:lmax_dlm)
         c=-two/sqrt(pi)
         cofr(0:lmax_dlm)=c*cofr(0:lmax_dlm)
      endif
   endif
!
   if (real_space_scheme) then
       write(6,'('' Real space method is enforeced'')')
   endif
!
   if(print_level.ge.2) then
      write(6,'(/,a,'':: energy ='',2d20.13)')sname,kappa_in*kappa_in
      write(6,'('' kx,ky,kz = '',3f10.5)')k_in(1),k_in(2),k_in(3)
   endif
!
   kvec(1)=k_in(1)*ScalingFactor
   kvec(2)=k_in(2)*ScalingFactor
   kvec(3)=k_in(3)*ScalingFactor
!  ===================================================================
!  generate the KKR lattice structure matrix for the energy and
!  (kx,ky,kz), and store result in dlke, in the unit of ScalingFactor
!  ===================================================================
!  na=0
!  do n=1,npij
   if (ia == ja) then
      n = 1
   else
      if (ia < ja) then
         n=(ja-1)*NumAtoms+ia-ja+2
      else
         n=(ja-1)*NumAtoms+ia-ja+1
      endif
      if (ia /= isub(n)) then
         call ErrorHandler(sname,'Inconsistent i sublattice index',ia,isub(n))
      else if (ja /= jsub(n)) then
         call ErrorHandler(sname,'Inconsistent j sublattice index',ja,jsub(n))
      endif
   endif
!  -------------------------------------------------------------------
   call caldlke(esame,ksame,n,kvec,kappa_new,energy_new,isub(n)-jsub(n))
!  -------------------------------------------------------------------
   dlke(1:kmax_dlm)=dlke(1:kmax_dlm)/ScalingFactor
!
   if (present(aij_opt)) then
      aij_opt(1) = aij(1,n)*ScalingFactor
      aij_opt(2) = aij(2,n)*ScalingFactor
      aij_opt(3) = aij(3,n)*ScalingFactor
   endif
!
   if(print_level >= 1) then
      write(6,'(''    l    m                   dlm(k,e)'')')
      do kl=1,kmax_dlm
         if(abs(dlke(kl)).gt.ten2m6) then
            write(6,'(2i5,2x,2d20.13)') lofk(kl),mofk(kl),dlke(kl)
         endif
      end do
   endif
!  na=na+1
!  strcon_matrix(1:kmaxi,1:kmaxj,na)=czero
   scm(1:kmaxi,1:kmaxj)=czero
   do klp=1,kmaxj
      lp=lofk(klp)
      do kl=1,kmaxi
         l=lofk(kl)
!        =============================================================
!        The following three lines are modified on May 16, 2014, and
!        the changes need to be understood algebracally. 
!        The changes are needed in order to make the charge density 
!        symmetry correct when the body-center atom in a BCC structure
!        is displaced along x- and y- directions, respectively.
!        -------------------------------------------------------------
!ywg     nj3 = getNumK3(kl,klp)
!ywg     kj3 => getK3(kl,klp)
!ywg     cgnt => getGauntFactor(kl,klp)
         nj3 = getNumK3(klp,kl)
         kj3 => getK3(klp,kl)
         cgnt => getGauntFactor(klp,kl)
!        =============================================================
         do j=1,nj3
!           strcon_matrix(kl,klp,na) = strcon_matrix(kl,klp,na)      &
!                                    + cgnt(j)*dlke(kj3(j))
            scm(kl,klp)=scm(kl,klp)+cgnt(j)*dlke(kj3(j))
         enddo
!        strcon_matrix(kl,klp,na)=pi4*illp(l+1,lp+1)*strcon_matrix(kl,klp,na)
         scm(kl,klp)=pi4*illp(l,lp)*scm(kl,klp)
      enddo
   enddo
!  if(jsub(n)-isub(n) .eq. 1) then
!     strcon_matrix(1:kmaxi,1:kmaxj,na+1)=strcon_matrix(1:kmaxi,1:kmaxj,1)
!     ----------------------------------------------------------------
!     call zcopy(kmaxi*kmaxj,strcon_matrix(1,1,na),1,                 &
!                strcon_matrix(1,1,na+1),1)
!     ----------------------------------------------------------------
!     na=na+1
!  endif
   if(print_level.ge.2) then
      write(6,'(a)')'    L    L''                  Bllp(k,e)'
      do klp=1,kmaxj
         do kl=1,kmaxi
!           if(abs(strcon_matrix(kl,klp,n)).gt.ten2m6) then
!              write(6,'(2i5,2x,2d20.13)')kl,klp,strcon_matrix(kl,klp,na)
            if(abs(scm(kl,klp)).gt.ten2m6) then
               write(6,'(2i5,2x,2d20.13)')kl,klp,scm(kl,klp)
            endif        
         enddo
      enddo
   endif        
!  enddo
!
   nullify( kj3, cgnt )
!
   end function getStrConstMatrix0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStrConstMatrix1(k_in,kappa_in,ni,nj,id,jd,lmaxi,lmaxj) result(scm)
!  ===================================================================
   use MathParamModule, only : ten2m10
!
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
!  *******************************************************************
!  *
!  *  calculate KKR structure constant matrix: strcon_matrix
!  *
!  *  Note:
!  *        * kx, ky, kz, kappa_in, RSpaceBravais, atom_posi must 
!  *          be in the same units (a.u. or d.u.)
!  *
!  *        * the returned value strcon_matrix is in the same units as
!  *          the input quantities.   
!  *        * ni     = num. of i-array elements, ni < num_atoms   *
!  *        * nj     = num. of j-array elements, nj < num_atoms   *
!  *        * id     = array of ith atoms, 0 < id <= num_atoms    *
!  *        * jd     = array of jth atoms, 0 < jd <= num_atoms    *
!  *        * lmaxi: the array of lmax for L
!  *        * lmaxj: the array of lmax for L'
!  *
!  *******************************************************************
!
   implicit   none
!
   character (len=17), parameter :: sname='getStrConstMatrix'
!
   logical :: esame,ksame
!
   integer (kind=IntKind), intent(in) :: ni
   integer (kind=IntKind), intent(in) :: nj
   integer (kind=IntKind), intent(in) :: id(1:ni)
   integer (kind=IntKind), intent(in) :: jd(1:nj)
   integer (kind=IntKind), intent(in) :: lmaxi(1:ni)
   integer (kind=IntKind), intent(in) :: lmaxj(1:nj)
!
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: lp
   integer (kind=IntKind) :: i, j, ip, jp, ia, ja
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: kl
   integer (kind=IntKind) :: klp
   integer (kind=IntKind) :: kmaxi, kmaxj
   integer (kind=IntKind) :: nj3
   integer (kind=IntKind), pointer :: kj3(:)
!
   real (kind=RealKind), intent(in) :: k_in(3)
!
   real (kind=RealKind) :: kvec(3)
   real (kind=RealKind) :: rs_test
   real (kind=RealKind), pointer :: cgnt(:)
!
   complex (kind=CmplxKind), intent(in) :: kappa_in
!
   complex (kind=CmplxKind) :: pinv
   complex (kind=CmplxKind) :: c
   complex (kind=CmplxKind) :: z
   complex (kind=CmplxKind) :: zinv
!
   complex (kind=CmplxKind), pointer :: hl(:)
   complex (kind=CmplxKind), pointer :: scmm(:,:)
   type (ScmBlockStruct), pointer :: scm(:,:)
!
   if(.not.Initialized) then
      call ErrorHandler(sname,'need to call initStrConst first')
   endif
!
   if(abs(kappa_in).lt.ten2m8) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'energy = 0.0')
!     ----------------------------------------------------------------
   else if(ni > NumAtoms .or. ni < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'ni < 1 or ni > num_atoms',ni,NumAtoms)
!     ----------------------------------------------------------------
   else if(nj > NumAtoms .or. nj < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'nj < 1 or nj > num_atoms',nj,NumAtoms)
!     ----------------------------------------------------------------
   endif
!
   do jp = 1, nj
      if(lmaxj(jp) > lmax_max) then
!        -------------------------------------------------------------
         call ErrorHandler(sname,'lmaxj > lmax_max',lmaxj(jp),lmax_max)
!        -------------------------------------------------------------
      endif
   enddo
   do ip = 1, ni
      if(lmaxi(ip) > lmax_max) then
!        -------------------------------------------------------------
         call ErrorHandler(sname,'lmaxi > lmax_max',lmaxi(ip),lmax_max)
!        -------------------------------------------------------------
      endif
   enddo
!
   if (ni_sav < 1 .or. nj_sav < 1) then
!     ----------------------------------------------------------------
      call allocateScmBlocks(ni,nj,lmaxi,lmaxj)
!     ----------------------------------------------------------------
   else if (ni_sav < ni .or. nj_sav < nj) then
!     ----------------------------------------------------------------
      call deleteScmBlocks()
      call allocateScmBlocks(ni,nj,lmaxi,lmaxj)
!     ----------------------------------------------------------------
   else
      LOOP_nj: do jp = 1, nj
         do ip = 1, ni
            if ( scm_blocks(ip,jp)%lmaxi < lmaxi(ip) .or.                &
                 scm_blocks(ip,jp)%lmaxj < lmaxj(jp) ) then
!              -------------------------------------------------------
               call deleteScmBlocks()
               call allocateScmBlocks(ni,nj,lmaxi,lmaxj)
!              -------------------------------------------------------
               exit LOOP_nj
            endif
         enddo
      enddo LOOP_nj
   endif
   scm => scm_blocks(1:ni,1:nj)
!
!do jp = 1, nj
!   ja = jd(jp)
!   do ip = 1, ni
!      ia = id(ip)
!   scmm => getStrConstMatrix0(k_in,kappa_in,ia,ja,lmaxi(ip),lmaxj(jp))
!   scm_blocks(ip,jp)%strcon_matrix = scmm
!   enddo
!enddo
!return
!
   if(abs(k_in(1)-k_sav(1))+abs(k_in(2)-k_sav(2))+abs(k_in(3)-k_sav(3)) &
      < ten2m6) then
      ksame=.true.
      do jp = 1, nj
         do ip = 1, ni
            if (lmax_dlm < lmaxi(ip)+lmaxj(jp)) then
               write(6,'(a,3f12.5)')'kvec input: ',k_in(1:3)
               write(6,'(a,3f12.5)')'kvec saved: ',k_sav(1:3)
               call ErrorHandler(sname,'lmaxi or lmaxj has increased',  &
                                 lmaxi(ip),lmaxj(jp))
            endif
         enddo
      enddo
   else
      k_sav(1)=k_in(1)
      k_sav(2)=k_in(2)
      k_sav(3)=k_in(3)
      ksame=.false.
      do n=1,npij
         kst_done(n)=.false.
      enddo
   endif
!
   if (abs(kappa_in-kappa_sav) < ten2m6) then
      esame=.true.
      do jp = 1, nj
         do ip = 1, ni
            if (lmax_dlm < lmaxi(ip)+lmaxj(jp)) then
               write(6,'(a,2f12.5)')'Kappa input: ',kappa_in
               write(6,'(a,2f12.5)')'Kappa saved: ',kappa_sav
               call ErrorHandler(sname,'lmaxi or lmaxj has increased',&
                                 lmaxi(ip),lmaxj(jp))
            endif
         enddo
      enddo
   else
      kappa_sav=kappa_in
      esame=.false.
      do n=1,npij
         est_done(n)=.false.
      enddo
   endif
!
!  if (esame .and. ksame) then
!     ----------------------------------------------------------------
!     call WarningHandler(sname,                                      &
!             'matrix has been calculated or invalid kappa or k vector')
!     ----------------------------------------------------------------
!     return
!  endif
!
   lmax_dlm = 0
   do jp = 1, nj
      do ip = 1, ni
         lmax_dlm = max(lmaxi(ip)+lmaxj(jp), lmax_dlm)
      enddo
   enddo
   kmax_dlm=(lmax_dlm+1)**2
   if (.not.esame) then
      kappa_new=kappa_in*ScalingFactor
      energy_new=kappa_new*kappa_new
!
!     ================================================================
!     determine if real space method can be used:
!               real_space_scheme = .false., use Ewald's method    
!                                   .true., use real space method
!     ================================================================
      real_space_scheme=.false.
      if (check_realspace) then
         hl=>cofr
         z=kappa_new*sqrt(rslatsq(nrslat))
         zinv=cone/z
         hl(0)=-sqrtm1*exp(sqrtm1*z)*zinv
         if(lmax_dlm.gt.0) then
            hl(1)=-sqrtm1*hl(0)*(cone+sqrtm1*zinv)
            rs_test=max(abs(hl(0)),abs(hl(1)))
         else
            rs_test=abs(hl(0))
         endif
         do i=2,lmax_dlm
            hl(i)=(2*i-1)*hl(i-1)*zinv-hl(i-2)
            rs_test=max(abs(hl(i)),rs_test)
         enddo
         nullify( hl )
!        =============================================================
!        Allows for real space approach, if rs_test < 10**(-10)
!        =============================================================
         if (rs_test.lt.ten2m10) then
            real_space_scheme=.true.
            if (print_level >= 0) then
               write(6,'(/,'' Real space method, energy ='',2d15.8,/)')kappa_in*kappa_in
            endif
         endif
      endif
!     real_space_scheme=.true.
!
      if (.not. real_space_scheme) then
!        =============================================================
!        Call cald3 to calculate the third term:
!                                  inf
!        d3term = -sqrt(eta)/pi2 * sum (energy/eta)**n/[(2n-1)*n!]
!                                  n=0
!        -------------------------------------------------------------
         call cald3(energy_new,eta)
!        -------------------------------------------------------------
!
!        =============================================================
!        cofk and cofr are factors of the D1 and D2 terms, repectively.
!
!                                                -l
!        cofk = -4pi/omegaws * exp[energy/eta] * p
!
!                   l+1                 l
!        cofr = (-2)   /sqrt(pi) * (i/p) ,     where p = sqrt(energy)
!        =============================================================
         pinv=cone/kappa_new
         c=-two*sqrtm1*pinv ! this is how it was programmed originally
!        c=two*sqrtm1*pinv  ! See PRB 55, 12946 (1997)
         cofk(0)=cone
         cofr(0)=cone
         do l=1,lmax_dlm
            cofk(l)=pinv**l
            cofr(l)=c**l
         end do
         c=-pi4/omegaws*exp(energy_new/eta)
         cofk(0:lmax_dlm)=c*cofk(0:lmax_dlm)
         c=-two/sqrt(pi)
         cofr(0:lmax_dlm)=c*cofr(0:lmax_dlm)
      endif
   endif
!
!  ===================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This is a temporary fix of the problem: we force to run real space method
!  We need to look into this problem  3-6-2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ===================================================================
   if (real_space_scheme) then
       write(6,'('' Real space method is enforeced'')')
   endif
!  ===================================================================
!
   if(print_level.ge.2) then
      write(6,'(/,a,'':: energy ='',2d20.13)')sname,kappa_in*kappa_in
      write(6,'('' kx,ky,kz = '',3f10.5)')k_in(1),k_in(2),k_in(3)
   endif
!
   kvec(1)=k_in(1)*ScalingFactor
   kvec(2)=k_in(2)*ScalingFactor
   kvec(3)=k_in(3)*ScalingFactor
!
!  ===================================================================
!  generate the KKR lattice structure matrix for the energy and
!  (kx,ky,kz), and store result in dlke, in the unit of ScalingFactor
!  ===================================================================
   do jp = 1, nj
      ja = jd(jp)
      do ip = 1, ni
         ia = id(ip)
         kmaxi=(lmaxi(ip)+1)**2
         kmaxj=(lmaxj(jp)+1)**2
         if (ia == ja) then
            n = 1
         else
            if (ia < ja) then
               n=(ja-1)*NumAtoms+ia-ja+2
            else
               n=(ja-1)*NumAtoms+ia-ja+1
            endif
            if (ia /= isub(n)) then
               call ErrorHandler(sname,'Inconsistent i sublattice index',ia,isub(n))
            else if (ja /= jsub(n)) then
               call ErrorHandler(sname,'Inconsistent j sublattice index',ja,jsub(n))
            endif
         endif
!        -------------------------------------------------------------
         call caldlke(esame,ksame,n,kvec,kappa_new,energy_new,isub(n)-jsub(n))
!        -------------------------------------------------------------
         dlke(1:kmax_dlm)=dlke(1:kmax_dlm)/ScalingFactor
!
         if(print_level >= 1) then
            write(6,'(''    l    m                   dlm(k,e)'')')
            do kl=1,kmax_dlm
               if(abs(dlke(kl)).gt.ten2m6) then
                  write(6,'(2i5,2x,2d20.13)') lofk(kl),mofk(kl),dlke(kl)
               endif
            end do
         endif
         scm_blocks(ip,jp)%strcon_matrix(1:kmaxi,1:kmaxj)=czero
         do klp=1,kmaxj
            lp=lofk(klp)
            do kl=1,kmaxi
               l=lofk(kl)
!              =======================================================
!              The following three lines are modified on May 16, 2014, and
!              the changes need to be understood algebracally. 
!              The changes are needed in order to make the charge density 
!              symmetry correct when the body-center atom in a BCC structure
!              is displaced along x- and y- directions, respectively.
!              -------------------------------------------------------
!              nj3 = getNumK3(kl,klp)
!              kj3 => getK3(kl,klp)
!              cgnt => getGauntFactor(kl,klp)
               nj3 = getNumK3(klp,kl)
               kj3 => getK3(klp,kl)
               cgnt => getGauntFactor(klp,kl)
!              =======================================================
               do j=1,nj3
                  scm_blocks(ip,jp)%strcon_matrix(kl,klp) =            &
                      scm_blocks(ip,jp)%strcon_matrix(kl,klp)+cgnt(j)*dlke(kj3(j))
               enddo
               scm_blocks(ip,jp)%strcon_matrix(kl,klp) =               &
                      pi4*illp(l,lp)*scm_blocks(ip,jp)%strcon_matrix(kl,klp)
            enddo
         enddo
         if(print_level.ge.1) then
            write(6,'(a)')'    L    L''                  Bllp(k,e)'
            do klp=1,kmaxj
               do kl=1,kmaxi
                  if(abs(scm_blocks(ip,jp)%strcon_matrix(kl,klp)).gt.ten2m6) then
                     write(6,'(2i5,2x,2d20.13)')kl,klp,                &
                           scm_blocks(ip,jp)%strcon_matrix(kl,klp)
                  endif        
               enddo
            enddo
         endif        
      enddo
   enddo
!
   nullify( kj3, cgnt )
!
   end function getStrConstMatrix1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocateScmBlocks(ni,nj,lmaxi,lmaxj)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ni, nj
   integer (kind=IntKind), intent(in) :: lmaxi(1:ni), lmaxj(1:nj)
   integer (kind=IntKind) :: i, j, kmaxi, kmaxj
!
   allocate( scm_blocks(ni,nj) )
   do j = 1, nj
      kmaxj = (lmaxj(j)+1)**2
      do i = 1, ni
         kmaxi = (lmaxi(i)+1)**2
         scm_blocks(i,j)%lmaxj = lmaxj(j)
         scm_blocks(i,j)%lmaxi = lmaxi(i)
         allocate( scm_blocks(i,j)%strcon_matrix(kmaxi,kmaxj) )
      enddo
   enddo
   ni_sav = ni
   nj_sav = nj
!
   end subroutine allocateScmBlocks
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteScmBlocks()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, j
!
   do j = 1, nj_sav
      do i = 1, ni_sav
         deallocate( scm_blocks(i,j)%strcon_matrix )
      enddo
   enddo
   deallocate( scm_blocks )
   ni_sav = 0
   nj_sav = 0
!
   end subroutine deleteScmBlocks
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine caldlke(esame,ksame,n,kvec,kappa,energy,imj)
!  ===================================================================
!
!  *******************************************************************
!  *                                                                 *
!  * The subroutine to calculate the lattice structure constant      *
!  *                                                                 *
!  *                      ij  _                                      *
!  *                     D   (k,e), l=0,1,...,lmax_dlm               *
!  *                      l,m                                        *
!  *                                                                 *
!  * for a general (simple or complex) lattice at energy e and       *
!  * the reciprocal space vector kvec, and store them in dlke.       *
!  *                                                                 *
!  * Notice that Dlm are calculated based on the spherical harmonics *
!  * In order to get the Dlm based on the real spherical harmonics   *
!  * one may use the following formular provide e being real:        *
!  *                                                                 *
!  *    ij  _                         ij  _          m               *
!  *   d   (k,e) = sqrt(2.0d0) * Re[ D   (k,e) * (-1)  ],  m > 0     *
!  *    l,m                           l,m                            *
!  *                                                                 *
!  *                                                                 *
!  *                                                                 *
!  *    ij  _                         ij   _          m              *
!  *   d   (k,e) = sqrt(2.0d0) * Im[ D    (k,e) * (-1)  ], m < 0     *
!  *    l,m                           l,|m|                          *
!  *                                                                 *
!  *                                                                 *
!  *                                                                 *
!  *    ij  _       ij  _                                            *
!  *   d   (k,e) = D   (k,e)                            m = 0        *
!  *    l,0         l,0                                              *
!  *                                                                 *
!  *                                                                 *
!  * where,                                                          *
!  *                       ij  _                                     *
!  *                      d   (k,e), l=0,1,...,lmax_dlm              *
!  *                       l,m                                       *
!  *                                                                 *
!  * stands for the Dlm's based on the real spherical harmonics.     *
!  *                                                                 *
!  *                                                                 *
!  * April 20th, 1992.   Yang Wang                                   *
!  *                                                                 *
!  * updated on Aug. 12th, 1993  by Yang Wang                        *
!  * updated on Nov, 3rd,  1995  by Yang Wang                        *
!  * updated on Apr, 10th, 1996  by Yang Wang                        *
!  * updated on Jun, 22nd, 1996  by Yang Wang                        *
!  *                                                                 *
!  *******************************************************************
!
   implicit   none
!
   character (len= 7), parameter :: sname='caldlke'
!
   logical, intent(in) :: ksame,esame
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: imj
!
   integer (kind=IntKind) :: j0
   integer (kind=IntKind) :: kl
!
   real (kind=RealKind), intent(in) :: kvec(3)
!
   complex (kind=CmplxKind), intent(in) :: kappa
   complex (kind=CmplxKind), intent(in) :: energy
!
   if(imj.eq.0) then
      j0=2
   else
      j0=1
   endif
!
   if(.not.real_space_scheme) then
!     ================================================================
!     calculate the structure constants using Ewald's method.
!
!     Call knterm to calculate the term related to summation
!     over knlat. The result is stored in d2lke.
!
!     Call rsterm to calculate the term related to summation
!     over rslat. The result is stored in d1lke.
!     ----------------------------------------------------------------
      call knterm(ksame,n,kvec,energy)
!     ----------------------------------------------------------------
      do kl=1,kmax_dlm
         dlke(kl) = cofk(lofk(kl))*dterm(kl)
      enddo
!     ----------------------------------------------------------------
      call rsterm(esame,n,j0,kvec,energy)
!     ----------------------------------------------------------------
      do kl=1,kmax_dlm
         dlke(kl) = dlke(kl) + cofr(lofk(kl))*dterm(kl)
      enddo
!     ================================================================
!     add the third term:
!                               inf             n
!     d3term = -sqrt(eta)/pi2 * sum (energy/eta) /[(2n-1)n!] 
!                               n=0
!     ================================================================
      if(imj.eq.0) then
         dlke(1)=dlke(1)+d3term
      endif
   else
!     ================================================================
!     calculate the structure constants using real space summation.
!     ----------------------------------------------------------------
      call rdlke(esame,n,j0,kvec,kappa)
!     ----------------------------------------------------------------
   endif
!
   if (stop_routine .eq. sname) then
      call StopHandler(sname)
   end if
!
   end subroutine caldlke
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cald3(energy,eta)
!  ===================================================================
   use MathParamModule, only : ten2m14
!
   implicit   none
!
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: eta
   real (kind=RealKind) :: rn
!
   complex (kind=CmplxKind), intent(in) :: energy
!
   complex (kind=CmplxKind) :: sum
   complex (kind=CmplxKind) :: sum0
   complex (kind=CmplxKind) :: term
!   integer (kind=IntKind) ::   n
!   complex (kind=CmplxKind) :: eoe
!
!  *******************************************************************
!  *                           inf                                   *
!  * d3term = -sqrt(eta)/pi2 * sum (energy/eta)**n/[(2n-1)*n!]       *
!  *                           n=0                                   *
!  *******************************************************************
!
!!!eoe = energy/eta
!!!term = CONE
!!!rn = ZERO
!!!do while (abs(term) > ten2m14)
!!!   rn = rn + one
!!!   term = term*eoe/rn 
!!!enddo
!!!sum = term/(two*rn-one)
!!!n = rn
!!!do i = n-1,1,-1
!!!   rn = i
!!!   term = term*rn/eoe
!!!   sum=sum+term/(two*rn-one)
!!!enddo
!!!sum = sum - CONE ! add n=0 term
   rn=zero
   sum0=czero
   sum=-cone
   term=cone
   do while (abs(sum-sum0).ge.ten2m8)
      sum0=sum
      do i=1,50
         rn=rn+one
         term=term*energy/(eta*rn)
         sum=sum+term/(two*rn-one)
      end do
   end do
   d3term=-sqrt(eta)/pi2*sum
!
   end subroutine cald3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getscut(func,eta,lmax,a1,a2,a3,scut,nm1,nm2,nm3)
!  ===================================================================
   use MathParamModule, only : ten2m14
   implicit   none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(out) :: nm1
   integer (kind=IntKind), intent(out) :: nm2
   integer (kind=IntKind), intent(out) :: nm3
!
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: k
!
   real (kind=RealKind), intent(in) ::  eta
   real (kind=RealKind), intent(in) ::  a1(3)
   real (kind=RealKind), intent(in) ::  a2(3)
   real (kind=RealKind), intent(in) ::  a3(3)
   real (kind=RealKind), intent(out) ::  scut
!
   real (kind=RealKind) ::  r(3)
   real (kind=RealKind) ::  rm
   real (kind=RealKind) ::  term
!
   interface
      function func(x,eta,l)
         use KindParamModule, only : IntKind
         use KindParamModule, only : RealKind
         real (kind=RealKind) :: func
         real (kind=RealKind) :: x
         real (kind=RealKind) :: eta
         integer (kind=IntKind)  :: l
      end function func
   end interface
!
!  ===================================================================
!  calculate nm1,nm2,nm3...........................................
!  ===================================================================
   r(1)=sqrt(a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3))
   term=one
   nm1=0
   do while(term.gt.ten2m14)
      nm1=nm1+1
      rm=nm1*r(1)
      term=func(rm,eta,lmax)
   enddo
   r(2)=sqrt(a2(1)*a2(1)+a2(2)*a2(2)+a2(3)*a2(3))
   term=one
   nm2=0
   do while(term.gt.ten2m14)
      nm2=nm2+1
      rm=nm2*r(2)
      term=func(rm,eta,lmax)
   enddo
   r(3)=sqrt(a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3))
   term=one
   nm3=0
   do while(term.gt.ten2m14)
      nm3=nm3+1
      rm=nm3*r(3)
      term=func(rm,eta,lmax)
   enddo
!
!  ===================================================================
!  calculate scut..................................................
!  ===================================================================
   scut=r(1)*nm1
   do i=-1,1
      r(1)=i*a1(1)*nm1
      r(2)=i*a1(2)*nm1
      r(3)=i*a1(3)*nm1
      do j=-1,1
         r(1)=r(1)+j*a2(1)*nm2
         r(2)=r(2)+j*a2(2)*nm2
         r(3)=r(3)+j*a2(3)*nm2
         do k=-1,1
            r(1)=r(1)+k*a3(1)*nm3
            r(2)=r(2)+k*a3(2)*nm3
            r(3)=r(3)+k*a3(3)*nm3
            rm=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
            scut=max(scut,rm)
         enddo
      enddo
   enddo
!
   end subroutine getscut
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getstruc(vbrar,lmax_dlm)
!  ===================================================================
!
   implicit   none
!
   character (len= 8), parameter ::  sname='getstruc'
!
   integer (kind=IntKind), intent(in) :: lmax_dlm
   integer (kind=IntKind) :: n1
   integer (kind=IntKind) :: nm1
   integer (kind=IntKind) :: nm2
   integer (kind=IntKind) :: nm3
!
   real (kind=RealKind), intent(in) :: vbrar(3,3)
!
   real (kind=RealKind) :: vbrak(3,3)
   real (kind=RealKind) :: rscut
   real (kind=RealKind) :: kncut
   real (kind=RealKind) :: factor
!
!  *******************************************************************
!  *  Sets up real space and receprocal space Bravais lattice vectors
!  *
!  *  input:  vbrar   = basic real space Bravais lattice vectors
!  *          print_level  = print instruction parameter
!  *          stop_routine   = routine name for stop
!  *******************************************************************
!
!  ===================================================================
!  calculate rscut, the radius of real space truncation sphere.....
!  -------------------------------------------------------------------
   call getscut(rsfunc,eta,lmax_dlm,                                  &
                vbrar(1:3,1),vbrar(1:3,2),vbrar(1:3,3),rscut,nm1,nm2,nm3)
   call numlat(vbrar,rscut,nm1,nm2,nm3,nrslat)
!  -------------------------------------------------------------------
   if (nrslat > iprslat) then
      call ErrorHandler(sname,'nrslat > iprslat',nrslat,iprslat)
   endif
!
!  ===================================================================
!  generate the real space lattice vectors.
!  -------------------------------------------------------------------
   call lattice(vbrar,rscut,nm1,nm2,nm3,                              &
                rslat_x,rslat_y,rslat_z,rslatsq,nrslat,iprslat)
!  -------------------------------------------------------------------
   if(print_level.ge.0) then
      write(6,'(/,'' GETSTRUC:: nm1, nm2, nm3 = '',3i5)')nm1,nm2,nm3
      write(6,'(  ''            Rs cut radius = '',1f10.5)') rscut
      write(6,'(  ''            Number of Rs  = '',i5)') nrslat
   endif
!  ===================================================================
!  calculate the bravais lattice cell volume.
!  ===================================================================
   omegaws=(vbrar(2,1)*vbrar(3,2)-vbrar(3,1)*vbrar(2,2))*vbrar(1,3)+   &
           (vbrar(3,1)*vbrar(1,2)-vbrar(1,1)*vbrar(3,2))*vbrar(2,3)+   &
           (vbrar(1,1)*vbrar(2,2)-vbrar(2,1)*vbrar(1,2))*vbrar(3,3)
!
!  ===================================================================
!  generate basis vectors for reciprocal space.
!  ===================================================================
   factor=PI2/omegaws
   vbrak(1,1)=factor*(vbrar(2,2)*vbrar(3,3)-vbrar(3,2)*vbrar(2,3))
   vbrak(2,1)=factor*(vbrar(3,2)*vbrar(1,3)-vbrar(1,2)*vbrar(3,3))
   vbrak(3,1)=factor*(vbrar(1,2)*vbrar(2,3)-vbrar(2,2)*vbrar(1,3))
   vbrak(1,2)=factor*(vbrar(2,3)*vbrar(3,1)-vbrar(3,3)*vbrar(2,1))
   vbrak(2,2)=factor*(vbrar(3,3)*vbrar(1,1)-vbrar(1,3)*vbrar(3,1))
   vbrak(3,2)=factor*(vbrar(1,3)*vbrar(2,1)-vbrar(2,3)*vbrar(1,1))
   vbrak(1,3)=factor*(vbrar(2,1)*vbrar(3,2)-vbrar(3,1)*vbrar(2,2))
   vbrak(2,3)=factor*(vbrar(3,1)*vbrar(1,2)-vbrar(1,1)*vbrar(3,2))
   vbrak(3,3)=factor*(vbrar(1,1)*vbrar(2,2)-vbrar(2,1)*vbrar(1,2))
!
   omegaws=abs(omegaws)
!
!  ===================================================================
!  calculate rscut, the radius of real space truncation sphere.
!  -------------------------------------------------------------------
   call getscut(knfunc,eta,lmax_dlm,                                  &
                vbrak(1:3,1),vbrak(1:3,2),vbrak(1:3,3),kncut,nm1,nm2,nm3)
   call numlat(vbrak,kncut,nm1,nm2,nm3,nknlat)
!  -------------------------------------------------------------------
   if (nknlat > ipknlat) then
      call ErrorHandler(sname,'nrslat > iprslat',nrslat,iprslat)
   endif
!
!  ===================================================================
!  generate the reciprocal space lattice vectors.
!  -------------------------------------------------------------------
   call lattice(vbrak,kncut,nm1,nm2,nm3,                               &
                knlat_x,knlat_y,knlat_z,knlatsq,nknlat,ipknlat)
!  -------------------------------------------------------------------
   if(print_level.ge.0) then
      write(6,'(/,''            nm1, nm2, nm3 = '',3i5)') nm1,nm2,nm3
      write(6,'(  ''            Kn cut radius = '',1f10.5)') kncut
      write(6,'(  ''            Number of Kn  = '',i5)') nknlat
      write(6,'(/)')
      write(6,'(  '' Real Space Unit Cell Vol = '',1f20.13)')omegaws
      write(6,'(  ''    k-Space Unit Cell Vol = '',1f20.13)')PI2**3/omegaws
      write(6,'(/)')
   endif
   if(print_level.ge.2) then
      write(6,'(/)')
      write(6,'(16x,''n'',20x,''rslat'',18x,''rslatsq'')')
      write(6,'(12x,56(''=''))')
      write(6,'(12x,1i5,2x,4f12.5)')                                   &
           (n1,rslat_x(n1),rslat_y(n1),rslat_z(n1),rslatsq(n1),n1=1,nrslat)
      write(6,'(/)')
      write(6,'(14x,''n'',20x,''knlat'',18x,''knlatsq'')')
      write(6,'(12x,56(''=''))')
      write(6,'(12x,1i5,2x,4f12.5)')                                   &
           (n1,knlat_x(n1),knlat_y(n1),knlat_z(n1),knlatsq(n1),n1=1,nknlat)
   endif
!
!  ===================================================================
   if(sname.eq.stop_routine) then
      call StopHandler(sname)
   endif
!
   end subroutine getstruc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine intfac(rij,e,eta,lmax,dqint)
!  ===================================================================
   use MathParamModule, only : fourth
   use MathParamModule, only : four
   implicit   none
!
   character (len=6), parameter :: sname='intfac'
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: n2
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: rij
   real (kind=RealKind), intent(in) :: eta
!
   real (kind=RealKind) :: a
   real (kind=RealKind) :: r
   real (kind=RealKind) :: rfac
   real (kind=RealKind) :: x1
   real (kind=RealKind) :: y1
   real (kind=RealKind) :: x2
   real (kind=RealKind) :: y2
   real (kind=RealKind) :: xy12
   real (kind=RealKind) :: xy22
   real (kind=RealKind) :: erfc
   real (kind=RealKind) :: fn1
   real (kind=RealKind) :: gn1
   real (kind=RealKind) :: fn2
   real (kind=RealKind) :: gn2
!
   complex (kind=CmplxKind), intent(out) :: dqint(0:lmax)
   complex (kind=CmplxKind), intent(in) :: e
!
   complex (kind=CmplxKind) :: efac, bfac
   complex (kind=CmplxKind) :: b
   complex (kind=CmplxKind) :: t1
   complex (kind=CmplxKind) :: t2
   complex (kind=CmplxKind) :: remain1
   complex (kind=CmplxKind) :: remain2
   complex (kind=CmplxKind) :: dqm1
!
!  *******************************************************************
!  *  call intfac to calculate the integral:
!  *
!  *          -> ->  -1/2  inf       2*l
!  *  I(l) = |Rs+aij|    * int dx * x   * exp(-x**2*r**2-b**2/x**2)
!  *                        a      
!  *
!  *  and store it in dqint(l)
!  *
!  *  l = 0, 1, ..., lmax.
!  *
!  *                  -> ->              -> ->                 -> ->
!  *  a=0.5*sqrt(eta*|Rs+aij|), r=sqrt(|(Rs+aij|), b=i*sqrt(e*|Rs+aij|)/2.
!  *
!  *  using: (2*l+1)*I(l) = 2*r*r*I(l+1) - 2*b*b*I(l-1) - c(l) ,
!  *
!  *                                                        
!  *         c(l) = a**(2*l+1) * exp(-a**2*r**2-b**2/a**2) / r
!  *
!  *         I(0) = sqrt(pi)/(4*r*r) * [+exp(+2*r*b)*erfc(r*a+b/a))
!  *                                    +exp(-2*r*b)*erfc(r*a-b/a)) ]
!  *
!  *         I(-1)= sqrt(pi)/(4*b*r) * [-exp(+2*r*b)*erfc(r*a+b/a)
!  *                                    +exp(-2*r*b)*erfc(r*a-b/a)]
!  *
!  *  erfc(z) = 1 - erf(z)
!  *                              2
!  *  erf(x+i*y) = erf(x) + exp(-x )/(2pi*x)*[(1-cos(2xy))+i*sin(2xy)]
!  *
!  *                                   2  inf       2      2   2
!  *                      + 2/pi*exp(-x ) sum exp(-n /4)/(n +4x )
!  *                                      n=1
!  *
!  *                        * [fn(x,y)+i*gn(x,y)] + epsi(x,y)
!  *
!  *  where, fn(x,y) = 2x-2x*cosh(n*y)*cos(2xy)+n*sinh(ny)*sin(2xy)
!  *
!  *         gn(x,y) = 2x*cosh(n*y)*sin(2xy)+n*sinh(ny)*cos(2xy)
!  *                   .
!  *         epsi(x,y) = 1.0d-16 * |erf(x+i*y)|
!  *
!  *  assuming r and b <> 0.0d0.
!  *
!  *  If b=0.0d0, using:
!  *
!  *         (2*l+1)*I(l) = r*r*I(l+1) - c(l) ,
!  *
!  *         c(l) = a**(2*l+1) * exp(-a**2*r**2) / r ,
!  *
!  *         I(0) = sqrt(pi)/(2*r*r) * erfc(r*a)
!  *
!  *  if r=0.0d0, return with dqint=0.0d0.
!  *
!  *  Ref: Chapter 7, "Handbook of Mathematical Functions"
!  *       Edited by Milton Abramowitz and Irene A. Stegun
!  *
!  *  written by Yang Wang, April, 1992.
!  *  modified by Yang Wang, April, 1999.
!  *******************************************************************
!
   dqint(0:lmax)=czero
!
   if (rij.le.ten2m6) then
      return
   end if
!
   a   = sqrt(eta*rij)*half
   b   = sqrtm1*sqrt(e*rij)*half
   r   = sqrt(rij)
   efac=exp(-rij*rij*eta*fourth+e/eta)/r
!
   x1  = rij*sqrt(eta)*half+real(sqrtm1*sqrt(e),RealKind)/sqrt(eta)
   y1  = real(sqrt(e),RealKind)/sqrt(eta)
   x2  = rij*sqrt(eta)*half-real(sqrtm1*sqrt(e),RealKind)/sqrt(eta)
   y2  =-real(sqrt(e),RealKind)/sqrt(eta)
!
   if (abs(b) .lt. ten2m6) then
      dqint(0) = sqrt(pi)*half*erfc(x1)/rij
      dqm1     = czero
   else
      xy12     = two*x1*y1
      xy22     = two*x2*y2
      t1       = czero
      t2       = czero
      do n=1,10        ! n up to 10 is plenty for the summation.
         n2 = n*n
         fn1=two*x1-two*x1*cosh(n*y1)*cos(two*x1*y1)+n*sinh(n*y1)*sin(two*x1*y1)
         gn1=two*x1*cosh(n*y1)*sin(two*x1*y1)+n*sinh(n*y1)*cos(two*x1*y1)
         t1 =t1+exp(-n2/four)*cmplx(fn1,gn1,CmplxKind)/(n2+four*x1*x1)
         fn2=two*x2-two*x2*cosh(n*y2)*cos(two*x2*y2)+n*sinh(n*y2)*sin(two*x2*y2)
         gn2=two*x2*cosh(n*y2)*sin(two*x2*y2)+n*sinh(n*y2)*cos(two*x2*y2)
         t2 =t2+exp(-n2/four)*cmplx(fn2,gn2,CmplxKind)/(n2+four*x2*x2)
      end do
      remain1 = exp(-x1*x1)/pi2*(cmplx(one-cos(xy12),sin(xy12),CmplxKind)/x1 &
                                 + four*t1)
      remain2 = exp(-x2*x2)/pi2*(cmplx(one-cos(xy22),sin(xy22),CmplxKind)/x2 &
                                 + four*t2)
      dqint(0)= sqrt(pi)/(four*rij)*( exp( two*r*b)*(erfc(x1)-remain1)  &
                                     +exp(-two*r*b)*(erfc(x2)-remain2) )
      dqm1    = sqrt(pi)/(four*b*r)*(-exp( two*r*b)*(erfc(x1)-remain1)  &
                                     +exp(-two*r*b)*(erfc(x2)-remain2) )
   end if
!
   if(lmax.eq.0) then
      return
   endif
!
   bfac = e*rij*half         !  = -2*b*b
   rfac = half/rij           !  = 1/(2*r*r)
   dqint(1)= (dqint(0)-bfac*dqm1+a*efac)*rfac
   do l=1,lmax-1
      dqint(l+1)= ((2*l+1)*dqint(l)-bfac*dqint(l-1)+a**(2*l+1)*efac)*rfac
   end do
!
   if (stop_routine.eq.sname) then
      call StopHandler(sname)
   end if
!
   end subroutine intfac
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine knterm(ksame,n,kvec,energy)
!  ===================================================================
   implicit   none
!
   character (len= 6), parameter :: sname='knterm'
!
   logical, intent(in) :: ksame
!
   integer (kind=IntKind), intent(in) :: n
!
   integer (kind=IntKind) :: jkn
   integer (kind=IntKind) :: kl
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: np_kfac
!
   real (kind=RealKind), intent(in) :: kvec(3)
!
   real (kind=RealKind) :: xk(3)
   real (kind=RealKind) :: xij
   real (kind=RealKind) :: xnk
   real (kind=RealKind) :: xnk2
!
   complex (kind=CmplxKind), intent(in) :: energy
!
   complex (kind=CmplxKind), pointer :: ylmc(:)
   complex (kind=CmplxKind) :: fpol
!
   logical :: redo_kfac
!
!  *******************************************************************
!  * calculate the D1 terms:
!  *
!  *              _  _ l        _  _ 2          _  _ 2
!  * dterm = sum |Kn+k| * exp[-(Kn+k) /eta] / (|Kn+k| - e)
!  *          Kn
!  *              m *_  _            _  _  _
!  *           * Y  (Kn+k) * exp[-i*(Kn+k)*aij]
!  *              l
!  *
!  *             _  _ l        _  _ 2                              m *
!  * kfac_sav = |Kn+k| * exp[-(Kn+k) /eta] * exp[-i*(Kn+k)*aij] * Y  (Kn+k)
!  *                                                               l
!  *******************************************************************
!
   if (kfac_allocated_full .or. n == 1) then
      np_kfac = n
      redo_kfac = .false.
   else
      np_kfac = 2
      redo_kfac = .true.
   endif
!
   if(.not.ksame .or. .not.kst_done(n) .or.redo_kfac) then 
!     ================================================================
!     kylm_sav depends only on k
!     if kvec is the same as the previous value, kylm_sav needs not 
!     to be re-calculated.
!     ================================================================
      ylmc=>dterm
      do jkn=1,nknlat
         xk(1)=knlat_x(jkn)+kvec(1)
         xk(2)=knlat_y(jkn)+kvec(2)
         xk(3)=knlat_z(jkn)+kvec(3)
!        =============================================================
!        call calYlmConjg to obtain the complex conjugate of spherical 
!        harmonics.
!        -------------------------------------------------------------
         call calYlmConjg(xk,lmax_dlm,ylmc)
!        -------------------------------------------------------------
         xij=xk(1)*aij(1,n)+xk(2)*aij(2,n)+xk(3)*aij(3,n)
         xnk2=xk(1)*xk(1)+xk(2)*xk(2)+xk(3)*xk(3)
         xnk=sqrt(xnk2)
         fac(0)=exp(-sqrtm1*xij-xnk2/eta)
         do l=1,lmax_dlm
            fac(l)=xnk**l*fac(0)
         enddo
         do kl=1,kmax_dlm
            kfac_sav(kl,jkn,np_kfac)=fac(lofk(kl))*ylmc(kl)
         enddo
      end do
!     ----------------------------------------------------------------
      nullify(ylmc)
!     ----------------------------------------------------------------
      kst_done(n) = .true.
   endif
!
   dterm(1:kmax_dlm)=czero
   do jkn=nknlat,1,-1
      xk(1)=knlat_x(jkn)+kvec(1)
      xk(2)=knlat_y(jkn)+kvec(2)
      xk(3)=knlat_z(jkn)+kvec(3)
      fpol=xk(1)*xk(1)+xk(2)*xk(2)+xk(3)*xk(3)-energy
!     ================================================================
!     Warning: fpol is related to the free electron poles.
!     ================================================================
      if(abs(fpol).lt.ten2m6) then
!        -------------------------------------------------------------
         call WarningHandler(sname,'free electron pole is found',energy)
!        -------------------------------------------------------------
         fpol=ten2m6
      endif
!     ================================================================
      dterm(1:kmax_dlm)=dterm(1:kmax_dlm)+kfac_sav(1:kmax_dlm,jkn,np_kfac)/fpol
!     ----------------------------------------------------------------
!     call zaxpy(kmax_dlm,cone/fpol,kfac_sav(1,jkn,np_kfac),1,dterm,1)
!     ----------------------------------------------------------------
   end do
!
   if (print_level.ge.2) then
      write(6,'(''    l    m                   d1l(k,e)'')')
      do kl=1,kmax_dlm
         write(6,'(2i5,2x,2d20.13)')lofk(kl),mofk(kl),dterm(kl)
      end do
   end if
!
   if (stop_routine.eq.sname) then
      call StopHandler(sname)
   end if
!
   end subroutine knterm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine rdlke(esame,n,j0,kvec,kappa)
!  ===================================================================
   implicit   none
!
   character (len= 5), parameter  :: sname='rdlke'
!
   logical, intent(in) :: esame
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: j0
!
   integer (kind=IntKind) :: jrs
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: kl
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: np_rylm, np_rfac
!
   real (kind=RealKind), intent(in) :: kvec(3)
!
   real (kind=RealKind) :: rsij(3)
   real (kind=RealKind) :: rsk
!
   complex (kind=CmplxKind), intent(in) :: kappa
!
   complex (kind=CmplxKind) :: xfac
   complex (kind=CmplxKind) :: z
   complex (kind=CmplxKind) :: zinv
!
   logical :: redo_rfac
!
!  *******************************************************************
!  *  calculates Dlm's in using real space sum........................
!  *
!  *                    - -
!  *                  i*k*Rn   -l+1       -  -        m *-  - 
!  *  Dlm = -p * sum e      * i   * h (p*|Rn+aij|) * Y  (Rn+aij)
!  *              n                  l                l
!  *
!  *        -i * p/sqrt(4*pi)      <--- this term exists if l = 0, and
!  *                                                        aij = 0.
!  *  where p=sqrt(E)
!  *******************************************************************
!
   if (rylm_allocated_full .or. n == 1) then
      np_rylm = n
   else
      np_rylm = 2
      do jrs=1,nrslat
         rsij(1)=rslat_x(jrs)+aij(1,n)
         rsij(2)=rslat_y(jrs)+aij(2,n)
         rsij(3)=rslat_z(jrs)+aij(3,n)
!        -------------------------------------------------------------
         call calYlmConjg(rsij,lmax_max2,                             &
                          rylm_sav(1:kmax_max2,jrs,np_rylm))
!        -------------------------------------------------------------
      enddo
   endif
!
   if(rfac_allocated_full .or. n == 1) then 
      np_rfac = n
      redo_rfac = .false.
   else
      np_rfac = 2
      redo_rfac = .true.
   endif
!
   if (.not.esame .or. .not.est_done(n) .or. redo_rfac) then
!     ================================================================
!     rfac_sav depends only on energy
!     if energy is the same as the previous value, rfac_sav needs
!     not to be re-calculated.
!     ================================================================
      k=lmax_dlm-1
      do jrs=nrslat,j0,-1
         rsij(1)=rslat_x(jrs)+aij(1,n)
         rsij(2)=rslat_y(jrs)+aij(2,n)
         rsij(3)=rslat_z(jrs)+aij(3,n)
!        =============================================================
!        calculate the hankel function h (z), l = 0, 1, ..., 2*lmax.
!                                       l    
!        rfac_sav = h (z)
!                    l
!        =============================================================
         z=kappa*sqrt(rsij(1)**2+rsij(2)**2+rsij(3)**2)
         zinv=cone/z
         rfac_sav(0,jrs,np_rfac)=-sqrtm1*exp(sqrtm1*z)*zinv
         if(k.ge.1) then
            rfac_sav(1,jrs,np_rfac)=-sqrtm1*rfac_sav(0,jrs,np_rfac)*(cone+sqrtm1*zinv)
            do l=1,k
               rfac_sav(l+1,jrs,np_rfac)=(2*l+1)*rfac_sav(l,jrs,np_rfac)*zinv-  &
                                         rfac_sav(l-1,jrs,np_rfac)
            enddo
         endif
      enddo
      est_done(n) = .true.
   endif
!
   dlke(1:kmax_dlm)=czero
   if(n == 1) then
      dlke(1)=cone/sqrt(pi4)
   endif
!
   do jrs=nrslat,j0,-1
      rsk = rslat_x(jrs)*kvec(1)+rslat_y(jrs)*kvec(2)+rslat_z(jrs)*kvec(3)
      xfac=exp(sqrtm1*rsk)
      do kl=1,kmax_dlm
         dlke(kl)=dlke(kl)+xfac*rfac_sav(lofk(kl),jrs,np_rfac)*rylm_sav(kl,jrs,np_rylm)
      enddo
   enddo
!
   do kl=1,kmax_dlm
      dlke(kl)=kappa*dlke(kl)/ilp1(lofk(kl))
!     dlke(kl)=-kappa*dlke(kl)*ilp1(lofk(kl))
   enddo
   if(sname.eq.stop_routine) then
      call StopHandler(sname)
   endif
!
   end subroutine rdlke
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine rsterm(esame,n,j0,kvec,energy)
!  ===================================================================
   implicit   none
!
   character (len= 6), parameter :: sname='rsterm'
!
   logical, intent(in) :: esame
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: j0
!
   integer (kind=IntKind) :: jrs
   integer (kind=IntKind) :: kl
   integer (kind=IntKind) :: np_rylm, np_rfac
!
   real (kind=RealKind), intent(in) :: kvec(3)
!
   real (kind=RealKind) :: rsij(3)
   real (kind=RealKind) :: rij
   real (kind=RealKind) :: rsk
!
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind) :: erk
!
   logical :: redo_rfac
!
!  *******************************************************************
!  *  calculate the D2 terms:
!  *
!  *                   ->   ->  l       ->  ->     m * -> ->
!  *  dterm(l) = sum' |Rs + aij| * exp(ik * Rs) * Y   (Rs+aij)
!  *             Rs                                l
!  *
!  *             inf        2*l        2   -> ->  2          2
!  *           * int {dx * x   * exp[-x * (Rs+aij) + 0.25*e/x ]}
!  *              a    
!  *
!  *  a = sqrt(eta)/2
!  *                                      -> ->
!  *  The summation excludes the case of |Rs+aij| = 0.
!  *
!  *
!  *              m * -> ->
!  *  rylm_sav = Y   (Rs+aij)
!  *              l
!  *
!  *              -> ->  l  inf        2*l        2   -> ->  2          2
!  *  rfac_sav = |Rs+aij| * int {dx * x   * exp[-x * |Rs+aij| + 0.25*e/x ]
!  *                         a
!  *
!  *              -> ->  -1/2   inf        2*l       -> ->      2     e  
!  *           = |Rs+aij|     * int {dx * x   * exp[|Rs+aij|*(-x + ------- )]
!  *                             b                                  4*x*x
!  *
!  *         -> ->  1/2
!  *  b = a*|Rs+aij|
!  *******************************************************************
!
   if (rylm_allocated_full .or. n == 1) then
      np_rylm = n
   else
      np_rylm = 2
      do jrs=1,nrslat
         rsij(1)=rslat_x(jrs)+aij(1,n)
         rsij(2)=rslat_y(jrs)+aij(2,n)
         rsij(3)=rslat_z(jrs)+aij(3,n)
!        -------------------------------------------------------------
         call calYlmConjg(rsij,lmax_max2,                             &
                          rylm_sav(1:kmax_max2,jrs,np_rylm))
!        -------------------------------------------------------------
      enddo
   endif
!
   if(rfac_allocated_full .or. n == 1) then 
      np_rfac = n
      redo_rfac = .false.
   else
      np_rfac = 2
      redo_rfac = .true.
   endif
!
   if(.not.esame .or. .not.est_done(n) .or. redo_rfac) then 
!     ================================================================
!     rfac_sav depends only on energy
!     if energy is the same as the previous value, rfac_sav needs 
!     not to be re-calculated.
!     ================================================================
      do jrs=j0,nrslat
         rsij(1)=rslat_x(jrs)+aij(1,n)
         rsij(2)=rslat_y(jrs)+aij(2,n)
         rsij(3)=rslat_z(jrs)+aij(3,n)
         rij = sqrt(rsij(1)**2+rsij(2)**2+rsij(3)**2)
!        =============================================================
!        call intfac to calculate the integral and store in rfac_sav:
!
!         -> ->  l  inf       2*l             -> ->  2
!        |Rs+aij| * int dx * x   * exp[-x**2*|Rs+aij| +e/(4*x*x)]
!                   a
!
!        a=0.5*sqrt(eta)
!        -------------------------------------------------------------
         call intfac(rij,energy,eta,lmax_dlm,rfac_sav(0:lmax_dlm,jrs,np_rfac))
!        -------------------------------------------------------------
      enddo
      est_done(n) = .true.
   endif
!
   dterm(1:kmax_dlm)=czero
   do jrs=nrslat,j0,-1
      rsk = rslat_x(jrs)*kvec(1)+rslat_y(jrs)*kvec(2)+rslat_z(jrs)*kvec(3)
      erk=exp(sqrtm1*rsk)
      do kl=1,kmax_dlm
         dterm(kl)=dterm(kl)+erk*rfac_sav(lofk(kl),jrs,np_rfac)*rylm_sav(kl,jrs,np_rylm)
      enddo
   enddo
!
   if (print_level.ge.2) then
      write(6,'(''    l    m                   dterm(k,e)'')')
      do kl=1,kmax_dlm
         write(6,'(2i5,2x,2d20.13)')lofk(kl),mofk(kl),dterm(kl)
      end do
   end if
!
   if (stop_routine.eq.sname) then
      call StopHandler(sname)
   end if
!
   end subroutine rsterm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function rsfunc(x,eta0,l)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: lp
!
   real (kind=RealKind) :: rsfunc
   real (kind=RealKind) :: x
   real (kind=RealKind) :: eta0
!
   if(x.lt.ten2m8) then
      rsfunc=1.0d+14
      return
   endif
!
!  -------------------------------------------------------------------
   call intfac(x,czero,eta0,l,facr(0))
!  -------------------------------------------------------------------
   rsfunc=abs(facr(0))
   do lp=1,l
      rsfunc=max(abs(facr(lp)),rsfunc)
   enddo
!
   end function rsfunc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function knfunc(x,eta0,l)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind) :: l
!
   real (kind=RealKind) :: knfunc
   real (kind=RealKind) :: x
   real (kind=RealKind) :: eta0
!
   real (kind=RealKind) :: x2
!
   if(x.lt.ten2m8) then
      knfunc=1.0d+14
      return
   endif
!
   x2=x*x
   if(x.lt.one) then
      knfunc=exp(-x2/eta0)/x2
   else
      knfunc=x**l*exp(-x2/eta0)/x2
   endif
!
   end function knfunc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFreeElectronPoleSum(k_in,energy_in) result(fep)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: nm1
   integer (kind=IntKind) :: nm2
   integer (kind=IntKind) :: nm3
   integer (kind=IntKind) :: ig
!
   real (kind=RealKind) :: vbrak(3,3)
   real (kind=RealKind) :: vol
   real (kind=RealKind) :: kncut, kn2cut
   real (kind=RealKind) :: factor, eimag, fr, fi
!
   real (kind=RealKind), intent(in) :: k_in(3)
!
   complex (kind=CmplxKind), intent(in) :: energy_in
   complex (kind=CmplxKind)             :: e_in
!
   complex (kind=CmplxKind) :: ksq
   complex (kind=CmplxKind) :: fep
!
   real(kind=RealKind),parameter :: eps=ten2m8
!
   eimag = real(-sqrtm1*energy_in)
   if (eimag.lt.eps) then
      e_in=cmplx(real(energy_in,kind=RealKind),eps,Kind=RealKind)
   else
      e_in=energy_in
   endif
   
   if (nknlat_fep == 0) then
      nknlat_fep = nknlat
      ipknlat_fep = nknlat_fep
      allocate(knlat_x_fep(1:nknlat_fep),knlat_y_fep(1:nknlat_fep),   &
               knlat_z_fep(1:nknlat_fep),knlatsq_fep(1:nknlat_fep))
      knlat_x_fep(1:nknlat_fep) = knlat_x(1:nknlat_fep)/ScalingFactor
      knlat_y_fep(1:nknlat_fep) = knlat_y(1:nknlat_fep)/ScalingFactor
      knlat_z_fep(1:nknlat_fep) = knlat_z(1:nknlat_fep)/ScalingFactor
      knlatsq_fep(1:nknlat_fep) = knlatsq(1:nknlat_fep)/ScalingFactor**2
   endif
!
   if ( knlatsq_fep(nknlat_fep) < real(e_in,kind=RealKind) ) then
      deallocate(knlat_x_fep,knlat_y_fep,knlat_z_fep,knlatsq_fep)
!
      vol=(Bravais(2,1)*Bravais(3,2)-Bravais(3,1)*Bravais(2,2))*Bravais(1,3)+   &
          (Bravais(3,1)*Bravais(1,2)-Bravais(1,1)*Bravais(3,2))*Bravais(2,3)+   &
          (Bravais(1,1)*Bravais(2,2)-Bravais(2,1)*Bravais(1,2))*Bravais(3,3)
!
!     ================================================================
!     generate the basis vectors for reciprocal space.
!     ================================================================
      factor=TWO*PI/vol
      vbrak(1,1)=factor*(Bravais(2,2)*Bravais(3,3)-Bravais(3,2)*Bravais(2,3))
      vbrak(2,1)=factor*(Bravais(3,2)*Bravais(1,3)-Bravais(1,2)*Bravais(3,3))
      vbrak(3,1)=factor*(Bravais(1,2)*Bravais(2,3)-Bravais(2,2)*Bravais(1,3))
      vbrak(1,2)=factor*(Bravais(2,3)*Bravais(3,1)-Bravais(3,3)*Bravais(2,1))
      vbrak(2,2)=factor*(Bravais(3,3)*Bravais(1,1)-Bravais(1,3)*Bravais(3,1))
      vbrak(3,2)=factor*(Bravais(1,3)*Bravais(2,1)-Bravais(2,3)*Bravais(1,1))
      vbrak(1,3)=factor*(Bravais(2,1)*Bravais(3,2)-Bravais(3,1)*Bravais(2,2))
      vbrak(2,3)=factor*(Bravais(3,1)*Bravais(1,2)-Bravais(1,1)*Bravais(3,2))
      vbrak(3,3)=factor*(Bravais(1,1)*Bravais(2,2)-Bravais(2,1)*Bravais(1,2))
!
      vol=abs(vol)
!
      kn2cut = real(e_in,kind=RealKind) + 0.01 
      nm1 = ceiling(kn2cut/(vbrak(1,1)**2+vbrak(2,1)**2+vbrak(3,1)**2))
      nm2 = ceiling(kn2cut/(vbrak(1,2)**2+vbrak(2,2)**2+vbrak(3,2)**2))
      nm3 = ceiling(kn2cut/(vbrak(1,3)**2+vbrak(2,3)**2+vbrak(3,3)**2))
      kncut = sqrt(kn2cut)
!
!     ================================================================
!     generate the reciprocal space lattice vectors.
!     ----------------------------------------------------------------
      call numlat(vbrak,kncut,nm1,nm2,nm3,ipknlat_fep)
!     ----------------------------------------------------------------
      allocate(knlat_x_fep(1:ipknlat_fep),knlat_y_fep(1:ipknlat_fep), &
               knlat_z_fep(1:ipknlat_fep),knlatsq_fep(1:ipknlat_fep))
!     ----------------------------------------------------------------
      call lattice(vbrak,kncut,nm1,nm2,nm3,                           &
                   knlat_x_fep,knlat_y_fep,knlat_z_fep,knlatsq_fep,   &
                   nknlat_fep, ipknlat_fep)
!     ----------------------------------------------------------------
      if(print_level.ge.0) then
         write(6,'(/,a)')                                             &
  ' A new set of reciprocal space lattice vectors is generated for Loyd formula'
         write(6,'(  ''            nm1, nm2, nm3 = '',3i5)') nm1,nm2,nm3
         write(6,'(  ''            Kn cut radius = '',1f10.5)') kncut
         write(6,'(  ''            Number of Kn  = '',i5)') nknlat_fep
         write(6,'(/)')
      endif
   endif
!
   fep = CZERO
   do ig = 1, nknlat_fep
      ksq = (k_in(1)+knlat_x_fep(ig))**2 + (k_in(2)+knlat_y_fep(ig))**2 + &
            (k_in(3)+knlat_z_fep(ig))**2
      fep = fep + log(ksq-e_in)
   enddo
   fep = -fep/PI
!
   if (eimag.lt.eps) then
      fr = real(fep,kind=RealKind)
      fi = real(-sqrtm1*fep,kind=RealKind)
      fi = real(nint(fi),kind=RealKind)
      fep=cmplx(fr,fi,kind=CmplxKind)
   endif
!
   end function getFreeElectronPoleSum
!  ===================================================================
end module StrConstModule
