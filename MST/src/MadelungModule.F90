! ********************************************************************
! *  MadelungModule                                                  *
! *  ==============                                                  *
! *     Purpose: compute Madelung matrix                             *
! *                                                                  *
! *     Public Functions:                                            *
! *            initMadelung -- initialize the module                 *
! *            endMadelung  -- clean up the arrays allocated in the  *
! *                            module                                *
! *       getMadelungMatrix -- returns M_{0,0}(1:Na,i)               *
! *             getDLMatrix -- returns alat^l * DL(1:Na,1:jmax,i)    *
! *             getDLFactor -- returns the prefactor of DL matrix.   *
! *                            The Madelung matrix is a product of   *
! *                            this factor and DL matrix.            *
! *                                                                  *
! ********************************************************************

module MadelungModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO, SQRTm1, TEN2m6, TEN2m8,  &
                               HALF, ONE, TWO, FOUR, PI, PI2, PI4, Y0inv
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler,      &
                                  StopHandler
   use IntegerFactorsModule, only : lofk, mofk, lofj, mofj, kofj
!
public :: initMadelung,         &
          endMadelung,          &
          printMadelungMatrix,  &
          getMadelungMatrix,    &
          getDLMatrix,          &
          getDLFactor,          &
          getRSpaceLattice,     &
          getKSpaceLattice,     &
          getNumRSVectors,      &
          getNumKSVectors,      &
          getJelliumPotential
!
private
   logical :: Initialized = .false.
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: PrintLevel
   integer (kind=IntKind) :: nrslat, nknlat
!
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
!
   real (kind=RealKind), allocatable, target :: madmat(:,:)
   real (kind=RealKind), allocatable, target :: DL_factor(:,:)
   real (kind=RealKind) :: a0, alat
!
   real (kind=RealKind), allocatable, target :: rspace(:)
   real (kind=RealKind), allocatable, target :: rslat(:,:)
   real (kind=RealKind), allocatable, target :: rslatsq(:)
   real (kind=RealKind), pointer :: rslat_x(:)
   real (kind=RealKind), pointer :: rslat_y(:)
   real (kind=RealKind), pointer :: rslat_z(:)
   real (kind=RealKind), allocatable, target :: knlat(:,:)
   real (kind=RealKind), allocatable, target :: knlatsq(:)
   real (kind=RealKind), pointer :: knlat_x(:)
   real (kind=RealKind), pointer :: knlat_y(:)
   real (kind=RealKind), pointer :: knlat_z(:)
   real (kind=RealKind), allocatable :: atom_position(:,:)
   real (kind=RealKind) :: eta, omegbra
!
   complex (kind=CmplxKind), allocatable, target :: DL_matrix(:,:,:)
   complex (kind=CmplxKind), allocatable, target :: cspace(:)
   complex (kind=CmplxKind), allocatable :: Ylm(:), ctmp(:)
!
   integer (kind=IntKind) :: lmax_rho = 0
   integer (kind=IntKind) :: jmax_rho = 1
   integer (kind=IntKind) :: kmax_rho = 1
   integer (kind=IntKind) :: lmax_pot = 0
   integer (kind=IntKind) :: jmax_pot = 1
   integer (kind=IntKind) :: kmax_pot = 1
   integer (kind=IntKind) :: lmax_mad = 0
   integer (kind=IntKind) :: jmax_mad = 1
   integer (kind=IntKind) :: kmax_mad = 1
!
contains
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMadelung(num_local_atoms,num_atoms,gindex,         &
                           lmax_rho_in,lmax_pot_in,bravais,posi,iprint)
!  ==================================================================
   use IntegerFactorsModule, only : initIntegerFactors, endIntegerFactors
!
   use MathParamModule, only : THIRD, THREE
!
   use GauntFactorsModule, only : isGauntInitialized => isInitialized
   use GauntFactorsModule, only : initGauntFactors
   use GauntFactorsModule, only : getGauntFactor
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: num_local_atoms
   integer (kind=IntKind), intent(in) :: gindex(num_local_atoms)
   integer (kind=IntKind), intent(in) :: iprint
!
   integer (kind=IntKind), intent(in) :: lmax_rho_in
   integer (kind=IntKind), intent(in) :: lmax_pot_in
!
   integer (kind=IntKind) :: i, n1, l, m, kl, jl, jl_pot
   integer (kind=IntKind) :: l_pot, kl_pot, l_rho, kl_rho
   integer (kind=IntKind) :: l_sum, m_dif
!
   real (kind=RealKind), intent(in) :: bravais(3,3)
   real (kind=RealKind), target, intent(in) :: posi(3,num_atoms)
!
   real (kind=RealKind), pointer :: factmat(:)
   real (kind=RealKind) :: gaunt
   real (kind=RealKind) :: sfac
   real (kind=RealKind) :: factor
   real (kind=RealKind) :: vatom
   real (kind=RealKind) :: vbrar(3,3)
   real (kind=RealKind) :: vbrak(3,3)
!
   if (num_local_atoms < 1) then
      call ErrorHandler('initMadelung','invalid local number of atoms', &
                        num_local_atoms)
   else if (num_atoms < 1) then
      call ErrorHandler('initMadelung','invalid total number of atoms', &
                        num_atoms)
   else if (lmax_rho_in < 0) then
      call ErrorHandler('initMadelung','invalid 1st lmax value',lmax_rho_in)
   else if (lmax_pot_in < 0) then
      call ErrorHandler('initMadelung','invalid 2nd lmax value',lmax_pot_in)
   endif
!
   LocalNumAtoms = num_local_atoms
   GlobalNumAtoms = num_atoms
!
   lmax_pot = lmax_pot_in
   jmax_pot = (lmax_pot+1)*(lmax_pot+2)/2
   kmax_pot = (lmax_pot+1)*(lmax_pot+1)
   lmax_rho = lmax_rho_in
   jmax_rho = (lmax_rho+1)*(lmax_rho+2)/2
   kmax_rho = (lmax_rho+1)*(lmax_rho+1)
   lmax_mad = 2*max(lmax_pot_in,lmax_rho_in)
   jmax_mad = (lmax_mad+1)*(lmax_mad+2)/2
   kmax_mad = (lmax_mad+1)*(lmax_mad+1)
!
   allocate(GlobalIndex(LocalNumAtoms))
   allocate(atom_position(3,GlobalNumAtoms))
   allocate(madmat(GlobalNumAtoms,LocalNumAtoms))
   allocate(rspace(1:lmax_mad+1))
!
   do i = 1, LocalNumAtoms
      GlobalIndex(i) = gindex(i)
   enddo
!
   PrintLevel = iprint
!
   atom_position = posi
!
!  -------------------------------------------------------------------
   call calScalingFactor(bravais,sfac)
!  -------------------------------------------------------------------
   a0 = sfac
!
   vbrar(1:3,1:3) = bravais(1:3,1:3)
!
!  ===================================================================
!  change units so that both vbrar and atom_posi_* are in
!  in the units of a0
!  -------------------------------------------------------------------
   call dscal(9,ONE/a0,vbrar,1)
   call dscal(GlobalNumAtoms,ONE/a0,atom_position(1,1),3)
   call dscal(GlobalNumAtoms,ONE/a0,atom_position(2,1),3)
   call dscal(GlobalNumAtoms,ONE/a0,atom_position(3,1),3)
!  -------------------------------------------------------------------
!
   omegbra=(vbrar(2,1)*vbrar(3,2)-vbrar(3,1)*vbrar(2,2))*vbrar(1,3)+ &
           (vbrar(3,1)*vbrar(1,2)-vbrar(1,1)*vbrar(3,2))*vbrar(2,3)+ &
           (vbrar(1,1)*vbrar(2,2)-vbrar(2,1)*vbrar(1,2))*vbrar(3,3)
   factor=PI2/omegbra
   omegbra=abs(omegbra)
   vatom = omegbra/real(GlobalNumAtoms,RealKind)
   alat=a0*(THREE*vatom/PI4)**THIRD
!
!  ===================================================================
!  generate basis vectors for reciprocal space....................
!  ===================================================================
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
!  ===================================================================
!  obtain the lattice vectors for the big cell.....................
!  rslat, rslatsq, knlat, and knlatsq are in the units of a0 ......
!  the first argument is l+1, which is 1 in this case.................
!  -------------------------------------------------------------------
   call genLattice(vbrar,vbrak)
!  -------------------------------------------------------------------
!
   if(PrintLevel > 1) then
      write(6,'(/)')
      write(6,*) "Madelung:: scaling factor: ", sfac
      write(6,'(12x,a)')                                              &
              '    n                    rslat                  rslatsq'
      write(6,'(12x,56(''=''))')
      write(6,'(12x,1i5,2x,4f12.5)')                                  &
            (n1,rslat_x(n1),rslat_y(n1),rslat_z(n1),rslatsq(n1),n1=1,nrslat)
      write(6,'(/)')
      write(6,'(12x,a)')                                              &
              '    n                    knlat                  knlatsq'
      write(6,'(12x,56(''=''))')
      write(6,'(12x,1i5,2x,4f12.5)')                                  &
            (n1,knlat_x(n1),knlat_y(n1),knlat_z(n1),knlatsq(n1),n1=1,nknlat)
   endif
!
!  ===================================================================
!  Note: rslatsq will be used as scramble space later.
!  ===================================================================
!
   if ( jmax_mad > 1 ) then
      if (.not.isGauntInitialized()) then
!        -------------------------------------------------------------
         call initGauntFactors(lmax_mad,'xxxx',iprint)
!        -------------------------------------------------------------
      endif
!
      call initIntegerFactors(lmax_mad)
!
      allocate(DL_matrix(GlobalNumAtoms,1:kmax_mad,LocalNumAtoms))
      allocate(DL_factor(kmax_mad,jmax_mad))
      allocate(Ylm(1:kmax_mad),cspace(1:kmax_mad), ctmp(1:kmax_mad))
      DL_matrix = CZERO
      DL_factor = CZERO
!
      factmat => rspace(1:lmax_mad+1)
!
      factmat(1) = ONE
      do l = 1, lmax_mad
         factmat(l+1) = factmat(l)/(2*l+ONE)
      enddo
!
      do jl_pot = 1, jmax_mad
         l_pot = lofj(jl_pot)
         kl_pot = kofj(jl_pot)
         do kl_rho = 1, kmax_mad
            l_rho = lofk(kl_rho)
            l_sum = l_pot+l_rho
            m_dif = mofk(kl_rho)-mofj(jl_pot)
            kl = (l_sum+1)*(l_sum+1)-l_sum+m_dif
            gaunt = getGauntFactor(kl_pot,kl_rho,kl)
            DL_factor(kl_rho,jl_pot)= gaunt*factmat(l_pot+1)*factmat(l_rho+1)
         enddo
      enddo
      nullify(factmat)
   endif
!
   do i=1,LocalNumAtoms
!     ================================================================
!     set up the generalized Madelung matrix..........................
!     ----------------------------------------------------------------
      call madewd(i,GlobalIndex(i))
!     ----------------------------------------------------------------
   enddo
!
   if ( jmax_mad > 1 ) then
      deallocate(Ylm, cspace, ctmp)
   endif
!
   deallocate(GlobalIndex, rspace)
   deallocate(rslat, rslatsq)
   nullify(rslat_x, rslat_y, rslat_z)
   deallocate(knlat, knlatsq)
   nullify(knlat_x, knlat_y, knlat_z)
   if (jmax_mad > 1) then
      call endIntegerFactors()
   endif
!
   Initialized = .true.
!
   end subroutine initMadelung
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endMadelung()
!  ===================================================================
   implicit none
!
   Initialized = .false.
!
   deallocate( madmat, atom_position )
   if (jmax_mad > 1) then
      deallocate( DL_matrix, DL_factor )
   endif
!
   end subroutine endMadelung
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calScalingFactor(Bravais,sfac)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: Bravais(3,3)
   real (kind=RealKind), intent(out) :: sfac
!
   integer (kind=IntKind) :: nm1, nm2, nm3, iter
   integer (kind=IntKind), parameter :: max_iter = 10000
!
   real (kind=RealKind) :: a1, a2, a3
   real (kind=RealKind) :: vbrar(3,3), vbrak(3,3)
   real (kind=RealKind) :: volr, vfac
   real (kind=RealKind) :: rscut, kncut
!
   real (kind=RealKind), parameter :: tfac = 1.2d0
   real (kind=RealKind), parameter :: fstep = 0.02d0
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
   sfac=min(a1,a2,a3)
   eta=HALF+0.1d0*max(a1,a2,a3)/sfac
!
   sfac = sfac/PI2
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
      vfac=PI2/abs(volr)
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
      call getrscut(vbrar(1:3,1),vbrar(1:3,2),vbrar(1:3,3),rscut,nm1,nm2,nm3)
      call numlat(vbrar,rscut,nm1,nm2,nm3,nrslat)
!     ----------------------------------------------------------------
!
!     ================================================================
!     calculate rscut, the radius of real space truncation sphere.
!     ----------------------------------------------------------------
      call getkncut(vbrak(1:3,1),vbrak(1:3,2),vbrak(1:3,3),kncut,nm1,nm2,nm3)
      call numlat(vbrak,kncut,nm1,nm2,nm3,nknlat)
!     ----------------------------------------------------------------
!     write(6,'(a,3i8)')'iter, nrslat, nknlat = ',iter,nrslat,nknlat
      if (iter > max_iter .or. sfac <= 0.1d0) then
         write(6,'(a,3i8)')'iter, nrslat, nknlat = ',iter,nrslat,nknlat
!        =============================================================
!        If this message shows up, reduce fstep value.
!        -------------------------------------------------------------
         call WarningHandler('calScalingFactor',                        &
                             'The scaling factor may not be optimal',sfac)
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
!
   end subroutine calScalingFactor
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMadelungMatrix(i) result(pm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
!
   real (kind=RealKind), pointer :: pm(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getMadelungMatrix','MadelungModule not initialized')
   else if (i < 1 .or. i > LocalNumAtoms) then
      call ErrorHandler('getMadelungMatrix','invalid local atom index',i)
   endif
!
   pm => madmat(1:GlobalNumAtoms,i)
!
   end function getMadelungMatrix
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDLMatrix(i,afac) result(dm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
!
   real (kind=RealKind), intent(out) :: afac
!
   complex (kind=CmplxKind), pointer :: dm(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getDLMatrix','MadelungModule not initialized')
   else if (i < 1 .or. i > LocalNumAtoms) then
      call ErrorHandler('getDLMatrix','invalid local atom index',i)
   else if (lmax_mad < 1) then
      call WarningHandler('getDLMatrix','invalid for lmax_mad < 1',lmax_mad)
      nullify( dm )
   else
      afac = alat
      dm => DL_matrix(1:GlobalNumAtoms,1:kmax_mad,i)
   endif
!
   end function getDLMatrix
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDLFactor(jl) result(dfac)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: jl
!
   real (kind=RealKind), pointer :: dfac(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getDLFactor','MadelungModule not initialized')
   else if (jl < 1) then
      call ErrorHandler('getDLFactor','jl < 1',jl)
   else if (jl > jmax_mad) then
      call ErrorHandler('getDLFactor','jl > jmax_mad',jl,jmax_mad)
   else if (jmax_mad < 1) then
      call WarningHandler('getDLFactor','invalid for lmax_mad < 1',lmax_mad)
      nullify( dfac )
   else
      dfac => DL_factor(1:kmax_mad,jl)
   endif
!
   end function getDLFactor
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumRSVectors() result(nrs)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: nrs
!
   if (.not.Initialized) then
      call ErrorHandler('getNumRSVectors','MadelungModule not initialized')
   endif
!
   nrs = nrslat
!
   end function getNumRSVectors
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumKSVectors() result(nkn)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: nkn
!
   if (.not.Initialized) then
      call ErrorHandler('getNumKSVectors','MadelungModule not initialized')
   endif
!
   nkn = nknlat
!
   end function getNumKSVectors
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRSpaceLattice(fac) result(rs)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(out) :: fac
   real (kind=RealKind), pointer :: rs(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getRSpaceLattice','MadelungModule not initialized')
   endif
!
   rs => rslat(1:nrslat,1:3)
   fac = a0
!
   end function getRSpaceLattice
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKSpaceLattice(fac) result(kn)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(out) :: fac
   real (kind=RealKind), pointer :: kn(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getKSpaceLattice','MadelungModule not initialized')
   endif
!
   kn => knlat(1:nknlat,1:3)
   fac = ONE/a0
!
   end function getKSpaceLattice
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printMadelungMatrix(iprint)
!  ===================================================================
   use Atom2ProcModule, only : getGlobalIndex
!
   implicit none
!
   integer (kind=IntKind), optional :: iprint
   integer (kind=IntKind) :: i, j, id, print_level
!
   real (kind=RealKind) :: sum
!
   if (.not.Initialized) then
      call WarningHandler('printMadelungMatrix',                      &
                          'MadelungModule is not initialized')
      return
   else if (present(iprint)) then
      print_level = iprint
   else
      print_level = 0
   endif
!
   write(6,'(/,a)')'----------------------------------------------------------'
   write(6,'( a )')'                         *******************************'
   write(6,'( a )')'                         *        PrintMadelung        *'
   write(6,'( a )')'                         *******************************'
   write(6,'( / )')
   do i=1,LocalNumAtoms
      id = getGlobalIndex(i)
      sum=ZERO
      do j=1,GlobalNumAtoms
         if (print_level > 0) then
            write(6,'(a,i5,a,d20.13 )')' Madelung Matrix(j,i) for atom ', &
                                       id,': ',madmat(j,i)
         endif
         sum=sum+madmat(j,i)
      enddo
      write(6,'(a,i5,a,d20.13 )')                                     &
              ' Sum over j of Madelung Matrix(j,i) for atom ',        &
              id,':  ',sum
   enddo
   write(6,'( a )')'----------------------------------------------------------'
!
   end subroutine printMadelungMatrix
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genLattice(vbrar,vbrak)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind) :: nm1, nm2, nm3
   integer (kind=IntKind) :: nr, nk, ipmax
!
   real (kind=RealKind), intent(in) :: vbrar(3,3)
   real (kind=RealKind), intent(in) :: vbrak(3,3)
!
   real (kind=RealKind), pointer :: vec_x(:)
   real (kind=RealKind), pointer :: vec_y(:)
   real (kind=RealKind), pointer :: vec_z(:)
   real (kind=RealKind), pointer :: vecsq(:)
!
   real (kind=RealKind) :: rscut
   real (kind=RealKind) :: kncut
!
!  *******************************************************************
!  Sets up real space Bravais lattice vectors
!  *******************************************************************
!
!  ===================================================================
!  calculate rscut, the radius of real space truncation sphere.....
!  -------------------------------------------------------------------
   call getrscut(vbrar(1:3,1),vbrar(1:3,2),vbrar(1:3,3),rscut,nm1,nm2,nm3)
   call numlat(vbrar,rscut,nm1,nm2,nm3,nr)
!  -------------------------------------------------------------------
   nrslat = nr
   allocate(rslat(1:nrslat,1:3), rslatsq(1:nrslat))
   rslat   = ZERO
   rslatsq = ZERO
   rslat_x => rslat(1:nrslat,1)
   rslat_y => rslat(1:nrslat,2)
   rslat_z => rslat(1:nrslat,3)
!
!  ===================================================================
!  generate the real space lattice vectors.........................
!  ===================================================================
   vec_x => rslat_x(1:nrslat)
   vec_y => rslat_y(1:nrslat)
   vec_z => rslat_z(1:nrslat)
   vecsq => rslatsq(1:nrslat)
!
   ipmax = nrslat
!  -------------------------------------------------------------------
   call lattice(vbrar,rscut,nm1,nm2,nm3,vec_x,vec_y,vec_z,vecsq,nr,ipmax)
!  -------------------------------------------------------------------
!
   if(PrintLevel.ge.0) then
      write(6,'(/,'' Real Space Lattice:: nm1,nm2,nm3   = '',3i5)')nm1,nm2,nm3
      write(6,'(  ''                      Rs cut radius = '',1f10.5)') rscut
      write(6,'(  ''                      Number of Rs  = '',i5)') nrslat
   endif
!
!  *******************************************************************
!  Sets up receprocal space Bravais lattice vectors
!  *******************************************************************
!
!  ===================================================================
!  calculate kncut, the radius of k-space truncation sphere........
!  -------------------------------------------------------------------
   call getkncut(vbrak(1:3,1),vbrak(1:3,2),vbrak(1:3,3),kncut,nm1,nm2,nm3)
   call numlat(vbrak,kncut,nm1,nm2,nm3,nr)
!  -------------------------------------------------------------------
   nknlat = nr
   allocate(knlat(1:nknlat,1:3), knlatsq(1:nknlat))
   knlat   = ZERO
   knlatsq = ZERO
   knlat_x => knlat(1:nknlat,1)
   knlat_y => knlat(1:nknlat,2)
   knlat_z => knlat(1:nknlat,3)
!
!  ===================================================================
!  generate the reciprocal space lattice vectors...................
!  ===================================================================
   vec_x => knlat_x(1:nknlat)
   vec_y => knlat_y(1:nknlat)
   vec_z => knlat_z(1:nknlat)
   vecsq => knlatsq(1:nknlat)
!
   ipmax = nknlat
!  -------------------------------------------------------------------
   call lattice(vbrak,kncut,nm1,nm2,nm3,vec_x,vec_y,vec_z,vecsq,nk,ipmax)
!  -------------------------------------------------------------------
!
   if(PrintLevel.ge.0) then
      write(6,'(/,'' Reciprocal Lattice:: nm1,nm2,nm3   = '',3i5)')nm1,nm2,nm3
      write(6,'(  ''                      Kn cut radius = '',1f10.5)') kncut
      write(6,'(  ''                      Number of Kn  = '',i5)') nknlat
   endif
!
   nullify(vec_x, vec_y, vec_z, vecsq)
!
   end subroutine genLattice
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getrscut(a1,a2,a3,rscut,nm1,nm2,nm3)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(out) :: nm1
   integer (kind=IntKind), intent(out) :: nm2
   integer (kind=IntKind), intent(out) :: nm3
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: k
!
   real (kind=RealKind), intent(in) :: a1(3)
   real (kind=RealKind), intent(in) :: a2(3)
   real (kind=RealKind), intent(in) :: a3(3)
   real (kind=RealKind), intent(out) :: rscut
   real (kind=RealKind) :: r(3)
   real (kind=RealKind) :: rm
   real (kind=RealKind) :: term
   real (kind=RealKind), pointer :: gamma_l(:)
   real (kind=RealKind), parameter :: epsi=1.0d-14
!
   gamma_l => rspace(1:lmax_mad+1)
!
!  ===================================================================
!  calculate nm1,nm2,nm3...........................................
!  ===================================================================
   r(1)=sqrt(a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3))
   term=ONE
   nm1=0
   do while(term.gt.HALF*epsi)
      nm1=nm1+1
      rm=nm1*r(1)
!     term=erfc(rm/eta)/rm**lp1
      call calGammaFunc(rm/eta,lmax_mad,gamma_l(1:lmax_mad+1))
      term=gamma_l(lmax_mad+1)/(rm/TWO)**(lmax_mad+1)
   enddo
!
   r(2)=sqrt(a2(1)*a2(1)+a2(2)*a2(2)+a2(3)*a2(3))
   term=ONE
   nm2=0
   do while(term.gt.HALF*epsi)
      nm2=nm2+1
      rm=nm2*r(2)
!     term=erfc(rm/eta)/rm**lp1
      call calGammaFunc(rm/eta,lmax_mad,gamma_l(1:lmax_mad+1))
      term=gamma_l(lmax_mad+1)/(rm/TWO)**(lmax_mad+1)
   enddo
!
   r(3)=sqrt(a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3))
   term=ONE
   nm3=0
   do while(term.gt.HALF*epsi)
      nm3=nm3+1
      rm=nm3*r(3)
!     term=erfc(rm/eta)/rm**lp1
      call calGammaFunc(rm/eta,lmax_mad,gamma_l(1:lmax_mad+1))
      term=gamma_l(lmax_mad+1)/(rm/TWO)**(lmax_mad+1)
   enddo
!
!  ===================================================================
!  calculate rscut.................................................
!  ===================================================================
   rscut=r(1)*nm1
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
            rscut=max(rscut,rm)
         enddo
      enddo
   enddo
!
   end subroutine getrscut
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getkncut(a1,a2,a3,kncut,nm1,nm2,nm3)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(out) :: nm1
   integer (kind=IntKind), intent(out) :: nm2
   integer (kind=IntKind), intent(out) :: nm3
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: k
!
   real (kind=RealKind), intent(in) :: a1(3)
   real (kind=RealKind), intent(in) :: a2(3)
   real (kind=RealKind), intent(in) :: a3(3)
   real (kind=RealKind), intent(out) :: kncut
   real (kind=RealKind) :: r(3)
   real (kind=RealKind) :: rm2
   real (kind=RealKind) :: term
   real (kind=RealKind) :: fac
   real (kind=RealKind), parameter :: epsi=1.0d-14
!
   fac=eta*eta/FOUR
!
!  ===================================================================
!  calculate nm1,nm2,nm3...........................................
!  ===================================================================
   r(1)=a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3)
   term=ONE
   nm1=0
   do while(term.gt.HALF*epsi)
      nm1=nm1+1
      rm2=nm1*nm1*r(1)
!     term=exp(-fac*rm2)/rm2
      term=exp(-fac*rm2)*(sqrt(rm2))**(lmax_mad-2)
   enddo
!
   r(2)=a2(1)*a2(1)+a2(2)*a2(2)+a2(3)*a2(3)
   term=ONE
   nm2=0
   do while(term.gt.HALF*epsi)
      nm2=nm2+1
      rm2=nm2*nm2*r(2)
!     term=exp(-fac*rm2)/rm2
      term=exp(-fac*rm2)*(sqrt(rm2))**(lmax_mad-2)
   enddo
!
   r(3)=a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3)
   term=ONE
   nm3=0
   do while(term.gt.HALF*epsi)
      nm3=nm3+1
      rm2=nm3*nm3*r(3)
!     term=exp(-fac*rm2)/rm2
      term=exp(-fac*rm2)*(sqrt(rm2))**(lmax_mad-2)
   enddo
!
!  ===================================================================
!  calculate kncut.................................................
!  ===================================================================
   kncut=sqrt(r(1))*nm1
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
            rm2=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
            kncut=max(kncut,rm2)
         enddo
      enddo
   enddo
!
   end subroutine getkncut
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine madewd(id,myatom)
!  ===================================================================
!
!  *******************************************************************
!     program for calculating the madelung constants
!
!     input :
!     =====
!             id = local atom index
!             myatom = global atom index
!
!  *******************************************************************
#ifdef DIRECT_SUM
   use SphericalHarmonicsModule, only : calYlm
#endif
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: myatom
!
   integer (kind=IntKind) :: ibegin
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: l, kl
!
   real (kind=RealKind) :: aij(3)
   real (kind=RealKind) :: r0tm
   real (kind=RealKind) :: term0
   real (kind=RealKind) :: term12
#ifdef DIRECT_SUM
   real (kind=RealKind) :: vec(3), vm
   real (kind=RealKind), pointer :: factmat(:)
   complex (kind=CmplxKind), pointer :: dirsum(:)
#endif
!
   complex (kind=CmplxKind), pointer :: dlm(:)
!
!  *******************************************************************
!  calculate Madelung constant matrix:
!
!     for i <> j,
!                                   2   2       -> ->
!            4*pi          1    -eta *Kq /4 - i*Kq*aij
!     M   =  ---- * sum  ----- e
!      ij    tau    q<>0    2
!                         Kq
!
!                          ->   ->                  2
!                 1 - erf(|Rn + aij|/eta)     pi*eta
!          + sum ------------------------- - ---------
!             n         ->   ->                 tau
!                      |Rn + aij|
!
!     for i = j,
!                                   2   2
!            4*pi          1    -eta *Kq /4
!     M   =  ---- * sum  ----- e
!      ii    tau    q<>0    2
!                         Kq
!
!                                             2
!                  1 - erf(Rn/eta)      pi*eta           2
!          + sum  ----------------- - ---------- - --------------
!            n<>0         Rn             tau        sqrt(pi)*eta
!
!     eta is the Ewald parameter;
!     tau, atom_posi_*, rslat_* and knlat_* are in the units of a0;  
!     madmat is in the units of a0=alat;
!  *******************************************************************
!
   term0=-PI*eta*eta/omegbra
!
!  ===================================================================
!  start madmat calculation
!  ===================================================================
   do n =1,GlobalNumAtoms
!     ================================================================
!     aij is in the units of a0 ...................................
!     ================================================================
      aij(1) = atom_position(1,myatom) - atom_position(1,n)
      aij(2) = atom_position(2,myatom) - atom_position(2,n)
      aij(3) = atom_position(3,myatom) - atom_position(3,n)
!
      if( n .eq. myatom ) then
         ibegin=2
         r0tm=-TWO/sqrt(PI)/eta
      else
         ibegin=1
         r0tm=ZERO
      endif
!     ================================================================
!     perform the reciprocal-space sum and real-space sum
!     ----------------------------------------------------------------
      call madsum(ibegin,aij,term12)
!     ----------------------------------------------------------------
      madmat(n,id)=term12+r0tm+term0
!     ================================================================
!     The unit of madmat is finally resumed to a0 = alat
!     ================================================================
      madmat(n,id)=madmat(n,id)/a0
!
      if (jmax_mad > 1) then
         DL_matrix(n,1,id) = madmat(n,id)*Y0inv
#ifdef DIRECT_SUM
         dirsum => cspace(1:kmax_mad)
         dirsum(:) = CZERO
         do i=nrslat,ibegin,-1
            vec(1) = aij(1) - rslat_x(i)
            vec(2) = aij(2) - rslat_y(i)
            vec(3) = aij(3) - rslat_z(i)
            vm = a0*sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
!           ----------------------------------------------------------
            call calYlm(vec,lmax_mad,Ylm)
!           ----------------------------------------------------------
            do kl = kmax_mad, 2, -1
               dirsum(kl) = dirsum(kl) + Ylm(kl)/vm**(lofk(kl)+1)
            enddo
         enddo
!
         factmat => rspace(1:lmax_mad+1)
         factmat(1) = ONE
         do l = 1, lmax_mad
            factmat(l+1) = factmat(l)*(2*l+ONE)
         enddo
!
         do kl = 2, kmax_mad
            DL_matrix(n,kl,id) =                                         &
                      dirsum(kl)*PI4*factmat(lofk(kl))*alat**(lofk(kl))
         enddo
!
         nullify(dirsum, factmat)
#else
         dlm => cspace(1:kmax_mad)
         dlm = CZERO
!        -------------------------------------------------------------
         call dlsum(ibegin,aij,dlm)
!        -------------------------------------------------------------
         do kl = 2, kmax_mad
            l = lofk(kl)
            DL_matrix(n,kl,id) = dlm(kl)*(alat/a0)**l/a0
         enddo
         nullify(dlm)
#endif
      endif
   enddo                                  ! end do loop over n
!
   end subroutine madewd
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine madsum(ibegin,aij,term12)
!  ===================================================================
!
!  *******************************************************************
!  performs ewald summation for Madelung constant.
!  requires 'lastkn,lastrs,eta,knlat,aij'
!  returns 'madelung sum': term12
!  *******************************************************************
   implicit   none
!
   integer (kind=IntKind), intent(in) :: ibegin
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: aij(3)
   real (kind=RealKind), intent(out) :: term12
!
   real (kind=RealKind), pointer :: rslatmd(:)
   real (kind=RealKind) :: fac, rterm
   real (kind=RealKind) :: erfc
!
   rslatmd => rslatsq(1:nrslat)
!
!  ===================================================================
!  subrtact aij from rslat and calculate rslatmd which is used in
!  calculating the real-space integral
!  rslatmd, and aij are in the units of a0 = 1
!  ===================================================================
   do i = 1,nrslat
      rslatmd(i)=sqrt( (rslat_x(i)-aij(1))*(rslat_x(i)-aij(1))        &
                      +(rslat_y(i)-aij(2))*(rslat_y(i)-aij(2))        &
                      +(rslat_z(i)-aij(3))*(rslat_z(i)-aij(3)) )
   enddo
!
!  ===================================================================
!                                     2   2       -> ->
!              4*pi          1    -eta *Kq /4 - i*Kq*aij
!     term1 =  ---- * sum  ----- e
!              tau    q<>0    2
!                           Kq
!
!     note sum starts at 2 since kn=0.0 of 1/(rs-aij) is
!     canceled by kn=0.0 of the 1/rs sum.
!  ===================================================================
   term12 = ZERO
   fac=eta*eta/FOUR
   do i = nknlat,2,-1
      term12=term12+exp(-fac*knlatsq(i))/knlatsq(i)                   &
                   *cos( knlat_x(i)*aij(1)+knlat_y(i)*aij(2)          &
                        +knlat_z(i)*aij(3) )
   enddo
   term12 = PI4/omegbra*term12
!  ===================================================================
!
!                           ->   ->
!                  1 - erf(|Rn + aij|/eta)
!     term2 = sum -------------------------
!              n         ->   ->
!                       |Rn + aij|
!
!     note for calculation of aij=0.0 term ibegin=2.
!  ===================================================================
   rterm = ZERO
   do i=nrslat,ibegin,-1
      rterm = rterm + erfc(rslatmd(i)/eta)/rslatmd(i)
   enddo
   term12 = term12 + rterm
!
   nullify(rslatmd)
!
   end subroutine madsum
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine dlsum(ibegin,aij,dlm)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ibegin
!
   integer (kind=IntKind) :: i, kl, l, lp1, m
!
   real (kind=RealKind), intent(in) :: aij(3)
!
   real (kind=RealKind) :: vec(3), vlen, vhalf, rfac
   real (kind=RealKind), pointer :: gamma_l(:)
   real (kind=RealKind) :: tfac, expfac, aij2, sintfac, costfac
!
   complex (kind=CmplxKind) :: cfac
   complex (kind=CmplxKind), intent(out) :: dlm(1:kmax_mad)
!
   gamma_l => rspace(1:lmax_mad+1)
   gamma_l = ZERO
!
   aij2 = aij(1)*aij(1) + aij(2)*aij(2) + aij(3)*aij(3)
!  ===================================================================
!  start DL_matrix calculation
!
!  vec and aij are in the units of a0 = 1
!  The unit of DL_matrix needs to be resumed to a0 = alat later
!  ===================================================================
   dlm(1:kmax_mad) = CZERO
   do i=nrslat, ibegin, -1
      vec(1) = aij(1) + rslat_x(i)
      vec(2) = aij(2) + rslat_y(i)
      vec(3) = aij(3) + rslat_z(i)
      vlen = sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
!     ----------------------------------------------------------------
      call calYlm(vec,lmax_mad,Ylm)
      call calGammaFunc(vlen/eta,lmax_mad,gamma_l(1:lmax_mad+1))
!     ----------------------------------------------------------------
      vhalf = HALF*vlen
      do kl = kmax_mad, 2, -1
         l = lofk(kl)
         lp1 = l+1
         dlm(kl) = dlm(kl) + gamma_l(lp1)*Ylm(kl)/vhalf**lp1
      enddo
   enddo
   rfac = FOUR*sqrt(PI)
   dlm(2:kmax_mad) = rfac*dlm(2:kmax_mad)
!
   rfac= -eta*eta/FOUR
   ctmp(1:kmax_mad) = CZERO
   do i = nknlat,2,-1
      vec(1) = knlat_x(i)
      vec(2) = knlat_y(i)
      vec(3) = knlat_z(i)
      vlen = sqrt(knlatsq(i))
!     ----------------------------------------------------------------
      call calYlm(vec,lmax_mad,Ylm)
!     ----------------------------------------------------------------
      expfac = exp( rfac*knlatsq(i) )/knlatsq(i)
!
      if (aij2 > TEN2m8) then
         tfac = -(knlat_x(i)*aij(1)+knlat_y(i)*aij(2)+knlat_z(i)*aij(3))
         sintfac = sin(tfac)
         costfac = cos(tfac)
      else
         tfac = ZERO
         sintfac = ZERO
         costfac = ONE
      endif
!
!ywg  cfac = vlen                  ! 12/3/14
      cfac = -vlen
      do l = 1, lmax_mad, 2
         do m = -l, l
            kl = (l+1)*(l+1)-l+m
            ctmp(kl)=ctmp(kl)+expfac*cfac*sintfac*Ylm(kl)
         enddo
         cfac = -cfac*knlatsq(i)
      enddo
!
      cfac = -knlatsq(i)
      do l = 2, lmax_mad, 2
         do m = -l, l
            kl = (l+1)*(l+1)-l+m
            ctmp(kl)=ctmp(kl)+expfac*cfac*costfac*Ylm(kl)
         enddo
         cfac = -cfac*knlatsq(i)
      enddo
!
   enddo
   rfac = PI4*PI4/omegbra
   dlm(2:kmax_mad) = dlm(2:kmax_mad)+rfac*ctmp(2:kmax_mad)
!
   nullify(gamma_l)
!
   end subroutine dlsum
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calGammaFunc(a,lmax,gamma_l)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: a
!
   real (kind=RealKind) :: rfac
   real (kind=RealKind) :: erfc
!
   real (kind=RealKind), intent(out) :: gamma_l(0:lmax)
!
!  *******************************************************************
!  *  call calGammaFunc to calculate the integral:
!  *
!  *         inf       2*l
!  *  I(l) = int dx * x   * exp(-x**2)
!  *          a
!  *
!  *  and store it in gamma_l(l)
!  *
!  *  l = 0, 1, ..., lmax.
!  *
!  *  using: (2*l+1)*I(l) = 2*I(l+1) - c(l) ,
!  *
!  *
!  *         c(l) = a**(2*l+1) * exp(-a**2)
!  *
!  *         I(0) = sqrt(pi)/2 * erfc(a)
!  *
!  *  erfc(z) = 1 - erf(z)
!  *******************************************************************
!
   gamma_l(0) = sqrt(PI)*HALF*erfc(a)
!
   rfac = HALF*exp(-a*a)
   do l=1,lmax
      gamma_l(l)= HALF*(2*l-1)*gamma_l(l-1)+a**(2*l-1)*rfac
   end do
!
   end subroutine calGammaFunc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getJelliumPotential(id,rvec) result(vr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: gid
!
   integer (kind=IntKind) :: ibegin
   integer (kind=IntKind) :: n
!
   real (kind=RealKind), intent(in) :: rvec(3)
   real (kind=RealKind) :: vr
!
   real (kind=RealKind) :: aij(3), aijm
   real (kind=RealKind) :: r0tm
   real (kind=RealKind) :: term0
   real (kind=RealKind) :: term12
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getJelliumPotential','Invalid local atom index',id)
   endif
!
   gid = GlobalIndex(id)
!
!  *******************************************************************
!  calculate Madelung constant matrix:
!
!     for i <> j,
!                                   2   2       -> ->
!            4*pi          1    -eta *Kq /4 - i*Kq*aij
!     M   =  ---- * sum  ----- e
!      ij    tau    q<>0    2
!                         Kq
!
!                          ->   ->                  2
!                 1 - erf(|Rn + aij|/eta)     pi*eta
!          + sum ------------------------- - ---------
!             n         ->   ->                 tau
!                      |Rn + aij|
!
!     for i = j,
!                                   2   2
!            4*pi          1    -eta *Kq /4
!     M   =  ---- * sum  ----- e
!      ii    tau    q<>0    2
!                         Kq
!
!                                             2
!                  1 - erf(Rn/eta)      pi*eta           2
!          + sum  ----------------- - ---------- - --------------
!            n<>0         Rn             tau        sqrt(pi)*eta
!
!     eta is the Ewald parameter;
!     tau, atom_posi_*, rslat_* and knlat_* are in the units of a0=1;
!     madmat is in the units of a0=alat;
!  *******************************************************************
!
   term0=-PI*eta*eta/omegbra
!
!  ===================================================================
!  start the potential calculation
!  ===================================================================
   vr = ZERO
   do n =1,GlobalNumAtoms
!     ================================================================
!     aij is in the units of a0 = 1................................
!     ================================================================
      aij(1) = rvec(1)/a0 + atom_position(1,gid) - atom_position(1,n)
      aij(2) = rvec(2)/a0 + atom_position(2,gid) - atom_position(2,n)
      aij(3) = rvec(3)/a0 + atom_position(3,gid) - atom_position(3,n)
      aijm = sqrt(aij(1)*aij(1) + aij(2)*aij(2) + aij(3)*aij(3))
!
      if( aijm <= TEN2m8 ) then
         ibegin=2
         r0tm=-TWO/sqrt(PI)/eta
      else
         ibegin=1
         r0tm=ZERO
      endif
!
!     ================================================================
!     perform the reciprocal-space sum and real-space sum
!     ----------------------------------------------------------------
      call madsum(ibegin,aij,term12)
!     ----------------------------------------------------------------
      vr = vr + term12+r0tm+term0
   enddo                                  ! end do loop over n
   vr = vr/a0                             ! reset the unit
!
   end function getJelliumPotential
!  ===================================================================
end module MadelungModule
