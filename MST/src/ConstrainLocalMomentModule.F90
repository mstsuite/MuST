module ConstrainLocalMomentModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, &
                                  WarningHandler
   use MathParamModule, only : ZERO, HALF, ONE, TWO, PI, PI4,     &
                               Ten2m8, Ten2m6, Ten2m2
!
public :: initConstrainLM,   &
          endConstrainLM,    &
          calConstrainLM,    &
          updateConstrainLM, &
          getFieldRms,       &
          printConstrainMoment
!
private
   logical :: Initialized = .false.
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: tstep
   integer (kind=IntKind) :: itstep_save
   integer (kind=IntKind) :: max_Nrpts
   integer (kind=IntKind), allocatable :: print_level(:)
   integer (kind=IntKind), allocatable :: Nrpts(:)
   integer (kind=IntKind), allocatable :: NumSpecies(:)
!
   real (kind=RealKind), parameter :: tol=Ten2m6
   real (kind=RealKind), target, allocatable :: PotSph_r(:,:,:,:)
   real (kind=RealKind), target :: field_rms(1:2)
!
   type ConstrainingStruct
      real (kind=RealKind) :: bcon_rms
      real (kind=RealKind) :: evec_rms
      real (kind=RealKind) :: edote
      real (kind=RealKind) :: b_con_mag
      real (kind=RealKind) :: b_con_g(3)
      real (kind=RealKind) :: b_basis(3,3)
   end type ConstrainingStruct
!
   type (ConstrainingStruct), allocatable :: SiteConstrains(:)
!
contains
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initConstrainLM( na, cant, pola, t_stp, iprint )
!  ===================================================================
   use PublicTypeDefinitionsModule, only : GridStruct
   use RadialGridModule, only : getGrid
   use AtomModule, only : getLocalEvecOld, getLocalNumSpecies
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, cant, pola, iprint(na)
!
   real (kind=RealKind), intent(in) :: t_stp
!
   integer (kind=IntKind) :: n, i, ia, MaxSpecies
!
   real (kind=RealKind) :: evec(3)
!
   type(GridStruct), pointer :: Grid
!
   LocalNumAtoms = na
   n_spin_pola = pola
   n_spin_cant = cant
   tstep       = t_stp
!
   allocate( print_level(LocalNumAtoms) )
   allocate( SiteConstrains(LocalNumAtoms) )
   allocate( Nrpts(LocalNumAtoms), NumSpecies(LocalNumAtoms) )
!
   max_Nrpts = 0; MaxSpecies = 0
   do n = 1,LocalNumAtoms
      Grid => getGrid(n)
      Nrpts(n)  = Grid%jend
!      Nrpts(n)  = Grid%jmt
      NumSpecies(n) = getLocalNumSpecies(n)
      max_Nrpts = max(max_Nrpts,Nrpts(n))
      MaxSpecies = max(MaxSpecies, NumSpecies(n))
      print_level(n) = iprint(n)
   enddo
!
   if ( n_spin_cant == 2 ) then
      allocate(PotSph_r(max_Nrpts,n_spin_cant,MaxSpecies,LocalNumAtoms))
      do n = 1,LocalNumAtoms
         SiteConstrains(n)%evec_rms = ZERO
         SiteConstrains(n)%bcon_rms = ZERO
         SiteConstrains(n)%edote = ZERO
         SiteConstrains(n)%b_con_mag = ZERO
         SiteConstrains(n)%b_con_g(1:3) = ZERO
         do ia = 1, NumSpecies(n)
            evec = getLocalEvecOld(n)
            call set_b_basis(n,ia,evec)
         enddo
      enddo
   else
      do n = 1,LocalNumAtoms
         SiteConstrains(n)%evec_rms = ZERO
         SiteConstrains(n)%bcon_rms = ZERO
         SiteConstrains(n)%edote = ZERO
         SiteConstrains(n)%b_con_mag = ZERO
         SiteConstrains(n)%b_basis(1:3,1:3) = ZERO
         SiteConstrains(n)%b_con_g(1:3) = ZERO
         do i = 1,3
            SiteConstrains(n)%b_basis(i,i) = ONE
         enddo
      enddo
   endif
!
   Initialized = .true.
   nullify(Grid)
!
   end subroutine initConstrainLM
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endConstrainLM()
!  ===================================================================
   implicit none
!
   deallocate( print_level )
   deallocate( SiteConstrains )
   deallocate( Nrpts, NumSpecies )
   if ( n_spin_cant == 2 ) then
      deallocate( PotSph_r )
   endif
!
   Initialized = .false.
!
   end subroutine endConstrainLM
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printConstrainMoment(id)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in), optional :: id
   integer (kind=IntKind) :: id_fst, id_lst, i
!
   write(6,'(/,80(''-''))')
   write(6,'(/,23x,a)')'************************************'
   write(6,'( 23x,a )')'* Output from printConstrainMoment *'
   write(6,'(23x,a,/)')'************************************'
!
   if (present(id)) then
      id_fst=id
      id_lst=id
   else
      id_fst=1
      id_lst=LocalNumAtoms
   endif
!
   do i = id_fst,id_lst
      write(6,'(a,i4)') "Site:: ",i
      write(6,'(3x,a,1x,f15.8)') "evec_rms :: ", SiteConstrains(i)%evec_rms
      write(6,'(3x,a,1x,f15.8)') "bcon_rms :: ", SiteConstrains(i)%bcon_rms
      write(6,'(3x,a,1x,f15.8)') "edote    :: ", SiteConstrains(i)%edote
      write(6,'(3x,a,1x,f15.8)') "B_con_mag:: ", SiteConstrains(i)%b_con_mag
      write(6,'(3x,a,3(1x,f15.8))') "B_con_g  :: ", SiteConstrains(i)%b_con_g(1:3)
      write(6,'(3x,a)') "Basis :: "
      write(6,'(3x,3(1x,f15.8))') SiteConstrains(i)%b_basis(1:3,1)
      write(6,'(3x,3(1x,f15.8))') SiteConstrains(i)%b_basis(1:3,2)
      write(6,'(3x,3(1x,f15.8))') SiteConstrains(i)%b_basis(1:3,3)
   enddo
!
   write(6,'(/,80(''-''))')
   end subroutine printConstrainMoment
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calConstrainLM(itstep)
!  ===================================================================
   use SpinRotationModule, only : resetSpinRotation
   use PotentialModule, only : getSphPotr, setSphPotr
   use AtomModule, only : getLocalEvecOld, setLocalEvecOld,               &
                          getLocalConstrainField, setLocalConstrainField, &
                          getMixingParam4Evec
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: itstep
!
   integer (kind=IntKind) :: j, i, n, ia
!
   real (kind=RealKind) :: evec_rms
   real (kind=RealKind) :: evec(3)
   real (kind=RealKind) :: evec_r(3)
!
   real (kind=RealKind) :: c
   real (kind=RealKind) :: b_mag
   real (kind=RealKind) :: b_con_mag
   real (kind=RealKind) :: b_con(3)
   real (kind=RealKind) :: b_con_g(3)
   real (kind=RealKind) :: b_basis(3,3)
!
   real (kind=RealKind), pointer :: pot_up(:), pot_down(:)
   real (kind=RealKind), allocatable :: vvr(:), br(:)
!
   allocate( vvr(max_Nrpts), br(max_Nrpts) )
!  ===================================================================
!  Check to see if B_con=0: if so then set evec_r=evec and return..
!  ===================================================================
!
   NumSitesLoop: do n = 1,LocalNumAtoms
!
      Loop_ia: do ia = 1, NumSpecies(n)
!
         evec = getLocalEvecOld(n)
         if ( itstep > 1 .and. itstep/=itstep_save ) then
            call resetSpinRotation( n, evec )
            b_con = getLocalConstrainField(n)
            b_basis = SiteConstrains(n)%b_basis
            b_con_mag = ZERO
            do i = 1,3
               b_con_g(i) = b_con(1)*b_basis(i,1) +                   &
                            b_con(2)*b_basis(i,2) +                   &
                            b_con(3)*b_basis(i,3)
               b_con_mag  = b_con_mag + b_con_g(i)*b_con_g(i)
               SiteConstrains(n)%b_con_g(i) = b_con_g(i)
            enddo
            b_con_mag = sqrt(b_con_mag)
            SiteConstrains(n)%b_con_mag  = b_con_mag
!
!           ----------------------------------------------------------
            call calNewMomentDirection( n, ia, evec )
!           ----------------------------------------------------------
            call set_b_basis( n, ia, evec )
!           ----------------------------------------------------------
            call setLocalEvecOld( n, evec )
!           ----------------------------------------------------------
            if ( n==LocalNumAtoms .and. ia==NumSpecies(n)) then
               itstep_save = itstep
            endif
         endif
!
         if ( n_spin_cant /= 2 ) then
            call resetSpinRotation( n, evec )
            cycle Loop_ia
         endif
!
         b_con = getLocalConstrainField(n)
         b_basis = SiteConstrains(n)%b_basis
         b_con_mag = ZERO
         do i = 1,3
            b_con_g(i) = b_con(1)*b_basis(i,1) + b_con(2)*b_basis(i,2) +  &
                         b_con(3)*b_basis(i,3)
            b_con_mag  = b_con_mag + b_con_g(i)*b_con_g(i)
            SiteConstrains(n)%b_con_g(i) = b_con_g(i)
         enddo
         b_con_mag = sqrt(b_con_mag)
         SiteConstrains(n)%b_con_mag  = b_con_mag
!
         if ( abs(b_con_mag) < tol ) then
            call resetSpinRotation( n, evec )
         else if ( abs(b_con_mag) > ONE ) then
            call WarningHandler('calConstrainLM','b_con_mag is .gt. one')
         endif
         pot_up   => getSphPotr(n,ia,1)
         pot_down => getSphPotr(n,ia,2)
!        =============================================================
!        save the spherical component of the potential before adding
!        constraint field. The potential will be restored when the
!        constraint field is updated(see updateConstrainLM )
!        =============================================================
         PotSph_r(1:Nrpts(n),1,ia,n) = pot_up(1:Nrpts(n))
         PotSph_r(1:Nrpts(n),2,ia,n) = pot_down(1:Nrpts(n))
!        =============================================================
!        set up V and B ..............................................
!        =============================================================
!        -------------------------------------------------------------
         call v_plus_minus(half,pot_down,pot_up,vvr,br,Nrpts(n))
!        -------------------------------------------------------------
!        =============================================================
!        calculate the new direction .................................
!        =============================================================
         b_mag = ZERO
         do i = 1,3
            evec_r(i) = evec(i) + b_con_g(i)
            b_mag     = b_mag + evec_r(i)*evec_r(i)
         end do
         b_mag = sqrt(b_mag)
         do i = 1,3
            evec_r(i) = evec_r(i)/b_mag
         end do
!
         if ( print_level(n) > 1 ) then
            write(6,'('' CONSTRAINT: vectors:'')')
            write(6,'(29x,'' Direction:'',18x,''Magnitude :'')')
            write(6,'('' B xcor      :'',3f12.5,5x)') (evec(i),i=1,3)
            if ( abs(b_con_mag) >= tol ) then
               write(6,'('' B constraint:'',3f12.5,5x,f10.5)')            &
                                  (b_con_g(i)/b_con_mag,i=1,3),b_con_mag
            else
               write(6,'('' B constraint:'',3f12.5,5x,f10.5)')            &
                                           (b_con_g(i),i=1,3),b_con_mag
            endif
            write(6,'('' B sum       :'',3f12.5,5x,f10.5)')               &
                                        (evec_r(i),i=1,3),b_mag
         endif
!        =============================================================
!        generate new B=sqrt(1+b_con_mag^2)*B_xc......................
!        =============================================================
         c = sqrt( ONE + b_con_mag*b_con_mag )
         do j = 1,Nrpts(n)
            br(j) = c*br(j)
         enddo
!        =============================================================
!        generate new v_up and v_down.................................
!        =============================================================
!        -------------------------------------------------------------
         call v_plus_minus( ONE, vvr, br, pot_down, pot_up, Nrpts(n) )
!        -------------------------------------------------------------
         call setSphPotr(n,ia,1,pot_up)
!        -------------------------------------------------------------
         call setSphPotr(n,ia,2,pot_down)
!        -------------------------------------------------------------
         call resetSpinRotation( n, evec_r )
!        -------------------------------------------------------------
         if ( print_level(n) > 1 ) then
            write(6,'(/,a,3(1x,f12.8))') "CONSTRAINT:: NewB = ",c,    &
                                   pot_down( Nrpts(n)), pot_up( Nrpts(n))
         endif
!
      enddo Loop_ia
   enddo NumSitesLoop
!
   nullify(pot_up, pot_down)
   deallocate(vvr,br)
!
   end subroutine calConstrainLM
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calNewMomentDirection( id, ia, mdir )
!  ===================================================================
   use RadialGridModule, only : getRmesh
   use PotentialModule, only : getSphPotr
   use DataServiceCenterModule, only : getDataStorage, RealMark
   use PotentialTypeModule, only : isFullPotential
   use AtomModule, only : getLocalNumSpecies
   use IntegrationModule, only : calIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), intent(inout) :: mdir(3)
!
   integer (kind=IntKind) ::   i
!
   real (kind=RealKind) :: rtpi4, moment, emach
   real (kind=RealKind) :: beff, dummy
   real (kind=RealKind) :: mnorm, rms
!
   real (kind=RealKind) :: m(3), b(3), mb(3), mmb(3)
   real (kind=RealKind), pointer :: pmom(:,:)
   real (kind=RealKind), pointer :: pot_u(:)
   real (kind=RealKind), pointer :: pot_d(:)
   real (kind=RealKind), pointer :: rad(:)
   real (kind=RealKind), allocatable :: rtmp(:)
   real (kind=RealKind), allocatable :: dr(:)
   real (kind=RealKind), allocatable :: dvr(:)
   real (kind=RealKind), allocatable :: w1(:)
!
   emach = Ten2m2
!
   rms    = ZERO
   rtpi4  = sqrt(PI4)
!
   allocate( rtmp(0:Nrpts(id)) )
   allocate( dr(0:Nrpts(id)) )
   allocate( dvr(0:Nrpts(id)) )
   allocate( w1(0:Nrpts(id)) )
!  ===================================================================
!  we assume that rho=rho(r)*r^2 and pot=v(r)*r
!  ===================================================================
   if ( isFullPotential() ) then
      call WarningHandler('calNewMomentDirection','Needs fix')
   endif
   rad => getRmesh(id)
   pot_u => getSphPotr(id,ia,1)
   pot_d => getSphPotr(id,ia,2)
   pmom  => getDataStorage(id,'OldSphericalMomentDensity', &
                           Nrpts(id),getLocalNumSpecies(id),RealMark)
!
!  ===================================================================
!  In full potential case, the interation should be carried over the
!  entire cell.
!  ===================================================================
   rtmp(0) = ZERO
   dr(0)  = ZERO
   dvr(0) = ZERO
   do i = 1,Nrpts(id)
      rtmp(i) = sqrt(rad(i))
      dr(i)  = pmom(i,ia)*rad(i)*rad(i)
!     dvr(i) = rtpi4*rad(i)*pmom(i,ia)*(pot_u(i) - pot_d(i))
      dvr(i) = rad(i)*pmom(i,ia)*(pot_u(i) - pot_d(i))
   enddo
!
!  ------------------------------------------------------------------
   call calIntegration( Nrpts(id)+1, rtmp(0:Nrpts(id)), dr(0:Nrpts(id)), &
                        w1(0:Nrpts(id)), 1 )
!  ------------------------------------------------------------------
   moment = w1(Nrpts(id))*PI4*TWO
!
!  ------------------------------------------------------------------
   call calIntegration( Nrpts(id)+1, rtmp(0:Nrpts(id)), dvr(0:Nrpts(id)),&
                        w1(0:Nrpts(id)), 1 )
!  ------------------------------------------------------------------
   beff = w1(Nrpts(id))*PI4*TWO
!
   if ( abs(moment) > emach ) then
      beff = beff / moment
   else
      beff = ZERO
   endif
!
   if ( print_level(id) > 1 ) then
      write(6,'('' LLeq: moment:'',f15.8)') moment
      write(6,'('' LLeq:   beff:'',f15.8)')   beff
   endif
!
   do i = 1,3
      b(i) = beff * SiteConstrains(id)%b_con_g(i)
      m(i) = mdir(i)
   enddo
!
   if ( print_level(id) > 1) then
      write(6,'('' LLeq: cfield:'',3f15.8)') SiteConstrains(id)%b_con_g(1:3)
      write(6,'('' LLeq:      b:'',3f15.8)')      b
      write(6,'('' LLeq:      m:'',3f15.8)')      m
   endif
!  -------------------------------------------------------------------
   call vcross(m,b,mb)
   call vcross(m,mb,mmb)
!  -------------------------------------------------------------------
   if ( print_level(id) > 1) then
      write(6,'('' LLeq:     mb:'',3f15.8)')   mb
      write(6,'('' LLeq:    mmb:'',3f15.8)')  mmb
   endif
   dummy = ZERO
   mnorm = ZERO
   do i = 1,3
      m(i) = m(i) - tstep * mmb(i)
      mnorm = mnorm + m(i)**2
      dummy = dummy + (m(i)-mdir(i))**2
   enddo
   if (mnorm > emach) then
      mnorm = ONE/sqrt(mnorm)
      mdir(1:3) = m(1:3) * mnorm
   else
      m(1:3) = ZERO
      mdir(1:2) = ZERO; mdir(3) = ONE
   endif
   dummy = sqrt(dummy)
!
   rms = rms + dummy
   SiteConstrains(id)%bcon_rms = rms
!
!  write(6,*) ' RMS_MOM=',rms
   nullify( pmom, pot_u, pot_d, rad)
   deallocate( rtmp, dr, dvr, w1 )
!
   end subroutine calNewMomentDirection
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_b_basis( id, ia, evec )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), intent(in) :: evec(3)
!
   integer (kind=IntKind) :: i, j, k
!
   real (kind=RealKind) :: b_basis(3,3)
   real (kind=RealKind) :: sint, cost
   real (kind=RealKind) :: sinp, cosp
   real (kind=RealKind) :: dot, norm
!
!  ===================================================================
!  determine axis system in which b_con is specified..................
!  ===================================================================
   norm = evec(1)*evec(1) + evec(2)*evec(2) + evec(3)*evec(3)
   if ( abs(norm-one) > tol ) then
      call WarningHandler("set_b_basis","Evec is not a unit vector" )
   endif
   norm = sqrt(norm)
!
   cost = evec(3)/norm
   sint = sqrt( ONE - cost*cost )
   if ( abs(sint) > tol ) then
      cosp = evec(1)/(sint*norm)
      sinp = evec(2)/(sint*norm)
   else
      cosp = ONE
      sinp = ZERO
      cost = ONE
      sint = ZERO
   endif
!
   b_basis(1,3) = evec(1)/norm
   b_basis(2,3) = evec(2)/norm
   b_basis(3,3) = evec(3)/norm
!
   b_basis(1,1) = cost*cosp
   b_basis(2,1) = cost*sinp
   b_basis(3,1) = -sint
!
   b_basis(1,2) = b_basis(2,3)*b_basis(3,1) - b_basis(3,3)*b_basis(2,1)
   b_basis(2,2) =-b_basis(1,3)*b_basis(3,1) + b_basis(3,3)*b_basis(1,1)
   b_basis(3,2) = b_basis(1,3)*b_basis(2,1) - b_basis(2,3)*b_basis(1,1)
!
!  ===================================================================
   if ( print_level(id) >= 0 ) then
      write(6,'(10x,''SET_B_BASIS: my_b_basis 1'',t40,''='',3f10.5)') &
               (b_basis(i,1),i=1,3)
      write(6,'(10x,''             my_b_basis 2'',t40,''='',3f10.5)') &
               (b_basis(i,2),i=1,3)
      write(6,'(10x,''             my_b_basis 3'',t40,''='',3f10.5)') &
               (b_basis(i,3),i=1,3)
   endif
!
!  ===================================================================
!  Check that B-axes are mutually othogonal...........................
!  ===================================================================
   do i = 1,3
      do j = i+1,3
         dot = ZERO
         do k = 1,3
            dot = dot + b_basis(k,i)*b_basis(k,j)
         enddo
         if ( abs(dot) > tol ) then
            write(6,'(''b_basis incorrect i,j,abs(dot)'',2i5,d16.8)')  &
                                                       i, j, abs(dot)
            call ErrorHandler("set_b_basis",'Error b_basis incorrect')
         endif
      enddo
   enddo
!
!  ===================================================================
   if ( print_level(id) >= 0 ) then
      write(6,'(10x,''SET_B_BASIS: my_b_basis 1'',t40,''='',3f10.5)') &
               (b_basis(i,1),i=1,3)
      write(6,'(10x,''             my_b_basis 2'',t40,''='',3f10.5)') &
               (b_basis(i,2),i=1,3)
      write(6,'(10x,''             my_b_basis 3'',t40,''='',3f10.5)') &
               (b_basis(i,3),i=1,3)
   endif
!
   SiteConstrains(id)%b_basis(1:3,1:3) = b_basis(1:3,1:3)
!
   end subroutine set_b_basis
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateConstrainLM()
!  ===================================================================
   use Atom2ProcModule, only : getGlobalIndex
!
   use AtomModule, only : getLocalConstrainField,                      &
                          setLocalConstrainField,                      &
                          getLocalEvecOut, getLocalEvecNew,            &
                          setLocalEvecOld, getLocalEvecOld,            &
                          getMixingParam4Evec
!
   use SystemModule, only : setConstrainField
   use PotentialModule, only : setSphPotr
!
   implicit   none
!
   integer (kind=IntKind) :: n,ng,j,i
!
   real (kind=RealKind) :: evec_rms
   real (kind=RealKind) :: bcon_rms
   real (kind=RealKind) :: evec(3)
   real (kind=RealKind) :: evec_new(3)
   real (kind=RealKind) :: b_basis(3,3)
   real (kind=RealKind) :: b_con_g(3)
   real (kind=RealKind) :: b_con_o(3)
   real (kind=RealKind) :: b_con(3)
   real (kind=RealKind) :: b_con_in(3)
   real (kind=RealKind) :: moment(3)
   real (kind=RealKind) :: moment_dir(3)
   real (kind=RealKind) :: moment_mag
   real (kind=RealKind) :: mproj, bproj, mbproj
   real (kind=RealKind) :: ctmp
   real (kind=RealKind) :: b_con_mag
   real (kind=RealKind) :: bb1,bb2
   real (kind=RealKind) :: alphev
   real (kind=RealKind), pointer :: pot_r(:)
!
!  ===================================================================
!  Generate B_con in global frame.....................................
!  ===================================================================
   LoopAtoms: do n = 1, LocalNumAtoms
!
      ng = getGlobalIndex(n)
      evec    = getLocalEvecOld(n)
      evec_new = getLocalEvecNew(n)
      alphev = getMixingParam4Evec(n)
      if ( abs(alphev) > tol ) then
!
         SiteConstrains(n)%edote = evec(1)*evec_new(1) +              &
                                   evec(2)*evec_new(2) +              &
                                   evec(3)*evec_new(3)
!        =============================================================
!        update the Evec for the new nstep iteration
!        also update the Constrain Field stored in System module
!        -------------------------------------------------------------
         call setLocalEvecOld(n,evec_new)
!        -------------------------------------------------------------
         evec_rms = ZERO
         do i=1,3
            evec_rms = evec_rms + (evec_new(i)-evec(i))**2
         enddo
         SiteConstrains(n)%evec_rms=sqrt(evec_rms)
         SiteConstrains(n)%bcon_rms=ZERO
         cycle LoopAtoms
      endif
!
      b_basis = SiteConstrains(n)%b_basis
      b_con   = getLocalConstrainField(n)
      moment  = getLocalEvecOut(n)
!
      b_con_mag=zero
      do i = 1,3
         b_con_g(i)  = (b_con(1)*b_basis(i,1) + b_con(2)*b_basis(i,2)) + &
                       b_con(3)*b_basis(i,3)
         b_con_o(i)  = b_con_g(i)
         b_con_in(i) = b_con(i)
         b_con_mag   = b_con_mag+b_con_g(i)*b_con_g(i)
      enddo
      b_con_mag = sqrt(b_con_mag)
!
!     ================================================================
!     Generate unit vector in the direction of moment.................
!     ================================================================
      moment_mag = sqrt(moment(1)**2 + moment(2)**2 + moment(3)**2)
      do i = 1,3
         moment_dir(i) = moment(i)/moment_mag
      enddo
!
      evec_rms = ZERO
      do i=1,3
         evec_rms = evec_rms + (moment(i)-evec(i))**2
      enddo
      SiteConstrains(n)%evec_rms=sqrt(evec_rms)
!
!     ================================================================
!     calculate moment_dir*evec.......................................
!     ================================================================
      mproj = moment_dir(1)*evec(1) + moment_dir(2)*evec(2) +          &
              moment_dir(3)*evec(3)
      bproj = moment_dir(1)*b_con_g(1) + moment_dir(2)*b_con_g(2) +    &
              moment_dir(3)*b_con_g(3)
!
      mbproj = ZERO
      do i=1,3
         mbproj = mbproj + moment_dir(i)*(evec(i)+b_con_g(i))
      enddo
!     ================================================================
!     generate new constraining field in global frame.................
!     ================================================================
      bb1 = zero
      do i = 1,3
!         b_con_g(i) = (b_con_g(i)+evec(i))-mbproj*moment_dir(i)
         b_con_g(i) = b_con_g(i)-(moment_dir(i)-mproj*evec(i))
         bb1 = bb1+b_con_g(i)**2
      enddo
      bb1 = sqrt(bb1)
!
!     ================================================================
!     project out components in b_basis frame ...................
!     ================================================================
!
      bcon_rms = ZERO
      do j = 1,3
         ctmp = ZERO
         do i = 1,3
            ctmp = ctmp + b_con_g(i)*b_basis(i,j)
         enddo
         bcon_rms = bcon_rms + (ctmp-b_con(j))**2
         b_con(j) = ctmp
         if (j==3) b_con(j) = ZERO
      enddo
      SiteConstrains(n)%bcon_rms=sqrt(bcon_rms)
      call setLocalConstrainField(n,b_con)
!
!     ====================
!     Check..........
!     ====================
      bb2 = ZERO
      do i = 1,3
         b_con_g(i) = b_con(1)*b_basis(i,1) + b_con(2)*b_basis(i,2) +  &
                      b_con(3)*b_basis(i,3)
         bb2 = bb2 + b_con_g(i)*b_con_g(i)
      enddo
      bb2 = sqrt(bb2)
      if ( abs(bb1-bb2) >= tol ) then
         write(6,*)"updateConstrains:: abs(bb1-bb2) > tol", bb2, bb2
         call WarningHandler("updateConstrains",'Wrong algorithm coded')
      endif
      SiteConstrains(n)%b_con_g = b_con_g
      SiteConstrains(n)%b_con_mag = bb2
!
!     ================================
!     Calculate errror parameter
!     ================================
      bb2 = ZERO
      do i = 1,3
         bb2 = bb2 + ( b_con(i) - b_con_in(i) )**2
      enddo
      bb2 = sqrt(bb2)
!     ================================================================
!     printout if needed ............................................
!     ================================================================
      if ( print_level(n) >= 0 ) then
         write(6,'('' NEW CONSTRAINT: vectors:'')')
         write(6,'(29x,'' Direction:'',18x,''Magnitude :'')')
         write(6,'('' evec,mproj      :'',3f12.5,5x,f10.5)')           &
                     (evec(i),i=1,3),mproj
         write(6,'('' New moment      :'',3f12.5,5x,f10.5)')           &
                     (moment(i),i=1,3),moment_mag
         if ( abs(b_con_mag) > tol ) then
            write(6,'('' Old B constraint:'',3f12.5,5x,f10.5)')        &
                                (b_con_o(i)/b_con_mag,i=1,3),b_con_mag
         else
            write(6,'('' Old B constraint:'',3f12.5,5x,f10.5)')        &
                                (b_con_o(i),i=1,3),b_con_mag
         endif
         if ( abs(bb1) > tol ) then
            write(6,'('' New B constraint:'',3f12.5,5x,f10.5)')        &
                                (b_con_g(i)/bb1,i=1,3),bb1
         else
            write(6,'('' New B constraint:'',3f12.5,5x,f10.5)')        &
                                (b_con_g(i),i=1,3),bb1
         endif
         write(6,'('' Change in b_con:'',f12.5)') bb2
      endif
!
      SiteConstrains(n)%edote = evec(1)*evec_new(1) +                  &
                                evec(2)*evec_new(2) +                  &
                                evec(3)*evec_new(3)
!     ================================================================
!     update the Evec for the new nstep iteration
!     also update the Constrain Field stored in System module
!     ----------------------------------------------------------------
      call setLocalEvecOld(n,evec_new)
      call setConstrainField(ng,b_con)
!     ----------------------------------------------------------------
      if ( n_spin_cant==2 ) then
         do i = 1,NumSpecies(n)
            pot_r => PotSph_r(1:Nrpts(n),1,i,n)
!           ----------------------------------------------------------
            call setSphPotr(n,i,1,pot_r)
!           ----------------------------------------------------------
            pot_r => PotSph_r(1:Nrpts(n),2,i,n)
!           ----------------------------------------------------------
            call setSphPotr(n,i,2,pot_r)
!           ----------------------------------------------------------
         enddo
      endif
   enddo LoopAtoms
!
   end subroutine updateConstrainLM
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine v_plus_minus( fac, vr1, vr2, vr, br, npts )
!  ===================================================================
!
   implicit   none
!
!     ****************************************************************
!     calculates either v=v_up+v_down & b=v_up-v_down (fac=1)
!     or                v_up=(v+b)/2 & v_down=(v-b)/2 (fac=1/2)
!     ****************************************************************
!
!
   integer (kind=IntKind), intent(in) :: npts
!
   real (kind=RealKind), intent(in)  :: fac
   real (kind=RealKind), target, intent(in) :: vr1(npts)
   real (kind=RealKind), target, intent(in) :: vr2(npts)
!
   real (kind=RealKind), target, intent(out) :: vr(npts)
   real (kind=RealKind), target, intent(out) :: br(npts)
!
   integer (kind=IntKind) :: n
!
   do n = 1,npts
      vr(n) = fac*( vr1(n) + vr2(n) )
      br(n) = fac*( vr1(n) - vr2(n) )
   enddo
!
   end subroutine v_plus_minus
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFieldRms(id)                             result(pf_rms)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), optional, intent(in) :: id
   integer (kind=IntKind) :: n
   real (kind=RealKind), pointer :: pf_rms(:)
!
   field_rms=ZERO
   if ( present(id) ) then
      field_rms(1) = SiteConstrains(id)%evec_rms
      field_rms(2) = SiteConstrains(id)%bcon_rms
   else
      do n=1,LocalNumAtoms
         field_rms(1) = max(field_rms(1),SiteConstrains(n)%evec_rms)
         field_rms(2) = max(field_rms(2),SiteConstrains(n)%bcon_rms)
      enddo
   endif
   pf_rms => field_rms(1:2)
!
   end function getFieldRms
!  ===================================================================
end module ConstrainLocalMomentModule
