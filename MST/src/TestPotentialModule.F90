module TestPotentialModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicTypeDefinitionsModule, only : GridStruct 
   use RadialGridModule, only : getGrid
   use MathParamModule, only : ZERO, ONE, HALF, SQRT_PI, PI4, TEN2m8
   use MathParamModule, only : CZERO, CONE, SQRTm1, Y0, TWO
!
public :: initTestPotential, &
          endTestPotential,  &
          getTestPotential,  &
          readTestPotential, &
          getTestPotEf,      &
          getTestValenceNum, &
          getTestV0,         &
          testPotential
!
private
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
!
   type PotentialStruct
      integer (kind=IntKind) :: ptype
      integer (kind=IntKind) :: lmax
      integer (kind=IntKind) :: jmax
      real (kind=RealKind) :: xval
      real (kind=RealKind) :: afac(3)
      real (kind=RealKind) :: ufac(3)
      type (GridStruct), pointer :: Grid
      complex (kind=CmplxKind), pointer :: pot_l(:,:)
   end type PotentialStruct
!
   type (PotentialStruct), allocatable :: TestPoten(:,:)
!
   logical :: Initialized = .false.
!
   real (kind=RealKind) :: Ef
   real (kind=RealKind) :: u0(2)
!
   integer (kind=IntKind), parameter :: MuffinTinTestPotential = 0
   integer (kind=IntKind), parameter :: EmptyLatticePotential = 1
   integer (kind=IntKind), parameter :: MathieuPotential  = 2
   integer (kind=IntKind), parameter :: CoulombPotential  = 3
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initTestPotential(na,ns)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, ns
   integer (kind=IntKind) :: id, is
!
   if (na < 1) then
      call ErrorHandler('initTestPotential','Invalid number of atoms',na)
   else if (ns < 1 .or. ns > 2) then
      call ErrorHandler('initTestPotential','Invalid number of spins',ns)
   endif
!
   LocalNumAtoms = na
   n_spin_pola = ns
!
   allocate( TestPoten(na,ns) )
!
   do is = 1, ns
      do id = 1, na
         TestPoten(id,is)%jmax = 0
         nullify( TestPoten(id,is)%pot_l )
      enddo
   enddo
!
   Ef = -10.0d0
!
   Initialized = .true.
!
   end subroutine initTestPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endTestPotential()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id, is
!
   do is = 1, n_spin_pola
      do id = 1, LocalNumAtoms
         if ( TestPoten(id,is)%jmax >= 1) then
            deallocate(TestPoten(id,is)%pot_l)
         endif
         nullify(TestPoten(id,is)%Grid)
      enddo
   enddo
!
   deallocate( TestPoten )
!
   Initialized = .false.
!
   end subroutine endTestPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readTestPotential(lmax)
!  ===================================================================
   use PotentialTypeModule, only : isMuffinTinTestPotential, &
                                   isEmptyLatticePotential,  &
                                   isMathieuPotential, isCoulombPotential
!
   use AtomModule, only : getInPotFileName
!
   implicit none
!
   character (len=80) :: text
   character (len=50) :: radius_info
   character (len=1) :: radius_type
   character (len=13), parameter :: sname='readTestPotential'
!
   integer (kind=IntKind), intent(in) :: lmax(*)
   integer (kind=IntKind) :: ia, is, js, ns, status
   integer (kind=IntKind), parameter :: funit=81
!
   real (kind=RealKind) :: v0(2), u(3,2), a(3,2), e, val(2), Za
!
   do ia = 1, LocalNumAtoms
      if (lmax(ia) < 0) then
!        -------------------------------------------------------------
         call ErrorHandler('readTestPotential','lmax < 0',lmax(ia))
!        -------------------------------------------------------------
      else if (lmax(ia) > 0 .and. (isMuffinTinTestPotential() .or.    &
                            isEmptyLatticePotential())) then
!        -------------------------------------------------------------
         call ErrorHandler('readTestPotential','lmax /= 0',lmax(ia))
!        -------------------------------------------------------------
      endif
!
!     ----------------------------------------------------------------
      open(unit=funit,file=getInPotFileName(ia),form='formatted',status='old')
!
      read(funit,'(a)')text
      do while(text(1:1) == '#' .or. text(1:1) == '!' .or. text(1:1) == ' ')
         read(funit,'(a)',iostat=status)text
      enddo
      read(text,*)Za
!
      read(funit,'(a)')text
      text = adjustl(text)
      do while(text(1:1) == '#' .or. text(1:1) == '!' .or. text(1:1) == ' ')
         read(funit,'(a)',iostat=status)text
      enddo
      read(text,*)ns,e,radius_info
!
      if (ns /= n_spin_pola) then
         call ErrorHandler('readTestPotential','invalid total spin index',ns)
      endif
!
      Ef = max(Ef,e)
!
      radius_info = adjustl(radius_info)
      radius_type = radius_info(1:1)
!
      do is = 1, n_spin_pola
         read(funit,'(a)')text
         text = adjustl(text)
         do while(text(1:1) == '#' .or. text(1:1) == '!' .or. text(1:1) == ' ')
            read(funit,'(a)',iostat=status)text
         enddo
         if (isMathieuPotential()) then
            read(text,*)js,val(is),v0(is),u(1:3,is),a(1:3,is)
         else
            read(text,*)js,val(is),v0(is)
         endif
         if (js /= is) then
            call ErrorHandler('readTestPotential','invalid spin index',js)
         endif
      enddo
!
      close(unit=funit)
!     ----------------------------------------------------------------
!
      do is = 1, n_spin_pola
         if (isMuffinTinTestPotential()) then
!           ----------------------------------------------------------
            call setMuffinTinTestPotential(ia,is,v0(is))
!           ----------------------------------------------------------
         else if (isEmptyLatticePotential()) then
!           ----------------------------------------------------------
            call setEmptyLatticePotential(ia,is,v0(is))
!           ----------------------------------------------------------
         else if (isMathieuPotential()) then
!           ----------------------------------------------------------
            call setMathieuPotential(ia,is,v0(is),u(1:3,is),a(1:3,is),lmax(ia))
!           ----------------------------------------------------------
         else if (isCoulombPotential()) then
!           ----------------------------------------------------------
            call setCoulombPotential(ia,is,v0(is),Za,radius_type)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call ErrorHandler('readTestPotential',                    &
                              'This type of test potential is not implemented')
!           ----------------------------------------------------------
         endif
         TestPoten(ia,is)%xval = val(is)
      enddo
   enddo
!
   end subroutine readTestPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setMuffinTinTestPotential(id,is,v0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is
   integer (kind=IntKind) :: ir, jmt, jend_plus_n
!
   real (kind=RealKind), intent(in) :: v0
!
   complex (kind=CmplxKind) :: c
!
   if (.not.Initialized) then
      call ErrorHandler('setMuffinTinTestPotential','Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setMuffinTinTestPotential',                   &
                        'Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('setMuffinTinTestPotential','Invalid spin index',is)
   endif
!
   TestPoten(id,is)%Grid => getGrid(id)
   TestPoten(id,is)%lmax = 0
   TestPoten(id,is)%jmax = 1
   TestPoten(id,is)%afac(1:3) = ZERO
   TestPoten(id,is)%ufac(1:3) = ZERO
   u0(is) = v0
!
   jmt = TestPoten(id,is)%Grid%jmt
   jend_plus_n = TestPoten(id,is)%Grid%jend_plus_n
   allocate( TestPoten(id,is)%pot_l(jend_plus_n,1) )
!
   c = cmplx(v0/Y0,ZERO,kind=CmplxKind)
!
   do ir = 1, jmt
      TestPoten(id,is)%pot_l(ir,1) = c
   enddo
!
   do ir = jmt+1, jend_plus_n
      TestPoten(id,is)%pot_l(ir,1) = CZERO
   enddo
!
   TestPoten(id,is)%ptype = MuffinTinTestPotential
!
   end subroutine setMuffinTinTestPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setEmptyLatticePotential(id,is,v0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is
   integer (kind=IntKind) :: ir, jend_plus_n
!
   real (kind=RealKind), intent(in) :: v0
!
   complex (kind=CmplxKind) :: c
!
   if (.not.Initialized) then
      call ErrorHandler('setEmptyLatticePotential','Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setEmptyLatticePotential',                   &
                        'Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('setEmptyLatticePotential','Invalid spin index',is)
   endif
!
   TestPoten(id,is)%Grid => getGrid(id)
   TestPoten(id,is)%lmax = 0
   TestPoten(id,is)%jmax = 1
   TestPoten(id,is)%afac(1:3) = ZERO
   TestPoten(id,is)%ufac(1:3) = ZERO
   u0(is) = v0
!
   jend_plus_n = TestPoten(id,is)%Grid%jend_plus_n
   allocate( TestPoten(id,is)%pot_l(jend_plus_n,1) )
!
   c = cmplx(v0/Y0,ZERO,kind=CmplxKind)
!
   do ir = 1, jend_plus_n
      TestPoten(id,is)%pot_l(ir,1) = c
   enddo
!
   TestPoten(id,is)%ptype = EmptyLatticePotential
!
   end subroutine setEmptyLatticePotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setCoulombPotential(id,is,v0,Za,rt)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is
   integer (kind=IntKind) :: ir, jend_plus_n, n
!
   character (len=1), intent(in) :: rt
!
   real (kind=RealKind), intent(in) :: v0, Za
!
   complex (kind=CmplxKind) :: c0, c1
!
   if (.not.Initialized) then
      call ErrorHandler('setEmptyLatticePotential','Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setEmptyLatticePotential',                   &
                        'Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('setEmptyLatticePotential','Invalid spin index',is)
   endif
!
   TestPoten(id,is)%Grid => getGrid(id)
   TestPoten(id,is)%lmax = 0
   TestPoten(id,is)%jmax = 1
   TestPoten(id,is)%afac(1:3) = ZERO
   TestPoten(id,is)%ufac(1:3) = ZERO
   u0(is) = v0
!
   jend_plus_n = TestPoten(id,is)%Grid%jend_plus_n
   allocate( TestPoten(id,is)%pot_l(jend_plus_n,1) )
!
   if (rt == 'M' .or. rt == 'm') then
      n = TestPoten(id,is)%Grid%jmt
      write(6,'(a,i5)')'The Coulomb potential is established up to the Mufifn-tin Sphere Radius: ',n
   else
      n = TestPoten(id,is)%Grid%jend
      write(6,'(a,i5)')'The Coulomb potential is established up to the Bounding Sphere Radius: ',n
   endif
!
!
   c0 = cmplx(v0/Y0,ZERO,kind=CmplxKind)
   c1 = cmplx(TWO*Za/Y0,ZERO,kind=CmplxKind)
   TestPoten(id,is)%pot_l = CZERO
!
   do ir = 1, n
      TestPoten(id,is)%pot_l(ir,1) = c0-c1/TestPoten(id,is)%Grid%r_mesh(ir)
   enddo
!
   TestPoten(id,is)%ptype = CoulombPotential
!
   end subroutine setCoulombPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setMathieuPotential(id,is,U0c,U,A,lmax)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors, endIntegerFactors, &
                                    lofj, kofj
   use SphericalHarmonicsModule, only : calYlmConjg
!
   use BesselModule, only : SphericalBessel
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, lmax
   integer (kind=IntKind) :: ir, iend, jmax, kmax, i, jl, l, m
!
   real (kind=RealKind), intent(in) :: U0c, U(3), A(3)
   real (kind=RealKind) :: unit_vec(3,3)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: cfac
   complex (kind=CmplxKind), allocatable :: itol(:), bjl(:)
   complex (kind=CmplxKind), allocatable :: ylmCC(:), ylm_fac(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('setMathieuPotential','Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setMathieuPotential','Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('setMathieuPotential','Invalid spin index',is)
   else if (lmax < 0) then
      call ErrorHandler('setMathieuPotential','Invalid lmax',lmax)
   endif
!
!  ===================================================================
!  Mathieu-potential is given below:
!
!       ->                              
!     v(r ) = U + U *cos(A *x) + U *cos(A *y) + U *cos(A *z)
!              0   1      1       2      2       3      3
!
!                           l         l                      * ->
!           = U + 2*pi*sum i * [1+(-1) ] * [ U * j (A *r) * Y (x )
!              0        L                     1   l  1       L
!                                                            * ->
!                                           +U * j (A *r) * Y (y )
!                                             2   l  2       L
!                                                            * ->         ->
!                                           +U * j (A *r) * Y (z ) ] * Y (r )
!                                             3   l  3       L          L
!
!  ===================================================================
!
   jmax = (lmax+1)*(lmax+2)/2
   kmax = (lmax+1)**2
   TestPoten(id,is)%Grid => getGrid(id)
   TestPoten(id,is)%lmax = lmax
   TestPoten(id,is)%jmax = jmax
   TestPoten(id,is)%afac(1:3) = A(1:3)
   TestPoten(id,is)%ufac(1:3) = U(1:3)
   u0(is) = U0c
!
   iend = TestPoten(id,is)%Grid%jend_plus_n
   allocate( TestPoten(id,is)%pot_l(iend,jmax) )
!
   allocate(ylmCC(kmax),ylm_fac(jmax,3),itol(0:lmax),bjl(0:lmax))
!  ------------------------------------------------------------------
   call initIntegerFactors(lmax)
!  ------------------------------------------------------------------
!
   itol(0)=CONE
   do l=1,lmax
      itol(l)=itol(l-1)*SQRTm1
   enddo
!
   unit_vec(1,1)=ONE
   unit_vec(2,1)=ZERO
   unit_vec(3,1)=ZERO
!
   unit_vec(1,2)=ZERO
   unit_vec(2,2)=ONE
   unit_vec(3,2)=ZERO
!
   unit_vec(1,3)=ZERO
   unit_vec(2,3)=ZERO
   unit_vec(3,3)=ONE
!
!  ===================================================================
!  Calculate ylm_fac:
!
!                l               l     * ->
!         PI2 * i * U * [1 + (-1) ] * Y (x )
!                    1                 L
!
!                l               l     * ->
!         PI2 * i * U * [1 + (-1) ] * Y (y )
!                    2                 L
!
!                l               l     * ->
!         PI2 * i * U * [1 + (-1) ] * Y (z )
!                    3                 L
!
!  for 0 <= m <= l, l <= lmax
!  ===================================================================
   do i=1,3
!     ----------------------------------------------------------------
      call calYlmConjg(unit_vec(1:3,i),lmax,ylmCC)
!     ----------------------------------------------------------------
      do jl=1,jmax 
         l=lofj(jl)
         ylm_fac(jl,i)=itol(l)*(1-mod(l,2))*ylmCC(kofj(jl))*PI4*U(i)
      enddo
   enddo
!
!  ===================================================================
!  Calculate pot_l(ir) (l>=0, m>=0)
!  ===================================================================
   do ir=1,iend
      TestPoten(id,is)%pot_l(ir,1) = cmplx(U0c/Y0,ZERO,kind=CmplxKind)
   enddo
!
   do jl=2,jmax
      do ir=1,iend
         TestPoten(id,is)%pot_l(ir,jl) = CZERO
      enddo
   enddo
!
   r_mesh => TestPoten(id,is)%Grid%r_mesh(:)
   do i=1,3
      do ir=1,iend
         cfac=cmplx(r_mesh(ir)*A(i),ZERO,kind=CmplxKind)
!        =============================================================
!        calculate jl for l=0->lmax,
!        -------------------------------------------------------------
         call SphericalBessel(lmax,cfac,bjl)
!        -------------------------------------------------------------
         do jl=1,jmax
            TestPoten(id,is)%pot_l(ir,jl) = TestPoten(id,is)%pot_l(ir,jl) &
                                           +ylm_fac(jl,i)*bjl(lofj(jl))
         enddo
      enddo
   enddo
!
!  ===================================================================
   deallocate( ylmCC, ylm_fac, itol, bjl )
!  ===================================================================
!  ------------------------------------------------------------------
   call endIntegerFactors()
!  ------------------------------------------------------------------
!
   TestPoten(id,is)%ptype = MathieuPotential
!
   end subroutine setMathieuPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTestPotential(id,is,jl) result(pot_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, jl
   integer (kind=IntKind) :: iend
!
   complex (kind=CmplxKind), pointer :: pot_l(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getTestPotential','Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTestPotential','Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getTestPotential','Invalid spin index',is)
   else if (jl < 1 .or. jl > TestPoten(id,is)%jmax) then
      call ErrorHandler('getTestPotential','Invalid jl index',jl)
   endif
!
   iend = TestPoten(id,is)%Grid%jend_plus_n
   pot_l => TestPoten(id,is)%pot_l(1:iend,jl)
!
   end function getTestPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine testPotential(id,is,x,y,z)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is
   integer (kind=IntKind) :: lmax, jmax, kmax, jl, kl, ir, m, l, iend, jmt
!
   real (kind=RealKind), intent(in) :: x, y, z
   real (kind=RealKind) :: vec(3), vr, vi, ve, x_0, y_0, z_0, r
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), allocatable :: ylm(:)
   complex (kind=CmplxKind) :: vc
!
   if (.not.Initialized) then
      call ErrorHandler('testPotential','Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('testPotential','Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('testPotential','Invalid spin index',is)
   endif
!
   lmax = TestPoten(id,is)%lmax
   jmax = TestPoten(id,is)%jmax
   kmax = (lmax+1)**2
!
   allocate(ylm(kmax))
!
   r = sqrt(x*x + y*y + z*z)
   if (r < TEN2m8) then
      vec(1)=ZERO
      vec(2)=ZERO
      vec(3)=ZERO
      ylm(1) = Y0
      do kl = 2,kmax
         ylm(kl) = CZERO
      enddo
      iend = 1
   else
      vec(1)=x/r
      vec(2)=y/r
      vec(3)=z/r
!     ----------------------------------------------------------------
      call calYlm(vec(1:3),lmax,ylm)
!     ----------------------------------------------------------------
      iend = TestPoten(id,is)%Grid%jend_plus_n
   endif
   jmt = TestPoten(id,is)%Grid%jmt
!
   if (TestPoten(id,is)%ptype == EmptyLatticePotential) then
      write(6,'(//,a)') &
' ====================    Empty Lattice Potential Test    ===================='
   else if (TestPoten(id,is)%ptype == MuffinTinTestPotential) then
      write(6,'(//,a)') &
' ======================    Muffin-tin Potential Test    ======================'
   else if (TestPoten(id,is)%ptype == MathieuPotential) then
      write(6,'(//,a)') &
' =======================    Mathieu Potential Test    ======================='
   endif
!
   write(6,'(/,13x,a)') &
           '(x, y, z)                   CALC. POTENTIAL    EXACT POTENTIAL'
!
   r_mesh => TestPoten(id,is)%Grid%r_mesh(:)
!
   do ir = 1, iend
      x_0 = vec(1)*TestPoten(id,is)%afac(1)*r_mesh(ir)
      y_0 = vec(2)*TestPoten(id,is)%afac(2)*r_mesh(ir)
      z_0 = vec(3)*TestPoten(id,is)%afac(3)*r_mesh(ir)
      if (TestPoten(id,is)%ptype == MuffinTinTestPotential .and. ir > jmt) then
         ve = ZERO
      else
         ve = u0(is)                                                  &
             +TestPoten(id,is)%ufac(1)*cos(x_0)                        &
             +TestPoten(id,is)%ufac(2)*cos(y_0)                        &
             +TestPoten(id,is)%ufac(3)*cos(z_0)
      endif
      vc = CZERO
      jl = 0
      kl = 0
      do l=0,lmax
         jl = jl + 1
         kl = kl + l + 1
         vc = vc + TestPoten(id,is)%pot_l(ir,jl)*ylm(kl)
         do m=1,l
            jl = jl + 1
            kl = kl + 1
            vc = vc + TestPoten(id,is)%pot_l(ir,jl)*ylm(kl)           &
                    + (1-2*mod(m,2))*ylm(kl-2*m)*                     &
                      conjg(TestPoten(id,is)%pot_l(ir,jl))
         enddo
      enddo
      vr = real(vc,kind=RealKind)
      vi = aimag(vc)
      write(6,'(1x,a1,f10.5,2x,f10.5,2x,f10.5,a1,2(4x,d15.8))')  &
           '(',vec(1)*r_mesh(ir),vec(2)*r_mesh(ir),vec(3)*r_mesh(ir),')',vr,ve
      if (abs(vi) > TEN2m8) then
         call WarningHandler('testPotential','Potential is complex',vc)
      endif
   enddo
!
   deallocate(ylm)
!
   end subroutine testPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTestPotEf() result(e)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: e
!
   if (.not.Initialized) then
      call ErrorHandler('getTestPotEf','Module is not initialized')
   endif
!
   e = Ef
!
   end function getTestPotEf
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTestValenceNum(id,is) result(v)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is
!
   real (kind=RealKind) :: v
!
   if (.not.Initialized) then
      call ErrorHandler('getTestValenceNum','Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTestValenceNum','Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getTestValenceNum','Invalid spin index',is)
   endif
!
   v = TestPoten(id,is)%xval
!
   end function getTestValenceNum
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTestV0(is) result(v)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind) :: v
!
   if (.not.Initialized) then
      call ErrorHandler('getTestV0','Module is not initialized')
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getTestV0','Invalid spin index',is)
   endif
!
   v = u0(is)
!
   end function getTestV0
!  ===================================================================
end module TestPotentialModule
