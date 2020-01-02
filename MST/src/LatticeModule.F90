module LatticeModule
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ZERO, HALF, THIRD, ONE, TWO, PI, PI2, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler
   use MPPModule, only : MyPE
!
public :: initLattice,           &
          endLattice,            &
          genLattice,            &
          getLatticeType,        &
          getNumLatticeVectors,  &
          getLatticeVectors
!
   interface initLattice
      module procedure initLatt0, initLatt1
   end interface
!
private
!
   character (len=12) :: lattice_name = 'Null'
!
   integer (kind=IntKind), target :: num_rlat
   integer (kind=IntKind), target :: num_klat
!
   real (kind=RealKind), target :: rbrav(3,3)
   real (kind=RealKind), target :: kbrav(3,3)
   real (kind=RealKind), allocatable, target :: rlat(:,:)
   real (kind=RealKind), allocatable, target :: klat(:,:)
   real (kind=RealKind), parameter :: tol = TEN2m8
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initLatt0(brav_in)
!  ===================================================================
   use VectorModule, only : getVecLength, getDotProduct
   implicit none
!
   real (kind=RealKind), intent(in) :: brav_in(3,3)
   real (kind=RealKind) :: a0, b0, c0, cosa, cosb, cosg
!
   rbrav = brav_in
!
!  -------------------------------------------------------------------
   call genKBravais()
!  -------------------------------------------------------------------
!
   a0 = getVecLength(3,rbrav(1:3,1))
   b0 = getVecLength(3,rbrav(1:3,2))
   c0 = getVecLength(3,rbrav(1:3,3))
   cosa = getDotProduct(3,rbrav(1:3,2),rbrav(1:3,3))/(b0*c0)
   cosb = getDotProduct(3,rbrav(1:3,1),rbrav(1:3,3))/(a0*c0)
   cosg = getDotProduct(3,rbrav(1:3,1),rbrav(1:3,2))/(a0*b0)
!
   if (abs(cosa) < tol) then
      cosa = ZERO
   endif
!
   if (abs(cosb) < tol) then
      cosb = ZERO
   endif
!
   if (abs(cosg) < tol) then
      cosg = ZERO
   endif
!
!  -------------------------------------------------------------------
   call calLatticeSystem(a0,b0,c0,cosa,cosb,cosg,lattice_name)
!  -------------------------------------------------------------------
!
   end subroutine initLatt0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initLatt1(a0,b0,c0,alpha,beta,gamma,isRadian)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: a0,b0,c0,alpha,beta,gamma
   real (kind=RealKind) :: rfac, cosa, cosb, cosg, sing
!
   logical, intent(in), optional :: isRadian
!
   rfac = PI/180.0d0
!
   if (present(isRadian)) then
      if (isRadian) then
         rfac = ONE
      endif
   endif
!
   if (alpha < tol .or. beta < tol .or. gamma < tol) then
      call ErrorHandler('initLattice','Invalid lattice parameters')
   else if (a0 < tol .or. b0 < tol .or. c0 < tol) then
      call ErrorHandler('initLattice','Invalid lattice parameters')
   endif
!
   cosa = cos(alpha*rfac)
   if (abs(cosa) < tol) then
      cosa = ZERO
   endif
!
   cosb = cos(beta*rfac)
   if (abs(cosb) < tol) then
      cosb = ZERO
   endif
!
   cosg = cos(gamma*rfac)
   if (abs(cosg) < tol) then
      cosg = ZERO
      sing = ONE
   else
      sing = sqrt(ONE-cosg*cosg)
   endif
!
!  -------------------------------------------------------------------
   call calLatticeSystem(a0,b0,c0,cosa,cosb,cosg,lattice_name)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  We set up our coordinates system in such a way that x-axis is along
!  rbrav(1:3,1), and rbrav(1:3,2) is on the x-y plane.
!  ===================================================================
   rbrav(1,1) = a0
   rbrav(2,1) = ZERO
   rbrav(3,1) = ZERO
   rbrav(1,2) = b0*cosg
   rbrav(2,2) = b0*sing
   rbrav(3,2) = ZERO
   rbrav(1,3) = c0*cosb
   rbrav(2,3) = c0*(cosa-cosg*cosb)/sing
   rbrav(3,3) = c0*sqrt(ONE-cosb*cosb-(cosa-cosg*cosb)**2/sing**2)
!
!  -------------------------------------------------------------------
   call genKBravais()
!  -------------------------------------------------------------------
!
   end subroutine initLatt1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endLattice()
!  ===================================================================
   implicit none
!
   if (allocated(rlat)) then
      deallocate( rlat )
   endif
   if (allocated(klat)) then
      deallocate( klat )
   endif
!
   num_rlat = 0; num_klat = 0
   lattice_name = 'Null'
!
   end subroutine endLattice
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLatticeType() result(t)
!  ===================================================================
   implicit none
!
   character (len=12) :: t
!
   t = lattice_name
!
   end function getLatticeType
!  =====================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumLatticeVectors(s) result(n)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: s
!
   integer (kind=IntKind) :: n
!
   if (s == 'R' .or. s == 'r') then
      n = num_rlat
   else if (s == 'K' .or. s == 'k') then
      n = num_klat
   else
      call ErrorHandler('getNumLatticeVectors','Invalid space type',s)
   endif
!
   end function getNumLatticeVectors
!  =====================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLatticeVectors(s) result(v)
!  =====================================================================
   implicit none
!
   character (len=1), intent(in) :: s
!
   real (kind=RealKind), pointer :: v(:,:)
!
   if (s == 'R' .or. s == 'r') then
      v => rlat
   else if (s == 'K' .or. s == 'k') then
      v => klat
   else
      call ErrorHandler('getNumLatticeVectors','Invalid space type',s)
   endif
!
   end function getLatticeVectors
!  =====================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calLatticeSystem(a0,b0,c0,cosa,cosb,cosg,latsys)
!  =====================================================================
   implicit none
!
   character (len=*), intent(out) :: latsys
!
   real (kind=RealKind), intent(in) :: a0, b0, c0, cosa, cosb, cosg
!
   if (abs(cosa) < tol .and. abs(cosb) < tol .and. abs(cosg) < tol) then
      if (abs(a0-b0) < tol .and. abs(a0-c0) < tol) then
         latsys = 'Cubic'
      else if (abs(a0-b0) < tol .or. abs(a0-c0) < tol .or. abs(b0-c0) < tol) then
         latsys = 'Tetragonal'
      else
         latsys = 'Orthorhombic'
      endif
   else if (abs(cosa) < tol .and. abs(cosb) < tol) then
      if (abs(a0-b0) < tol .and.                                        &
          (abs(cosg-HALF) < tol .or. abs(cosg+HALF) < tol)) then
         latsys = 'Hexagonal'
      else 
         latsys = 'Monoclinic'
      endif
   else if (abs(cosa) < tol .and. abs(cosg) < tol) then
      if (abs(a0-c0) < tol .and.                                        &
          (abs(cosb-HALF) < tol .or. abs(cosb+HALF) < tol)) then
         latsys = 'Hexagonal'
      else 
         latsys = 'Monoclinic'
      endif
   else if (abs(cosb) < tol .and. abs(cosg) < tol) then
      if (abs(b0-c0) < tol .and.                                        &
          (abs(cosa-HALF) < tol .or. abs(cosa+HALF) < tol)) then
         latsys = 'Hexagonal'
      else 
         latsys = 'Monoclinic'
      endif
   else if (abs(a0-b0) < tol .and. abs(a0-c0) < tol) then
      if (abs(cosa-cosb) < tol .and. abs(cosa-cosg) < tol) then
         if (abs(abs(cosa)-HALF) < tol .or. abs(cosa+THIRD) < tol) then
            latsys = 'Cubic'
         else
            latsys = 'Rhombohedral'
         endif
      else if (abs(cosa-cosb) < tol .or. abs(cosa-cosg) < tol .or. abs(cosb-cosg) < tol) then
         latsys = 'Tetragonal'
      else
         latsys = 'Rhombohedral'
      endif
   else if (abs(cosa-cosb) < tol .and. abs(cosa-cosg) < tol) then
      if (abs(a0-b0) < tol .or. abs(a0-c0) < tol .or. abs(b0-c0) < tol) then
         latsys = 'Tetragonal'
      else
         latsys = 'Orthorhombic'
      endif
   else
      latsys = 'Triclinic'
   endif
!
   if (MyPE == 0) then
      write(6,'(/,80(''=''),/)')
      write(6,'(a,a)')' The lattice system is: ',latsys
   endif
!
   end subroutine calLatticeSystem
!  =====================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genKBravais()
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: vol, fac
!
   vol=(rbrav(2,1)*rbrav(3,2)-rbrav(3,1)*rbrav(2,2))*rbrav(1,3)+          &
       (rbrav(3,1)*rbrav(1,2)-rbrav(1,1)*rbrav(3,2))*rbrav(2,3)+          &
       (rbrav(1,1)*rbrav(2,2)-rbrav(2,1)*rbrav(1,2))*rbrav(3,3)
   vol=abs(vol)

   fac=PI2/vol
   kbrav(1,1)=fac*(rbrav(2,2)*rbrav(3,3)-rbrav(3,2)*rbrav(2,3))
   kbrav(2,1)=fac*(rbrav(3,2)*rbrav(1,3)-rbrav(1,2)*rbrav(3,3))
   kbrav(3,1)=fac*(rbrav(1,2)*rbrav(2,3)-rbrav(2,2)*rbrav(1,3))
   kbrav(1,2)=fac*(rbrav(2,3)*rbrav(3,1)-rbrav(3,3)*rbrav(2,1))
   kbrav(2,2)=fac*(rbrav(3,3)*rbrav(1,1)-rbrav(1,3)*rbrav(3,1))
   kbrav(3,2)=fac*(rbrav(1,3)*rbrav(2,1)-rbrav(2,3)*rbrav(1,1))
   kbrav(1,3)=fac*(rbrav(2,1)*rbrav(3,2)-rbrav(3,1)*rbrav(2,2))
   kbrav(2,3)=fac*(rbrav(3,1)*rbrav(1,2)-rbrav(1,1)*rbrav(3,2))
   kbrav(3,3)=fac*(rbrav(1,1)*rbrav(2,2)-rbrav(2,1)*rbrav(1,2))
!
   end subroutine genKBravais
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genLattice(s,vcut,nm1_in,nm2_in,nm3_in)
!  ===================================================================
   use VectorModule, only : getVecLength
   implicit none
!
   character (len=1), intent(in) :: s
!
   real (kind=RealKind), intent(in) :: vcut
!
   integer (kind=IntKind), intent(in), optional :: nm1_in
   integer (kind=IntKind), intent(in), optional :: nm2_in
   integer (kind=IntKind), intent(in), optional :: nm3_in
!
   real (kind=RealKind), pointer :: vbra(:,:)
   real (kind=RealKind), pointer :: vlat(:,:)
   real (kind=RealKind), allocatable :: vlat_1(:)
   real (kind=RealKind), allocatable :: vlat_2(:)
   real (kind=RealKind), allocatable :: vlat_3(:)
   real (kind=RealKind), allocatable :: vsq(:)
   real (kind=RealKind) :: vn(3)
   real (kind=RealKind) :: vnsq
   real (kind=RealKind) :: vcut2
!
   integer (kind=IntKind), pointer :: nv
   integer (kind=IntKind) :: nm1, nm2, nm3, n1, n2, n3, i
   integer (kind=IntKind) :: ipmax
   integer (kind=IntKind) :: tn1p1, tn2p1, tn3p1
   integer (kind=IntKind) :: nt12, nt23, nt123
!
   if (s == 'R' .or. s == 'r') then
      vbra => rbrav
      vlat => rlat
      nv => num_rlat
   else if (s == 'K' .or. s == 'k') then
      vbra => kbrav
      vlat => klat
      nv => num_klat
   else
      call ErrorHandler('getNumLatticeVectors','Invalid space type',s)
   endif
!
   if (.not.present(nm1_in)) then
      nm1 = ceiling(vcut/getVecLength(3,vbra(1:3,1)))
      nm2 = ceiling(vcut/getVecLength(3,vbra(1:3,2)))
      nm3 = ceiling(vcut/getVecLength(3,vbra(1:3,3)))
   else if (.not.present(nm2_in)) then
      nm1 = nm1_in
      nm2 = ceiling(vcut/getVecLength(3,vbra(1:3,2)))
      nm3 = ceiling(vcut/getVecLength(3,vbra(1:3,3)))
   else if (.not.present(nm3_in)) then
      nm1 = nm1_in
      nm2 = nm2_in
      nm3 = ceiling(vcut/getVecLength(3,vbra(1:3,3)))
   else
      nm1 = nm1_in
      nm2 = nm2_in
      nm3 = nm3_in
   endif
!
   ipmax = (2*nm1+1)*(2*nm2+1)*(2*nm3+1)
   allocate(vlat_1(ipmax), vlat_2(ipmax), vlat_3(ipmax), vsq(ipmax))
!
!  ===================================================================
!  generate lattice vectors........................................
!  ===================================================================
   vlat_1 = ZERO
   vlat_2 = ZERO
   vlat_3 = ZERO
   vsq = ZERO
!
   nv=0
   vcut2=vcut*vcut+tol
   tn1p1=2*nm1+1
   tn2p1=2*nm2+1
   tn3p1=2*nm3+1
   nt12=tn1p1*tn2p1
   nt23=tn2p1*tn3p1
   nt123=nt12*tn3p1
   do i=1,nt123
      n1=i-1
      n1=mod(n1,tn1p1)-nm1
      n2=(i-1)/tn1p1
      n2=mod(n2,tn2p1)-nm2
      n3=(i-1)/nt12
      n3=mod(n3,tn3p1)-nm3
      vn(1) = n1*vbra(1,1)+n2*vbra(1,2)+n3*vbra(1,3)
      vn(2) = n1*vbra(2,1)+n2*vbra(2,2)+n3*vbra(2,3)
      vn(3) = n1*vbra(3,1)+n2*vbra(3,2)+n3*vbra(3,3)
      vnsq=vn(1)*vn(1)+vn(2)*vn(2)+vn(3)*vn(3)
      if(vnsq.le.vcut2) then
         if(nv+1.gt.ipmax) then
            call ErrorHandler('genLattice','nv > ipmax',nv,ipmax)
         endif
!        -------------------------------------------------------------
         call ord3v(vlat_1,vlat_2,vlat_3,vsq,nv,vn,vnsq)
!        -------------------------------------------------------------
      endif
   enddo
!
   do i = 1, nv
     vlat(1,i) = vlat_1(i)
     vlat(2,i) = vlat_2(i)
     vlat(3,i) = vlat_3(i)
   enddo
!
   deallocate(vlat_1, vlat_2, vlat_3, vsq)
!
   end subroutine genLattice
!  ===================================================================
end module LatticeModule
