!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine numlat(vbra,vcut,nm1,nm2,nm3,nv)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : TEN2m6
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: nm1
   integer (kind=IntKind), intent(in) :: nm2
   integer (kind=IntKind), intent(in) :: nm3
   integer (kind=IntKind), intent(out) :: nv
   integer (kind=IntKind) :: tn1p1
   integer (kind=IntKind) :: tn2p1
   integer (kind=IntKind) :: tn3p1
   integer (kind=IntKind) :: nt12
   integer (kind=IntKind) :: nt23
   integer (kind=IntKind) :: nt123
   integer (kind=IntKind) :: n1
   integer (kind=IntKind) :: n2
   integer (kind=IntKind) :: n3
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: vbra(3,3)
   real (kind=RealKind) :: vn(3)
   real (kind=RealKind) :: vnsq
   real (kind=RealKind) :: vcut
   real (kind=RealKind) :: vcut2
!
!  *******************************************************************
!  Calculates number of lattice vectors, for the given cutoff radius..
!  *******************************************************************
   nv=0
   vcut2=vcut*vcut+TEN2m6
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
         nv = nv + 1
      endif
   enddo
   end subroutine numlat
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine lattice(vbra,vcut,nm1,nm2,nm3,vlat_1,vlat_2,vlat_3,vsq,nv, &
                      ipmax)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, TEN2m6
!
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: ipmax
   integer (kind=IntKind), intent(in) :: nm1
   integer (kind=IntKind), intent(in) :: nm2
   integer (kind=IntKind), intent(in) :: nm3
   integer (kind=IntKind), intent(out) :: nv
   integer (kind=IntKind) :: tn1p1
   integer (kind=IntKind) :: tn2p1
   integer (kind=IntKind) :: tn3p1
   integer (kind=IntKind) :: nt12
   integer (kind=IntKind) :: nt23
   integer (kind=IntKind) :: nt123
   integer (kind=IntKind) :: n1
   integer (kind=IntKind) :: n2
   integer (kind=IntKind) :: n3
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: vbra(3,3)
   real (kind=RealKind), intent(out) :: vlat_1(ipmax)
   real (kind=RealKind), intent(out) :: vlat_2(ipmax)
   real (kind=RealKind), intent(out) :: vlat_3(ipmax)
   real (kind=RealKind), intent(out) :: vsq(ipmax)
   real (kind=RealKind) :: vn(3)
   real (kind=RealKind) :: vnsq
   real (kind=RealKind) :: vcut
   real (kind=RealKind) :: vcut2
!
!  *******************************************************************
!     Calculates lattice vectors .....................................
!     input:
!       vbra     : basis vectors for bravais lattice
!       vcut     : cut off radius for lattice vectors
!       nm1      : number of repeats of vbra(1)
!       nm2      : number of repeats of vbra(2)
!       nm3      : number of repeats of vbra(3)
!       ipmax    : dimension of array that will contain vectors
!     output:
!       vlat     : bravais lattice vectors
!       vsq      : square of length of lattice vectors
!       nv       : number of lattice vectors
!
!  *******************************************************************
!
!  ===================================================================
!  generate lattice vectors........................................
!  ===================================================================
   vlat_1(1:ipmax) = ZERO
   vlat_2(1:ipmax) = ZERO
   vlat_3(1:ipmax) = ZERO
   vsq(1:ipmax) = ZERO
!
   nv=0
   vcut2=vcut*vcut+TEN2m6
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
            write(6,'('' lattice:: nv.gt.ipmax:'',2i5)')nv,ipmax
            call ErrorHandler('lattice','nv > ipmax',nv,ipmax)
         endif
!        -------------------------------------------------------------
         call ord3v(vlat_1,vlat_2,vlat_3,vsq,nv,vn,vnsq)
!        -------------------------------------------------------------
      endif
   enddo
!
   end subroutine lattice
