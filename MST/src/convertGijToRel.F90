!====================================================================
subroutine convertGijToRel(gij, bgij, kkr1, kkr2, ce)
!====================================================================

! Made from setgij and relmtrx from LSMS 1.9

    use KindParamModule, only : IntKind, RealKind, CmplxKind
    use PhysParamModule, only : LightSpeed
    use MathParamModule, only : CZERO

    implicit none

    character (len=32), parameter :: sname = 'convertGijToRel'

    integer (kind=IntKind) :: im, in
    integer (kind=IntKind) :: i, j
    integer (kind=IntKind) :: i1, i2, j1, j2
    integer (kind=IntKind) :: is
    integer (kind=IntKind), intent(in) :: kkr1
    integer (kind=IntKind) :: kkr1_ns
    integer (kind=IntKind), intent(in) :: kkr2
    integer (kind=IntKind) :: kkr2_ns
    integer (kind=IntKind) :: ind1(500), ind2(500)

    real (kind=RealKind) :: u1(500), u2(500)

    complex (kind=CmplxKind), intent(in) :: gij(:,:)
    complex (kind=CmplxKind), intent(out) :: bgij(:,:)
    complex (kind=CmplxKind), intent(in) :: ce
    complex (kind=CmplxKind) :: fac
    complex (kind=CmplxKind) :: psq

    kkr1_ns = 2*kkr1
    kkr2_ns = 2*kkr2

    psq = ce + ce*ce/(LightSpeed*LightSpeed)

    call clebsch(u1,u2,ind1,ind2,500)

    bgij = CZERO

! Begin relmtrx

    do j=1,kkr2_ns
       j1=ind1(j)
       j2=ind2(j)
       do i=1,kkr1_ns
          i1=ind1(i)
          i2=ind2(i)
          bgij(i,j)=u1(i)*gij(i1,j1)*u1(j)+ &
                    u2(i)*gij(i2,j2)*u2(j)
       end do
    end do

! End relmtrx

    fac=psq/ce
    do j=1,kkr2_ns
       do i=1,kkr1_ns
          bgij(i,j)=fac*bgij(i,j)
       end do
    end do

end subroutine convertGijToRel

!=================================================================
   subroutine clebsch(u1,u2,ind1,ind2,n)
!=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO
   implicit none

   integer (kind=IntKind) :: i, l, m, inr, kap, ir, n

   real (kind=RealKind), parameter :: tiny = 1.0d-15
   real (kind=RealKind) :: twolp1
   real (kind=RealKind):: u1(n),u2(n)
   integer (kind=IntKind) :: ind1(n),ind2(n)

!  clebsch-gordan rectangular matrices to transform from (lm) to
! (kappa,my) basis

   do 101 i=1,500
      u1(i)=ZERO
      u2(i)=ZERO
      ind1(i)=1
      ind2(i)=1
   101 continue

   inr=0
   do 103 l=0,12
      twolp1=dfloat(2*l+1)
      do m=-l,l
         inr=inr+1

! j=l-1/2
         kap=l
         if(kap.eq.0) goto 102

! ms=-1/2
         ir=2*kap*kap+kap+m
         u1(ir)=sqrt((l+m)/twolp1)
         ind1(ir)=inr

! ms=+1/2
         ir=2*kap*kap+kap+m+1
         u2(ir)=-sqrt((l-m)/twolp1)
         ind2(ir)=inr
    102  continue

! j=l+1/2
         kap=-l-1

! ms=-1/2
         ir=2*kap*kap+kap+m
         u1(ir)=sqrt((l-m+1)/twolp1)
         ind1(ir)=inr

! ms=+1/2
         ir=2*kap*kap+kap+m+1
         u2(ir)=sqrt((l+m+1)/twolp1)
         ind2(ir)=inr
      enddo
  103 continue

  do ir=1,inr
     if(abs(u1(ir)).lt.tiny) ind1(ir)=1
  end do
  do ir=1,inr
     if(abs(u2(ir)).lt.tiny) ind2(ir)=1
  end do

   end subroutine clebsch
