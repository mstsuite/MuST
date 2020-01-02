!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValueAtPosi(r_mesh,posi,val_l) result(val)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TWO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use IntegerFactorsModule, only : lofj, kofj, mofj
!
   use SphericalHarmonicsModule, only : calYlm
!
   use InterpolationModule, only : PolyInterp
!
   implicit none
!
   real (kind=RealKind), intent(in) :: r_mesh(:)
   real (kind=RealKind), intent(in) :: posi(3)
!
   complex (kind=CmplxKind), intent(in) :: val_l(:,:)
!
   integer (kind=IntKind) :: jl, ir, irp, l, kl, iend, kmax, lmax, jmax
   integer (kind=IntKind), parameter :: n_inter = 5 ! order of polynomial
!
   real (kind=RealKind) :: val, err, r
!
   complex (kind=CmplxKind) :: val_in
   complex (kind=CmplxKind), allocatable :: ylm(:)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   iend = size(val_l,dim=1)
   jmax = size(val_l,dim=2)
   lmax = lofj(jmax)
   kmax = kofj(jmax)
   r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
!
   if (iend > size(r_mesh)) then
      call ErrorHandler('getValueAtPosi','iend > size of r_mesh',iend,size(r_mesh))
   else if (r > r_mesh(iend)) then
      val = ZERO
      return
   endif
!
   allocate( ylm(kmax) )
!  -------------------------------------------------------------------
   call calYlm(posi,lmax,ylm)
   call hunt(iend,r_mesh,r,ir)
!  -------------------------------------------------------------------
   if (ir > iend-(n_inter-1)/2) then
      irp=iend-n_inter+1
   else if (2*ir+1 > n_inter) then
      irp=ir-(n_inter-1)/2
   else
      irp=1
   endif
!
   val = ZERO
!
   do jl = 1, jmax
      l = lofj(jl)
!     ----------------------------------------------------------------
      call PolyInterp(n_inter, r_mesh(irp:irp+n_inter-1),             &
                      val_l(irp:irp+n_inter-1,jl), r, val_in, err)
!     ----------------------------------------------------------------
      kl = (l+1)*(l+1)-l+mofj(jl)
      if (mofj(jl) == 0) then
         val = val + real(val_in*ylm(kl),RealKind)
      else
         val = val + TWO*real(val_in*ylm(kl),RealKind)
      endif
   enddo
!
   deallocate( ylm )
!
   end function getValueAtPosi
