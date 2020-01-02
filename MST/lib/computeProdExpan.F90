!  *******************************************************************
!  *                                                                 *
!  *  This subroutine aims to compute the spherical harmonic         *
!  *  expansion of the product of two real functions, f and g, which *
!  *  are expanded in spherical harmonics up to lmax_f and lmax_g,   *
!  *  respectively. The product function h = f*g is expanded up to   *
!  *  lmax_h.                                                        *
!  *                                                                 *
!  *  Note: jmax = (lmax+1)*(lmax+2)/2. It is easy to show that      *
!  *                                                                 *
!  *                L'                          L"     *             *
!  *  h (r) = sum  C    * f (r) * g (r) = sum  C    * f (r) * g (r)  *
!  *   L"     L,L'  L,L"   L       L'     L,L'  L,L'   L       L'    *
!  *                                                                 *
!  *******************************************************************
subroutine computeProdExpan(nr,lmax_f,f,lmax_g,g,lmax_h,h)
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : CZERO, CONE
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
   use IntegerFactorsModule, only : lofk, mofk, jofk, kofj, m1m
!
   implicit none
!
   integer (kind=intKind), intent(in) :: nr, lmax_f, lmax_g, lmax_h
   integer (kind=intKind) :: kmax_f, kmax_g, kmax_h, jl_f, jl_g, jl_h, ir, j, m_f, m_g, kl_h, kl_f, kl_g
   integer (kind=IntKind), pointer :: nj3(:,:)
   integer (kind=IntKind), pointer :: kj3(:,:,:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
!
   complex (kind=CmplxKind), intent(in) :: f(:,:), g(:,:)
   complex (kind=CmplxKind), intent(out) :: h(:,:)
   complex (kind=CmplxKind) :: cfac
!
   if (size(f,1) < nr) then
      call ErrorHandler('computeProdExpan','size(f,1) < nr',size(f,1),nr) 
   else if (size(f,2) < (lmax_f+1)*(lmax_f+2)/2) then
      call ErrorHandler('computeProdExpan','size(f,2) < jmax_f',size(f,2),(lmax_f+1)*(lmax_f+2)/2) 
   else if (size(g,1) < nr) then
      call ErrorHandler('computeProdExpan','size(g,1) < nr',size(g,1),nr) 
   else if (size(g,2) < (lmax_g+1)*(lmax_g+2)/2) then
      call ErrorHandler('computeProdExpan','size(g,2) < jmax_g',size(g,2),(lmax_g+1)*(lmax_g+2)/2) 
   else if (size(h,1) < nr) then
      call ErrorHandler('computeProdExpan','size(h,1) < nr',size(h,1),nr) 
   else if (size(h,2) < (lmax_h+1)*(lmax_h+2)/2) then
      call ErrorHandler('computeProdExpan','size(h,2) < jmax_h',size(h,2),(lmax_h+1)*(lmax_h+2)/2) 
   endif 
!
   kmax_f = (lmax_f+1)**2
   kmax_g = (lmax_g+1)**2
   kmax_h = (lmax_h+1)**2
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   h = CZERO
   do kl_g = 1, kmax_g
      m_g = mofk(kl_g)
      jl_g = jofk(kl_g)
      do kl_f = 1, kmax_f
         m_f = mofk(kl_f)
         jl_f = jofk(kl_f)
         do j = 1, nj3(kl_f,kl_g)
            kl_h = kj3(j,kl_f,kl_g)
            if (mofk(kl_h) >= 0 .and. kl_h <= kmax_h) then
               jl_h = jofk(kl_h)
               if (m_f >= 0 .and. m_g >= 0) then
                  cfac = cgnt(j,kl_f,kl_g)
                  do ir = 1, nr
                     h(ir,jl_h) = h(ir,jl_h) + cfac*conjg(f(ir,jl_f))*g(ir,jl_g)
                  enddo
               else if (m_f < 0 .and. m_g >= 0) then
                  cfac = cgnt(j,kl_f,kl_g)*m1m(m_f)
                  do ir = 1, nr
                     h(ir,jl_h) = h(ir,jl_h) + cfac*f(ir,jl_f)*g(ir,jl_g)
                  enddo
               else if (m_f >= 0 .and. m_g < 0) then
                  cfac = cgnt(j,kl_f,kl_g)*m1m(m_g)
                  do ir = 1, nr
                     h(ir,jl_h) = h(ir,jl_h) + cfac*conjg(f(ir,jl_f)*g(ir,jl_g))
                  enddo
               else
                  cfac = cgnt(j,kl_f,kl_g)*m1m(m_f+m_g)
                  do ir = 1, nr
                     h(ir,jl_h) = h(ir,jl_h) + cfac*f(ir,jl_f)*conjg(g(ir,jl_g))
                  enddo
               endif
            endif
         enddo
      enddo
   enddo
!
end subroutine computeProdExpan
