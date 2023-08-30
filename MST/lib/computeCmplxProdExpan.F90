!  **********************************************************************
!  *                                                                    *
!  *  This subroutine aims to compute the spherical harmonic            *
!  *  expansion of the product of two complex functions, f and g,       *
!  *  which are expanded in spherical harmonics up to lmax_f and lmax_g,*
!  *  respectively. The product function h = f*g is expanded up to      *
!  *  lmax_h. There are three situations, depending on act symbol:      *
!  *                                                                    *
!  *  If act = 'n' or 'N', the product h = f*g, and its expansion is    *
!  *     calculated as follows                                          *
!  *                            L"                                      *
!  *             h (r) = sum   C    * f (r) * g (r)                     *
!  *              L      L',L"  L',L   L'      L"                       *
!  *                                                                    *
!  *  If act = 's' or 'S', the product h = f^{s}*g, and its expansion   *
!  *     is calculated as follows                                       *
!  *                            L       s                               *
!  *             h (r) = sum   C     * f (r) * g (r)                    *
!  *              L      L',L"  L',L"   L'      L"                      *
!  *  Note:                                                             *
!  *         s ->         s     * ^     ->                 ^            *
!  *        f (r ) = sum f (r)*Y (r), g(r ) = sum g (r)*Y (r)           *
!  *                  L   L     L              L   L     L              *
!  *                                                                    *
!  *  If act = 'c' or 'C', the product h = f^{*}*g, where {*} is only   *
!  *     allied to the spherical harmonics, and the expansion of h is   *
!  *     calculated as follows                                          *
!  *                            L                                       *
!  *             h (r) = sum   C     * f (r) * g (r)                    *
!  *              L      L',L"  L',L"   L'      L"                      *
!  *                                                                    *
!  **********************************************************************
subroutine computeCmplxProdExpan(nr,kmax_f,f,act,kmax_g,g,kmax_h,h)
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
   use IntegerFactorsModule, only : lofk, mofk, kofj, m1m
!
   implicit none
!
   character (len=1), intent(in) :: act
   integer (kind=intKind), intent(in) :: nr, kmax_f, kmax_g, kmax_h
   integer (kind=intKind) :: ir, j, m_f, m_g, kl_h, kl_f, kl_g, kl_f_bar
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
      call ErrorHandler('computeCmplxProdExpan','size(f,1) < nr',size(f,1),nr) 
   else if (size(g,1) < nr) then
      call ErrorHandler('computeCmplxProdExpan','size(g,1) < nr',size(g,1),nr) 
   else if (size(h,1) < nr) then
      call ErrorHandler('computeCmplxProdExpan','size(h,1) < nr',size(h,1),nr) 
   endif 
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   h = CZERO
   if (act == 'n' .or. act == 'N') then
      do kl_h = 1, kmax_h
         do kl_f = 1, kmax_f
            do j = 1, nj3(kl_f,kl_h)
               kl_g = kj3(j,kl_f,kl_h)
               if (kl_g <= kmax_g) then
                  cfac = cgnt(j,kl_f,kl_h)
                  do ir = 1, nr
                     h(ir,kl_h) = h(ir,kl_h) + cfac*f(ir,kl_f)*g(ir,kl_g)
                  enddo
               endif
            enddo
         enddo
      enddo
   else if (act == 's' .or. act == 'S') then
      do kl_g = 1, kmax_g
         do kl_f = 1, kmax_f
            m_f = mofk(kl_f)
            kl_f_bar = kl_f - 2*m_f
            do j = 1, nj3(kl_f,kl_g)
               kl_h = kj3(j,kl_f,kl_g)
               if (kl_h <= kmax_h) then
                  cfac = m1m(m_f)*cgnt(j,kl_f,kl_g)
                  do ir = 1, nr
                     h(ir,kl_h) = h(ir,kl_h) + cfac*f(ir,kl_f_bar)*g(ir,kl_g)
                  enddo
               endif
            enddo
         enddo
      enddo
   else if (act == 'c' .or. act == 'C') then
      do kl_g = 1, kmax_g
         do kl_f = 1, kmax_f
            do j = 1, nj3(kl_f,kl_g)
               kl_h = kj3(j,kl_f,kl_g)
               if (kl_h <= kmax_h) then
                  cfac = cgnt(j,kl_f,kl_g)
                  do ir = 1, nr
                     h(ir,kl_h) = h(ir,kl_h) + cfac*f(ir,kl_f)*g(ir,kl_g)
                  enddo
               endif
            enddo
         enddo
      enddo
   else
      call ErrorHandler('computeCmplxProdExpan','Invaid action symbol',act)
   endif
!
end subroutine computeCmplxProdExpan
