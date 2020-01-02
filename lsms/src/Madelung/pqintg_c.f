c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pqintg(eta,rssq,lastrs,dqint,ibegin,iprslat_mad)
c     ================================================================
c
c     ****************************************************************
c     input: eta,rssq,lastrs,ibegin
c     output: dqint
c     ****************************************************************
c
      implicit   none
c
      integer    iprslat_mad
      integer    ibegin
      integer    lastrs
      integer    igns
      integer    i
c
      real*8     dqint(iprslat_mad)
      real*8     rssq(iprslat_mad)
      real*8     eta
      real*8     delta
      real*8     xi
      real*8     xisq
      real*8     arg
      real*8     sum
      real*8     term
      real*8     zero
      real*8     half
      real*8     three
      real*8     factor
      real*8     tol
c
      parameter (zero=0.0d0)
      parameter (factor=0.005d0)
      parameter (half=0.5d0)
      parameter (three=3.0d0)
      parameter (tol=0.000000001d0)
c
c     ================================================================
      do i = ibegin,lastrs
c        delta = 5.e-3*eta
         delta = factor*eta
         xi = sqrt(eta)*half
         xisq = xi*xi
         arg = - rssq(i)*xisq
         sum = exp(arg)
c
         if( sum .ne. zero ) then
            igns = -1
 1          continue
            xi = xi + delta
            xisq = xi * xi
            arg = - rssq(i)*xisq
            igns = -igns
            term = ( 3 + igns )*exp(arg)
            sum = sum + term
            if( term/sum .gt. tol ) go to 1
            dqint(i) = ( sum*delta )/three
         else
            dqint(i) = zero
         endif
      enddo
c
      return
      end
