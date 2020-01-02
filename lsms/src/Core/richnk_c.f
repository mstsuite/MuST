c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine richnk (nn,y,bh,dh)
c     ================================================================
c
c     ****************************************************************
c     calculates the riccati-bessel functions hl.
c     this version calculates the derivatives.
c     ****************************************************************
c 
      implicit   none
c
      integer    l
      integer    nn
c
      complex*16 bh(nn)
      complex*16 dh(nn)
      complex*16 y
c
      complex*16 x
      complex*16 xinv
      complex*16 xpon
      complex*16 one
      complex*16 sqrtm1
      complex*16 zero
c
      parameter (zero=(0.d0,0.d0))
      parameter (one=(1.d0,0.d0))
      parameter (sqrtm1=(0.d0,1.d0))
c
      x=y
      if (abs(x).eq.0.d00) then
         write(6,'('' trouble in richnk argument ='',2f10.5)') x
	 call fstop('richnk')
      endif
c
c     ================================================================
c     recursion relations
c     ================================================================
      xpon=exp( sqrtm1*x )
      xinv=one/x
      bh(1)=-sqrtm1
      bh(2)= - one - sqrtm1*xinv 
      dh(1)= one
      dh(2)= sqrtm1 * ( xinv - ( sqrtm1 + x ) ) * xinv
c
c     ================================================================
c     flm=2l+1 for real l=0,1,2,3, ...  quantum number
c     ================================================================
      do l=3,nn
         bh(l) = (2*l-3)*bh(l-1)*xinv - bh(l-2)
         dh(l) = bh(l-1) - (l-1)*bh(l)*xinv
      enddo 
      do l=1,nn
         bh(l)=bh(l)*xpon
         dh(l)=dh(l)*xpon
      enddo
c
      if( abs(x) .le. 0.01d0 ) then
c        =============================================================
c        power-series for j if abs(x) is smaller than 0.9. trouble!!!
c        =============================================================
         write(6,'('' trouble in richnk: small argument'')')
         call fstop('richnk')
      endif
c
      return
      end
