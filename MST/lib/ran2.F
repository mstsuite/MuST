c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function  ran2(idum)
c     ================================================================
c
c     ****************************************************************
c
c     this is function ran2(idum)  NUMERICAL RECIPES 2nd Ed.
c     copied by hand !
c
c     ****************************************************************
c
      implicit  none
c
      integer   im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv,idum
      real*8    zam,eps,rnmx,ran2          ! Change to double precision
      parameter (im1=2147483563,im2=2147483399,zam=(1./im1),imm1=im1-1,
     >           ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     >           ir2=3791,ntab=32,ndiv=1+imm1/ntab,
     >           eps=1.2e-14,rnmx=1.-eps)
      integer   idum2,j,k,iv(ntab),iy
      save      iv,iy,idum2
      data      idum2/123456789/, iv/ntab*0/, iy/0/

      if (idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if (idum.lt.0) idum=idum+im1
            if (j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
       endif
       k=idum/iq1
       idum=ia1*(idum-k*iq1)-k*ir1
       if (idum.lt.0) idum=idum+im1
       k=idum2/iq2
       idum2=ia2*(idum2-k*iq2)-k*ir2
       if (idum2.lt.0) idum2=idum2+im2
       j=1+iy/ndiv
       iy=iv(j)-idum2
       iv(j)=idum
       if (iy.lt.1) iy=iy+imm1
       ran2=min(zam*iy,rnmx)
c
       return
       end
