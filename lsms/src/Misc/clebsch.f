      subroutine clebsch
c=======================
      implicit real*8 (a-h,o-z)
c
      include 'cgc.h'
!     dimension u1(50),ind1(50)
!     dimension u2(50),ind2(50)
!     common/cgc/u1,u2,ind1,ind2
      data tiny/1.0d-15/
c
c clebsch-gordan rectangular matrices to transform from (lm) to
c (kappa,my) basis
c
      do 1 i=1,500
      u1(i)=0.d0
      u2(i)=0.d0
      ind1(i)=1
      ind2(i)=1
    1 continue

      inr=0
      do 33 l=0,12
         twolp1=dfloat(2*l+1)
         do 3 m=-l,l
            inr=inr+1
!           write(6,*) 'inr',inr
c
c j=l-1/2
c
            kap=l
            if(kap.eq.0) goto 2
c
c ms=-1/2
c
            ir=2*kap*kap+kap+m
            u1(ir)=dsqrt((l+m)/twolp1)
            ind1(ir)=inr
c
c ms=+1/2
c
            ir=2*kap*kap+kap+m+1
            u2(ir)=-dsqrt((l-m)/twolp1)
            ind2(ir)=inr
    2       continue
c
c j=l+1/2
c
            kap=-l-1
c
c ms=-1/2
c
            ir=2*kap*kap+kap+m
            u1(ir)=dsqrt((l-m+1)/twolp1)
            ind1(ir)=inr
c
c ms=+1/2
c
            ir=2*kap*kap+kap+m+1
            u2(ir)=dsqrt((l+m+1)/twolp1)
            ind2(ir)=inr
c
    3 continue
   33 continue
c
c      write(6,*)
       do ir=1,inr
         if(dabs(u1(ir)).lt.tiny) ind1(ir)=1
c        write(6,'(i3,2x,i3,2x,f20.10)') ir,ind1(ir),u1(ir)
       end do
c      write(6,*)
       do ir=1,inr
         if(dabs(u2(ir)).lt.tiny) ind2(ir)=1
c        write(6,'(i3,2x,i3,2x,f20.10)') ir,ind2(ir),u2(ir)
       end do
c
      return
      end
      subroutine ruthi(kap,my,ii,norder)
c=====================
c     find index for (kapa-my)-representation
c
      dimension kapdex(50),mydex(50)
      data kapdex/-1,-1,1,1,-2,-2,-2,-2,2,2,2,2,-3,-3,-3,-3,-3,-3
     * ,3,3,3,3,3,3,-4,-4,-4,-4,-4,-4,-4,-4,
     * 4,4,4,4,4,4,4,4,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5/
      data mydex/-1,1,-1,1,-3,-1,1,3,-3,-1,1,3,-5,-3,-1,1,3,5,
     * -5,-3,-1,1,3,5,-7,-5,-3,-1,1,3,5,7,
     * -7,-5,-3,-1,1,3,5,7,-9,-7,-5,-3,-1,1,3,5,7,9/
c
      do 100 i=1,norder
      if(kapdex(i).eq.kap) go to 150
  100 continue
      go to 500
  150 do 170 j=i,norder
      if(mydex(j).eq.my) go to 200
  170 continue
      go to 500
c
  200 ii=j
      return
  500 stop
      end
