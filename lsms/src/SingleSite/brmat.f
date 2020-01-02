      subroutine brmat(l,my,n,bspr,bopr,
     >                 ba11,ba12,ba22,bb11,bb22,iprpts)
c=====================
c
      implicit real*8 (a-h,o-z)
!     include '../param.h'                 
!     include 'atom_param.h'
      include '../Misc/cgc.h'
c
      dimension bspr(iprpts),bopr(iprpts,2)
      dimension ba11(iprpts),ba12(iprpts),ba22(iprpts)
      dimension bb11(iprpts),bb22(iprpts),br(iprpts,2)
c
!     dimension u1(50),ind1(50)
!     dimension u2(50),ind2(50)
!     common/cgc/u1,u2,ind1,ind2
c
      call zeroout(ba11,iprpts)
      call zeroout(ba12,iprpts)
      call zeroout(ba22,iprpts)
      call zeroout(bb11,iprpts)
      call zeroout(bb22,iprpts)
c
      do i=1,n
        br(i,1)=-bspr(i)
        br(i,2)=bspr(i)
      end do
      if(l.eq.2) then
        do i=1,n
          br(i,1)=br(i,1)+bopr(i,1)*(my+1.0d0)/2.0d0
          br(i,2)=br(i,2)+bopr(i,2)*(my-1.0d0)/2.0d0
        end do
      end if
c
      IF(iabs(my).eq.2*l+1) THEN
c
        kap2=-l-1
        kapb2=l+1
        kap2my=2*kap2*kap2+kap2+(my+1)/2
        kapb2my=2*kapb2*kapb2+kapb2+(my+1)/2
        ca22_d=u1(kap2my)*u1(kap2my)
        ca22_u=u2(kap2my)*u2(kap2my)
        cb22_d=u1(kapb2my)*u1(kapb2my)
        cb22_u=u2(kapb2my)*u2(kapb2my)
        do i=1,n
          ba22(i)=ca22_d*br(i,1)+ca22_u*br(i,2)
          bb22(i)=cb22_d*br(i,1)+cb22_u*br(i,2)
        end do
c
      ELSE 
c
        kap1=l
        kap2=-l-1
        kapb1=-l 
        kapb2=l+1
        kap1my=2*kap1*kap1+kap1+(my+1)/2
        kapb1my=2*kapb1*kapb1+kapb1+(my+1)/2
        kap2my=2*kap2*kap2+kap2+(my+1)/2
        kapb2my=2*kapb2*kapb2+kapb2+(my+1)/2
        ca11_d=u1(kap1my)*u1(kap1my)
        ca11_u=u2(kap1my)*u2(kap1my)
        ca12_d=u1(kap1my)*u1(kap2my)
        ca12_u=u2(kap1my)*u2(kap2my)
        ca22_d=u1(kap2my)*u1(kap2my)
        ca22_u=u2(kap2my)*u2(kap2my)
        cb11_d=u1(kapb1my)*u1(kapb1my)
        cb11_u=u2(kapb1my)*u2(kapb1my)
        cb22_d=u1(kapb2my)*u1(kapb2my)
        cb22_u=u2(kapb2my)*u2(kapb2my)
        do i=1,n
          ba11(i)=ca11_d*br(i,1)+ca11_u*br(i,2)
          ba12(i)=ca12_d*br(i,1)+ca12_u*br(i,2)
          ba22(i)=ca22_d*br(i,1)+ca22_u*br(i,2)
          bb11(i)=cb11_d*br(i,1)+cb11_u*br(i,2)
          bb22(i)=cb22_d*br(i,1)+cb22_u*br(i,2)
        end do
c
      END IF
c
      return
      end
