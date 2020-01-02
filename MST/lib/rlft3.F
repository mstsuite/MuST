      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      use KindParamModule, only : IntKind, RealKind, CmplxKind
      use MathParamModule, only : PI2, HALF, ZERO, TWO, ONE
      INTEGER (kind=IntKind) :: isign,nn1,nn2,nn3
      COMPLEX (kind=CmplxKind), target :: data(nn1/2,nn2,nn3)
      COMPLEX (kind=CmplxKind) ::speq(nn2,nn3)
CU    USES fourn
      INTEGER (kind=IntKind) :: i1,i2,i3,j1,j2,j3,nn(3)
      REAL (kind=RealKind) :: theta,wi,wpi,wpr,wr,wtemp
      COMPLEX (kind=CmplxKind) :: c1,c2,h1,h2,w
c
      c1=cmplx(HALF,ZERO,kind=CmplxKind)
      c2=cmplx(ZERO,-HALF*isign,kind=CmplxKind)
c     theta=6.28318530717959d0/real(isign*nn1,kind=RealKind)
      theta=PI2/real(isign*nn1,kind=RealKind)
      wpr=-TWO*sin(HALF*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=ONE
        wi=ZERO
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=(wr*wpr-wi*wpi)+wr
          wi=(wi*wpr+wtemp*wpi)+wi
          w=cmplx(wr,wi,kind=CmplxKind)
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END
