      SUBROUTINE fourn(data,nn,ndim,isign)
      use KindParamModule, only : IntKind, RealKind, CmplxKind
      use MathParamModule, only : PI2, ZERO, HALF, ONE, TWO
      INTEGER (kind=IntKind) :: isign,ndim,nn(ndim)
      REAL (kind=RealKind) :: data(*)
      INTEGER (kind=IntKind) :: i1,i2,i2rev,i3,i3rev,ibit,idim
      INTEGER (kind=IntKind) :: ifp1,ifp2,ip1,ip2,ip3,k1
      INTEGER (kind=IntKind) :: k2,n,nprev,nrem,ntot
      REAL (kind=RealKind) :: tempi,tempr
      REAL (kind=RealKind) :: theta,wi,wpi,wpr,wr,wtemp
c
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
c         theta=isign*6.28318530717959d0/(ifp2/ip1)
          theta=isign*PI2/(ifp2/ip1)
          wpr=-TWO*sin(HALF*theta)**2
          wpi=sin(theta)
          wr=ONE
          wi=ZERO
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=wr*data(k2)-wi*data(k2+1)
                tempi=wr*data(k2+1)+wi*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=(wr*wpr-wi*wpi)+wr
            wi=(wi*wpr+wtemp*wpi)+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END
