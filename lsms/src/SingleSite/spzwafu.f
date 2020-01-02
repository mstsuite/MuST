      subroutine spzwafu(socsc,ce,psq,l,my,vr,br,bopr,dx,ns,rs,
     >                   tminv,gz,fz,gj,fj,iflag,iprpts,iplmax)
c=======================
c
c input:  e - any complex energy
c         l,my - angular momentum indices
c         vr,br,bopr,dx,ns,rs,nws,rws - as in 'readpot'
c output: tminv - (l,my)-like 2x2 block of inverse t-matrix
c         gz - large component of regular radial solutions * r
c         fz - small component of regular radial solutions * r
c         gj - large component of irregular radial solutions * r
c         fj - small component of irregular radial solutions * r
c
      implicit real*8 (a-h,o-z)
c
      logical scale
c
!     include '../param.h'
!     include 'atom_param.h'

      character*32 sname
      character*32 istop

c
      integer iprpts, iplmax
      dimension vr(iprpts),br(iprpts),bopr(iprpts,2)
      complex*16 tminv(2,2)
      complex*16 fz(iprpts,2,2),gz(iprpts,2,2)
      complex*16 fj(iprpts,2,2),gj(iprpts,2,2)
      complex*16 f1(2,2,iprpts),g1(2,2,iprpts),gp1(2,2)
      complex*16 f2(2,2,iprpts),g2(2,2,iprpts),gp2(2,2)
      complex*16 f11(iprpts),g11(iprpts),f12(iprpts),g12(iprpts)
      complex*16 gp11,gp12
      complex*16 fb(0:iplmax+1),fn(0:iplmax+1),fh(0:iplmax+1)
      complex*16 sqrtm1,ce,psq,p,pl,sl,smlm1,detl
      complex*16 aa1,aa2,bb1,bb2,gg1,ff1,ggam1,gg2,ff2,ggam2,ggamb   
      complex*16 a1(2,2),a2(2,2),b1(2,2),b2(2,2)
      complex*16 g1mat(2,2),f1mat(2,2),gam1(2,2)
      complex*16 g2mat(2,2),f2mat(2,2),gam2(2,2)
      complex*16 jmat(2,2),jbmat(2,2),nmat(2,2),nbmat(2,2),gamb(2,2)
      complex*16 x1(2,2),x2(2,2),x3(2,2),x4(2,2)
c
      parameter (sname='spzwafu')

      data sqrtm1/(0.d0,1.d0)/
c
c     write(6,*) 
c     write(6,*) ' L,MU=:',l,my
c     write(6,*) 
c
!     write(6,*) 'entering ',sname
!     write(6,*) 'ns=',ns

      scale=dabs(1.0d0-socsc).gt.1.0d-3
!     write(6,*) 'SPZWAFU scale=',scale
c
      call zeroout(gz,8*iprpts)
      call zeroout(fz,8*iprpts)
      call zeroout(gj,8*iprpts)
      call zeroout(fj,8*iprpts)
      call zeroout(tminv,8)
      call zeroout(jmat,8)
      call zeroout(jbmat,8)
      call zeroout(nmat,8)
      call zeroout(nbmat,8)
      call zeroout(gamz,8)
      call zeroout(gamj,8)
      call zeroout(gamb,8)
c
      xnot=dlog(rs)-(ns-1)*dx
c
c---> c in rydberg units:
      c=274.072d0
c
      p=cdsqrt(psq)
      xl=dfloat(l)
      sl=1.d0
      smlm1=-1.d0
      sl=sl*ce/psq
      smlm1=smlm1*ce/psq
c
      call csbf(l+1,p,rs,fb,fn,fh)
c
      if(iabs(my) .gt. (2*l+1)) then
        stop ' spzwafu: my is out of range!'
      end if

!     if(iabs(my)-(2*l+1)) 1,2,3
!   3 stop ' spzwafu: my is out of range!'
c
!   2 continue
      if(iabs(my).eq.(2*l+1)) then

!     write(6,*) sname,2,1
c
      call dirmago1op(socsc,ce,l,my,vr,br,bopr,dx,xnot,rs,ns,
     >              g11,f11,gp11,iprpts)
!     write(6,*) sname,2,2
c     call dirmago1(socsc,ce,l,my,vr,br,dx,xnot,rs,ns,g11,f11,gp11)
!     write(6,*) sname,2,3
c
      if(scale) then
        gg1=g11(ns)/rs
        ff1=gp11/rs
c       write(6,'('' gmt reg  '',2d15.8)') gg1*rs
c       write(6,'('' gpmt reg '',2d15.8)') ff1*rs
        ggam1=ff1/gg1 
        aa2=p*((ggam1-xl/rs)*fn(l)+p*fn(l+1))/
     >        ((ggam1-xl/rs)*fb(l)+p*fb(l+1))   
      else
        gg1=g11(ns)/rs
        ff1=f11(ns)/rs
c       write(6,'('' gmt reg  '',2d15.8)') gg1*rs
c       write(6,'('' fmt reg  '',2d15.8)') ff1*rs
        ggam1=ff1/gg1 
        aa2=p*(ggam1*fn(l)-smlm1*p*fn(l+1)/c)/
     >        (ggam1*fb(l)-smlm1*p*fb(l+1)/c)   
      end if
c
c This is now the inverse reactance!
      tminv(2,2)=-aa2
c This is now the inverse t-matrix!
      tminv(2,2)=tminv(2,2)+sqrtm1*p
c
!     write(6,*) sname,2,3
      if(iflag.eq.0) return
!     write(6,*) sname,2,4
      call dirmagi1op(socsc,ce,l,my,vr,br,bopr,dx,xnot,rs,ns,
     >              g12,f12,gp12,iprpts)
!     write(6,*) sname,2,5
c     call dirmagi1(socsc,ce,l,my,vr,br,dx,xnot,rs,ns,g12,f12,gp12)
!     write(6,*) sname,2,6
c
      if(scale) then
        gg2=g12(ns)/rs
        ff2=gp12/rs
c       write(6,'('' gmt irr  '',2d15.8)') gg2*rs
c       write(6,'('' gpmt irr  '',2d15.8)') ff2*rs
        ggam2=ff2/gg2 
        ggamb=(xl*fb(l)/rs-p*fb(l+1))/fb(l)
        aa1=p*((ggamb-xl/rs)*fn(l)+p*fn(l+1))/(ggamb*gg1-ff1)   
        bb1=  ((ggam2-xl/rs)*fb(l)+p*fb(l+1))/(ggam2*gg1-ff1)   
        bb2=  ((ggam1-xl/rs)*fb(l)+p*fb(l+1))/(ggam1*gg2-ff2)   
        do i=1,ns
          gz(i,2,2)=g11(i)*aa1
          gj(i,2,2)=g11(i)*bb1+g12(i)*bb2
        end do
      else
        gg2=g12(ns)/rs
        ff2=f12(ns)/rs
!       write(6,'('' gmt irr  '',2d15.8)') gg2*rs
!       write(6,'('' fmt irr  '',2d15.8)') ff2*rs*c
        ggam2=ff2/gg2 
        ggamb=smlm1*p*fb(l+1)/(c*fb(l))
        aa1=p*(ggamb*fn(l)-smlm1*p*fn(l+1)/c)/(ggamb*gg1-ff1)   
        bb1=  (ggam2*fb(l)-smlm1*p*fb(l+1)/c)/(ggam2*gg1-ff1)   
        bb2=  (ggam1*fb(l)-smlm1*p*fb(l+1)/c)/(ggam1*gg2-ff2)   
        do i=1,ns
!         write(6,*) '>',i
          gz(i,2,2)=g11(i)*aa1
          fz(i,2,2)=f11(i)*aa1
          gj(i,2,2)=g11(i)*bb1+g12(i)*bb2
          fj(i,2,2)=f11(i)*bb1+f12(i)*bb2
        end do
      end if

!     write(6,*)sname,2,7
c
      return 
c
!   1 continue
      end if
!     write(6,*) sname,1,1
      call dirmago2op(socsc,ce,l,my,vr,br,bopr,dx,xnot,rs,ns,
     >                g1,f1,gp1,iprpts)
!     write(6,*) sname,1,2
c   1 call dirmago2(socsc,ce,l,my,vr,br,dx,xnot,rs,ns,g1,f1,gp1)
c
      jmat(1,1)=fb(l)
      jmat(2,2)=fb(l)
      nmat(1,1)=p*fn(l)
      nmat(2,2)=p*fn(l)
c
      if(scale) then
        jbmat(1,1)=xl*fb(l)/rs-p*fb(l+1)
        jbmat(2,2)=jbmat(1,1)
        nbmat(1,1)=xl*p*fn(l)/rs-p*p*fn(l+1)
        nbmat(2,2)=nbmat(1,1)
        do ii=1,2
        do jj=1,2
          g1mat(ii,jj)=g1(ii,jj,ns)/rs 
          f1mat(ii,jj)=gp1(ii,jj)/rs 
c         write(6,'('' gmt reg  '',2d15.8)') g1(ii,jj,ns)
c         write(6,'('' gpmt reg  '',2d15.8)') gp1(ii,jj) 
        end do
        end do
      else
        jbmat(1,1)=sl*p*fb(l-1)/c
        jbmat(2,2)=smlm1*p*fb(l+1)/c
        nbmat(1,1)=sl*p*p*fn(l-1)/c
        nbmat(2,2)=smlm1*p*p*fn(l+1)/c
        do ii=1,2
        do jj=1,2
          g1mat(ii,jj)=g1(ii,jj,ns)/rs 
          f1mat(ii,jj)=f1(ii,jj,ns)/rs 
c         write(6,'('' gmt reg  '',2d15.8)') g1(ii,jj,ns)
c         write(6,'('' fmt reg  '',2d15.8)') f1(ii,jj,ns)
        end do
        end do
      end if
c
!     write(6,*) sname,1,3
      call repl(x1,g1mat,2,2)
      call gjinv(x1,2,2,detl)
      call repl(gam1,f1mat,2,2)
      call doubmt(gam1,x1,2,2)
c
      call repl(x1,gam1,2,2) 
      call doubmt(x1,jmat,2,2) 
      call submat(x1,jbmat,2,2) 
      call repl(a2,x1,2,2) 
      call gjinv(a2,2,2,detl)
      call repl(x2,gam1,2,2) 
      call doubmt(x2,nmat,2,2) 
      call submat(x2,nbmat,2,2) 
      call doubmt(a2,x2,2,2)
!     write(6,*) sname,1,4 
c     
c This is now the inverse reactance matrix!
      do ii=1,2                 
      do jj=1,2                 
        tminv(ii,jj)=-a2(ii,jj)
      end do
      end do
c This is now the inverse t-matrix!
      tminv(1,1)=tminv(1,1)+sqrtm1*p
      tminv(2,2)=tminv(2,2)+sqrtm1*p
c
!     write(6,*) sname,1,5
      if(iflag.eq.0) return
!     write(6,*) sname,1,6
      call dirmagi2op(socsc,ce,l,my,vr,br,bopr,dx,xnot,rs,ns,
     >                g2,f2,gp2,iprpts)
!     write(6,*) sname,1,7
c     call dirmagi2(socsc,ce,l,my,vr,br,dx,xnot,rs,ns,g2,f2,gp2)
c
      if(scale) then
        do ii=1,2
        do jj=1,2
          g2mat(ii,jj)=g2(ii,jj,ns)/rs 
          f2mat(ii,jj)=gp2(ii,jj)/rs 
        end do
        end do
      else
        do ii=1,2
        do jj=1,2
          g2mat(ii,jj)=g2(ii,jj,ns)/rs 
          f2mat(ii,jj)=f2(ii,jj,ns)/rs 
        end do
        end do
      end if
c
      call repl(x2,g2mat,2,2)
      call gjinv(x2,2,2,detl)
      call repl(gam2,f2mat,2,2)
      call doubmt(gam2,x2,2,2)
      gamb(1,1)=jbmat(1,1)/jmat(1,1)
      gamb(2,2)=jbmat(2,2)/jmat(2,2)
c
      call repl(b2,gam1,2,2) 
      call doubmt(b2,g2mat,2,2) 
      call submat(b2,f2mat,2,2) 
      call gjinv(b2,2,2,detl)
      call doubmt(b2,x1,2,2) 
c
      call repl(a1,gamb,2,2)
      call doubmt(a1,g1mat,2,2)
      call submat(a1,f1mat,2,2) 
      call gjinv(a1,2,2,detl)
      call repl(x1,gamb,2,2)
      call doubmt(x1,nmat,2,2)
      call submat(x1,nbmat,2,2) 
      call doubmt(a1,x1,2,2) 
c
      call repl(b1,gam2,2,2)
      call doubmt(b1,g1mat,2,2)
      call submat(b1,f1mat,2,2) 
      call gjinv(b1,2,2,detl)
      call repl(x1,gam2,2,2)
      call doubmt(x1,jmat,2,2)
      call submat(x1,jbmat,2,2) 
      call doubmt(b1,x1,2,2)
!     write(6,*) sname,1,8 
c
      do i=1,ns
        call repl(x1,g1(1,1,i),2,2)
        call repl(x2,f1(1,1,i),2,2)
        call doubmt(x1,a1,2,2)
        call doubmt(x2,a1,2,2)
        do ii=1,2
        do jj=1,2
          gz(i,ii,jj)=x1(ii,jj)
          fz(i,ii,jj)=x2(ii,jj)
        end do
        end do
        call repl(x1,g1(1,1,i),2,2)
        call repl(x2,f1(1,1,i),2,2)
        call repl(x3,g2(1,1,i),2,2)
        call repl(x4,f2(1,1,i),2,2)
        call doubmt(x1,b1,2,2)
        call doubmt(x2,b1,2,2)
        call doubmt(x3,b2,2,2)
        call doubmt(x4,b2,2,2)
        call addmat(x1,x3,2,2)
        call addmat(x2,x4,2,2)
        do ii=1,2
        do jj=1,2
          gj(i,ii,jj)=x1(ii,jj)
          fj(i,ii,jj)=x2(ii,jj)
        end do
        end do
      end do
c
      if(scale) then
        call zeroout(fz,8*iprpts)
        call zeroout(fj,8*iprpts)
      endif
c
      return
      end
