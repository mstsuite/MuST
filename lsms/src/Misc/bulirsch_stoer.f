!           call bsstep (y,dy,nv,nr,rfrom,htry,watol,yp,hdid,hnext,idid,
!     &      rad(jp),einside,pot(1,jp),
!     &      vso,1,isymm,eunit,dfv,b,allp1)

      subroutine bulirsch_stoer(y,dy,nv,x0,dx,
     &                          rn,ecomp,vr,eunit,b,allp1,nextrp)
      implicit none

      integer nv,nextrp
      real*8 y(nv),dy(nv),x0,dx

      real*8 rn,vr(*),eunit,b,allp1
      complex*16 ecomp

      integer max_steps
!     parameter (max_steps=8)
      parameter (max_steps=12)
      real*8 Tik(nv,max_steps),xi(max_steps),df(nv)
      real*8 Qik(nv,max_steps),Dik(nv,max_steps)
      real*8 y_out(nv)
      real*8 h,hh,err,den

      integer i,ni,j,nj,ii
      PARAMETER (nj=12)

      hh=dx/nj
      do j=1,nj
      do i=1,max_steps
        ni=2*i
        h=dx/(ni*nj)
        call mod_midpoint(y,dy,nv,x0+(j-1)*hh,hh,ni,y_out,
     &                    rn,ecomp,vr,eunit,b,allp1,nextrp)
        call extrapolateT0(h*h,y_out,Tik,xi,df,nv,max_steps,i)
!       call extrapolateQD0(h*h,y_out,Tik,Qik,Dik,xi,df,nv,max_steps,i)

!        err=0.0d0
!        do ii=1,nv,2
!          den=abs(yscal(ii))+abs(yscal(ii+1))
!          if(den.gt.1.0e-10) then 
!            err=err+(abs(df(ii))+abs(df(ii+1)))/den
!          endif
!        enddo

      end do
      do i=1,nv
        y(i)=Tik(i,1)
      end do
      call dfv_m(x0+dx,y,dy,nv,nextrp,rn,ecomp,vr,eunit,b,allp1)
      end do

      end

      subroutine extrapolateQD0(x,f,fout,Qik,Dik,xi,df,n_val,imax,i)
      implicit none

! calculate the diagonal ik for a fixed i, k=1..k of the interpolation Tableau
! store Qik as Qik(*,i-k+1)
! the new extrapolation to h=0 will be in fout

      real*8 fout(n_val),Qik(n_val,imax),Dik(n_val,imax)
      real*8 xi(imax),f(n_val),c(n_val)
      real*8 df(n_val)
      real*8 d,x,delta
      integer k,j,n_val,imax,i

      xi(i)=x
! set Ti1
      do j=1,n_val
        Qik(j,i)=f(j)
        Dik(j,i)=f(j)
        c(j)=0.0
!        df(j)=Tik(j,1)
      end do
! calculate remaining entries in Qik
      do k=i-1,1,-1
        d=1.0/(xi(k)-xi(i))
        do j=1,n_val
          delta=Dik(j,k+1)-Qik(j,k)
          Qik(j,k)=delta*xi(i)*d
          Dik(j,k)=delta*xi(k)*d
          c(j)=c(j)+Qik(j,k)
        end do
      end do
      do j=1,n_val
        fout(j)=c(j)+f(j)
! is this the right error expression?
        df(j)=Qik(j,1)
      end do
      end

      subroutine extrapolateT0(x,f,Tik,xi,df,n_val,imax,i)
      implicit none

! calculate the diagonal ik for a fixed i, k=1..k of the interpolation Tableau
! store Tik as Tik(*,i-k+1)
! the new extrapolation to h=0 will be in Tik(*,1)

      real*8 Tik(n_val,imax)
      real*8 xi(imax),f(n_val)
      real*8 df(n_val)
      real*8 d,x
      integer k,j,n_val,imax,i

      xi(i)=x
! set Ti1
      do j=1,n_val
        Tik(j,i)=f(j)
        df(j)=Tik(j,1)
      end do
! calculate remaining entries in Tik
      do k=i-1,1,-1
        d=1.0/(xi(k)/xi(i) - 1.0)
        do j=1,n_val
          Tik(j,k)=Tik(j,k+1)+(Tik(j,k+1)-Tik(j,k))*d
        end do
      end do
      do j=1,n_val
        df(j)=2.0d0*Tik(j,2)-df(j)-Tik(j,1)
      end do
      end
