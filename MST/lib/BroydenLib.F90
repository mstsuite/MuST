!***


!****q* math/giterat
!
! NAME
!   giterat
!
! COPYRIGHT
! 
! DESCRIPTION
!
! REFERENCES
!
! INPUTS
!
! OUTPUT
!
! NOTES
!
! WARNING
!
! PARENTS
!
! CHILDREN
!   dglei
! 
! SOURCE
!
!---------------------------------------------------------------------


!..............................................................................
!  09.03.95  ******************************************************************
!******************************************************************************
!------------------------------------------------------------------------------
subroutine giterat(niter,mdim,pmix,x,n,t,u,m,k,p,ppq,pq,w0,tp,tolerance)
!------------------------------------------------------------------------------
  use KindParamModule, only : IntKind, RealKind
  implicit none
!io............................................................................
  integer (kind=IntKind) :: niter,mdim,n,m,k
  real (kind=RealKind) :: x(:),t(:),pmix,u,p,ppq,pq,w0,tp,tolerance
!local.........................................................................
  integer (kind=IntKind) :: i,m1,mm1,nu,nu2,nn,na,nb,i1,i2,ju,j,ju2,indef,iu
  real (kind=RealKind) :: err,tm,q,s,uu,ux
!******************************************************************************
!******************************************************************************
! iterative solution of the n-dimensional non-linear equation system
!
!    x = vector(x,n,func)
!
! by searching in m-dimensional subspaces (m.le.mdim: mdim.ge.1).
!
! convergence for u**2 = (x - vector(x,n))**2 .le. tolerance.
! permitted relative range of linear extrapolation: 1/w0.
! iteration steps dx = p*u from the interpolated minimum of u**2
! in the actual subspace are used (p self-controlled).
! the maximum number of iterations is kmax.
! if(print) the iteration is protocolled.
! 4integers and
! 1+(2*n+2*mdim)*(mdim+1) real elements of common/terat/ are used.
!
! storage in /terat/:
! common/terat/k,m,m1,m1d,u,
!              ((x (i),i=1,n),j=0,mdim),((u (i),i=1,n),j=0,mdim),
!                 j                        j
!              1                     nu                      nu2
!
!*             ((u *u ,i=0,mdim),j=0,mdim),
!                 i  j
!               (i.ge.j)                na
!
!*             ((a(i,j),i=1,mdim-1),j=1,mdim-1),
!
!                                            nb
!
!*             (b(i),i=1,mdim-1),(mu(i),i=1,mdim-1)
!
!                             nn
!
! (mu solution of the linear equation system a*mu=b for the
! minimum of u**2).
!
!***********************************************************************
!
! storage x ==> x  and calculation of u  and of u *u , i=1,...,m
!                m                     m         i  m
!..............................................................................
  mm1=mdim+1
  nu=n*mm1
  nu2=2*nu
  na=nu2+mm1**2
  nb=na+(mdim-1)**2
  nn=nb+mdim-1

  m1=m+1
  k=k+1
  if(k.gt.niter) return

  i1=n*m
  i2=nu+i1
  err=0.d0
  do i=1,n
     t(i2+i)=x(i)-t(i1+i)
     err=err+(x(i)-t(i1+i))**2
  enddo
  ju2=nu2+m1-mm1
  do j=1,m1
     ju2=ju2+mm1
     i1=nu+m*n
     i2=nu+(j-1)*n
     tp=0.d0
     do i=1,n
        tp=tp+t(i1+i)*t(i2+i)
     enddo
     t(ju2)=tp
  enddo
  tp=err
!..............................................................................
! convergence: restoration of x and return
!..............................................................................
  if(sqrt(tp).lt.tolerance) then
     do i=1,n
        x(i)=t(n*m+i)
     enddo
     return
!..............................................................................
! case m=0 simple mixing.......................................................
!..............................................................................
  elseif(mdim.eq.1.or.m.eq.0) then
     do i=1,n
        x(i)=t(i)+p*t(nu+i)
     enddo
     u=t(nu2+1)
     m=1
     if(mdim.eq.1) m=0
     do i=1,n
        t(n*m+i)=x(i)
     enddo
     ppq=max(ppq,w0)
     return
  endif
!
! exclusionn of useless points
!
  if(m.eq.1) goto 32
  ju=0
  if(m.eq.mdim) ju=1
  if(u/tp.lt.w0.or.t(nu2+1)/tp.lt.1.d0) ju=m/2
  pq=min(1.d0,sqrt(u/tp))
  if(ju.eq.0) go to 32
38 i1=(m-ju+1)*n
  i2=ju*n
  do i=1,i1
     t(i)=t(i+i2)
     t(nu+i)=t(nu+i+i2)
  enddo
  i1=(m-ju+1)*mm1
  i2=ju*(mm1+1)
  do i=1,i1
     t(nu2+i)=t(nu2+i+i2)
  enddo
  m=m-ju
  if(m.gt.0) go to 32
  
  m=1
  do i=1,n
     t(n*m+i)=x(i)
  enddo
  ppq=max(ppq,w0)
  return
!
! rearranging of x , u , u *u  with descending order of u **2;
!                 i   i   i  j                           i
!
32 do ju=1,m
     j=m-ju
     ju2=nu2+mm1*j+j+1
     uu=t(ju2)
     ux=t(ju2+mm1+1)
     if(ux.le.uu) exit
     i1=j*n
     i2=i1+n
     do i=1,n
        i1=i1+1
        i2=i2+1
        s=t(i1)
        t(i1)=t(i2)
        t(i2)=s
        s=t(nu+i1)
        t(nu+i1)=t(nu+i2)
        t(nu+i2)=s
     enddo
     t(ju2)=ux
     t(ju2+mm1+1)=uu
     if(j.ne.0) then
        i1=nu2+j+1
        i2=i1+1
        do i=1,j
           s=t(i1)
           t(i1)=t(i2)
           t(i2)=s
           i1=i1+mm1
           i2=i2+mm1
        enddo
     endif
     if(ju.eq.1) cycle
     i1=ju2+1
     j=j+3
     do i=j,m1
        i1=i1+1
        i2=i1+mm1
        s=t(i1)
        t(i1)=t(i2)
        t(i2)=s
     enddo
  enddo
!
! finding the minimum of u**2 for u(mu) = u  + sum mu *(u -u ).
!                                          m   i=0   i   i  m
  u=t(nu2+mm1*m+m+1)
  i1=nu2-mm1
  iu=i1+m+1
  do i=1,m
     i1=i1+mm1
     iu=iu+mm1
     t(nb+i)=1.d0-t(iu)/u
     i2=iu
     ju=na+i+m*(i-1)
     do j=i,m
        t(ju)=(u+t(i1+j)-t(iu)-t(i2))/u
        i2=i2+mm1
        ju=ju+m
     enddo
  enddo
  if(m.gt.1) go to 48
  uu=t(na+1)
  if(uu.eq.0.d0) go to 42
  t(nn+1)=t(nb+1)/uu
  goto 44
48 call dglei(t(na+1),t(nb+1),t(nn+1),m,m,indef)
  if(indef.eq.0) goto 44
!
!     creation of a stochastic step, if no minimum of u**2 is found
!
42 do i=1,n
      u=t(nu+n*m+i)
      x(i)=t(n*m+i)+ppq*u*sin(2.d0*i*k)
   enddo
   u=1.d29
   m=m+1
!..............................................................................
   do i=1,n
      t(n*m+i)=x(i)
   enddo
   ppq=max(ppq,w0)
   return
!..............................................................................
!
!                  m-1
! setting x = x  + sum mu *(x -x ) + ppq*u   ,
!              m   i=0   i   i  m         min
!
!               m-1
! where p = abs(sum mu *(x -x ) )/abs(u ).
!               i=0   i   i  m         m
!
! and ppq is the mean of this p and the previous one,
! eventually reduced by a factor pp.
!..............................................................................
44 q=p
   p=0.d0
   iu=n*m
   ju=nu+iu
   do i=1,n
      ux=t(iu+i)
      i1=i
      s=0.d0
      do j=1,m
         s=s+t(nn+j)*(t(i1)-ux)
         i1=i1+n
      enddo
      p=p+s**2
   enddo
   p=(sqrt(p/u)+q)/2.d0
   p=min(pmix/pq,p)
   ppq=p*pq
   s=1.d0
   do i=1,m
      s=s-t(nn+i)
   enddo
   u=0.d0
   do i=1,n
      ux=s*t(iu+i)
      uu=s*t(ju+i)
      i1=i
      do j=1,m
         tm=t(nn+j)
         ux=ux+tm*t(i1)
         uu=uu+tm*t(nu+i1)
         i1=i1+n
      enddo
      u=u+uu**2
      x(i)=ux+ppq*uu
   enddo
   ju=m
   if(u/tp.lt.w0) go to 38
   m=m+1
!==============================================================================
  do i=1,n
     t(n*m+i)=x(i)
  enddo
  ppq=max(ppq,w0)
  return

!101 format(/'dimension',i5,12x,'last deviation u=',1pd14.7,          &
!         &     /'new pmix=',0pf12.6,5x,                              &
!         &     'new deviation  u=',1pd14.7)

end subroutine giterat

!***


!****q* math/giterbr
!
! NAME
!   giterbr
!
! COPYRIGHT
! 
! DESCRIPTION
!
! REFERENCES
!
! INPUTS
!
! OUTPUT
!
! NOTES
!
! WARNING
!
! PARENTS
!
! CHILDREN
!   gelim, subs
! 
! SOURCE
!-------------------------------------------------------------------------
subroutine giterbr(pmix,w0,mdim,x,nx,m,nm,nml,fm,fml,delta,bkni,tp)
  use KindParamModule, only : IntKind, RealKind
  implicit none
!i/o...........................................................................
  integer (kind=IntKind) :: nx,mdim,m
  real (kind=RealKind) :: x(nx),pmix,w0,nm(nx),nml(nx),fm(nx),fml(nx),&
       & delta(nx,mdim,2),bkni(mdim,mdim),tp
!local.........................................................................
  integer (kind=IntKind) :: i,j,k,n
  real (kind=RealKind) :: fnorm,s,sm
  integer (kind=IntKind), allocatable :: isc(:)
  real (kind=RealKind), allocatable :: df(:),deltan(:),deltaf(:),w(:),bkn(:,:),&
     & fkm(:),gm(:),bk(:,:)
!simple mixing................................................................
  m=m+1
  if(mdim.eq.1.or.m.eq.1) then
     nml=nm
     fml=x-nm

     fm=fml

     x(1:nx)=x(1:nx)*pmix+(1.d0-pmix)*nml(1:nx)
     nm(1:nx)=x(1:nx)
     
     tp=sum(fm*fm)
     return
  endif

  fm=x-nm
  tp=sum(fm*fm)

  allocate(isc(mdim),df(nx),deltan(nx),deltaf(nx),w(mdim),bkn(mdim,mdim),&
       & fkm(mdim),gm(mdim),bk(mdim,mdim))
  isc=0
  df=0.d0
  w=0.d0
  bkn=0.d0
  fkm=0.d0
  gm=0.d0
  m=min(m,mdim)
!============================================================================
  s=0.d0
  do i=1,nx
     s=s+fm(i)*fm(i)
  enddo
  w(m-1)=1.d0/sqrt(s)
  do i=1,nx                                 
     df(i)=fm(i)-fml(i)
  enddo
  fnorm=0.d0
  do i=1,nx
     fnorm=fnorm+df(i)*df(i)
  enddo
  fnorm = 1.d0/sqrt(fnorm)                
  do i=1,nx                             
     deltan(i)=fnorm*(nm(i)-nml(i) )  
     deltaf(i)=fnorm*df(i)
  enddo
!------------------------------------------------------------------------------
!     inverse beta(k,n) (13a)
!------------------------------------------------------------------------------
  bkni(m-1,m-1)=w0*w0+w(m-1)*w(m-1)   
  do j=1, m-2                    
     w(j)=sqrt(bkni(j,j)-w0*w0)
  enddo
!==============================================================================
  s=0.d0
  do i=1,nx
     s=s+deltaf(i)*fm(i)
  enddo
  fkm(m-1)=s
!===========================================================================
  do n=1,m-2                     
     df(1:nx)=delta(1:nx,n,1)
     s=0.d0
     sm=0.d0
     do i=1,nx
        s=s+df(i)*fm(i)
        sm=sm+df(i)*deltaf(i)
     enddo

     fkm(n)=s
     bkni(n,m-1)=sm*w(m-1)*w(n) 
  enddo

  delta(1:nx,m-1,1)=deltaf(1:nx)

  nml(1:nx)=nm(1:nx)
  fml(1:nx)=fm(1:nx)

  do i=1,m-2            
     do j =i+1,m-1       
        bkni(j,i) =bkni(i,j)  
     enddo
  enddo
!==============================================================================
  bk=bkni
  if (m.eq.2) then  
     bkn(1,1)=1.0d0/bkni(1,1)
  else                     
     do i=1,m-1
        do j=1,m-1
           bkn(j,i)=0.0d0
        enddo
        bkn(i,i)=1.0d0
     enddo
     call gelim(bk,isc,mdim,m-1,1.d-14)
     do i=1,m-1
        call subs(bk,isc,bkn(1,i),mdim,m-1)
     enddo
  endif
!==============================================================================
  do i=1,m-2                                
     do j=i+1,m-1                           
        bkn(i,j)=0.5d0*(bkn(i,j)+bkn(j,i)) 
        bkn(j,i)=bkn(i,j)                
     enddo
  enddo
!==============================================================================
!     gamma(m,nx)
!===========================================================================
  do n=1,m-1                         
     s=0.0d0                          
     do k=1,m-1                      
        s=s+w(k)*bkn(n,k)*fkm(k)
     enddo
     gm(n)=s
  enddo
!===========================================================================
!     (15a)
!===========================================================================
  do i=1,nx
     deltan(i)=deltan(i)+deltaf(i)*pmix
  enddo
  do i=1,nx
     deltaf(i)=0.d0
  enddo

  do n=1,m-2                        
     df(1:nx)=delta(1:nx,n,2)
     do i=1,nx
        df(i)=df(i)*w(n)*gm(n)
        deltaf(i)=deltaf(i)+df(i)
     enddo
  enddo

  delta(1:nx,m-1,2)=deltan(1:nx)    

  do i=1,nx
     deltaf(i)=deltaf(i)+deltan(i)*w(m-1)*gm(m-1)
  enddo

  do i=1,nx
     nm(i)=nm(i)+fm(i)*pmix
  enddo

  do i=1,nx
     nm(i)=nm(i)-deltaf(i)
  enddo
  deallocate(isc,df,deltan,deltaf,w,bkn,fkm,gm,bk)
  
  x=nm
  return
end subroutine giterbr

!***

!****q* math/dglei
!
! NAME
!   dglei
!
! COPYRIGHT
! 
! DESCRIPTION
!
! REFERENCES
!
! INPUTS
!
! OUTPUT
!
! NOTES
!
! WARNING
!
! PARENTS
!
! CHILDREN
!   dchol
! 
! SOURCE
!-----------------------------------------------------------------------------



!***********************************************************************
subroutine dglei (a,b,x,n,m,indef)
!***********************************************************************
!***********************************************************************
!     solution of a linear inhomogenous equation system with a
!     symmetric positive definite real matrix.
!     double precision algebra.
!
!     input:
!       a     : matrix of the problem
!       b     : inhomogenity vector
!       n     : order of the problem
!       m     : specified dimension of a(m,m)
!
!     output:
!       x     : solution vector
!       indef : =0 if a is positive definite, =1 otherwise (no solution)
!***********************************************************************
  use KindParamModule, only : IntKind, RealKind
  implicit none
  integer (kind=IntKind) :: m,i,n,indef,k,k1,i1
  real (kind=RealKind) ::  a(m,m),b(m),x(m),s
!***********************************************************************
  call dchol(a,n,m,indef)
  if(indef.eq.1) return
  do i=1,n
     x(i)=b(i)
  enddo !i
  s=x(1)/a(1,1)
  x(1)=s
  do i=2,n
     i1=i-1
     do k=i,n
        x(k)=x(k)-s*a(i1,k)
     enddo !k
     s=x(i)/a(i,i)
     x(i)=s
  enddo !i
  s=x(n)/a(n,n)
  x(n)=s
  do i1=2,n
     i=n+2-i1
     do k1=i1,n
        k=n+1-k1
        x(k)=x(k)-a(k,i)*s
     enddo !k1
     i=i-1
     s=x(i)/a(i,i)
     x(i)=s
  enddo !i1
  return
end subroutine dglei

!***


!****q* math/dchol
!
! NAME
!   dchol
!
! COPYRIGHT
! 
! DESCRIPTION
!
! REFERENCES
!
! INPUTS
!
! OUTPUT
!
! NOTES
!
! WARNING
!
! PARENTS
!
! CHILDREN
! 
! SOURCE
!---------------------------------------------------------------------



!  10.09.93  ***********************************************************
!***********************************************************************
subroutine dchol (a,n,m,indef)
!***********************************************************************
!***********************************************************************
!     cholesky procedure for the reduction of a positive definite
!     symmetric real matrix into a trigonal matrix.
!     double precision algebra.
!
!     input:
!       a     : symmetric matrix
!       n     : order of a
!       m     : specified dimension of a(m,m)
!
!     output:
!       a     : upper trigonal matrix
!       indef : =0 if a is positive definite, =1 otherwise
!***********************************************************************
  use KindParamModule, only : IntKind, RealKind
  implicit none
  integer (kind=IntKind) :: n,m,ip,indef,lunit,k,kp,i
  real (kind=RealKind) :: a(m,m),ap
!***********************************************************************
  indef=0
  do ip=1,n
     ap=a(ip,ip)
     if (ap.le.0.d0) then
        indef=1
!       write(lunit(),*)' dchol: ',n,ip,ap
        write(6,*)' dchol: ',n,ip,ap
        return
     endif
     ap=sqrt(ap)
     a(ip,ip)=ap
     if (ip.eq.n) return
     kp=ip+1
     do k=kp,n
        a(ip,k)=a(ip,k)/ap
     enddo
     do i=kp,n
        do k=i,n
           a(i,k)=a(i,k)-a(ip,i)*a(ip,k)
        enddo
     enddo
  enddo
  return
end subroutine dchol

!***

!==============================================================================
subroutine gelim(ar,nt,np,n,emach)
  use KindParamModule, only : IntKind, RealKind
  implicit none
  integer (kind=IntKind) :: n
  integer (kind=IntKind) :: nt(n),np,ii,i,in,j,k
  real (kind=RealKind) :: ar(np,n),yrr,dum,emach

  if(n.lt.2) goto 15
  do ii=2,n
     i=ii-1
     yrr=ar(i,i)
     in=i
     do j=ii,n
        if(abs(yrr).lt.abs(ar(j,i))) then 
           yrr=ar(j,i)
           in=j
        endif
     enddo
     nt(i)=in
     if(in.lt.i.or.in.gt.i) then
        do j=i,n
           dum=ar(i,j)
           ar(i,j)=ar(in,j)
           ar(in,j)=dum
        enddo
     endif
     if(abs(yrr).lt.emach) then
        ar(i,i)=emach*emach
        goto 12
     endif
     do j=ii,n
        if(abs(ar(j,i)).gt.emach) then
           ar(j,i)=ar(j,i)/yrr
           do k=ii,n
              ar(j,k)=ar(j,k)-ar(i,k)*ar(j,i)
           enddo
        endif
     enddo
12   continue
  enddo
15 if(abs(ar(n,n)).le.emach) then
     ar(n,n)=emach*emach
  endif
  return
end subroutine gelim
!==============================================================================
subroutine subs(ar,nt,xr,np,n)
  use KindParamModule, only : IntKind, RealKind
  implicit none
  integer (kind=IntKind) :: n
  integer (kind=IntKind) :: nt(n),np,i,ii,in,ij,j
  real (kind=RealKind) :: ar(np,n),xr(n),dum
  if(n.lt.2) goto 18
  do ii=2,n
     i=ii-1
     if(nt(i).lt.i) then
        goto 16
     elseif(nt(i).eq.i) then
        goto 17
     else
        goto 16
     endif
16   in=nt(i)
     dum=xr(in)
     xr(in)=xr(i)
     xr(i)=dum
17   do j=ii,n
        xr(j)=xr(j)-ar(j,i)*xr(i)
     enddo
  enddo
!  back substitution
18 do ii=1,n
     i=n-ii+1
     ij=i+1
     if(i.lt.n) then
        goto 21
     elseif(i.eq.n) then
        goto 25      
     else
        goto 21
     endif
21   do j=ij,n
        xr(i)=xr(i)-ar(i,j)*xr(j)
     enddo
25   xr(i)=xr(i)/ar(i,i)
  enddo
  return
end subroutine subs

!***
