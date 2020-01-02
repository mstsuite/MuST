c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makegij(lmaxi,kkri,lmaxj,kkrj,
     >                   lmax,kkrsz,ndlj,ndlm,
     >                   prel,rij,sinmp,cosmp,
     >                   clm,plm,cgnt,lmax_cg,lofk,mofk,
     >                   ilp1,illp,
     >                   hfn,dlm,gij,
     >                   pi4,iprint,istop)
c     ==============================================================
c
      implicit   none
c
      character  sname*32
      character  istop*32
c
      integer    lmaxi
      integer    lmaxj
      integer    lend
      integer    kkri
      integer    kkrj
      integer    lmax,lmax_cg
      integer    kkrsz
      integer    ndlj
      integer    ndlm,ndlm_local
      integer    lofk(ndlj)
      integer    mofk(ndlj)
      integer    l,ll
      integer    l1
      integer    l2
      integer    m
      integer    ma
      integer    m1
      integer    m2
      integer    lm1
      integer    lm2
      integer    l3,m3,llow
      integer    j
      integer    iprint
c
      real*8     rij(3)
      real*8     sinmp(0:lmaxi+lmaxj)
      real*8     cosmp(0:lmaxi+lmaxj)
      real*8     clm(ndlm)
      real*8     plm(ndlm)
!     original definition of cgnt in LSMS_1
!     real*8     cgnt(lmax+1,kkrsz,kkrsz)
      real*8     cgnt(lmax_cg+1,(lmax_cg+1)**2,(lmax_cg+1)**2)
      real*8     pi4
      real*8     rmag
      real*8     pmag
      real*8     costheta
      real*8     m1m
      real*8     ptol
      real*8     zero
      real*8     one
c
      complex*16 prel
      complex*16 ilp1(0:lmaxi+lmaxj)
!     original definition of illp in LSMS_1
!     complex*16 illp(kkrsz,kkrsz)
      complex*16 illp((lmax+1)**2,kkrsz)
      complex*16 hfn(0:lmaxi+lmaxj)
      complex*16 dlm(ndlj)
      complex*16 gij(kkri,kkrj)
      complex*16 z
      complex*16 fac
      complex*16 sqrtm1
      complex*16 czero
      complex*16 cone
c
      parameter (sqrtm1=(0.0d0,1.0d0))
      parameter (czero =(0.0d0,0.0d0))
      parameter (cone  =(1.0d0,0.0d0))
      parameter (ptol  =1.0d-06)
      parameter (zero  =0.0d0)
      parameter (one   =1.0d0)
      parameter (sname ='makegij')
c
c     *****************************************************************
c    
c
c      ij         l+1                                m ->  *
c     D  [E]  = -i    * sqrt(E) * h [sqrt(E)*R  ] * Y [R  ]
c      L                           l          ij     l  ij
c    
c
c                 -l+1                               m  ->  *
c             = -i    * sqrt(E) * h [sqrt(E)*R  ] * Y [-R  ]
c                                  l          ij     l   ij
c    
c The extra factor (-1)^m removed by xgz July 97 from plm (plglmax_c.f)
c to make Ylm consistent with Numerical Recipes.
c
c      m       [ (2*l+1)*(l-|m|)!]  m
c     Y  = sqrt[-----------------]*P [cos(theta)]*exp(i*m*phi)
c      l       [   4*pi*(l+|m|)! ]  l
c
c     for m>=0 
c
c      m                  m  -m           *
c     Y [theta,phi] = (-1)  Y  [theta,phi]      for m<0
c      l                     l
c
c     ->    ->   ->
c     R   = R  - R  ==> [theta,phi]
c      ij    j    i
c
c     cos(theta)=Rij(3)/sqrt(Rij(1)**2+Rij(2)**2+Rij(3)**2)
c
c     cos(phi)  =Rij(1)/sqrt(Rij(1)**2+Rij(2)**2)
c
c                m     m
c     Indexing: P  :  P(l*(l+1)/2+m+1)  only +m are calculated
c                l     l
c
c                m     m
c     Indexing: C  :  C(l*(l+1)/2+m+1)  only +m are calculated
c                l     l
c
c                m     m
c     Indexing: D  :  D(l*(l+1)+m+1)    all   m are calculated
c                l     l
c                    
c     Now for the real space structure contant we have :
c                    
c      ij             l-l'         L       ij
c     G   (E) = 4*pi*i    * SUM { C     * D  (E) }
c      L,L'                  L"    L',L"   L"
c
c     *****************************************************************
c
c     =================================================================
      if(lmax.lt.0) then
         write(6,'('' makegij:: bad arguments: lmax='',i5)') lmax
         call fstop(sname)
      endif
      lend=lmaxi+lmaxj
c     =================================================================
c     calculate the hankel function.[dangerous code if z is close to 0]
c     hankel function is hfn*exp(i*z)/z
      rmag=sqrt(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
      if(rmag.lt.ptol) stop 'Error in makegij: rmag = 0.'
      if(prel.eq.czero) then
c if prel is zero then we calculate Gij for multipole Coulomb interaction.
	hfn(0)=cone/rmag
	do l=1,lend
	  hfn(l)=sqrtm1*(2*l-1)*hfn(l-1)/rmag
	enddo
      else
      z=prel*rmag
      hfn(0)=-sqrtm1
      hfn(1)=-cone-sqrtm1/z
      do l=1,lend-1
         hfn(l+1)=(2*l+1)*hfn(l)/z - hfn(l-1)
      enddo
c     =================================================================
c             l+1
c     hfn = -i   *h (k*R  )*sqrt(E)
c                  l    ij
c     =================================================================
      z=exp(sqrtm1*z)/rmag
      do l=0,lend
         hfn(l)=-hfn(l)*z*ilp1(l)     
      enddo
      endif
c
c     =================================================================
c     calculate p(l,m)'s...............................................
      costheta=rij(3)/rmag
c     -----------------------------------------------------------------
!     call plglmax(lend,costheta,plm)
!     replace with normalized associated legendre functions:
      call plm_normalized(lend,costheta,plm)
c     -----------------------------------------------------------------
c     =================================================================
c     multiply be the normalization constant...........................
      ndlm_local=(lend+1)*(lend+2)/2
      if(ndlm_local.gt.ndlm) then
	write(6,'(''MAKEGIJ:: ndlm incorrect!'')')
        write(6,*) 'ndlm=',ndlm
        write(6,*) 'ndlm_local=',ndlm_local
	call fstop(sname)
      endif
      do j=1,ndlm_local
         plm(j)=clm(j)*plm(j)
      enddo
c     =================================================================
c     calculate cos(phi) and sin(phi) .................................
      pmag=sqrt(rij(1)*rij(1)+rij(2)*rij(2))
      cosmp(0)=one
      sinmp(0)=zero
      if(pmag.ge.ptol) then
         cosmp(1)=rij(1)/pmag
         sinmp(1)=rij(2)/pmag
      else
         cosmp(1)=zero
         sinmp(1)=zero
      endif
      do m=2,lend
         cosmp(m)=cosmp(m-1)*cosmp(1)-sinmp(m-1)*sinmp(1)
         sinmp(m)=sinmp(m-1)*cosmp(1)+cosmp(m-1)*sinmp(1)
      enddo
c        
      j=0
      do l=0,lend
	 ll=l*(l+1)
	 j=ll+1
	 ll=ll/2+1
	 m1m=one
         dlm(j)= hfn(l)*plm(ll)
         do m=1,l
	    m1m=-m1m
	    fac=plm(ll+m)*dcmplx(cosmp(m),sinmp(m))
            dlm(j-m)= hfn(l)*m1m*fac
            dlm(j+m)= hfn(l)*conjg(fac)
         enddo
      enddo
c     =================================================================
      if(iprint.ge.3) then
         write(6,'(/)')
         write(6,*) "Rij=",rij(1),rij(2),rij(3)
         write(6,*) "cos(theta)=",costheta
         write(6,'(/)')
         write(6,'('' makegij:: l,m,dlm(l,m):'')')
	 do j=1,ndlj
            write(6,'(2i3,2x,f10.5,1x,d16.8)')
     >      lofk(j),mofk(j),dlm(j)
         enddo
         write(6,*) "i, plm(i)"
         do j=1,ndlm
           write(6,*) j,plm(j)
         end do
      endif
c
c     ================================================================
c     calculate g(R_ij)...............................................
c     ----------------------------------------------------------------
!     call zeroout(gij,2*kkri*kkrj)
      gij=0.0d0
c     ----------------------------------------------------------------
c     ================================================================
c     loop over l1,m1............................................
      do lm1=1,kkrj
         l1=lofk(lm1)
         m1=mofk(lm1)
c        =============================================================
c        loop over l2,m2..............................................
         do lm2=1,kkri
            l2=lofk(lm2)
            m2=mofk(lm2)
c           ==========================================================
c                            l2-l1
c           illp(lm2,lm1) = i
c
c           perform sum over l3 with gaunt # ......................
c           ==========================================================
	    m3=m2-m1
	    llow=max(abs(m3),abs(l1-l2))
	    if(prel.eq.czero) llow=l1+l2
            do l3=l1+l2,llow,-2
               j=l3*(l3+1)+m3+1
               gij(lm2,lm1) = gij(lm2,lm1)+cgnt(l3/2+1,lm1,lm2)*dlm(j)
            enddo
            gij(lm2,lm1)=pi4*illp(lm2,lm1)*gij(lm2,lm1)
         enddo
      enddo
c
c     ================================================================
c     printout if needed..............................................
      if(iprint.ge.3) then
c        =============================================================
c        loop over l1,m1..............................................
         do lm1=1,kkrj
            l1=lofk(lm1)
            m1=mofk(lm1)
c           ==========================================================
c           loop over l2,m2...........................................
            do lm2=1,kkri
               l2=lofk(lm2)
               m2=mofk(lm2)
               write(6,'('' lm1,lm2,l1,m1,l2,m2:'',6i3,2d12.5)')
     >         lm1,lm2,l1,m1,l2,m2,gij(lm2,lm1)
            enddo
         enddo
      endif
c     =================================================================
      if (istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
