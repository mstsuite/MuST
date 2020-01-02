      subroutine caltnode(r,xp,nbnd,
     >                    edge,edgp,nedge,
     >                    begth,endth,
     >                    np_index,tnode,node)
c     ================================================================
c
      implicit   none
c
      integer    nbnd
      integer    nedge
      integer    node
      integer    np_index(2)
      integer    sigma
      integer    i,n,i1
c
      real*8     r
      real*8     edge(3,nedge)
      real*8     edgp(3,nedge)
      real*8     dn2
      real*8     tnode(nbnd+2*nedge+2)
      real*8     dp2,rhop
      real*8     xp(3,nbnd)
      real*8     begth,endth
      real*8     t,cost2,e2
      real*8     one,tol
      real*8     x0
      real*8     y0
      real*8     z0
      real*8     a,r2
c
      parameter  (one=1.0d0)
      parameter  (tol=1.0d-14)
c     parameter  (tol=1.0d-8)
c
c     ================================================================
      if(begth+one.gt.-tol) then
         n=1
         tnode(1)=begth
      else
         n=0
      endif
      r2=r*r
      do i=1,nbnd
	 rhop=xp(1,i)*xp(1,i)+xp(2,i)*xp(2,i)
	 dp2=rhop+xp(3,i)*xp(3,i)
         if (r2.gt.dp2) then
           if(np_index(1).ne.i .and. np_index(2).ne.i) then
            cost2 = sqrt((r2-dp2)*rhop/dp2)
c
	    do i1=1,2
            z0=xp(3,i)+cost2
	    if(z0.gt.r) z0=r
	    if(z0.lt.-r) z0=-r
            a=sqrt((r2-z0*z0)/rhop)
            x0=xp(1,i)*a
            y0=xp(2,i)*a
	    if(sigma(x0,y0,z0,xp,nbnd,1).eq.1.or.
     &             sigma(-x0,-y0,z0,xp,nbnd,1).eq.1) then
               n=n+1
               tnode(n)=z0/r
            endif
	    cost2=-cost2
	    enddo
c
           end if
         end if
      enddo
      if (n.gt.1) then
	call cleanup(tnode,n,tol,endth)
      end if
c
c     ================================================================
c     look for tnode by solving:
c
c            ->       ->  2    2
c          ( p  + t * e  )  = r
c             i        i
c
c                        ->       ->
c          r * tnode = ( p  + t * e  )
c                         i        i  z
c     ================================================================
      do i=1,nedge
         dn2=edgp(1,i)*edgp(1,i)+edgp(2,i)*edgp(2,i)+edgp(3,i)*edgp(3,i)
         if (r2.gt.dn2) then
            e2=edge(1,i)*edge(1,i)+edge(2,i)*edge(2,i)+
     &         edge(3,i)*edge(3,i)
            t=sqrt((r2-dn2)/e2)
c
	    do i1=1,2
            x0=edgp(1,i)+t*edge(1,i)
            y0=edgp(2,i)+t*edge(2,i)
            z0=edgp(3,i)+t*edge(3,i)
            if(sigma(x0,y0,z0,xp,nbnd,1).eq.1) then
               n=n+1
               tnode(n)=z0/r
            endif
	    t=-t
	    enddo
         end if
      enddo
c
      if (n.gt.1) then
	call cleanup(tnode,n,tol,endth)
      end if
      if(endth-one.gt.tol) then
	 node=n
      else
         n=n+1
         tnode(n)=endth
	 node=n
      endif
c        write(6,'(/,''  the number of tnode = '',i4)')node
c        write(6,'(  ''  tnode:: '',i4,2x,1e15.8)')
c    >                             (i,tnode(i),i=1,node)
c
      return
      end

      subroutine cleanup(tnode,node,tol,endth)
      implicit real*8 (a-h,o-z)
      real*8 tnode(node)
      n=1
c        -------------------------------------------------------------
         call sort(node,tnode)
c        -------------------------------------------------------------
      do i=2,node
         if ( abs(tnode(i)-tnode(n)).gt.tol .and.
     >        tnode(i).lt.endth) then
            n=n+1
            tnode(n)=tnode(i)
         end if
      end do
      node=n

      return
      end
