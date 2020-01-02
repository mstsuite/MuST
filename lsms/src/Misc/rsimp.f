      function rsimp(f,r,irn,dx)
c===================
c     radial integration via simpson
c
      implicit real*8 (a-h,o-z)
c
      dimension f(*),r(*)
c
      isw=0
      rsimp=0.d0
      if(irn.le.2) return
      if(irn/2*2.eq.irn) isw=1
      np=irn-isw
      s=f(1)*r(1)+f(np)*r(np)
      nl=np-1
      do 5 i=2,nl,2
    5 s=s+4.d0*f(i)*r(i)
      nl=nl-1
      if(nl.lt.3) goto 15
      do 10 i=3,nl,2
   10 s=s+2.d0*f(i)*r(i)
   15 s=s*dx/3.d0
      if(isw.eq.1) goto 30
      rsimp=s
      return
   30 rsimp=s+(f(irn)*r(irn)+f(irn-1)*r(irn-1))*0.5d0*dx
      return
      end
