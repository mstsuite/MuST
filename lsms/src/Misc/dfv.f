      subroutine dfv(r,y,dy,nv,n,rn,ecomp,vr,vso,
     &           nsizm,isymm,eunit,b,allp1)
c     =============
c
c------ compute derivatives
c
      implicit real*8 (a-h,o-z)
c
      complex*16 ecomp
c
!     common /dfbits/ b,allp1
c
      real*8 y(nv),dy(nv),rn,vr(n)
      real*8 eunit
c     ==================================================
c     dimensioned to stop run time array bounds problems
c     ==================================================
      integer isymm(nsizm)
      real*8  vso(nsizm,nsizm)
c
      s1or=1.0d0/r
      s1or2=s1or*s1or
      ereal=dreal(ecomp)
      eimag=dimag(ecomp)
c
      r1=r-rn
      v=vr(1)*r1+vr(2)
      if(n.ge.3) then
	v=v/(1.d0+vr(3)*r1*r1)
	if(n.eq.4) v=v+vr(n)
      endif
c
      s2vmer=eunit*(v*s1or-ereal)
      frelr=1.0d0-b*s2vmer
      if(nv.eq.2) then
	if(b.eq.0.d0) then
	if(y(1).eq.0.d0.and.y(2).eq.0.d0) then
	  dy(1)=frelr
	  dy(2)=allp1*s1or2/frelr+s2vmer
        else
          facr=allp1*s1or2/frelr+s2vmer
          dy(1)=frelr*y(2)+s1or*y(1)
          dy(2)=facr*y(1)-y(2)*s1or
	endif
        else
	  escale=0.5d0*eunit
	  facr=allp1*s1or
          s2vmer=escale*(v*s1or-ereal)/b
	if(y(1).eq.0.d0.and.y(2).eq.0.d0) then
	  dy(1)=-s2vmer+2.d0*b
	  dy(2)=facr*facr-s2vmer*(s2vmer-2.d0*b)
	else
          dy(1)=(-s2vmer+2.d0*b)*y(2)-facr*y(1)
          dy(2)=s2vmer*y(1)+facr*y(2)
	endif
	endif
      else
      s2vmei=       -eunit*eimag
      freli=     -b*s2vmei
      den=allp1/(frelr*frelr+freli*freli)
      den=den*s1or2
      facr=+den*frelr+s2vmer
      faci=-den*freli+s2vmei
      dy(1)=frelr*y(3)-freli*y(4)+s1or*y(1)
      dy(2)=freli*y(3)+frelr*y(4)+s1or*y(2)
      dy(3)=facr*y(1)-faci*y(2)-y(3)*s1or
      dy(4)=faci*y(1)+facr*y(2)-y(4)*s1or
      endif
c
      return
      end
