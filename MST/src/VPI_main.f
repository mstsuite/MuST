!*****************************************************************************
!*****             Isoparametric integration                             *****
!*****************************************************************************
      subroutine VPI_main(iVP,Nqp,tnAtom,lattice,r,rmt,nsublat,iType,
     &                    t,nVertex1,nVertex2,Ivdex2,nfaces,
     &                    Iface,vertices,fcount,weight,rmag,vj)
!=======================================================================
! [Reference]
!  "The Finite Elements Method (linear static and dynamic finite element
!   analysis" by Thomas J.R. Hughes
!   ps: check chapter-3
!***********************************************************************
      implicit none
      include 'imp_inp.h' 
      integer*4 two,three,four
      parameter(two=2,three=3,four=4)
!--[General veriables]-------------------------------------------------
      integer*4 nAtom,tnAtom,nsublat,iType(tnAtom),n1,ITYP,iVP
      integer*4 NEA(nsublat)
      logical wholeVP,ErrTest,ETDone,ConvTest,UseSym
!--[input from Bernal]-------------------------------------------------
      real*8 vertices(tnAtom,three,MNV),r(tnAtom),rmt(nsublat)
      real*8 t(3,3),tt(3),lattice(3)
      integer*4 Ivdex2(tnAtom,MNV,two)
      integer*4 Iface(tnAtom,MNF,MNV+1),nfaces(tnAtom),fcount(nsublat)
      Integer*8 nVertex1(tnAtom),nVertex2(tnAtom)
!--[for face analysis]-------------------------------------------------
      real*8 vertices1(three,MNV)
      integer*4 vf(tnAtom,MNF,four),cut_count,fdim
!--[for abscissas & weight]--------------------------------------------     
      integer*4 Nqp,Nmin
      real*8 al,be,Qxfix(2),abscissas(MNqp),weight(MNqp)
!--[for poly_fit]------------------------------------------------------
      integer*4 jpf,error
      parameter(error=-2)
      real*8 a1pf(npol+1,tnAtom),SD
C      real*8  basis(npol+1,Nints,tnAtom)
!--[for Jacobian & print result]---------------------------------------
      integer*4 i,j,k,I_face,Nfix,IS,L2,m
CDEBUG
      integer*4 iii,jj,kptr

      real*8 TM(8,8),T1(8,8),Rmt_chk(tnAtom,MNF),diffmin
      real*8 X1(8),Y1(8),Z1(8),s1234(4),s567(4),alpha(0:7,3),rX1(3)
      real*8 vjj,func,volume,rec_volume(MNqp,4)
      real*8 rmagn,GQSum(2,tnAtom),VPInt(2)
      real*8 vj(MNqp,MNqp,MNqp,MNF,nsublat) 
      real*8 rmag(3,MNqp,MNqp,MNqp,MNF,nsublat)
!-------------------------------------------------------------------
      if(iVP.eq.0)then
        wholeVP=.false.
      else 
        wholeVP=.true.
      endif 
      ErrTest=.false.
      ConvTest=.true.
      ETDone=.false.
!
!     ETDone: An internal use flag, for error-testing
!     nsublat: Total # of different sub-lattice.
!     NEA: total # of equivalent atoms for certain sub-lattice. 
!          e.g.: NEA(2)=3 indicate the sub-lattice 2 has 3 equivalent
!                atoms.
!
      data TM/1.0D0, 1.0D0,-1.0D0, 1.0D0,  -1.0D0,-1.0D0, 1.0D0,-1.0D0,
     &        1.0D0,-1.0D0,-1.0D0, 1.0D0,   1.0D0,-1.0D0,-1.0D0, 1.0D0,
     &        1.0D0,-1.0D0, 1.0D0, 1.0D0,  -1.0D0, 1.0D0,-1.0D0,-1.0D0,
     &        1.0D0, 1.0D0, 1.0D0, 1.0D0,   1.0D0, 1.0D0, 1.0D0, 1.0D0,
     &	      1.0D0, 1.0D0,-1.0D0,-1.0D0,  -1.0D0, 1.0D0,-1.0D0, 1.0D0,
     &	      1.0D0,-1.0D0,-1.0D0,-1.0D0,   1.0D0, 1.0D0, 1.0D0,-1.0D0,
     &	      1.0D0,-1.0D0, 1.0D0,-1.0D0,  -1.0D0,-1.0D0, 1.0D0, 1.0D0,
     &	      1.0D0, 1.0D0, 1.0D0,-1.0D0,   1.0D0,-1.0D0,-1.0D0,-1.0D0/
      T1=TM/8.0D0
!
!      do iii =2,nVertex2(tnAtom)
!         write(6,'(a,i5,3(1x,f18.10))') "Vertices::",iii, 
!     &                   vertices(1,1:3,iii)
!      enddo
!
      call find_Rmt(MNV,MNF,tnAtom,nAtom,Iface,nfaces,Ivdex2,
     &                    nVertex2,vertices,Rmt_chk)
        do iii=1,tnAtom
           diffmin=100.0d0
          do jj=1,nfaces(iii)
           diffmin=min(diffmin,Rmt_chk(iii,jj)) 
           write(6,*) "Rmt_chk::",iii,jj, Rmt_chk(iii,jj)
          enddo
          do jj=1,nfaces(iii)
            if(diffmin.eq.Rmt_chk(iii,jj))then
              kptr=jj
              go to 131
            endif
          enddo  
131         r(iii)=Rmt_chk(iii,kptr)
            write(6,*) "Rmt_chk r::", r(iii)
        enddo  

	!---[put scale back]---
        t=t*lattice(1)
	r=r*lattice(1)
         do i=1,tnAtom
           rmt(iType(i))=r(i)
         enddo
	vertices=vertices*lattice(1)
	!----------------------
	if(wholeVP) r=0.0d0
      if(ConvTest) Nmin=1
       vf=0
       fcount=0
	 NEA=0
       al=0.0D0
       be=0.0D0
       Nfix=0
c       call absc_weight(Nqp,abscissas,weight)

        call gauleg(-1.0d0,+1.0d0,abscissas,weight,Nqp)       

       do 230 nAtom=1,tnAtom 
        !---[use Sym.]---

          ITYP=iType(nAtom)
          if(NEA(ITYP)>0)then
           NEA(ITYP)=NEA(ITYP)+1
	   goto 230
	  end if
          NEA(ITYP)=NEA(ITYP)+1

        call face_analysis(MNV,MNF,nAtom,tnAtom,nfaces,Iface,
     &    nVertex2,Ivdex2,nVertex1,vertices,vertices1,vf,fcount(ITYP))
        do i = 1,nfaces(nAtom)
           do j=1,Iface(nAtom,i,1)
              print* , "Iface::",i,j+1,Iface(nAtom,i,1+j)
           enddo
        enddo       
        do I_face=1,fcount(ITYP)
         call drive_jacobian(tnAtom,MNF,MNV,vertices,nAtom,vf,
     &                 I_face,r,X1,Y1,Z1,s1234,s567,T1,alpha,wholeVP)

         do i=1,Nqp
         do j=1,Nqp
         do k=1,Nqp 
          call jacobian(tnAtom,nAtom,Nqp,i,j,k,X1,Y1,Z1,s1234,s567,
     &                  r,abscissas,alpha,rX1,vjj)
C           rmag(i,j,k,I_face,ITYP)=sqrt(rX1(1)**2+rX1(2)**2+rX1(3)**2)
           rmag(1,i,j,k,I_face,ITYP)=rX1(1)
           rmag(2,i,j,k,I_face,ITYP)=rX1(2)
           rmag(3,i,j,k,I_face,ITYP)=rX1(3)           
           vj(i,j,k,I_face,ITYP) = vjj
         enddo
         enddo
         enddo
        enddo

  230  continue

      return
      end 

!===================================================================
!     face_analysis(MNV,MNF,nAtom,nfaces,Iface,nVertex2,    !
!                   Ivdex2,nVertex1,vertices,vertices1,vf,fcount)  !
!===================================================================
!**********************************************************************
! nVertex1: number of vertices for "non-degenerate" part
! vertices: coordinate of the vertices for "non-degenerate" part
! Ivdex1: index of the vertices for vertices
!-----------
! nVertex2: number of vertices for "degenerate" part
! Ivdex2: index of the vertices for degenerate vertices
!-----------
! nfaces: number of faces
! Iface(32,nfaces): face information includes
!           ==>1st row: face index
!              2dn row: number of vertices with possible duplication; 
!              3rd-32th row: all vertex-index of this face.
!              32 rows is just some large number, here, 32 means we could have 
!              maximum (32-2) vertices for each face.
! vertices: vertices{x,y,z} for each plane with counter-clockwise sequence
!**********************************************************************

        subroutine face_analysis(MNV,MNF,nAtom,tnAtom,nfaces,Iface,
     &         nVertex2,Ivdex2,nVertex1,vertices,vertices1,vf,fcount)

      implicit none
       integer*4 two,three,four
       parameter(two=2,three=3,four=4)
      integer*4 i,j,k,MNV,MNF,nAtom,tnAtom,L2,n1
      integer*4 Ivdex2(tnAtom,MNV,two)
      integer*4 Iface(tnAtom,MNF,MNV+1),nfaces(tnAtom)
      integer*8 nVertex1(tnAtom)
      integer*8 nVertex2(tnAtom)
      integer*4 vf(tnAtom,MNF,four),Vdex0(MNV)
      integer*4 cut_count,fcount,FI,faceDex(MNF,four)
      real*8 vertices(tnAtom,three,MNV),vertices1(three,MNV)

      do 403 i=1,nfaces(nAtom) !faces
       do 404 j=2,31 !vertices index
        if(Iface(nAtom,i,j)==0) goto 403
        L2=Iface(nAtom,i,j)
        do 405 k=1,nVertex2(nAtom)
          if(L2==Ivdex2(nAtom,k,1)) Iface(nAtom,i,j)=Ivdex2(nAtom,k,2)
          !replace the degenerate vertex to non-degenerate vertex
  405   continue
  404  continue
  403 continue

      !elliminate the duplicate vertices
      do 406 i=1,nfaces(nAtom) !faces
       do 407 j=3,31 !start from 3 because the first column is not vertex
  409   do 408 k=j-1,2,-1
          if(Iface(nAtom,i,j)==Iface(nAtom,i,k))then
            Iface(nAtom,i,j:30)=Iface(nAtom,i,j+1:31)
            !if there is one vertex repeat, over write it by the remain column
            goto 409                !then re-compare it again==> goto 407
          else if(Iface(nAtom,i,j)==0)then
            Iface(nAtom,i,1)=j-2
            !substract 2 is because the first column is "number of 
            ! vertices with possible duplication",	so if the element 
            !of 5th-vertex(which is the 6th-column of Iface) is 0, 
            !then 6(6th)-2=4 ==> 4 sides for this plane
            goto 406 !when it hit "0", change to next face
          end if
  408   continue
  407  continue
  406 continue
!=======================================================================

      fcount=0
      do FI=1,nfaces(nAtom) !Number of Face
       Vdex0=0
       fcount=fcount+1
       if(Iface(nAtom,FI,1)==3)then
        vf(nAtom,fcount,1:3)=Iface(nAtom,FI,2:4)
	  !made the 3rd- & 4th-vertex be degenerate
        vf(nAtom,fcount,4)=Iface(nAtom,FI,4)     
	  !ps: triangle don't need to re-arrange the order
       else if(Iface(nAtom,FI,1)==4)then
        Vdex0(1:Iface(nAtom,FI,1))=Iface(nAtom,FI,2:Iface(nAtom,FI,1)+1)
        call determineSQ(Iface(nAtom,FI,1),MNV,MNF,Vdex0, !re-arrange the order
     &            vertices(nAtom,:,:),Iface(nAtom,:,:),FI,vertices1)
        vf(nAtom,fcount,1:4)=Vdex0(1:4)
       else if(Iface(nAtom,FI,1)>4)then
        Vdex0(1:Iface(nAtom,FI,1))=Iface(nAtom,FI,2:Iface(nAtom,FI,1)+1)
        call determineSQ(Iface(nAtom,FI,1),MNV,MNF,Vdex0, !re-arrange the order
     &            vertices(nAtom,:,:),Iface(nAtom,:,:),FI,vertices1)
      !call bestcut(Iface(nAtom,FI,1),MNV,MNF,recVdex,cut_count,Vdex0,vertices1)
        call simple_cut(Iface(nAtom,FI,1),MNV,MNF,vdex0,cut_count,
     &                  faceDex)
        !vf(nAtom,fcount:fcount+cut_count-1,:)=recVdex(1:cut_count,cut_count,1:4)
        vf(nAtom,fcount:fcount+cut_count-1,:)=faceDex(1:cut_count,:)
        fcount=fcount+cut_count-1
       end if
      enddo

      return
      end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine determineSQ(NN,MNV,MNF,Vdex0,vertices,Iface,
     &                       FI,vertices1)
      implicit none
      integer*4 NN,MNV,MNF !NN-->nVertex
      integer*4 I,J,K,L,Nap,N1,N2,Ibri,Iface(MNF,MNV+1),FI,Vdex0(MNV)
      real*8 vertices(3,MNV),vertices1(3,MNV),vector(3,NN),origin(3)
      real*8 theta(NN),rx(3),ry(3),rz(3),zN,yN,Rot(3,3),x,y,pi,bri

      pi=4.0D0*datan(1.0D0)
      vertices1=0.0D0
      origin=vertices(:,Iface(FI,2)) !define 1st-column(1st vertex) to be origin
      do i=1,NN !define a relative-position vector
      vertices1(:,i)=vertices(:,Iface(FI,i+1))-origin
      end do

      call cross(vertices1(:,2),vertices1(:,3),rz)
      ry=vertices1(:,2) 
	!define vector-ab to be y, so x can't be zero, then we could use tangent to define angle
      zN=sqrt(dot_product(rz,rz))
      yN=sqrt(dot_product(ry,ry))
      rz=rz/zN
      ry=ry/yN
      call cross(ry,rz,rx)
      Rot(1,1:3)=rx
      Rot(2,1:3)=ry
      Rot(3,1:3)=rz
      call dMRRRR(3,3,Rot,3,3,NN,vertices1(:,1:NN),3,3,NN,vector,3)
      !!!call dwrrrl('   Rot   ',3,3,Rot,3,0,'(f6.3)','number','number')
      vertices1(:,1:NN)=vector !here, we turn the vertices into x-y plane
      ! call dwrrrl('   vertices(after)   ',3,NN,vertices1,3,0,'(f16.13)','
      !number','number')
      !define the angle between each vector(:,i) and x-axis
      theta(1)=0.0d0
      theta(2)=pi/2.0d0
      do i=3,NN
       x=vertices1(1,i)
       y=vertices1(2,i)
       if(y==0)then
        if(x>0) theta(i)=0.0d0
        if(x<0) theta(i)=pi
       else if(y>0)then
        if(x>0) theta(i)=datan(y/x)
        if(x<0) theta(i)=pi+datan(y/x) !in this region, datan(y/x)<0
       else if(y<0)then
        if(x>0) theta(i)=2.0d0*pi+datan(y/x) !in this region, datan(y/x)<0
        if(x<0) theta(i)=pi+datan(y/x)
       end if
      enddo

                   !determine the sequence of the vertices
         do i=2,NN !the first vertex had set to be first in the sequence
  100     do j=i+1,NN
           if(theta(j)<theta(i))then
            !exchange the sequence of theta
            bri=theta(i)
            theta(i)=theta(j)
            theta(j)=bri
            !exchange the sequence of Vdex0
            Ibri=Vdex0(i)
            Vdex0(i)=Vdex0(j)
            Vdex0(j)=Ibri
  	        goto 100
	       end if
	      enddo
         enddo
            !!!call dwrrrl('   theta(after)   ',NN,1,theta,NN,0,'(f8.3)','
            ! number','number')
            !!!call wrirl('   vdex in sequence   ',NN,1,Vdex0,NN,0,'(I2)','
            !none','none')
      return
      end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cross(a,b,c) !a (cross) b = c
      implicit none
      real*8 a(3),b(3),c(3)

      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)

      return
      end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine simple_cut(NN,MNV,MNF,vdex0,cut_count,faceDex)
      implicit none
      integer*4 i,j,NN,MNV,MNF,Vdex0(MNV),faceDex(MNF,4)
	integer*4 cut_count,temp(4)
      faceDex=0
      cut_count=0

      j=NN
      temp(1)=Vdex0(1)
      do 100 i=1,NN-3
       if(j>=4)then
        temp(2:4)=Vdex0(2*i:2*i+2)
        faceDex(i,:)=temp
        j=j-2
        if(j==2) goto 10
       else if(j==3)then
        temp(2:3)=Vdex0(2*i:2*i+1)
        temp(4)=Vdex0(2*i+1) !degenerate point
        faceDex(i,:)=temp
        goto 10
       endif
  100 continue
   10 cut_count=i
      ! write(6,*)  cut_count
      return
      end 
!===================================================================
!     drive_jacobian(MNF,MNV,vertices,nAtom,vf,I_face,r,           !
!                    X1,Y1,Z1,s1234,s567,T1,alpha,,wholeVP)        !
!===================================================================
      subroutine drive_jacobian(tnAtom,MNF,MNV,vertices,nAtom,vf,
     &                 I_face,r,X1,Y1,Z1,s1234,s567,T1,alpha,wholeVP)
      implicit none
      integer*4 tnAtom,MNF,MNV,nAtom,vf(tnAtom,MNF,4)
	integer*4 I_face,IS,iRmt,iDis,pwr1
      real*8 vertices(tnAtom,3,MNV)
	real*8 X1(8),Y1(8),Z1(8),X2(8),Y2(8),Z2(8)
      real*8 L_scale(4),r(tnAtom),s1234(4),s567(4),t,T1(8,8)
      real*8 XYZ2(8,3),alpha(8,3),distance,shift
	logical wholeVP

      X1(1)=vertices(nAtom,1,vf(nAtom,I_face,1))
      X1(2)=vertices(nAtom,1,vf(nAtom,I_face,2))
      X1(3)=vertices(nAtom,1,vf(nAtom,I_face,3))
      X1(4)=vertices(nAtom,1,vf(nAtom,I_face,4))

      Y1(1)=vertices(nAtom,2,vf(nAtom,I_face,1))
      Y1(2)=vertices(nAtom,2,vf(nAtom,I_face,2))
      Y1(3)=vertices(nAtom,2,vf(nAtom,I_face,3))
      Y1(4)=vertices(nAtom,2,vf(nAtom,I_face,4))

      Z1(1)=vertices(nAtom,3,vf(nAtom,I_face,1))
      Z1(2)=vertices(nAtom,3,vf(nAtom,I_face,2))
      Z1(3)=vertices(nAtom,3,vf(nAtom,I_face,3))
      Z1(4)=vertices(nAtom,3,vf(nAtom,I_face,4))


      L_scale(1)=r(nAtom)/sqrt(x1(1)**2+y1(1)**2+z1(1)**2)
      L_scale(2)=r(nAtom)/sqrt(x1(2)**2+y1(2)**2+z1(2)**2)
      L_scale(3)=r(nAtom)/sqrt(x1(3)**2+y1(3)**2+z1(3)**2)
      L_scale(4)=r(nAtom)/sqrt(x1(4)**2+y1(4)**2+z1(4)**2)
      do IS=5,8
      X1(IS)=L_scale(IS-4)*x1(IS-4)
      Y1(IS)=L_scale(IS-4)*y1(IS-4)
      Z1(IS)=L_scale(IS-4)*z1(IS-4)
      end do

!--------[generate the coefficient of plane equation]--------
      s1234(1)=(Y1(1)-Y1(2))*(Z1(2)-Z1(3))-(Z1(1)-Z1(2))*(Y1(2)-Y1(3))
      s1234(2)=(Z1(1)-Z1(2))*(X1(2)-X1(3))-(X1(1)-X1(2))*(Z1(2)-Z1(3))
      s1234(3)=(X1(1)-X1(2))*(Y1(2)-Y1(3))-(Y1(1)-Y1(2))*(X1(2)-X1(3))
      s1234(4)=-(s1234(1)*X1(4)+s1234(2)*Y1(4)+s1234(3)*Z1(4))

      s567(1)=(Y1(5)-Y1(6))*(Z1(6)-Z1(7))-(Z1(5)-Z1(6))*(Y1(6)-Y1(7))
      s567(2)=(Z1(5)-Z1(6))*(X1(6)-X1(7))-(X1(5)-X1(6))*(Z1(6)-Z1(7))
      s567(3)=(X1(5)-X1(6))*(Y1(6)-Y1(7))-(Y1(5)-Y1(6))*(X1(6)-X1(7))
      s567(4)=-(s567(1)*X1(7)+s567(2)*Y1(7)+s567(3)*Z1(7))
!-----[check the muffin]
      distance=-s1234(4)/sqrt(dot_product(s1234(1:3),s1234(1:3)))
	distance=dabs(distance)
        pwr1=6
	shift=dble(10**pwr1)
C	iRmt=kidint(r(nAtom)*shift)
C	iDis=kidint(distance*shift)
C       DEBUG
C       if(iRmt>iDis) stop 'muffin-tin radius exceed the VP boundary!'
	!If lattice>>1, we might need to decrease 'pwr1' to save digits
	!for integer part.
!---[Generate X2]-----------------------
      !if(r(nAtom)==0.0D0)then
	if(wholeVP)then
       X2=X1
       Y2=Y1
       Z2=Z1
       goto 210 ! skip the modify
      end if
      t=-s567(4)/(s567(1)*X1(8)+s567(2)*Y1(8)+s567(3)*Z1(8))
      !assume t[X1(8),Y1(8),Z1(8)] is the solution for ax+by+cz+d=0
      !ps: t[X1(8),Y1(8),Z1(8)] --> [X2(8),Y2(8),Z2(8)] which is what we need
      !-->t[a*X1(8)+b*Y1(8)+c*Z1(8)]+d=0
      !so t=-d/[a*X1(8)+b*Y1(8)+c*Z1(8)]
      X2(8)=t*X1(8)
      Y2(8)=t*Y1(8)
      Z2(8)=t*Z1(8)    ! we choose point-5, -6 and -7 to forming a plane, so,
      X2(1:7)=X1(1:7)  ! point 8 may not be on plane 567, thus we modify point 8
      Y2(1:7)=Y1(1:7)
      Z2(1:7)=Z1(1:7)
  210 continue

      XYZ2(:,1)=X2
      XYZ2(:,2)=Y2
      XYZ2(:,3)=Z2
      !call dwrrrn('XYZ2=',8,3,XYZ2,8,0)
      !call dwrrrl('XYZ2=',8,3,XYZ2,8,0,'(f20.16)','number','number')

      call dMRRRR(8,8,T1,8,8,3,XYZ2,8,8,3,alpha,8) !attain the coefficients
      !call dwrrrl('alpha=',8,3,alpha,8,0,'(f20.16)','number','number')

      return
      end 

!===================================================================
!    jacobian(Nqp,Nqp,i,j,k,X1,Y1,Z1,s1234,s567,r,abscissas,        !
!             alpha,rX1,vj)                                        !
!===================================================================
      subroutine jacobian(tnAtom,nAtom,Nqp,i,j,k,X1,Y1,Z1,
     &                    s1234,s567,r,abscissas,alpha,rX1,vj)
      implicit none
      integer*4 Nqp,i,j,k,tnAtom,nAtom
	real*8 pi,zero
      real*8 abscissas(Nqp),Xi,Eta,Zeta,Xi_vector(8),alpha(0:7,3)
      real*8 X1(8),Y1(8),Z1(8),rX1(3),rX2(3),r1,r2,r(tnAtom)
      real*8 vj1,vj2,vj3,vj,s1234(4),s567(4),dir(3),vl1,vl2
	real*8 theta1,phi1,theta2,phi2
      real*8 al1467,be1467,ga1467,al2457,be2457,ga2457
	real*8 al3567,be3567,ga3567

      pi=4.0D0*datan(1.0D0)
      zero=1.0D-14 ! 14th decimal place....

      Xi=abscissas(i)
      Eta=abscissas(j)
      Zeta=abscissas(k)
!---------------------------------------
!---[get coordinate rX2]----------------
      Xi_vector(1)=1.0D0
      Xi_vector(2)=Xi
      Xi_vector(3)=Eta
      Xi_vector(4)=Zeta
      Xi_vector(5)=Xi*Eta
      Xi_vector(6)=Eta*Zeta
      Xi_vector(7)=Zeta*Xi
      Xi_vector(8)=Xi*Eta*Zeta

      call dMRRRR(1,8,Xi_vector,1,8,3,alpha,8,1,3,rX2,1)

      if(r(nAtom)==0.0D0)then 
      vJ1=1.0D0
      vJ2=1.0D0
      rX1=rX2
      goto 90  ! skip Transformation(2)
      end if

!---[get coordinate r2]----------------
      r2=sqrt(rX2(1)**2+rX2(2)**2+rX2(3)**2)

      if(dabs(rX2(1))<zero .and. dabs(rX2(2))<zero)then
	 !x=y=0 ==> theta2=phi2=0.0 ==> dir=(0,0,+1) or (0,0,-1)
       dir(3)=1.0D0
       if(rX2(3)<0.0D0)dir(3)=-1.0D0
       vl1=-s567(4)/(s567(3)*dir(3))
       vl2=-s1234(4)/(s1234(3)*dir(3))
       r1=(vl2*(r(nAtom)-vl1)+r2*(vl2-r(nAtom)))/(vl2-vl1)
       vJ1=r1**2/r2**2
       rX1(1:2)=0.0D0
       rX1(3)=r1*dir(3)
       goto 85
      else if(dabs(rX2(1))<zero .and. rX2(2)>zero)then !x=0, y>0
       theta2=dacos(rX2(3)/r2)
       phi2=pi/2.0D0
       goto 80
      else if(dabs(rX2(1))<zero .and. rX2(2)<zero)then !x=0, y<0
       theta2=dacos(rX2(3)/r2)
       phi2=pi*3.0D0/2.0D0
       goto 80
      end if

      theta2=dacos(rX2(3)/r2)

      if(rX2(1)>0.0D0)then !x>0, normal
       phi2=datan(rX2(2)/rX2(1))
      else if(rX2(1)<0.0D0)then !x<0, +pi
       phi2=datan(rX2(2)/rX2(1))+pi
      end if


   80 continue
!---[get coordinate r1]----------------
      dir(1)=dsin(theta2)*dcos(phi2)
      dir(2)=dsin(theta2)*dsin(phi2)
      dir(3)=dcos(theta2) !direction vector
      !r=sqrt(X1(5)**2+Y1(5)**2+Z1(5)**2) !radius of inner sphere

      vl1=-s567(4)/(s567(1)*dir(1)+s567(2)*dir(2)+s567(3)*dir(3))
      vl2=-s1234(4)/(s1234(1)*dir(1)+s1234(2)*dir(2)+s1234(3)*dir(3))
      r1=(vl2*(r(nAtom)-vl1)+r2*(vl2-r(nAtom)))/(vl2-vl1)
      theta1=theta2
      phi1=phi2

!---[get coordinate rX1]----------------
      rX1(1)=r1*dir(1)
      rX1(2)=r1*dir(2)
      rX1(3)=r1*dir(3) !for calculate function value

      !call dwrrrl('rX1',3,1,rX1,3,0,'(f20.16)','number','number')

!---------------------------------------
!---------------------------------------
      vJ1=r1**2*dsin(theta1)/(r2**2*dsin(theta2))
   85 continue
      vJ2=(vl2-r(nAtom))/(vl2-vl1)

      !because the notation, we set the row of alpha from 0-7;==> alpha(0:7,3) not alpha(8,3)
   90 continue
      al1467=alpha(1,1)+alpha(4,1)*eta+alpha(6,1)*zeta+
     &       alpha(7,1)*eta*zeta
      be1467=alpha(1,2)+alpha(4,2)*eta+alpha(6,2)*zeta+
     &       alpha(7,2)*eta*zeta
      ga1467=alpha(1,3)+alpha(4,3)*eta+alpha(6,3)*zeta+
     &       alpha(7,3)*eta*zeta

      al2457=alpha(2,1)+alpha(4,1)*xi+alpha(5,1)*zeta+alpha(7,1)*xi*zeta
      be2457=alpha(2,2)+alpha(4,2)*xi+alpha(5,2)*zeta+alpha(7,2)*xi*zeta
      ga2457=alpha(2,3)+alpha(4,3)*xi+alpha(5,3)*zeta+alpha(7,3)*xi*zeta

      al3567=alpha(3,1)+alpha(5,1)*eta+alpha(6,1)*xi+alpha(7,1)*eta*xi
      be3567=alpha(3,2)+alpha(5,2)*eta+alpha(6,2)*xi+alpha(7,2)*eta*xi
      ga3567=alpha(3,3)+alpha(5,3)*eta+alpha(6,3)*xi+alpha(7,3)*eta*xi

      vJ3=(al1467*be2457*ga3567+al2457*be3567*ga1467+
     &     al3567*be1467*ga2457)-(al3567*be2457*ga1467+
     &     al2457*be1467*ga3567+al1467*be3567*ga2457)

      vJ=vJ1*vJ2*vJ3
      vJ=dabs(vJ)

      return
      end 

!===================================================================
!    subroutine find_Rmt(MNV,MNF,tnAtom,nAtom,Iface,nfaces,Ivdex2, !
!     &                    nVertex2,vertices)                      !
!===================================================================
      subroutine find_Rmt(MNV,MNF,tnAtom,nAtom,Iface,nfaces,Ivdex2,
     &                    nVertex2,vertices,Rmt_chk)
	implicit none
	integer*4 i,j,MNV,MNF,tnAtom,nAtom,NF,Iface(tnAtom,MNF,MNV+1)
	integer*4 v1,v2,v3,v4,nfaces(tnAtom),Ivdex2(tnAtom,MNV,2)
	integer*4 Iface_dup(tnAtom,MNF,MNV+1)
	integer*8 nVertex2(tnAtom)
	real*8 vertices(tnAtom,3,MNV),Rmt_chk(tnAtom,MNF)
        real*8 x1(4),y1(4),z1(4)
	real*8 s1234(4)

      Rmt_chk=0.0d0
      Iface_dup=Iface ! cause we need to make some change that may cause
                      ! problem after this stage
      do 10 nAtom=1,tnAtom
	 do 20 NF=1,nfaces(nAtom)
	  do 30 i=2,Iface(nAtom,NF,1)+1
	   do 40 j=1,nVertex2(nAtom)
	    if(Iface(nAtom,NF,i)==Ivdex2(nAtom,j,1))then
	     Iface(nAtom,NF,i)=Ivdex2(nAtom,j,2)
	    end if
   40    continue
   30	  continue
   20  continue
   10	continue

      do 100 nAtom=1,tnAtom
       do 200 NF=1,nfaces(nAtom)
         v1=Iface_dup(nAtom,NF,1+1) !+1 is because the 1st column is # of vertex
         v2=Iface_dup(nAtom,NF,2+1) !not the vdex(vertices-index)
         v3=Iface_dup(nAtom,NF,3+1)
         v4=Iface_dup(nAtom,NF,4+1)
        
        write(6,'(a,4(i5))') "v1-4::",v1,v2,v3,v4
        write(6,'(3(f22.14))') vertices(nAtom,1,v1)
        write(6,'(3(f22.14))') vertices(nAtom,1,v2)
        X1(1)=vertices(nAtom,1,v1)
        X1(2)=vertices(nAtom,1,v2)
        X1(3)=vertices(nAtom,1,v3)
        X1(4)=vertices(nAtom,1,v4)
        
        Y1(1)=vertices(nAtom,2,v1)
        Y1(2)=vertices(nAtom,2,v2)
        Y1(3)=vertices(nAtom,2,v3)
        Y1(4)=vertices(nAtom,2,v4)

        Z1(1)=vertices(nAtom,3,v1)
        Z1(2)=vertices(nAtom,3,v2)
        Z1(3)=vertices(nAtom,3,v3)
        Z1(4)=vertices(nAtom,3,v4)
        write(6,*) "XYZ::", X1(1),Y1(1),Z1(1)
        write(6,*) "XYZ::", X1(4),Y1(4),Z1(4)
        s1234(1)=(Y1(1)-Y1(2))*(Z1(2)-Z1(3))-(Z1(1)-Z1(2))*(Y1(2)-Y1(3))
        s1234(2)=(Z1(1)-Z1(2))*(X1(2)-X1(3))-(X1(1)-X1(2))*(Z1(2)-Z1(3))
        s1234(3)=(X1(1)-X1(2))*(Y1(2)-Y1(3))-(Y1(1)-Y1(2))*(X1(2)-X1(3))
        s1234(4)=-(s1234(1)*X1(4)+s1234(2)*Y1(4)+s1234(3)*Z1(4))
        write(6,*) "s1234::", s1234(1:4)
       Rmt_chk(nAtom,NF)=-s1234(4)/
     *                    sqrt(dot_product(s1234(1:3),s1234(1:3)))
  200  continue
  100 continue
      Rmt_chk=dabs(Rmt_chk)

	return
	end 

!=======================================================================
C    Calculate the abscissa and weights for Gauss Quadr. Integration
!=======================================================================

       SUBROUTINE gauleg(x1,x2,x,w,n)
c     ------------------------------------------------------------------
c     given the lower and upper limits of integration x1 and x2, and
c     given n, this routine returns arrays x(1:n) and w(1:n) containing
c     the abscissas and weights of the gauss-legendre quadrature
c     formulae
c          int(x1,x2) dx f(x) = sum(j=1,n) w[j] f(x[j])
c     ------------------------------------------------------------------
       INTEGER n
       REAL*8 x1,x2,x(n),w(n)
       DOUBLE PRECISION EPS
       PARAMETER (EPS=3.d-14)
       INTEGER i,j,m
       DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1

       m=(n+1)/2
       xm=0.5d0*(x2+x1)
       xl=0.5d0*(x2-x1)
       do 12 i=1,m
         z=cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
1       continue
           p1=1.0d0
           p2=0.0d0
           do 11 j=1,n
             p3=p2
             p2=p1
             p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
11        continue
           pp=n*(z*p1-p2)/(z*z-1.0d0)
           z1=z
           z=z1-p1/pp
         if(abs(z-z1).gt.EPS)goto 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.0d0*xl/((1.0d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
12    continue
       return
       END

!===================================================================
       SUBROUTINE DMRRRR(NRA,NCA,A,LDA,NRB,NCB,B,LDB,
     +               NRC,NCC,C,LDC)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION A(LDA,*),B(LDB,*),C(LDC,*)
       DO 10 I=1,NRA
          DO 20 J=1,NCB
             SUM=0D0
             DO 30 K=1,NCA
                SUM=SUM+A(I,K)*B(K,J)
 30          CONTINUE
             C(I,J)=SUM
 20       CONTINUE
 10    CONTINUE
       RETURN
       END





