c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnewtc(diamtx,tcpa,tab,
     >                  xab,wspace,
c    >                  w1,w2,w3,w4,
c    >                  v1,v2,v3,
     >                  kkrsz,
     >                  atcon,
     >                  komp,mxcomp,
     >                  errcpa,
     >                  iprint,istop)
c     ================================================================
c
c     ****************************************************************
c
c     called by: 
c     calls    : 
c
c     input:
c                tau00   (tauc)
c                tcpa    (tc)
c                tab     (t-matrix for atoms a & b)
c                conc    (concentration of species a)
c                kkrsz   (size of KKR-matrix)
c                komp    (number of components on sublattice)
c                istop   (index of subroutine prog. stops in)
c                roote   (energy in ryd.)
c     output:
c                tcpa    (next iterate for tc)
c                errcpa  (error in cpa solution (tcold-tcnew) )
c
c     method:
c                iteration scheme uses
c                                         -1                -1
c            mc(n+1) = mc(n) - [ xc(mc(n))   + tauc(mc(n)) ]
c
c                where,
c                                                -1             -1  
c            xc = sum over i { c(i) * [ (mc-m(i))   - tauc(mc) ]  }
c
c     mc(n+1) is obtained by assuming that it can obtain the same
c     xc that was obtained by using m(i). it is like an "equivalent
c     ata." this solution must always be started with an ata to
c     ensure convergence, as shown by bob mills and co-workers.
c     ****************************************************************
c
      implicit   none
c
      character  istop*10
      character  sname*10
c
      integer    komp
      integer    mxcomp,kkrsz
      integer    i,ic
      integer    iprint
c
      real*8     atcon(mxcomp)
      real*8     errcpa
      real*8     zero
      real*8     one
      real*8     onem
c
      complex*16 diamtx(kkrsz,kkrsz)
      complex*16 tcpa(kkrsz*kkrsz)
      complex*16 tab(kkrsz*kkrsz,mxcomp)
      complex*16 xab(kkrsz*kkrsz,mxcomp)
      complex*16, target :: wspace(3*kkrsz+4*kkrsz*kkrsz)
      complex*16, pointer :: v1(:)
      complex*16, pointer :: v2(:)
      complex*16, pointer :: v3(:)
      complex*16, pointer :: w1(:)
      complex*16, pointer :: w2(:)
      complex*16, pointer :: w3(:)
      complex*16, pointer :: w4(:)
      complex*16 det
c
      parameter  (zero=0.0)
      parameter  (one=1.0)
      parameter  (onem=-1.0)
      parameter  (sname='mnewtc')
c
      v1 => wspace(1:kkrsz)
      v2 => wspace(kkrsz+1:2*kkrsz)
      v3 => wspace(2*kkrsz+1:3*kkrsz)
      w1 => wspace(3*kkrsz+1:3*kkrsz+kkrsz*kkrsz)
      w2 => wspace(3*kkrsz+kkrsz*kkrsz+1:3*kkrsz+2*kkrsz*kkrsz)
      w3 => wspace(3*kkrsz+2*kkrsz*kkrsz+1:3*kkrsz+3*kkrsz*kkrsz)
      w4 => wspace(3*kkrsz+3*kkrsz*kkrsz+1:3*kkrsz+4*kkrsz*kkrsz)
c
c     *****************************************************************
c
      errcpa=zero
c
c     ================================================================
c     tcpa => t_{c}(old) \equiv t_{c}.................................
c     w4   => t_{c}^{-1}..............................................
c     ================================================================
c     ----------------------------------------------------------------
c     call mbeqa(tcpa,w3,2*kkrsz*kkrsz)
      w3 = tcpa
      call mtxinv(w3,w4,det,kkrsz,v1,v2,v3,iprint,istop)
c     ----------------------------------------------------------------
c
      if(iprint.ge.1) then
         write(6,'('' mcpait:: t_c(old) ::'')')
         call wrtmtx(tcpa,kkrsz,istop)
         write(6,'('' mcpait:: m_c(old) ::'')')
         call wrtmtx(w4,kkrsz,istop)
      endif
c
      do ic=1,komp
c        =============================================================
c        w2 => t_{\alpha}^{-1}........................................
c        =============================================================
c        -------------------------------------------------------------
c        call mbeqa(tab(1,ic),w3,2*kkrsz*kkrsz)
         w3 = tab(:,ic)
         call mtxinv(w3,w2,det,kkrsz,v1,v2,v3,iprint,istop)
c        -------------------------------------------------------------
c
c        =============================================================
c        w1 => t_c^{-1} - t_{\alpha}^{-1}..........................
c        =============================================================
c        -------------------------------------------------------------
         call madd(w4,onem,w2,w1,kkrsz,iprint,istop)
c        -------------------------------------------------------------
c
c        =============================================================
c        w2 => [t_c^{-1} - t_{\alpha}^{-1}]^{-1}...................
c        =============================================================
c        -------------------------------------------------------------
c        call mbeqa(w1,w3,2*kkrsz*kkrsz)
         w3 = w1
         call mtxinv(w3,w2,det,kkrsz,v1,v2,v3,iprint,istop)
c        -------------------------------------------------------------
c
c        =============================================================
c        w1 => [[t_c^{-1} - t_{\alpha}^{-1}]^{-1} - \tau_c]...........
c        =============================================================
c        -------------------------------------------------------------
         call madd(w2,onem,diamtx,w1,kkrsz,iprint,istop)
c        -------------------------------------------------------------
c
c        =============================================================
c        w1 => [[t_c^{-1} - t_{\alpha}^{-1}]^{-1} - \tau_c]^{-1}...
c        w1 =>  X_\alpha...........................................
c        =============================================================
c        -------------------------------------------------------------
c        call mbeqa(w1,w3,2*kkrsz*kkrsz)
         w3 = w1
         call mtxinv(w3,xab(1,ic),det,kkrsz,v1,v2,v3,iprint,istop)
c        -------------------------------------------------------------
      enddo
c
c     ================================================================
c     w1 => X_c = \Sum_\alpha{C_\alpha X_\alpha}......................
c     ================================================================
c     ----------------------------------------------------------------
      call mcav(atcon,xab,w1,kkrsz,komp,mxcomp,iprint,istop)
c     ----------------------------------------------------------------
c
c     ================================================================
c     w2 => X_{c}^{-1}................................................
c     ================================================================
c     ----------------------------------------------------------------
c     call mbeqa(w1,w3,2*kkrsz*kkrsz)
      w3 = w1
      call mtxinv(w3,w2,det,kkrsz,v1,v2,v3,iprint,istop)
c     ----------------------------------------------------------------
c
c     ================================================================
c     w1 => X_{c}^{-1} + \tau_c.......................................
c     ================================================================
c     ----------------------------------------------------------------
      call madd(w2,one,diamtx,w1,kkrsz,iprint,istop)
c     ----------------------------------------------------------------
c
c     ================================================================
c     w2 => [X_{c}^{-1} + \tau_c]^{-1}................................
c     ================================================================
c     ----------------------------------------------------------------
c     call mbeqa(w1,w3,2*kkrsz*kkrsz)
      w3 = w1
      call mtxinv(w3,w2,det,kkrsz,v1,v2,v3,iprint,istop)
c     ----------------------------------------------------------------
c
c     ================================================================
c     w1 =>  t_{c}^{-1} - [X_{c}^{-1} + \tau_c]^{-1}..................
c     w1 =>  t_{c}^{-1}(new).......................................... 
c     ================================================================
c     ----------------------------------------------------------------
      call madd(w4,onem,w2,w1,kkrsz,iprint,istop)
c     ----------------------------------------------------------------
c
c     ================================================================
c     w2 =>  t_{c}(new)............................................... 
c     ================================================================
c     ----------------------------------------------------------------
c     call mbeqa(w1,w3,2*kkrsz*kkrsz)
      w3 = w1
      call mtxinv(w3,w2,det,kkrsz,v1,v2,v3,iprint,istop)
c     ----------------------------------------------------------------
c
      if(iprint.ge.1) then
         write(6,'('' mcpait:: t_c(new) ::'')')
         call wrtmtx(w2,kkrsz,istop)
      endif
c
c     ================================================================
c     w1 => \delta t_{c} = t_{c}(new) - t_{c}(old)....................
c     ================================================================
c     ----------------------------------------------------------------
      call madd(w2,onem,tcpa,w1,kkrsz,iprint,istop)
c     ----------------------------------------------------------------
c
c     ================================================================
c     check for convergence...........................................
c     ================================================================
      do i=1,kkrsz*kkrsz
         errcpa=max(errcpa,abs(w1(i)))
         tcpa(i)=w2(i)
      enddo
c
c     ****************************************************************
      if(istop.eq.sname) then
        call fstop(sname)
      endif
c
      return
      end
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mbeqa(a,b,n)
c     ================================================================
c
      integer    i,n
c
      real*8     a(n)
      real*8     b(n)
c
      do i=1,n
         b(i)=a(i)
      enddo
c
      return
      end
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine madd(amt,const,bmt,cmt,kkrsz,iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  istop*10,sname*10
c
      integer    kkrsz
      integer    i
      integer    n
      integer    iprint
c
      real*8     const
c
      complex*16 amt(kkrsz*kkrsz)
      complex*16 bmt(kkrsz*kkrsz)
      complex*16 cmt(kkrsz*kkrsz)
c
      parameter (sname='madd')
c
c     ****************************************************************
c     calculates C=A+cB     : A,B,C are matrices and c is a constant
c     used in iteration of CPA equation...........
c     bg & gms [Nov 1991].....
c     ****************************************************************
c
      n=kkrsz*kkrsz
      do i=1,n
         cmt(i)=amt(i)+const*bmt(i)
      enddo
c
c     ================================================================
c     print if needed.................................................
c     ================================================================
      if (iprint.ge.1) then
         write(6,'('' madd:: kkrsz,const :'',i5,d12.4)') kkrsz,const
         call wrtmtx(cmt,kkrsz,istop)
      endif
c
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mcav(atcon,xab,xav,kkrsz,komp,mxcomp,
     >                iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  istop*10
      character  sname*10
c
      integer    mxcomp
      integer    kkrsz
      integer    komp
      integer    ic
      integer    i
      integer    iprint
c
      real*8     atcon(mxcomp)
c
      complex*16 xab(kkrsz*kkrsz,mxcomp)
      complex*16 xav(kkrsz*kkrsz)
c
      parameter  (sname='mcav')
c
c     ****************************************************************
c     calculates concentration average of xab and puts in xav.........
c     used in iteration of CPA equation...........
c     bg & gms [Nov 1991].....
c     ****************************************************************
c
c     ================================================================
c     zeroout array that average is stored in.........................
c     ================================================================
c     ----------------------------------------------------------------
      call zeroout(xav,2*kkrsz*kkrsz)
c     ----------------------------------------------------------------
c
c     ================================================================
c     average over components.........................................
c     ================================================================
      do ic=1,komp
         do i=1,kkrsz*kkrsz
            xav(i)=xav(i)+atcon(ic)*xab(i,ic)
         enddo
      enddo
c     ================================================================
c     print if needed.................................................
      if (iprint.ge.1) then
         write(6,'('' mcav:: kkrsz,komp :'',2i5)') kkrsz,komp
         call wrtmtx(xav,kkrsz,istop)
      endif
c
c     ================================================================
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mtxinv(amt,bmt,det,nrmat,
     >                  td,ad,bd,
     >                  iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  istop*10
      character  sname*10
c
      integer    i,j,k,ik,jk
      integer    nrmat
      integer    iprint
c
      complex*16 amt(nrmat,nrmat)
      complex*16 bmt(nrmat,nrmat)
      complex*16 td(nrmat)
      complex*16 ad(nrmat)
      complex*16 bd(nrmat)
      complex*16 det
      complex*16 czero
      complex*16 cone
      complex*16 amtinv
c
      parameter  (czero=(0.0d0,0.0d0))
      parameter  (cone=(1.0d0,0.0d0))
      parameter  (sname='mtxinv')
c
c     ****************************************************************
c     full matrix storage version...........bg & gms [Nov 1991].
c     version for scalar-relativistic code..gms [Apr 1992].
c
c     output:     bmt = amt**(-1)
c                 det = Det|amt|
c     ****************************************************************
c
c     ================================================================
c     set up unit matrix..............................................
c     ================================================================
c
c     ----------------------------------------------------------------
      call zeroout(bmt,2*nrmat*nrmat)
c     ----------------------------------------------------------------
c
      do i=1,nrmat
         bmt(i,i)=cone
      enddo
c
c     ================================================================
c     transform amt and bmt...........................................
c     ================================================================
c
      do i=1,nrmat
         amtinv=cone/amt(i,i)
         do j=1,nrmat
            td(j)=amtinv*amt(j,i)
         enddo
         td(i)=czero
         do k=1,nrmat
            ad(k)=amt(i,k)
            bd(k)=bmt(i,k)
         enddo
         do k=i,nrmat
            do j=1,nrmat
               amt(j,k)=amt(j,k)-(td(j)*ad(k))
            enddo
         enddo
         do k=1,i
            do j=1,nrmat
               bmt(j,k)=bmt(j,k)-(td(j)*bd(k))
            enddo
         enddo
      enddo
c
c     ================================================================
c     calculate determinant ..........................................
c     ================================================================
c
      det=cone
      do ik=1,nrmat
         det=det*amt(ik,ik)
      enddo
c
      do jk=1,nrmat
         do ik=1,nrmat
            bmt(ik,jk)=bmt(ik,jk)/amt(ik,ik)
         end do
      end do
c
c     ================================================================
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wrtmtx(x,n,istop)
c     ================================================================
c
      implicit   none
c
      character  sname*10
      character  istop*10
c
      integer    n
      integer    i
      integer    j
c
      real*8     tol
c
      complex*16 x(n,n)
c
      parameter  (sname='wrtmtx')
      parameter  (tol=0.1**6)
c
c     ****************************************************************
c     writes out the non-zero elements of a NxN complex matrix........
c     ****************************************************************
c
      do j=1,n
         do i=1,n
            if(abs(x(i,j)).gt.tol) then
               write(6,'(2i4,2d15.8)') i,j,x(i,j)
            endif
         enddo
      enddo
c
      if (istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fstop(name)
c     =================================================================
c
c     this routine will stop the program when an error has occured
c     in a subroutine
c
      character (len=*), intent(in) :: name
c
      if(name .ne. 'null' .and. name .ne. '     ') then
         write(6,'(///)')
         write(6,'(3x,''***************************************'')')
         write(6,'(3x,''* fstop called with index '',a10,t42,''*'')')
     >         name
         write(6,'(3x,''***************************************'')')
      endif
c
c     call kill(0,9)
      stop
      end
