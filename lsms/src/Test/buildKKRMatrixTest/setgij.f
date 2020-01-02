      subroutine setgij(gij, bgij,kkr1, kkr1_ns,kkr2,kkr2_ns,
     >                  n_spin_cant,nrel_rel,psq,ce)
c     ==============================================================
c
      implicit   none
c
      character  istop*32
      character  sname*32
      integer    im,in
      integer    i,j
      integer    is
      integer    kkr1
      integer    kkr1_ns
      integer    kkr2
      integer    kkr2_ns
      integer    n_spin_cant
      integer    nrel_rel
c
      complex*16 gij(kkr1,kkr2)
      complex*16 bgij(kkr1_ns,kkr2_ns)
      complex*16 fac
      complex*16 psq
      complex*16 ce

      parameter (sname='setgij')
c
      call zeroout(bgij,2*kkr1_ns*kkr2_ns)
c
      if(nrel_rel.eq.0) then
        do is=1,n_spin_cant
       	    im=(is-1)*kkr1
       	    in=(is-1)*kkr2
            do i=1,kkr1
            do j=1,kkr2
                 bgij(im+i,in+j)=gij(i,j)
            end do
            end do
        end do
      else
!       if(kkr1.ne.kkr2) then
!         write(6,*) 'SETGIJ: nonsquare rel. Gij not yet supported'
!         call fstop(sname)
!       end if
!       call relmtrx(gij,bgij,lmax)

         write(6,*) 'SETGIJ: relativity currently not supported'
         call fstop(sname)
c$$$        call relmtrx(gij,bgij,kkr1,kkr2)
c$$$        fac=psq/ce
c$$$        do i=1,kkr1_ns
c$$$          do j=1,kkr2_ns
c$$$            bgij(i,j)=fac*bgij(i,j)
c$$$          end do
c$$$        end do
      end if
c
      return
      end
