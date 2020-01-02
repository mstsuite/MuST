c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ord3v(v3out_1,v3out_2,v3out_3,vsqout,nv3,v3in,vsqin,
     >                 istop)
c     ================================================================
c
      implicit  none
c
      character sname*32
      character istop*32
c
      integer   nv3
      integer   nv
      integer   nvm
c
      real*8    v3out_1(nv3+1)
      real*8    v3out_2(nv3+1)
      real*8    v3out_3(nv3+1)
      real*8    vsqout(nv3+1)
      real*8    v3in(3)
      real*8    vsqin
c
      parameter (sname='ord3v')
c
c     ****************************************************************
c     inserts a vector in a list of vectors such that they are in 
c     order of increasing length......................................
c     ****************************************************************
c
c     ================================================================
      do nv=1,nv3
         if(vsqout(nv).ge.vsqin) then
            do nvm=nv3,nv,-1
               vsqout(nvm+1)=vsqout(nvm)
            enddo
            do nvm=nv3,nv,-1
               v3out_1(nvm+1) = v3out_1(nvm)
            enddo
            do nvm=nv3,nv,-1
               v3out_2(nvm+1) = v3out_2(nvm)
            enddo
            do nvm=nv3,nv,-1
               v3out_3(nvm+1) = v3out_3(nvm)
            enddo
            vsqout(nv)=vsqin
            v3out_1(nv) = v3in(1)
            v3out_2(nv) = v3in(2)
            v3out_3(nv) = v3in(3)
            nv3=nv3+1
            return
         endif
      enddo
      vsqout(nv3+1)=vsqin
      v3out_1(nv3+1) = v3in(1)
      v3out_2(nv3+1) = v3in(2)
      v3out_3(nv3+1) = v3in(3)
      nv3=nv3+1
c     ================================================================
c     return condition................................................
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
