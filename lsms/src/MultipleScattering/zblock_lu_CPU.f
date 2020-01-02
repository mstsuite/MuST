
      subroutine zblock_lu(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
c does a partitioning of a matrix and a partial inversion to
c get the upper subblock of the inverse
c   a -- complex*16 of (lda,size) the matrix to be inverted
c        where size is the sum of the all elements of blk_sz
c        on return, the upper subblock of a is untouched, and
c        postprocessing is needed to obtain the inverse
c   blk_sz -- integer of (nblk), the size of each subblock of matrix a
c   nblk -- number of blocks (number of elements in blk_sz)
c   ipvt -- integer of (mp), work array
c   idcol -- integer of (blk_sz(1)), if idcol(i)=idcol(j), then
c            the two columns are equivalent by symmetry, and only
c            the min(i,j) th column of the inverse is calculated.
c   k -- returns actual number of columns in the calculated inverse
      implicit none
      integer lda,na,mp,nblk
      integer ipvt(mp),blk_sz(*)
      integer i,j,k,ioff,m,joff,n,iblk,info
      integer idcol(blk_sz(1))
      complex*16 a(lda,*)
      complex*16 cone,cmone
      parameter (cone=(1.d0,0.d0))
      parameter (cmone=(-1.d0,0.d0))

      na=0
      do i=1,abs(nblk)
        na=na+blk_sz(i)
      enddo

!       write(6,*) 'idcol(1)=',idcol(1)
      if(idcol(1).eq.0) then
        k=1
      else
c eliminate columns that are equiv due to symmetry
      k=blk_sz(1)+1
      do i=blk_sz(1),1,-1
          if(idcol(1).eq.0.or.idcol(i).eq.i) then
          k=k-1
          if(k.ne.i) then
          call zcopy(na-blk_sz(1),a(blk_sz(1)+1,i),1,a(blk_sz(1)+1,k),1)
          endif
        endif
      enddo
      endif

      if(nblk.gt.0) then
c Do block LU
      n=blk_sz(nblk)
      joff=na-n
      do iblk=nblk,2,-1
      m=n
      ioff=joff
      n=blk_sz(iblk-1)
      joff=joff-n
c invert the diagonal blk_sz(iblk) x blk_sz(iblk) block
      call zgetrf(m,m,a(ioff+1,ioff+1),lda,ipvt,info)
      if(info.ne.0) then
        write(*,*) "zgetrf info=",info," ioff=",ioff
      end if
c calculate the inverse of above multiplying the row block
c blk_sz(iblk) x ioff
      call zgetrs('n',m,ioff,a(ioff+1,ioff+1),lda,ipvt,
     &     a(ioff+1,1),lda,info)
      if(info.ne.0) then
         write(*,*) "zgetrs info=",info," ioff=",ioff
      end if
      if(iblk.gt.2) then
      call zgemm('n','n',n,ioff-k+1,na-ioff,cmone,a(joff+1,ioff+1),lda,
     &     a(ioff+1,k),lda,cone,a(joff+1,k),lda)
      call zgemm('n','n',joff,n,na-ioff,cmone,a(1,ioff+1),lda,
     &     a(ioff+1,joff+1),lda,cone,a(1,joff+1),lda)
      endif
      enddo
      call zgemm('n','n',blk_sz(1),blk_sz(1)-k+1,na-blk_sz(1),cmone,
     &     a(1,blk_sz(1)+1),lda,a(blk_sz(1)+1,k),lda,cone,a,lda)
      endif ! nblk.gt.0

      k=blk_sz(1)-k+1
!     write(6,*) "k out =",k
      return
      end
