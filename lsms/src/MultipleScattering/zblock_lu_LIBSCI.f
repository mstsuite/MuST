
      integer function idx2f(i,j,ld)
      implicit none
      integer i,j,ld

      idx2f=((((j)-1)*(ld))+((i)-1))
      end function

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
      use iso_c_binding
      implicit none

!     integer*8  dev_ptr
      integer sizeof_Z,sizeof_I,ierr
      parameter (sizeof_Z=16)
      parameter (sizeof_I=8)


      integer lda,na,mp,nblk
      integer ipvt(mp),blk_sz(*)
      integer i,j,k,ioff,m,joff,n,iblk,info
      integer idcol(blk_sz(1))
      complex*16 a(lda,*)
      complex*16 cone,cmone
      parameter (cone=(1.d0,0.d0))
      parameter (cmone=(-1.d0,0.d0))
      external libsci_acc_alloc, libsci_acc_free, idx2f,
     &         libsci_acc_getmatrix
      integer  libsci_acc_alloc,libsci_acc_free
      integer  idx2f,  libsci_acc_getmatrix


      na=0
      do i=1,abs(nblk)
        na=na+blk_sz(i)
      enddo
c eliminate columns that are equiv due to symmetry
      k=blk_sz(1)+1
      do i=blk_sz(1),1,-1
        if(idcol(1).eq.0.or.idcol(i).eq.i) then
          k=k-1
          if(k.ne.i) then
            call zcopy(na-blk_sz(1),a(blk_sz(1)+1,i),1,
     &                 a(blk_sz(1)+1,k),1)
          endif
        endif
      enddo



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
c calculate the inverse of above multiplying the row block
c blk_sz(iblk) x ioff
      call zgesv_lsms( m, ioff, a(ioff+1,ioff+1),lda, ipvt,
     &            a(ioff+1,1),lda,info) 
!     call zgetrf(m,m,a(ioff+1,ioff+1),lda, ipvt,info)
!     call zgetrs('n',m,ioff,a(ioff+1,ioff+1),lda,ipvt,
!    &     a(ioff+1,1),lda,info)

!
!     Not executed in libsci_acc version     
!
      if(iblk.gt.2) then
      call zgemm('n','n',n,ioff-k+1,na-ioff,cmone,a(joff+1,ioff+1),lda,
     &     a(ioff+1,k),lda,cone,a(joff+1,k),lda)
      call zgemm('n','n',joff,n,na-ioff,cmone,a(1,ioff+1),lda,
     &     a(ioff+1,joff+1),lda,cone,a(1,joff+1),lda)
      endif
!
! 
!
      enddo

!     call zgemm_acc('n','n',blk_sz(1),blk_sz(1)-k+1,na-blk_sz(1),cmone,
!    &     dev_ptr + idx2f(1,blk_sz(1)+1,na)*sizeof_Z,na,
!    &     dev_ptr + idx2f(blk_sz(1)+1,k,na)*sizeof_Z,na,
!    &     cone, dev_ptr, na)
!     ierr = libsci_acc_getmatrix( blk_sz(1), blk_sz(1), sizeof_Z,
!    &                             dev_ptr, na, a, lda )
      call zgemm_cpu('n','n',blk_sz(1),blk_sz(1)-k+1,na-blk_sz(1),cmone,
     &     a(1,blk_sz(1)+1),lda,a(blk_sz(1)+1,k),lda,cone,a,lda)
      endif ! nblk.gt.0

      k=blk_sz(1)-k+1

      return
      end
