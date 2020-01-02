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


      implicit none
      integer lda,na,mp,nblk
      integer ipvt(mp),blk_sz(*)
      integer i,j,k,ioff,m,joff,n,iblk,info
      integer idcol(blk_sz(1))
      complex*16 a(lda,*)
      complex*16 cone,cmone
      parameter (cone=(1.d0,0.d0))
      parameter (cmone=(-1.d0,0.d0))


      external idx2f
      integer idx2f
      external cublas_alloc,cublas_set_matrix,cublas_get_matrix
      external cublas_free
      integer cublas_alloc,cublas_set_matrix,cublas_get_matrix
      integer cublas_free

      integer sizeof_Z,sizeof_I
      parameter (sizeof_Z=16)
      parameter (sizeof_I=8)
      integer*8 devA, get_dev_m
!     include 'CULA_Common.h'
!      include '../Accelerator/CULA_Common.h'

      devA = get_dev_m()

      na=0
      do i=1,abs(nblk)
        na=na+blk_sz(i)
      enddo
      k=1

!     copy matrix to device
!     info = CUBLAS_ALLOC(lda*na,  sizeof_Z , devA)
!      info = cublas_set_matrix(lda, na, sizeof_Z, a, lda, devA, lda)
!      write(*,*) 'cublas_set_matrix of A :',info, na

!! #endif

      if(nblk.gt.0) then
c Do block LU
      n=blk_sz(nblk)
      joff=na-n
      do iblk=nblk,2,-1
!      write(*,*) 'iblk = ',iblk
      m=n
      ioff=joff
      n=blk_sz(iblk-1)
      joff=joff-n
c invert the diagonal blk_sz(iblk) x blk_sz(iblk) block
!      write(*,*) m,ioff,lda,idx2f(ioff+1,ioff+1,lda)
!     call zgetrf(m,m,a(ioff+1,ioff+1),lda,ipvt,info)
!      write(*,*) 'm',m,'m',m,'lda',lda
      call zgetrf_acc2(m,m,
     &       devA+idx2f(ioff+1,ioff+1,lda)*sizeof_Z,lda,
     &       IPVT,info)
!     write(*,*) 'zgetrf info=',info
c calculate the inverse of above multiplying the row block
c blk_sz(iblk) x ioff
!     call zgetrs('n',m,ioff,a(ioff+1,ioff+1),lda,ipvt,
!    &     a(ioff+1,1),lda,info)
      call zgetrs_acc2('n',m,ioff,
     &     devA+idx2f(ioff+1,ioff+1,lda)*sizeof_Z,lda,IPVT,
     &     devA+idx2f(ioff+1,1,lda)*sizeof_Z,lda,info)
      if(iblk.gt.2) then
!     call zgemm('n','n',n,ioff-k+1,na-ioff,cmone,a(joff+1,ioff+1),lda,
!    &     a(ioff+1,k),lda,cone,a(joff+1,k),lda)
      call zgemm_acc2('n','n',n,ioff-k+1,na-ioff,cmone,
     &     devA+idx2f(joff+1,ioff+1,lda)*sizeof_Z,lda,
     &     devA+idx2f(ioff+1,k,lda)*sizeof_Z,lda,cone,
     &     devA+idx2f(joff+1,k,lda)*sizeof_Z,lda)
!     call zgemm('n','n',joff,n,na-ioff,cmone,a(1,ioff+1),lda,
!    &     a(ioff+1,joff+1),lda,cone,a(1,joff+1),lda)
      call zgemm_acc2('n','n',joff,n,na-ioff,cmone,
     &     devA+idx2f(1,ioff+1,lda)*sizeof_Z,lda,
     &     devA+idx2f(ioff+1,joff+1,lda)*sizeof_Z,lda,cone,
     &     devA+idx2f(1,joff+1,lda)*sizeof_Z,lda)
!     info = cublas_get_matrix(lda, na, sizeof_Z, devA, lda, a, lda)
      endif
      enddo
!     call zgemm('n','n',blk_sz(1),blk_sz(1)-k+1,na-blk_sz(1),cmone,
!    &     a(1,blk_sz(1)+1),lda,a(blk_sz(1)+1,k),lda,cone,a,lda)
      call zgemm_acc2('n','n',blk_sz(1),blk_sz(1)-k+1,na-blk_sz(1),
     &     cmone,devA+idx2f(1,blk_sz(1)+1,lda)*sizeof_Z,lda,
     &     devA+idx2f(blk_sz(1)+1,k,lda)*sizeof_Z,lda,cone,devA,lda)
      info = cublas_get_matrix(lda,blk_sz(1),sizeof_Z,devA,lda,a,lda)
!     info = cublas_free( devA )
      endif ! nblk.gt.0

      k=blk_sz(1)-k+1
      return
      end
