c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lmfacts(lmax,ndlj,ndlm,lofj,mofj,lofk,mofk)
c     ================================================================
c
      implicit   none
c
      integer    lmax
      integer    ndlj
      integer    ndlm
      integer    lofj((lmax+1)*(lmax+2)/2)
      integer    mofj((lmax+1)*(lmax+2)/2)
      integer    lofk((lmax+1)*(lmax+1))
      integer    mofk((lmax+1)*(lmax+1))
      integer    j
      integer    k
      integer    l
      integer    m
      integer    ma
c
      ndlj=(lmax+1)**2
      ndlm=(lmax+1)*(lmax+2)/2
      j=0
      k=0
      do l=0,lmax
         do m=0,l
            j=j+1
            lofj(j)=l
            mofj(j)=m
         enddo
         do m=-l,l
            k=k+1
            lofk(k)=l
            mofk(k)=m
         enddo
      enddo           
c
      return
      end
