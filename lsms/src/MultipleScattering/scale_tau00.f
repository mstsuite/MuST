c
      subroutine scale_tau00(tau00,kkrsz1,kkrsz2,lofk,n_spin_cant,
     &                       kappa_rmt)
      implicit none

      integer kkrsz1,kkrsz2,n_spin_cant
      integer lofk(kkrsz1+kkrsz2)

      complex*16 tau00(kkrsz1,n_spin_cant,kkrsz2,n_spin_cant)

      integer i
      integer j
      integer ispin
      integer jspin

      complex*16 kappa_rmt
      complex*16 factor(0:100)

      factor(0)=1.d0
      do i=1,lofk(max(kkrsz1,kkrsz2))
        factor(i)=factor(i-1)*(2.d0*i+1.d0)/kappa_rmt
      enddo

      do ispin=1,n_spin_cant
        do i=1,kkrsz2
	  do jspin=1,n_spin_cant
            do j=1,kkrsz1
	      tau00(j,jspin,i,ispin)=tau00(j,jspin,i,ispin)
     &                               *factor(lofk(i))*factor(lofk(j))
	    enddo
          enddo
        enddo
      enddo

      return
      end
