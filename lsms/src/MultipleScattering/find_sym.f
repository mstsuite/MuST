c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_sym(delta,sym_ops,nlen,idcol,ipvt,iprint)
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ********************************************************************
c     find the symmetry relations between the columns of delta
c     The symmetry matrices are stored in sym_ops
c     ********************************************************************
c
      implicit none
c
      integer nlen
      integer ipvt(nlen),idcol(nlen)
      integer i,j,k,ioff,iprint,k0,itmp
      integer ks,ip(2)
      complex*16 delta(nlen,nlen)
c     complex*16 sym_ops(nlen,nlen,nlen-1)
      complex*16 sym_ops(nlen,nlen,(nlen-1)/2)
      complex*16 a,b,v1,v2,s
      real*8 c
      complex*16 zone,zero
      parameter (zone=(1.d0,0.d0))
      parameter (zero=(0.d0,0.d0))
      real*8 sum1,sum2,absq
      complex*16 x
      absq(x)=dreal(x*conjg(x))

c For each column, check to see if the norm is the same as a previous
c column
      ioff=1
      idcol(1)=1
      do i=2,nlen
c Symmetry is turned off by the following goto statement
c       goto 11
c if ioff exceeds the dimension of sym_op we forget about symmetry
	if(ioff.gt.(nlen-1)/2) goto 11
      do j=1,i-1
      sum1=0.d0
      sum2=0.d0
      do k=1,nlen
        sum1=sum1+absq(delta(k,j))
        sum2=sum2+absq(delta(k,i))
      enddo
      if(abs(sum1-sum2).lt.1.d-10*sum1) then
c If the norms are the same they are possibly equivelent
	do k=1,nlen
	  ipvt(k)=k
	  sym_ops(k,1,ioff)=0.d0
	enddo
	ks=0
	do k=1,nlen
c Don't do any zero's since they'll mess up the matching
	  if(absq(delta(k,i)).lt.1.d-20*sum1) goto 20
	  do k0=1,nlen
	      if(sym_ops(ipvt(k0),1,ioff).eq.0.d0) then
	      if(absq(delta(k,i)-delta(ipvt(k0),j)).lt.1.d-20*sum1) then
		itmp=ipvt(k)
		ipvt(k)=ipvt(k0)
		ipvt(k0)=itmp
		sym_ops(ipvt(k),1,ioff)=1.d0
		goto 20
	      endif
	      endif
	  enddo  ! k0
c We haven't found any match
	  ks=ks+1
c If there are more than two elements not matched we give up and quit
	  if(ks.gt.2) goto 11
	  ip(ks)=k
 20       continue
	enddo  ! k
c Match the ks elements with the nonzero elements in jth col
	do k=1,ks
	  if(absq(delta(ipvt(ip(k)),j)).lt.1.d-20*sum1) then
	    do k0=1,nlen
	      if(sym_ops(ipvt(k0),1,ioff).eq.0.d0.and.
     &          absq(delta(ipvt(k0),j)).gt.1.d-20*sum1) then
		itmp=ipvt(ip(k))
		ipvt(ip(k))=ipvt(k0)
		ipvt(k0)=itmp
	      endif
	    enddo  ! k0
	  endif
	enddo  ! ks
c Now match the zero's
        do k=1,nlen
	  if(absq(delta(k,i)).lt.1.d-20*sum1) then
	    if(absq(delta(ipvt(k),j)).lt.1.d-20*sum1) then
c The zero's match
	      sym_ops(ipvt(k),1,ioff)=1.d0
	    else
c The zero's don't match
	      ks=ks+1
	      if(ks.gt.2) goto 11
	      ip(ks)=k
	    endif
	  endif
	enddo  ! k

	call zeroout(sym_ops(1,1,ioff),2*nlen*nlen)
	do k=1,nlen
	  sym_ops(ipvt(k),k,ioff)=1.d0
	enddo

c now find the rotation matrices
	if(ks.eq.2) then
	  k=ipvt(ip(1))
	  k0=ipvt(ip(2))
	  v2=delta(k,j)
	  b=delta(k0,j)
	  call zrotg(v2,b,c,s)
	  sym_ops(ip(1),k,ioff)=c
	  sym_ops(ip(2),k,ioff)=-conjg(s)
	  sym_ops(ip(1),k0,ioff)=s
	  sym_ops(ip(2),k0,ioff)=c

	  v1=delta(ip(1),i)
	  b=delta(ip(2),i)
	  call zrotg(v1,b,c,s)
	  if(absq(v1+v2).lt.1.d-20*sum1) then
	    c=-c
	    s=-s
	  elseif(absq(v1-v2).gt.1.d-20*sum1) then
	    goto 11
c           write(6,'(''find_sym:: A bug in find_sym!'')')
c           call fstop('find_sym')
	  endif
	  a=c*sym_ops(ip(1),k,ioff)-conjg(s)*sym_ops(ip(2),k,ioff)
	  b=s*sym_ops(ip(1),k,ioff)+c*sym_ops(ip(2),k,ioff)
	  sym_ops(ip(1),k,ioff)=a
	  sym_ops(ip(2),k,ioff)=b
	  a=c*sym_ops(ip(1),k0,ioff)-conjg(s)*sym_ops(ip(2),k0,ioff)
	  b=s*sym_ops(ip(1),k0,ioff)+c*sym_ops(ip(2),k0,ioff)
	  sym_ops(ip(1),k0,ioff)=a
	  sym_ops(ip(2),k0,ioff)=b
        else if(ks.eq.1) then
	  goto 11
c         ks cannot be 1
c         write(6,'(''FIND_SYM::Code doesnt work for the symmetry'')')
c         write(6,'(''FIND_SYM::ks is 1'')')
c         call fstop('find_sym')
	endif  ! ks.eq.2

	      if(iprint.ge.0) then
		write(6,'(''col '',i3,'' is equiv to col '',i3)') i,j
c        do k=1,nlen
c          if(absq(delta(k,j)).gt.1.d-20.or.absq(delta(k,i)).gt.1.d-20)
c    &       then
c          write(6,'(i3,4e15.6)') k,delta(k,j),delta(k,i)
c          endif
c        enddo
c	 call wrtmtx(sym_ops(1,1,ioff),nlen,' ')
	      endif
	      idcol(i)=j
	      ioff=ioff+1
	      goto 10
      endif
      enddo
 11     idcol(i)=i
 10   continue
      enddo

      return
      end
