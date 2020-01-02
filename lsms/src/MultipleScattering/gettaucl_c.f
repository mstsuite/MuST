c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gettaucl(mynod,n_spin_cant,nrel_rel,
     >                    lmax,kkrsz,ndlj,ndlm,
     >                    clm,cgnt,
     >                    lofk,mofk,lofj,mofj,
     >                    ilp1,illp,
     >                    prel,energy,kappa_rmt,
     >                    rsclu_1,rsclu_2,rsclu_3,nrsclu,nrmat,lmaxclu,
     >                    tmat,tmat_n,tau00,delta,wbig,ialg,
     >                    iswtch,ied,nume,idcol,sym_ops,
     >                    nbortyp,nsend,msg_send,nsnd_tot,
     >                    pi4,iprint,istop)
c     ================================================================
      use MPPModule, only : sendMessage, recvMessage
      use MPPModule, only : nbsendMessage, nbrecvMessage
      use MPPModule, only : waitMessage
c
      implicit   none
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include    'atom_param.h'
      include    'param_rsp.h'
c     include    'mpif.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      character  istop*32
      character  sname*32
      logical    received(iprsclu)
c
      integer    mynod
      integer    n_spin_cant
      integer    nrel_rel
      integer    lmax
      integer    kkrsz
      integer    ndlj
      integer    ndlm
      integer    nrmat
      integer    lofk(ndlj)
      integer    mofk(ndlj)
      integer    lofj(ndlm)
      integer    mofj(ndlm)
      integer    nrsclu
      integer    nsnd_tot
      integer    lmaxclu(nrsclu)
      integer    nbortyp(nrsclu)
      integer    nsend(nsnd_tot)
      integer    msg_send(nsnd_tot)
      integer, allocatable ::  msgids(:)
      integer    msgidr
      integer    iprint
      integer    info
      integer    ipvt(ipcmat*ip_spin_cant)
      integer    blk_sz(iprsclu*ip_spin_cant),nblk
      integer    min_sz
c
      integer    kkrsz_ns
      integer    nrmat_ns
      integer    mtxsize
      integer    msgsize
      integer    is
      integer    i
      integer    im
      integer    js
      integer    jsm
      integer    j
      integer    jm
      integer    ii
c*********************************TMP****************
      integer    ir
      integer    l
c*********************************TMP****************
      integer    ir1
      integer    ir2
      integer    kkr1
      integer    kkr1_ns
      integer    kkr2
      integer    kkr2_ns
      integer    kkr12
      integer    nrst
      integer    ncst
      integer    ialg
      integer    ied
      integer    nume
      integer    iswtch(nume)
      integer    nm
      integer    nmsg_sent
      integer    nrsclu_max
      integer    idcol(kkrsz*n_spin_cant)
      integer, allocatable :: iwork(:)
c
      real*8     clm(ndlm)
      real*8     cgnt((lmax+1)*kkrsz*kkrsz)
      real*8     rsclu_1(nrsclu)
      real*8     rsclu_2(nrsclu)
      real*8     rsclu_3(nrsclu)
      real*8     pi4
      real*8     one
      real*8     rij(3)
      real*8     plm(ipdlm)
      real*8     sinmp(2*iplmax+1)
      real*8     cosmp(2*iplmax+1)
      real*8     time_direct
      real*8     time_to_wait
      real*8     t0
      real*8     clock_time
      integer    block_size
      integer    block_dir
      real*8     time_old,time_new
      real*8, allocatable :: rwork(:)
c
c     complex*16 sym_ops(kkrsz*kkrsz*n_spin_cant*n_spin_cant*
c    &          (kkrsz*n_spin_cant-lmax-1))
      complex*16 sym_ops(kkrsz*kkrsz*n_spin_cant*n_spin_cant*
     &          (kkrsz*n_spin_cant-1)/2)
      complex*16 ilp1(0:2*lmax)
      complex*16 illp(kkrsz*kkrsz)
      complex*16 prel
      complex*16 energy
      complex*16 kappa_rmt
      complex*16 tmat(kkrsz*kkrsz*n_spin_cant*n_spin_cant)
      complex*16 tmat_n(kkrsz*kkrsz*n_spin_cant*n_spin_cant)
      complex*16 tau00(kkrsz*kkrsz*n_spin_cant*n_spin_cant)
      complex*16 delta(kkrsz*kkrsz*n_spin_cant*n_spin_cant)
      complex*16 wbig(nrmat*nrmat*n_spin_cant*n_spin_cant)
      complex*16, allocatable ::  work1(:)
      complex*16, allocatable, target :: vbig(:)
      complex*16, pointer :: gij(:)
      complex*16, pointer :: hfn(:)
      complex*16, pointer :: dlm(:)
c     **************************************************
c
      complex*16 cone
      complex*16 cmone
      complex*16 czero
      complex*16 sqrtm1
c
      parameter (one=1.0d0)
      parameter (cone=(1.0d0,0.0d0))
      parameter (cmone=(-1.0d0,0.0d0))
      parameter (czero=(0.0d0,0.0d0))
      parameter (sqrtm1=(0.0d0,1.0d0))
      parameter (sname='gettaucl')
      data block_size/3/,block_dir/1/
      data time_old/0.d0/
      data time_direct/0.d0/
      save block_size,block_dir,time_old,time_direct
c
c     ****************************************************************
c     full matrix storage version.........bg & gms [Nov 1991].........
c     ****************************************************************
c
      if(iprint.ge.1) then
         write(6,'('' gettaucl:: mynod ='',i3)') mynod
         write(6,'('' gettaucl:: lmax  ='',i3)') lmax
         write(6,'('' gettaucl:: ndlj  ='',i3)') ndlj
         write(6,'('' gettaucl:: ndlm  ='',i3)') ndlm
         write(6,'('' gettaucl:: kkrsz ='',i3)') kkrsz
         write(6,'('' gettaucl:: nrsclu='',i3)') nrsclu
         do ir1=1,nrsclu
            write(6,'(''gettaucl:: mynod,rsclu(i,ir1):'',2i5,3d15.5)') 
     >            mynod,ir1,rsclu_1(ir1),rsclu_2(ir1),rsclu_3(ir1)
         enddo
         call flush(6)
      endif
c
c     call MPI_BARRIER( MPI_COMM_WORLD, info )
c
      kkrsz_ns=kkrsz*n_spin_cant
      nrmat_ns=nrmat*n_spin_cant
      allocate(msgids(nsnd_tot))
c
c     ================================================================
c     setup unit matrix...............................................
c     ----------------------------------------------------------------
      call cmtruni(wbig,nrmat_ns)
c     ----------------------------------------------------------------
c
      mtxsize=kkrsz_ns*kkrsz_ns
      msgsize=16*mtxsize
c     ================================================================
      allocate( iwork(1:kkrsz*n_spin_cant*nrsclu) )
      allocate( rwork(1:kkrsz*n_spin_cant*nrsclu) )
      allocate( work1(1:kkrsz*n_spin_cant*nrsclu) )
      j=max(nrmat*(kkrsz*n_spin_cant*6+6)*n_spin_cant,
     &      nrmat*kkrsz*n_spin_cant*n_spin_cant+kkrsz*kkrsz+
     &      2*lmax+1+ndlj)
      allocate( vbig(1:j) )
      gij=>vbig(nrmat*kkrsz*n_spin_cant*n_spin_cant+1:
     >          nrmat*kkrsz*n_spin_cant*n_spin_cant+kkrsz*kkrsz)
      hfn=>vbig(nrmat*kkrsz*n_spin_cant*n_spin_cant+kkrsz*kkrsz+1:
     >          nrmat*kkrsz*n_spin_cant*n_spin_cant+kkrsz*kkrsz+
     >          2*lmax+1)
      dlm=>vbig(nrmat*kkrsz*n_spin_cant*n_spin_cant+kkrsz*kkrsz+
     >          2*lmax+1+1:
     >          nrmat*kkrsz*n_spin_cant*n_spin_cant+kkrsz*kkrsz+
     >          2*lmax+1+ndlj)
c     ================================================================
c
c     call scale_tau00(tmat,kkrsz,kkrsz,lofk,n_spin_cant,kappa_rmt)
c     ================================================================
c     for each Rij coordinate to set up Tau**(-1)=[1-tG] matrix.......
c     ================================================================
      nrst=0
      nmsg_sent=0
      do ir1=1,nrsclu
	received(ir1)=.false.
      enddo
c     ----------------------------------------------------------------
      call zcopy(mtxsize,tmat,1,vbig,1)
c     ----------------------------------------------------------------
      if(ialg.ge.10) then
      nblk=2
      blk_sz(1)=kkrsz_ns
      blk_sz(2)=nrmat_ns-kkrsz_ns
      else
c The minimum block size is 16
      blk_sz(1)=kkrsz_ns
        if(ialg.eq.1) then
	 nblk=nrsclu
         do ir1=2,nrsclu
         kkr1_ns=(lmaxclu(ir1)+1)**2*n_spin_cant
	 blk_sz(ir1)=kkr1_ns
	 enddo
	else
        if(block_size.eq.0) block_size=1
        nblk=min(block_size,(nrmat_ns-blk_sz(1))/16)+1
c       write(6,'('' GETTAUCL:: mynod,nblk,nrmat_ns,block_size''
c    >  ,5i5)')                 mynod,nblk,nrmat_ns,block_size
	min_sz=(nrmat_ns-blk_sz(1))/(nblk-1)
	 do i=2,nblk-1
	 blk_sz(i)=min_sz
	 enddo
	 blk_sz(nblk)=nrmat_ns-blk_sz(1)-(nblk-2)*min_sz
	endif
      endif
      do ir1=1,nrsclu
         kkr1=(lmaxclu(ir1)+1)**2
         kkr1_ns=kkr1*n_spin_cant
         ncst=0
c        =============================================================
c        download the t-matrix from the temperory area (vbig) to tmat_n
c        and re-arrange tmat_n so that tmat_n is in the form 
c           tmat_n(kkr1*n_spin_cant,kkr1*n_spin_cant)
c        =============================================================
         im=0
	 do js=1,n_spin_cant
	    jsm=kkrsz*kkrsz_ns*(js-1)
	    do j=1,kkr1
	 if(.not.received(ir1)) then
	       do is=1,n_spin_cant
c                 ----------------------------------------------------
	          jm=jsm+kkrsz_ns*(j-1)+kkrsz*(is-1)
		  call zcopy(kkr1,vbig(jm+1),1,tmat_n(im+1),1)
c                 ----------------------------------------------------
		  im=im+kkr1
	       enddo
	 else
c the t-matrix has been stored in an appropriate place in wbig
c                 ----------------------------------------------------
	          jm=nrst+nrmat_ns*(kkr1*(js-1)+j-1)
		  call zcopy(kkr1_ns,wbig(jm+1),1,tmat_n(im+1),1)
c restore wbig for use later
		  call zeroout(wbig(jm+1),2*kkr1_ns)
c                 ----------------------------------------------------
		  im=im+kkr1_ns
	 endif
	    enddo
	 enddo
c
         if(ir1.lt.nrsclu) then
c           ==========================================================
c           post send command to send t-matrix to other nodes.........
c           ==========================================================
            nm=0
c           write(6,'(''GETTAUCL: nsnd_tot '',i5)') nsnd_tot
            call flush(6)
            do i=1,nsnd_tot
              if(msg_send(i).eq.ir1+1) then
                  nm=nm+1
c                 write(6,'('' SENDING: ir1,i,nm,msg_send,nsend '',
c    >            5i5)') ir1,i,nm,msg_send(i),nsend(i)
c                 write(6,'('' MPI_COMM_WORLD '',i10)') MPI_COMM_WORLD
c                 do ii = 1,mtxsize,32+1
c                   write(6,'('' i,tmat '',i5,2d14.6)') ii,tmat(ii)
c                 enddo
                  call flush(6)
c                 ----------------------------------------------------
c                 call sendMessage(tmat,mtxsize,msg_send(i),nsend(i))
c                 ----------------------------------------------------
                  msgids(nm)=nbsendMessage(tmat,mtxsize,msg_send(i),
     &                                     nsend(i))
c                 ----------------------------------------------------
c                 call MPI_send(tmat,mtxsize,MPI_DOUBLE_COMPLEX,
c    >            nsend(i),msg_send(i), MPI_COMM_WORLD,info)
c                 ----------------------------------------------------
c                 write(6,'('' EXITING SEND: ir1,i,nm '',3i5)') 
c    >            ir1,i,nm
c                 call flush(6)
              endif
            enddo
         else
            t0=clock_time()
         endif
c
c        =============================================================
         do ir2=1,nrsclu
c           write(6,'('' IR1,IR2 '',2i5)') ir1,ir2
c           call flush(6)
            kkr2=(lmaxclu(ir2)+1)**2
            kkr2_ns=kkr2*n_spin_cant
            if (ir1.ne.ir2) then
               rij(1)=rsclu_1(ir1)-rsclu_1(ir2)
               rij(2)=rsclu_2(ir1)-rsclu_2(ir2)
               rij(3)=rsclu_3(ir1)-rsclu_3(ir2)
c              =======================================================
c              g(Rij) calculation.....................................
c              -------------------------------------------------------
               call makegij(lmaxclu(ir1),kkr1,lmaxclu(ir2),kkr2,
     >                      lmax,kkrsz,ndlj,ndlm,
     >                      prel,rij,sinmp,cosmp,
     >                      clm,plm,cgnt,lofk,mofk,
     >                      ilp1,illp,
     >                      hfn,dlm,gij,
     >                      pi4,iprint,istop)
c              -------------------------------------------------------
c              call inv_scale_tau00(gij,kkr1,kkr2,lofk,1,
c    &                              kappa_rmt)
c              -------------------------------------------------------
	       do is=1,n_spin_cant
c                 ====================================================
c                    s,s'  ij
c              form t   * g  (e ) for current Rij and store the result
c                   -i    -    s
c
c              into wbig => 1-t*G .........................
c              -------------------------------------------------------
	          im=(is-1)*kkr1_ns
                  call zgemm('n','n',kkr1_ns,kkr2,kkr1,cmone,
     >                       tmat_n(im*kkr1+1),kkr1_ns,gij,kkr1,czero,
     &                       wbig(nrmat_ns*(ncst+(is-1)*kkr2)+nrst+1),
     &                       nrmat_ns)
c                 ----------------------------------------------------
               enddo
            endif
	 if(nbortyp(ir2).eq.nbortyp(ir1).and.ir1.gt.1.and.
     &     ir2.gt.ir1.and..not.received(ir2)) then
	 received(ir2)=.true.
         im=0
               do js=1,n_spin_cant
	    do j=1,kkr2
c                 ----------------------------------------------------
c store the t-matrix in the first column of wbig corresponding to atom ir2.
               do is=1,n_spin_cant
	          jm=ncst+nrmat_ns*(j-1+kkr2*(js-1))+kkr2*(is-1)
c                 ----------------------------------------------------
		  call zcopy(kkr2,tmat_n(im+1),1,wbig(jm+1),1)
c                 ----------------------------------------------------
		  im=im+kkr1
               enddo
	    enddo
                  im=im+kkr1_ns*(kkr1-kkr2)
               enddo
	 endif
            ncst=ncst+kkr2_ns
         enddo
         nrst=nrst+kkr1_ns
c
         if(ir1.lt.nrsclu) then
c           ==========================================================
c           post receive command to receive t-matrix from other node..
c           ==========================================================
            if(nbortyp(ir1+1).ne.mynod+1) then
c              -------------------------------------------------------
               if(.not.received(ir1+1)) then
c                write(6,'('' RECEIVING: ir1,ir1+1, nbortyp '',3i5)')
c    >           ir1,ir1+1,nbortyp(ir1+1)-1
c                call flush(6)
c                 ----------------------------------------------------
c                 call recvMessage(vbig,mtxsize,ir1+1,nbortyp(ir1+1)-1)
c                 ----------------------------------------------------
                  msgidr=nbrecvMessage(vbig,mtxsize,ir1+1,
     >            nbortyp(ir1+1)-1)
c                 ----------------------------------------------------
                  call waitMessage(msgidr)
c                 ----------------------------------------------------
               endif
c              -------------------------------------------------------
               do i=1,nm
c                -----------------------------------------------------
                  call waitMessage(msgids(i))
c                -----------------------------------------------------
               enddo
c              -------------------------------------------------------
c              write(6,'('' Exiting WAIT ON MSGIDS '')')
c              call flush(6)
c              -------------------------------------------------------
            else
c              -------------------------------------------------------
               call zcopy(mtxsize,tmat,1,vbig,1)
c              -------------------------------------------------------
            endif
            nmsg_sent=nm+nmsg_sent
         else
            time_to_wait=clock_time()-t0
         endif
      enddo
c
      if(nmsg_sent.lt.nsnd_tot) then
         write(6,'('' nmsg_sent < nsnd_tot ='',2i5)')nmsg_sent,nsnd_tot
         nrsclu_max=0
         do i=1,nsnd_tot
            nrsclu_max=max(nrsclu_max,msg_send(i))
         enddo
         do ir1=nrsclu+1,nrsclu_max
c           ----------------------------------------------------------
            nm=0
            do i=1,nsnd_tot
               if(msg_send(i).eq.ir1) then
                  nm=nm+1
c                 ----------------------------------------------------
                  call sendMessage(tmat,mtxsize,msg_send(i),nsend(i))
c                 ----------------------------------------------------
               endif
            enddo
            nmsg_sent=nmsg_sent+nm
         enddo
         write(6,'('' Finally, nmsg_sent ='',i5)')nmsg_sent
      else if(nmsg_sent.gt.nsnd_tot) then
         write(6,'('' nmsg_sent > nsnd_tot ='',2i5)')nmsg_sent,nsnd_tot
      endif
c
c     ================================================================
      if(nrst.ne.ncst) then
         write(6,'('' GETTAUCL:: nrst <> ncst:'',2i5)')nrst,ncst
         call fstop(sname)
      else if(nrst.ne.nrmat_ns) then
         write(6,'('' GETTAUCL:: nrst <> nrmat_ns:'',2i5)')nrst,nrmat_ns
         call fstop(sname)
      endif
c     ================================================================
c     print if needed.................................................
      if(iprint.ge.3) then
c        write(6,'('' gettaucl:: 1-t*g(Rij) :: mynod ='',1i5)') mynod
c        call wrtmtx(wbig,nrmat_ns,istop)
      endif
c     ================================================================
c     invert [1-t*G] : obtain the kkrsz_ns*kkrsz_ns block only........
c     returned in delta is 1-t*tau00^{-1}
c     tau00=(1-delta)^{-1}*t
c     ----------------------------------------------------------------
      time_new=time_direct
      call block_inv(wbig,vbig,nrmat_ns,nrmat_ns,kkrsz_ns,ipvt,
     &      blk_sz,nblk,delta,
     &      iwork,rwork,work1,iswtch,ied,
     &      nume,time_new,idcol,sym_ops,iprint)
      if(time_direct.eq.0.d0) then
	time_direct=time_new
	block_size=block_size+block_dir
      else if(time_old.eq.0.d0) then
c Make sure we are using LU at this energy and at the previous energy
	if(iswtch(ied-1).eq.2.and.iswtch(ied).eq.2) then
	if(time_new.le.time_direct) then
	time_old=time_direct
	time_direct=time_new
	block_size=block_size+block_dir
	else
	time_old=time_new
	block_dir=-block_dir
	block_size=block_size+2*block_dir
	endif
	endif
      else if(time_direct.lt.time_new) then
	block_size=block_size-block_dir
	block_dir=0
      else if(block_dir.ne.0) then
	time_old=time_direct
	time_direct=time_new
	block_size=block_size+block_dir
      endif
c     setup unit matrix...............................................
c     ----------------------------------------------------------------
      call cmtruni(wbig,kkrsz_ns)
c     ----------------------------------------------------------------
c     get 1-delta and put it in wbig
      call zaxpy(mtxsize,cmone,delta,1,wbig,1)
c     ================================================================
c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
c     ----------------------------------------------------------------
      call zgetrf(kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,info)
      call zcopy(kkrsz_ns*kkrsz_ns,tmat,1,tau00,1)
      call zgetrs('n',kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,tau00,
     &           kkrsz_ns,info)
c     ----------------------------------------------------------------
c*********************************TMP****************
!     if(istop.eq.sname) ied=2
!     if(ied.eq.2)then
      if((iprint.gt.0).or.(istop.eq.sname)) then
        write(6,'(''real part gettau_cl tau00'')')
        do ir=1,kkrsz
          write(6,'(d10.3,2x,3d10.3,2x,5d10.3,2x,7d10.3)')
     >(dble(tau00(ir+(l-1)*kkrsz)),l=1,kkrsz)
          if(ir.eq.1.or.ir.eq.4.or.ir.eq.9)write(6,'('' '')')
         enddo
         write(6,'(''imaginary part tau00'')')
         do ir=1,kkrsz
           write(6,'(d10.3,2x,3d10.3,2x,5d10.3,2x,7d10.3)')
     >(dimag(tau00(ir+(l-1)*kkrsz)),l=1,kkrsz)
         if(ir.eq.1.or.ir.eq.4.or.ir.eq.9)write(6,'('' '')')
         enddo
         endif
c*********************************TMP****************
c     ----------------------------------------------------------------
c     call inv_scale_tau00(tau00,kkrsz,kkrsz,lofk,n_spin_cant,kappa_rmt)
c        -------------------------------------------------------------
c     call inv_scale_tau00(tmat,kkrsz,kkrsz,lofk,n_spin_cant,kappa_rmt)
c        -------------------------------------------------------------
c     ================================================================
c     print if needed.................................................
      if(iprint.ge.1) then
         write(6,'('' gettaucl:: tau00:: mynod ='',1i5)') mynod
         call wrtmtx(tau00,kkrsz_ns,istop)
         call flush(6)
      endif
c     ================================================================
      nullify(gij, hfn, dlm)
      deallocate( iwork, rwork, work1, vbig )
      deallocate( msgids )
c     ================================================================
      if (istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
