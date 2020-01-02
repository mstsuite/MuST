program extractPotential  
!
!   use MPPModule, only : initMPP, endMPP
!   use MPPModule, only : MyPE, NumPEs
!   use MPPModule, only : bcastMessage
!
    use ChemElementModule, only : getName, getZcor, getZsem, getZval
!
    implicit   none
!
    character (len=(160+40*30)*2) :: cmsgbuf
    character (len=2) :: atname
    character (len=5) :: stmp
    character (len=5), allocatable :: lst(:)
    character (len=80) :: text
    character (len=80) :: header, jtitle
    character (len=80) :: wdir
    character (len=80) :: vpname
    character (len=80) :: vpfile
    character (len=80) :: outfile
!
    integer, parameter :: MyPE = 0
    integer, parameter :: NUmPEs = 1
!
    integer :: iparam(2)
    integer :: imsgbuf(10)
    integer :: num_atoms, atom_id
    integer :: i0, i1, istep, id, ns, ir, slen
    integer :: zt, zs, zc, zv
    integer :: vunit, fp_pos
    integer :: integer4_size, real8_size
    integer :: nspin, n_spin_pola, n_spin_cant
    integer :: jmtmax,jwsmax,numcmax
    integer :: msg_bytes
    integer :: jmt,jws,numc, str_len
    integer :: cms, fms, j, ncd, n
    integer, allocatable :: nc(:)
    integer, allocatable :: lc(:)
    integer, allocatable :: kc(:)
    integer, parameter :: ounit = 10
!
    real*8, allocatable  :: r_mesh(:)
    real*8, allocatable  :: vr(:)
    real*8, allocatable  :: rhotot(:)
    real*8, allocatable  :: ec(:)
    real*8, allocatable  :: fmsgbuf(:)
    real*8 :: evec(3)
    real*8 :: efermi
    real*8 :: xmt, xstart, alat
    real*8 :: vdif, vzero
    real*8 :: xvalws
    real*8 :: ztotss,zcorss
    real*8, parameter :: ZERO = 0.0d0
    real*8, parameter :: ONE = 1.0d0
!
!   ==================================================================
!   Initilize MPP parameters..........................................
!   ------------------------------------------------------------------
!ywgcall InitMPP()
!   ------------------------------------------------------------------
    call c_dtsize(integer4_size,real8_size)
!   ------------------------------------------------------------------
!
    if (MyPE == 0) then
        read(5,'(a)')text
        read(5,*)iparam(1)
        read(5,'(a)')text
        read(5,*)iparam(2)
        read(5,'(a)')text
        read(5,'(a)')wdir
        read(5,'(a)')text
        read(5,'(a)')vpname
    endif
!
!   ------------------------------------------------------------------
!ywgcall bcastMessage(iparam,2,0)
!ywgcall bcastMessage(wdir,0)
!ywgcall bcastMessage(vpname,0)
!   ------------------------------------------------------------------
!
    num_atoms=iparam(1)
    atom_id=iparam(2)
!
    print *,'num_atoms = ',num_atoms
    print *,'atom_id = ',atom_id
!
    if (atom_id <= 0) then
        i0 = MyPE+1
        i1 = num_atoms
        istep = NumPEs
    else
        i0 = atom_id
        i1 = atom_id-MyPE
        istep = 1
    endif
!
    vpfile = trim(adjustl(wdir))//adjustl(vpname)
    str_len = len(trim(vpfile)) 
!
    if (i1 >= i0) then
        call c_gopen(vunit,trim(vpfile),str_len,'old',3)

        fp_pos=MyPE*integer4_size*10+1
!       --------------------------------------------------------------
        call c_fseek(vunit,fp_pos,0)
!       --------------------------------------------------------------
        call c_read_integer(vunit,imsgbuf,integer4_size*10)
!       --------------------------------------------------------------
        if (imsgbuf(1) /= MyPE+1) then
            write(6,'(a,i5,a,i5)')'ERROR: MyPE+1 = ',MyPE+1,          &
                                  ', imsgbuf(1) = ',imsgbuf(1)
        endif

        nspin=imsgbuf(4)
        jmtmax=imsgbuf(8)
        jwsmax=imsgbuf(9)
        numcmax=imsgbuf(10)

        print *,'nspin = ',nspin
        print *,'jmtmax = ',jmtmax
        print *,'jwsmax = ',jwsmax
        print *,'numcmax = ',numcmax
!
        if (nspin == 1) then
            n_spin_pola = 1
            n_spin_cant = 1
        else if (nspin == 2) then
            n_spin_pola = 2
            n_spin_cant = 1
        else if (nspin == 3) then
            n_spin_pola = 2
            n_spin_cant = 2
        else
            write(6,'(a,i5)')'ERROR: nspin = ',nspin
        endif

        allocate( fmsgbuf((9+jmtmax+jwsmax)*n_spin_pola) )
        allocate( vr(1:jwsmax) )
        allocate( rhotot(1:jwsmax) )
        allocate( ec(1:numcmax) )
        allocate( lc(1:numcmax) )
        allocate( kc(1:numcmax) )
        allocate( nc(1:numcmax) )
        allocate( lst(1:numcmax) )
!
        do id = i0, i1, istep
            write(stmp,'(i5)')10000+id
            outfile=trim(wdir)//'pot_'//stmp(2:5)
            print *,outfile

            open(unit=ounit,file=outfile,form='formatted',status='unknown')

            fp_pos=(id-1)*integer4_size*10+1
!           ----------------------------------------------------------
            call c_fseek(vunit,fp_pos,0)
!           ----------------------------------------------------------
            call c_read_integer(vunit,imsgbuf,integer4_size*10)
!           ----------------------------------------------------------
            jmt=imsgbuf(5)
            jws=imsgbuf(6)
            numc=imsgbuf(7)
            print *,'jmt = ',jmt
            print *,'jws = ',jws
            print *,'numc = ',numc
!
!           ==========================================================
!           reading cmsgbuf, fmsgbuf and evec arrays from vfile........
!           ==========================================================
            msg_bytes=imsgbuf(2)+(imsgbuf(3)+3)*real8_size
            fp_pos=num_atoms*integer4_size*10+(id-1)*msg_bytes+1
!           ----------------------------------------------------------
            call c_fseek(vunit,fp_pos,0)
!           ----------------------------------------------------------
        print *,'imsgbuf(3) =',imsgbuf(3)
        print *,'integer4_size =',integer4_size
        print *,'real8_size =',real8_size
        print *,'msg_bytes =',msg_bytes
            print *,'fp_pos =',fp_pos
            print *,'imsgbuf(2) =',imsgbuf(2)
            call c_read_string(vunit,cmsgbuf,imsgbuf(2),slen)
!           ----------------------------------------------------------
            call c_read_double(vunit,fmsgbuf,imsgbuf(3)*real8_size)
!           ----------------------------------------------------------
            if (nspin == 3) then
!               ------------------------------------------------------
   	        call c_read_double(vunit,evec,3*real8_size)
!               ------------------------------------------------------
            else
                evec(1)=ZERO
                evec(2)=ZERO
                evec(3)=ONE
            endif
!
!           ==========================================================
!           decoding cmsgbuf and fmsgbuf..............................
!           ==========================================================
            do ns=1,n_spin_pola
                cms = (160+40*numcmax)*(ns-1)
                header=cmsgbuf(cms+1:cms+80)
                jtitle=cmsgbuf(cms+81:cms+160)
                print *,header
                print *,jtitle

                do j=1,numc
!                   --------------------------------------------------
                    text=cmsgbuf(cms+160+(j-1)*40+1:cms+160+j*40)
!                   --------------------------------------------------
                    read(text,'(3i5,1e20.13,1a5)')nc(j),lc(j),kc(j),ec(j),lst(j)
                    print *,text
                enddo
!
!               ======================================================
!               decoding fmsgbuf......................................
!               ======================================================
                fms    = (9+jmtmax+jwsmax)*(ns-1)
                alat   = fmsgbuf(fms+1)
                efermi = fmsgbuf(fms+2)
                vdif   = fmsgbuf(fms+3)
                ztotss = fmsgbuf(fms+4)
                zcorss = fmsgbuf(fms+5)
                xvalws = fmsgbuf(fms+6)
                vzero  = fmsgbuf(fms+7)
                xstart = fmsgbuf(fms+8)
                xmt    = fmsgbuf(fms+9)
!
                if (ns == 1) then
                    write(ounit,'(a)') header
                    write(ounit,'(i5,2x,d20.13)') n_spin_pola,vdif
                    zt = ztotss
!                   call getaninfo(atname,zt,zc,zs,zv)
                    atname = getName(zt)
                    zc = getZcor(zt)
                    zs = getZsem(zt)
                    zv = getZval(zt)
                endif
                write(ounit,'(a,t18,a2,t25,a,f4.0,t35,a,f10.5)')      &
     &                      ' LSMS:',atname,'z=',ztotss,'xvalws=',xvalws
                write(ounit,'(f5.0,17x,f12.5,f5.0,d20.13)')           &
     &                        ztotss,alat,zcorss,efermi
                write(ounit,'(17x,2d20.13,2i5)') xstart,xmt,jmt,jws
!
                n=9+fms
                do ir=1,jmt
                    vr(ir)=fmsgbuf(n+ir)
                enddo
                write(ounit,'(4d20.13)') (vr(ir),ir=1,jmt)
                write(ounit,'(35x,d20.13)') vzero
!
                n=9+jmtmax+fms
                do ir=1,jws
                   rhotot(ir)=fmsgbuf(n+ir)
                enddo
                write(ounit,'(i5,d20.13)') jws,xvalws
                write(ounit,'(4d20.13)') (rhotot(ir),ir=1,jws)

                ncd = 0
                write(ounit,'(2i5)') numc,ncd
                do j=1,numc
                   write(ounit,'(3i5,f12.5,2x,a5)')                   &
     &                   nc(j),lc(j),kc(j),ec(j),lst(j)
                enddo
            enddo
            close(unit=ounit)
        enddo
        call c_close(vunit)
    endif
!
!   ------------------------------------------------------------------
!ywgcall endMPP()
!   ------------------------------------------------------------------
    stop 'Ok'
end program extractPotential
