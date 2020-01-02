!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getpotg(local_id,num_atoms,num_species,MaxNumSpecies,   &
                      nspin,n_spin_pola,efpot,evec,                   &
                      nr,nrmax,jmt,jws,r_mesh,jmax_pot,               &
                      vr,vdif,rhotot,xvalws,                          &
                      ztotss,zcorss,numc,nc,lc,kc,ec,lst,numcmax,     &
                      pot_l,isSphericalData,header,vunit,iprint,istop)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, TEN2m6, HALF, TWO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use MPPModule, only : sendMessage, recvMessage, waitMessage
   use MPPModule, only : setCommunicator, resetCommunicator, syncAllPEs
!
   use ParallelIOModule, only : isInputProc, getMyInputProc,          &
                                getNumInputClients, getInputClient,   &
                                getIOCommunicator, getMyPEinIOGroup,  &
                                getNumPEsInIOGroup
!
   use SystemModule, only : getAlloyTableSize,                        &
                            getNumAlloyElements, getAlloySpeciesIndex
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use LdaCorrectionModule, only : checkLdaCorrection, insertDataPackFromInput
!
   use InterpolationModule, only : getInterpolation
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: numcmax
   integer (kind=IntKind), intent(in) :: n_spin_pola
!
   character (len=*), intent(in) :: istop
   character (len=*), intent(out) :: header
   character(len=5), intent(out) :: lst(numcmax,n_spin_pola,num_species)
!
   character (len=5) :: ctmp
   character(len=1) :: cmsgbuf((160+40*numcmax)*2,MaxNumSpecies)
!
   logical, intent(out) :: isSphericalData
   logical :: checking_atoms = .true.
   logical :: checking_spin = .true.
!
   integer (kind=IntKind), intent(in) :: local_id
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: num_species
   integer (kind=IntKind), intent(in) :: maxNumSpecies
   integer (kind=IntKind), intent(in) :: nspin
   integer (kind=IntKind), intent(in) :: jmt
   integer (kind=IntKind), intent(in) :: jws
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: nrmax
   integer (kind=IntKind), intent(in) :: vunit
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind), intent(in) :: jmax_pot
   integer (kind=IntKind), intent(out) :: numc(num_species)
   integer (kind=IntKind), intent(out) :: nc(numcmax,n_spin_pola,num_species)
   integer (kind=IntKind), intent(out) :: lc(numcmax,n_spin_pola,num_species)
   integer (kind=IntKind), intent(out) :: kc(numcmax,n_spin_pola,num_species)
!
   integer (kind=IntKind) :: jmtmax
   integer (kind=IntKind) :: jwsmax
   integer (kind=IntKind) :: numcmax_in
   integer (kind=IntKind) :: is, ns, nspin_in, n_spin_pola_in
   integer (kind=IntKind) :: present_atom
   integer (kind=IntKind) :: fp_pos
   integer (kind=IntKind) :: msg_bytes
   integer (kind=IntKind) :: itmp
   integer (kind=IntKind) :: i, ig, ia
   integer (kind=IntKind) :: jl, nr_old, ir, jmt0, nj, lenc, lenf
   integer (kind=IntKind) :: slen, pad_bytes, table_size
   integer (kind=IntKind) :: imsgbuf(20,MaxNumSpecies)
!
   integer (kind=IntKind) :: num_clients
   integer (kind=IntKind) :: proc_client
   integer (kind=IntKind) :: integer4_size
   integer (kind=IntKind) :: real8_size
   integer (kind=IntKind) :: dsize_ldapu, dsize_nspot, dsize_min
   integer (kind=IntKind) :: comm, MyPEinGroup, NumPEsInGroup
   integer (kind=IntKind) :: na_in
!
   integer (kind=IntKind), parameter :: fsize = 10000 ! > (9+2*nrmax)*n_spin_pola
!
   real (kind=RealKind), intent(in) :: r_mesh(nr)
   real (kind=RealKind), intent(in) :: ztotss(num_species)
   real (kind=RealKind), intent(in) :: zcorss(num_species)
   real (kind=RealKind), intent(out) :: evec(3)
   real (kind=RealKind), intent(out) :: efpot(num_species)
   real (kind=RealKind), intent(out) :: vr(nr,n_spin_pola,num_species)
   real (kind=RealKind), intent(out) :: rhotot(nr,n_spin_pola,num_species)
   real (kind=RealKind), intent(out) :: xvalws(n_spin_pola,num_species)
   real (kind=RealKind), intent(out) :: vdif
   real (kind=RealKind), intent(out) :: ec(numcmax,n_spin_pola,num_species)
   real (kind=RealKind) :: fspace(fsize,MaxNumSpecies)
!  real (kind=RealKind), intent(out) :: fspace((9+2*nrmax)*n_spin_pola)
!
   real (kind=RealKind) :: alat
   real (kind=RealKind) :: rtmp, xst, xmt, hh, dvr
   real (kind=RealKind), allocatable :: data_ldapu(:,:)
   real (kind=RealKind), allocatable, target :: data_nspot(:,:)
   real (kind=RealKind), allocatable :: r_mesh_old(:)
!
   complex (kind=CmplxKind), target :: pot_l(nr,jmax_pot,n_spin_pola,num_species)
   complex (kind=CmplxKind), allocatable, target :: c_data_nspot(:)
   complex (kind=CmplxKind), pointer :: p_data_nspot(:)
!
   interface
      function getAtomPointer(vunit,ig,ia,table_size,integer4_size,real8_size) result(p)
         use KindParamModule, only : IntKind, RealKind
         implicit none
!
         integer (kind=IntKind), intent(in) :: vunit, ig, ia
         integer (kind=IntKind), intent(in) :: table_size
         integer (kind=IntKind), intent(in) :: integer4_size,real8_size
         integer (kind=IntKind) :: p
      end function getAtomPointer
   end interface
!
!  *******************************************************************
!  cccccccccc read in the potentials for current sublattice  ccccccccc
!  *******************************************************************
!
   MyPEinGroup = getMyPEinIOGroup()
   NumPEsInGroup = getNumPEsInIOGroup()
   comm = getIOCommunicator()
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!  -------------------------------------------------------------------
   rhotot = ZERO
   imsgbuf = 0
   call c_dtsize(integer4_size,real8_size)
   dsize_ldapu = 0
   dsize_nspot = 0
!  -------------------------------------------------------------------
!
!  ===================================================================
!  using gopen to perform the global openning......................
!  -------------------------------------------------------------------
   if( isInputProc() ) then
!     ================================================================
!     check the number of atoms in the potential file.................
!     ================================================================
      if (checking_atoms) then
         call c_fseek(vunit,1,0)
         call c_read_integer(vunit,imsgbuf(1:4,1),4)
         if (imsgbuf(4,1) > 10) then
            table_size = (imsgbuf(4,1)-mod(imsgbuf(4,1),10))/10
         else
            table_size = 0 ! number of sites
            fp_pos = 1
            LOOP_do: do
               call c_fseek(vunit,fp_pos,0)
               call c_read_integer(vunit,imsgbuf(1:10,1),10)
               fp_pos = fp_pos + integer4_size*10
               if (imsgbuf(1,1) < 0) then
!                 write(6,'(a)')'The potential contains modified data format.........!'
                  imsgbuf(1,1) = -imsgbuf(1,1)
               endif
               if (imsgbuf(1,1) /= table_size + 1 .and.               &
                   imsgbuf(1,1) /= table_size) then
                  exit LOOP_do
               endif
               table_size = imsgbuf(1,1)
            enddo LOOP_do
            if (table_size /= num_atoms) then
               write(6,'(a,i6)')'The number atoms in the potential file = ',table_size
            else if (table_size < 1) then
               call ErrorHandler('getpot','table_size < 1',table_size)
            endif
         endif
      else ! Assume the potential is consistent with the position data.
         table_size = getAlloyTableSize()
      endif
!     ================================================================
!     read in the potentials..........................................
!     ================================================================
      num_clients = getNumInputClients()
      do i = 1,num_clients
         proc_client = getInputClient(i)
         ig = getGlobalIndex(local_id,proc_client)
         do ia = 1, getNumAlloyElements(ig)
!           ----------------------------------------------------------
            present_atom = getAtomPointer(vunit,ig,ia,table_size,     &
                                          integer4_size,real8_size)
!           ----------------------------------------------------------
            if (present_atom < 1) then
               call ErrorHandler('getpotg','My atom potential is not found in the potential data',ig)
            endif
!           ==========================================================
            fp_pos=(present_atom-1)*integer4_size*10+1
!           ----------------------------------------------------------
            call c_fseek(vunit,fp_pos,0)
            call c_read_integer(vunit,imsgbuf(1:10,ia),10)
!           ----------------------------------------------------------
            if (imsgbuf(3,ia) > fsize) then
               call ErrorHandler('getpotg','imsgbuf(3) > size(fspace): ', &
                                 present_atom,imsgbuf(3,ia),fsize)
            else if (imsgbuf(10,ia) > numcmax) then
               call ErrorHandler('getpotg','imsgbuf(10) > numcmax',imsgbuf(10,ia),numcmax)
            endif
            if (imsgbuf(1,ia) < 0) then
               fp_pos=table_size*integer4_size*10+fp_pos
               call c_fseek(vunit,fp_pos,0)
               call c_read_integer(vunit,imsgbuf(11:20,ia),10)
               if (dsize_ldapu < imsgbuf(13,ia)) then
                  if ( allocated(data_ldapu) ) then
                     deallocate( data_ldapu )
                  endif
                  allocate( data_ldapu(imsgbuf(13,ia),MaxNumSpecies) )
                  dsize_ldapu = imsgbuf(13,ia)
               endif
               if (dsize_nspot < imsgbuf(14,ia)) then
                  if ( allocated(data_nspot) ) then
                     deallocate( data_nspot )
                  endif
                  allocate( data_nspot(imsgbuf(14,ia),MaxNumSpecies) )
                  dsize_nspot = imsgbuf(14,ia)
               endif
            endif
!           ==========================================================
!           reading cmsgbuf, fspace and evec arrays from vfile........
!           ==========================================================
            call c_string_padsize(vunit,imsgbuf(2,ia),pad_bytes)
!                print *,'imsgbuf = ',imsgbuf(2,ia),', pad_bytes = ',pad_bytes
            msg_bytes=imsgbuf(2,ia)+pad_bytes+(imsgbuf(3,ia)+3)*real8_size
            if (imsgbuf(1,ia) > 0) then
               fp_pos=table_size*integer4_size*10+(present_atom-1)*msg_bytes+1
            else
               fp_pos=imsgbuf(12,ia)*real8_size +                     &
                      table_size*integer4_size*20+(present_atom-1)*msg_bytes+1
            endif
!           ----------------------------------------------------------
            call c_fseek(vunit,fp_pos,0)
            call c_read_string(vunit,cmsgbuf(:,ia),imsgbuf(2,ia),slen)
            call c_read_double(vunit,fspace(:,ia),imsgbuf(3,ia))
            call c_read_double(vunit,evec,3)
!           ----------------------------------------------------------
            if (imsgbuf(13,ia) > 0) then  
               call c_read_double(vunit,data_ldapu(:,ia),imsgbuf(13,ia))
            endif
            if (imsgbuf(14,ia) > 0) then
               call c_read_double(vunit,data_nspot(:,ia),imsgbuf(14,ia))
            endif
!           ----------------------------------------------------------
            if (mod(imsgbuf(4,ia),10) < 3) then
               evec(1)=0.0d0
               evec(2)=0.0d0
               evec(3)=1.0d0
            endif
!           ----------------------------------------------------------
            call sendMessage(imsgbuf(:,ia),20,23457,proc_client)
            call sendMessage(cmsgbuf(:,ia),imsgbuf(2,ia),23458,proc_client)
            call sendMessage(fspace(:,ia),imsgbuf(3,ia),23459,proc_client)
            call sendMessage(evec,3,23465,proc_client)
!           ----------------------------------------------------------
            if (imsgbuf(13,ia) > 0) then  
!              -------------------------------------------------------
               call sendMessage(data_ldapu(:,ia),imsgbuf(13,ia),23466,proc_client)
!              -------------------------------------------------------
            endif
            if (imsgbuf(14,ia) > 0) then  
!              -------------------------------------------------------
               call sendMessage(data_nspot(:,ia),imsgbuf(14,ia),23467,proc_client)
!              -------------------------------------------------------
            endif
         enddo ! ia
      enddo ! i
!
!     ================================================================
!     Read potential data for local_id site
!     ================================================================
      ig = getGlobalIndex(local_id)
      do ia = 1, num_species
!        -------------------------------------------------------------
         present_atom = getAtomPointer(vunit,ig,ia,table_size,        &
                                       integer4_size,real8_size)
!        -------------------------------------------------------------
         fp_pos=(present_atom-1)*integer4_size*10+1
!        -------------------------------------------------------------
         call c_fseek(vunit,fp_pos,0)
         call c_read_integer(vunit,imsgbuf(:,ia),10)
         if (imsgbuf(3,ia) > fsize) then
            call ErrorHandler('getpotg','imsgbuf(3) > size(fspace): ',&
                              present_atom,imsgbuf(3,ia),fsize)
         else if (imsgbuf(10,ia) > numcmax) then
            call ErrorHandler('getpotg','imsgbuf(10) > numcmax',imsgbuf(10,ia),numcmax)
         endif
         if (imsgbuf(1,ia) < 0) then
            fp_pos=table_size*integer4_size*10+fp_pos
!           ----------------------------------------------------------
            call c_fseek(vunit,fp_pos,0)
            call c_read_integer(vunit,imsgbuf(11:20,ia),10)
!           ----------------------------------------------------------
            if (dsize_ldapu < imsgbuf(13,ia)) then
               if ( allocated(data_ldapu) ) then
                  deallocate( data_ldapu )
               endif
               allocate( data_ldapu(imsgbuf(13,ia),MaxNumSpecies) )
               dsize_ldapu = imsgbuf(13,ia)
            endif
            if (dsize_nspot < imsgbuf(14,ia)) then
               if ( allocated(data_nspot) ) then
                  deallocate( data_nspot )
               endif
               allocate( data_nspot(imsgbuf(14,ia),MaxNumSpecies) )
               dsize_nspot = imsgbuf(14,ia)
            endif
         endif
!        =============================================================
!        reading cmsgbuf, fspace and evec arrays from vfile...........
!        =============================================================
         call c_string_padsize(vunit,imsgbuf(2,ia),pad_bytes)
!             print *,'pad_bytes = ',pad_bytes
         msg_bytes=imsgbuf(2,ia)+pad_bytes+(imsgbuf(3,ia)+3)*real8_size
         if (imsgbuf(1,ia) > 0) then
            fp_pos=table_size*integer4_size*10+(present_atom-1)*msg_bytes+1
         else
            fp_pos=imsgbuf(12,ia)*real8_size +                        &
                   table_size*integer4_size*20+(present_atom-1)*msg_bytes+1
         endif
!        -------------------------------------------------------------
         call c_fseek(vunit,fp_pos,0)
         call c_read_string(vunit,cmsgbuf(:,ia),imsgbuf(2,ia),slen)
         call c_read_double(vunit,fspace(:,ia),imsgbuf(3,ia))
         call c_read_double(vunit,evec,3)
!        -------------------------------------------------------------
         if(mod(imsgbuf(4,ia),10) < 3) then
            evec(1)=0.0d0
            evec(2)=0.0d0
            evec(3)=1.0d0
         endif
         if (imsgbuf(13,ia) > 0) then  
            call c_read_double(vunit,data_ldapu(:,ia),imsgbuf(13,ia))
         endif
         if (imsgbuf(14,ia) > 0) then
            call c_read_double(vunit,data_nspot(:,ia),imsgbuf(14,ia))
         endif
!
         if (checkLdaCorrection(local_id,1) .and. imsgbuf(13,ia) > 1) then
            call insertDataPackFromInput(local_id,1,imsgbuf(13,ia),data_ldapu(:,ia))
         endif
      enddo
   else
      do ia = 1, num_species
!        -------------------------------------------------------------
         call recvMessage(imsgbuf(:,ia),20,23457,getMyInputProc())
         call recvMessage(cmsgbuf(:,ia),imsgbuf(2,ia),23458,getMyInputProc())
         call recvMessage(fspace(:,ia),imsgbuf(3,ia),23459,getMyInputProc())
         call recvMessage(evec,3,23465,getMyInputProc())
!        -------------------------------------------------------------
         if (imsgbuf(1,ia) < 0) then
            if (dsize_ldapu < imsgbuf(13,ia)) then
               if ( allocated(data_ldapu) ) then
                  deallocate( data_ldapu )
               endif
               allocate( data_ldapu(imsgbuf(13,ia),MaxNumSpecies) )
               dsize_ldapu = imsgbuf(13,ia)
            endif
            if (dsize_nspot < imsgbuf(14,ia)) then
               if ( allocated(data_nspot) ) then
                  deallocate( data_nspot )
               endif
               allocate( data_nspot(imsgbuf(14,ia),MaxNumSpecies) )
               dsize_nspot = imsgbuf(14,ia)
            endif
         endif
         if (imsgbuf(13,ia) > 0) then  
!           ----------------------------------------------------------
            call recvMessage(data_ldapu(:,ia),imsgbuf(13,ia),23466,getMyInputProc())
!           ----------------------------------------------------------
         endif
         if (imsgbuf(14,ia) > 0) then  
!           ----------------------------------------------------------
            call recvMessage(data_nspot(:,ia),imsgbuf(14,ia),23467,getMyInputProc())
!           ----------------------------------------------------------
         endif
         if (checkLdaCorrection(local_id,1) .and. imsgbuf(13,ia) > 1) then
            call insertDataPackFromInput(local_id,1,imsgbuf(13,ia),data_ldapu(:,ia))
         endif
!
      enddo
   endif
!
   if (checking_spin) then
      nspin_in = mod(imsgbuf(4,1),10)
      if (nspin_in /= nspin) then
         write(6,'(a,i6)')'The spin parameter in the potential file = ',nspin_in
         n_spin_pola_in = min(2,nspin_in)
      else
         n_spin_pola_in = n_spin_pola
      endif
   else
      n_spin_pola_in = n_spin_pola
   endif
!
   if (imsgbuf(1,1) < 0) then
      isSphericalData = .false.
   else
      isSphericalData = .true.
   endif
!
!  ===================================================================
!  decoding cmsgbuf and fspace.....................................
!  Note: fspace will be scratched, we need to store xst and xmt here in
!        case they are needed later
!  ===================================================================
   xst=fspace(8,1)
   xmt=fspace(9,1)
   do ia = 1, num_species
      jmtmax=imsgbuf(8,ia)
      jwsmax=imsgbuf(9,ia)
      numcmax_in=imsgbuf(10,ia)
      lenc = 160+40*numcmax_in
      lenf = 9+jmtmax+jwsmax
      do is = 1, n_spin_pola_in
!        In case of n_spin_pola_in /= n_spin_pola in the potential data file
         if (n_spin_pola_in /= n_spin_pola) then
            ns = 1
         else
            ns = is
         endif
!        -------------------------------------------------------------
         call potredg(jmtmax, jwsmax, nr, numcmax_in, imsgbuf(5:7,ia),&
                      cmsgbuf(lenc*(ns-1)+1:lenc*ns,ia),              &
                      fspace(lenf*(ns-1)+1:lenf*ns,ia),               &
                      alat,efpot(ia),                                 &
                      jmt,jws,r_mesh,vr(1:nr,ns,ia),vdif,             &
                      rhotot(1:nr,ns,ia),                             &
                      xvalws(ns,ia),ztotss(ia),zcorss(ia),numc(ia),   &
                      nc(1:numcmax_in,ns,ia),lc(1:numcmax_in,ns,ia),  &
                      kc(1:numcmax_in,ns,ia),ec(1:numcmax_in,ns,ia),  &
                      lst(1:numcmax_in,ns,ia),header,iprint,istop)
!        -------------------------------------------------------------
      enddo
!
      if (n_spin_pola_in < n_spin_pola) then
         cmsgbuf(lenc+1:lenc*2,ia) = cmsgbuf(1:lenc,ia)
         fspace(lenf+1:lenf*2,ia)= fspace(1:lenf,ia)
         vr(1:nr,2,ia) = vr(1:nr,1,ia)
         rhotot(1:nr,1,ia) = HALF*rhotot(1:nr,1,ia)
         rhotot(1:nr,2,ia) = rhotot(1:nr,1,ia)
         xvalws(1,ia) = HALF*xvalws(1,ia)
         xvalws(2,ia) = xvalws(1,ia)
         nc(1:numcmax_in,2,ia) = nc(1:numcmax_in,1,ia)
         lc(1:numcmax_in,2,ia) = lc(1:numcmax_in,1,ia)
         kc(1:numcmax_in,2,ia) = kc(1:numcmax_in,1,ia)
         ec(1:numcmax_in,2,ia) = ec(1:numcmax_in,1,ia)
      else if (n_spin_pola_in > n_spin_pola) then
         rhotot(1:nr,1,ia) = TWO*rhotot(1:nr,1,ia)
         xvalws(1,ia) = TWO*xvalws(1,ia)
      endif
!
!     ================================================================
!     check if the moment is negative.................................
!     if necessary, make the moment positive, put evec to its opposite 
!     direction, and rearrange the potential and density..............
!     fspace is used as a scratch space now.
!     ================================================================
      if(xvalws(1,ia)-xvalws(n_spin_pola,ia).lt.-0.001 .and. nspin.eq.3)then
         vdif=-vdif
         rtmp=xvalws(1,ia)
         xvalws(1,ia)=xvalws(n_spin_pola,ia)
         xvalws(n_spin_pola,ia)=rtmp
         evec(1)=-evec(1)
         evec(2)=-evec(2)
         evec(3)=-evec(3)
!        -------------------------------------------------------------
!        call dcopy(jmt,vr(1,1,ia),1,fspace(1,ia),1)
         fspace(1:jmt,ia) = vr(1:jmt,1,ia)
!        -------------------------------------------------------------
!        call dcopy(jmt,vr(1,n_spin_pola,ia),1,vr(1,1,ia),1)
         vr(1:jmt,1,ia) = vr(1:jmt,n_spin_pola,ia)
!        -------------------------------------------------------------
!        call dcopy(jmt,fspace(1,ia),1,vr(1,n_spin_pola,ia),1)
         vr(1:jmt,n_spin_pola,ia) = fspace(1:jmt,ia)
!        -------------------------------------------------------------
!        call dcopy(jws,rhotot(1,1,ia),1,fspace(1,ia),1)
         fspace(1:jws,ia) = rhotot(1:jws,1,ia)
!        -------------------------------------------------------------
!        call dcopy(jws,rhotot(1,n_spin_pola,ia),1,rhotot(1,1,ia),1)
         rhotot(1:jws,1,ia) = rhotot(1:jws,n_spin_pola,ia)
!        -------------------------------------------------------------
!        call dcopy(jws,fspace(1,ia),1,rhotot(1,n_spin_pola,ia),1)
         rhotot(1:jws,n_spin_pola,ia) = fspace(1:jws,ia)
!        -------------------------------------------------------------
         do i = 1, numc(ia)
            itmp=nc(i,1,ia)
            nc(i,1,ia)=nc(i,n_spin_pola,ia)
            nc(i,n_spin_pola,ia)=itmp
            itmp=lc(i,1,ia)
            lc(i,1,ia)=lc(i,n_spin_pola,ia)
            lc(i,n_spin_pola,ia)=itmp
            itmp=kc(i,1,ia)
            kc(i,1,ia)=kc(i,n_spin_pola,ia)
            kc(i,n_spin_pola,ia)=itmp
            rtmp=ec(i,1,ia)
            ec(i,1,ia)=ec(i,n_spin_pola,ia)
            ec(i,n_spin_pola,ia)=rtmp
            ctmp=lst(i,1,ia)
            lst(i,1,ia)=lst(i,n_spin_pola,ia)
            lst(i,n_spin_pola,ia)=ctmp
         enddo
      endif
   enddo
!
   call resetCommunicator()
!
!  ===================================================================
!  Needs to take care the radial mesh interpolation...
!  ===================================================================
   if ( dsize_nspot > 0) then
!     ================================================================
!     perform interpolation: added by Yang @08/06/14
!     Note : imsgbuf(4) = nspin_old; imsgbuf(17) = jmax_pot_old
!     ================================================================
      nr_old = dsize_nspot/(2*imsgbuf(17,1)*n_spin_pola*num_species)
      allocate(r_mesh_old(nr_old), c_data_nspot(dsize_nspot/2))
      jmt0=imsgbuf(5,1)
      hh=(xmt-xst)/real(jmt0-1,kind=RealKind) ! Assuming hin = hout
      do ir=1,nr_old
         r_mesh_old(ir)=exp(xst+(ir-1)*hh)
      enddo
      c_data_nspot = transfer(data_nspot,c_data_nspot)
      nj = nr_old*imsgbuf(17,1)
      do ia = 1, num_species
         do ns = 1, n_spin_pola
            if (nr /= nr_old .or. abs(r_mesh(jmt)-r_mesh_old(jmt0)) > TEN2m6) then
               do jl = 1, min(imsgbuf(17,ia),jmax_pot)
                  p_data_nspot => c_data_nspot(nj*n_spin_pola*(ia-1)+          &
                                               nj*(ns-1)+nr_old*(jl-1)+1:      &
                                               nj*n_spin_pola*(ia-1)+nj*(ns-1)+nr_old*jl)
                  do ir = 1, nr
!                    -------------------------------------------------
                     pot_l(ir,jl,ns,ia) = getInterpolation(nr_old,r_mesh_old,  &
                                                           p_data_nspot,r_mesh(ir),dvr)      
!                    -------------------------------------------------
                  enddo
               enddo
            else
               p_data_nspot => c_data_nspot(nj*n_spin_pola*(ia-1)+nj*(ns-1)+1: &
                                            nj*n_spin_pola*(ia-1)+nj*ns)
!              -------------------------------------------------------
               call zcopy(nr*min(imsgbuf(17,ia),jmax_pot),p_data_nspot,1,pot_l(1,1,ns,ia),1)
!              -------------------------------------------------------
            endif
         enddo
      enddo
      deallocate(r_mesh_old, c_data_nspot)
   endif
!
   if ( allocated(data_ldapu) ) then
      deallocate( data_ldapu )
   endif
   if ( allocated(data_nspot) ) then
      deallocate( data_nspot )
   endif
!
   end subroutine getpotg
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomPointer(vunit,ig,ia,table_size,integer4_size,real8_size) result(p)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use SystemModule, only : getAlloyTableSize, getAlloySpeciesIndex,  &
                            getAtomicNumber
   implicit none
!
   integer (kind=IntKind), intent(in) :: vunit, ig, ia
   integer (kind=IntKind), intent(in) :: table_size
   integer (kind=IntKind), intent(in) :: integer4_size,real8_size
   integer (kind=IntKind) :: p
   integer (kind=IntKind) :: it, msg_bytes, imsgbuf(20), fp_pos, pad_bytes
!
   real (kind=RealKind) :: za
!
   if (table_size /= getAlloyTableSize()) then
      p = -1
!     ================================================================
!     Loop over the atoms in the potential file to find the matching atom
!     ================================================================
      LOOP_it: do it = 1, table_size
         fp_pos=(it-1)*integer4_size*10+1
!        -------------------------------------------------------------
         call c_fseek(vunit,fp_pos,0)
         call c_read_integer(vunit,imsgbuf,10)
         call c_string_padsize(vunit,imsgbuf(2),pad_bytes)
!        -------------------------------------------------------------
         msg_bytes=imsgbuf(2)+pad_bytes+(imsgbuf(3)+3)*real8_size
         if (imsgbuf(1) > 0) then
            fp_pos=table_size*integer4_size*10+(it-1)*msg_bytes+1
         else
            fp_pos=table_size*integer4_size*10+fp_pos
!           ----------------------------------------------------------
            call c_fseek(vunit,fp_pos,0)
            call c_read_integer(vunit,imsgbuf(11:20),10)
!           ----------------------------------------------------------
            fp_pos=imsgbuf(12)*real8_size+table_size*integer4_size*20+(it-1)*msg_bytes+1
         endif
         fp_pos = fp_pos + imsgbuf(2) + 3*real8_size
!        -------------------------------------------------------------
         call c_fseek(vunit,fp_pos,0)
         call c_read_double(vunit,za,1)
!        -------------------------------------------------------------
         if (int(za) /= getAtomicNumber(ig,ia)) then
            write(6,'(a,2i5)')'Za in file /= Za of my atom: ',int(za),getAtomicNumber(ig,ia)
         else
            p = it
            exit LOOP_it
         endif
      enddo LOOP_it
   else
      p = getAlloySpeciesIndex(ig,ia)
   endif
!
   end function getAtomPointer
