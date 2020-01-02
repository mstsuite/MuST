!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine potredg(jmtmax,jwsmax,nr,numcmax,imsgbuf,                   &
                      cmsgbuf,fmsgbuf,alat,efpot,                         &
                      jmt,jws,r_mesh,vr,vdif,rhotot,xvalws,ztotss,zcorss, &
                      numc,nc,lc,kc,ec,lst,header,iprint,istop)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ZERO, TWO, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: numcmax
!
   character (len=*), intent(in) :: istop
   character (len=1), intent(in) :: cmsgbuf(160+40*numcmax)
   character (len=5), intent(out) :: lst(numcmax)
   character (len=*), intent(out) :: header
!
   character (len=40) :: temp_buf
   character (len=80) :: jtitle
!
   integer (kind=IntKind), intent(in) :: jmtmax
   integer (kind=IntKind), intent(in) :: jwsmax
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: imsgbuf(3)
   integer (kind=IntKind), intent(in) :: jmt
   integer (kind=IntKind), intent(in) :: jws
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind), intent(out) :: numc
   integer (kind=IntKind), intent(out) :: nc(numcmax)
   integer (kind=IntKind), intent(out) :: lc(numcmax)
   integer (kind=IntKind), intent(out) :: kc(numcmax)
!
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: jmt0
   integer (kind=IntKind) :: ntchg
   integer (kind=IntKind) :: ir
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: nr_old
!
   real (kind=RealKind), intent(in) :: r_mesh(nr)
   real (kind=RealKind), intent(in) :: fmsgbuf(9+jmtmax+jwsmax)
   real (kind=RealKind), intent(in) :: ztotss
   real (kind=RealKind), intent(in) :: zcorss
   real (kind=RealKind), intent(out) :: vr(nr)
   real (kind=RealKind), intent(out) :: rhotot(nr)
   real (kind=RealKind), intent(out) :: vdif
   real (kind=RealKind), intent(out) :: xvalws
   real (kind=RealKind), intent(out) :: ec(numcmax)
   real (kind=RealKind), intent(out) :: efpot
!
   real (kind=RealKind) :: alat
   real (kind=RealKind) :: azed
   real (kind=RealKind) :: zcd
   real (kind=RealKind) :: xst
   real (kind=RealKind) :: xmt
   real (kind=RealKind) :: hh
!  real (kind=RealKind) :: ylag
   real (kind=RealKind) :: dvr, rfac
   real (kind=RealKind) :: vzero
!
   real (kind=RealKind), allocatable :: x_mesh_old(:)
   real (kind=RealKind), allocatable :: r_mesh_old(:)
   real (kind=RealKind), allocatable :: vrold(:)
   real (kind=RealKind), allocatable :: chgold(:)
!
   logical :: performed_Interpolation
!
   interface
      subroutine copyCharArray2String(a,b,n)
         use KindParamModule, only : IntKind
         integer (kind=IntKind), intent(in) :: n
         character (len=1), intent(in) :: a(:)
         character (len=*), intent(out) :: b
      end subroutine copyCharArray2String
   end interface
!
!  ===================================================================
!  decoding imsgbuf................................................
!  ===================================================================
   jmt0   =imsgbuf(1)
   ntchg  =imsgbuf(2)
   numc   =imsgbuf(3)
!
   nr_old = max(jmt0,ntchg)
   allocate( x_mesh_old(nr_old), r_mesh_old(nr_old) )
   allocate( vrold(nr_old) )
   if ( ntchg /=0 ) then
      allocate( chgold(nr_old) )
      chgold(1:ntchg) = ZERO
   endif
!
   vr(1:nr) = ZERO
   vrold(1:nr_old) = ZERO
!
!  ===================================================================
!  decoding cmsgbuf................................................
!  -------------------------------------------------------------------
   call copyCharArray2String(cmsgbuf,header,80)
   call copyCharArray2String(cmsgbuf(81:160),jtitle,80)
!  -------------------------------------------------------------------
   do j=1,numc
!     ----------------------------------------------------------------
      call copyCharArray2String(cmsgbuf(160+(j-1)*40+1:160+j*40),temp_buf,40)
!     ----------------------------------------------------------------
      read(temp_buf,'(3i5,1e20.13,1a5)')nc(j),lc(j),kc(j),ec(j),lst(j)
   enddo
!
!  ===================================================================
!  decoding fmsgbuf................................................
!  ===================================================================
   alat  =fmsgbuf(1)
   efpot =fmsgbuf(2)
   vdif  =fmsgbuf(3)
   azed  =fmsgbuf(4)
   zcd   =fmsgbuf(5)
   xvalws=fmsgbuf(6)
   vzero =fmsgbuf(7)
   xst   =fmsgbuf(8)
   xmt   =fmsgbuf(9)
   do ir=1,jmt0
      if (ztotss > TEN2m8) then
         vrold(ir)=max(fmsgbuf(ir+9),-two*ztotss)
      else
         vrold(ir)=fmsgbuf(ir+9)
      endif
   enddo
   n=9+jmtmax
   do ir=1,ntchg
      chgold(ir)=fmsgbuf(ir+n)
   enddo
!
!  ===================================================================
   if(iprint.ge.0) then
      write(6,'(/,'' POTREDG:: jmt0''   ,t40,''='',i5)') jmt0
      write(6,'(''           ntchg''  ,t40,''='',i5)') ntchg
      write(6,'(''           numc''   ,t40,''='',i5)') numc
      write(6,'(''           xvalws'',t40,''='',e20.13)') xvalws
      write(6,'(''           vdif''   ,t40,''='',e20.13)') vdif
      write(6,'(''           alat''   ,t40,''='',e20.13)') alat
      write(6,'(''           azed''   ,t40,''='',e20.13)') azed
      write(6,'(''           zcd''    ,t40,''='',e20.13)') zcd
      write(6,'(''           efpot''  ,t40,''='',e20.13)') efpot
      write(6,'(''           xmt''    ,t40,''='',e20.13)') xmt
      write(6,'(''           xst''    ,t40,''='',e20.13)') xst
      write(6,'(''           vzero''  ,t40,''='',e20.13)') vzero
      do j=1,numc
         write(6,'(11x,''nc,lc,kc,ec,lst'',3i3,e12.4,1x,a5)')         &
                         nc(j),lc(j),kc(j),ec(j),lst(j)
      enddo
   endif
!
!  ===================================================================
!  check the parameters............................................
!  ===================================================================
   if (abs(ztotss-azed) > TEN2m6) then
!     ----------------------------------------------------------------
      call ErrorHandler('potredg','ztotss <> azed',ztotss,azed)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  generate grid for potential read-in ............................
!  ===================================================================
   hh=(xmt-xst)/real(jmt0-1,kind=RealKind)
   do j=1,nr_old
      x_mesh_old(j)=xst+(j-1)*hh
   enddo
   do j=1,nr_old
      r_mesh_old(j)=exp(x_mesh_old(j))
   enddo
   do j=1,jmt0
      vrold(j)=vrold(j)-vzero*r_mesh_old(j)
   enddo
!  ===================================================================
!  Print out of old grid and potential if needed...................
!  ===================================================================
   if(iprint.ge.0) then
      write(6,'(/,''  POTREDG: Old Potential Grid Parameters'')')
      write(6,'(''  POTREDG:    jmt0:'',i4)') jmt0
      write(6,'(''  POTREDG:    alat:'',f10.5)') alat
      write(6,'(''  POTREDG:       h:'',e23.14)') hh
      write(6,'(''  POTREDG: xst,xmt:'',2e23.14)')x_mesh_old(1),x_mesh_old(jmt0)
      write(6,'(''  POTREDG: rst,rmt:'',2e23.14)')r_mesh_old(1),r_mesh_old(jmt0)
      write(6,'(''  POTREDG: vst,vmt:'',2e23.14)') vrold(1),vrold(jmt0)
   endif
!
!  ===================================================================
!  interpolate chg density and cor density onto rsw mesh
!  ===================================================================
   performed_Interpolation = .false.
   if (jmt /= jmt0 .or. abs(log(r_mesh(jmt))-x_mesh_old(jmt0)) > TEN2m8 .or. &
       abs(log(r_mesh(1))-x_mesh_old(1)) > TEN2m8) then 
      do j=1,jmt
!        vr(j)=ylag(x_mesh(j),x_mesh_old,vrold,0,3,jmt0,iex)
         call interp(r_mesh_old,vrold,nr_old,r_mesh(j),vr(j),dvr)
      enddo
      if (ztotss > TEN2m8) then
         j = 1
         rfac = vrold(1)/vr(1)
         do while (r_mesh(j) < r_mesh_old(1))
            if (vr(j) > vrold(1) .or. int(vr(j)) /= int(vrold(1))) then
               vr(j) = rfac * vr(j) ! This is to ensure vr -> -2*Z when r -> 0.
            endif
            j = j + 1
         enddo
      endif
      performed_Interpolation = .true.
   else
      vr = vrold
   endif
   if(ntchg.gt.0) then
      if (jws /= ntchg .or. abs(r_mesh(jws)-r_mesh_old(ntchg)) > TEN2m6) then 
!         do j=1,jmt
!           rhotot(j)=ylag(x_mesh(j),x_mesh_old,chgold,0,3,ntchg,iex)
!            call interp(r_mesh_old,chgold,jmt0,r_mesh(j),rhotot(j),dvr)
!         enddo
         do j=1,jws
!           rhotot(j)=ylag(x_mesh(j),x_mesh_old,chgold,0,3,ntchg,iex)
!           print*,'potredg:: ir=',j, ntchg, jws 
            call interp(r_mesh_old,chgold,ntchg,r_mesh(j),rhotot(j),dvr)
            if(rhotot(j) < ZERO) rhotot(j)=ZERO
         enddo
      else
         rhotot = chgold
      endif
   endif
!  ===================================================================
!  Print out of new grid and potential if needed...................
!  ===================================================================
   if(iprint.ge.0) then
      if (performed_Interpolation) then
         write(6,'(/,''  POTREDG: Interpolation is performed...'')')
         write(6,'(  ''  ======================================'')')
      else
         write(6,'(/,''  ======================================'')')
      endif
      write(6,'(''  POTREDG: New Potential Grid Parameters'')')
      write(6,'(''  POTREDG:     jmt:'',i4)') jmt
      write(6,'(''  POTREDG:    alat:'',f10.5)') alat
      write(6,'(''  POTREDG: rst,rmt:'',2e23.14)') r_mesh(1),r_mesh(jmt)
      write(6,'(''  POTREDG: vst,vmt:'',2e23.14)') vr(1),vr(jmt)
   endif
!
   deallocate( x_mesh_old, r_mesh_old, vrold )
!
   if(ntchg.gt.0) then
      deallocate(chgold)
   endif
!
   end subroutine potredg
