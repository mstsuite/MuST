program tst_findFit
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : CZERO, ZERO, TEN2m6, TEN2m4, HALF, ONE
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   character (len=1) :: text
!
   integer (kind=IntKind), parameter :: lmax_rho = 8
   integer (kind=IntKind), parameter :: jmax_rho = (lmax_rho+1)*(lmax_rho+2)/2
   integer (kind=IntKind) :: jl_flag(jmax_rho)
   integer (kind=IntKind) :: status
   integer (kind=IntKind) :: l, m, jl, NumRs, ir, ic, nc, atype, ir0(jmax_rho)
!
   real (kind=RealKind) :: r0, r, r2, r3, rl
   real (kind=RealKind), parameter :: epsilon = TEN2m6
   real (kind=RealKind), allocatable :: r_mesh(:), sqrtr(:)
   real (kind=RealKind), allocatable :: rhor_l(:,:), rhoi_l(:,:)
   real (kind=RealKind), allocatable :: rhort(:), rhoit(:)
   real (kind=RealKind), allocatable :: drhor(:), ddrhor(:)
   real (kind=RealKind), allocatable :: drhoi(:), ddrhoi(:)
!
   complex (kind=CmplxKind) :: rhol_r0,drhol_r0,ddrhol_r0,rhoc_fit,rhoq_fit
   complex (kind=CmplxKind) :: ac(jmax_rho),bc(jmax_rho),cc(jmax_rho), &
                               dc(jmax_rho)
   complex (kind=CmplxKind) :: aq(jmax_rho),bq(jmax_rho),cq(jmax_rho)
!
   open(unit=10,file='rhol.dat',status='old',form='formatted')
!
   NumRs = -1
   LOOP_do: do
      read(10,'(a)',iostat=status)text
      if (status  < 0) then
         exit LOOP_do
      else
         NumRs = NumRs + 1
      endif
   enddo LOOP_do
!
   rewind(unit=10)
!
   read(10,'(153i4)')(jl_flag(jl),jl=1,jmax_rho)
   nc = 0
   do jl = 1, jmax_rho
      if (jl_flag(jl) == 1) then
         nc = nc + 1
      endif
   enddo
!
   print *,'NumRs = ',NumRs,',  nc = ',nc
!
   allocate( r_mesh(NumRs), sqrtr(NumRs) )
   allocate( rhor_l(NumRs,jmax_rho), rhoi_l(NumRs,jmax_rho) )
   allocate( drhor(NumRs), ddrhor(NumRs), rhort(nc) )
   allocate( drhoi(NumRs), ddrhoi(NumRs), rhoit(nc) )
   do ir = 1, NumRs
      read(10,'(i4,1x,d16.8,153(2x,d16.8,1x,d16.8))') &
                atype,r_mesh(ir),(rhort(ic),rhoit(ic),ic=1,nc)
      sqrtr(ir) = sqrt(r_mesh(ir))
      ic = 0
      do jl = 1, jmax_rho
         if (jl_flag(jl) == 1) then
            ic = ic + 1
            rhor_l(ir,jl) = rhort(ic)
            rhoi_l(ir,jl) = rhoit(ic)
         else
            rhor_l(ir,jl) = ZERO
            rhoi_l(ir,jl) = ZERO
         endif
      enddo
   enddo
!
   close(unit=10)
!
   r0 = ZERO
   do while (r0 <= TEN2m6)
      write(6,'(a,$)')'Enter fitting radius (> 0) = '
      read(5,*)r0
   enddo
!
!  -------------------------------------------------------------------
   call hunt(NumRs,r_mesh,r0,ir)
!  -------------------------------------------------------------------
   print *,'actual r0 = ',r_mesh(ir)
   ir0(1:jmax_rho) = ir
!
   jl = 0
   do l = 0, lmax_rho
      do m = 0, l
         jl = jl + 1
         if (jl_flag(jl) == 1) then
!           ----------------------------------------------------------
            call newder(rhor_l(1:NumRs,jl),drhor,sqrtr,NumRs)
!           ----------------------------------------------------------
            do ir = 1, NumRs
               drhor(ir) = HALF*drhor(ir)/sqrtr(ir)
            enddo
!           ----------------------------------------------------------
            call newder(drhor,ddrhor,sqrtr,NumRs)
!           ----------------------------------------------------------
            do ir = 1, NumRs
               ddrhor(ir) = HALF*ddrhor(ir)/sqrtr(ir)
            enddo
!           ----------------------------------------------------------
            call newder(rhoi_l(1:NumRs,jl),drhoi,sqrtr,NumRs)
!           ----------------------------------------------------------
            do ir = 1, NumRs
               drhoi(ir) = HALF*drhoi(ir)/sqrtr(ir)
            enddo
!           ----------------------------------------------------------
            call newder(drhoi,ddrhoi,sqrtr,NumRs)
!           ----------------------------------------------------------
            do ir = 1, NumRs
               ddrhoi(ir) = HALF*ddrhoi(ir)/sqrtr(ir)
            enddo
            ir = ir0(jl)
            do while (abs(drhor(ir))+abs(drhoi(ir)) <= epsilon) 
               ir = ir + 1
               if (ir >= NumRs) then
                  exit
               endif
            enddo
            if (ir >= NumRs) then
               ir = ir0(jl) - 1
               do while (abs(drhor(ir))+abs(drhoi(ir)) <= epsilon) 
                  ir = ir - 1
                  if (ir <= 0) then
                     call ErrorHandler('tst_findFit','no idea where to fit')
                  endif
               enddo
            endif
! 
            ir0(jl) = ir
            r0 = r_mesh(ir)
            print *,'ir0 = ',ir,',  r0 = ',r0,',  drho = ',drhor(ir),drhoi(ir)
            rhol_r0   = cmplx(rhor_l(ir,jl),rhoi_l(ir,jl),kind=CmplxKind)
            drhol_r0  = cmplx(drhor(ir),drhoi(ir),kind=CmplxKind)
            ddrhol_r0 = cmplx(ddrhor(ir),ddrhoi(ir),kind=CmplxKind)
!
!           ----------------------------------------------------------
            call findQuadraticFit(l,r0,rhol_r0,drhol_r0,aq(jl),bq(jl),cq(jl))
!           ----------------------------------------------------------
!
!           ----------------------------------------------------------
            call findCubicFit(l,r0,rhol_r0,drhol_r0,ddrhol_r0,        &
                              ac(jl),bc(jl),cc(jl),dc(jl))
!           ----------------------------------------------------------
         endif
      enddo
   enddo
!
   open(unit=11,file='rhol_fit.dat',status='unknown',form='formatted')
!
   write(11,'(1x,a14,$)')'r-mesh        '
   jl = 0
   ic = 0
   do l = 0, lmax_rho
      do m = 0, l
         jl = jl + 1
         if (jl_flag(jl) == 1) then
            ic = ic + 1
            write(11,'(2x,a1,i2,a1,i2,a1,1x,a7,$)')                   &
                  '(',l,',',m,')','Re Part'
            write(11,'(2x,a1,i2,a1,i2,a1,1x,a7,$)')                   &
                  '(',l,',',m,')','Im Part'
            write(11,'(2x,a1,i2,a1,i2,a1,1x,a7,$)')                   &
                  '(',l,',',m,')','Re Cfit'
            write(11,'(2x,a1,i2,a1,i2,a1,1x,a7,$)')                   &
                  '(',l,',',m,')','Im Cfit'
            write(11,'(2x,a1,i2,a1,i2,a1,1x,a7,$)')                   &
                  '(',l,',',m,')','Re Qfit'
            if (ic < nc) then
               write(11,'(2x,a1,i2,a1,i2,a1,1x,a7,$)')                &
                     '(',l,',',m,')','Im Qfit'
            else
               write(11,'(2x,a1,i2,a1,i2,a1,1x,a7)')                  &
                     '(',l,',',m,')','Im Qfit'
            endif
         endif
      enddo
   enddo
!
   do ir = 1, NumRs
      write(11,'(d15.8,''  '',$)')r_mesh(ir)
      r = r_mesh(ir); r2 = r*r; r3 = r*r*r
      jl = 0
      ic = 0
      rl = ONE
      do l = 0, lmax_rho
         do m = 0, l
            jl = jl + 1
            if (jl_flag(jl) == 1) then
               ic = ic + 1
               write(11,'(2(2x,d15.8),$)')rhor_l(ir,jl),rhoi_l(ir,jl)
               if (ir <= ir0(jl)) then
                  rhoc_fit = rl*(ac(jl)+bc(jl)*r+cc(jl)*r2+dc(jl)*r3)
                  rhoq_fit = rl*(aq(jl)+bq(jl)*r+cq(jl)*r2)
               else
                  rhoc_fit = cmplx(rhor_l(ir,jl),rhoi_l(ir,jl),kind=CmplxKind)
                  rhoq_fit = cmplx(rhor_l(ir,jl),rhoi_l(ir,jl),kind=CmplxKind)
               endif
               if (ic < nc) then
                  write(11,'(4(2x,d15.8),$)')rhoc_fit,rhoq_fit
               else
                  write(11,'(4(2x,d15.8))')  rhoc_fit,rhoq_fit
               endif
            endif
         enddo
         rl = rl*r
      enddo
   enddo
!
   close(unit=11)
!
end program tst_findFit
