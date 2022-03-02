program testAdaptiveIntegration
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use MathParamModule, only : ZERO, TEN2m12, HALF, ONE, PI
!
   use MPPModule, only : initMPP, endMPP, MyPE, NumPEs, bcastMessage
!
   use GroupCommModule
!
   use AdaptIntegrationModule, only : initAdaptIntegration
   use AdaptIntegrationModule, only : endAdaptIntegration
   use AdaptIntegrationModule, only : getUniFormIntegration
   use AdaptIntegrationModule, only : getUniMeshIntegration
   use AdaptIntegrationModule, only : getGaussianIntegration
   use AdaptIntegrationModule, only : getRombergIntegration
!
   implicit   none
!
   integer (kind=IntKind) :: proc_dim(1)
   integer (kind=IntKind) :: info(42), ms(2), nx, nm, ns, nlor, i, j
!
   real (kind=RealKind) :: a, b, gamma, x0, y(22), ua, ub, ya, yb, ex, r
   real (kind=RealKind) :: s1, s2, s3, s4, s5, err1, err2
!
   complex (kind=CmplxKind) :: aux(1)
!
!  -------------------------------------------------------------------
   call initMPP()
!  -------------------------------------------------------------------
!
   proc_dim(1) = NumPEs
!  -------------------------------------------------------------------
   call initGroupComm()
   call createProcGrid(1,'MPI World',proc_dim)
!  -------------------------------------------------------------------
   write(6,'(a,i3)')'Number of groups = ',getNumGroups()  
   write(6,'(a,i3)')'Group ID = ',getGroupID('MPI World')
   write(6,'(3a,i3)')'For group title: ',trim(getGroupLabel(1)),', size = ',getNumPEsInGroup(1)
   write(6,'(2(a,i3))')'MyPE = ',MyPE,', My rank in the group = ',getMyPEinGroup(1)
!
   info = 0
   ms = 0
   y = ZERO
!
   if (MyPE == 0) then
      write(6,'(a,$)')'Enter the number of Lorentzian compositions for integration: '
      read(5,*)ms(1)
      if (ms(1) < 1) then
         call ErrorHandler('main','The number of Lorentzians < 1',ms(1))
      else if (ms(1) > 10) then
         call ErrorHandler('main','The number of Lorentzians > 10',ms(1))
      endif
!
      write(6,'(a,$)')'Enter the number of mesh points for integration: '
      read(5,*)ms(2)
      if (ms(2) < 3) then
         call ErrorHandler('main','The mesh number < 3',ms(2))
      endif
!
      write(6,'(a,$)')'Enter integration domain a, b: '
      read(5,*)y(1:2)
      if (y(1) >= y(2)) then
         call ErrorHandler('main','a >= b',y(1),y(2))
      endif
!
      j = 0
      do i = 1, 2*ms(1), 2
         j = j + 1
         write(6,'(a,i2,$)')'Enter gamma and x0 for Lorenzian #',j,': '
         read(5,*)y(i+2),y(i+3)
         if (y(2+i) < TEN2m12) then
            call ErrorHandler('main','gamma < 10^(-12)',y(2+i))
         endif
      enddo
   endif
!  -------------------------------------------------------------------
   call bcastMessage(ms,2,0)
   call bcastMessage(y,22,0)
!  -------------------------------------------------------------------
   nlor = ms(1)
   nx = ms(2)
   a = y(1); b = y(2)
!
!  -------------------------------------------------------------------
   call initAdaptIntegration(n=1,gid=1)
!  -------------------------------------------------------------------
!
   if (MyPE == 0) then
      info(1) = nlor
      j = 0
      do i = 1, 4*nlor, 4
         j = j + 2
         gamma = y(j+1)
         call getIntegerExpresion(gamma,info(i+1),info(i+2))
         r = info(i+2)
         write(6,'(a,f12.8)')'Reconstructed gamma = ',info(i+1)/r
         if (info(i+1)/r < TEN2m12) then
            call ErrorHandler('main','After reconstruction, gamma < 10^(-12)')
         endif
!        info(i+1) = floor(log10(gamma)); print *,'info(i+1) = ',info(i+1)
!        info(i+2) = int(gamma/10.0d0**info(i+1)); print *,'info(i+2) = ',info(i+2)
!        write(6,'(a,i2,a,i3,a)')'Reconstructed gamma = ',info(i+2),'.0 X 10^(',info(i+1),')'
!        write(6,'(a,d9.2)')     '                    = ',info(i+2)*10.0d0**info(i+1)
!        if (info(i+2)*10.0d0**info(i+1) < TEN2m12) then
!           call ErrorHandler('main','After reconstruction, gamma < 10^(-12)')
!        endif
         x0 = y(j+2)
         call getIntegerExpresion(x0,info(i+3),info(i+4))
         r = info(i+4)
         write(6,'(a,f10.5)')'Reconstructed x0    = ',info(i+3)/r
!        info(i+3) = int(x0); print *,'info(i+3) = ',info(i+3)
!        info(i+4) = int(abs(x0-info(i+3))*10.0d0); print *,'info(i+4) = ',info(i+4)
!        write(6,'(a,i5,a,i1)')'Reconstructed x0    = ',info(i+3),'.',info(i+4)
!        write(6,'(a,f10.2)')  '                    = ',info(i+3)+info(i+4)/10.0d0
      enddo
   endif
!  -------------------------------------------------------------------
   call bcastMessage(info,42,0)
!  -------------------------------------------------------------------
!
   ex = ZERO
   do i = 1, 2*nlor, 2
      gamma = y(i+2)
      x0 = y(i+3)
      if (abs(a-x0) < TEN2m12) then
         ya = HALF
      else
         ua = gamma*HALF/(a-x0)
         ya = atan(ua)/PI
      endif
      if (abs(b-x0) < TEN2m12) then
         yb = HALF
      else
         ub = gamma*HALF/(b-x0)
         yb = atan(ub)/PI
      endif
      if (ya < ZERO) then
         ya = ya + ONE
      endif
      if (yb < ZERO) then
         yb = yb + ONE
      endif
      ex = ex + (ya-yb)
   enddo
!
   s1 = getUniformIntegration(nx,a,b,info,Lorentzian,nm)
   s2 = getUniMeshIntegration(nx,a,b,info,Lorentzian,nm)
   s3 = getGaussianIntegration(nx,a,b,info,Lorentzian)
   s4 = getRombIntegral(nx,a,b,info,Lorentzian,err1,ns)
   s5 = getRombergIntegration(nx,a,b,info,Lorentzian,err2,nm)
!
   if (MyPE == 0) then
      write(6,'(/,a,d15.8)')'The exact result    = ', ex
      write(6,'(a,d15.8)')'Uniform integration = ', s1
      write(6,'(a,d15.8)')'UniMesh integration = ', s2
      write(6,'(a,d15.8)')'Gaussian integration= ', s3
      write(6,'(a,d15.8)')'Local Romb integral = ', s4
      write(6,'(a,i5)')   '  integration steps = ', ns
      write(6,'(a,i5)')   '  integration points= ', 2**(ns-1)+1
      write(6,'(a,d15.8)')'  integration error = ', err1
      write(6,'(a,d15.8)')'Romberg integration = ', s5
      write(6,'(a,i5)')   '  integration points= ', nm
      write(6,'(a,d15.8)')'  integration error = ', err2
   endif
!
!  -------------------------------------------------------------------
   call endAdaptIntegration()
!  -------------------------------------------------------------------
   call endGroupComm()
   call endMPP()
!  -------------------------------------------------------------------
!
!  ===================================================================
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function Lorentzian(info,x,aux,fac,t) result(Lz)
!  ===================================================================
!  return Lorentzian function
!         1/PI * (gamma*w) / [x^2+(gamma*w)^2]
!  if info(1) = 1, w = 1;
!             = 2, w = 1/2
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, HALF, ONE, PI, CZERO
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: i, n
!
   real (kind=RealKind), intent(in) :: x
   real (kind=RealKind), intent(in), optional :: fac
   real (kind=RealKind) :: Lz, g, w, r
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
!
   logical, intent(in), optional :: t
!
   w = HALF
   n = info(1)
!
   Lz = ZERO
   do i = 1, 4*n, 4
      r = info(i+2)
      g = info(i+1)/r*w
      r = info(i+4)
      x0 = info(i+3)/r
!     g = info(i+2)*10.0d0**(info(i+1))*w
      if (g < TEN2m12) then
         call ErrorHandler('Lorentzian','g < 10.0^(-12)')
      endif
      Lz = Lz + g/((x-x0)**2+g**2)
   enddo
   Lz = Lz/PI
!
   aux = CZERO
!
   end function Lorentzian
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRombIntegral(nstep,a,b,info,func,err,nm) result(romb)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : WarningHandler
!
   use MathParamModule, only : ZERO, TEN2m6, TEN2m8, FOURTH, ONE
!
   use InterpolationModule, only : PolyInterp
!
   integer (kind=IntKind), intent(in) :: info(:)
!
   INTEGER (kind=IntKind), intent(in) :: nstep
   INTEGER (kind=IntKind), intent(out) :: nm
   INTEGER (kind=IntKind), parameter :: inter = 5
   integer (kind=IntKind) :: j
!
   REAL (kind=RealKind), parameter :: EPS = TEN2m6
!
   REAL (kind=RealKind), intent(in) :: a,b
   REAL (kind=RealKind), intent(out) :: err
   REAL (kind=RealKind) :: romb
   REAL (kind=RealKind) :: h(nstep+1),s(nstep+1)
!
   complex (kind=CmplxKind) :: aux(1)
!
   interface
      function func(info,x,aux,xi,redundant) result(y)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: x
         real (kind=RealKind), intent(in), optional :: xi
         real (kind=RealKind) :: y
         complex (kind=CmplxKind), intent(out), target :: aux(:)
         logical, intent(in), optional :: redundant
      end function func
   end interface
!
   h(1) = ONE
   do j = 1, nstep
      nm = j
!     ----------------------------------------------------------------
      call trapzd(info,func,aux,a,b,s(j),j)
!     ----------------------------------------------------------------
      if (j .ge. inter) then
!        -------------------------------------------------------------
         call PolyInterp(inter,h(j-inter+1:),s(j-inter+1:),ZERO,romb,err)
!        -------------------------------------------------------------
         if (abs(err) .le. EPS*abs(romb)) then
            return
         endif
      endif
      s(j+1)=s(j)
      h(j+1)=FOURTH*h(j)
   enddo
!
!  -------------------------------------------------------------------
   call WarningHandler('getRombIntegral',                             &
                       'Number of Romberg steps > nstep',nstep)
!  -------------------------------------------------------------------
!
   END function getRombIntegral
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   SUBROUTINE trapzd(info,func,aux,a,b,s,n)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, HALF
!
   integer (kind=IntKind), intent(in) :: info(:)
!
   INTEGER (kind=IntKind), intent(in) :: n
   INTEGER (kind=IntKind) :: it,j
!
   REAL (kind=RealKind), intent(in) :: a,b
   REAL (kind=RealKind), intent(inout) :: s
   REAL (kind=RealKind) :: del,sum,tnm,x
!
   complex (kind=CmplxKind), intent(out) :: aux(:)
!
   interface
      function func(info,x,aux,xi,redundant) result(y)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: x
         real (kind=RealKind), intent(in), optional :: xi
         real (kind=RealKind) :: y
         complex (kind=CmplxKind), intent(out), target :: aux(:)
         logical, intent(in), optional :: redundant
      end function func
   end interface
!
   if (n == 1) then
      s=HALF*(b-a)*(func(info,a,aux)+func(info,b,aux))
   else
      it=2**(n-2)
      tnm=it
      del=(b-a)/tnm
      x=a+HALF*del
      sum=ZERO
      do j=1,it
         sum=sum+func(info,x,aux)
         x=x+del
      enddo
      s=HALF*(s+(b-a)*sum/tnm)
   endif
!
   END SUBROUTINE trapzd
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getIntegerExpresion(r,i1,i2)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, TEN2m10, ONE
!
   implicit none
!
   integer (kind=IntKind), intent(out) :: i1, i2
   integer (kind=IntKind) :: n, k, n10
!
   real (kind=RealKind), intent(in) :: r
   real (kind=RealKind) :: a, b, c, d
!
   a = abs(r)
!
   if (a < TEN2m12) then
      i1 = 0; i2 = 1
      return
   endif
!
   i1 = int(a)
   a = a - i1
!
   b = a
   c = ZERO
   d = ONE
   k = 0
   n10 = 0
   do while (abs(a-c) > TEN2m10 .and. d > TEN2m10)  
      b = b*10.0d0
      n = int(b)
      b = b - n
      d = d/10.0d0
!
      k = k*10 + n
      c = c + n*d
      n10 = n10 + 1
   enddo
!
   i2 = 1
   if (k > 0) then
      do while (mod(k,5) == 0)
         k = k/5
         i2 = i2*2
         n10 = n10 - 1
      enddo
      do while (mod(k,2) == 0)
         k = k/2
         i2 = i2*5
         n10 = n10 - 1
      enddo
      do n = 1, n10
         i2 = i2*10
      enddo
   endif
!
   i1 = k + i1*i2
!
   if (r < ZERO) then
      i1 = -i1
   endif
!
   end subroutine getIntegerExpresion
!  ===================================================================
end program testAdaptiveIntegration
