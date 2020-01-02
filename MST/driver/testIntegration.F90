!  *******************************************************************
!  *   TEST INTEGRATION SCHEMES                                      *
!  *******************************************************************
program testIntegration
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO, ONE, CONE, SQRTm1, TWO, THIRD, HALF
   use IntegrationModule, only : calIntegration
!
   integer, parameter :: n = 2000
   integer :: i
!
   real (kind=RealKind) :: h
   real (kind=RealKind) :: r(0:n), f(0:n), fint(0:n)
   real (kind=RealKind) :: rtmp(0:n), ftmp(0:n), ftint(0:n)
   real (kind=RealKind) :: s, st, exr, err
!
   complex (kind=CmplxKind) :: cf(0:n), cfint(0:n)
   complex (kind=CmplxKind) :: cftmp(0:n), cftint(0:n)
   complex (kind=CmplxKind) :: cs, cst, exc
!
   h = ONE/real(n,RealKind)
   do i=0,n
      r(i) = h*i                   ! define mesh between [0, 1]
      f(i) = r(i)*exp(r(i)*r(i))   ! define a real function: r*exp(r^2)
      cs = SQRTm1*r(i)
      cf(i) = r(i)*exp(cs)      ! define a complex function: r*exp(i*r)
   enddo
!
!  ===================================================================
!  use our traditional routine
!  ===================================================================
   call calIntegration(0,n+1,r(0:n),f(0:n),fint(0:n))
   call calIntegration(0,n+1,r(0:n),cf(0:n),cfint(0:n))
   s = fint(n)
   cs = cfint(n)
!
!  ===================================================================
!  define rtmp, ftmp, and cftmp for testing XGZ's integration code
!  ===================================================================
   do i=0,n
      rtmp(i) = r(i)
      ftmp(i) = f(i)
      cftmp(i) = cf(i)
   enddo
   call calIntegration(n+1,rtmp(0:n),ftmp(0:n),ftint(0:n),0)
   call calIntegration(n+1,rtmp(0:n),cftmp(0:n),cftint(0:n),0)
   st = ftint(n)
   cst = cftint(n)
!
   write(6,'(/,a)')'TEST 1'
   exr = HALF*(exp(ONE)-ONE)
   write(6,'(a,d19.11)')'exact result       = ',exr
   write(6,'(a,d19.11,a,d15.8)')'traditional scheme = ',s,', Error = ',abs(s-exr)
   write(6,'(a,d19.11,a,d15.8)')'XYG scheme         = ',st,', Error = ',abs(st-exr)
!
   exc = (CONE-SQRTm1)*exp(SQRTm1)-CONE
   write(6,'(a,2d19.11)')'exact result       = ',exc
   write(6,'(a,2d19.11,a,d15.8)')'traditional scheme = ',cs,', Error = ',abs(cs-exc)
   write(6,'(a,2d19.11,a,d15.8)')'XGZ scheme         = ',cst,', Error = ',abs(cst-exc)
!
!  ===================================================================
!  use our traditional routine
   call calIntegration(0,n,r(1:n),f(1:n),fint(1:n))
   call calIntegration(0,n,r(1:n),cf(1:n),cfint(1:n))
   s = fint(n) + f(1)*h
   cs = cfint(n) + cf(1)*h
!
!  ===================================================================
!  call XGZ's routine
!  ===================================================================
   do i=1,n
      rtmp(i) = sqrt(r(i))
      ftmp(i) = TWO*f(i)/r(i)
      cftmp(i) = TWO*cf(i)/r(i)
   enddo
   call calIntegration(n,rtmp(1:n),ftmp(1:n),ftint(1:n),3)
   call calIntegration(n,rtmp(1:n),cftmp(1:n),cftint(1:n),3)
   st = ftint(n) + f(1)*h
   cst = cftint(n) + cf(1)*h
!
   write(6,'(/,a)')'TEST 2'
   exr = HALF*(exp(ONE)-ONE)
   write(6,'(a,d19.11)')'exact result       = ',exr
   write(6,'(a,d19.11,a,d15.8)')'traditional scheme = ',s,', Error = ',abs(s-exr)
   write(6,'(a,d19.11,a,d15.8)')'XYG scheme         = ',st,', Error = ',abs(st-exr)
!
   exc = (CONE-SQRTm1)*exp(SQRTm1)-CONE
   write(6,'(a,2d19.11)')'exact result       = ',exc
   write(6,'(a,2d19.11,a,d15.8)')'traditional scheme = ',cs,', Error = ',abs(cs-exc)
   write(6,'(a,2d19.11,a,d15.8)')'XGZ scheme         = ',cst,', Error = ',abs(cst-exc)
!
!  ===================================================================
!  Another test
!  ===================================================================
   do i=1,n
      f(i) = ONE/sqrt(r(i))          ! define a real function: 1/sqrt(r)
      cf(i) = CONE/sqrt(r(i)+SQRTm1) ! define a complex function: 1/sqrt(r+i)
   enddo
!
   call calIntegration(0,n,r(1:n),f(1:n),fint(1:n))
   call calIntegration(0,n,r(1:n),cf(1:n),cfint(1:n))
   s = fint(n)
   cs = cfint(n)
!
   do i=1,n
      rtmp(i) = sqrt(r(i))
      ftmp(i) = TWO*f(i)/r(i)**2
      cftmp(i) = TWO*cf(i)/r(i)**2
   enddo
   call calIntegration(n,rtmp(1:n),ftmp(1:n),ftint(1:n),5)
   call calIntegration(n,rtmp(1:n),cftmp(1:n),cftint(1:n),5)
   st = ftint(n)
   cst = cftint(n)
!
   write(6,'(/,a)')'TEST 3'
   exr = TWO*(ONE-sqrt(r(1)))
   write(6,'(a,d19.11)')'exact result       = ',exr
   write(6,'(a,d19.11,a,d15.8)')'traditional scheme = ',s,', Error = ',abs(s-exr)
   write(6,'(a,d19.11,a,d15.8)')'XYG scheme         = ',st,', Error = ',abs(st-exr)
!
   exc = TWO*(sqrt(CONE+SQRTm1)-sqrt(r(1)+SQRTm1))
   write(6,'(a,2d19.11)')'exact result       = ',exc
   write(6,'(a,2d19.11,a,d15.8)')'traditional scheme = ',cs,', Error = ',abs(cs-exc)
   write(6,'(a,2d19.11,a,d15.8)')'XGZ scheme         = ',cst,', Error = ',abs(cst-exc)
!
!  ===================================================================
!  One more test
!  ===================================================================
   do i=0,n
      f(i) = r(i)*r(i)          ! define a real function: r**2
      cf(i) = r(i)*r(i)         ! define a complex function r**2
   enddo
!
   call calIntegration(0,n+1,r(0:n),f(0:n),fint(0:n))
   call calIntegration(0,n+1,r(0:n),cf(0:n),cfint(0:n))
   s = fint(n)
   cs = cfint(n)
!
   do i=0,n
      f(i) = ONE          ! define a real function: 1/sqrt(r)
      cf(i) = CONE ! define a complex function: 1/sqrt(r+i)
   enddo
!
   call calIntegration(n+1,r(0:n),f(0:n),ftint(0:n),2)
   call calIntegration(n+1,r(0:n),cf(0:n),cftint(0:n),2)
   st = ftint(n)
   cst = cftint(n)
!
   write(6,'(/,a)')'TEST 4'
   exr = THIRD
   write(6,'(a,d19.11)')'exact result       = ',exr
   write(6,'(a,d19.11,a,d15.8)')'traditional scheme = ',s,', Error = ',abs(s-exr)
   write(6,'(a,d19.11,a,d15.8)')'XYG scheme         = ',st,', Error = ',abs(st-exr)
!
   exc = THIRD
   write(6,'(a,2d19.11)')'exact result       = ',exc
   write(6,'(a,2d19.11,a,d15.8)')'traditional scheme = ',cs,', Error = ',abs(cs-exc)
   write(6,'(a,2d19.11,a,d15.8)')'XGZ scheme         = ',cst,', Error = ',abs(cst-exc)
!
   stop 'OK'
end program testIntegration
