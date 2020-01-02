!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function erfc(xx)
! ====================================================================
  use KindParamModule, only : IntKind, RealKind
!
  implicit none
!     
! sandia mathematical program library
! applied mathematics division 2613
! sandia laboratories
! albuquerque, new mexico  87185
! control data 6600/7600  version 7.2  may 1978
! $Id: erfc.x,v 1.1 89/03/11 19:09:54 sverre Exp $
!
! $Log:	erfc.x,v $
! Revision 1.1  89/03/11  19:09:54  sverre
! Initial revision
!     
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *                 issued by sandia laboratories                     *
! *                   a prime contractor to the                       *
! *                united states department of energy                 *
! * * * * * * * * * * * * * * * notice  * * * * * * * * * * * * * * * *
! * this report was prepared as an account of work sponsored by the   *
! * united states government.  neither the united states nor the      *
! * united states department of energy nor any of their employees,    *
! * nor any of their contractors, subcontractors, or their employees  *
! * makes any warranty, express or implied, or assumes any legal      *
! * liability or responsibility for the accuracy, completeness or     *
! * usefulness of any information, apparatus, product or process      *
! * disclosed, or represents that its use would not infringe          *
! * owned rights.                                                     *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * the primary document for the library of which this routine is     *
! * part is sand77-1441.                                              *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     
! written by j.e. vogel from approximations derived by w.j. cody .
!     
! abstract
!     
! .  erfc(x) computes 2.0/sqrt(pi) times the integral from x to
! .  infinity of exp(-x**2). this is done using rational approx-
! .  imations.  eleven correct significant figures are provided.
!     
! description of parameters
!     
! .  x may be any real value
!     
! erfc is documented completely in sc-m-70-275.
!     
      integer (kind=IntKind) :: i
!
      real (kind=RealKind), intent(in) :: xx
      real (kind=RealKind) :: erfc
!
      real (kind=RealKind) :: x, x2, xi2, r, a
      real (kind=RealKind) :: p1(4),q1(4),p2(6),q2(6),p3(4),q3(4)
!
      data (p1(i),i=1,4)/ 242.6679552305318D0,  21.97926161829415D0,  &
                          6.996383488619136D0, -3.560984370181539D-2/
      data (q1(i),i=1,4)/ 215.0588758698612D0,  91.16490540451490D0,  &
                          15.08279763040779D0,  1.0D0/
      data (p2(i),i=1,6)/   22.898992851659D0,    26.094746956075D0,  &
                            14.571898596926D0,    4.2677201070898D0,  &
                           0.56437160686381D0,  -6.0858151959688D-6/
      data (q2(i),i=1,6)/   22.898985749891D0,   51.933570687552D0,   &
                            50.273202863803D0,   26.288795758761D0,   &
                            7.5688482293618D0,  1.0D0/
      data (p3(i),i=1,4)/-1.21308276389978D-2, -0.1199039552681460D0, &
                         -0.243911029488626D0, -3.24319519277746D-2/
      data (q3(i),i=1,4)/ 4.30026643452770D-2,  0.489552441961437D0,  &
                           1.43771227937118D0,  1.0D0/
      real (kind=RealKind), parameter :: sqpi=0.564189583547756D0
!     
      x=abs(xx)
      x2=x*x
!     
      if (xx .lt. -6.D0) then
         erfc = 2.D0
         return
      else if (xx .gt. 25.8D0) then
         erfc = 0.D0
         return
      else if (x .gt. 4.D0) then
         xi2 = 1.D0/x2
         r = xi2*(p3(1)+xi2*(p3(2)+xi2*(p3(3)+xi2*p3(4))))/         &
              (q3(1)+xi2*(q3(2)+xi2*(q3(3)+xi2*q3(4))))
         a = exp(-x2)*(sqpi+r)/x
         if (xx .lt. 0.D0) then
            erfc = 2.D0-a
         else
            erfc = a
         end if
         return
      else if (x .gt. 0.46875D0) then
         a = exp(-x2)
         a = a*(p2(1)+x*(p2(2)+x*(p2(3)+x*(p2(4)+x*(p2(5)+x*p2(6))))))
         a = a/(q2(1)+x*(q2(2)+x*(q2(3)+x*(q2(4)+x*(q2(5)+x*q2(6))))))
         if (xx .lt. 0.D0) then
            erfc = 2.D0-a
         else
            erfc = a
         end if
         return
      else
         a= x*(p1(1)+x2*(p1(2)+x2*(p1(3)+x2*p1(4))))
         a = a/(q1(1)+x2*(q1(2)+x2*(q1(3)+x2*q1(4))))
         if (xx .lt. 0.D0) then
            erfc = 1.D0+a
         else
            erfc = 1.D0-a
         end if
         return
      end if
  end function erfc
