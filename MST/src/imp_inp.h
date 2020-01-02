!==========================================================================
!   * ==============       INPUT PARAM for VP-iteg     ==================*
!==========================================================================
!  #1-->MNV ; integer*4
!        maximum number of vertices for atom
!
!  #2-->MNF ; integer*4
!        maximum number of faces for atom

!  #3-->MNqp ; integer*4
!        Maximum number of quadrature points.
!        e.g.: For volume, it takes 2 points to get converge.
!              For volume-integral it takes at least 13 points
!              to get converge usually.

!  #4-->Nf ; integer*4
!        Number of data points. Assume each atom have the same number
!        of data points, no matter the size of atom.

!  #5-->npol ; integer*4
!   The order of the polynomial, which use in polynomial fitting

!  #6-->poly_type ; integer*4
!   Type of the polynomial which use in polynomial fitting.
!   poly_type=1 --> simple polynomial
!   poly_type=2 --> Laguerre polynomial
!   poly_type=3 --> Legendre polynomial
!   poly_type=4 --> Chebyshev polynomial
!   poly_type=5 --> ?

!  #7-->Ngauss ; integer*4
!        Number of Gauss quadrature points used in the calc.

!  #8-->NNN ; integer*4
!        Number of nearest neighbors used (i.e. R in Ewald sum), system dependent

!  #9-->Npt ; integer*4
!        Number of radial points along each high symm direc. while 
!        calculating V(r).
!==========================================================================
	 integer*4 MNV,MNF,MNqp,Nf,npol,poly_type
	 integer*4 Ngauss,NNN,Npt
! These param are fixed 
         parameter (MNV=100,MNF=50,MNqp=25,Nf=4000,npol=80,poly_type=2)
! These param might need to vary for accuracy purpose, and with system 
         parameter (Ngauss=10,NNN=86,Npt=30)
