SUBROUTINE ZGETF2_nopivot( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
  INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
  INTEGER             IPIV( * )
  COMPLEX(kind(0.d0)) A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
  COMPLEX (KIND(0.D0)) ONE, ZERO
  PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
  DOUBLE PRECISION   SFMIN
  INTEGER            I, J, JP
!     ..
!     .. External Functions ..
  DOUBLE PRECISION   DLAMCH
  INTEGER            IZAMAX
  EXTERNAL           DLAMCH, IZAMAX
!     ..
!     .. External Subroutines ..
  EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZGETF2', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 )  RETURN
!
!     Compute machine safe minimum
!
  SFMIN = DLAMCH('S') 
!
  DO J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
!         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
     JP = J
     IPIV( J ) = JP
     IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
!            IF( JP.NE.J )
!     $         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
        IF( J.LT.M ) THEN
           IF( ABS(A( J, J )) .GE. SFMIN ) THEN
              CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
           ELSE
              DO I = 1, M-J
                 A( J+I, J ) = A( J+I, J ) / A( J, J )
              END DO
           END IF
        END IF
!
     ELSE IF( INFO.EQ.0 ) THEN
!
        INFO = J
     END IF
!
     IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
        CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),LDA, A( J+1, J+1 ), LDA )
     END IF
  END DO
  RETURN
!
!     End of ZGETF2
!
END SUBROUTINE ZGETF2_nopivot
