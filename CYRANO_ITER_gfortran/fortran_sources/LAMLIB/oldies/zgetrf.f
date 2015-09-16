      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
c
c     Alias for Cray single precision
c
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
      call cgetrf( M, N, A, LDA, IPIV, INFO )
      RETURN
*
*     End of ZGETRF
*
      END
