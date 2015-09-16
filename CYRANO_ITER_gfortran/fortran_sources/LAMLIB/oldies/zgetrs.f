      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c
c     Alias for Cray single precision
c
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
      call CGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      RETURN
*
*     End of ZGETRS
*
      END
