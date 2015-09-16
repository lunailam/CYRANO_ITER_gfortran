        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:38:28 2012
        MODULE F06PAF__genmod
          INTERFACE 
            SUBROUTINE F06PAF(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE F06PAF
          END INTERFACE 
        END MODULE F06PAF__genmod
