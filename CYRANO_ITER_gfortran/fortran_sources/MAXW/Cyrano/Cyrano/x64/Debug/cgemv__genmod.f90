        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:38:45 2012
        MODULE CGEMV__genmod
          INTERFACE 
            SUBROUTINE CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: A(LDA,*)
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=4) :: BETA
              COMPLEX(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE CGEMV
          END INTERFACE 
        END MODULE CGEMV__genmod
