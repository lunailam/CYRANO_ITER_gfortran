        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:33:50 2012
        MODULE C06FRF__genmod
          INTERFACE 
            SUBROUTINE C06FRF(M,N,X,Y,INIT,TRIG,WORK,IFAIL)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: X(M*N)
              REAL(KIND=8) :: Y(M*N)
              CHARACTER(LEN=1) :: INIT
              REAL(KIND=8) :: TRIG(2*N)
              REAL(KIND=8) :: WORK(2*M*N)
              INTEGER(KIND=4) :: IFAIL
            END SUBROUTINE C06FRF
          END INTERFACE 
        END MODULE C06FRF__genmod
