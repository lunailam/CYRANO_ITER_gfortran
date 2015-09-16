        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:37:56 2012
        MODULE DL2CCG__genmod
          INTERFACE 
            SUBROUTINE DL2CCG(N,A,LDA,FAC,LDFAC,IPVT,RCOND,Z)
              INTEGER(KIND=4) :: LDFAC
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: FAC(LDFAC,*)
              INTEGER(KIND=4) :: IPVT(*)
              REAL(KIND=8) :: RCOND
              COMPLEX(KIND=8) :: Z(*)
            END SUBROUTINE DL2CCG
          END INTERFACE 
        END MODULE DL2CCG__genmod
