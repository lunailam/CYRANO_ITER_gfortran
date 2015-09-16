        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:37:49 2012
        MODULE DL2TCG__genmod
          INTERFACE 
            SUBROUTINE DL2TCG(N,A,LDA,FAC,LDFAC,IPVT,SCALE)
              INTEGER(KIND=4) :: LDFAC
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: FAC(LDFAC,*)
              INTEGER(KIND=4) :: IPVT(*)
              COMPLEX(KIND=8) :: SCALE(*)
            END SUBROUTINE DL2TCG
          END INTERFACE 
        END MODULE DL2TCG__genmod
