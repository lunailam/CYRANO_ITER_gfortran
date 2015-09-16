        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:36:41 2012
        MODULE MCRCR__genmod
          INTERFACE 
            SUBROUTINE MCRCR(NRA,NCA,A,LDA,NRB,NCB,B,LDB,NRC,NCC,C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: NRA
              INTEGER(KIND=4) :: NCA
              COMPLEX(KIND=4) :: A(LDA,*)
              INTEGER(KIND=4) :: NRB
              INTEGER(KIND=4) :: NCB
              COMPLEX(KIND=4) :: B(LDB,*)
              INTEGER(KIND=4) :: NRC
              INTEGER(KIND=4) :: NCC
              COMPLEX(KIND=4) :: C(LDC,*)
            END SUBROUTINE MCRCR
          END INTERFACE 
        END MODULE MCRCR__genmod
