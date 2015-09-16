        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:32:26 2012
        MODULE PLASMB__genmod
          INTERFACE 
            SUBROUTINE PLASMB(PLA,PLAS,IPOINT,K,V,IV1,LDV,SWITCH,ADD)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: IV1
              COMPLEX(KIND=8) :: PLA(6,6)
              COMPLEX(KIND=8) :: PLAS(6,6)
              INTEGER(KIND=4) :: IPOINT
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: V(IV1:IV1+LDV-1,28)
              INTEGER(KIND=4) :: SWITCH
              LOGICAL(KIND=4) :: ADD
            END SUBROUTINE PLASMB
          END INTERFACE 
        END MODULE PLASMB__genmod
