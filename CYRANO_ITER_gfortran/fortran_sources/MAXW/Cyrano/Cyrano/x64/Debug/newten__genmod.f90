        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:33:00 2012
        MODULE NEWTEN__genmod
          INTERFACE 
            SUBROUTINE NEWTEN(CASE,OMC,OMP,OM,KPERP,RKPAR,VELEC,VPERP,  &
     &VPAR,V0,UDRIFT,KNEW)
              CHARACTER(LEN=5) :: CASE
              REAL(KIND=8) :: OMC
              REAL(KIND=8) :: OMP
              REAL(KIND=8) :: OM
              COMPLEX(KIND=8) :: KPERP
              REAL(KIND=8) :: RKPAR
              REAL(KIND=8) :: VELEC
              REAL(KIND=8) :: VPERP
              REAL(KIND=8) :: VPAR
              REAL(KIND=8) :: V0
              REAL(KIND=8) :: UDRIFT
              COMPLEX(KIND=8) :: KNEW(3,3)
            END SUBROUTINE NEWTEN
          END INTERFACE 
        END MODULE NEWTEN__genmod
