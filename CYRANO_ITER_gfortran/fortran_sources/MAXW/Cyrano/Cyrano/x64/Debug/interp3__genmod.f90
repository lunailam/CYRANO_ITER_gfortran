        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:36:47 2012
        MODULE INTERP3__genmod
          INTERFACE 
            SUBROUTINE INTERP3(XTAB,INCX,REFTAB,INCREF,NREF,XABS,TAB,   &
     &INCTAB,NTAB)
              INTEGER(KIND=4) :: NTAB
              INTEGER(KIND=4) :: INCTAB
              INTEGER(KIND=4) :: NREF
              INTEGER(KIND=4) :: INCREF
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: XTAB(1+(NREF-1)*INCX)
              REAL(KIND=8) :: REFTAB(1+(NREF-1)*INCREF)
              REAL(KIND=8) :: XABS(NTAB)
              REAL(KIND=8) :: TAB(1+(NTAB-1)*INCTAB)
            END SUBROUTINE INTERP3
          END INTERFACE 
        END MODULE INTERP3__genmod
