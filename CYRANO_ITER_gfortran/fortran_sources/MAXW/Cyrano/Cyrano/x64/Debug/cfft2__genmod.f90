        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:38:54 2012
        MODULE CFFT2__genmod
          INTERFACE 
            SUBROUTINE CFFT2(INIT,IDIR,NPFFT,CR,WORK,CVC)
              INTEGER(KIND=4) :: NPFFT
              INTEGER(KIND=4) :: INIT
              INTEGER(KIND=4) :: IDIR
              COMPLEX(KIND=8) :: CR(0:NPFFT-1)
              REAL(KIND=8) :: WORK(4*NPFFT)
              COMPLEX(KIND=8) :: CVC(0:NPFFT-1)
            END SUBROUTINE CFFT2
          END INTERFACE 
        END MODULE CFFT2__genmod
