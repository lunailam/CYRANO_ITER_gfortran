        !COMPILER-GENERATED INTERFACE MODULE: Sat Feb 04 20:32:10 2012
        MODULE WRP2D__genmod
          INTERFACE 
            SUBROUTINE WRP2D(NOFILE,IPLOT,NX,X,NY,Y,Z,LDZ,INCX,INCY,    &
     &TITLE,XUNIT,XLABEL,YUNIT,YLABEL,ZUNIT,ZLABEL,NTEXT,TEXT)
              INTEGER(KIND=4) :: NTEXT
              INTEGER(KIND=4) :: LDZ
              INTEGER(KIND=4) :: NY
              INTEGER(KIND=4) :: NX
              INTEGER(KIND=4) :: NOFILE
              INTEGER(KIND=4) :: IPLOT
              REAL(KIND=8) :: X(NX)
              REAL(KIND=8) :: Y(NY)
              REAL(KIND=8) :: Z(LDZ,*)
              INTEGER(KIND=4) :: INCX
              INTEGER(KIND=4) :: INCY
              CHARACTER(LEN=48) :: TITLE
              CHARACTER(LEN=10) :: XUNIT
              CHARACTER(LEN=48) :: XLABEL
              CHARACTER(LEN=10) :: YUNIT
              CHARACTER(LEN=48) :: YLABEL
              CHARACTER(LEN=10) :: ZUNIT
              CHARACTER(LEN=48) :: ZLABEL
              CHARACTER(LEN=48) :: TEXT(NTEXT)
            END SUBROUTINE WRP2D
          END INTERFACE 
        END MODULE WRP2D__genmod
