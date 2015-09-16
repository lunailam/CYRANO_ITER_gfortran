      SUBROUTINE LOCTRA(IELE, ELEM)

      IMPLICIT NONE
      LOGICAL ELEM
      INTEGER IELE
      
C     Computes matrix 'BC' of rot.rot term at gauss points of element
C     (ELEM=.T.) or at local reduced abscissa Y (ELEM=.F.)
C     For ELEM=.T., also at left boundary of element (useful at magn. axis)
C     The pointer INTAB is assumed set to the index of first Gauss point of
C     current element.
C     Trigono. data of magnetic field angle provided by comma2.
C     Only valid for circular concentric equilibrium.

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'comma2.copy'
      include 'commag.copy'
      include 'comrot.copy'
      include 'comfin.copy'
      include 'comin2.copy'
      include 'comswe.copy'
      include 'comphy.copy'

      INTEGER IG1, IG2, I, J, INTABL
      DOUBLE PRECISION T
      COMPLEX*16 CT
C     RESULTS ARE SENT THROUGH COMROT

      T = SQRT2I
      CT = DCMPLX(T,0.D0)

      CALL ZSET( 360, CZERO, BC1, 1 )

      IF (ELEM) THEN
      IG1 = 0
      IG2 = NGAUSS
      ELSE
      IG1 = 1
      IG2 = 1
      END IF

      DO 100 IG = IG1, IG2
      INTABL = INTAB + IG - 1
      SI = SITAB(INTABL)
      CO = COTAB(INTABL)
      SI1 = SI1TAB(INTABL)
      CO1 = CO1TAB(INTABL)

      BC(1,1,IG) = CT
      BC(1,3,IG) = CT
      BC(2,2,IG) = CT
      BC(2,4,IG) = CT

      BC(3,1,IG) =    DCMPLX(0.D0,-T*CO)
      BC(3,3,IG) =  - BC(3,1,IG)
      BC(3,5,IG) =    DCMPLX(SI,0.D0)

      BC(4,1,IG) = DCMPLX(0.D0,-T*CO1)
      BC(4,2,IG)=BC(3,1,IG)
      BC(4,3,IG)=-BC(4,1,IG)
      BC(4,4,IG)=-BC(4,2,IG)
      BC(4,5,IG)=DCMPLX(SI1,0.D0)
      BC(4,6,IG)=DCMPLX(SI,0.D0)

      BC(5,1,IG) = DCMPLX(0.D0,T*SI)
      BC(5,3,IG) =-BC(5,1,IG)
      BC(5,5,IG) = DCMPLX(CO,0.D0)

      BC(6,1,IG) = DCMPLX(0.D0,T*SI1)
      BC(6,2,IG) = BC(5,1,IG)
      BC(6,3,IG) =-BC(6,1,IG)
      BC(6,4,IG) =-BC(6,2,IG)
      BC(6,5,IG) = DCMPLX(CO1,0.D0)
      BC(6,6,IG) = DCMPLX(CO,0.D0)

      DO 1 I=1,6
      DO 1 J=1,6
      BCH(I,J,IG)=DCONJG(BC(J,I,IG))
  1   CONTINUE

  100 CONTINUE

      RETURN
      END

C
C****************************************
C
C         TEST LOCTRA
C
CC      INCLUDE 'COMMA2.COPY'
CC      INCLUDE 'COMMAG.COPY'
CC      INCLUDE 'COMFIN.COPY'
CC      INCLUDE 'COMROT.COPY'
CC      INCLUDE 'COMIN2.COPY'
CC      INCLUDE 'NAMMAG.COPY'
C
CC      DOUBLE PRECISION H
CC      NGAUSS=4
CC      FX0(1)=0.D0
CC      H=0.3D0
CC      READ(51,NAMMAG)
C
CC      CALL BAFAGA
CC      DO 1 I=1,5
CC      FL(I)=H
CC      CALL LOCTRA(I)
CC      FX0(I+1)=H*I
CC  1   CONTINUE
CC      STOP
CC      END
