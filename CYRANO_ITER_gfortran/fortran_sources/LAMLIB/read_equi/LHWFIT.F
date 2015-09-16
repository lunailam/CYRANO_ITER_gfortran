      SUBROUTINE LHWFIT (NPPF,RPPF,ZPPF,OPAR)
* Determination of Lao-Hirshman-Wieland representation of flux surface
*  Input:  RPPF(i), ZPPF(i) - flux surface (R,Z) points
*  Output: OPAR(1) - R0
*              (2) - Z0
*              (3) - r  (minor radius)
*              (4) - E  (elongation)
*              (5) - R2 (triangularity)
*              (6) - xi
*              (7) - s2
*              (8) - Error
cpl	  indicates P.L. comment (tentative interpretation of the code)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RPPF(*),ZPPF(*),OPAR(8)
      DIMENSION RBND(1000),ZBND(1000)
      DIMENSION AN(20),BN(20)
      DIMENSION THETA(1000),SINT(1000),COST(1000)
      COMMON /COMAB1/ THETA,SINT,COST,RBND,ZBND
C
      DATA TWOPI/6.28318 53071 79586 5D0/
C
      DO I=1,8
       OPAR(I) = -1.
      END DO
      DO I=1,3
         AN(I)=0.0
         BN(I)=0.0
      END DO
      NOUT=3
      E=0.0
C
C  COPY SLICE
C
      IF(RPPF(1).EQ.RPPF(NPPF) .AND. ZPPF(1).EQ.ZPPF(NPPF)) THEN
         NBND=NPPF-1
      ELSE
         NBND=NPPF
      END IF
C
      DO I=1,NBND
         RBND(I)=RPPF(I)
         ZBND(I)=ZPPF(I)
      END DO
      RBND(NBND+1)=RPPF(1)
      ZBND(NBND+1)=ZPPF(1)
C
C  CALC R0,Z0 FOR INITIAL GUESS
C
      TOTR=0.0
      TOTZ=0.0
      TOTL=0.0
      RMIN=RBND(1)
      ZMIN=ZBND(1)
      RMAX=RBND(1)
      ZMAX=ZBND(1)
      DO I=1,NBND
cpl	     element of poloidal length:
         XL=SQRT((RBND(I)-RBND(I+1))**2+(ZBND(I)-ZBND(I+1))**2)
cpl	     perimeter:
         TOTL=TOTL+XL
         TOTR=TOTR+(RBND(I)+RBND(I+1))*XL
         TOTZ=TOTZ+(ZBND(I)+ZBND(I+1))*XL
         RMIN=MIN(RMIN,RBND(I+1))
         RMAX=MAX(RMAX,RBND(I+1))
         ZMIN=MIN(ZMIN,ZBND(I+1))
         ZMAX=MAX(ZMAX,ZBND(I+1))
      END DO
      IF(TOTL.EQ.0.0) GOTO 9
cpl	  geom. centre of surface:
      R0=TOTR*0.5/TOTL
      Z0=TOTZ*0.5/TOTL
      E=(ZMAX-ZMIN)/(RMAX-RMIN)
CCC   WRITE(6,*) 'R0=',R0,' Z0=',Z0,' ELLIPTICITY=',E
C
C  CALC THETA FOR INITIAL GUESS
C
      DO I=1,NBND
         THETA(I)=ARCTAN(RBND(I)-R0,(ZBND(I)-Z0)/E)
      END DO
C
      ERRL=1E10
      ERR3=0.0
      DO J=1,NOUT
         CALL LSFE(NBND,RBND,ZBND,THETA,SINT,COST,AN,BN,E,ERR)
         CALL NTHN(NBND,RBND,ZBND,THETA,SINT,COST,AN,BN,E,ERR3)
CCC      WRITE(6,*) 'ITERATION ',J,'  ERR ',ERR,' ERR3 ',ERR3
      END DO
      OPAR(1) = AN(1)
      OPAR(2) = BN(1)*E
      OPAR(3) = AN(2)
      OPAR(4) = E
      OPAR(5) = AN(3)
      OPAR(6) = BN(2)
      OPAR(7) = BN(3)
      OPAR(8) = ERR3
C
CCC   WRITE(8,*) 'E = ',E
CCC   DO I=1,3
CCC      WRITE(8,*) 'AN',I,' = ',AN(I)
CCC   END DO
CCC   DO I=1,3
CCC      WRITE(8,*) 'BN',I,' = ',BN(I)
CCC   END DO
C
    9 CONTINUE
      RETURN
      END
*-----------------------------------------------------------------------
      FUNCTION ARCTAN(X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA TWOPI/6.28318 53071 79586 5D0/
      DATA THPIO2/4.71238 89803 84689 9D0/
      DATA PI/3.14159 26535 89793 2D0/
      DATA PIO2/1.57079 63267 94896 6D0/
      IF(DABS(X).LE.DABS(Y*1.0D-20)) GOTO 20
      T=DATAN(Y/X)
      T=DABS(T)
      IF(X.LT.0.0D0) T=PI-T
      IF(Y.LT.0.0D0) T=TWOPI-T
      GOTO 30
   20 CONTINUE
      T=0.0D0
      IF(Y) 10,30,40
   40 CONTINUE
      T=PIO2
      GOTO 30
   10 CONTINUE
      T=THPIO2
   30 CONTINUE
      ARCTAN=T
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE LSFE(NBND,RBND,ZBND,THETA,SINT,COST,AN,BN,E,ERR)
C
C  USE LEAST SQUARES FIT TO UPDATE AN(*),BN(*),E
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RBND(*),ZBND(*),THETA(*),SINT(*),COST(*),AN(*),BN(*)
      DIMENSION ABE(7)
      DIMENSION XMAT(49),COS2T(1000),SIN2T(1000), 
     ;          AMAT(7,7), FAC(7,7), IPVT(7)
      DATA TWOPI/6.28318 53071 79586 5E0/
	integer i,j,k
C
      Z=0.0
      ZZ=0.0
      ZS1=0.0
      ZS2=0.0
      ZC1=0.0
      ZC2=0.0
      R=0.0
      RS1=0.0
      RS2=0.0
      RC1=0.0
      RC2=0.0
      S1=0.0
      S2=0.0
      S3=0.0
      C1=0.0
      C2=0.0
      C3=0.0
      DO I=1,NBND
         S=SIN(THETA(I))
         C=COS(THETA(I))
         SINT(I)=S
         COST(I)=C
         S2T=2.*S*C
         C2T=2.*C*C-1.0
         SIN2T(I)=S2T
         COS2T(I)=C2T
         S3T=S2T*C+C2T*S
         C3T=C*C2T-S*S2T
C
         S1=S1+S
         S2=S2+S2T
         S3=S3+S3T
         C1=C1+C
         C2=C2+C2T
         C3=C3+C3T
C
         ZI=ZBND(I)
         Z=Z+ZI
         ZZ=ZZ+ZI*ZI
         ZC1=ZC1+ZI*C
         ZC2=ZC2+ZI*C2T
         ZS1=ZS1+ZI*S
         ZS2=ZS2+ZI*S2T
C
         RI=RBND(I)
         R=R+RI
         RC1=RC1+RI*C
         RC2=RC2+RI*C2T
         RS1=RS1+RI*S
         RS2=RS2+RI*S2T
      END DO
C BUILD LSF MATRIX AND RHS
C
      XNBND=NBND
C>>>  DO K=1,NBND
C  E PARAM
         ABE(1)=0.0
         XMAT(1)=ZZ
         XMAT(2)=0.0
         XMAT(3)=-ZS1
         XMAT(4)=ZS2
         XMAT(5)=-Z
         XMAT(6)=-ZC1
         XMAT(7)=-ZC2
C
C  AN(1) PARAM: R0
         ABE(2)=R
         XMAT(8)=XNBND
         XMAT(9)=C1
         XMAT(10)=C2
         XMAT(11)=0.0
         XMAT(12)=S1
         XMAT(13)=S2
C
C  AN(2) PARAM: CR0
         ABE(3)=RC1
         XMAT(14)=XNBND
         XMAT(15)=C3
         XMAT(16)=S1
         XMAT(17)=S2
         XMAT(18)=S3
C
C  AN(3)
         ABE(4)=RC2
         XMAT(19)=XNBND
         XMAT(20)=-S2
         XMAT(21)=-S1
         XMAT(22)=0.0
C
C  BN(1) PARAM: Z0
         ABE(5)=0.0
         XMAT(23)=XNBND
         XMAT(24)=C1
         XMAT(25)=C2
C
C  BN(2) PARAM
         ABE(6)=RS1
         XMAT(26)=XNBND
         XMAT(27)=C1
C
C  BN(3) PARAM
         ABE(7)=RS2
         XMAT(28)=XNBND
C
C>>>  END DO
C
      NV=7


cERN      CALL DPPF(XMAT,NV,1)
cERN      CALL DPPS(XMAT,NV,ABE,1)

	k = 0
       do j=1,NV
          do i=j,NV
             k = k+1
             AMAT(i,j)=XMAT(k)
		   AMAT(j,i)=XMAT(k)
         end do
      end do


	call DLFTSF(NV, AMAT, NV, FAC, NV, IPVT)
	call DLFSSF(NV, FAC, NV, IPVT, ABE, ABE)

C
C  UNPACK RESULTS
C
      E=1./ABE(1)
      DO I=1,3
         AN(I)=ABE(I+1)
      END DO
      DO I=1,3
         BN(I)=ABE(I+4)
      END DO
C
C  CALC RMS ERROR
C
      ACC=0.0
      DO J=1,NBND
         RACC=AN(1)+AN(2)*COST(J)+AN(3)*COS2T(J)+
     1        BN(2)*SINT(J)+BN(3)*SIN2T(J)
         ZACC=BN(1)+AN(2)*SINT(J)-AN(3)*SIN2T(J)+
     1        BN(2)*COST(J)+BN(3)*COS2T(J)
         ACCI=(RACC-RBND(J))**2+(ZACC*E-ZBND(J))**2
         ACC=ACC+ACCI
      END DO
      ERR=SQRT(ACC/NBND)
C
    9 CONTINUE
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE NTHN(NBND,RBND,ZBND,THETA,SINT,COST,AN,BN,E,ERR)
C
C  UPDATE THETA COORDS OF BOUNDARY POINTS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RBND(*),ZBND(*),THETA(*),SINT(*),COST(*),AN(*),BN(*)
      DATA TWOPI/6.28318 53071 79586 5E0/
C
      RE=1./E
      A1=AN(1)
      A2=AN(2)
      A3=AN(3)
      B1=BN(1)
      B2=BN(2)
      B3=BN(3)
      ANP3=A3*2.
      ANPP3=A3*4.
      BNP3=B3*2.
      BNPP3=B3*4.
C
      ACC=0.0
      DO I=1,NBND
         T=THETA(I)
         S1=SINT(I)
         C1=COST(I)
         DO 10 NRIT=1,3
            S2=2.*S1*C1
            C2=2.*C1**2-1.0
            R=A1+A2*C1+A3*C2+BN(2)*S1+B3*S2
            DRDT=-A2*S1-ANP3*S2+BN(2)*C1+BNP3*C2
            D2RDT2=-A2*C1-ANPP3*C2-BN(2)*S1-BNPP3*S2
            Z=B1+A2*S1-A3*S2+BN(2)*C1+B3*C2
            DZDT=A2*C1-ANP3*C2-BN(2)*S1-BNP3*S2
            D2ZDT2=-A2*S1+ANPP3*S2-BN(2)*C1-BNPP3*C2
            DR=R-RBND(I)
            DZ=Z-ZBND(I)*RE
            F=DR*DRDT+DZ*DZDT
            FP=D2RDT2*DR+DRDT**2+D2ZDT2*DZ+DZDT**2
            DT=F/FP
            T=T-DT
            S1N=S1-DT*C1
            C1=C1+DT*S1
            S1=S1N
   10    CONTINUE
CCC      WRITE(8,*) I,THETA(I),T,R-RBND(I),Z*E-ZBND(I)
         THETA(I)=T
         ACC=ACC+((R-RBND(I))**2+(Z*E-ZBND(I))**2)
      END DO
      ERR=SQRT(ACC/NBND)
C
    9 CONTINUE
      RETURN
      END

