      SUBROUTINE PA03A(AA,R,N)

C######CALLS FD05A

      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION AA(4), R(3)

C**************************************************************
C*  PURPOSE....                                               *
C*  TO FIND THE ROOTS OF THE REAL CUBIC ......                *
C*          A(4)*X**3+A(3)*X**2+A(2)*X+A(1)                   *
C*                                                            *
C*  ARGUMENT LIST....                                         *
C*  AA  IS A REAL ARRAY OF LENGTH 4 WHICH WILL CONTAIN THE    *
C*      COEFFICIENTS A.                                       *
C*  R   IS A REAL ARRAY WHICH WILL HAVE ITS COMPONENTS        *
C*      SET TO THE ROOTS.IF THERE ARE THREE REAL ROOTS,THEN   *
C*      R(1).LE.R(2).LE.R(3). FOR ONE REAL ROOT,R(1) IS SET   *
C*      TO IT,R(2) IS SET TO THE REAL PART OF BOTH COMPLEX    *
C*      ROOTS AND R(3) IS SET TO THE IMAGINARY PART WHICH IS  *
C*      IS POSITIVE.THE DUMMY VALUE OF 1D+70 IS RETURNED FOR  *
C*      EACH INFINITE ROOT.THIS CORRESPONDS TO A ZERO         *
C*      LEADING COEFFICIENT.                                  *
C*  N   IS AN INTEGER WHOSE VALUE WILL BE SET BY THE ROUTINE  *
C*      TO THE NUMBER OF REAL ROOTS.                          *
C**************************************************************
      DOUBLE PRECISION A(3),B(2),C1,C5,C3,C0,C4,C23, XINF, X, Y, Z, TA, TB
     ;, TC, TE, TF
      REAL*8 FD05A
      EXTERNAL FD05A 
      DATA   C1/1./,C5/.5/,C3/3./,C0/0./,C4/4./
      DATA   C23/0.666666666666666667/
C     DATA   XINF/1.D+70/
C
C     IS LEADING COEFFICIENT ZERO,IE INFINITE ROOT.
      XINF=FD05A(5)*.1
      IF(AA(4).EQ.C0)GO TO 80
      IF(AA(1).EQ.C0)GO TO 70
      A(3)=AA(3)/(C3*AA(4))
      A(2)=AA(2)/(C3*AA(4))
      A(1)=AA(1)/AA(4)
      X=A(2)-A(3)*A(3)
      Y=A(1)-A(3)*(X+X+A(2))
      Z=Y**2+C4*X**3
      IF(Z.GE.C0)GOTO100
C
C     THERE ARE THREE REAL ROOTS.
      N=3
      R(1)=-2.0*SQRT(-X)
      Y=Y/(R(1)*X)
      X=R(1)
      Y =ATAN2(SQRT(C1-Y),SQRT(C1+Y))*C23
      IF(A(3).LT.C0)Y=Y+2.094395102393195
C
C     CALCULATE ROOT WHICH DOES NOT INVOLVE CANCELLATION
      R(1)=X*COS(Y)-A(3)
C
C     DEFLATE CUBIC FROM OPTIMAL END.
   10 B(1)=-A(1)/R(1)
      B(2)=(B(1)-C3*A(2))/R(1)
      IF(ABS(R(1)**3).LE.ABS(A(1)))B(2)=R(1)+C3*A(3)
   20 X=B(2)*B(2)-C4*B(1)
C
C     IS THE PAIR OF ROOTS REAL OR COMPLEX.
      IF(X.LT.C0)GOTO 60
      R(3)=-SIGN(C5,B(2))*(SQRT(X)+ABS(B(2)))
      R(2)=C0
      IF(R(3).NE.C0)R(2)=B(1)/R(3)
      IF(R(1).LE.R(2))GOTO30
      TA=R(2)
      R(2)=R(1)
      R(1)=TA
   30 IF(R(2).LE.R(3))GO TO 50
      TA=R(3)
      IF(R(1).LE.R(3))GO TO 40
      TA=R(1)
      R(1)=R(3)
   40 R(3)=R(2)
      R(2)=TA
   50 N=3
      GO TO 150
   60 R(2)=-C5*B(2)
      R(3)=C5*SQRT(-X)
      N=1
      GO TO 150
   70 R(1)=C0
      B(1)=AA(2)/AA(4)
      B(2)=AA(3)/AA(4)
      GO TO 20
C
C     THE CUBIC HAS LEADING COEFFICIENT ZERO,IE QUADRATIC
   80 R(1)=XINF
      IF(AA(3).EQ.C0)GO TO 90
      B(1)=AA(1)/AA(3)
      B(2)=AA(2)/AA(3)
      GO TO 20
C
C     CUBIC HAS FIRST TWO LEADING COEFFICIENTS ZERO
   90 IF(AA(2).NE.C0)R(1)=-AA(1)/AA(2)
      R(2)=XINF
      R(3)=XINF
      N=3
      GO TO 150
C
C     THERE IS ONE REAL ROOT.
  100 N=1
      TA=SQRT(Z)
      TB=(ABS(Y)+TA)*C5
      TC=TB**0.333333333333333333
      IF(TC.GT.0.0)GO TO 110
      R(1)=-A(3)
      R(2)=-A(3)
      R(3)=-A(3)
      N=3
      GO TO 150
  110 TC=TC-(TC**3-TB)/(C3*TC*TC)
      TE=TC*TC+ABS(X)
      TF=C1/((X/TC)**2+TE)
      IF(X.LT.C0)GOTO120
      X=TE/TC
      Z=Y*TF
      GOTO130
  120 X=TA*TF
      Z=SIGN(C1,Y)*TE/TC
  130 IF(Z*A(3).LT.C0)GOTO140
      R(1)=-Z-A(3)
      GO TO 10
  140 R(2)=C5*Z-A(3)
      R(3)=C5*SQRT(C3)*ABS(X)
      R(1)=-A(1)/(R(2)*R(2)+R(3)*R(3))
  150 RETURN
      END
