C-----------------------------------------------------------------------
C  IMSL Name:  IZAMAX (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Find the smallest index of the component of a
C              double-complex vector having maximum magnitude.
C
C  Usage:      IZAMAX(N, ZX, INCX)
C
C  Arguments:
C     N      - Length of vector X.  (Input)
C     ZX     - Double complex vector of length N*INCX.  (Input)
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be ZX(1+(I-1)*INCX). INCX must be
C              greater than zero.
C     IZAMAX - The smallest index I such that DCABS(X(I)) is the maximum
C              of DCABS(X(J)) for J=1 to N.  (Output)
C              X(I) refers to a specific element of ZX.
C
C  GAMS:       D1a2
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION IZAMAX (N, ZX, INCX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX
      COMPLEX    *16 ZX(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX
      DOUBLE PRECISION SMAX
      COMPLEX    *16 ZDUM
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DABS,DBLE,DIMAG
      INTRINSIC  DABS, DBLE, DIMAG
      DOUBLE PRECISION DABS, DBLE, DIMAG
      DOUBLE PRECISION ZABS1
C
      ZABS1(ZDUM) = DABS(DBLE(ZDUM)) + DABS(DIMAG(ZDUM))
C
      IZAMAX = 0
      IF (N .GE. 1) THEN
         IZAMAX = 1
         IF (N .NE. 1) THEN
            IF (INCX .NE. 1) THEN
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1
               IX = 1
               IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
               SMAX = ZABS1(ZX(IX))
               IX = IX + INCX
               DO 10  I=2, N
                  IF (ZABS1(ZX(IX)) .GT. SMAX) THEN
                     IZAMAX = I
                     SMAX = ZABS1(ZX(IX))
                  END IF
                  IX = IX + INCX
   10          CONTINUE
            ELSE
C                                  CODE FOR INCREMENTS EQUAL TO 1
               SMAX = ZABS1(ZX(1))
               DO 20  I=2, N
                  IF (ZABS1(ZX(I)) .GT. SMAX) THEN
                     IZAMAX = I
                     SMAX = ZABS1(ZX(I))
                  END IF
   20          CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  N1RTY
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 6, 1984
C
C  Purpose:    Retrieve an error type.
C
C  Usage:      N1RTY(IOPT)
C
C  Arguments:
C     IOPT   - Integer specifying the level.  (Input)
C              If IOPT=0 the error type for the current level is
C              returned.  If IOPT=1 the error type for the most
C              recently called routine (last pop) is returned.
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION N1RTY (IOPT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IOPT
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                              SPECIFICATIONS FOR COMMON /ERCOM1/
      INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51),
     &           PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7,
     &           IALLOC(51), HDRFMT(7), TRACON(7)
      COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE,
     &           PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT,
     &           TRACON
      SAVE       /ERCOM1/
C                              SPECIFICATIONS FOR COMMON /ERCOM2/
      CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
      COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
      SAVE       /ERCOM2/
C                              SPECIFICATIONS FOR COMMON /ERCOM3/
      DOUBLE PRECISION ERCKSM
      COMMON     /ERCOM3/ ERCKSM
      SAVE       /ERCOM3/
C                              SPECIFICATIONS FOR COMMON /ERCOM4/
      LOGICAL    ISUSER(51)
      COMMON     /ERCOM4/ ISUSER
      SAVE       /ERCOM4/
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1PRT, M1VECH
C
      IF (IOPT.NE.0 .AND. IOPT.NE.1) THEN
         ERTYPE(CALLVL) = 5
         ERCODE(CALLVL) = 1
         MSGLEN = 47
         CALL M1VECH ('.  The argument passed to N1RTY must be 0 or '//
     &                '1. ', MSGLEN, MSGSAV, MSGLEN)
         CALL E1PRT
         STOP
      ELSE
         N1RTY = ERTYPE(CALLVL+IOPT)
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1STD
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 6, 1984
C
C  Purpose:    To store a real number for subsequent use within an error
C              message.
C
C  Usage:      CALL E1STD(ID, DVALUE)
C
C  Arguments:
C     ID     - Integer specifying the substitution index.  ID must be
C              between 1 and 9.  (Input)
C     DVALUE - The double precision number to be stored.  (Input)
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1STD (ID, DVALUE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    ID
      DOUBLE PRECISION DVALUE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IBEG, ILEN
      CHARACTER  ARRAY(24), SAVE*24
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    IFINIT
      CHARACTER  BLANK(1)
      SAVE       BLANK, IFINIT
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                              SPECIFICATIONS FOR COMMON /ERCOM1/
      INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51),
     &           PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7,
     &           IALLOC(51), HDRFMT(7), TRACON(7)
      COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE,
     &           PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT,
     &           TRACON
      SAVE       /ERCOM1/
C                              SPECIFICATIONS FOR COMMON /ERCOM2/
      CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
      COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
      SAVE       /ERCOM2/
C                              SPECIFICATIONS FOR COMMON /ERCOM3/
      DOUBLE PRECISION ERCKSM
      COMMON     /ERCOM3/ ERCKSM
      SAVE       /ERCOM3/
C                              SPECIFICATIONS FOR COMMON /ERCOM4/
      LOGICAL    ISUSER(51)
      COMMON     /ERCOM4/ ISUSER
      SAVE       /ERCOM4/
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1INIT, E1INPL
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1ERIF
      INTEGER    I1ERIF
C
      DATA BLANK/' '/, IFINIT/0/
C                                  INITIALIZE IF NECESSARY
      IF (IFINIT .EQ. 0) THEN
         CALL E1INIT
         IFINIT = 1
      END IF
      IF (DVALUE .EQ. 0.0D0) THEN
         WRITE (SAVE,'(D24.15)') DVALUE
      ELSE
         WRITE (SAVE,'(1PD24.15)') DVALUE
      END IF
      READ (SAVE,'(24A1)') (ARRAY(I),I=1,24)
      IBEG = I1ERIF(ARRAY,24,BLANK,1)
      IF (ID.GE.1 .AND. ID.LE.9) THEN
         ILEN = 25 - IBEG
         CALL E1INPL ('D', ID, ILEN, ARRAY(IBEG))
      END IF
C
      RETURN
      END

C-----------------------------------------------------------------------
C  IMSL Name:  CCGCG/DCCGCG (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    June 5, 1985
C
C  Purpose:    Copy a complex general matrix.
C
C  Usage:      CALL CCGCG (N, A, LDA, B, LDB)
C
C  Arguments:
C     N      - Order of the matrices A and B.  (Input)
C     A      - Complex matrix of order N.  (Input)
C     LDA    - Leading dimension of A exactly as specified in the
C              dimension statement of the calling program.  (Input)
C     B      - Complex matrix of order N containing a copy of A.
C              (Output)
C     LDB    - Leading dimension of B exactly as specified in the
C              dimension statement of the calling program.  (Input)
C
C  GAMS:       D1b8
C
C  Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DCCGCG (N, A, LDA, B, LDB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, LDA, LDB
      COMPLEX    *16 A(LDA,*), B(LDB,*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    J
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   ZCOPY, E1MES, E1POP, E1PSH, E1STI
C
      CALL E1PSH ('DCCGCG ')
C                                  Check N
      IF (N .LT. 1) THEN
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'The argument N = %(I1).  It must be at '//
     &               'least 1.')
         GO TO 9000
      END IF
C                                  Check LDA
      IF (LDA .LT. N) THEN
         CALL E1STI (1, LDA)
         CALL E1STI (2, N)
         CALL E1MES (5, 2, 'The argument LDA = %(I1).  It must be at '//
     &               'least as large as N = %(I2).')
         GO TO 9000
      END IF
C                                  Check LDB
      IF (LDB .LT. N) THEN
         CALL E1STI (1, LDB)
         CALL E1STI (2, N)
         CALL E1MES (5, 3, 'The argument LDB = %(I1).  It must be at '//
     &               'least as large as N = %(I2).')
         GO TO 9000
      END IF
C                                  Copy
      IF (LDA.EQ.N .AND. LDB.EQ.N) THEN
         CALL ZCOPY (N*N, A, 1, B, 1)
      ELSE IF (LDA .GE. LDB) THEN
         DO 10  J=1, N
            CALL ZCOPY (N, A(1,J), 1, B(1,J), 1)
   10    CONTINUE
      ELSE
         DO 20  J=N, 1, -1
            CALL ZCOPY (N, A(1,J), -1, B(1,J), -1)
   20    CONTINUE
      END IF
C
 9000 CONTINUE
      CALL E1POP ('DCCGCG ')
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  DZASUM (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Sum the absolute values of the real part together with
C              the absolute values of the imaginary part of the
C              components of a double-complex vector.
C
C  Usage:      DZASUM(N, ZX, INCX)
C
C  Arguments:
C     N      - Length of vectors X.  (Input)
C     ZX     - Double complex vector of length N*INCX.  (Input)
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be ZX(1+(I-1)*INCX).  INCX must be
C              greater than 0.
C     DZASUM - Sum from I=1 to N of DABS(REAL(X(I)))*DABS(AIMAG(X(I)))).
C              (Output)
C              X(I) refers to a specific element of ZX.
C
C  GAMS:       D1a
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DZASUM (N, ZX, INCX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX
      COMPLEX    *16 ZX(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX
      DOUBLE PRECISION STEMP
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DABS,DBLE,DIMAG
      INTRINSIC  DABS, DBLE, DIMAG
      DOUBLE PRECISION DABS, DBLE, DIMAG
C
      DZASUM = 0.0D0
      STEMP = 0.0D0
      IF (N .GT. 0) THEN
         IF (INCX .NE. 1) THEN
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1
            IX = 1
            IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
            DO 10  I=1, N
               STEMP = STEMP + DABS(DIMAG(ZX(IX))) + DABS(DBLE(ZX(IX)))
               IX = IX + INCX
   10       CONTINUE
            DZASUM = STEMP
         ELSE
C                                  CODE FOR INCREMENTS EQUAL TO 1
            DO 20  I=1, N
               STEMP = STEMP + DABS(DIMAG(ZX(I))) + DABS(DBLE(ZX(I)))
   20       CONTINUE
            DZASUM = STEMP
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ZDSCAL (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Multiply a double-complex vector by a double-precision
C              scalar, y = ay.
C
C  Usage:      CALL ZDSCAL (N, DA, ZX, INCX)
C
C  Arguments:
C     N      - Length of vectors X.  (Input)
C     DA     - Double precison scalar.  (Input)
C     ZX     - Double complex vector of length MAX(N*IABS(INCX),1).
C                 (Input/Output)
C              ZDSCAL replaces X(I) with DA*X(I) for I = 1,...,N.
C              X(I) refers to a specific element of ZX.
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be ZX(1+(I-1)*INCX). INCX must be
C              greater than 0.
C
C  GAMS:       D1a6
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZDSCAL (N, DA, ZX, INCX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX
      DOUBLE PRECISION DA
      COMPLEX    *16 ZX(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DCMPLX
      INTRINSIC  DCMPLX
      COMPLEX    *16 DCMPLX
C
      IF (N .GT. 0) THEN
         IF (INCX .NE. 1) THEN
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1
            IX = 1
            IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
            DO 10  I=1, N
               ZX(IX) = DCMPLX(DA,0.0D0)*ZX(IX)
               IX = IX + INCX
   10       CONTINUE
         ELSE
C                                  CODE FOR INCREMENTS EQUAL TO 1
            DO 20  I=1, N
               ZX(I) = DCMPLX(DA,0.0D0)*ZX(I)
   20       CONTINUE
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ZGERU (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    July 24, 1986
C
C  Purpose:    Perform the rank-one matrix update A = A + alpha*x*y',
C              where A is a full storage mode matrix, all double
C              complex.
C
C  Usage:      CALL ZGERU (M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
C
C  Arguments:
C     M      - Number of rows in A.  (Input)
C     N      - Number of columns in A.  (Input)
C     ALPHA  - Double complex scalar.  (Input)
C     X      - Double complex vector of length (M-1)*IABS(INCX)+1.
C              (Input)
C     INCX   - Displacement between elements of X.  (Input)
C     Y      - Double complex vector of length (N-1)*IABS(INCY)+1.
C              (Input)
C     INCY   - Displacement between elements of Y.  (Input)
C     A      - Double complex array of size M by N.  (Input/Output)
C              On input, A contains the matrix to be updated.
C              On output, A contains the updated matrix.
C     LDA    - Leading dimension of A exactly as specified in the
C              calling routine.  (Input)
C
C  GAMS:       D1b
C
C  Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZGERU (M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    M, N, INCX, INCY, LDA
      COMPLEX    *16 ALPHA, X(*), Y(*)
      COMPLEX    *16 A(*)
      INTEGER    I1X
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IY, J
C                                  SPECIFICATIONS FOR SPECIAL CASES
      EXTERNAL   ZAXPY
C                                  Quick return if possible
      IF (M.EQ.0 .OR. N.EQ.0 .OR. ALPHA.EQ.(0.0D0,0.0D0)) GO TO 9000
C
      IY = 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
C
      I1X = 1
      DO 10  J=1, N
         CALL ZAXPY (M, ALPHA*Y(IY), X, INCX, A(I1X), 1)
         IY = IY + INCY
         I1X = I1X + LDA
   10 CONTINUE
C
 9000 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ZSCAL (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Multiply a vector by a scalar, y = ay, both double
C              complex.
C
C  Usage:      CALL ZSCAL (N, ZA, ZX, INCX)
C
C  Arguments:
C     N      - Length of vectors X.  (Input)
C     ZA     - Double complex scalar.  (Input)
C     ZX     - Double complex vector of length MAX(N*IABS(INCX),1).
C                 (Input/Output)
C              ZSCAL replaces X(I) with ZA*X(I) for I = 1,...,N.
C              X(I) refers to a specific element of ZX.
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be ZX(1+(I-1)*INCX). INCX must be
C              greater than 0.
C
C  GAMS:       D1a6
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSCAL (N, ZA, ZX, INCX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX
      COMPLEX    *16 ZA, ZX(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX
C
      IF (N .GT. 0) THEN
         IF (INCX .NE. 1) THEN
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1
            IX = 1
            IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
            DO 10  I=1, N
               ZX(IX) = ZA*ZX(IX)
               IX = IX + INCX
   10       CONTINUE
         ELSE
C                                  CODE FOR INCREMENTS EQUAL TO 1
            DO 20  I=1, N
               ZX(I) = ZA*ZX(I)
   20       CONTINUE
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ZSWAP (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Interchange vectors X and Y, both double complex.
C
C  Usage:      CALL ZSWAP (N, ZX, INCX, ZY, INCY)
C
C  Arguments:
C     N      - Length of vectors X and Y.  (Input)
C     ZX     - Double complex vector of length MAX(N*IABS(INCX),1).
C              (Input/Output)
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be
C                ZX(1+(I-1)*INCX) if INCX.GE.0  or
C                ZX(1+(I-N)*INCX) if INCX.LT.0.
C     ZY     - Double complex vector of length MAX(N*IABS(INCY),1).
C              (Input/Output)
C     INCY   - Displacement between elements of ZY.  (Input)
C              Y(I) is defined to be
C                ZY(1+(I-1)*INCY) if INCY.GE.0  or
C                ZY(1+(I-N)*INCY) if INCY.LT.0.
C
C  GAMS:       D1a5
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSWAP (N, ZX, INCX, ZY, INCY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX, INCY
      COMPLEX    *16 ZX(*), ZY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX, IY
      COMPLEX    *16 ZTEMP
C
      IF (N .GT. 0) THEN
         IF (INCX.NE.1 .OR. INCY.NE.1) THEN
C                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL
C                                    INCREMENTS NOT EQUAL TO 1
            IX = 1
            IY = 1
            IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
            IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
            DO 10  I=1, N
               ZTEMP = ZX(IX)
               ZX(IX) = ZY(IY)
               ZY(IY) = ZTEMP
               IX = IX + INCX
               IY = IY + INCY
   10       CONTINUE
         ELSE
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
            DO 20  I=1, N
               ZTEMP = ZX(I)
               ZX(I) = ZY(I)
               ZY(I) = ZTEMP
   20       CONTINUE
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ZAXPY (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Compute the scalar times a vector plus a vector,
C              y = ax + y, all double complex.
C
C  Usage:      CALL ZAXPY (N, ZA, ZX, INCX, ZY, INCY)
C
C  Arguments:
C     N      - Length of vectors X and Y.  (Input)
C     ZA     - Complex scalar.  (Input)
C     ZX     - Complex vector of length MAX(N*IABS(INCX),1).  (Input)
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be
C                 ZX(1+(I-1)*INCX) if INCX.GE.0  or
C                 ZX(1+(I-N)*INCX) if INCX.LT.0.
C     ZY     - Complex vector of length MAX(N*IABS(INCY),1).
C                 (Input/Output)
C              ZAXPY replaces Y(I) with ZA*X(I) + Y(I) for I = 1,...N.
C              X(I) and Y(I) refer to specific elements of ZX and ZY.
C     INCY   - Displacement between elements of ZY.  (Input)
C              Y(I) is defined to be
C                 ZY(1+(I-1)*INCY) if INCY.GE.0  or
C                 ZY(1+(I-N)*INCY) if INCY.LT.0.
C
C  GAMS:       D1a7
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZAXPY (N, ZA, ZX, INCX, ZY, INCY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX, INCY
      COMPLEX    *16 ZA, ZX(*), ZY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX, IY
      DOUBLE PRECISION ZZABS
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DABS,DBLE,DIMAG
      INTRINSIC  DABS, DBLE, DIMAG
      DOUBLE PRECISION DABS, DBLE, DIMAG
C
      IF (N .GT. 0) THEN
         ZZABS = DABS(DBLE(ZA)) + DABS(DIMAG(ZA))
         IF (ZZABS .NE. 0.0D0) THEN
            IF (INCX.NE.1 .OR. INCY.NE.1) THEN
C                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL
C                                    INCREMENTS NOT EQUAL TO 1
               IX = 1
               IY = 1
               IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
               IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
               DO 10  I=1, N
                  ZY(IY) = ZY(IY) + ZA*ZX(IX)
                  IX = IX + INCX
                  IY = IY + INCY
   10          CONTINUE
            ELSE
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
               DO 20  I=1, N
                  ZY(I) = ZY(I) + ZA*ZX(I)
   20          CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ZCOPY (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Copy a vector X to a vector Y, both double complex.
C
C  Usage:      CALL ZCOPY (N, ZX, INCX, ZY, INCY)
C
C  Arguments:
C     N      - Length of vectors X and Y.  (Input)
C     ZX     - Double complex vector of length MAX(N*IABS(INCX),1).
C              (Input)
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be
C                 ZX(1+(I-1)*INCX) if INCX.GE.0  or
C                 ZX(1+(I-N)*INCX) if INCX.LT.0.
C     ZY     - Double complex vector of length MAX(N*IABS(INCY),1).
C              (Output)
C              ZCOPY copies X(I) to Y(I) for I = 1,...N. X(I) and Y(I)
C              refer to specific elements of ZX and ZY.
C     INCY   - Displacement between elements of ZY.  (Input)
C              Y(I) is defined to be
C                 ZY(1+(I-1)*INCY) if INCY.GE.0  or
C                 ZY(1+(I-N)*INCY) if INCY.LT.0.
C
C  GAMS:       D1a
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCOPY (N, ZX, INCX, ZY, INCY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX, INCY
      COMPLEX    *16 ZX(*), ZY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX, IY
C
      IF (N .GT. 0) THEN
         IF (INCX.NE.1 .OR. INCY.NE.1) THEN
C                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL
C                                    INCREMENTS NOT EQUAL TO 1
            IX = 1
            IY = 1
            IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
            IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
            DO 10  I=1, N
               ZY(IY) = ZX(IX)
               IX = IX + INCX
               IY = IY + INCY
   10       CONTINUE
         ELSE
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
            DO 20  I=1, N
               ZY(I) = ZX(I)
   20       CONTINUE
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ZSET (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Set the components of a vector to a scalar, all double
C              complex.
C
C  Usage:      CALL ZSET (N, ZA, ZX, INCX)
C
C  Arguments:
C     N      - Length of vector X.  (Input)
C     ZA     - Double complex scalar.  (Input)
C     ZX     - Double complex vector of length N*INCX.  (Input/Output)
C              ZSET replaces X(I) with ZA for I=1,...,N. X(I) refers to
C              a specific element of ZX. See INCX argument description.
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be ZX(1+(I-1)*INCX). INCX must be
C              greater than zero.
C
C  GAMS:       D1a1
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSET (N, ZA, ZX, INCX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX
      COMPLEX    *16 ZA, ZX(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, NINCX
C
      IF (N .GT. 0) THEN
         IF (INCX .NE. 1) THEN
C                                  CODE FOR INCREMENT NOT EQUAL TO 1
            NINCX = N*INCX
            DO 10  I=1, NINCX, INCX
               ZX(I) = ZA
   10       CONTINUE
         ELSE
C                                  CODE FOR INCREMENT EQUAL TO 1
            DO 20  I=1, N
               ZX(I) = ZA
   20       CONTINUE
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ZDOTC (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Compute the double-complex conjugate dot product,
C              conjg(x)*y.
C
C  Usage:      ZDOTC(N, ZX, INCX, ZY, INCY)
C
C  Arguments:
C     N      - Length of vectors X and Y.  (Input)
C     ZX     - Double complex vector of length MAX(N*IABS(INCX),1).
C              (Input)
C     INCX   - Displacement between elements of ZX.  (Input)
C              X(I) is defined to be
C                 ZX(1+(I-1)*INCX) if INCX.GE.0  or
C                 ZX(1+(I-N)*INCX) if INCX.LT.0.
C     ZY     - Double complex vector of length MAX(N*IABS(INCY),1).
C              (Input)
C     INCY   - Displacement between elements of ZY.  (Input)
C              Y(I) is defined to be
C                 ZY(1+(I-1)*INCY) if INCY.GE.0  or
C                 ZY(1+(I-N)*INCY) if INCY.LT.0.
C     ZDOTC  - Double complex sum from I=1 to N of CONJ(X(I))*Y(I).
C              (Output)
C              X(I) and Y(I) refer to specific elements of ZX and ZY
C              respectively.
C
C  GAMS:       D1a4
C
C  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      COMPLEX *16FUNCTION ZDOTC (N, ZX, INCX, ZY, INCY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX, INCY
      COMPLEX    *16 ZX(*), ZY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX, IY
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DCONJG
      INTRINSIC  DCONJG
      COMPLEX    *16 DCONJG
C
      ZDOTC = (0.0D0,0.0D0)
      IF (N .GT. 0) THEN
         IF (INCX.NE.1 .OR. INCY.NE.1) THEN
C                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL
C                                    INCREMENTS NOT EQUAL TO 1
            IX = 1
            IY = 1
            IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
            IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
            DO 10  I=1, N
               ZDOTC = ZDOTC + DCONJG(ZX(IX))*ZY(IY)
               IX = IX + INCX
               IY = IY + INCY
   10       CONTINUE
         ELSE
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
            DO 20  I=1, N
               ZDOTC = ZDOTC + DCONJG(ZX(I))*ZY(I)
   20       CONTINUE
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  DMACH (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    March 15, 1984
C
C  Purpose:    Generate double precision machine constants.
C
C  Usage:      DMACH(N)
C
C  Arguments:
C     N      - Index of desired constant.  (Input)
C     DMACH  - Machine constant.  (Output)
C              DMACH(1) = B**(EMIN-1), the smallest positive magnitude.
C              DMACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C              DMACH(3) = B**(-T), the smallest relative spacing.
C              DMACH(4) = B**(1-T), the largest relative spacing.
C              DMACH(5) = LOG10(B), the log, base 10, of the radix.
C              DMACH(6) = not-a-number.
C              DMACH(7) = positive machine infinity.
C              DMACH(8) = negative machine infinity.
C
C  GAMS:       R1
C
C  Chapters:   MATH/LIBRARY Reference Material
C              STAT/LIBRARY Reference Material
C              SFUN/LIBRARY Reference Material
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DMACH (N)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      DOUBLE PRECISION RMACH(8)
      SAVE       RMACH
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI
C                                  DEFINE CONSTANTS
      DATA RMACH(1)/2.93941D-39/
      DATA RMACH(2)/1.70102D38/
      DATA RMACH(3)/1.38810D-17/
      DATA RMACH(4)/2.77620D-17/
      DATA RMACH(5)/.3010299956639811952137388947245D0/
      DATA RMACH(6)/1.70110D38/
      DATA RMACH(7)/1.70102D38/
      DATA RMACH(8)/-1.70102D38/
C
      IF (N.LT.1 .OR. N.GT.8) THEN
         CALL E1PSH ('DMACH ')
         DMACH = RMACH(6)
         CALL E1STI (1, N)
         CALL E1MES (5, 5, 'The argument must be between 1 '//
     &               'and 8 inclusive. N = %(I1)')
         CALL E1POP ('DMACH ')
      ELSE
         DMACH = RMACH(N)
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  L2CCG/DL2CCG  (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the LU factorization of a complex general matrix
C              and estimate its L1 condition number.
C
C  Usage:      CALL L2CCG (N, A, LDA, FAC, LDFAC, IPVT, RCOND, Z)
C
C  Arguments:  See LFCCG/DLFCCG.
C
C  Remarks:    See LFCCG/DLFCCG.
C
C  Chapter:    MATH/LIBRARY Linear Systems
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DL2CCG (N, A, LDA, FAC, LDFAC, IPVT, RCOND, Z)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, LDA, LDFAC, IPVT(*)
      DOUBLE PRECISION RCOND
      COMPLEX    *16 A(LDA,*), FAC(LDFAC,*), Z(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    J, K, KP1, L
      DOUBLE PRECISION ANORM, S, SM, YNORM
      COMPLEX    *16 EK, T, WK, WKM, ZDUM, ZDUM1, ZDUM2, ZERO
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DABS,DIMAG,DMAX1,DCMPLX,DCONJG,DBLE
      INTRINSIC  DABS, DIMAG, DMAX1, DCMPLX, DCONJG, DBLE
      DOUBLE PRECISION DABS, DIMAG, DMAX1, DBLE
      COMPLEX    *16 DCMPLX, DCONJG
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   ZAXPY, ZSET, ZDSCAL, E1MES, E1POP, E1PSH, E1STI,
     &           E1STD, DL2TCG
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   DMACH, ZDOTC, N1RTY, DZASUM
      INTEGER    N1RTY
      DOUBLE PRECISION DMACH, DZASUM
      COMPLEX    *16 ZDOTC
      DOUBLE PRECISION CABS1
      COMPLEX    *16 CSIGN1
C
      CABS1(ZDUM) = DABS(DBLE(ZDUM)) + DABS(DIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C
      CALL E1PSH ('DL2CCG ')
C
      IF (N .LE. 0) THEN
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'The order of the matrix must be '//
     &               'positive while N = %(I1) is given.')
         GO TO 9000
      END IF
C
      IF (N .GT. LDA) THEN
         CALL E1STI (1, N)
         CALL E1STI (2, LDA)
         CALL E1MES (5, 2, 'The order of the matrix must be '//
     &               'less than or equal to its leading dimension '//
     &               'while N = %(I1) and LDA = %(I2) are given.')
         GO TO 9000
      END IF
C
      IF (N .GT. LDFAC) THEN
         CALL E1STI (1, N)
         CALL E1STI (2, LDFAC)
         CALL E1MES (5, 3, 'The order of the matrix must be '//
     &               'less than or equal to its leading dimension '//
     &               'while N = %(I1) and LDFAC = %(I2) are given.')
         GO TO 9000
      END IF
C                                  COMPUTE 1-NORM OF A
      RCOND = 0.0D0
      ANORM = 0.0D0
      DO 10  J=1, N
         ANORM = DMAX1(ANORM,DZASUM(N,A(1,J),1))
   10 CONTINUE
C                                  FACTORIZATION STEP
C
      CALL DL2TCG (N, A, LDA, FAC, LDFAC, IPVT, Z)
      IF (N1RTY(1) .EQ. 4) GO TO 9000
C                                  RCOND = 1/(NORM(A)*(ESTIMATE OF
C                                  NORM(INVERSE(A)))).  ESTIMATE =
C                                  NORM(Z)/NORM(Y) WHERE A*Z = Y AND
C                                  CTRANS(A)*Y = E.  CTRANS(A) IS THE
C                                  TRANSPOSE OF A.  THE COMPONENTS OF E
C                                  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C                                  GROWTH IN THE ELEMENTS OF W WHERE
C                                  CTRANS(U)*W = E.  THE VECTORS ARE
C                                  FREQUENTLY RESCALED TO AVOID
C                                  OVERFLOW.
C                                  SOLVE CTRANS(U)*W = E
      ZERO = DCMPLX(0.0D0,0.0D0)
      EK = DCMPLX(1.0D0,0.0D0)
      CALL ZSET (N, ZERO, Z, 1)
      DO 40  K=1, N
         IF (CABS1(Z(K)) .NE. 0.0D0) EK = CSIGN1(EK,-Z(K))
         IF (CABS1(EK-Z(K)) .GT. CABS1(FAC(K,K))) THEN
            S = CABS1(FAC(K,K))/CABS1(EK-Z(K))
            CALL ZDSCAL (N, S, Z, 1)
            EK = DCMPLX(S,0.0D0)*EK
         END IF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(FAC(K,K)) .NE. 0.0D0) THEN
            WK = WK/DCONJG(FAC(K,K))
            WKM = WKM/DCONJG(FAC(K,K))
         ELSE
            WK = DCMPLX(1.0D0,0.0D0)
            WKM = DCMPLX(1.0D0,0.0D0)
         END IF
         KP1 = K + 1
         IF (KP1 .LE. N) THEN
            DO 20  J=KP1, N
               SM = SM + CABS1(Z(J)+WKM*DCONJG(FAC(K,J)))
               Z(J) = Z(J) + WK*DCONJG(FAC(K,J))
               S = S + CABS1(Z(J))
   20       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               DO 30  J=KP1, N
                  Z(J) = Z(J) + T*DCONJG(FAC(K,J))
   30          CONTINUE
            END IF
         END IF
         Z(K) = WK
   40 CONTINUE
      S = 1.0D0/DZASUM(N,Z,1)
      CALL ZDSCAL (N, S, Z, 1)
C                                  SOLVE CTRANS(L)*Y = W
      DO 50  K=N, 1, -1
         IF (K .LT. N) Z(K) = Z(K) + ZDOTC(N-K,FAC(K+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .GT. 1.0D0) THEN
            S = 1.0D0/CABS1(Z(K))
            CALL ZDSCAL (N, S, Z, 1)
         END IF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
   50 CONTINUE
      S = 1.0D0/DZASUM(N,Z,1)
      CALL ZDSCAL (N, S, Z, 1)
C
      YNORM = 1.0D0
C                                  SOLVE L*V = Y
      DO 60  K=1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL ZAXPY (N-K, T, FAC(K+1,K), 1, Z(K+1), 1)
         IF (CABS1(Z(K)) .GT. 1.0D0) THEN
            S = 1.0D0/CABS1(Z(K))
            CALL ZDSCAL (N, S, Z, 1)
            YNORM = S*YNORM
         END IF
   60 CONTINUE
      S = 1.0D0/DZASUM(N,Z,1)
      CALL ZDSCAL (N, S, Z, 1)
      YNORM = S*YNORM
C                                  SOLVE U*Z = V
      DO 70  K=N, 1, -1
         IF (CABS1(Z(K)) .GT. CABS1(FAC(K,K))) THEN
            S = CABS1(FAC(K,K))/CABS1(Z(K))
            CALL ZDSCAL (N, S, Z, 1)
            YNORM = S*YNORM
         END IF
         IF (CABS1(FAC(K,K)) .NE. 0.0D0) THEN
            Z(K) = Z(K)/FAC(K,K)
         ELSE
            Z(K) = DCMPLX(1.0D0,0.0D0)
         END IF
         T = -Z(K)
         CALL ZAXPY (K-1, T, FAC(1,K), 1, Z(1), 1)
   70 CONTINUE
C                                  MAKE ZNORM = 1.0
      S = 1.0D0/DZASUM(N,Z,1)
      CALL ZDSCAL (N, S, Z, 1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (RCOND .LE. DMACH(4)) THEN
         CALL E1STD (1, RCOND)
         CALL E1MES (3, 1, 'The matrix is algorithmically '//
     &               'singular.  An estimate of the reciprocal '//
     &               'of its L1 condition number is RCOND = %(D1).')
      END IF
C
 9000 CALL E1POP ('DL2CCG ')
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  L2TCG/DL2TCG (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the LU factorization of a complex general matrix.
C
C  Usage:      CALL L2TCG (N, A, LDA, FAC, LDFAC, IPVT, SCALE)
C
C  Arguments:  See LFTCG/DLFTCG.
C
C  Remarks:    See LFTCG/DLFTCG.
C
C  Chapter:    MATH/LIBRARY Linear Systems
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DL2TCG (N, A, LDA, FAC, LDFAC, IPVT, SCALE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, LDA, LDFAC, IPVT(*)
      COMPLEX    *16 A(LDA,*), FAC(LDFAC,*), SCALE(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, INDJ, INFO, J, K, L
      DOUBLE PRECISION BIG, CURMAX, SMALL, VALUE
      COMPLEX    *16 T, ZDUM
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DABS,DIMAG,DCMPLX,DBLE
      INTRINSIC  DABS, DIMAG, DCMPLX, DBLE
      DOUBLE PRECISION DABS, DIMAG, DBLE
      COMPLEX    *16 DCMPLX
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   ZCOPY, ZGERU, ZSCAL, ZSWAP, E1MES, E1POP, E1PSH, E1STI
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   DMACH, IZAMAX
      INTEGER    IZAMAX
      DOUBLE PRECISION DMACH
      DOUBLE PRECISION CABS1
C
      CABS1(ZDUM) = DABS(DBLE(ZDUM)) + DABS(DIMAG(ZDUM))
C
      CALL E1PSH ('DL2TCG ')
C
      IF (N .LE. 0) THEN
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'The order of the matrix must be '//
     &               'positive while N = %(I1) is given.')
         GO TO 9000
      END IF
C
      IF (N .GT. LDA) THEN
         CALL E1STI (1, N)
         CALL E1STI (2, LDA)
         CALL E1MES (5, 2, 'The order of the matrix must be '//
     &               'less than or equal to its leading dimension '//
     &               'while N = %(I1) and LDA = %(I2) are given.')
         GO TO 9000
      END IF
C
      IF (N .GT. LDFAC) THEN
         CALL E1STI (1, N)
         CALL E1STI (2, LDFAC)
         CALL E1MES (5, 3, 'The order of the matrix must be '//
     &               'less than or equal to its leading dimension '//
     &               'while N = %(I1) and LDFAC = %(I2) are given.')
         GO TO 9000
      END IF
C                                  PRESERVE A COPY OF THE INPUT MATRIX
      DO 10  J=1, N
         CALL ZCOPY (N, A(1,J), 1, FAC(1,J), 1)
   10 CONTINUE
C                                  COMPUTE THE INFINITY NORM OF EACH
C                                  ROW OF A FOR SCALING PURPOSE
      DO 20  I=1, N
         INDJ = IZAMAX(N,FAC(I,1),LDFAC)
         SCALE(I) = DCMPLX(CABS1(FAC(I,INDJ)),0.0D0)
   20 CONTINUE
C                                  GAUSSIAN ELIMINATION WITH PARTIAL
C                                  PIVOTING
      INFO = 0
      SMALL = DMACH(1)
      BIG = DMACH(2)
      IF (SMALL*BIG .LT. 1.0D0) SMALL = 1.0D0/BIG
      DO 40  K=1, N - 1
C                                  FIND L = PIVOT INDEX
         L = K
         CURMAX = 0.0D0
         DO 30  I=K, N
            IF (CABS1(SCALE(I)) .NE. 0.0D0) THEN
               VALUE = CABS1(FAC(I,K)/SCALE(I))
            ELSE
               VALUE = CABS1(FAC(I,K))
            END IF
            IF (VALUE .GT. CURMAX) THEN
               CURMAX = VALUE
               L = I
            END IF
   30    CONTINUE
         IPVT(K) = L
C                                  ZERO PIVOT IMPLIES THIS COLUMN
C                                  ALREADY TRIANGULARIZED
         IF (CABS1(FAC(L,K)) .NE. 0.0D0) THEN
C                                  INTERCHANGE IF NECESSARY
            IF (L .NE. K) THEN
               T = FAC(L,K)
               FAC(L,K) = FAC(K,K)
               FAC(K,K) = T
               SCALE(L) = SCALE(K)
            END IF
C                                  COMPUTE MULTIPLIERS
            IF (CABS1(FAC(K,K)) .GT. SMALL) THEN
               T = DCMPLX(-1.0D0,0.0D0)/FAC(K,K)
               CALL ZSCAL (N-K, T, FAC(K+1,K), 1)
            END IF
C                                  ROW ELIMINATION WITH COLUMN INDEXING
            CALL ZSWAP (N-K, FAC(K,K+1), LDFAC, FAC(L,K+1), LDFAC)
            CALL ZGERU (N-K, N-K, DCMPLX(1.0D0,0.0D0), FAC(K+1,K), 1,
     &                  FAC(K,K+1), LDFAC, FAC(K+1,K+1), LDFAC)
         ELSE
            INFO = K
         END IF
   40 CONTINUE
      IPVT(N) = N
      IF (CABS1(FAC(N,N)) .LT. SMALL) INFO = N
C
      IF (INFO .NE. 0) THEN
         CALL E1MES (4, 2, 'The input matrix is singular.  '//
     &               'Some of the diagonal elements of the upper '//
     &               'triangular matrix U of the LU factorization '//
     &               'are close to zero.')
      END IF
C
 9000 CALL E1POP ('DL2TCG ')
      RETURN
      END
