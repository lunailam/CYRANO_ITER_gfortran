C-----------------------------------------------------------------------
C  IMSL Name:  E1POP
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 13, 1984
C
C  Purpose:    To pop a subroutine name from the error control stack.
C
C  Usage:      CALL E1POP(NAME)
C
C  Arguments:
C     NAME   - A character string of length six specifying the name
C              of the subroutine.  (Input)
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1POP (NAME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      CHARACTER  NAME*(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IERTYP, IR
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
      EXTERNAL   E1MES, E1PRT, E1PSH, E1STI, E1STL, I1KRL
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1KST
      INTEGER    I1KST
C
      IF (CALLVL .LE. 1) THEN
         CALL E1PSH ('E1POP ')
         CALL E1STL (1, NAME)
         CALL E1MES (5, 1, 'Error condition in E1POP.  Cannot pop '//
     &               'from %(L1) because stack is empty.')
         STOP
      ELSE IF (NAME .NE. RNAME(CALLVL)) THEN
         CALL E1STL (1, NAME)
         CALL E1STL (2, RNAME(CALLVL))
         CALL E1MES (5, 2, 'Error condition in E1POP.  %(L1) does '//
     &               'not match the name %(L2) in the stack.')
         STOP
      ELSE
         IERTYP = ERTYPE(CALLVL)
         IF (IERTYP .NE. 0) THEN
C                                  M1VE ERROR TYPE AND ERROR CODE TO
C                                    PREVIOUS LEVEL FOR ERROR TYPES 2-7
            IF (IERTYP.GE.2 .AND. IERTYP.LE.7) THEN
               ERTYPE(CALLVL-1) = ERTYPE(CALLVL)
               ERCODE(CALLVL-1) = ERCODE(CALLVL)
            END IF
C                                  CHECK PRINT TABLE TO DETERMINE
C                                    WHETHER TO PRINT STORED MESSAGE
            IF (IERTYP .LE. 4) THEN
               IF (ISUSER(CALLVL-1) .AND. PRINTB(IERTYP).EQ.1)
     &             CALL E1PRT
            ELSE
               IF (PRINTB(IERTYP) .EQ. 1) CALL E1PRT
            END IF
C                                  CHECK STOP TABLE AND ERROR TYPE TO
C                                    DETERMINE WHETHER TO STOP
            IF (IERTYP .LE. 4) THEN
               IF (ISUSER(CALLVL-1) .AND. STOPTB(IERTYP).EQ.1) THEN
                  STOP
               END IF
            ELSE IF (IERTYP .EQ. 5) THEN
               IF (STOPTB(IERTYP) .EQ. 1) THEN
                  STOP
               END IF
            ELSE IF (HDRFMT(IERTYP) .EQ. 1) THEN
               IF (ISUSER(CALLVL-1)) THEN
                  IF (N1RGB(0) .NE. 0) THEN
                     STOP
                  END IF
               END IF
            END IF
         END IF
C                                  SET ERROR TYPE AND CODE
         IF (CALLVL .LT. MAXLEV) THEN
            ERTYPE(CALLVL+1) = -1
            ERCODE(CALLVL+1) = -1
         END IF
C                                  SET IR = AMOUNT OF WORKSPACE
C                                  ALLOCATED AT THIS LEVEL
         IR = I1KST(1) - IALLOC(CALLVL-1)
         IF (IR .GT. 0) THEN
C                                  RELEASE WORKSPACE
            CALL I1KRL (IR)
            IALLOC(CALLVL) = 0
         ELSE IF (IR .LT. 0) THEN
            CALL E1STI (1, CALLVL)
            CALL E1STI (2, IALLOC(CALLVL-1))
            CALL E1STI (3, I1KST(1))
            CALL E1MES (5, 3, 'Error condition in E1POP. '//
     &                  ' The number of workspace allocations at '//
     &                  'level %(I1) is %(I2).  However, the total '//
     &                  'number of workspace allocations is %(I3).')
            STOP
         END IF
C                                  DECREASE THE STACK POINTER BY ONE
         CALLVL = CALLVL - 1
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1PSH
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 2, 1984
C
C  Purpose:    To push a subroutine name onto the error control stack.
C
C  Usage:      CALL E1PSH(NAME)
C
C  Arguments:
C     NAME   - A character string of length six specifing the name of
C              the subroutine.  (Input)
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1PSH (NAME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      CHARACTER  NAME*(*)
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    IFINIT
      SAVE       IFINIT
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
      EXTERNAL   E1INIT, E1MES, E1STI
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1KST
      INTEGER    I1KST
C
      DATA IFINIT/0/
C                                  INITIALIZE ERROR TABLE IF NECESSARY
      IF (IFINIT .EQ. 0) THEN
         CALL E1INIT
         IFINIT = 1
      END IF
      IF (CALLVL .GE. MAXLEV) THEN
         CALL E1STI (1, MAXLEV)
         CALL E1MES (5, 1, 'Error condition in E1PSH.  Push would '//
     &               'cause stack level to exceed %(I1). ')
         STOP
      ELSE
C                                  STORE ALLOCATION LEVEL
         IALLOC(CALLVL) = I1KST(1)
C                                  INCREMENT THE STACK POINTER BY ONE
         CALLVL = CALLVL + 1
C                                  PUT SUBROUTINE NAME INTO STACK
         RNAME(CALLVL) = NAME
C                                  SET ERROR TYPE AND ERROR CODE
         ERTYPE(CALLVL) = 0
         ERCODE(CALLVL) = 0
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  I1KRL
C
C  Computer:   VAX/SINGLE
C
C  Revised:    August 9, 1983
C
C  Purpose:    Deallocate the last N allocations made in the workspace.
C              stack by I1KGT
C
C  Usage:      CALL I1KRL(N)
C
C  Arguments:
C     N      - Number of allocations to be released top down (Input)
C
C  Copyright:  1983 by IMSL, Inc.  All Rights Reserved
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE I1KRL (N)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IN, LALC, LBND, LBOOK, LMAX, LNEED, LNOW, LOUT,
     &           LUSED, NDX, NEXT
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      LOGICAL    FIRST
      SAVE       FIRST
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                                  SPECIFICATIONS FOR COMMON /WORKSP/
      REAL       RWKSP(5000)
      REAL       RDWKSP(5000)
      DOUBLE PRECISION DWKSP(2500)
      COMPLEX    CWKSP(2500)
      COMPLEX    CZWKSP(2500)
      COMPLEX    *16 ZWKSP(1250)
      INTEGER    IWKSP(5000)
      LOGICAL    LWKSP(5000)
      EQUIVALENCE (DWKSP(1), RWKSP(1))
      EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
      EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
      EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
      COMMON     /WORKSP/ DWKSP
C                                  SPECIFICATIONS FOR EQUIVALENCE
      EQUIVALENCE (LOUT, IWKSP(1))
      EQUIVALENCE (LNOW, IWKSP(2))
      EQUIVALENCE (LUSED, IWKSP(3))
      EQUIVALENCE (LBND, IWKSP(4))
      EQUIVALENCE (LMAX, IWKSP(5))
      EQUIVALENCE (LALC, IWKSP(6))
      EQUIVALENCE (LNEED, IWKSP(7))
      EQUIVALENCE (LBOOK, IWKSP(8))
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1STI, IWKIN
C
      DATA FIRST/.TRUE./
C
      IF (FIRST) THEN
C                                  INITIALIZE WORKSPACE IF NEEDED
         FIRST = .FALSE.
         CALL IWKIN (0)
      END IF
C                                  CALLING I1KRL(0) WILL CONFIRM
C                                  INTEGRITY OF SYSTEM AND RETURN
      IF (N .LT. 0) THEN
         CALL E1MES (5, 10, 'Error from subroutine I1KRL:  Attempt'//
     &               ' to release a negative number of workspace'//
     &               ' allocations. ')
         GO TO 9000
      END IF
C                                  BOOKKEEPING OVERWRITTEN
      IF (LNOW.LT.LBOOK .OR. LNOW.GT.LUSED .OR. LUSED.GT.LMAX .OR.
     &    LNOW.GE.LBND .OR. LOUT.GT.LALC) THEN
         CALL E1MES (5, 11, 'Error from subroutine I1KRL:  One or '//
     &               'more of the first eight bookkeeping locations '//
     &               'in IWKSP have been overwritten.  ')
         GO TO 9000
      END IF
C                                  CHECK ALL THE POINTERS IN THE
C                                  PERMANENT STORAGE AREA.  THEY MUST
C                                  BE MONOTONE INCREASING AND LESS THAN
C                                  OR EQUAL TO LMAX, AND THE INDEX OF
C                                  THE LAST POINTER MUST BE LMAX+1.
      NDX = LBND
      IF (NDX .NE. LMAX+1) THEN
         DO 10  I=1, LALC
            NEXT = IWKSP(NDX)
            IF (NEXT .EQ. LMAX+1) GO TO 20
C
            IF (NEXT.LE.NDX .OR. NEXT.GT.LMAX) THEN
               CALL E1MES (5, 12, 'Error from subroutine I1KRL:  '//
     &                     'A pointer in permanent storage has been '//
     &                     ' overwritten. ')
               GO TO 9000
            END IF
            NDX = NEXT
   10    CONTINUE
         CALL E1MES (5, 13, 'Error from subroutine I1KRL:  A '//
     &               'pointer in permanent storage has been '//
     &               'overwritten. ')
         GO TO 9000
      END IF
   20 IF (N .GT. 0) THEN
         DO 30  IN=1, N
            IF (LNOW .LE. LBOOK) THEN
               CALL E1MES (5, 14, 'Error from subroutine I1KRL:  '//
     &                     'Attempt to release a nonexistant '//
     &                     'workspace  allocation. ')
               GO TO 9000
            ELSE IF (IWKSP(LNOW).LT.LBOOK .OR. IWKSP(LNOW).GE.LNOW-1)
     &              THEN
C                                  CHECK TO MAKE SURE THE BACK POINTERS
C                                  ARE MONOTONE.
               CALL E1STI (1, LNOW)
               CALL E1MES (5, 15, 'Error from subroutine I1KRL:  '//
     &                     'The pointer at IWKSP(%(I1)) has been '//
     &                     'overwritten.  ')
               GO TO 9000
            ELSE
               LOUT = LOUT - 1
               LNOW = IWKSP(LNOW)
            END IF
   30    CONTINUE
      END IF
C
 9000 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  N1RGB
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 2, 1984
C
C  Purpose:    Return a positive number as a flag to indicated that a
C              stop should occur due to one or more global errors.
C
C  Usage:      N1RGB(IDUMMY)
C
C  Arguments:
C     IDUMMY - Integer scalar dummy argument.
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION N1RGB (IDUMMY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDUMMY
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
C                                  INITIALIZE FUNCTION
      N1RGB = 0
C                                  CHECK FOR GLOBAL ERROR TYPE 6
      IF (IFERR6 .GT. 0) THEN
         N1RGB = STOPTB(6)
         IFERR6 = 0
      END IF
C                                  CHECK FOR GLOBAL ERROR TYPE 7
      IF (IFERR7 .GT. 0) THEN
         N1RGB = STOPTB(7)
         IFERR7 = 0
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1STL
C
C  Computer:   VAX/SINGLE
C
C  Revised:    November 8, 1985
C
C  Purpose:    To store a string for subsequent use within an error
C              message.
C
C  Usage:      CALL E1STL(IL,STRING)
C
C  Arguments:
C     IL     - Integer specifying the substitution index.  IL must be
C              between 1 and 9.  (Input)
C     STRING - A character string.  (Input)
C
C  Copyright:  1985 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1STL (IL, STRING)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IL
      CHARACTER  STRING*(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, LEN2
      CHARACTER  STRGUP(255)
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    IFINIT
      SAVE       IFINIT
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
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  IABS,LEN,MIN0
      INTRINSIC  IABS, LEN, MIN0
      INTEGER    IABS, LEN, MIN0
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1INIT, E1INPL
C
      DATA IFINIT/0/
C                                  INITIALIZE IF NECESSARY
      IF (IFINIT .EQ. 0) THEN
         CALL E1INIT
         IFINIT = 1
      END IF
      LEN2 = LEN(STRING)
      LEN2 = MIN0(LEN2,255)
      DO 10  I=1, LEN2
         STRGUP(I) = STRING(I:I)
   10 CONTINUE
      IF (IABS(IL).GE.1 .AND. IABS(IL).LE.9) THEN
         CALL E1INPL ('L', IL, LEN2, STRGUP)
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  I1KST
C
C  Computer:   VAX/SINGLE
C
C  Revised:    August 9, 1983
C
C  Purpose:    Return control information about the workspace stack.
C
C  Usage:      I1KST(NFACT)
C
C  Arguments:
C     NFACT  - Integer value between 1 and 6 inclusive returns the
C                 following information: (Input)
C                   NFACT = 1 - LOUT: number of current allocations
C                               excluding permanent storage. At the
C                               end of a run, there should be no
C                               active allocations.
C                   NFACT = 2 - LNOW: current active length
C                   NFACT = 3 - LTOTAL: total storage used thus far
C                   NFACT = 4 - LMAX: maximum storage allowed
C                   NFACT = 5 - LALC: total number of allocations made
C                               by I1KGT thus far
C                   NFACT = 6 - LNEED: number of numeric storage units
C                               by which the stack size must be
C                               increased for all past allocations
C                               to succeed
C     I1KST  - Integer function. (Output) Returns a workspace stack
C              statistic according to value of NFACT.
C
C  Copyright:  1983 by IMSL, Inc.  All Rights Reserved
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION I1KST (NFACT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NFACT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    ISTATS(7)
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      LOGICAL    FIRST
      SAVE       FIRST
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                                  SPECIFICATIONS FOR COMMON /WORKSP/
      REAL       RWKSP(5000)
      REAL       RDWKSP(5000)
      DOUBLE PRECISION DWKSP(2500)
      COMPLEX    CWKSP(2500)
      COMPLEX    CZWKSP(2500)
      COMPLEX    *16 ZWKSP(1250)
      INTEGER    IWKSP(5000)
      LOGICAL    LWKSP(5000)
      EQUIVALENCE (DWKSP(1), RWKSP(1))
      EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
      EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
      EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
      COMMON     /WORKSP/ DWKSP
C                                  SPECIFICATIONS FOR EQUIVALENCE
      EQUIVALENCE (ISTATS(1), IWKSP(1))
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, IWKIN
C
      DATA FIRST/.TRUE./
C
      IF (FIRST) THEN
C                                  INITIALIZE WORKSPACE IF NEEDED
         FIRST = .FALSE.
         CALL IWKIN (0)
      END IF
C
      IF (NFACT.LE.0 .OR. NFACT.GE.7) THEN
         CALL E1MES (5, 9, 'Error from subroutine I1KST:  Argument'//
     &               ' for I1KST must be between 1 and 6 inclusive.')
      ELSE IF (NFACT .EQ. 1) THEN
C                                  LOUT
         I1KST = ISTATS(1)
      ELSE IF (NFACT .EQ. 2) THEN
C                                  LNOW + PERMANENT
         I1KST = ISTATS(2) + (ISTATS(5)-ISTATS(4)+1)
      ELSE IF (NFACT .EQ. 3) THEN
C                                  LUSED + PERMANENT
         I1KST = ISTATS(3) + (ISTATS(5)-ISTATS(4)+1)
      ELSE IF (NFACT .EQ. 4) THEN
C                                  LMAX
         I1KST = ISTATS(5)
      ELSE IF (NFACT .EQ. 5) THEN
C                                  LALC
         I1KST = ISTATS(6)
      ELSE IF (NFACT .EQ. 6) THEN
C                                  LNEED
         I1KST = ISTATS(7)
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  IWKIN (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    January 17, 1984
C
C  Purpose:    Initialize bookkeeping locations describing the
C              workspace stack.
C
C  Usage:      CALL IWKIN (NSU)
C
C  Argument:
C     NSU    - Number of numeric storage units to which the workspace
C              stack is to be initialized
C
C  GAMS:       N4
C
C  Chapters:   MATH/LIBRARY Reference Material
C              STAT/LIBRARY Reference Material
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IWKIN (NSU)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NSU
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    ISIZE(6), LALC, LBND, LBOOK, LMAX, LNEED, LNOW, LOUT,
     &           LUSED, MELMTS, MTYPE
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      LOGICAL    FIRST
      SAVE       FIRST
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                                  SPECIFICATIONS FOR COMMON /WORKSP/
      REAL       RWKSP(5000)
      REAL       RDWKSP(5000)
      DOUBLE PRECISION DWKSP(2500)
      COMPLEX    CWKSP(2500)
      COMPLEX    CZWKSP(2500)
      COMPLEX    *16 ZWKSP(1250)
      INTEGER    IWKSP(5000)
      LOGICAL    LWKSP(5000)
      EQUIVALENCE (DWKSP(1), RWKSP(1))
      EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
      EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
      EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
      COMMON     /WORKSP/ DWKSP
C                                  SPECIFICATIONS FOR EQUIVALENCE
      EQUIVALENCE (LOUT, IWKSP(1))
      EQUIVALENCE (LNOW, IWKSP(2))
      EQUIVALENCE (LUSED, IWKSP(3))
      EQUIVALENCE (LBND, IWKSP(4))
      EQUIVALENCE (LMAX, IWKSP(5))
      EQUIVALENCE (LALC, IWKSP(6))
      EQUIVALENCE (LNEED, IWKSP(7))
      EQUIVALENCE (LBOOK, IWKSP(8))
      EQUIVALENCE (ISIZE(1), IWKSP(11))
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  MAX0
      INTRINSIC  MAX0
      INTEGER    MAX0
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1STI
C
      DATA FIRST/.TRUE./
C
      IF (.NOT.FIRST) THEN
         IF (NSU .NE. 0) THEN
            CALL E1STI (1, LMAX)
            CALL E1MES (5, 100, 'Error from subroutine IWKIN:  '//
     &                  'Workspace stack has previously been '//
     &                  'initialized to %(I1). Correct by making the '//
     &                  'call to IWKIN the first executable '//
     &                  'statement in the main program.  ')
C
            STOP
C
         ELSE
            RETURN
         END IF
      END IF
C
      IF (NSU .EQ. 0) THEN
C                                  IF NSU=0 USE DEFAULT SIZE 5000
         MELMTS = 5000
      ELSE
         MELMTS = NSU
      END IF
C                                  NUMBER OF ITEMS .LT. 0
      IF (MELMTS .LE. 0) THEN
         CALL E1STI (1, MELMTS)
         CALL E1MES (5, 1, 'Error from subroutine IWKIN:  Number '//
     &               'of numeric storage units is not positive. NSU '//
     &               '= %(I1) ')
      ELSE
C
         FIRST = .FALSE.
C                                  HERE TO INITIALIZE
C
C                                  SET DATA SIZES APPROPRIATE FOR A
C                                  STANDARD CONFORMING FORTRAN SYSTEM
C                                  USING THE FORTRAN
C                                  *NUMERIC STORAGE UNIT* AS THE
C                                  MEASURE OF SIZE.
C
C                                  TYPE IS REAL
         MTYPE = 3
C                                  LOGICAL
         ISIZE(1) = 1
C                                  INTEGER
         ISIZE(2) = 1
C                                  REAL
         ISIZE(3) = 1
C                                  DOUBLE PRECISION
         ISIZE(4) = 2
C                                  COMPLEX
         ISIZE(5) = 2
C                                  DOUBLE COMPLEX
         ISIZE(6) = 4
C                                  NUMBER OF WORDS USED FOR BOOKKEEPING
         LBOOK = 16
C                                  CURRENT ACTIVE LENGTH OF THE STACK
         LNOW = LBOOK
C                                  MAXIMUM VALUE OF LNOW ACHIEVED THUS
C                                  FAR
         LUSED = LBOOK
C                                  MAXIMUM LENGTH OF THE STORAGE ARRAY
         LMAX = MAX0(MELMTS,((LBOOK+2)*ISIZE(2)+ISIZE(3)-1)/ISIZE(3))
C                                  LOWER BOUND OF THE PERMANENT STORAGE
C                                  WHICH IS ONE WORD MORE THAN THE
C                                  MAXIMUM ALLOWED LENGTH OF THE STACK
         LBND = LMAX + 1
C                                  NUMBER OF CURRENT ALLOCATIONS
         LOUT = 0
C                                  TOTAL NUMBER OF ALLOCATIONS MADE
         LALC = 0
C                                  NUMBER OF WORDS BY WHICH THE ARRAY
C                                  SIZE MUST BE INCREASED FOR ALL PAST
C                                  ALLOCATIONS TO SUCCEED
         LNEED = 0
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1MES
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 2, 1984
C
C  Purpose:    Set an error state for the current level in the stack.
C              The message is printed immediately if the error type is
C              5, 6, or 7 and the print attribute for that type is YES.
C
C  Usage:      CALL E1MES(IERTYP,IERCOD,MSGPKD)
C
C  Arguments:
C     IERTYP - Integer specifying the error type.  (Input)
C                IERTYP=1,  informational/note
C                IERTYP=2,  informational/alert
C                IERTYP=3,  informational/warning
C                IERTYP=4,  informational/fatal
C                IERTYP=5,  terminal
C                IERTYP=6,  PROTRAN/warning
C                IERTYP=7,  PROTRAN/fatal
C     IERCOD - Integer specifying the error code.  (Input)
C     MSGPKD - A character string containing the message.
C              (Input)  Within the message, any of following may appear
C                %(A1),%(A2),...,%(A9) for character arrays
C                %(C1),%(C2),...,%(C9) for complex numbers
C                %(D1),%(D2),...,%(D9) for double precision numbers
C                %(I1),%(I2),...,%(I9) for integer numbers
C                %(K1),%(K2),...,%(K9) for keywords
C                %(L1),%(L2),...,%(L9) for literals (strings)
C                %(R1),%(R2),...,%(R9) for real numbers
C                %(Z1),%(Z2),...,%(Z9) for double complex numbers
C              This provides a way to insert character arrays, strings,
C              numbers, and keywords into the message.  See remarks
C              below.
C
C  Remarks:
C     The number of characters in the message after the insertion of
C     the corresponding strings, etc. should not exceed 255.  If the
C     limit is exceeded, only the first 255 characters will be used.
C     The appropriate strings, etc. need to have been previously stored
C     in common via calls to E1STA, E1STD, etc.  Line breaks may be
C     specified by inserting the two characters '%/' into the message
C     at the desired locations.
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1MES (IERTYP, IERCOD, MSGPKD)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IERTYP, IERCOD
      CHARACTER  MSGPKD*(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    ERTYP2, I, IER, IPLEN, ISUB, LAST, LEN2, LOC, M, MS,
     &           NLOC, NUM, PBEG
      CHARACTER  MSGTMP(255)
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    IFINIT, NFORMS
      CHARACTER  BLNK, DBB(3), FIND(4), FORMS(9), INREF(25), LPAR,
     &           NCHECK(3), PERCNT, RPAR
      SAVE       BLNK, DBB, FIND, FORMS, IFINIT, INREF, LPAR, NCHECK,
     &           NFORMS, PERCNT, RPAR
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
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  LEN,MIN0
      INTRINSIC  LEN, MIN0
      INTEGER    LEN, MIN0
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   C1TCI, E1INIT, E1PRT, E1UCS, M1VE, M1VECH
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1DX
      INTEGER    I1DX
C
      DATA FORMS/'A', 'C', 'D', 'I', 'K', 'L', 'R', 'S', 'Z'/,
     &     NFORMS/9/
      DATA PERCNT/'%'/, LPAR/'('/, RPAR/')'/, BLNK/' '/
      DATA INREF/' ', 'i', 'n', ' ', 'r', 'e', 'f', 'e', 'r',
     &     'e', 'n', 'c', 'e', ' ', 't', 'o', ' ', 'k', 'e',
     &     'y', 'w', 'o', 'r', 'd', ' '/
      DATA NCHECK/'N', '1', '*'/, DBB/'.', ' ', ' '/
      DATA FIND/'*', ' ', ' ', '*'/
      DATA IFINIT/0/
C                                  INITIALIZE ERROR TABLE IF NECESSARY
      IF (IFINIT .EQ. 0) THEN
         CALL E1INIT
         IFINIT = 1
      END IF
C                                  CHECK AND SET ERROR TYPE IF NECESSARY
      IF (IERTYP .NE. -1) THEN
         ERTYPE(CALLVL) = IERTYP
      ELSE IF (IERTYP.LT.-1 .OR. IERTYP.GT.7) THEN
         MSGLEN = 51
         CALL M1VECH ('.  Error from E1MES.  Illegal error type'//
     &                ' specified. ', MSGLEN, MSGSAV, MSGLEN)
         CALL E1PRT
         STOP
      END IF
C
      ERTYP2 = ERTYPE(CALLVL)
C                                  SET ERROR CODE IF NECESSARY
      IF (IERCOD .GT. -1) ERCODE(CALLVL) = IERCOD
      LEN2 = LEN(MSGPKD)
C
      IF (IERTYP.EQ.0 .OR. IERCOD.EQ.0) THEN
C                                  REMOVE THE ERROR STATE
         MSGLEN = 0
      ELSE IF (LEN2.EQ.0 .OR. (LEN2.EQ.1.AND.MSGPKD(1:1).EQ.BLNK)) THEN
         IF (ERTYP2 .EQ. 6) IFERR6 = 1
         IF (ERTYP2 .EQ. 7) IFERR7 = 1
C                                  UPDATE CHECKSUM PARAMETER ERCKSM
         CALL E1UCS
C                                  PRINT MESSAGE IF NECESSARY
         IF (ERTYP2.GE.5 .AND. PRINTB(ERTYP2).EQ.1) CALL E1PRT
      ELSE
C                                  FILL UP MSGSAV WITH EXPANDED MESSAGE
         LEN2 = MIN0(LEN2,255)
         DO 10  I=1, LEN2
            MSGTMP(I) = MSGPKD(I:I)
   10    CONTINUE
         MS = 0
         M = 0
C                                  CHECK PLIST FOR KEYWORD NAME
         NLOC = I1DX(PLIST,PLEN,NCHECK,3)
         IF (NLOC.GT.0 .AND. HDRFMT(ERTYP2).EQ.3) THEN
C                                  M1VE INREF INTO MSGSAV
            CALL M1VE (INREF, 1, 25, 25, MSGSAV, 1, 25, 25, IER)
C                                  GET LENGTH OF KEYWORD NAME
            CALL C1TCI (PLIST(NLOC+3), 3, IPLEN, IER)
            PBEG = NLOC + 3 + IER
C                                  M1VE KEYWORD NAME INTO MSGSAV
            CALL M1VE (PLIST, PBEG, PBEG+IPLEN-1, PLEN, MSGSAV, 26,
     &                 IPLEN+25, 255, IER)
C                                  UPDATE POINTER
            MS = IPLEN + 25
         END IF
C                                  INSERT DOT, BLANK, BLANK
         CALL M1VE (DBB, 1, 3, 3, MSGSAV, MS+1, MS+3, 255, IER)
         MS = MS + 3
C                                  LOOK AT NEXT CHARACTER
   20    M = M + 1
         ISUB = 0
         IF (M .GT. LEN2-4) THEN
            LAST = LEN2 - M + 1
            DO 30  I=1, LAST
   30       MSGSAV(MS+I) = MSGTMP(M+I-1)
            MSGLEN = MS + LAST
            GO TO 40
         ELSE IF (MSGTMP(M).EQ.PERCNT .AND. MSGTMP(M+1).EQ.LPAR .AND.
     &           MSGTMP(M+4).EQ.RPAR) THEN
            CALL C1TCI (MSGTMP(M+3), 1, NUM, IER)
            IF (IER.EQ.0 .AND. NUM.NE.0 .AND. I1DX(FORMS,NFORMS,
     &          MSGTMP(M+2),1).NE.0) THEN
C                                  LOCATE THE ITEM IN THE PARAMETER LIST
               CALL M1VE (MSGTMP(M+2), 1, 2, 2, FIND, 2, 3, 4, IER)
               LOC = I1DX(PLIST,PLEN,FIND,4)
               IF (LOC .GT. 0) THEN
C                                  SET IPLEN = LENGTH OF STRING
                  CALL C1TCI (PLIST(LOC+4), 4, IPLEN, IER)
                  PBEG = LOC + 4 + IER
C                                  ADJUST IPLEN IF IT IS TOO BIG
                  IPLEN = MIN0(IPLEN,255-MS)
C                                  M1VE STRING FROM PLIST INTO MSGSAV
                  CALL M1VE (PLIST, PBEG, PBEG+IPLEN-1, PLEN, MSGSAV,
     &                       MS+1, MS+IPLEN, 255, IER)
                  IF (IER.GE.0 .AND. IER.LT.IPLEN) THEN
C                                  UPDATE POINTERS
                     M = M + 4
                     MS = MS + IPLEN - IER
C                                  BAIL OUT IF NO MORE ROOM
                     IF (MS .GE. 255) THEN
                        MSGLEN = 255
                        GO TO 40
                     END IF
C                                  SET FLAG TO SHOW SUBSTITION WAS MADE
                     ISUB = 1
                  END IF
               END IF
            END IF
         END IF
         IF (ISUB .EQ. 0) THEN
            MS = MS + 1
            MSGSAV(MS) = MSGTMP(M)
         END IF
         GO TO 20
   40    ERTYP2 = ERTYPE(CALLVL)
         IF (ERTYP2 .EQ. 6) IFERR6 = 1
         IF (ERTYP2 .EQ. 7) IFERR7 = 1
C                                  UPDATE CHECKSUM PARAMETER ERCKSM
         CALL E1UCS
C                                  PRINT MESSAGE IF NECESSARY
         IF (ERTYP2.GE.5 .AND. PRINTB(ERTYP2).EQ.1) CALL E1PRT
      END IF
C                                  CLEAR PARAMETER LIST
      PLEN = 1
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1STI
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 6, 1984
C
C  Purpose:    To store an integer for subsequent use within an error
C              message.
C
C  Usage:      CALL E1STI(II, IVALUE)
C
C  Arguments:
C     II     - Integer specifying the substitution index.  II must be
C              between 1 and 9.  (Input)
C     IVALUE - The integer to be stored.  (Input)
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1STI (II, IVALUE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    II, IVALUE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IBEG, IER, ILEN
      CHARACTER  ARRAY(14)
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
      EXTERNAL   C1TIC, E1INIT, E1INPL
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
      CALL C1TIC (IVALUE, ARRAY, 14, IER)
      IBEG = I1ERIF(ARRAY,14,BLANK,1)
      IF (II.GE.1 .AND. II.LE.9 .AND. IER.EQ.0) THEN
         ILEN = 15 - IBEG
         CALL E1INPL ('I', II, ILEN, ARRAY(IBEG))
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1INIT
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 13, 1984
C
C  Purpose:    Initialization.
C
C  Usage:      CALL E1INIT
C
C  Arguments:  None
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1INIT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    L
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    ISINIT
      SAVE       ISINIT
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
C                              SPECIFICATIONS FOR COMMON /ERCOM8/
      INTEGER    PROLVL, XXLINE(10), XXPLEN(10), ICALOC(10), INALOC(10)
      COMMON     /ERCOM8/ PROLVL, XXLINE, XXPLEN, ICALOC, INALOC
      SAVE       /ERCOM8/
C                              SPECIFICATIONS FOR COMMON /ERCOM9/
      CHARACTER  XXPROC(10)*31
      COMMON     /ERCOM9/ XXPROC
      SAVE       /ERCOM9/
C
      DATA ISINIT/0/
C
      IF (ISINIT .EQ. 0) THEN
C                                  INITIALIZE
         CALLVL = 1
         ERCODE(1) = 0
         ERTYPE(1) = 0
         IALLOC(1) = 0
         ISUSER(1) = .TRUE.
         IFERR6 = 0
         IFERR7 = 0
         PLEN = 1
         MAXLEV = 50
         DO 10  L=2, 51
            ERTYPE(L) = -1
            ERCODE(L) = -1
            IALLOC(L) = 0
            ISUSER(L) = .FALSE.
   10    CONTINUE
         DO 20  L=1, 7
            HDRFMT(L) = 1
            TRACON(L) = 1
   20    CONTINUE
         PROLVL = 1
         DO 30  L=1, 10
   30    ICALOC(L) = 0
         XXLINE(1) = 0
         XXPLEN(1) = 1
         XXPROC(1) = '?'
         RNAME(1) = 'USER'
         PRINTB(1) = 0
         PRINTB(2) = 0
         DO 40  L=3, 7
   40    PRINTB(L) = 1
         STOPTB(1) = 0
         STOPTB(2) = 0
         STOPTB(3) = 0
         STOPTB(4) = 1
         STOPTB(5) = 1
         STOPTB(6) = 0
         STOPTB(7) = 1
         ERCKSM = 0.0D0
C                                  SET FLAG TO INDICATE THAT
C                                    INITIALIZATION HAS OCCURRED
         ISINIT = 1
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1PRT
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 14, 1984
C
C  Purpose:    To print an error message.
C
C  Usage:      CALL E1PRT
C
C  Arguments:  None
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1PRT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    ALL, I, IBEG, IBLOC, IBLOC2, IEND, IER, IHDR, J,
     &           LERTYP, LOC, LOCM1, LOCX, MAXLOC, MAXTMP, MLOC, MOD,
     &           NCBEG, NLOC, NOUT
      CHARACTER  MSGTMP(70), STRING(10)
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      CHARACTER  ATLINE(9), BLANK(1), DBB(3), FROM(6), MSGTYP(8,7),
     &           PERSLA(2), QMARK, UNKNOW(8)
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
C                              SPECIFICATIONS FOR COMMON /ERCOM8/
      INTEGER    PROLVL, XXLINE(10), XXPLEN(10), ICALOC(10), INALOC(10)
      COMMON     /ERCOM8/ PROLVL, XXLINE, XXPLEN, ICALOC, INALOC
      SAVE       /ERCOM8/
C                              SPECIFICATIONS FOR COMMON /ERCOM9/
      CHARACTER  XXPROC(10)*31
      COMMON     /ERCOM9/ XXPROC
      SAVE       /ERCOM9/
      SAVE       ATLINE, BLANK, DBB, FROM, MSGTYP, PERSLA, QMARK,
     &           UNKNOW
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  MIN0
      INTRINSIC  MIN0
      INTEGER    MIN0
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   C1TIC, M1VE, UMACH
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1DX, I1ERIF
      INTEGER    I1DX, I1ERIF
C
      DATA MSGTYP/'N', 'O', 'T', 'E', ' ', ' ', ' ', ' ', 'A',
     &     'L', 'E', 'R', 'T', ' ', ' ', ' ', 'W', 'A', 'R',
     &     'N', 'I', 'N', 'G', ' ', 'F', 'A', 'T', 'A', 'L',
     &     ' ', ' ', ' ', 'T', 'E', 'R', 'M', 'I', 'N', 'A',
     &     'L', 'W', 'A', 'R', 'N', 'I', 'N', 'G', ' ', 'F',
     &     'A', 'T', 'A', 'L', ' ', ' ', ' '/
      DATA UNKNOW/'U', 'N', 'K', 'N', 'O', 'W', 'N', ' '/
      DATA ATLINE/' ', 'a', 't', ' ', 'l', 'i', 'n', 'e', ' '/
      DATA BLANK/' '/, FROM/' ', 'f', 'r', 'o', 'm', ' '/
      DATA DBB/'.', ' ', ' '/, PERSLA/'%', '/'/
      DATA QMARK/'?'/
C
      IF (MSGLEN .LE. 0) RETURN
      CALL UMACH (2, NOUT)
      MAXTMP = 70
      MOD = 0
      LERTYP = ERTYPE(CALLVL)
      IHDR = HDRFMT(LERTYP)
      IF (IHDR .EQ. 3) THEN
         IF (XXPROC(PROLVL)(1:1).EQ.QMARK .AND. XXLINE(PROLVL).EQ.0)
     &       THEN
            IHDR = 1
         END IF
      END IF
      IEND = 0
      IF (IHDR.EQ.1 .AND. ERTYPE(CALLVL).LE.4) THEN
         MSGTMP(1) = BLANK(1)
         IEND = 1
C                                  CONVERT ERROR CODE INTO CHAR STRING
         CALL C1TIC (ERCODE(CALLVL), STRING, 10, IER)
C                                  LOCATE START OF NON-BLANK CHARACTERS
         IBEG = I1ERIF(STRING,10,BLANK,1)
C                                  M1VE IT TO MSGTMP
         CALL M1VE (STRING, IBEG, 10, 10, MSGTMP, IEND+1,
     &              IEND+11-IBEG, MAXTMP, IER)
         IEND = IEND + 11 - IBEG
      END IF
      IF (IHDR .NE. 2) THEN
         CALL M1VE (FROM, 1, 6, 6, MSGTMP, IEND+1, IEND+6, MAXTMP, IER)
         IEND = IEND + 6
      END IF
      IF (IHDR .EQ. 3) THEN
C                                  THIS IS A PROTRAN RUN TIME ERROR MSG.
C                                  RETRIEVE THE PROCEDURE NAME
         CALL M1VE (XXPROC(PROLVL), 1, XXPLEN(PROLVL), 31, MSGTMP,
     &              IEND+1, IEND+XXPLEN(PROLVL), MAXTMP, IER)
         MLOC = IEND + XXPLEN(PROLVL) + 1
         MSGTMP(MLOC) = BLANK(1)
         IEND = IEND + I1DX(MSGTMP(IEND+1),XXPLEN(PROLVL)+1,BLANK,1) -
     &          1
         IF (XXLINE(PROLVL) .GT. 0) THEN
C                                  INSERT ATLINE
            CALL M1VE (ATLINE, 1, 9, 9, MSGTMP, IEND+1, IEND+9,
     &                 MAXTMP, IER)
            IEND = IEND + 9
C                                  CONVERT PROTRAN GLOBAL LINE NUMBER
            CALL C1TIC (XXLINE(PROLVL), STRING, 10, IER)
C                                  LOCATE START OF NON-BLANK CHARACTERS
            IBEG = I1ERIF(STRING,10,BLANK,1)
C                                  M1VE GLOBAL LINE NUMBER TO MSGTMP
            CALL M1VE (STRING, IBEG, 10, 10, MSGTMP, IEND+1,
     &                 IEND+11-IBEG, MAXTMP, IER)
            IEND = IEND + 11 - IBEG
         END IF
      ELSE
C                                  THIS IS EITHER A LIBRARY ERROR MSG
C                                  OR A PROTRAN PREPROCESSOR ERROR MSG
         IF (IHDR .EQ. 1) THEN
C                                  THIS IS A LIBRARY ERROR MESSAGE.
C                                  RETRIEVE ROUTINE NAME
            CALL M1VE (RNAME(CALLVL), 1, 6, 6, MSGTMP, IEND+1, IEND+6,
     &                 MAXTMP, IER)
            MSGTMP(IEND+7) = BLANK(1)
            IEND = IEND + I1DX(MSGTMP(IEND+1),7,BLANK,1) - 1
         END IF
C                                  ADD DOT, BLANK, BLANK IF NEEDED
         IF (I1DX(MSGSAV,3,DBB,3) .NE. 1) THEN
            CALL M1VE (DBB, 1, 3, 3, MSGTMP, IEND+1, IEND+3, MAXTMP,
     &                 IER)
            IEND = IEND + 3
            MOD = 3
         END IF
      END IF
C                                  MSGTMP AND MSGSAV NOW CONTAIN THE
C                                   ERROR MESSAGE IN FINAL FORM.
      NCBEG = 59 - IEND - MOD
      ALL = 0
      IBLOC = I1DX(MSGSAV,MSGLEN,PERSLA,2)
      IF (IBLOC.NE.0 .AND. IBLOC.LT.NCBEG) THEN
         LOCM1 = IBLOC - 1
         LOC = IBLOC + 1
      ELSE IF (MSGLEN .LE. NCBEG) THEN
         LOCM1 = MSGLEN
         ALL = 1
      ELSE
         LOC = NCBEG
C                                  CHECK FOR APPROPRIATE PLACE TO SPLIT
   10    CONTINUE
         IF (MSGSAV(LOC) .NE. BLANK(1)) THEN
            LOC = LOC - 1
            IF (LOC .GT. 1) GO TO 10
            LOC = NCBEG + 1
         END IF
         LOCM1 = LOC - 1
      END IF
C                                  NO BLANKS FOUND IN FIRST NCBEG CHARS
      IF (LERTYP.GE.1 .AND. LERTYP.LE.7) THEN
         WRITE (NOUT,99995) (MSGTYP(I,LERTYP),I=1,8),
     &                     (MSGTMP(I),I=1,IEND), (MSGSAV(I),I=1,LOCM1)
      ELSE
         WRITE (NOUT,99995) (UNKNOW(I),I=1,8), (MSGTMP(I),I=1,IEND),
     &                     (MSGSAV(I),I=1,LOCM1)
      END IF
      IF (ALL .EQ. 0) THEN
C                                  PREPARE TO WRITE CONTINUATION OF
C                                    MESSAGE
C
C                                  FIND WHERE TO BREAK MESSAGE
C                                    LOC = NUMBER OF CHARACTERS OF
C                                          MESSAGE WRITTEN SO FAR
   20    LOCX = LOC + 64
         NLOC = LOC + 1
         IBLOC2 = IBLOC
         MAXLOC = MIN0(MSGLEN-LOC,64)
         IBLOC = I1DX(MSGSAV(NLOC),MAXLOC,PERSLA,2)
         IF (MSGSAV(NLOC).EQ.BLANK(1) .AND. IBLOC2.EQ.0) NLOC = NLOC +
     &       1
         IF (IBLOC .GT. 0) THEN
C                                  PAGE BREAK FOUND AT IBLOC
            LOCX = NLOC + IBLOC - 2
            WRITE (NOUT,99996) (MSGSAV(I),I=NLOC,LOCX)
            LOC = NLOC + IBLOC
            GO TO 20
C                                  DON'T BOTHER LOOKING FOR BLANK TO
C                                    BREAK AT IF LOCX .GE. MSGLEN
         ELSE IF (LOCX .LT. MSGLEN) THEN
C                                  CHECK FOR BLANK TO BREAK THE LINE
   30       CONTINUE
            IF (MSGSAV(LOCX) .EQ. BLANK(1)) THEN
C                                  BLANK FOUND AT LOCX
               WRITE (NOUT,99996) (MSGSAV(I),I=NLOC,LOCX)
               LOC = LOCX
               GO TO 20
            END IF
            LOCX = LOCX - 1
            IF (LOCX .GT. NLOC) GO TO 30
            LOCX = LOC + 64
C                                  NO BLANKS FOUND IN NEXT 64 CHARS
            WRITE (NOUT,99996) (MSGSAV(I),I=NLOC,LOCX)
            LOC = LOCX
            GO TO 20
         ELSE
C                                  ALL THE REST WILL FIT ON 1 LINE
            LOCX = MSGLEN
            WRITE (NOUT,99996) (MSGSAV(I),I=NLOC,LOCX)
         END IF
      END IF
C                                  SET LENGTH OF MSGSAV AND PLEN
C                                    TO SHOW THAT MESSAGE HAS
C                                    ALREADY BEEN PRINTED
 9000 MSGLEN = 0
      PLEN = 1
      IF (TRACON(LERTYP).EQ.1 .AND. CALLVL.GT.2) THEN
C                                  INITIATE TRACEBACK
         WRITE (NOUT,99997)
         DO 9005  J=CALLVL, 1, -1
            IF (J .GT. 1) THEN
               IF (ISUSER(J-1)) THEN
                  WRITE (NOUT,99998) RNAME(J), ERTYPE(J), ERCODE(J)
               ELSE
                  WRITE (NOUT,99999) RNAME(J), ERTYPE(J), ERCODE(J)
               END IF
            ELSE
               WRITE (NOUT,99998) RNAME(J), ERTYPE(J), ERCODE(J)
            END IF
 9005    CONTINUE
      END IF
C
      RETURN
99995 FORMAT (/, ' *** ', 8A1, ' ERROR', 59A1)
99996 FORMAT (' *** ', 9X, 64A1)
99997 FORMAT (14X, 'Here is a traceback of subprogram calls',
     &       ' in reverse order:', /, 14X, '      Routine    Error ',
     &       'type    Error code', /, 14X, '      -------    ',
     &       '----------    ----------')
99998 FORMAT (20X, A6, 5X, I6, 8X, I6)
99999 FORMAT (20X, A6, 5X, I6, 8X, I6, 4X, '(Called internally)')
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1INPL
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 2, 1984
C
C  Purpose:    To store a character string in the parameter list PLIST
C              for use by the error message handler.
C
C  Usage:      CALL E1INPL(FORM,NUM,SLEN,STRUP)
C
C  Arguments:
C     FORM   - A character string of length one to be inserted into
C              PLIST which specifies the form of the string.  (Input)
C              For example, 'L' for string, 'A' for character array,
C              'I' for integer, 'K' for keyword (PROTRAN only).  An
C              asterisk is inserted into PLIST preceding FORM.
C     NUM    - Integer to be inserted as a character into PLIST
C              immediately following FORM.  (Input)  NUM must be between
C              1 and 9.
C     SLEN   - The number of characters in STRUP.  (Input)  LEN must be
C              less than or equal to 255.  The character representation
C              of SLEN is inserted into PLIST after NUM and an asterisk.
C     STRUP  - A character string of length LEN which is to be inserted
C              into PLIST.  (Input)  Trailing blanks are ignored.
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1INPL (FORM, NUM, SLEN, STRUP)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NUM, SLEN
      CHARACTER  FORM, STRUP(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IER, L, LEN2, LENCK, LOC, NLEN, NNUM
      CHARACTER  STRNCH(3)
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      CHARACTER  BLANK, PRCNT(1), TEMP(4)
      SAVE       BLANK, PRCNT, TEMP
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
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  IABS
      INTRINSIC  IABS
      INTEGER    IABS
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   C1TIC, M1VE
C
      DATA TEMP/'*', ' ', ' ', '*'/, PRCNT/'%'/, BLANK/' '/
C
      NNUM = IABS(NUM)
      LENCK = PLEN + SLEN + 8
      IF (NNUM.GE.1 .AND. NNUM.LE.9 .AND. LENCK.LE.300) THEN
         TEMP(2) = FORM
         CALL C1TIC (NNUM, TEMP(3), 1, IER)
         LOC = PLEN + 1
         IF (LOC .EQ. 2) LOC = 1
         CALL M1VE (TEMP, 1, 4, 4, PLIST(LOC), 1, 4, 262, IER)
         LOC = LOC + 4
         IF (NUM .LT. 0) THEN
            LEN2 = SLEN
         ELSE
            DO 10  L=1, SLEN
               LEN2 = SLEN - L + 1
               IF (STRUP(LEN2) .NE. BLANK) GO TO 20
   10       CONTINUE
            LEN2 = 1
   20       CONTINUE
         END IF
         NLEN = 1
         IF (LEN2 .GE. 10) NLEN = 2
         IF (LEN2 .GE. 100) NLEN = 3
         CALL C1TIC (LEN2, STRNCH, NLEN, IER)
         CALL M1VE (STRNCH, 1, NLEN, 3, PLIST(LOC), 1, NLEN, 262, IER)
         LOC = LOC + NLEN
         CALL M1VE (PRCNT, 1, 1, 1, PLIST(LOC), 1, 1, 262, IER)
         LOC = LOC + 1
         CALL M1VE (STRUP, 1, LEN2, LEN2, PLIST(LOC), 1, LEN2, 262,
     &              IER)
         PLEN = LOC + LEN2 - 1
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1UCS
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 8, 1984
C
C  Purpose:    To update the checksum number for error messages.
C
C  Usage:      CALL E1UCS
C
C  Arguments:  None
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1UCS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IBEG, IBEG2, IEND, ILOC, IPOS, JLOC, NCODE, NLEN
      DOUBLE PRECISION DNUM
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      DOUBLE PRECISION DMAX
      CHARACTER  BLANK(1), COMMA(1), EQUAL(1), LPAR(1)
      SAVE       BLANK, COMMA, DMAX, EQUAL, LPAR
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
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DMOD
      INTRINSIC  DMOD
      DOUBLE PRECISION DMOD
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   S1ANUM
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   ICASE, I1X
      INTEGER    ICASE, I1X
C
      DATA BLANK(1)/' '/, COMMA(1)/','/, LPAR(1)/'('/
      DATA EQUAL(1)/'='/, DMAX/1.0D+9/
C
      IF (MSGLEN .GT. 1) THEN
         IPOS = 0
         IBEG2 = 1
   10    IBEG = IBEG2
         IEND = MSGLEN
C                                  LOOK FOR BLANK, COMMA, LEFT PAREN.,
C                                  OR EQUAL SIGN
         ILOC = I1X(MSGSAV(IBEG),IEND-IBEG+1,BLANK,1)
         JLOC = I1X(MSGSAV(IBEG),IEND-IBEG+1,COMMA,1)
         IF (ILOC.EQ.0 .OR. (JLOC.GT.0.AND.JLOC.LT.ILOC)) ILOC = JLOC
         JLOC = I1X(MSGSAV(IBEG),IEND-IBEG+1,LPAR,1)
         IF (ILOC.EQ.0 .OR. (JLOC.GT.0.AND.JLOC.LT.ILOC)) ILOC = JLOC
         JLOC = I1X(MSGSAV(IBEG),IEND-IBEG+1,EQUAL,1)
         IF (ILOC.EQ.0 .OR. (JLOC.GT.0.AND.JLOC.LT.ILOC)) ILOC = JLOC
         IF (ILOC .GE. 1) THEN
            CALL S1ANUM (MSGSAV(IBEG+ILOC), IEND-IBEG-ILOC+1, NCODE,
     &                   NLEN)
            IF (NCODE.EQ.2 .OR. NCODE.EQ.3) THEN
C                                  FLOATING POINT NUMBER FOUND.
C                                  SET POINTERS TO SKIP OVER IT
               IBEG2 = IBEG + ILOC + NLEN
               IF (IBEG2 .LE. MSGLEN) THEN
                  CALL S1ANUM (MSGSAV(IBEG2), IEND-IBEG2+1, NCODE,
     &                         NLEN)
                  IF ((MSGSAV(IBEG2).EQ.'+'.OR.MSGSAV(IBEG2).EQ.
     &                '-') .AND. NCODE.EQ.1) THEN
C                                  INTEGER IMMEDIATELY FOLLOWS A REAL AS
C                                  WITH SOME CDC NOS. LIKE 1.2345678+123
C                                  SET POINTERS TO SKIP OVER IT
                     IBEG2 = IBEG2 + NLEN
                  END IF
               END IF
            ELSE
               IBEG2 = IBEG + ILOC
            END IF
            IEND = IBEG + ILOC - 1
         END IF
C                                  UPDATE CKSUM USING PART OF MESSAGE
         DO 20  I=IBEG, IEND
            IPOS = IPOS + 1
            DNUM = ICASE(MSGSAV(I))
            ERCKSM = DMOD(ERCKSM+DNUM*IPOS,DMAX)
   20    CONTINUE
C                                  GO BACK FOR MORE IF NEEDED
         IF (IEND.LT.MSGLEN .AND. IBEG2.LT.MSGLEN) GO TO 10
C                                  UPDATE CKSUM USING ERROR TYPE
         DNUM = ERTYPE(CALLVL)
         ERCKSM = DMOD(ERCKSM+DNUM*(IPOS+1),DMAX)
C                                  UPDATE CKSUM USING ERROR CODE
         DNUM = ERCODE(CALLVL)
         ERCKSM = DMOD(ERCKSM+DNUM*(IPOS+2),DMAX)
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  M1VECH
C
C  Computer:   VAX/SINGLE
C
C  Revised:    December 31, 1984
C
C  Purpose:    Character substring assignment.
C
C  Usage:      CALL M1VECH (STR1, LEN1, STR2, LEN2)
C
C  Arguments:
C     STR1   - Source substring.  (Input)
C              The source substring is STR1(1:LEN1).
C     LEN1   - Length of STR1.  (Input)
C     STR2   - Destination substring.  (Output)
C              The destination substring is STR2(1:LEN2).
C     LEN2   - Length of STR2.  (Input)
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE M1VECH (STR1, LEN1, STR2, LEN2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    LEN1, LEN2
      CHARACTER  STR1*(*), STR2*(*)
C
      STR2(1:LEN2) = STR1(1:LEN1)
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  C1TCI
C
C  Computer:   VAX/SINGLE
C
C  Revised:    August 13, 1984
C
C  Purpose:    Convert character string into corresponding integer
C              form.
C
C  Usage:      CALL C1TCI (CHRSTR, SLEN, NUM, IER)
C
C  Arguments:
C   CHRSTR  - Character array that contains the number description.
C             (Input)
C   SLEN    - Length of the character array.  (Input)
C   NUM     - The answer.  (Output)
C   IER     - Completion code.  (Output)  Where
C                IER =-2  indicates that the number is too large to
C                         be converted;
C                IER =-1  indicates that SLEN <= 0;
C                IER = 0  indicates normal completion;
C                IER > 0  indicates that the input string contains a
C                         nonnumeric character.  IER is the index of
C                         the first nonnumeric character in CHRSTR.
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE C1TCI (CHRSTR, SLEN, NUM, IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    SLEN, NUM, IER
      CHARACTER  CHRSTR(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    COUNT, I, IMACH5, J, N, S, SIGN
      CHARACTER  ZERO
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      CHARACTER  BLANK, DIGIT*10, MINUS, PLUS
      SAVE       BLANK, DIGIT, MINUS, PLUS
C                                  SPECIFICATIONS FOR EQUIVALENCE
      EQUIVALENCE (DIGIT, ZERO)
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  INDEX
      INTRINSIC  INDEX
      INTEGER    INDEX
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   IMACH
      INTEGER    IMACH
C
      DATA DIGIT/'0123456789'/
      DATA BLANK/' '/, MINUS/'-'/, PLUS/'+'/
C
C                                  CHECK SLEN
      NUM = 0
      IF (SLEN .LE. 0) THEN
         IER = -1
         GO TO 50
      END IF
C                                  HANDLE LEADING BLANKS
      SIGN = 1
      I = 1
   10 IF (I .LE. SLEN) THEN
         IF (CHRSTR(I) .EQ. BLANK) THEN
            I = I + 1
            GO TO 10
         END IF
      ELSE
         IER = 1
         GO TO 50
      END IF
C                                  CHECK FOR SIGN, IF ANY
      S = I
      IF (CHRSTR(I) .EQ. MINUS) THEN
         SIGN = -1
         I = I + 1
      ELSE IF (CHRSTR(I) .EQ. PLUS) THEN
         I = I + 1
      END IF
   20 IF (I .LE. SLEN) THEN
         IF (CHRSTR(I) .EQ. BLANK) THEN
            I = I + 1
            GO TO 20
         END IF
      ELSE
         IER = S
         GO TO 50
      END IF
C                                  SKIP LEADING ZERO
      J = I
   30 IF (I .LE. SLEN) THEN
         IF (CHRSTR(I) .EQ. ZERO) THEN
            I = I + 1
            GO TO 30
         END IF
      ELSE
         IER = 0
         GO TO 50
      END IF
C                                  CHECK FIRST NONBLANK CHARACTER
      COUNT = 0
C                                  CHECK NUMERIC CHARACTERS
      IMACH5 = IMACH(5)
   40 N = INDEX(DIGIT,CHRSTR(I))
      IF (N .NE. 0) THEN
         COUNT = COUNT + 1
         IF (NUM .GT. ((IMACH5-N)+1)/10) THEN
            IER = -2
            GO TO 50
         ELSE
            NUM = NUM*10 - 1 + N
            I = I + 1
            IF (I .LE. SLEN) GO TO 40
         END IF
      END IF
C
      IF (COUNT .EQ. 0) THEN
         IF (I .GT. J) THEN
            IER = I
         ELSE
            IER = S
         END IF
      ELSE IF (I .GT. SLEN) THEN
         NUM = SIGN*NUM
         IER = 0
      ELSE
         NUM = SIGN*NUM
         IER = I
      END IF
C
   50 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  C1TIC
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 9, 1984
C
C  Purpose:    Convert an integer to its corresponding character form.
C              (Right justified)
C
C  Usage:      CALL C1TIC(NUM, CHRSTR, SLEN, IER)
C
C  Arguments:
C     NUM    - Integer number.  (Input)
C     CHRSTR - Character array that receives the result.  (Output)
C     SLEN   - Length of the character array.  (Input)
C     IER    - Completion code.  (Output) Where
C                 IER < 0  indicates that SLEN <= 0,
C                 IER = 0  indicates normal completion,
C                 IER > 0  indicates that the character array is too
C                       small to hold the complete number.  IER
C                       indicates how many significant digits are
C                       being truncated.
C
C  Remarks:
C  1. The character array is filled in a right justified manner.
C  2. Leading zeros are replaced by blanks.
C  3. Sign is inserted only for negative number.
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE C1TIC (NUM, CHRSTR, SLEN, IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NUM, SLEN, IER
      CHARACTER  CHRSTR(SLEN)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, J, K, L
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      CHARACTER  BLANK(1), DIGIT(10), MINUS(1)
      SAVE       BLANK, DIGIT, MINUS
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  IABS
      INTRINSIC  IABS
      INTEGER    IABS
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   M1VE
C
      DATA DIGIT/'0', '1', '2', '3', '4', '5', '6', '7', '8',
     &     '9'/
      DATA BLANK/' '/, MINUS/'-'/
C                                  CHECK SLEN
      IF (SLEN .LE. 0) THEN
         IER = -1
         RETURN
      END IF
C                                  THE NUMBER IS ZERO
      IF (NUM .EQ. 0) THEN
         CALL M1VE (BLANK, 1, 1, 1, CHRSTR, 1, SLEN-1, SLEN, I)
         CHRSTR(SLEN) = DIGIT(1)
         IER = 0
         RETURN
      END IF
C                                  CONVERT NUMBER DIGIT BY DIGIT TO
C                                  CHARACTER FORM
      J = SLEN
      K = IABS(NUM)
   10 IF (K.GT.0 .AND. J.GE.1) THEN
         L = K
         K = K/10
         L = L - K*10
         CHRSTR(J) = DIGIT(L+1)
         J = J - 1
         GO TO 10
      END IF
C
   20 IF (K .EQ. 0) THEN
         IF (NUM .LT. 0) THEN
            CALL M1VE (MINUS, 1, 1, 1, CHRSTR, J, J, SLEN, I)
            IF (I .NE. 0) THEN
               IER = 1
               RETURN
            END IF
            J = J - 1
         END IF
         IER = 0
         CALL M1VE (BLANK, 1, 1, 1, CHRSTR, 1, J, SLEN, I)
         RETURN
      END IF
C                                  DETERMINE THE NUMBER OF SIGNIFICANT
C                                  DIGITS BEING TRUNCATED
      I = 0
   30 IF (K .GT. 0) THEN
         K = K/10
         I = I + 1
         GO TO 30
      END IF
C
      IF (NUM .LT. 0) I = I + 1
      IER = I
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  I1ERIF
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 13, 1984
C
C  Purpose:    Return the position of the first element of a given
C              character array which is not an element of another
C              character array.
C
C  Usage:      I1ERIF(STR1, LEN1, STR2, LEN2)
C
C  Arguments:
C     STR1   - Character array to be searched.  (Input)
C     LEN1   - Length of STR1.  (Input)
C     STR2   - Character array to be searched for.  (Input)
C     LEN2   - Length of STR2.  (Input)
C     I1ERIF - Integer function.  (Output)
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION I1ERIF (STR1, LEN1, STR2, LEN2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    LEN1, LEN2
      CHARACTER  STR1(LEN1), STR2(LEN2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1X
      INTEGER    I1X
C                              FIRST EXECUTABLE STATEMENT
      IF (LEN1.LE.0 .OR. LEN2.LE.0) THEN
         I1ERIF = 1
      ELSE
         DO 10  I=1, LEN1
            IF (I1X(STR2,LEN2,STR1(I),1) .EQ. 0) THEN
               I1ERIF = I
               RETURN
            END IF
   10    CONTINUE
         I1ERIF = 0
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  IMACH (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 26, 1984
C
C  Purpose:    Retrieve integer machine constants.
C
C  Usage:      IMACH(N)
C
C  Arguments:
C     N      - Index of desired constant.  (Input)
C     IMACH  - Machine constant.  (Output)
C
C  Remark:
C     Following is a description of the assorted integer machine
C     constants.
C
C     Words
C
C        IMACH( 1) = Number of bits per integer storage unit.
C        IMACH( 2) = Number of characters per integer storage unit.
C
C     Integers
C
C        Assume integers are represented in the S-DIGIT, BASE-A form
C        SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C        where 0 .LE. X(I) .LT. A for I=0,...,S-1.  Then
C
C        IMACH( 3) = A, the base.
C        IMACH( 4) = S, number of BASE-A digits.
C        IMACH( 5) = A**S - 1, largest magnitude.
C
C     Floating-point numbers
C
C        Assume floating-point numbers are represented in the T-DIGIT,
C        BASE-B form SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C        where 0 .LE. X(I) .LT. B for I=1,...,T,
C        0 .LT. X(1), and EMIN .LE. E .LE. EMAX.  Then
C
C        IMACH( 6) = B, the base.
C
C        Single precision
C
C           IMACH( 7) = T, number of BASE-B digits.
C           IMACH( 8) = EMIN, smallest exponent E.
C           IMACH( 9) = EMAX, largest exponent E.
C
C        Double precision
C
C           IMACH(10) = T, number of BASE-B digits.
C           IMACH(11) = EMIN, smallest exponent E.
C           IMACH(12) = EMAX, largest exponent E.
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
      INTEGER FUNCTION IMACH (N)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    NOUT
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    IMACHV(12)
      SAVE       IMACHV
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   UMACH
C                                  DEFINE CONSTANTS
      DATA IMACHV(1)/32/
      DATA IMACHV(2)/4/
      DATA IMACHV(3)/2/
      DATA IMACHV(4)/31/
      DATA IMACHV(5)/2147483647/
      DATA IMACHV(6)/2/
      DATA IMACHV(7)/24/
      DATA IMACHV(8)/-127/
      DATA IMACHV(9)/127/
      DATA IMACHV(10)/56/
      DATA IMACHV(11)/-127/
      DATA IMACHV(12)/127/
C
      IF (N.LT.1 .OR. N.GT.12) THEN
C                                  ERROR.  INVALID RANGE FOR N.
         CALL UMACH (2, NOUT)
         WRITE (NOUT,99999) N
99999    FORMAT (/, ' *** TERMINAL ERROR 5 from IMACH.  The argument',
     &          /, ' ***          must be between 1 and 12 inclusive.'
     &          , /, ' ***          N = ', I6, '.', /)
         IMACH = 0
         STOP
C
      ELSE
         IMACH = IMACHV(N)
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  I1DX (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    September 9, 1985
C
C  Purpose:    Determine the array subscript indicating the starting
C              element at which a key character sequence begins.
C              (Case-insensitive version)
C
C  Usage:      I1DX(CHRSTR, I1LEN, KEY, KLEN)
C
C  Arguments:
C     CHRSTR - Character array to be searched.  (Input)
C     I1LEN  - Length of CHRSTR.  (Input)
C     KEY    - Character array that contains the key sequence.  (Input)
C     KLEN   - Length of KEY.  (Input)
C     I1DX   - Integer function.  (Output)
C
C  Remarks:
C  1. Returns zero when there is no match.
C
C  2. Returns zero if KLEN is longer than ISLEN.
C
C  3. Returns zero when any of the character arrays has a negative or
C     zero length.
C
C  GAMS:       N5c
C
C  Chapter:    MATH/LIBRARY Utilities
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION I1DX (CHRSTR, I1LEN, KEY, KLEN)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    I1LEN, KLEN
      CHARACTER  CHRSTR(*), KEY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, II, J
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   ICASE, I1CSTR
      INTEGER    ICASE, I1CSTR
C
      I1DX = 0
      IF (KLEN.LE.0 .OR. I1LEN.LE.0) GO TO 9000
      IF (KLEN .GT. I1LEN) GO TO 9000
C
      I = 1
      II = I1LEN - KLEN + 1
   10 IF (I .LE. II) THEN
         IF (ICASE(CHRSTR(I)) .EQ. ICASE(KEY(1))) THEN
            IF (KLEN .NE. 1) THEN
               J = KLEN - 1
               IF (I1CSTR(CHRSTR(I+1),J,KEY(2),J) .EQ. 0) THEN
                  I1DX = I
                  GO TO 9000
               END IF
            ELSE
               I1DX = I
               GO TO 9000
            END IF
         END IF
         I = I + 1
         GO TO 10
      END IF
C
 9000 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  S1ANUM
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 28, 1984
C
C  Purpose:    Scan a token and identify it as follows: integer, real
C              number (single/double), FORTRAN relational operator,
C              FORTRAN logical operator, or FORTRAN logical constant.
C
C  Usage:      CALL S1ANUM(INSTR, SLEN, CODE, OLEN)
C
C  Arguments:
C     INSTR  - Character string to be scanned.  (Input)
C     SLEN   - Length of INSTR.  (Input)
C     CODE   - Token code.  (Output)  Where
C                 CODE =  0  indicates an unknown token,
C                 CODE =  1  indicates an integer number,
C                 CODE =  2  indicates a (single precision) real number,
C                 CODE =  3  indicates a (double precision) real number,
C                 CODE =  4  indicates a logical constant (.TRUE. or
C                               .FALSE.),
C                 CODE =  5  indicates the relational operator .EQ.,
C                 CODE =  6  indicates the relational operator .NE.,
C                 CODE =  7  indicates the relational operator .LT.,
C                 CODE =  8  indicates the relational operator .LE.,
C                 CODE =  9  indicates the relational operator .GT.,
C                 CODE = 10  indicates the relational operator .GE.,
C                 CODE = 11  indicates the logical operator .AND.,
C                 CODE = 12  indicates the logical operator .OR.,
C                 CODE = 13  indicates the logical operator .EQV.,
C                 CODE = 14  indicates the logical operator .NEQV.,
C                 CODE = 15  indicates the logical operator .NOT..
C     OLEN   - Length of the token as counted from the first character
C              in INSTR.  (Output)  OLEN returns a zero for an unknown
C              token (CODE = 0).
C
C  Remarks:
C  1. Blanks are considered significant.
C  2. Lower and upper case letters are not significant.
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE S1ANUM (INSTR, SLEN, CODE, OLEN)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    SLEN, CODE, OLEN
      CHARACTER  INSTR(SLEN)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IBEG, IIBEG, J
      LOGICAL    FLAG
      CHARACTER  CHRSTR(6)
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    TABPTR(16), TDCNST, TICNST, TOKEN(13), TRCNST, TZERR
      CHARACTER  DIGIT(10), LETTER(52), MINUS, PERIOD, PLUS, TABLE(38)
      SAVE       DIGIT, LETTER, MINUS, PERIOD, PLUS, TABLE, TABPTR,
     &           TDCNST, TICNST, TOKEN, TRCNST, TZERR
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1X, I1CSTR
      INTEGER    I1X, I1CSTR
C
      DATA TOKEN/5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 4/
      DATA TABLE/'D', 'E', 'E', 'Q', 'N', 'E', 'L', 'T', 'L',
     &     'E', 'G', 'T', 'G', 'E', 'A', 'N', 'D', 'O', 'R',
     &     'E', 'Q', 'V', 'N', 'E', 'Q', 'V', 'N', 'O', 'T',
     &     'T', 'R', 'U', 'E', 'F', 'A', 'L', 'S', 'E'/
      DATA TABPTR/1, 2, 3, 5, 7, 9, 11, 13, 15, 18, 20, 23, 27, 30,
     &     34, 39/
      DATA DIGIT/'0', '1', '2', '3', '4', '5', '6', '7', '8',
     &     '9'/
      DATA LETTER/'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
     &     'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
     &     'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c',
     &     'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
     &     'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w',
     &     'x', 'y', 'z'/
      DATA PERIOD/'.'/, PLUS/'+'/, MINUS/'-'/
      DATA TZERR/0/, TICNST/1/
      DATA TRCNST/2/, TDCNST/3/
C
      IF (SLEN .LE. 0) THEN
         CODE = 0
         OLEN = 0
         RETURN
      END IF
C                                  STATE 0 - ASSUME ERROR TOKEN
      IBEG = 1
      CODE = TZERR
C                                  CHECK SIGN
      IF (INSTR(IBEG).EQ.MINUS .OR. INSTR(IBEG).EQ.PLUS) THEN
         FLAG = .TRUE.
         IIBEG = IBEG
         IBEG = IBEG + 1
      ELSE
         FLAG = .FALSE.
      END IF
C                                  STATE 1 - ASSUME INTEGER CONSTANT
      IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
         CODE = TICNST
         IIBEG = IBEG
         IBEG = IBEG + 1
C
   10    IF (IBEG .LE. SLEN) THEN
C
            IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
               IIBEG = IBEG
               IBEG = IBEG + 1
               GO TO 10
C
            END IF
C
         ELSE
            GO TO 80
C
         END IF
C
         IF (INSTR(IBEG) .NE. PERIOD) GO TO 80
      END IF
C                                  STATE 2 - ASSUME REAL CONSTANT
      IF (CODE .EQ. TICNST) THEN
         CODE = TRCNST
         IIBEG = IBEG
         IBEG = IBEG + 1
         IF (IBEG .GT. SLEN) GO TO 80
      ELSE IF (INSTR(IBEG).EQ.PERIOD .AND. SLEN.GE.2) THEN
         IF (I1X(DIGIT,10,INSTR(IBEG+1),1) .NE. 0) THEN
            CODE = TRCNST
            IIBEG = IBEG + 1
            IBEG = IBEG + 2
            IF (IBEG .GT. SLEN) GO TO 80
         END IF
      END IF
C
      IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
         CODE = TRCNST
         IIBEG = IBEG
         IBEG = IBEG + 1
C
   20    IF (IBEG .LE. SLEN) THEN
C
            IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
               IIBEG = IBEG
               IBEG = IBEG + 1
               GO TO 20
C
            END IF
C
         ELSE
            GO TO 80
C
         END IF
C
      END IF
C
      IF (CODE .EQ. TZERR) THEN
         IF (INSTR(IBEG) .NE. PERIOD) GO TO 80
         IBEG = IBEG + 1
         IF (IBEG .GT. SLEN) GO TO 80
      END IF
C
      IF (I1X(LETTER,52,INSTR(IBEG),1) .EQ. 0) GO TO 80
      CHRSTR(1) = INSTR(IBEG)
C
      DO 30  I=2, 6
         IBEG = IBEG + 1
         IF (IBEG .GT. SLEN) GO TO 80
         IF (I1X(LETTER,52,INSTR(IBEG),1) .EQ. 0) GO TO 40
         CHRSTR(I) = INSTR(IBEG)
   30 CONTINUE
C
      GO TO 80
C
   40 CONTINUE
C
      DO 50  J=1, 15
         IF (I1CSTR(CHRSTR,I-1,TABLE(TABPTR(J)),TABPTR(J+1)-TABPTR(J))
     &        .EQ. 0) GO TO 60
   50 CONTINUE
C
      GO TO 80
C                                  STATE 4 - LOGICAL OPERATOR
   60 IF (J .GT. 2) THEN
C
         IF (CODE .EQ. TRCNST) THEN
C
            IF (INSTR(IBEG) .EQ. PERIOD) THEN
               CODE = TICNST
               IIBEG = IIBEG - 1
            END IF
C
            GO TO 80
C
         ELSE IF (INSTR(IBEG) .NE. PERIOD) THEN
            GO TO 80
C
         ELSE IF (FLAG) THEN
            GO TO 80
C
         ELSE
            CODE = TOKEN(J-2)
            IIBEG = IBEG
            GO TO 80
C
         END IF
C
      END IF
C                                  STATE 5 - DOUBLE PRECISION CONSTANT
      IF (CODE .NE. TRCNST) GO TO 80
      IF (INSTR(IBEG).EQ.MINUS .OR. INSTR(IBEG).EQ.PLUS) IBEG = IBEG +
     &    1
      IF (IBEG .GT. SLEN) GO TO 80
C
      IF (I1X(DIGIT,10,INSTR(IBEG),1) .EQ. 0) THEN
         GO TO 80
C
      ELSE
         IIBEG = IBEG
         IBEG = IBEG + 1
C
   70    IF (IBEG .LE. SLEN) THEN
C
            IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
               IIBEG = IBEG
               IBEG = IBEG + 1
               GO TO 70
C
            END IF
C
         END IF
C
      END IF
C
      IF (J .EQ. 1) CODE = TDCNST
C
   80 CONTINUE
C
      IF (CODE .EQ. TZERR) THEN
         OLEN = 0
C
      ELSE
         OLEN = IIBEG
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  UMACH (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 21, 1984
C
C  Purpose:    Set or retrieve input or output device unit numbers.
C
C  Usage:      CALL UMACH (N, NUNIT)
C
C  Arguments:
C     N      - Index of desired unit.  (Input)
C              The values of N are defined as follows:
C              N = 1, corresponds to the standard input unit.
C              N = 2, corresponds to the standard output unit.
C     NUNIT  - I/O unit.  (Input or Output)
C              If the value of N is negative, the unit corresponding
C              to the index is reset to the value given in NUNIT.
C              Otherwise, the value corresponding to the index is
C              returned in NUNIT.
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
      SUBROUTINE UMACH (N, NUNIT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, NUNIT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    NN, NOUT
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    UNIT(2)
      SAVE       UNIT
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  IABS
      INTRINSIC  IABS
      INTEGER    IABS
C
      DATA UNIT(1)/5/
      DATA UNIT(2)/6/
C
      NN = IABS(N)
      IF (NN.NE.1 .AND. NN.NE.2) THEN
C                                  ERROR.  INVALID RANGE FOR N.
         NOUT = UNIT(2)
         WRITE (NOUT,99999) NN
99999    FORMAT (/, ' *** TERMINAL ERROR 5 from UMACH.  The absolute',
     &          /, ' ***          value of the index variable must be'
     &          , /, ' ***          1 or 2.  IABS(N) = ', I6,
     &          '.', /)
         STOP
C                                  CHECK FOR RESET OR RETRIEVAL
      ELSE IF (N .LT. 0) THEN
C                                  RESET
         UNIT(NN) = NUNIT
      ELSE
C                                  RETRIEVE
         NUNIT = UNIT(N)
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  I1CSTR (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    September 10, 1985
C
C  Purpose:    Case insensitive comparison of two character arrays.
C
C  Usage:      I1CSTR(STR1, LEN1, STR2, LEN2)
C
C  Arguments:
C     STR1   - First character array.  (Input)
C     LEN1   - Length of STR1.  (Input)
C     STR2   - Second character array.  (Input)
C     LEN2   - Length of STR2.  (Input)
C     I1CSTR - Integer function.  (Output) Where
C              I1CSTR = -1  if STR1 .LT. STR2,
C              I1CSTR =  0  if STR1 .EQ. STR2,
C              I1CSTR =  1  if STR1 .GT. STR2.
C
C  Remarks:
C  1. If the two arrays, STR1 and STR2,  are of unequal length, the
C     shorter array is considered as if it were extended with blanks
C     to the length of the longer array.
C
C  2. If one or both lengths are zero or negative the I1CSTR output is
C     based on comparison of the lengths.
C
C  GAMS:       N5c
C
C  Chapter:    MATH/LIBRARY Utilities
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION I1CSTR (STR1, LEN1, STR2, LEN2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    LEN1, LEN2
      CHARACTER  STR1(LEN1), STR2(LEN2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IC1, IC2, ICB, IS, L, LENM
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  ISIGN,MIN0
      INTRINSIC  ISIGN, MIN0
      INTEGER    ISIGN, MIN0
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   ICASE
      INTEGER    ICASE
C
      IF (LEN1.GT.0 .AND. LEN2.GT.0) THEN
C                                  COMPARE FIRST LENM CHARACTERS
         LENM = MIN0(LEN1,LEN2)
         DO 10  L=1, LENM
            IC1 = ICASE(STR1(L))
            IC2 = ICASE(STR2(L))
            IF (IC1 .NE. IC2) THEN
               I1CSTR = ISIGN(1,IC1-IC2)
               RETURN
            END IF
   10    CONTINUE
      END IF
C                                  COMPARISON BASED ON LENGTH OR
C                                  TRAILING BLANKS
      IS = LEN1 - LEN2
      IF (IS .EQ. 0) THEN
         I1CSTR = 0
      ELSE
         IF (LEN1.LE.0 .OR. LEN2.LE.0) THEN
C                                  COMPARISON BASED ON LENGTH
            I1CSTR = ISIGN(1,IS)
         ELSE
C                                  COMPARISON BASED ON TRAILING BLANKS
C                                  TO EXTEND SHORTER ARRAY
            LENM = LENM + 1
            ICB = ICASE(' ')
            IF (IS .GT. 0) THEN
C                                  EXTEND STR2 WITH BLANKS
               DO 20  L=LENM, LEN1
                  IC1 = ICASE(STR1(L))
                  IF (IC1 .NE. ICB) THEN
                     I1CSTR = ISIGN(1,IC1-ICB)
                     RETURN
                  END IF
   20          CONTINUE
            ELSE
C                                  EXTEND STR1 WITH BLANKS
               DO 30  L=LENM, LEN2
                  IC2 = ICASE(STR2(L))
                  IF (ICB .NE. IC2) THEN
                     I1CSTR = ISIGN(1,ICB-IC2)
                     RETURN
                  END IF
   30          CONTINUE
            END IF
C
            I1CSTR = 0
         END IF
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  M1VE
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 5, 1984
C
C  Purpose:    Move a subset of one character array to another.
C
C  Usage:      CALL M1VE(INSTR, INBEG, INEND, INLEN, OUTSTR, OUTBEG,
C                         OUTEND, OUTLEN, IER)
C
C  Arguments:
C     INSTR  - Source character array.  (Input)
C     INBEG  - First element of INSTR to be moved.  (Input)
C     INEND  - Last element of INSTR to be moved.  (Input)
C              The source subset is INSTR(INBEG),...,INSTR(INEND).
C     INLEN  - Length of INSTR.  (Input)
C     OUTSTR - Destination character array.  (Output)
C     IUTBEG - First element of OUTSTR destination.  (Input)
C     IUTEND - Last element of OUTSTR  destination.  (Input)
C              The destination subset is OUTSRT(IUTBEG),...,
C              OUTSTR(IUTEND).
C     IUTLEN - Length of OUTSTR.  (Input)
C     IER    - Completion code.  (Output)
C              IER = -2  indicates that the input parameters, INBEG,
C                        INEND, INLEN, IUTBEG, IUTEND are not
C                        consistent.  One of the conditions
C                        INBEG.GT.0, INEND.GE.INBEG, INLEN.GE.INEND,
C                        IUTBEG.GT.0, or IUTEND.GE.IUTBEG is not
C                        satisfied.
C              IER = -1  indicates that the length of OUTSTR is
C                        insufficient to hold the subset of INSTR.
C                        That is, IUTLEN is less than IUTEND.
C              IER =  0  indicates normal completion
C              IER >  0  indicates that the specified subset of OUTSTR,
C                        OUTSTR(IUTBEG),...,OUTSTR(IUTEND) is not long
C                        enough to hold the subset INSTR(INBEG),...,
C                        INSTR(INEND) of INSTR.  IER is set to the
C                        number of characters that were not moved.
C
C  Remarks:
C  1. If the subset of OUTSTR is longer than the subset of INSTR,
C     trailing blanks are moved to OUTSTR.
C  2. If the subset of INSTR is longer than the subset of OUTSTR,
C     the shorter subset is moved to OUTSTR and IER is set to the number
C     of characters that were not moved to OUTSTR.
C  3. If the length of OUTSTR is insufficient to hold the subset,
C     IER is set to -2 and nothing is moved.
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE M1VE (INSTR, INBEG, INEND, INLEN, OUTSTR, IUTBEG,
     &                 IUTEND, IUTLEN, IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    INBEG, INEND, INLEN, IUTBEG, IUTEND, IUTLEN, IER
      CHARACTER  INSTR(*), OUTSTR(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IUTLAS, KI, KO
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      CHARACTER  BLANK
      SAVE       BLANK
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  MIN0
      INTRINSIC  MIN0
      INTEGER    MIN0
C
      DATA BLANK/' '/
C                                  CHECK INBEG, INEND, INLEN, IUTBEG,
C                                  AND IUTEND
C
      IF (INBEG.LE.0 .OR. INEND.LT.INBEG .OR. INLEN.LT.INEND .OR.
     &    IUTBEG.LE.0 .OR. IUTEND.LT.IUTBEG) THEN
         IER = -2
         RETURN
      ELSE IF (IUTLEN .LT. IUTEND) THEN
         IER = -1
         RETURN
      END IF
C                                  DETERMINE LAST CHARACTER TO M1VE
      IUTLAS = IUTBEG + MIN0(INEND-INBEG,IUTEND-IUTBEG)
C                                  M1VE CHARACTERS
      KI = INBEG
      DO 10  KO=IUTBEG, IUTLAS
         OUTSTR(KO) = INSTR(KI)
         KI = KI + 1
   10 CONTINUE
C                                   SET IER TO NUMBER OF CHARACTERS THAT
C                                   WHERE NOT MOVED
      IER = KI - INEND - 1
C                                   APPEND BLANKS IF NECESSARY
      DO 20  KO=IUTLAS + 1, IUTEND
         OUTSTR(KO) = BLANK
   20 CONTINUE
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  I1X (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    August 30, 1985
C
C  Purpose:    Determine the array subscript indicating the starting
C              element at which a key character sequence begins.
C              (Case-sensitive version)
C
C  Usage:      I1X(CHRSTR, I1LEN, KEY, KLEN)
C
C  Arguments:
C     CHRSTR - Character array to be searched.  (Input)
C     I1LEN  - Length of CHRSTR.  (Input)
C     KEY    - Character array that contains the key sequence.  (Input)
C     KLEN   - Length of KEY.  (Input)
C     I1X    - Integer function.  (Output)
C
C  Remarks:
C  1. Returns zero when there is no match.
C
C  2. Returns zero if KLEN is longer than ISLEN.
C
C  3. Returns zero when any of the character arrays has a negative or
C     zero length.
C
C  GAMS:       N5c
C
C  Chapter:    MATH/LIBRARY Utilities
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION I1X (CHRSTR, I1LEN, KEY, KLEN)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    I1LEN, KLEN
      CHARACTER  CHRSTR(*), KEY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, II, J
C
      I1X = 0
      IF (KLEN.LE.0 .OR. I1LEN.LE.0) GO TO 9000
      IF (KLEN .GT. I1LEN) GO TO 9000
C
      I = 1
      II = I1LEN - KLEN + 1
   10 IF (I .LE. II) THEN
         IF (CHRSTR(I) .EQ. KEY(1)) THEN
            DO 20  J=2, KLEN
               IF (CHRSTR(I+J-1) .NE. KEY(J)) GO TO 30
   20       CONTINUE
            I1X = I
            GO TO 9000
   30       CONTINUE
         END IF
         I = I + 1
         GO TO 10
      END IF
C
 9000 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  ICASE (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    September 9, 1985
C
C  Purpose:    Convert from character to the integer ASCII value without
C              regard to case.
C
C  Usage:      ICASE(CH)
C
C  Arguments:
C     CH     - Character to be converted.  (Input)
C     ICASE  - Integer ASCII value for CH without regard to the case
C              of CH.  (Output)
C              ICASE returns the same value as IMSL routine IACHAR for
C              all but lowercase letters.  For these, it returns the
C              IACHAR value for the corresponding uppercase letter.
C
C  GAMS:       N3
C
C  Chapter:    MATH/LIBRARY Utilities
C              STAT/LIBRARY Utilities
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION ICASE (CH)
C                                  SPECIFICATIONS FOR ARGUMENTS
      CHARACTER  CH
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   IACHAR
      INTEGER    IACHAR
C
      ICASE = IACHAR(CH)
      IF (ICASE.GE.97 .AND. ICASE.LE.122) ICASE = ICASE - 32
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  IACHAR (Single precision version)
C
C  Computer:   VAX/SINGLE
C
C  Revised:    September 9, 1985
C
C  Purpose:    Return the integer ASCII value of a character argument.
C
C  Usage:      IACHAR(CH)
C
C  Arguments:
C     CH     - Character argument for which the integer ASCII value
C              is desired.  (Input)
C     IACHAR - Integer ASCII value for CH.  (Output)
C              The character CH is in the IACHAR-th position of the
C              ASCII collating sequence.
C
C  GAMS:       N3
C
C  Chapter:    MATH/LIBRARY Utilities
C              STAT/LIBRARY Utilities
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION IACHAR (CH)
C                                  SPECIFICATIONS FOR ARGUMENTS
      CHARACTER  CH
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      IACHAR = ICHAR(CH)
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  FFTCI/DFFTCI (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute parameters needed by FFTCF and FFTCB.
C
C  Usage:      CALL FFTCI (N, WFFTC)
C
C  Arguments:
C     N      - Length of the sequence to be transformed.  (Input)
C     WFFTC  - Array of length 4*N+15 containing parameters needed by
C              FFTCF and FFTCB.  (Output)
C
C  Remark:
C     Different WFFTC arrays are needed for different values of N.
C
C  GAMS:       J1a2
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DFFTCI (N, WFFTC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
      DOUBLE PRECISION WFFTC(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DF3TCI
C
C                                  CHECK ARGUMENT N
      IF (N .LT. 1) THEN
         CALL E1PSH ('DFFTCI ')
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'The length of the sequence N = %(I1).  '//
     &               'It must be at least 1.')
         CALL E1POP ('DFFTCI ')
         GO TO 9000
      END IF
C
      IF (N .GT. 1) THEN
         CALL DF3TCI (N, WFFTC(2*N+1), WFFTC(4*N+1))
      END IF
C
 9000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F2TCB/DF2TCB (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the complex periodic sequence from its Fourier
C              coefficients.
C
C  Usage:      CALL F2TCB (N, COEF, SEQ, WFFTC, CPY)
C
C  Arguments:  (See FFTCB)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF2TCB (N, COEF, SEQ, WFFTC, CPY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
      COMPLEX    *16 COEF(*), SEQ(*)
      DOUBLE PRECISION WFFTC(*), CPY(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DCOPY, DF3TCB
C
C                                  CHECK ARGUMENT N
      IF (N .LT. 1) THEN
         CALL E1PSH ('DF2TCB ')
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'The length of the sequence N = %(I1).  '//
     &               'It must be at least 1.')
         CALL E1POP ('DF2TCB ')
         GO TO 9000
      END IF
C                                  COPY COEF TO CPY
      DO 10  I=1, N
         CPY(2*I-1) = DBLE(COEF(I))
         CPY(2*I) = DIMAG(COEF(I))
   10 CONTINUE
C
      IF (N .GT. 1) THEN
         CALL DF3TCB (N, CPY, WFFTC, WFFTC(2*N+1), WFFTC(4*N+1))
      END IF
C                                  COPY CPY TO SEQ
      DO 20  I=1, N
         SEQ(I) = DCMPLX(CPY(2*I-1),CPY(2*I))
   20 CONTINUE
C
 9000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F3TCB/DF3TCB (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the complex periodic sequence from its Fourier
C              coefficients.
C
C  Usage:      CALL F3TCB (N, C, CH, WA, FAC)
C
C  Arguments:
C     N      - Length of the sequence to be transformed.  (Input)
C     C      - Real array of length 2*N, on input containing SEQ,
C              on output containing the Fourier coefficients.
C              (Input/Output)
C     CH     - Real array of length 2*N, needed as workspace.  Will
C              contain the Fourier coefficients, which get copied
C              into C.  (Workspace)
C     WA     - Real array of length 2*N from FFTCI containing powers of
C              e**(2*PI*i/N).  (Input)
C     FAC    - Real array of length 15 from FFTCI containing N (the
C              length of the sequence), NF (the number of factors of
C              N), and the prime factors of N.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF3TCB (N, C, CH, WA, FAC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
      DOUBLE PRECISION C(*), CH(*), WA(*), FAC(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IDL1, IDO, IDOT, IP, IW, IX2, IX3, IX4, K1, L1, L2,
     &           NA, NAC, NF
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  NINT
      INTRINSIC  NINT
      INTEGER    NINT
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   DCOPY, DF4TCB, DF5TCB, DF6TCB, DF7TCB, DF8TCB
C
      NF = NINT(FAC(2))
      NA = 0
      L1 = 1
      IW = 1
      DO 60  K1=1, NF
         IP = NINT(FAC(K1+2))
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO + IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 10
         IX2 = IW + IDOT
         IX3 = IX2 + IDOT
         IF (NA .NE. 0) THEN
            CALL DF6TCB (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
         ELSE
            CALL DF6TCB (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
         END IF
         NA = 1 - NA
         GO TO 50
   10    IF (IP .NE. 2) GO TO 20
         IF (NA .NE. 0) THEN
            CALL DF4TCB (IDOT, L1, CH, C, WA(IW))
         ELSE
            CALL DF4TCB (IDOT, L1, C, CH, WA(IW))
         END IF
         NA = 1 - NA
         GO TO 50
   20    IF (IP .NE. 3) GO TO 30
         IX2 = IW + IDOT
         IF (NA .NE. 0) THEN
            CALL DF5TCB (IDOT, L1, CH, C, WA(IW), WA(IX2))
         ELSE
            CALL DF5TCB (IDOT, L1, C, CH, WA(IW), WA(IX2))
         END IF
         NA = 1 - NA
         GO TO 50
   30    IF (IP .NE. 5) GO TO 40
         IX2 = IW + IDOT
         IX3 = IX2 + IDOT
         IX4 = IX3 + IDOT
         IF (NA .NE. 0) THEN
            CALL DF7TCB (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3),
     &                   WA(IX4))
         ELSE
            CALL DF7TCB (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3),
     &                   WA(IX4))
         END IF
         NA = 1 - NA
         GO TO 50
   40    IF (NA .NE. 0) THEN
            CALL DF8TCB (NAC, IDOT, IP, L1, IDL1, CH, CH, CH, C, C,
     &                   WA(IW))
         ELSE
            CALL DF8TCB (NAC, IDOT, IP, L1, IDL1, C, C, C, CH, CH,
     &                   WA(IW))
         END IF
         IF (NAC .NE. 0) NA = 1 - NA
   50    L1 = L2
         IW = IW + (IP-1)*IDOT
   60 CONTINUE
      IF (NA .NE. 0) CALL DCOPY (2*N, CH, 1, C, 1)
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F3TCI/DF3TCI (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute parameters needed by FFTCI and FFTCB.
C
C  Usage:      CALL F3TCI (N, WA, FAC)
C
C  Arguments:
C     N      - Length of the sequence to be transformed.  (Input)
C     WA     - Vector of length 2*N.  (Output)
C     FAC    - Vector of length 14.  (Output)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF3TCI (N, WA, FAC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
      DOUBLE PRECISION WA(*), FAC(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, I1, IB, IDO, IDOT, II, IP, IPM, J, K1, L1, L2, LD,
     &           NF, NL, NQ, NR, NTRY
      DOUBLE PRECISION ARG, ARGH, ARGLD, FI, TPI
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    NTRYH(4)
      SAVE       NTRYH
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DCOS,DBLE,DSIN
      INTRINSIC  DCOS, DBLE, DSIN
      DOUBLE PRECISION DCOS, DBLE, DSIN
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   DCOPY
C
      DATA NTRYH/3, 4, 2, 5/
C                                  FAC WILL BE A VECTOR OF LENGTH
C                                  14. FAC(1) WILL CONTAIN N,
C                                  FAC(2) THRU FAC(14) WILL HAVE
C                                  THE PRIME FACTORS OF N, IN ASCENDING
C                                  ORDER.
C                                  WA WILL BE A VECTOR OF LENGTH 2*N,
C                                  CONTAINING THE REAL AND IMAGINARY
C                                  PARTS OF CEXP(2*PI*?/N)
      NL = N
      NF = 0
      J = 0
   10 J = J + 1
      IF (J-4) 20,20,30
   20 NTRY = NTRYH(J)
      GO TO 40
   30 NTRY = NTRY + 2
   40 NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR) 10,50,10
   50 NF = NF + 1
      FAC(NF+2) = DBLE(NTRY)
      NL = NQ
      IF (NTRY .NE. 2) GO TO 60
      IF (NF .EQ. 1) GO TO 60
      CALL DCOPY (NF-1, FAC(3), -1, FAC(4), -1)
      FAC(3) = 2.0D0
   60 IF (NL .NE. 1) GO TO 40
      FAC(1) = DBLE(N)
      FAC(2) = DBLE(NF)
      TPI = 2.0D0*3.1415926535897932384626433831D0
      ARGH = TPI/DBLE(N)
      I = 2
      L1 = 1
      DO 90  K1=1, NF
         IP = NINT(FAC(K1+2))
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO + IDO + 2
         IPM = IP - 1
         DO 80  J=1, IPM
            I1 = I
            WA(I-1) = 1.0D0
            WA(I) = 0.0D0
            LD = LD + L1
            FI = 0.0D0
            ARGLD = DBLE(LD)*ARGH
            DO 70  II=4, IDOT, 2
               I = I + 2
               FI = FI + 1.0D0
               ARG = FI*ARGLD
               WA(I-1) = DCOS(ARG)
               WA(I) = DSIN(ARG)
   70       CONTINUE
            IF (IP .GT. 5) THEN
               WA(I1-1) = WA(I-1)
               WA(I1) = WA(I)
            END IF
   80    CONTINUE
         L1 = L2
   90 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F7TCB/DF7TCB (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the complex periodic sequence from its Fourier
C              coefficients.
C
C  Usage:      CALL F7TCB (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C
C  Arguments:
C     IDO    - 2*N divided by present and previous factors.  (Input)
C     L1     - Product of previous factors used.  (Input)
C     CC     - SEQ or CH from F3TCB.  (Input)
C     CH     - A better estimate for COEF.  (Output)  Note:  The roles
C              of CC and CH are reversed when NA=1.
C     WA1    - Real array containing the needed elements of array WA
C              om F3TCB.  (Input)
C     WA2    - Real array containing other needed elements of array WA
C              om F3TCB.  (Input)
C     WA3    - Real array containing other needed elements of array WA
C              om F3TCB.  (Input)
C     WA4    - Real array containing other needed elements of array WA
C              om F3TCB.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF7TCB (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, L1
      DOUBLE PRECISION CC(IDO,5,*), CH(IDO,L1,*), WA1(*), WA2(*),
     &           WA3(*), WA4(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, K
      DOUBLE PRECISION CI2, CI3, CI4, CI5, CR2, CR3, CR4, CR5, DI2,
     &           DI3, DI4, DI5, DR2, DR3, DR4, DR5, TI2, TI3, TI4,
     &           TI5, TR2, TR3, TR4, TR5
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      DOUBLE PRECISION TI11, TI12, TR11, TR12
      SAVE       TI11, TI12, TR11, TR12
C
      DATA TR11/0.309016994374947424102293417183D0/
      DATA TI11/0.951056516295153572116439333379D0/
      DATA TR12/-0.809016994374947424102293417183D0/
      DATA TI12/0.587785252292473129168705954639D0/
C
      IF (IDO .EQ. 2) THEN
         DO 10  K=1, L1
            TI5 = CC(2,2,K) - CC(2,5,K)
            TI2 = CC(2,2,K) + CC(2,5,K)
            TI4 = CC(2,3,K) - CC(2,4,K)
            TI3 = CC(2,3,K) + CC(2,4,K)
            TR5 = CC(1,2,K) - CC(1,5,K)
            TR2 = CC(1,2,K) + CC(1,5,K)
            TR4 = CC(1,3,K) - CC(1,4,K)
            TR3 = CC(1,3,K) + CC(1,4,K)
            CH(1,K,1) = CC(1,1,K) + TR2 + TR3
            CH(2,K,1) = CC(2,1,K) + TI2 + TI3
            CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(2,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(2,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            CH(1,K,2) = CR2 - CI5
            CH(1,K,5) = CR2 + CI5
            CH(2,K,2) = CI2 + CR5
            CH(2,K,3) = CI3 + CR4
            CH(1,K,3) = CR3 - CI4
            CH(1,K,4) = CR3 + CI4
            CH(2,K,4) = CI3 - CR4
            CH(2,K,5) = CI2 - CR5
   10    CONTINUE
      ELSE
         DO 30  K=1, L1
            DO 20  I=2, IDO, 2
               TI5 = CC(I,2,K) - CC(I,5,K)
               TI2 = CC(I,2,K) + CC(I,5,K)
               TI4 = CC(I,3,K) - CC(I,4,K)
               TI3 = CC(I,3,K) + CC(I,4,K)
               TR5 = CC(I-1,2,K) - CC(I-1,5,K)
               TR2 = CC(I-1,2,K) + CC(I-1,5,K)
               TR4 = CC(I-1,3,K) - CC(I-1,4,K)
               TR3 = CC(I-1,3,K) + CC(I-1,4,K)
               CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
               CH(I,K,1) = CC(I,1,K) + TI2 + TI3
               CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
               CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
               CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
               CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
               CR5 = TI11*TR5 + TI12*TR4
               CI5 = TI11*TI5 + TI12*TI4
               CR4 = TI12*TR5 - TI11*TR4
               CI4 = TI12*TI5 - TI11*TI4
               DR3 = CR3 - CI4
               DR4 = CR3 + CI4
               DI3 = CI3 + CR4
               DI4 = CI3 - CR4
               DR5 = CR2 + CI5
               DR2 = CR2 - CI5
               DI5 = CI2 - CR5
               DI2 = CI2 + CR5
               CH(I-1,K,2) = WA1(I-1)*DR2 - WA1(I)*DI2
               CH(I,K,2) = WA1(I-1)*DI2 + WA1(I)*DR2
               CH(I-1,K,3) = WA2(I-1)*DR3 - WA2(I)*DI3
               CH(I,K,3) = WA2(I-1)*DI3 + WA2(I)*DR3
               CH(I-1,K,4) = WA3(I-1)*DR4 - WA3(I)*DI4
               CH(I,K,4) = WA3(I-1)*DI4 + WA3(I)*DR4
               CH(I-1,K,5) = WA4(I-1)*DR5 - WA4(I)*DI5
               CH(I,K,5) = WA4(I-1)*DI5 + WA4(I)*DR5
   20       CONTINUE
   30    CONTINUE
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F5TCB/DF5TCB (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the complex periodic sequence from its Fourier
C              coefficients.
C
C  Usage:      CALL F5TCB (IDO, L1, CC, CH, WA1, WA2)
C
C  Arguments:
C     IDO    - 2*N divided by present and previous factors.  (Input)
C     L1     - Product of previous factors used.  (Input)
C     CC     - SEQ or CH from F3TCB.  (Input)
C     CH     - A better estimate for COEF.  (Output)  Note:  The roles
C              of CC and CH are reversed when NA=1.
C     WA1    - Real array containing needed elements of array WA
C              in F3TCB.  (Input)
C     WA2    - Real array containing other needed elements of array
C              WA in F3TCB.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF5TCB (IDO, L1, CC, CH, WA1, WA2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, L1
      DOUBLE PRECISION CC(IDO,3,*), CH(IDO,L1,*), WA1(*), WA2(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, K
      DOUBLE PRECISION CI2, CI3, CR2, CR3, DI2, DI3, DR2, DR3, TI2, TR2
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      DOUBLE PRECISION TAUI, TAUR
      SAVE       TAUI, TAUR
C
      DATA TAUR/-0.5D0/
      DATA TAUI/0.866025403784438646763723170753D0/
C
      IF (IDO .EQ. 2) THEN
         DO 10  K=1, L1
            TR2 = CC(1,2,K) + CC(1,3,K)
            CR2 = CC(1,1,K) + TAUR*TR2
            CH(1,K,1) = CC(1,1,K) + TR2
            TI2 = CC(2,2,K) + CC(2,3,K)
            CI2 = CC(2,1,K) + TAUR*TI2
            CH(2,K,1) = CC(2,1,K) + TI2
            CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
            CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
            CH(1,K,2) = CR2 - CI3
            CH(1,K,3) = CR2 + CI3
            CH(2,K,2) = CI2 + CR3
            CH(2,K,3) = CI2 - CR3
   10    CONTINUE
      ELSE
         DO 30  K=1, L1
            DO 20  I=2, IDO, 2
               TR2 = CC(I-1,2,K) + CC(I-1,3,K)
               CR2 = CC(I-1,1,K) + TAUR*TR2
               CH(I-1,K,1) = CC(I-1,1,K) + TR2
               TI2 = CC(I,2,K) + CC(I,3,K)
               CI2 = CC(I,1,K) + TAUR*TI2
               CH(I,K,1) = CC(I,1,K) + TI2
               CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
               CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
               DR2 = CR2 - CI3
               DR3 = CR2 + CI3
               DI2 = CI2 + CR3
               DI3 = CI2 - CR3
               CH(I,K,2) = WA1(I-1)*DI2 + WA1(I)*DR2
               CH(I-1,K,2) = WA1(I-1)*DR2 - WA1(I)*DI2
               CH(I,K,3) = WA2(I-1)*DI3 + WA2(I)*DR3
               CH(I-1,K,3) = WA2(I-1)*DR3 - WA2(I)*DI3
   20       CONTINUE
   30    CONTINUE
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F4TCB/DF4TCB (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the complex periodic sequence from its Fourier
C              coefficients.
C
C  Usage:      CALL F4TCB (IDO, L1, CC, CH, WA1)
C
C  Arguments:
C     IDO    - 2*N divided by present and previous factors.  (Input)
C     L1     - Product of previous factors used.  (Input)
C     CC     - SEQ or CH from F3TCB.  (Input)
C     CH     - A better estimate for COEF.  (Output)  Note:  The roles
C              of CC and CH are reversed when NA=1.
C     WA1    - Real array containing the needed elements of array WA
C              in F3TCB.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF4TCB (IDO, L1, CC, CH, WA1)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, L1
      DOUBLE PRECISION CC(IDO,2,*), CH(IDO,L1,*), WA1(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, K
      DOUBLE PRECISION TI2, TR2
C
      IF (IDO .LE. 2) THEN
         DO 10  K=1, L1
            CH(1,K,1) = CC(1,1,K) + CC(1,2,K)
            CH(1,K,2) = CC(1,1,K) - CC(1,2,K)
            CH(2,K,1) = CC(2,1,K) + CC(2,2,K)
            CH(2,K,2) = CC(2,1,K) - CC(2,2,K)
   10    CONTINUE
      ELSE
         DO 30  K=1, L1
            DO 20  I=2, IDO, 2
               CH(I-1,K,1) = CC(I-1,1,K) + CC(I-1,2,K)
               TR2 = CC(I-1,1,K) - CC(I-1,2,K)
               CH(I,K,1) = CC(I,1,K) + CC(I,2,K)
               TI2 = CC(I,1,K) - CC(I,2,K)
               CH(I,K,2) = WA1(I-1)*TI2 + WA1(I)*TR2
               CH(I-1,K,2) = WA1(I-1)*TR2 - WA1(I)*TI2
   20       CONTINUE
   30    CONTINUE
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F8TCB/DF8TCB (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the complex periodic sequence from its Fourier
C              coefficients.
C
C  Usage:      CALL F8TCB (NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2,
C                          WA)
C
C  Arguments:
C     NAC    - Integer set to 1 if IDO=2, set to 0 otherwise.  (Output)
C     IDO    - 2*N divided by present and previous factors.  (Input)
C     IP     - The factor now being used.  (Input)
C     L1     - Product of previous factors used.  (Input)
C     IDL1   - ID0*L1
C     CC     - Matrix containing C or CH from F3TCB.  (Input/Output)
C     C1     - Matrix containing C or CH from F3TCB.  (Input/Output)
C     C2     - Matrix containing C or CH from F3TCB.  (Input/Output)
C     CH     - Matrix that contains a better estimate for
C              COEF.  (Output)
C     CH2    - Matrix that contains a better estimate for
C              COEF.  (Output)  Note:  The roles of CC and CH are
C              reversed when NA=1.
C     WA     - Real array containing the needed elements of array WA
C              in F3TCB.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF8TCB (NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2,
     &                   WA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NAC, IDO, IP, L1, IDL1
      DOUBLE PRECISION CC(IDO,IP,*), C1(IDO,L1,*), C2(IDL1,*),
     &           CH(IDO,L1,*), CH2(IDL1,*), WA(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IDIJ, IDJ, IDL, IDLJ, IDOT, IDP, IK, INC, IPP2,
     &           IPPH, J, JC, K, L, LC, NT
      DOUBLE PRECISION WAI, WAR
C
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP + 2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 60
      DO 30  J=2, IPPH
         JC = IPP2 - J
         DO 20  K=1, L1
            DO 10  I=1, IDO
               CH(I,K,J) = CC(I,J,K) + CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K) - CC(I,JC,K)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      DO 50  K=1, L1
         DO 40  I=1, IDO
            CH(I,K,1) = CC(I,1,K)
   40    CONTINUE
   50 CONTINUE
      GO TO 120
   60 DO 90  J=2, IPPH
         JC = IPP2 - J
         DO 80  I=1, IDO
            DO 70  K=1, L1
               CH(I,K,J) = CC(I,J,K) + CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K) - CC(I,JC,K)
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
      DO 110  I=1, IDO
         DO 100  K=1, L1
            CH(I,K,1) = CC(I,1,K)
  100    CONTINUE
  110 CONTINUE
  120 IDL = 2 - IDO
      INC = 0
      DO 160  L=2, IPPH
         LC = IPP2 - L
         IDL = IDL + IDO
         DO 130  IK=1, IDL1
            C2(IK,L) = CH2(IK,1) + WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  130    CONTINUE
         IDLJ = IDL
         INC = INC + IDO
         DO 150  J=3, IPPH
            JC = IPP2 - J
            IDLJ = IDLJ + INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ - IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 140  IK=1, IDL1
               C2(IK,L) = C2(IK,L) + WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC) + WAI*CH2(IK,JC)
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
      DO 180  J=2, IPPH
         DO 170  IK=1, IDL1
            CH2(IK,1) = CH2(IK,1) + CH2(IK,J)
  170    CONTINUE
  180 CONTINUE
      DO 200  J=2, IPPH
         JC = IPP2 - J
         DO 190  IK=2, IDL1, 2
            CH2(IK-1,J) = C2(IK-1,J) - C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J) + C2(IK,JC)
            CH2(IK,J) = C2(IK,J) + C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J) - C2(IK-1,JC)
  190    CONTINUE
  200 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) GO TO 310
      NAC = 0
      DO 210  IK=1, IDL1
         C2(IK,1) = CH2(IK,1)
  210 CONTINUE
      DO 230  J=2, IP
         DO 220  K=1, L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  220    CONTINUE
  230 CONTINUE
      IF (IDOT .GT. L1) GO TO 270
      IDIJ = 0
      DO 260  J=2, IP
         IDIJ = IDIJ + 2
         DO 250  I=4, IDO, 2
            IDIJ = IDIJ + 2
            DO 240  K=1, L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J) -
     &                       WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J) + WA(IDIJ)*CH(I-1,K,J)
  240       CONTINUE
  250    CONTINUE
  260 CONTINUE
      GO TO 310
  270 IDJ = 2 - IDO
      DO 300  J=2, IP
         IDJ = IDJ + IDO
         DO 290  K=1, L1
            IDIJ = IDJ
            DO 280  I=4, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J) -
     &                       WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J) + WA(IDIJ)*CH(I-1,K,J)
  280       CONTINUE
  290    CONTINUE
  300 CONTINUE
  310 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F6TCB/DF6TCB (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the complex periodic sequence from its Fourier
C              coefficients.
C
C  Usage:      CALL F6TCB (IDO, L1, CC, CH, WA1, WA2, WA3)
C
C  Arguments:
C     IDO    - 2*N divided by present and previous factors.  (Input)
C     L1     - Product of previous factors used.  (Input)
C     CC     - SEQ or CH from F3TCB.  (Input)
C     CH     - A better estimate for COEF.  (Output)  Note:  The roles
C              of CC and CH are reversed when NA=1.
C     WA1    - Real array containing the needed elements of array WA
C              in F3TCB.  (Input)
C     WA2    - Real array containing other needed elements of array
C              WA in F3TCB.  (Input)
C     WA3    - Real array containing other needed elements of array
C              WA in F3TCB.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF6TCB (IDO, L1, CC, CH, WA1, WA2, WA3)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, L1
      DOUBLE PRECISION CC(IDO,4,*), CH(IDO,L1,*), WA1(*), WA2(*),
     &           WA3(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, K
      DOUBLE PRECISION CI2, CI3, CI4, CR2, CR3, CR4, TI1, TI2, TI3,
     &           TI4, TR1, TR2, TR3, TR4
C
      IF (IDO .EQ. 2) THEN
         DO 10  K=1, L1
            TI1 = CC(2,1,K) - CC(2,3,K)
            TI2 = CC(2,1,K) + CC(2,3,K)
            TR4 = CC(2,4,K) - CC(2,2,K)
            TI3 = CC(2,2,K) + CC(2,4,K)
            TR1 = CC(1,1,K) - CC(1,3,K)
            TR2 = CC(1,1,K) + CC(1,3,K)
            TI4 = CC(1,2,K) - CC(1,4,K)
            TR3 = CC(1,2,K) + CC(1,4,K)
            CH(1,K,1) = TR2 + TR3
            CH(1,K,3) = TR2 - TR3
            CH(2,K,1) = TI2 + TI3
            CH(2,K,3) = TI2 - TI3
            CH(1,K,2) = TR1 + TR4
            CH(1,K,4) = TR1 - TR4
            CH(2,K,2) = TI1 + TI4
            CH(2,K,4) = TI1 - TI4
   10    CONTINUE
      ELSE
         DO 30  K=1, L1
            DO 20  I=2, IDO, 2
               TI1 = CC(I,1,K) - CC(I,3,K)
               TI2 = CC(I,1,K) + CC(I,3,K)
               TI3 = CC(I,2,K) + CC(I,4,K)
               TR4 = CC(I,4,K) - CC(I,2,K)
               TR1 = CC(I-1,1,K) - CC(I-1,3,K)
               TR2 = CC(I-1,1,K) + CC(I-1,3,K)
               TI4 = CC(I-1,2,K) - CC(I-1,4,K)
               TR3 = CC(I-1,2,K) + CC(I-1,4,K)
               CH(I-1,K,1) = TR2 + TR3
               CR3 = TR2 - TR3
               CH(I,K,1) = TI2 + TI3
               CI3 = TI2 - TI3
               CR2 = TR1 + TR4
               CR4 = TR1 - TR4
               CI2 = TI1 + TI4
               CI4 = TI1 - TI4
               CH(I-1,K,2) = WA1(I-1)*CR2 - WA1(I)*CI2
               CH(I,K,2) = WA1(I-1)*CI2 + WA1(I)*CR2
               CH(I-1,K,3) = WA2(I-1)*CR3 - WA2(I)*CI3
               CH(I,K,3) = WA2(I-1)*CI3 + WA2(I)*CR3
               CH(I-1,K,4) = WA3(I-1)*CR4 - WA3(I)*CI4
               CH(I,K,4) = WA3(I-1)*CI4 + WA3(I)*CR4
   20       CONTINUE
   30    CONTINUE
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F2TRF/DF2TRF (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the Fourier coefficients of a real periodic
C              sequence.
C
C  Usage:      CALL F2TRF (N, SEQ, COEF, WFFTR)
C
C  Arguments:  (See FFTRF)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF2TRF (N, SEQ, COEF, WFFTR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
      DOUBLE PRECISION SEQ(*), COEF(*), WFFTR(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DCOPY, DF3TRF
C
C                                  CHECK ARGUMENT N
      IF (N .LT. 1) THEN
         CALL E1PSH ('DF2TRF ')
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'The length of the sequence N = %(I1).  '//
     &               'It must be at least 1.')
         CALL E1POP ('DF2TRF ')
         GO TO 9000
      END IF
C                                  COPY SEQ TO COEF
      CALL DCOPY (N, SEQ, 1, COEF, 1)
C
      IF (N .GT. 1) THEN
         CALL DF3TRF (N, COEF, WFFTR, WFFTR(N+1), WFFTR(2*N+1))
      END IF
C
 9000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  FFTRI/DFFTRI (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute parameters needed by FFTRF and FFTRB.
C
C  Usage:      CALL FFTRI (N, WFFTR)
C
C  Arguments:
C     N      - Length of the sequence to be transformed.  (Input)
C     WFFTR  - Array of length 2*N+15 containing parameters needed by
C              FFTRF and FFTRB.  (Output)
C
C  Remark:
C     Different WFFTR arrays are needed for different values of N.
C
C  GAMS:       J1a1
C
C  Chapters:   MATH/LIBRARY Transforms
C              STAT/LIBRARY Mathematical Support
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DFFTRI (N, WFFTR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
      DOUBLE PRECISION WFFTR(*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DF3TRI
C
C                                  CHECK ARGUMENT N
      IF (N .LT. 1) THEN
         CALL E1PSH ('DFFTRI ')
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'The length of the sequence N = %(I1).  '//
     &               'It must be at least 1.')
         CALL E1POP ('DFFTRI ')
         GO TO 9000
      END IF
C
      IF (N .GT. 1) THEN
         CALL DF3TRI (N, WFFTR(N+1), WFFTR(2*N+1))
      END IF
C
 9000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F3TRF/DF3TRF (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the Fourier coefficients of a real periodic
C              sequence.
C
C  Usage:      CALL F3TRF (N, C, CH, WA, FAC)
C
C  Arguments:
C     N      - Length of the sequence to be transformed.  (Input)
C     C      - Real array.  (Input/Output)
C     CH     - Real array.  (Workspace)
C     WA     - Real array from FFTRI.  (Input)
C     FAC    - Real array from FFTRI.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF3TRF (N, C, CH, WA, FAC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
      DOUBLE PRECISION C(*), CH(*), WA(*), FAC(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IDL1, IDO, IP, IW, IX2, IX3, IX4, K1, KH, L1, L2, NA,
     &           NF
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  NINT
      INTRINSIC  NINT
      INTEGER    NINT
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   DCOPY, DF4TRF, DF5TRF, DF6TRF, DF7TRF, DF8TRF
C
      IF (N .EQ. 1) GO TO 30
      NF = NINT(FAC(2))
      NA = 1
      L2 = N
      IW = N
      DO 20  K1=1, NF
         KH = NF - K1
         IP = NINT(FAC(KH+3))
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW - (IP-1)*IDO
         NA = 1 - NA
         IF (IP .EQ. 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA .NE. 0) THEN
               CALL DF6TRF (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL DF6TRF (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            END IF
            GO TO 10
         END IF
         IF (IP .EQ. 2) THEN
            IF (NA .NE. 0) THEN
               CALL DF4TRF (IDO, L1, CH, C, WA(IW))
            ELSE
               CALL DF4TRF (IDO, L1, C, CH, WA(IW))
            END IF
            GO TO 10
         END IF
         IF (IP .EQ. 3) THEN
            IX2 = IW + IDO
            IF (NA .NE. 0) THEN
               CALL DF5TRF (IDO, L1, CH, C, WA(IW), WA(IX2))
            ELSE
               CALL DF5TRF (IDO, L1, C, CH, WA(IW), WA(IX2))
            END IF
            GO TO 10
         END IF
         IF (IP .EQ. 5) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IX4 = IX3 + IDO
            IF (NA .NE. 0) THEN
               CALL DF7TRF (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3),
     &                      WA(IX4))
            ELSE
               CALL DF7TRF (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3),
     &                      WA(IX4))
            END IF
            GO TO 10
         END IF
         IF (IDO .EQ. 1) NA = 1 - NA
         IF (NA .NE. 0) THEN
            CALL DF8TRF (IDO, IP, L1, IDL1, CH, CH, CH, C, C, WA(IW))
            NA = 0
         ELSE
            CALL DF8TRF (IDO, IP, L1, IDL1, C, C, C, CH, CH, WA(IW))
            NA = 1
         END IF
   10    L2 = L1
   20 CONTINUE
      IF (NA .NE. 1) CALL DCOPY (N, CH, 1, C, 1)
   30 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F8TRF/DF8TRF (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the Fourier coefficients of a real periodic
C              sequence.
C
C  Usage:      CALL F8TRF (IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
C
C  Arguments:
C     IDO    - Real scalar.  (Input)
C     IP     - Real scalar.  (Input)
C     L1     - Real scalar.  (Input)
C     IDL1   - Real scalar.  (Input)
C     CC     - Real array.  (Input/Output)
C     C1     - Real array.  (Input/Output)
C     C2     - Real array.  (Input/Output)
C     CH     - Real array.  (Output)
C     CH2    - Real array.  (Output)
C     WA     - Real vector.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF8TRF (IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, IP, L1, IDL1
      DOUBLE PRECISION CC(IDO,IP,*), C1(IDO,L1,*), C2(IDL1,*),
     &           CH(IDO,L1,*), CH2(IDL1,*), WA(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IC, IDIJ, IDP2, IK, IPP2, IPPH, IS, J, J2, JC, K,
     &           L, LC, NBD
      DOUBLE PRECISION AI1, AI2, AR1, AR1H, AR2, AR2H, ARG, DC2, DCP,
     &           DS2, DSP
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DCOS,DBLE,DSIN
      INTRINSIC  DCOS, DBLE, DSIN
      DOUBLE PRECISION DCOS, DBLE, DSIN
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   DCOPY
C
      ARG = 2.0D0*3.1415926535897932384626433831D0/DBLE(IP)
      DCP = DCOS(ARG)
      DSP = DSIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP + 2
      IDP2 = IDO + 2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 180
      CALL DCOPY (IDL1, C2, 1, CH2, 1)
      DO 20  J=2, IP
         DO 10  K=1, L1
            CH(1,K,J) = C1(1,K,J)
   10    CONTINUE
   20 CONTINUE
      IF (NBD .GT. L1) GO TO 60
      IS = -IDO
      DO 50  J=2, IP
         IS = IS + IDO
         IDIJ = IS
         DO 40  I=3, IDO, 2
            IDIJ = IDIJ + 2
            DO 30  K=1, L1
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J) +
     &                       WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J) - WA(IDIJ)*C1(I-1,K,J)
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      GO TO 100
   60 IS = -IDO
      DO 90  J=2, IP
         IS = IS + IDO
         DO 80  K=1, L1
            IDIJ = IS
            DO 70  I=3, IDO, 2
               IDIJ = IDIJ + 2
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J) +
     &                       WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J) - WA(IDIJ)*C1(I-1,K,J)
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
  100 IF (NBD .LT. L1) GO TO 140
      DO 130  J=2, IPPH
         JC = IPP2 - J
         DO 120  K=1, L1
            DO 110  I=3, IDO, 2
               C1(I-1,K,J) = CH(I-1,K,J) + CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J) - CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J) + CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC) - CH(I-1,K,J)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
      GO TO 190
  140 DO 170  J=2, IPPH
         JC = IPP2 - J
         DO 160  I=3, IDO, 2
            DO 150  K=1, L1
               C1(I-1,K,J) = CH(I-1,K,J) + CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J) - CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J) + CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC) - CH(I-1,K,J)
  150       CONTINUE
  160    CONTINUE
  170 CONTINUE
      GO TO 190
  180 CALL DCOPY (IDL1, CH2, 1, C2, 1)
  190 DO 210  J=2, IPPH
         JC = IPP2 - J
         DO 200  K=1, L1
            C1(1,K,J) = CH(1,K,J) + CH(1,K,JC)
            C1(1,K,JC) = CH(1,K,JC) - CH(1,K,J)
  200    CONTINUE
  210 CONTINUE
C
      AR1 = 1.0D0
      AI1 = 0.0D0
      DO 250  L=2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         DO 220  IK=1, IDL1
            CH2(IK,L) = C2(IK,1) + AR1*C2(IK,2)
            CH2(IK,LC) = AI1*C2(IK,IP)
  220    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 240  J=3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            DO 230  IK=1, IDL1
               CH2(IK,L) = CH2(IK,L) + AR2*C2(IK,J)
               CH2(IK,LC) = CH2(IK,LC) + AI2*C2(IK,JC)
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      DO 270  J=2, IPPH
         DO 260  IK=1, IDL1
            CH2(IK,1) = CH2(IK,1) + C2(IK,J)
  260    CONTINUE
  270 CONTINUE
C
      IF (IDO .LT. L1) GO TO 300
      DO 290  K=1, L1
         DO 280  I=1, IDO
            CC(I,1,K) = CH(I,K,1)
  280    CONTINUE
  290 CONTINUE
      GO TO 330
  300 DO 320  I=1, IDO
         DO 310  K=1, L1
            CC(I,1,K) = CH(I,K,1)
  310    CONTINUE
  320 CONTINUE
  330 DO 350  J=2, IPPH
         JC = IPP2 - J
         J2 = J + J
         DO 340  K=1, L1
            CC(IDO,J2-2,K) = CH(1,K,J)
            CC(1,J2-1,K) = CH(1,K,JC)
  340    CONTINUE
  350 CONTINUE
      IF (IDO .EQ. 1) GO TO 430
      IF (NBD .LT. L1) GO TO 390
      DO 380  J=2, IPPH
         JC = IPP2 - J
         J2 = J + J
         DO 370  K=1, L1
            DO 360  I=3, IDO, 2
               IC = IDP2 - I
               CC(I-1,J2-1,K) = CH(I-1,K,J) + CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J) - CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J) + CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC) - CH(I,K,J)
  360       CONTINUE
  370    CONTINUE
  380 CONTINUE
      GO TO 430
  390 DO 420  J=2, IPPH
         JC = IPP2 - J
         J2 = J + J
         DO 410  I=3, IDO, 2
            IC = IDP2 - I
            DO 400  K=1, L1
               CC(I-1,J2-1,K) = CH(I-1,K,J) + CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J) - CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J) + CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC) - CH(I,K,J)
  400       CONTINUE
  410    CONTINUE
  420 CONTINUE
  430 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F7TRF/DF7TRF (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the Fourier coefficients of a real periodic
C              sequence.
C
C  Usage:      CALL F7TRF (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C
C  Arguments:
C     IDO    - Real scalar.  (Input)
C     L1     - Real scalar.  (Input)
C     CC     - Real array.  (Input/Output)
C     CH     - Real array.  (Input/Output)
C     WA1    - Real vector.  (Input)
C     WA2    - Real vector.  (Input)
C     WA3    - Real vector.  (Input)
C     WA4    - Real vector.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF7TRF (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, L1
      DOUBLE PRECISION CC(IDO,L1,*), CH(IDO,5,*), WA1(*), WA2(*),
     &           WA3(*), WA4(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IC, IDP2, K
      DOUBLE PRECISION CI2, CI3, CI4, CI5, CR2, CR3, CR4, CR5, DI2,
     &           DI3, DI4, DI5, DR2, DR3, DR4, DR5, TI2, TI3, TI4,
     &           TI5, TR2, TR3, TR4, TR5
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      DOUBLE PRECISION TI11, TI12, TR11, TR12
      SAVE       TI11, TI12, TR11, TR12
C
      DATA TR11/0.309016994374947424102293417183D0/
      DATA TI11/0.951056516295153572116439333379D0/
      DATA TR12/-0.809016994374947424102293417183D0/
      DATA TI12/0.587785252292473129168705954639D0/
C
      DO 10  K=1, L1
         CR2 = CC(1,K,5) + CC(1,K,2)
         CI5 = CC(1,K,5) - CC(1,K,2)
         CR3 = CC(1,K,4) + CC(1,K,3)
         CI4 = CC(1,K,4) - CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2 + CR3
         CH(IDO,2,K) = CC(1,K,1) + TR11*CR2 + TR12*CR3
         CH(1,3,K) = TI11*CI5 + TI12*CI4
         CH(IDO,4,K) = CC(1,K,1) + TR12*CR2 + TR11*CR3
         CH(1,5,K) = TI12*CI5 - TI11*CI4
   10 CONTINUE
      IF (IDO .EQ. 1) GO TO 40
      IDP2 = IDO + 2
      DO 30  K=1, L1
         DO 20  I=3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5) + WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5) - WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2 + DR5
            CI5 = DR5 - DR2
            CR5 = DI2 - DI5
            CI2 = DI2 + DI5
            CR3 = DR3 + DR4
            CI4 = DR4 - DR3
            CR4 = DI3 - DI4
            CI3 = DI3 + DI4
            CH(I-1,1,K) = CC(I-1,K,1) + CR2 + CR3
            CH(I,1,K) = CC(I,K,1) + CI2 + CI3
            TR2 = CC(I-1,K,1) + TR11*CR2 + TR12*CR3
            TI2 = CC(I,K,1) + TR11*CI2 + TR12*CI3
            TR3 = CC(I-1,K,1) + TR12*CR2 + TR11*CR3
            TI3 = CC(I,K,1) + TR12*CI2 + TR11*CI3
            TR5 = TI11*CR5 + TI12*CR4
            TI5 = TI11*CI5 + TI12*CI4
            TR4 = TI12*CR5 - TI11*CR4
            TI4 = TI12*CI5 - TI11*CI4
            CH(I-1,3,K) = TR2 + TR5
            CH(IC-1,2,K) = TR2 - TR5
            CH(I,3,K) = TI2 + TI5
            CH(IC,2,K) = TI5 - TI2
            CH(I-1,5,K) = TR3 + TR4
            CH(IC-1,4,K) = TR3 - TR4
            CH(I,5,K) = TI3 + TI4
            CH(IC,4,K) = TI4 - TI3
   20    CONTINUE
   30 CONTINUE
   40 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F5TRF/DF5TRF (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the Fourier coefficients of a real periodic
C              sequence.
C
C  Usage:      CALL F5TRF (IDO, L1, CC, CH, WA1, WA2)
C
C  Arguments:
C     IDO    - Real scalar.  (Input)
C     L1     - Real scalar.  (Input)
C     CC     - Real array.  (Input/Output)
C     CH     - Real array.  (Input/Output)
C     WA1    - Real vector.  (Input)
C     WA2    - Real vector.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF5TRF (IDO, L1, CC, CH, WA1, WA2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, L1
      DOUBLE PRECISION CC(IDO,L1,*), CH(IDO,3,*), WA1(*), WA2(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IC, IDP2, K
      DOUBLE PRECISION CI2, CR2, DI2, DI3, DR2, DR3, TI2, TI3, TR2, TR3
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      DOUBLE PRECISION TAUI, TAUR
      SAVE       TAUI, TAUR
C
      DATA TAUR/-0.5D0/
      DATA TAUI/0.866025403784438646763723170753D0/
C
      DO 10  K=1, L1
         CR2 = CC(1,K,2) + CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1) + TAUR*CR2
   10 CONTINUE
      IF (IDO .EQ. 1) GO TO 40
      IDP2 = IDO + 2
      DO 30  K=1, L1
         DO 20  I=3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2 + DR3
            CI2 = DI2 + DI3
            CH(I-1,1,K) = CC(I-1,K,1) + CR2
            CH(I,1,K) = CC(I,K,1) + CI2
            TR2 = CC(I-1,K,1) + TAUR*CR2
            TI2 = CC(I,K,1) + TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2 + TR3
            CH(IC-1,2,K) = TR2 - TR3
            CH(I,3,K) = TI2 + TI3
            CH(IC,2,K) = TI3 - TI2
   20    CONTINUE
   30 CONTINUE
   40 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F6TRF/DF6TRF (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the Fourier coefficients of a real periodic
C              sequence.
C
C  Usage:      CALL F6TRF (IDO, L1, CC, CH, WA1, WA2, WA3)
C
C  Arguments:
C     IDO    - Real scalar.  (Input)
C     L1     - Real scalar.  (Input)
C     CC     - Real array.  (Input/Output)
C     CH     - Real array.  (Input/Output)
C     WA1    - Real vector.  (Input)
C     WA2    - Real vector.  (Input)
C     WA3    - Real vector.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF6TRF (IDO, L1, CC, CH, WA1, WA2, WA3)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, L1
      DOUBLE PRECISION CC(IDO,L1,*), CH(IDO,4,*), WA1(*), WA2(*),
     &           WA3(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IC, IDP2, K
      DOUBLE PRECISION CI2, CI3, CI4, CR2, CR3, CR4, TI1, TI2, TI3,
     &           TI4, TR1, TR2, TR3, TR4
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      DOUBLE PRECISION HSQT2
      SAVE       HSQT2
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  MOD
      INTRINSIC  MOD
      INTEGER    MOD
C
      DATA HSQT2/0.707106781186547524400844362105D0/
C
      DO 10  K=1, L1
         TR1 = CC(1,K,2) + CC(1,K,4)
         TR2 = CC(1,K,1) + CC(1,K,3)
         CH(1,1,K) = TR1 + TR2
         CH(IDO,4,K) = TR2 - TR1
         CH(IDO,2,K) = CC(1,K,1) - CC(1,K,3)
         CH(1,3,K) = CC(1,K,4) - CC(1,K,2)
   10 CONTINUE
      IF (IDO-2) 70,50,20
   20 IDP2 = IDO + 2
      DO 40  K=1, L1
         DO 30  I=3, IDO, 2
            IC = IDP2 - I
            CR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2 + CR4
            TR4 = CR4 - CR2
            TI1 = CI2 + CI4
            TI4 = CI2 - CI4
            TI2 = CC(I,K,1) + CI3
            TI3 = CC(I,K,1) - CI3
            TR2 = CC(I-1,K,1) + CR3
            TR3 = CC(I-1,K,1) - CR3
            CH(I-1,1,K) = TR1 + TR2
            CH(IC-1,4,K) = TR2 - TR1
            CH(I,1,K) = TI1 + TI2
            CH(IC,4,K) = TI1 - TI2
            CH(I-1,3,K) = TI4 + TR3
            CH(IC-1,2,K) = TR3 - TI4
            CH(I,3,K) = TR4 + TI3
            CH(IC,2,K) = TR4 - TI3
   30    CONTINUE
   40 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) GO TO 70
   50 CONTINUE
      DO 60  K=1, L1
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
         CH(IDO,1,K) = TR1 + CC(IDO,K,1)
         CH(IDO,3,K) = CC(IDO,K,1) - TR1
         CH(1,2,K) = TI1 - CC(IDO,K,3)
         CH(1,4,K) = TI1 + CC(IDO,K,3)
   60 CONTINUE
   70 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F3TRI/DF3TRI (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute parameters needed by RFFTF and RFFTB.
C
C  Usage:      CALL F3TRI (N, WA, FAC)
C
C  Arguments:
C     N      - Length of the sequence to be transformed.  (Input)
C     WA     - Real vector.  (Output)
C     FAC    - Real vector.  (Output)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF3TRI (N, WA, FAC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N
      DOUBLE PRECISION WA(*), FAC(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IB, IDO, II, IP, IPM, IS, J, K1, L1, L2, LD, NF,
     &           NFM1, NL, NQ, NR, NTRY
      DOUBLE PRECISION ARG, ARGH, ARGLD, FI, TPI
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    NTRYH(4)
      SAVE       NTRYH
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  DCOS,DBLE,DSIN
      INTRINSIC  DCOS, DBLE, DSIN
      DOUBLE PRECISION DCOS, DBLE, DSIN
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   DCOPY
C
      DATA NTRYH/4, 2, 3, 5/
C
      NL = N
      NF = 0
      J = 0
   10 J = J + 1
      IF (J-4) 20,20,30
   20 NTRY = NTRYH(J)
      GO TO 40
   30 NTRY = NTRY + 2
   40 NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR) 10,50,10
   50 NF = NF + 1
      FAC(NF+2) = DBLE(NTRY)
      NL = NQ
      IF (NTRY .NE. 2) GO TO 60
      IF (NF .EQ. 1) GO TO 60
      CALL DCOPY (NF-1, FAC(3), -1, FAC(4), -1)
      FAC(3) = 2.0D0
   60 IF (NL .NE. 1) GO TO 40
      FAC(1) = DBLE(N)
      FAC(2) = DBLE(NF)
      TPI = 2.0D0*3.1415926535897932384626433831D0
      ARGH = TPI/DBLE(N)
      IS = 0
      NFM1 = NF - 1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 90  K1=1, NFM1
         IP = NINT(FAC(K1+2))
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP - 1
         DO 80  J=1, IPM
            LD = LD + L1
            I = IS
            ARGLD = DBLE(LD)*ARGH
            FI = 0.0D0
            DO 70  II=3, IDO, 2
               I = I + 2
               FI = FI + 1.0D0
               ARG = FI*ARGLD
               WA(I-1) = DCOS(ARG)
               WA(I) = DSIN(ARG)
   70       CONTINUE
            IS = IS + IDO
   80    CONTINUE
         L1 = L2
   90 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  F4TRF/DF4TRF (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    January 1, 1985
C
C  Purpose:    Compute the Fourier coefficients of a real periodic
C              sequence.
C
C  Usage:      CALL F4TRF (IDO, L1, CC, CH, WA1)
C
C  Arguments:
C     IDO    - Real scalar.  (Input)
C     L1     - Real scalar.  (Input)
C     CC     - Real array.  (Input/Output)
C     CH     - Real array.  (Input/Output)
C     WA1    - Real vector.  (Input)
C
C  Chapter:    MATH/LIBRARY Transforms
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DF4TRF (IDO, L1, CC, CH, WA1)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDO, L1
      DOUBLE PRECISION CC(IDO,L1,*), CH(IDO,2,*), WA1(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IC, IDP2, K
      DOUBLE PRECISION TI2, TR2
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  MOD
      INTRINSIC  MOD
      INTEGER    MOD
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   DCOPY, DVCAL
C
      DO 10  K=1, L1
         CH(1,1,K) = CC(1,K,1) + CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1) - CC(1,K,2)
   10 CONTINUE
      IF (IDO-2) 60,50,20
   20 IDP2 = IDO + 2
      DO 40  K=1, L1
         DO 30  I=3, IDO, 2
            IC = IDP2 - I
            TR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1) + TI2
            CH(IC,2,K) = TI2 - CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1) + TR2
            CH(IC-1,2,K) = CC(I-1,K,1) - TR2
   30    CONTINUE
   40 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) GO TO 60
   50 CALL DVCAL (L1, -1.0D0, CC(IDO,1,2), IDO, CH(1,2,1), 2*IDO)
      CALL DCOPY (L1, CC(IDO,1,1), IDO, CH(IDO,1,1), 2*IDO)
   60 RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  DVCAL (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Multiply a vector by a scalar and store the result in
C              another vector, y = ax, all double precision.
C
C  Usage:      CALL DVCAL (N, DA, DX, INCX, DY, INCY)
C
C  Arguments:
C     N      - Length of vectors X.  (Input)
C     DA     - Double precision scalar.  (Input)
C     DX     - Double precision vector of length MAX(N*IABS(INCX),1).
C                 (Input)
C              DVCAL computes DA*X(I) for I = 1,...,N. X(I) refers
C              to a specific element of DX.
C     INCX   - Displacement between elements of DX.  (Input)
C              X(I) is defined to be DX(1+(I-1)*INCX). INCX must be
C              greater than 0.
C     DY     - Double precision vector of length MAX(N*IABS(INCY),1).
C                 (Output)
C              DVCAL sets Y(I) equal to DA*X(I) for I = 1,...,N.
C              Y(I) refers to a specific element of DY.
C     INCY   - Displacement between elements of DY.  (Input)
C              Y(I) is defined to be DY(1+(I-1)*INCY). INCY must be
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
      SUBROUTINE DVCAL (N, DA, DX, INCX, DY, INCY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX, INCY
      DOUBLE PRECISION DA, DX(*), DY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX, IY, M, MP1
C                                  SPECIFICATIONS FOR SPECIAL CASES
C     INTRINSIC  MOD
      INTRINSIC  MOD
      INTEGER    MOD
C
      IF (N .GT. 0) THEN
         IF (INCX.NE.1 .OR. INCY.NE.1) THEN
C                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL
C                                  INCREMENTS NOT EQUAL TO 1
            IX = 1
            IY = 1
            DO 10  I=1, N
               DY(IY) = DA*DX(IX)
               IX = IX + INCX
               IY = IY + INCY
   10       CONTINUE
         ELSE
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
            M = MOD(N,4)
C                                  CLEAN-UP LOOP
            DO 30  I=1, M
               DY(I) = DA*DX(I)
   30       CONTINUE
            MP1 = M + 1
            DO 40  I=MP1, N, 4
               DY(I) = DA*DX(I)
               DY(I+1) = DA*DX(I+1)
               DY(I+2) = DA*DX(I+2)
               DY(I+3) = DA*DX(I+3)
   40       CONTINUE
         END IF
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  DCOPY (Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 9, 1986
C
C  Purpose:    Copy a vector X to a vector Y, both double precision.
C
C  Usage:      CALL DCOPY (N, DX, INCX, DY, INCY)
C
C  Arguments:
C     N      - Length of vectors X and Y.  (Input)
C     DX     - Double precision vector of length MAX(N*IABS(INCX),1).
C              (Input)
C     INCX   - Displacement between elements of DX.  (Input)
C              X(I) is defined to be.. DX(1+(I-1)*INCX) if INCX .GE. 0
C              or DX(1+(I-N)*INCX) if INCX .LT. 0.
C     DY     - Double precision vector of length MAX(N*IABS(INCY),1).
C              (Output)
C              DCOPY copies X(I) to Y(I) for I=1,...,N. X(I) and Y(I)
C              refer to specific elements of DX and DY, respectively.
C              See INCX and INCY argument descriptions.
C     INCY   - Displacement between elements of DY.  (Input)
C              Y(I) is defined to be.. DY(1+(I-1)*INCY) if INCY .GE. 0
C              or DY(1+(I-N)*INCY) if INCY .LT. 0.
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
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX, INCY
      DOUBLE PRECISION DX(*), DY(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IX, IY, M, MP1
C                                  SPECIFICATIONS FOR SPECIAL CASES
C     INTRINSIC  MOD
      INTRINSIC  MOD
      INTEGER    MOD
C
      IF (N .GT. 0) THEN
         IF (INCX.NE.1 .OR. INCY.NE.1) THEN
C                                  CODE FOR UNEQUAL INCREMENTS.
            IX = 1
            IY = 1
            IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
            IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
            DO 10  I=1, N
               DY(IY) = DX(IX)
               IX = IX + INCX
               IY = IY + INCY
   10       CONTINUE
         ELSE
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                  CLEAN-UP LOOP SO REMAINING VECTOR
C                                  LENGTH IS A MULTIPLE OF 7.
            M = MOD(N,7)
            DO 30  I=1, M
               DY(I) = DX(I)
   30       CONTINUE
            MP1 = M + 1
            DO 40  I=MP1, N, 7
               DY(I) = DX(I)
               DY(I+1) = DX(I+1)
               DY(I+2) = DX(I+2)
               DY(I+3) = DX(I+3)
               DY(I+4) = DX(I+4)
               DY(I+5) = DX(I+5)
               DY(I+6) = DX(I+6)
   40       CONTINUE
         END IF
      END IF
      RETURN
      END
