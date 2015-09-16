C-----------------------------------------------------------------------
C  IMSL Name:  I1KGT
C
C  Computer:   VAX/SINGLE
C
C  Revised:    January 17, 1984
C
C  Purpose:    Allocate numerical workspace.
C
C  Usage:      I1KGT(NELMTS,ITYPE)
C
C  Arguments:
C     NELMTS - Number of elements of data type ITYPE to be
C              allocated.  (Input)
C     ITYPE  - Data type of array to be allocated.  (Input)
C                 1 - logical
C                 2 - integer
C                 3 - real
C                 4 - double precision
C                 5 - complex
C                 6 - double complex
C     I1KGT  - Integer function.  (Output)  Returns the index of the
C              first element in the current allocation.
C
C  Remarks:
C  1. On return, the array will occupy
C     WKSP(I1KGT), WKSP(I1KGT+1), ..., WKSP(I1KGT+NELMTS-1) where
C     WKSP is an array of data type ITYPE equivalenced to RWKSP.
C
C  2. If I1KGT is negative, the absolute value of I1KGT is the
C     additional workspace needed for the current allocation.
C
C  3. The allocator reserves the first sixteen integer locations of
C     the stack for its own internal bookkeeping.  These are initialized
C     by the function IWKIN upon the first call to the allocation
C     package.
C
C  4. The use of the first ten integer locations is as follows:
C      WKSP( 1) - LOUT    The number of current allocations
C      WKSP( 2) - LNOW    The current active length of the stack
C      WKSP( 3) - LUSED   The maximum value of WKSP(2) achieved
C                         thus far
C      WKSP( 4) - LBND    The lower bound of permanent storage which
C                         is one numeric storage unit more than the
C                         maximum allowed length of the stack.
C      WKSP( 5) - LMAX    The maximum length of the storage array
C      WKSP( 6) - LALC    The total number of allocations handled by
C                         I1KGT
C      WKSP( 7) - LNEED   The number of numeric storage units by which
C                         the array size must be increased for all past
C                         allocations to succeed
C      WKSP( 8) - LBOOK   The number of numeric storage units used for
C                         bookkeeping
C      WKSP( 9) - LCHAR   The pointer to the portion of the permanent
C                         stack which contains the bookkeeping and
C                         pointers for the character workspace
C                         allocation.
C      WKSP(10) - LLCHAR  The length of the array beginning at LCHAR
C                         set aside for character workspace bookkeeping
C                         and pointers.
C                 NOTE -  If character workspace is not being used,
C                         LCHAR and LLCHAR can be ignored.
C  5. The next six integer locations contain values describing the
C     amount of storage allocated by the allocation system to the
C     various data types.
C      WKSP(11) - Numeric storage units allocated to LOGICAL
C      WKSP(12) - Numeric storage units allocated to INTEGER
C      WKSP(13) - Numeric storage units allocated to REAL
C      WKSP(14) - Numeric storage units allocated to DOUBLE PRECISION
C      WKSP(15) - Numeric storage units allocated to COMPLEX
C      WKSP(16) - Numeric storage units allocated to DOUBLE COMPLEX
C
C  Copyright:  1984 by IMSL, Inc. All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION I1KGT (NELMTS, ITYPE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NELMTS, ITYPE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IDUMAL, IGAP, ILEFT, IPA, IPA7, ISA, ISA7,
     &           ISIZE(6), JTYPE, LALC, LBND, LBOOK, LMAX, LNEED,
     &           LNEED1, LNOW, LOUT, LUSED
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      LOGICAL    FIRST
      SAVE       FIRST
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                              SPECIFICATIONS FOR COMMON /ERCOM8/
      INTEGER    PROLVL, XXLINE(10), XXPLEN(10), ICALOC(10), INALOC(10)
      COMMON     /ERCOM8/ PROLVL, XXLINE, XXPLEN, ICALOC, INALOC
      SAVE       /ERCOM8/
C                              SPECIFICATIONS FOR COMMON /ERCOM9/
      CHARACTER  XXPROC(10)*31
      COMMON     /ERCOM9/ XXPROC
      SAVE       /ERCOM9/
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
C     INTRINSIC  IABS,MAX0,MOD
      INTRINSIC  IABS, MAX0, MOD
      INTEGER    IABS, MAX0, MOD
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1POS, E1PSH, E1STI, IWKIN
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1KQU
      INTEGER    I1KQU
C
      DATA FIRST/.TRUE./
C
      CALL E1PSH ('I1KGT ')
C
      IF (FIRST) THEN
C                                  INITIALIZE WORKSPACE IF NEEDED
         FIRST = .FALSE.
         CALL IWKIN (0)
      END IF
C                                  NUMBER OF ELEMENTS LESS THAN 0
      IF (NELMTS .LT. 0) THEN
         CALL E1STI (1, NELMTS)
         CALL E1MES (5, 2, 'Number of elements is not positive.%/'//
     &               'NELMTS = %(I1).')
         CALL E1POP ('I1KGT ')
         GO TO 9000
      END IF
C                                  ILLEGAL DATA TYPE REQUESTED
      IF (ITYPE.EQ.0 .OR. IABS(ITYPE).GE.7) THEN
         CALL E1MES (5, 3, 'Illegal data type requested.')
         CALL E1POP ('I1KGT ')
         GO TO 9000
      END IF
C                                  BOOKKEEPING OVERWRITTEN
      IF (LNOW.LT.LBOOK .OR. LNOW.GT.LUSED .OR. LUSED.GT.LMAX .OR.
     &    LNOW.GE.LBND .OR. LOUT.GT.LALC) THEN
         CALL E1MES (5, 4, 'One or more of the first eight '//
     &               'bookkeeping locations in IWKSP have been '//
     &               'overwritten.')
         CALL E1POP ('I1KGT ')
         GO TO 9000
      END IF
C
      CALL E1POP ('I1KGT ')
C                                  DETERMINE NUMBER OF LOCATIONS STILL
C                                  AVAILABLE FOR DATA TYPE ITYPE
C                                  NOTE: I1KQU ALLOWS FOR 2 INTEGER
C                                        POINTERS WHICH MUST BE HANDLED
C                                        ARTIFICIALLY IF ILEFT = 0.
      ILEFT = I1KQU(IABS(ITYPE))
C
      IF (ITYPE .GT. 0) THEN
C                                  RELEASABLE STORAGE
         IF (ILEFT .GE. NELMTS) THEN
            I1KGT = (LNOW*ISIZE(2)-1)/ISIZE(ITYPE) + 2
            I = ((I1KGT-1+NELMTS)*ISIZE(ITYPE)-1)/ISIZE(2) + 3
C                                  IWKSP(I-1) CONTAINS THE DATA TYPE FOR
C                                  THIS ALLOCATION. IWKSP(I) CONTAINS
C                                  LNOW FOR THE PREVIOUS ALLOCATION.
            IWKSP(I-1) = ITYPE
            IWKSP(I) = LNOW
            LOUT = LOUT + 1
            LALC = LALC + 1
            LNOW = I
            LUSED = MAX0(LUSED,LNOW)
            LNEED = 0
         ELSE
C                                  RELEASABLE STORAGE WAS REQUESTED
C                                  BUT THE STACK WOULD OVERFLOW.
C                                  THEREFORE, ALLOCATE RELEASABLE
C                                  SPACE THROUGH THE END OF THE STACK
            IF (LNEED .EQ. 0) THEN
               IDUMAL = (LNOW*ISIZE(2)-1)/ISIZE(ITYPE) + 2
               I = ((IDUMAL-1+ILEFT)*ISIZE(ITYPE)-1)/ISIZE(2) + 3
C                                  ADVANCE COUNTERS AND STORE POINTERS
C                                  IF THERE IS ROOM TO DO SO
               IF (I .LT. LBND) THEN
C                                  IWKSP(I-1) CONTAINS THE DATA TYPE FOR
C                                  THIS ALLOCATION. IWKSP(I) CONTAINS
C                                  LNOW FOR THE PREVIOUS ALLOCATION.
                  IWKSP(I-1) = ITYPE
                  IWKSP(I) = LNOW
                  LOUT = LOUT + 1
                  LALC = LALC + 1
                  LNOW = I
                  LUSED = MAX0(LUSED,LNOW)
               END IF
            END IF
C                                  CALCULATE AMOUNT NEEDED TO ACCOMODATE
C                                  THIS ALLOCATION REQUEST
            LNEED1 = (NELMTS-ILEFT)*ISIZE(ITYPE)
            IF (ILEFT .EQ. 0) THEN
               IGAP = ISIZE(ITYPE) - MOD(LNOW+LNEED,ISIZE(ITYPE))
               IF (IGAP .EQ. ISIZE(ITYPE)) IGAP = 0
               LNEED1 = LNEED1 + 2*ISIZE(2) + IGAP
            END IF
C                                  MODIFY LNEED ACCORDING TO THE SIZE
C                                  OF THE BASE BEING USED (D.P. HERE)
            LNEED = LNEED + ((LNEED1+ISIZE(3)-1)/ISIZE(3))
C                                  SINCE CURRENT ALLOCATION IS ILLEGAL,
C                                  RETURN THE NEGATIVE OF THE ADDITIONAL
C                                  AMOUNT NEEDED TO MAKE IT LEGAL
            I1KGT = -LNEED
         END IF
      ELSE
C                                  PERMANENT STORAGE
         IF (ILEFT .GE. NELMTS) THEN
            JTYPE = -ITYPE
            I1KGT = (LBND*ISIZE(2)-1)/ISIZE(JTYPE) + 1 - NELMTS
            I = ((I1KGT-1)*ISIZE(JTYPE))/ISIZE(2) - 1
C                                  IWKSP(I) CONTAINS LBND FOR PREVIOUS
C                                  PERMANENT STORAGE ALLOCATION.
C                                  IWKSP(I+1) CONTAINS THE DATA TYPE FOR
C                                  THIS ALLOCATION.
            IWKSP(I) = LBND
            IWKSP(I+1) = JTYPE
            LALC = LALC + 1
            LBND = I
            LNEED = 0
         ELSE
C                                  PERMANENT STORAGE WAS REQUESTED
C                                  BUT THE STACK WOULD OVERFLOW,
C                                  THEREFORE, ALLOCATE RELEASABLE
C                                  SPACE THROUGH THE END OF THE STACK
            IF (LNEED .EQ. 0) THEN
               JTYPE = -ITYPE
               IDUMAL = (LNOW*ISIZE(2)-1)/ISIZE(JTYPE) + 2
               I = ((IDUMAL-1+ILEFT)*ISIZE(JTYPE)-1)/ISIZE(2) + 3
C                                  ADVANCE COUNTERS AND STORE POINTERS
C                                  IF THERE IS ROOM TO DO SO
               IF (I .LT. LBND) THEN
C                                  IWKSP(I-1) CONTAINS THE DATA TYPE FOR
C                                  THIS ALLOCATION. IWKSP(I) CONTAINS
C                                  LNOW FOR THE PREVIOUS ALLOCATION.
                  IWKSP(I-1) = JTYPE
                  IWKSP(I) = LNOW
                  LOUT = LOUT + 1
                  LALC = LALC + 1
                  LNOW = I
                  LUSED = MAX0(LUSED,LNOW)
               END IF
            END IF
C                                  CALCULATE AMOUNT NEEDED TO ACCOMODATE
C                                  THIS ALLOCATION REQUEST
            LNEED1 = (NELMTS-ILEFT)*ISIZE(-ITYPE)
            IF (ILEFT .EQ. 0) THEN
               IGAP = ISIZE(-ITYPE) - MOD(LNOW+LNEED,ISIZE(-ITYPE))
               IF (IGAP .EQ. ISIZE(-ITYPE)) IGAP = 0
               LNEED1 = LNEED1 + 2*ISIZE(2) + IGAP
            END IF
C                                  MODIFY LNEED ACCORDING TO THE SIZE
C                                  OF THE BASE BEING USED (D.P. HERE)
            LNEED = LNEED + ((LNEED1+ISIZE(3)-1)/ISIZE(3))
C                                  SINCE CURRENT ALLOCATION IS ILLEGAL,
C                                  RETURN THE NEGATIVE OF THE ADDITIONAL
C                                  AMOUNT NEEDED TO MAKE IT LEGAL
            I1KGT = -LNEED
         END IF
      END IF
C                                  STACK OVERFLOW - UNRECOVERABLE ERROR
 9000 IF (LNEED .GT. 0) THEN
         CALL E1POS (-5, IPA, ISA)
         CALL E1POS (5, 0, 0)
         CALL E1POS (-7, IPA7, ISA7)
         CALL E1POS (7, 0, 0)
         CALL E1PSH ('I1KGT ')
         CALL E1STI (1, LNEED+(LMAX/ISIZE(3)))
         IF (XXLINE(PROLVL).GE.1 .AND. XXLINE(PROLVL).LE.999) THEN
            CALL E1MES (7, 1, 'Insufficient workspace for current '//
     &                  'allocation(s).  Correct by inserting the '//
     &                  'following PROTRAN line: $OPTIONS;WORKSPACE=%'//
     &                  '(I1)')
         ELSE
            CALL E1MES (5, 5, 'Insufficient workspace for current '//
     &                  'allocation(s). Correct by calling IWKIN '//
     &                  'from main program with the three following '//
     &                  'statements:  (REGARDLESS OF PRECISION)%/'//
     &                  '      COMMON /WORKSP/  RWKSP%/      REAL '//
     &                  'RWKSP(%(I1))%/      CALL IWKIN(%(I1))')
         END IF
         CALL E1POP ('I1KGT ')
         CALL E1POS (5, IPA, ISA)
         CALL E1POS (7, IPA7, ISA7)
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  E1POS
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 2, 1984
C
C  Purpose:    Set or retrieve print and stop attributes.
C
C  Usage:      CALL E1POS(IERTYP,IPATT,ISATT)
C
C  Arguments:
C     IERTYP - Integer specifying the error type for which print and
C              stop attributes are to be set or retrieved.  (Input)  If
C              IERTYP is 0 then the settings apply to all error types.
C              If IERTYP is between 1 and 7, then the settings only
C              apply to that specified error type.  If IERTYP is
C              negative then the current print and stop attributes will
C              be returned in IPATT and ISATT.
C     IPATT  - If IERTYP is positive, IPATT is an integer specifying the
C              desired print attribute as follows: -1 means no change,
C              0 means NO, 1 means YES, and 2 means assign the default
C              setting.  (Input)  If IERTYP is negative, IPATT is
C              returned as 1 if print is YES or 0 if print is NO for
C              error type IABS(IERTYP).  (Output)
C     ISATT  - If IERTYP is positive, ISATT is an integer specifying the
C              desired stop attribute as follows: -1 means no change,
C              0 means NO, 1 means YES, and 2 means assign the default
C              setting.  (Input)  If IERTYP is negative, ISATT is
C              returned as 1 if print is YES or 0 if print is NO for
C              error type IABS(IERTYP).  (Output)
C
C  Copyright:  1984 by IMSL, Inc.  All rights reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE E1POS (IERTYP, IPATT, ISATT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IERTYP, IPATT, ISATT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IER
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    DEFLTP(7), DEFLTS(7), IFINIT
      SAVE       DEFLTP, DEFLTS, IFINIT
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
      EXTERNAL   E1INIT, E1MES, E1STI
C
      DATA IFINIT/0/
      DATA DEFLTP/0, 0, 1, 1, 1, 1, 1/, DEFLTS/0, 0, 0, 1, 1, 0, 1/
C                                  INITIALIZE ERROR TABLE IF NECESSARY
      IF (IFINIT .EQ. 0) THEN
         CALL E1INIT
         IFINIT = 1
      END IF
      IER = 0
      IF (IERTYP .GE. 0) THEN
         IF (IPATT.LT.-1 .OR. IPATT.GT.2) THEN
            CALL E1STI (1, IPATT)
            CALL E1MES (5, 1, 'Invalid value specified for print '//
     &                  'table attribute.  IPATT must be -1, 0, 1, '//
     &                  'or 2.  IPATT = %(I1)')
            IER = 1
         END IF
         IF (ISATT.LT.-1 .OR. ISATT.GT.2) THEN
            CALL E1STI (1, ISATT)
            CALL E1MES (5, 1, 'Invalid value specified for stop '//
     &                  'table attribute.  ISATT must be -1, 0, 1, '//
     &                  'or 2.  ISATT = %(I1)')
            IER = 1
         END IF
      END IF
      IF (IER .EQ. 0) THEN
         IF (IERTYP .EQ. 0) THEN
            IF (IPATT.EQ.0 .OR. IPATT.EQ.1) THEN
               DO 10  I=1, 7
   10          PRINTB(I) = IPATT
            ELSE IF (IPATT .EQ. 2) THEN
C                                  ASSIGN DEFAULT SETTINGS
               DO 20  I=1, 7
   20          PRINTB(I) = DEFLTP(I)
            END IF
            IF (ISATT.EQ.0 .OR. ISATT.EQ.1) THEN
               DO 30  I=1, 7
   30          STOPTB(I) = ISATT
            ELSE IF (ISATT .EQ. 2) THEN
C                                  ASSIGN DEFAULT SETTINGS
               DO 40  I=1, 7
   40          STOPTB(I) = DEFLTS(I)
            END IF
         ELSE IF (IERTYP.GE.1 .AND. IERTYP.LE.7) THEN
            IF (IPATT.EQ.0 .OR. IPATT.EQ.1) THEN
               PRINTB(IERTYP) = IPATT
            ELSE IF (IPATT .EQ. 2) THEN
C                                  ASSIGN DEFAULT SETTING
               PRINTB(IERTYP) = DEFLTP(IERTYP)
            END IF
            IF (ISATT.EQ.0 .OR. ISATT.EQ.1) THEN
               STOPTB(IERTYP) = ISATT
            ELSE IF (ISATT .EQ. 2) THEN
C                                  ASSIGN DEFAULT SETTING
               STOPTB(IERTYP) = DEFLTS(IERTYP)
            END IF
         ELSE IF (IERTYP.LE.-1 .AND. IERTYP.GE.-7) THEN
            I = IABS(IERTYP)
            IPATT = PRINTB(I)
            ISATT = STOPTB(I)
         END IF
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  I1KQU
C
C  Computer:   VAX/SINGLE
C
C  Revised:    January 17, 1984
C
C  Purpose:    Return number of elements of data type ITYPE that
C              remain to be allocated in one request.
C
C  Usage:      I1KQU(ITYPE)
C
C  Arguments:
C     ITYPE  - Type of storage to be checked (Input)
C                 1 - logical
C                 2 - integer
C                 3 - real
C                 4 - double precision
C                 5 - complex
C                 6 - double complex
C     I1KQU  - Integer function. (Output) Returns number of elements
C              of data type ITYPE remaining in the stack.
C
C  Copyright:  1983 by IMSL, Inc.  All Rights Reserved
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION I1KQU (ITYPE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    ITYPE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    ISIZE(6), LALC, LBND, LBOOK, LMAX, LNEED, LNOW, LOUT,
     &           LUSED
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
      EXTERNAL   E1MES, E1POP, E1PSH, IWKIN
C
      DATA FIRST/.TRUE./
C
      CALL E1PSH ('I1KQU ')
C
      IF (FIRST) THEN
C                                  INITIALIZE WORKSPACE IF NEEDED
         FIRST = .FALSE.
         CALL IWKIN (0)
      END IF
C                                  BOOKKEEPING OVERWRITTEN
      IF (LNOW.LT.LBOOK .OR. LNOW.GT.LUSED .OR. LUSED.GT.LMAX .OR.
     &    LNOW.GE.LBND .OR. LOUT.GT.LALC) THEN
         CALL E1MES (5, 7, 'One or more of the first eight '//
     &               'bookkeeping locations in IWKSP have been '//
     &               'overwritten.')
      ELSE IF (ITYPE.LE.0 .OR. ITYPE.GE.7) THEN
C                                  ILLEGAL DATA TYPE REQUESTED
         CALL E1MES (5, 8, 'Illegal data type requested.')
      ELSE
C                                  THIS CALCULATION ALLOWS FOR THE
C                                  TWO POINTER LOCATIONS IN THE STACK
C                                  WHICH ARE ASSIGNED TO EACH ALLOCATION
         I1KQU = MAX0(((LBND-3)*ISIZE(2))/ISIZE(ITYPE)-(LNOW*ISIZE(2)-
     &           1)/ISIZE(ITYPE)-1,0)
      END IF
C
      CALL E1POP ('I1KQU ')
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

C-----------------------------------------------------------------------
C  IMSL Name:  N1RCD
C
C  Computer:   VAX/SINGLE
C
C  Revised:    March 6, 1984
C
C  Purpose:    Retrieve an error code.
C
C  Usage:      N1RCD(IOPT)
C
C  Arguments:
C     IOPT   - Integer specifying the level.  (Input)
C              If IOPT=0 the error code for the current level is
C              returned.  If IOPT=1 the error code for the most
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
      INTEGER FUNCTION N1RCD (IOPT)
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
         CALL M1VECH ('.  The argument passed to N1RCD must be 0 or '//
     &                '1. ', MSGLEN, MSGSAV, MSGLEN)
         CALL E1PRT
         STOP
      ELSE
         N1RCD = ERCODE(CALLVL+IOPT)
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

C------------------------------------------------
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

C------------------------------------------
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
C-----------------------------------------------------------------------
C  IMSL Name:  PPDER/DPPDER (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 11, 1986
C
C  Purpose:    Evaluate the derivative of a piecewise polynomial.
C
C  Usage:      PPDER(IDERIV, X, KORDER, NINTV, BREAK, PPCOEF)
C
C  Arguments:
C     IDERIV - Order of the derivative to be evaluated.  (Input)
C              In particular, IDERIV = 0 returns the value of the
C              polynomial.
C     X      - Point at which the polynomial is to be evaluated.
C              (Input)
C     KORDER - Order of the polynomial.  (Input)
C     NINTV  - Number of polynomial pieces.  (Input)
C     BREAK  - Array of length NINTV+1 containing the breakpoints of
C              the piecewise polynomial representation.  (Input)
C              BREAK must be strictly increasing.
C     PPCOEF - Array of size KORDER*NINTV containing the
C              local coefficients of the piecewise polynomial pieces.
C              (Input)
C              PPCOEF is treated internally as a matrix of size
C              KORDER by NINTV.
C     PPDER  - Value of the IDERIV-th derivative of the piecewise
C              polynomial at X.  (Output)
C
C  GAMS:       E3; K6
C
C  Chapter:    MATH/LIBRARY Interpolation and Approximation
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DPPDER (IDERIV, X, KORDER, NINTV,
     &                 BREAK, PPCOEF)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDERIV, KORDER, NINTV
      DOUBLE PRECISION X, BREAK(*), PPCOEF(KORDER,*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    J, LEFT
      DOUBLE PRECISION FMM, H, VALUE
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DP3DER
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   N1RTY
      INTEGER    N1RTY
C
      CALL E1PSH ('DPPDER ')
C                                  IN CASE OF ERRORS
      VALUE = 0.0D0
C                                  Check argument NINTV
      IF (NINTV .LT. 1) THEN
         CALL E1STI (1, NINTV)
         CALL E1MES (5, 1, 'The number of intervals must be '//
     &               'at least 1 while NINTV = %(I1) is given. ')
      END IF
C                                  Check argument IDERIV
      IF (IDERIV .LT. 0) THEN
         CALL E1STI (1, IDERIV)
         CALL E1MES (5, 2, 'The order of the derivative must '//
     &               'be positive while IDERIV = %(I1) is given.')
      END IF
C                                  Check argument KORDER
      IF (KORDER .LE. 0) THEN
         CALL E1STI (1, KORDER)
         CALL E1MES (5, 3, 'The order of the interpolating '//
     &               'polynomial must be positive while KORDER = '//
     &               '%(I1) is given.')
      END IF
C                                  Check for errors
      IF (N1RTY(0) .NE. 0) GO TO 9000
C                                  Derivatives of order KORDER or
C                                  higher are identically zero
      IF (IDERIV .GE. KORDER) GO TO 9000
C                                  Find index I of largest breakpoint
C                                  to the left of X
      CALL DP3DER (KORDER, NINTV, BREAK, X, LEFT)
C                                  Evaluate jderiv-th derivative of
C                                  I-th polynomial piece at X
      FMM = KORDER - IDERIV
      H = X - BREAK(LEFT)
      DO 10  J=KORDER, IDERIV + 1, -1
         VALUE = (VALUE/FMM)*H + PPCOEF(J,LEFT)
         FMM = FMM - 1.0D0
   10 CONTINUE
C
 9000 CONTINUE
      DPPDER = VALUE
      CALL E1POP ('DPPDER ')
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  C2DEC/DC2DEC (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    October 1, 1984
C
C  Purpose:    Compute the cubic spline interpolant with specified
C              derivative endpoint conditions.
C
C  Usage:      CALL C2DEC (NDATA, XDATA, FDATA, ILEFT, DLEFT, IRIGHT,
C                          DRIGHT, BREAK, CSCOEF, IPVT)
C
C  Arguments:  (See CSDEC)
C
C  GAMS:       E1a
C
C  Chapter:    MATH/LIBRARY Interpolation and Approximation
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DC2DEC (NDATA, XDATA, FDATA, ILEFT, DLEFT, IRIGHT,
     &                   DRIGHT, BREAK, CSCOEF, IPVT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NDATA, ILEFT, IRIGHT, IPVT(*)
      DOUBLE PRECISION DLEFT, DRIGHT, XDATA(*), FDATA(*), BREAK(*),
     &           CSCOEF(4,*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, J, L, M
      DOUBLE PRECISION DIVDF1, DIVDF3, DTAU, G
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DC1SOR
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   N1RTY
      INTEGER    N1RTY
C
      CALL E1PSH ('DC2DEC')
C                                  CHECK ARGUMENT NDATA
      IF (NDATA .LT. 2) THEN
         CALL E1STI (1, NDATA)
         CALL E1MES (5, 1, 'The number of data points must be '//
     &               '2 or more while NDATA = %(I1) is given.')
      END IF
C                                  CHECK ARGUMENT ILEFT
      IF (ILEFT.LT.0 .OR. ILEFT.GT.2) THEN
         CALL E1STI (1, ILEFT)
         CALL E1MES (5, 2, 'The argument ILEFT = %(I1).  '//
     &               'It must be 0, 1 or 2.')
      END IF
C                                  CHECK ARGUMENT IRIGHT
      IF (IRIGHT.LT.0 .OR. IRIGHT.GT.2) THEN
         CALL E1STI (1, IRIGHT)
         CALL E1MES (5, 3, 'The argument IRIGHT = %(I1).  '//
     &               'It must be 0, 1 or 2.')
      END IF
      IF (N1RTY(0) .NE. 0) GO TO 9000
C                                  COPY AND SORT INPUT DATA
      CALL DC1SOR (NDATA, XDATA, FDATA, BREAK, CSCOEF, 4, IPVT)
      IF (N1RTY(0) .NE. 0) GO TO 9000
C                                  A TRIDIAGONAL LINEAR SYSTEM FOR THE
C                                    UNKNOWN SLOPES S(I) OF F AT
C                                    BREAK(I), I=1,...,N, IS GENERATED
C                                    AND THEN SOLVED BY GAUSS
C                                    ELIMINATION, WITH S(I) ENDING UP
C                                    IN CSCOEF(2,I), ALL I. CSCOEF(3,.)
C                                    AND CSCOEF(4,.) ARE USED INITIALLY
C                                    FOR TEMPORARY STORAGE.
      L = NDATA - 1
C                                  COMPUTE FIRST DIFFERENCES OF TAU
C                                    SEQUENCE AND STORE IN CSCOEF(3,.).
C                                    ALSO COMPUTE FIRST DIVIDED
C                                    DIFFERENCE OF DATA AND STORE IN
C                                    CSCOEF(4,.)
      DO 10  M=2, NDATA
         CSCOEF(3,M) = BREAK(M) - BREAK(M-1)
         CSCOEF(4,M) = (CSCOEF(1,M)-CSCOEF(1,M-1))/CSCOEF(3,M)
   10 CONTINUE
C                                  CONSTRUCT FIRST EQUATION FROM THE
C                                    BOUNDARY CONDITION, OF THE FORM
C                                    CSCOEF(4,1)*S(1) +
C                                    CSCOEF(3,1)*S(2) = DLEFT
      IF (ILEFT .EQ. 0) THEN
         IF (NDATA .EQ. 2) THEN
C                                  NO CONDITION AT LEFT END AND NDATA =
C                                    2
            CSCOEF(4,1) = 1.0D0
            CSCOEF(3,1) = 1.0D0
            CSCOEF(2,1) = 2.0D0*CSCOEF(4,2)
         ELSE
C                                  NOT-A-KNOT CONDITION AT LEFT END AND
C                                    NDATA .GT. 2
            CSCOEF(4,1) = CSCOEF(3,3)
            CSCOEF(3,1) = CSCOEF(3,2) + CSCOEF(3,3)
            CSCOEF(2,1) = ((CSCOEF(3,2)+2.0D0*CSCOEF(3,1))*CSCOEF(4,2)*
     &                    CSCOEF(3,3)+CSCOEF(3,2)**2*CSCOEF(4,3))/
     &                    CSCOEF(3,1)
         END IF
      ELSE IF (ILEFT .EQ. 1) THEN
C                                  SLOPE PRESCRIBED AT LEFT END
         CSCOEF(4,1) = 1.0D0
         CSCOEF(3,1) = 0.0D0
         CSCOEF(2,1) = DLEFT
      ELSE IF (ILEFT .EQ. 2) THEN
C                                  SECOND DERIVATIVE PRESCRIBED AT LEFT
C                                    END
         CSCOEF(4,1) = 2.0D0
         CSCOEF(3,1) = 1.0D0
         CSCOEF(2,1) = 3.0D0*CSCOEF(4,2) - CSCOEF(3,2)/2.0D0*DLEFT
      END IF
      IF (NDATA .GT. 2) THEN
C                                  IF THERE ARE INTERIOR KNOTS,
C                                    GENERATE THE CORRESP. EQUATIONS
C                                    AND CARRY OUT THE FORWARD PASS OF
C                                    GAUSS ELIMINATION, AFTER WHICH THE
C                                    M-TH EQUATION READS
C                                    CSCOEF(4,M)*S(M) +
C                                    CSCOEF(3,M)*S(M+1) = CSCOEF(2,M)
         DO 20  M=2, L
            G = -CSCOEF(3,M+1)/CSCOEF(4,M-1)
            CSCOEF(2,M) = G*CSCOEF(2,M-1) + 3.0D0*(CSCOEF(3,M)*
     &                    CSCOEF(4,M+1)+CSCOEF(3,M+1)*CSCOEF(4,M))
            CSCOEF(4,M) = G*CSCOEF(3,M-1) + 2.0D0*(CSCOEF(3,M)+
     &                    CSCOEF(3,M+1))
   20    CONTINUE
C                                  CONSTRUCT LAST EQUATION FROM THE
C                                    SECOND BOUNDARY CONDITION, OF THE
C                                    FORM (-G*CSCOEF(4,NDATA-1))*S(NDATA
C                                    A-1)+CSCOEF(5
         IF (IRIGHT .EQ. 0) THEN
            IF (NDATA.EQ.3 .AND. ILEFT.EQ.0) THEN
C                                  NDATA=3 AND NOT-A-KNOT ALSO AT LEFT
C
               CSCOEF(2,NDATA) = 2.0D0*CSCOEF(4,NDATA)
               CSCOEF(4,NDATA) = 1.0D0
               G = -1.0D0/CSCOEF(4,NDATA-1)
            END IF
         ELSE IF (IRIGHT .EQ. 1) THEN
C                                  SLOPE IS PRESCRIBED AT RIGHT END, GO
C                                    DIRECTLY TO BACKSUBSTITUTION,
C                                    SINCE C ARRAY HAPPENS TO BE SET UP
C                                    JUST RIGHT FOR IT AT THIS POINT
            CSCOEF(2,NDATA) = DRIGHT
         ELSE IF (IRIGHT .EQ. 2) THEN
C                                  SECOND DERIVATIVE PRESCRIBED AT
C                                    RIGHT
            CSCOEF(2,NDATA) = 3.0D0*CSCOEF(4,NDATA) +
     &                        CSCOEF(3,NDATA)/2.0D0*DRIGHT
            CSCOEF(4,NDATA) = 2.0D0
            G = -1.0D0/CSCOEF(4,NDATA-1)
         END IF
C
         IF (IRIGHT.EQ.0 .AND. .NOT.(NDATA.EQ.3.AND.ILEFT.EQ.0)) THEN
            G = CSCOEF(3,NDATA-1) + CSCOEF(3,NDATA)
            CSCOEF(2,NDATA) = ((CSCOEF(3,NDATA)+2.0D0*G)*CSCOEF(4,
     &                        NDATA)*CSCOEF(3,NDATA-1)+
     &                        CSCOEF(3,NDATA)**2*
     &                        (CSCOEF(1,NDATA-1)-CSCOEF(1,NDATA-2))/
     &                        CSCOEF(3,NDATA-1))/G
            G = -G/CSCOEF(4,NDATA-1)
            CSCOEF(4,NDATA) = CSCOEF(3,NDATA-1)
         END IF
      ELSE
C                                  NDATA = 2
         IF (IRIGHT .EQ. 0) THEN
            IF (ILEFT .GT. 0) THEN
C                                  NDATA=2 AND NOT NOT-A-KNOT AT LEFT
C                                    END
C
               CSCOEF(2,NDATA) = 2.0D0*CSCOEF(4,NDATA)
               CSCOEF(4,NDATA) = 1.0D0
            ELSE
C                                  NDATA=2 AND NOT-A-KNOT AT BOTH
C                                    ENDPOINTS
               CSCOEF(2,NDATA) = CSCOEF(4,NDATA)
            END IF
         ELSE IF (IRIGHT .EQ. 1) THEN
            CSCOEF(2,NDATA) = DRIGHT
         ELSE IF (IRIGHT .EQ. 2) THEN
C                                  SECOND DERIVATIVE PRESCRIBED AT
C                                    RIGHT
            CSCOEF(2,NDATA) = 3.0D0*CSCOEF(4,NDATA) +
     &                        CSCOEF(3,NDATA)/2.0D0*DRIGHT
            CSCOEF(4,NDATA) = 2.0D0
         END IF
         G = -1.0D0/CSCOEF(4,NDATA-1)
      END IF
C
      IF (IRIGHT.NE.1 .AND. (NDATA.GT.2.OR.IRIGHT.NE.0.OR.ILEFT.NE.0))
     &    THEN
C                                  COMPLETE FORWARD PASS OF GAUSS
C                                    ELIMINATION
C
         CSCOEF(4,NDATA) = G*CSCOEF(3,NDATA-1) + CSCOEF(4,NDATA)
         CSCOEF(2,NDATA) = (G*CSCOEF(2,NDATA-1)+CSCOEF(2,NDATA))/
     &                     CSCOEF(4,NDATA)
      END IF
C                                  CARRY OUT BACK SUBSTITUTION
      DO 30  J=L, 1, -1
         CSCOEF(2,J) = (CSCOEF(2,J)-CSCOEF(3,J)*CSCOEF(2,J+1))/
     &                 CSCOEF(4,J)
   30 CONTINUE
C                                  GENERATE CUBIC COEFFICIENTS IN EACH
C                                    INTERVAL, I.E., THE DERIV.S AT ITS
C                                    LEFT ENDPOINT, FROM VALUE AND
C                                    SLOPE AT ITS ENDPOINTS.
      DO 40  I=2, NDATA
         DTAU = CSCOEF(3,I)
         DIVDF1 = (CSCOEF(1,I)-CSCOEF(1,I-1))/DTAU
         DIVDF3 = CSCOEF(2,I-1) + CSCOEF(2,I) - 2.0D0*DIVDF1
         CSCOEF(3,I-1) = 2.0D0*(DIVDF1-CSCOEF(2,I-1)-DIVDF3)/DTAU
         CSCOEF(4,I-1) = (DIVDF3/DTAU)*(6.0D0/DTAU)
   40 CONTINUE
C
 9000 CONTINUE
      CALL E1POP ('DC2DEC ')
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  C1SOR/DC1SOR (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 5, 1986
C
C  Purpose:    Copy and sort XDATA into BREAK and FDATA into the first
C              row of CSCOEF.
C
C  Usage:      CALL C1SOR (NDATA, XDATA, FDATA, BREAK, CSCOEF, LDC,
C                          IWORK)
C
C  Arguments:
C     NDATA  - Number of data points (must be at least 2).  (Input)
C     XDATA  - Array of length NDATA containing the data point
C              abscissas.  (Input)
C     FDATA  - Array of length NDATA containing the data point
C              ordinates.  (Input)
C     BREAK  - Array of length NDATA containing the breakpoints
C              for the piecewise cubic representation.  (Output)
C     CSCOEF - Matrix of size 4 by NDATA containing the cubic spline
C              (piecewise polynomial) representation in the first
C              NDATA-1 columns.  (Output)
C     LDC    - Leading dimension of CSCOEF exactly as specified in the
C              dimension statement of the calling program.  (Input)
C     IWORK  - Work array containing the permutation vector.  (Output)
C
C  GAMS:       E1a
C
C  Chapter:    MATH/LIBRARY Interpolation and Approximation
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DC1SOR (NDATA, XDATA, FDATA, BREAK, CSCOEF, LDC,
     &                   IWORK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NDATA, LDC, IWORK(*)
      DOUBLE PRECISION XDATA(*), FDATA(*), BREAK(*), CSCOEF(LDC,*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, J
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, E1STD, DCOPY, DSVRGP
C                                  CHECK FOR SORTED ARRAY
      CALL E1PSH ('DC1SOR ')
      DO 10  I=2, NDATA
         IF (XDATA(I-1) .GE. XDATA(I)) THEN
C                                  CHECK THAT XDATA VALUES ARE DISTINCT
            IF (XDATA(I-1) .EQ. XDATA(I)) THEN
               J = I - 1
               CALL E1STI (1, J)
               CALL E1STI (2, I)
               CALL E1STD (1, XDATA(I))
               CALL E1MES (5, 2, 'Points in the data point '//
     &                     'abscissas array, XDATA, must be '//
     &                     'distinct, but XDATA(%(I1)) = '//
     &                     'XDATA(%(I2)) = %(D1).')
               GO TO 9000
            ELSE
               GO TO 20
            END IF
         END IF
   10 CONTINUE
C                                  DATA IS ALREADY SORTED.  M1VE TO
C                                  BREAK AND CSCOEF.
      CALL DCOPY (NDATA, XDATA, 1, BREAK, 1)
      CALL DCOPY (NDATA, FDATA, 1, CSCOEF, LDC)
      GO TO 9000
C                                  SET INITIAL PERMUTATION
   20 DO 30  I=1, NDATA
         IWORK(I) = I
   30 CONTINUE
C                                  FIND SORTING PERMUTATION
      CALL DSVRGP (NDATA, XDATA, BREAK, IWORK)
C                                  APPLY PERMUTATION
      DO 40  I=1, NDATA
         CSCOEF(1,I) = FDATA(IWORK(I))
   40 CONTINUE
C                                  CHECK THE XDATA VALUES ARE DISTINCT
      DO 50  I=2, NDATA
         IF (BREAK(I-1) .EQ. BREAK(I)) THEN
            CALL E1STI (1, IWORK(I-1))
            CALL E1STI (2, IWORK(I))
            CALL E1STD (1, BREAK(I))
            CALL E1MES (5, 2, 'Points in the data point abscissas '//
     &                  'array, XDATA, must be distinct, '//
     &                  'but XDATA(%(I1)) = XDATA(%(I2)) = %(D1).')
            GO TO 9000
         END IF
   50 CONTINUE
C
 9000 CONTINUE
      CALL E1POP ('DC1SOR')
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  P3DER/DP3DER (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 11, 1986
C
C  Purpose:    Compute MAX (I, 1 .LE. NINTV. .AND. BREAK(I) .LE. X)
C
C  Usage:      CALL P3DER (KORD, NINTV, BREAK, X, LEFT)
C
C  Arguments:
C     KORD   - Order of the polynomial.  (Input)
C     NINTV  - Number of polynomial pieces.  (Input)
C     BREAK  - Vector of length NINTV+1 containing the breakpoints
C              of the piecewise polynomial representation.  (Input)
C     X      - The point whose location in BREAK is to be found.
C     LEFT   - Integer whose value is
C                LEFT
C                  1      IF                       X .LT.  BREAK(1)
C                  I      IF         BREAK(I).LE. X .LT. BREAK(I+1)
C                 NINTV   IF                    BREAK(NINTV) .LE. X
C              The asymmetric treatment of the interval is due to the
C              decision to make all PP functions continuous from the
C              right.  (Output)
C
C  Remark:
C     This routine is based in INTERV in de Boor, p92-93.
C
C  Chapter:    MATH/LIBRARY Interpolation and Approximation
C
C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DP3DER (KORD, NINTV, BREAK, X, LEFT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    KORD, NINTV, LEFT
      DOUBLE PRECISION X, BREAK(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IHI, ISTEP, MIDDLE
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    ILO
      SAVE       ILO
C
      DATA ILO/1/
C
      IHI = ILO + 1
      IF (IHI .GE. NINTV) THEN
         IF (X .GE. BREAK(NINTV)) THEN
            LEFT = NINTV
            GO TO 9000
         ELSE IF (NINTV .LE. 1) THEN
            LEFT = 1
            GO TO 9000
         END IF
         ILO = NINTV - 1
         IHI = NINTV
      END IF
C
      IF (X .LT. BREAK(IHI)) THEN
         IF (X .GE. BREAK(ILO)) THEN
            LEFT = ILO
            GO TO 9000
         END IF
C                                  Now X .LT. BREAK(ILO) - decrease ILO
C                                  to capture X
         ISTEP = 1
   10    CONTINUE
         IHI = ILO
         ILO = IHI - ISTEP
         IF (ILO .GT. 1) THEN
            IF (X .GE. BREAK(ILO)) GO TO 30
            ISTEP = ISTEP*2
            GO TO 10
         END IF
         ILO = 1
         IF (X .LT. BREAK(1)) THEN
            LEFT = 1
            GO TO 9000
         END IF
         GO TO 30
      END IF
C                                  Now X .GE. BREAK(IHI) - increase IHI
C                                  to capture X
      ISTEP = 1
   20 CONTINUE
      ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI .LT. NINTV) THEN
         IF (X .LT. BREAK(IHI)) GO TO 30
         ISTEP = ISTEP*2
         GO TO 20
      END IF
      IF (X .GE. BREAK(NINTV)) THEN
         LEFT = NINTV
         GO TO 9000
      END IF
      IHI = NINTV
C                                  Now BREAK(ILO) .LE. X .LT.
C                                  BREAK(IHI) - narrow the inteval
   30 CONTINUE
      MIDDLE = (ILO+IHI)/2
      IF (MIDDLE .EQ. ILO) THEN
         LEFT = ILO
         GO TO 9000
      END IF
C                                  It is assumed that MIDDLE = ILO in
C                                  case IHI = ILO+1
      IF (X .LT. BREAK(MIDDLE)) THEN
         IHI = MIDDLE
      ELSE
         ILO = MIDDLE
      END IF
      GO TO 30
 9000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  SVRGP/DSVRGP (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    February 19, 1985
C
C  Purpose:    Sort a real array by algebraic values and return the
C              permutations.
C
C  Usage:      CALL SVRGP (N, RA, RB, IPERM)
C
C  Arguments:
C     N      - Number of elements in the array to be sorted.  (Input)
C     RA     - Vector of length N containing the array to be sorted.
C              (Input)
C     RB     - Vector of length N containing the sorted array.  (Output)
C              If RA is not needed, RA and RB can share the same
C              storage locations.
C     IPERM  - Vector of length N.  (Input/Output)
C              On input IPERM should be initialized to the values
C              1, 2, ... N.  On output IPERM contains a record of
C              permutations made on the vector RA.
C
C  Remark:
C     For wider applicability, integers (1, 2, ... N) that are to be
C     associated with RA(I) for I = 1, 2, ... N may be entered into
C     IPERM(I) in any order.  Note that these integers must be unique.
C
C  GAMS:       N6a1b
C
C  Chapters:   MATH/LIBRARY Utilities
C              STAT/LIBRARY Utilities
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DSVRGP (N, RA, RB, IPERM)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, IPERM(*)
      DOUBLE PRECISION RA(*), RB(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IJ, IL(21), IT, ITT, IU(21), J, K, L, M
      DOUBLE PRECISION R, T, TT
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DCOPY
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   N1RCD
      INTEGER    N1RCD
C
      CALL E1PSH ('DSVRGP ')
      CALL DCOPY (N, RA, 1, RB, 1)
C                                  CHECK FOR INPUT ERRORS
      IF (N .LE. 0) THEN
         CALL E1STI (1, N)
         CALL E1MES (5, 1, 'The length of vectors RA, RB and '//
     &               'IPERM must be greater than or equal to one. '//
     &               ' On input the length, N, is given as %(I1).')
      END IF
      IF (N1RCD(0) .NE. 0) GO TO 9000
C
      M = 1
      I = 1
      J = N
      R = .375D0
   10 IF (I .EQ. J) GO TO 70
      IF (R .LE. .5898437D0) THEN
         R = R + 3.90625D-2
      ELSE
         R = R - .21875D0
      END IF
   20 K = I
C                                  SELECT A CENTRAL ELEMENT OF THE
C                                  ARRAY AND SAVE IT IN LOCATION T
      IJ = I + (J-I)*R
      T = RB(IJ)
      IT = IPERM(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (RB(I) .GT. T) THEN
         RB(IJ) = RB(I)
         RB(I) = T
         T = RB(IJ)
         IPERM(IJ) = IPERM(I)
         IPERM(I) = IT
         IT = IPERM(IJ)
      END IF
      L = J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  T, INTERCHANGE WITH T
      IF (RB(J) .GE. T) GO TO 40
      RB(IJ) = RB(J)
      RB(J) = T
      T = RB(IJ)
      IPERM(IJ) = IPERM(J)
      IPERM(J) = IT
      IT = IPERM(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (RB(I) .LE. T) GO TO 40
      RB(IJ) = RB(I)
      RB(I) = T
      T = RB(IJ)
      IPERM(IJ) = IPERM(I)
      IPERM(I) = IT
      IT = IPERM(IJ)
      GO TO 40
   30 IF (RB(L) .EQ. RB(K)) GO TO 40
      TT = RB(L)
      RB(L) = RB(K)
      RB(K) = TT
      ITT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = ITT
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN T
   40 L = L - 1
      IF (RB(L) .GT. T) GO TO 40
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN T
   50 K = K + 1
      IF (RB(K) .LT. T) GO TO 50
C                                  INTERCHANGE THESE ELEMENTS
      IF (K .LE. L) GO TO 30
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I .LE. J-K) GO TO 60
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 80
   60 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 80
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
   70 M = M - 1
      IF (M .EQ. 0) GO TO 9000
      I = IL(M)
      J = IU(M)
   80 IF (J-I .GE. 11) GO TO 20
      IF (I .EQ. 1) GO TO 10
      I = I - 1
   90 I = I + 1
      IF (I .EQ. J) GO TO 70
      T = RB(I+1)
      IT = IPERM(I+1)
      IF (RB(I) .LE. T) GO TO 90
      K = I
  100 RB(K+1) = RB(K)
      IPERM(K+1) = IPERM(K)
      K = K - 1
      IF (T .LT. RB(K)) GO TO 100
      RB(K+1) = T
      IPERM(K+1) = IT
      GO TO 90
C
 9000 CALL E1POP ('DSVRGP ')
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

C-----------------------------------------------------------------------
C  IMSL Name:  CSDEC/DCSDEC (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 1, 1986
C
C  Purpose:    Compute the cubic spline interpolant with specified
C              derivative endpoint conditions.
C
C  Usage:      CALL CSDEC (NDATA, XDATA, FDATA, ILEFT, DLEFT, IRIGHT,
C                          DRIGHT, BREAK, CSCOEF)
C
C  Arguments:
C     NDATA  - Number of data points.  (Input)
C     XDATA  - Array of length NDATA containing the data point
C              abscissas.  (Input)
C              The data point abscissas must be distinct.
C     FDATA  - Array of length NDATA containing the data point
C              ordinates.  (Input)
C     ILEFT  - Type of end condition at the left endpoint.  (Input)
C                   ILEFT    Condition
C                    0       'Not-a-knot' condition
C                    1       First derivative specified by DLEFT
C                    2       Second derivative specified by DLEFT
C     DLEFT  - Derivative at left endpoint if ILEFT is equal to
C              1 or 2.  (Input)
C              If ILEFT = 0 then DLEFT is ignored.
C     IRIGHT - Type of end condition at the right endpoint.  (Input)
C                   IRIGHT    Condition
C                    0       'Not-a-knot' condition
C                    1       First derivative specified by DRIGHT
C                    2       Second derivative specified by DRIGHT
C     DRIGHT - Derivative at right endpoint if IRIGHT is equal to
C              1 or 2.  (Input)
C              If IRIGHT = 0 then DRIGHT is ignored.
C     BREAK  - Array of length NDATA containing the breakpoints for the
C              piecewise cubic representation.  (Output)
C     CSCOEF - Matrix of size 4 by NDATA containing the local
C              coefficients of the cubic pieces.  (Output)
C
C  Remarks:
C  1. Automatic workspace usage is
C              CSDEC    NDATA units, or
C              DCSDEC   NDATA units.
C     Workspace may be explicitly provided, if desired, by use of
C     C2DEC/DC2DEC.  The reference is
C              CALL C2DEC (NDATA, XDATA, FDATA, ILEFT, DLEFT, IRIGHT,
C                          DRIGHT, BREAK, CSCOEF, IWK)
C     The additional argument is
C     IWK    - Work array of length NDATA.
C
C  2. The cubic spline can be evaluated using CSVAL; its derivative can
C     be evaluated using CSDER.
C
C  3. Note that column NDATA of CSCOEF is used as workspace.
C
C  GAMS:       E1a
C
C  Chapter:    MATH/LIBRARY Interpolation and Approximation
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DCSDEC (NDATA, XDATA, FDATA, ILEFT, DLEFT, IRIGHT,
     &                   DRIGHT, BREAK, CSCOEF)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NDATA, ILEFT, IRIGHT
      DOUBLE PRECISION DLEFT, DRIGHT, XDATA(*), FDATA(*), BREAK(*),
     &           CSCOEF(4,*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    INDATA
C                                  SPECIFICATIONS FOR SPECIAL CASES
C                                  SPECIFICATIONS FOR COMMON /WORKSP/
      REAL       RWKSP(5000)
      DOUBLE PRECISION RDWKSP(2500)
      DOUBLE PRECISION DWKSP(2500)
      COMPLEX    CWKSP(2500)
      COMPLEX    *16 CZWKSP(1250)
      COMPLEX    *16 ZWKSP(1250)
      INTEGER    IWKSP(5000)
      LOGICAL    LWKSP(5000)
      EQUIVALENCE (DWKSP(1), RWKSP(1))
      EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
      EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
      EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
      COMMON     /WORKSP/ DWKSP
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DC2DEC
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   I1KGT, N1RTY
      INTEGER    I1KGT, N1RTY
C
      CALL E1PSH ('DCSDEC')
C                                  CHECK ARGUMENT NDATA
      IF (NDATA .LT. 2) THEN
         CALL E1STI (1, NDATA)
         CALL E1MES (5, 1, 'The number of data points must be '//
     &               '2 or more while NDATA = %(I1) is given.')
      END IF
C                                  CHECK ARGUMENT ILEFT
      IF (ILEFT.LT.0 .OR. ILEFT.GT.2) THEN
         CALL E1STI (1, ILEFT)
         CALL E1MES (5, 2, 'The argument ILEFT = %(I1).  '//
     &               'It must be 0, 1 or 2.')
      END IF
C                                  CHECK ARGUMENT IRIGHT
      IF (IRIGHT.LT.0 .OR. IRIGHT.GT.2) THEN
         CALL E1STI (1, IRIGHT)
         CALL E1MES (5, 3, 'The argument IRIGHT = %(I1).  '//
     &               'It must be 0, 1 or 2.')
      END IF
      IF (N1RTY(0) .NE. 0) GO TO 9000
      INDATA = I1KGT(NDATA,2)
      IF (N1RTY(0) .NE. 0) THEN
         CALL E1MES (5, 5, ' ')
         CALL E1STI (1, NDATA)
         CALL E1MES (5, 5, 'Workspace allocation is based on '//
     &               'NDATA = %(I1).')
      ELSE
C
         CALL DC2DEC (NDATA, XDATA, FDATA, ILEFT, DLEFT, IRIGHT,
     &                DRIGHT, BREAK, CSCOEF, IWKSP(INDATA))
      END IF
C
 9000 CONTINUE
      CALL E1POP ('DCSDEC ')
      RETURN
      END
C-----------------------------------------------------------------------
C  IMSL Name:  CSVAL/DCSVAL (Single/Double precision version)
C
C  Computer:   VAX/DOUBLE
C
C  Revised:    August 1, 1986
C
C  Purpose:    Evaluate a cubic spline.
C
C  Usage:      CSVAL(X, NINTV, BREAK, CSCOEF)
C
C  Arguments:
C     X      - Point at which the spline is to be evaluated.
C              (Input)
C     NINTV  - Number of polynomial pieces.  (Input)
C     BREAK  - Array of length NINTV+1 containing the breakpoints
C              for the piecewise cubic representation.  (Input)
C              BREAK must be strictly increasing.
C     CSCOEF - Matrix of size 4 by NINTV+1 containing the local
C              coefficients of the cubic pieces.  (Input)
C     CSVAL  - Value of the polynomial at X.  (Output)
C
C  GAMS:       E3; K6
C
C  Chapter:    MATH/LIBRARY Interpolation and Approximation
C
C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DCSVAL (X, NINTV, BREAK, CSCOEF)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NINTV
      DOUBLE PRECISION X, BREAK(*), CSCOEF(4,*)
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   E1POP, E1PSH
C                                  SPECIFICATIONS FOR FUNCTIONS
      EXTERNAL   DPPDER
      DOUBLE PRECISION DPPDER
C
      CALL E1PSH ('DCSVAL ')
C                                  CUBIC SPLINE IS A SPECIAL PP
      DCSVAL = DPPDER(0,X,4,NINTV,BREAK,CSCOEF)
C
 9000 CONTINUE
      CALL E1POP ('DCSVAL ')
      RETURN
      END
