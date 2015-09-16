      SUBROUTINE CGEFA3(A, LDA, N, M, NPVT, IPVT, SCALE, SCIN, INFO)

      IMPLICIT NONE

      INTEGER LDA, N, M, NPVT, IPVT(M), INFO, SCIN
      DOUBLE PRECISION SCALE(N)
      COMPLEX*16 A(LDA,M)
C
C     CGEFA3 PARTIALLY FACTORS A N BY M COMPLEX MATRIX BY GAUSSIAN ELIMINATION.
C     (M<=N)
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, M)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE NUMBER OF ROWS OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF COLUMNS OF THE MATRIX  A .
C
C        NPVT    INTEGER
C                THE MAX.NUMBER OF ROWS TO BE USED FOR PIVOTING .
C
C        SCIN    INTEGER
C                A SWITCH FOR THE ROW INVERSE SCALING FACTOR SCALE:
C                SCIN=0: SCALE IS OUTPUT
C                SCIN NONZERO: SCALE IS AN INPUT.
c                In this case one must make sure that no element of
c                SCALE is zero (this produces errors in 
c                pivot selection)
C
C        SCALE   REAL(NPVT)
C                INPUT WHEN SCIN IS NONZERO. SCALE MUST THEN CONTAIN
C                THE INVERSE SCALING FACTOR OF EACH EQUATION TO BE USED
C                IN SEARCH OF PIVOTS.
C
C
C     ON RETURN
C
C        A       A 'PARTIALLY UPPER TRIANGULAR' MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(M)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        SCALE   REAL(N)
C                WHEN SCIN=0, THE MAXIMUM NORM OF EACH ROW OF A
C                USING THE CABS1 FUNCTION FOR ABSOLUTE VALUES .
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT CGESL OR CGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN CGECO2 FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     MODIFIED VERSION OF LINPACK CGEFA. THIS VERSION DATED 1/10/92 .
C     PURPOSE: PARTIAL FACTORIZATION;
C              SCALED ROW PIVOTING. THE SCALING FACTORS MAY BE INPUT
C              OR OUTPUT;
C              CHOICE OF PIVOTS MAY BE RESTRICTED TO LESS THAN N ROWS.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS ZAXPY,CSCAL,izamax
C     FORTRAN ABS,AIMAG,REAL
C
C     INTERNAL VARIABLES
C
      INTEGER izamax, I, J, K, KP1, L, NM1, MM1, NP
      DOUBLE PRECISION THR, THR2
      COMPLEX*16 T, ONE
      EXTERNAL izamax
C
      COMPLEX*16 ZDUM
      DOUBLE PRECISION CABS1
      DATA ONE/(1.0D0,0.0D0)/
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))

      NP = MIN0(NPVT, N)
C
C     ROW MAXIMUM NORMS ARE COMPUTED WHEN SCIN=0:
C
      IF(SCIN.EQ.0)THEN
      DO 1 I = 1, NP
      J = IZAMAX(M, A(I,1), LDA)
      SCALE(I) = CABS1( A(I,J) )
c     When partial factoriz. is performed, zero rows can occur.
c     They mustn't spoil the choice of pivots, hence
      if(scale(i) .eq. 0.d0)scale(i) = 1.d0
  1   CONTINUE
      END IF
C
C     GAUSSIAN ELIMINATION WITH RESTRICTED SCALED PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      MM1 = M 
      IF(M.GE.N)MM1 = NM1
      
      IF (MM1 .LT. 1) GO TO 70
      DO 60 K = 1, MM1
      KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
      L = K
      THR = CABS1(A(K,K))
        DO 2 I = KP1, NP
        THR2 = CABS1(A(I,K))
          IF(THR2*SCALE(L).GT.THR*SCALE(I))THEN
          L = I
          THR = THR2
          END IF
  2     CONTINUE
      IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (THR .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .NE. K) THEN
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
            END IF
C
C           COMPUTE MULTIPLIERS
C
            T = -ONE/A(K,K)
            CALL zscal(N-K,T,A(KP1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, M
               T = A(L,J)
               IF (L .NE. K) THEN
                  A(L,J) = A(K,J)
                  A(K,J) = T
               END IF
               CALL ZAXPY(N-K,T,A(KP1,K),1,A(KP1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE

      IF(M.EQ.N .and. NPVT.GE.N)THEN
        IPVT(N) = N
        IF (CABS1(A(N,N)) .EQ. 0.0D0) INFO = N
      END IF
      
      RETURN
      END
