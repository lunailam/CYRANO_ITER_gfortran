      INTEGER FUNCTION IIMAX(N, IX, INCX)

      IMPLICIT NONE
C
C     Returns the index of the largest component in integer vector IX(N)
C
      INTEGER N, IX(N), INCX

      INTEGER I, j, k
      
      IF(N.LE.0)THEN
      j = 0
      ELSE IF(N.EQ.1)THEN
      j = 1
      ELSE
      j = 1
      k = 1
      DO 1 I = 2, N
      k = k + incx
      IF(IX(k).GT.IX(j))j=I
  1   CONTINUE
      END IF
      IIMAX = j
      RETURN
      END
