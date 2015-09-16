      INTEGER FUNCTION IIMIN(N, IX, incx)

      IMPLICIT NONE

C
C     Returns the index of the smallest component in integer vector
C     IX(N)
C
      INTEGER N, IX(N), incx
      
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
      IF(IX(k).LT.IX(j))j=I
  1   CONTINUE
      END IF
      IIMIN = j
      RETURN
      END
