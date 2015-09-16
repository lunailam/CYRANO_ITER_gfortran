      INTEGER FUNCTION ISUM(N, IVEC, INCV)
      INTEGER N, IVEC(*), INCV, I, J
      ISUM = 0
      IF(N.LE.0)RETURN
      J = 1
      IF(INCV.LT.0)J = (-N+1)*INCV + 1
       DO 1 I = 1, N
      ISUM = ISUM + IVEC(J)
      J = J + INCV
    1 CONTINUE
      RETURN
      END
