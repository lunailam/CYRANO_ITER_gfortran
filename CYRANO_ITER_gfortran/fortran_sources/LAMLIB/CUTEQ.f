      SUBROUTINE CUTEQ(X0, X1, N, X)

      IMPLICIT NONE

      INTEGER N

      DOUBLE PRECISION X0, X1, X(0:N)

C  DIVIDES X0-X1 INTO N EQUAL ELEMENTS

c      include 'COMFIC.COPY'

      INTEGER I

      DOUBLE PRECISION H

      X(0) = X0
      H = (X1 - X0) / DFLOAT(N)

      DO 1 I = 1, N - 1
      X(I) = X(I-1) + H
  1   CONTINUE

      X(N) = X1

      RETURN
      END
