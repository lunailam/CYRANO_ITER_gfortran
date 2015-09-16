      SUBROUTINE WRP1D(NOFILE, IPLOT, N, X, Y, LDY, INCY, NCURV, TITLE
     ;, XUNIT, XLABEL, YUNIT, YLABEL, NTEXT, TEXT, HEADER)
      
      IMPLICIT NONE
      
      LOGICAL HEADER
      
      INTEGER NOFILE, IPLOT, N, LDY, INCY, NCURV, NTEXT
      
      DOUBLE PRECISION X(N), Y(LDY,1+(NCURV-1)*incy)
      
      CHARACTER*48 TITLE, XLABEL, YLABEL(NCURV), TEXT(NTEXT)
      
      CHARACTER*10 XUNIT, YUNIT
C
C     Writes 1d plot data to output file. 
C     There are NCURV curves and N abscissae.
C
      INTEGER I, j

      WRITE(NOFILE,1001)TITLE
      WRITE(NOFILE,*)N, NCURV+1, NTEXT
      
      if(header)then
c      WRITE(NOFILE,1001)IPLOT, TITLE
c      WRITE(NOFILE,*)N, NCURV, NTEXT
      WRITE(NOFILE,1002)XLABEL, (YLABEL(j),j=1,NCURV)
      WRITE(NOFILE,1003)XUNIT, YUNIT
      IF(NTEXT.GT.0)WRITE(NOFILE,*)(TEXT(j),j=1,NTEXT)
      end if
      
      DO I = 1, N
      WRITE(NOFILE,1000) X(I), (Y(I,1+(j-1)*INCY),j=1,NCURV)
      END DO
      
      RETURN
 1000 format(1h ,10(1x, g15.7))
 1001 format(1h ,a48)
 1002 format(1h ,10(1x,a48))
 1003 format(1h ,2(1x,a10))
      END
