      SUBROUTINE WRP2D(NOFILE, IPLOT, Nx, X, Ny, Y, Z, LDZ, INCX, INCY, TITLE
     ;, XUNIT, XLABEL, YUNIT, YLABEL, ZUNIT, ZLABEL, NTEXT, TEXT)

c     stores data for plot of z(x,y). Data on a rectangular mesh.
      
      IMPLICIT NONE
      
      INTEGER NOFILE, IPLOT, Nx, Ny, LDZ, INCX, INCY, NTEXT
      DOUBLE PRECISION X(Nx), Y(Ny), Z(ldz,*)
      CHARACTER*48 TITLE, XLABEL, YLABEL, ZLABEL, TEXT(NTEXT)
      CHARACTER*10 XUNIT, YUNIT, ZUNIT

      INTEGER i, j
      
      WRITE(NOFILE,*)IPLOT, TITLE
      WRITE(NOFILE,*)Nx, Ny, NTEXT
      WRITE(NOFILE,*)XLABEL, YLABEL, ZLABEL
      WRITE(NOFILE,*)XUNIT, YUNIT, ZUNIT
      IF(NTEXT.GT.0)WRITE(NOFILE,*)(TEXT(j),j=1,NTEXT)
      
      WRITE(NOFILE,*) (X(1+(i-1)*INCX),i=1,Nx)
      WRITE(NOFILE,*) (Y(1+(j-1)*INCY),j=1,Ny)
      DO i = 1, Nx
      WRITE(NOFILE,*) (Z(1+(i-1)*INCX,1+(j-1)*INCY),j=1,Ny)
      END DO
      
      RETURN
      END
