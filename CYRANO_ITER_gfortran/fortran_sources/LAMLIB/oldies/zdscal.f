      subroutine zdscal(n,da,zx,incx)
c
c     scales a vector by a constant.
c
c     Alias for Cray single precision
c
      double complex zx(*)
      double precision da
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      call csscal(n,da,zx,incx)
      return
      end

