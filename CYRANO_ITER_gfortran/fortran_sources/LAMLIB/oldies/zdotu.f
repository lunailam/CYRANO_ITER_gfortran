      double complex function zdotu(n,zx,incx,zy,incy)
c
c     Alias for Cray single precision
c
      double complex zx(*),zy(*)
      complex cdotu
      integer incx,incy,n
      zdotu = cdotu(n,zx,incx,zy,incy)
      return
      end
