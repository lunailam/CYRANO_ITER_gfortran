      integer function isrchfge(n, x, incx, targ)
      
      integer n, incx, i, j
      double precision x(1+(n-1)*incx), targ
      
c     search for first element of ordered array x >= targ
c     return its index

      if(targ.gt.x(1+(n-1)*incx))then
      isrchfge = 1 + n
      return
      end if
      
      i = 1
      j = 1
   1  continue
      if(x(i).ge.targ)then
      isrchfge = j
      return
      else if(j.lt.n)then
      j = j + 1
      i = i + incx
      go to 1
      end if
      
      end