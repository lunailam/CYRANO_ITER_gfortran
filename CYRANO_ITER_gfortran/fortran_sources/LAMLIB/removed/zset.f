      subroutine zset(n,a,sy,incy)
c
c     Set all components of complex*16 array sy to complex*16 scalar s
c
      double complex a,sy(*)
      integer incy,n

      integer i
c
c      call cset(n,a,sy,incy)
      if (incy .gt. 0) then
        do i = 1, n
        sy(1+(i-1)*incy) = a  
        end do
      end if

      return
      end
