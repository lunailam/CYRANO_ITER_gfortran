      subroutine iset(n,ia,sy,incy)
c
c     sets components of a vector, y, to a value ia.
c     uses unrolled loops for increments equal to 1.
c     based on
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      integer ia,sy(*)
      integer i,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incy.eq.1)go to 20
c
c        code for increment not equal to 1
c
      iy = 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = ia
        iy = iy + incy
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = ia
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = ia
        sy(i + 1) = ia
        sy(i + 2) = ia
        sy(i + 3) = ia
        sy(i + 4) = ia
        sy(i + 5) = ia
        sy(i + 6) = ia
   50 continue
      return
      end
