      subroutine cubeq(a, b, c, roots)

      implicit none
      complex*16 a, b, c, roots(3)

c     Solves the cubic equation with complex coefficients
c     x**3 + a x**2 + b x + c = 0
c     Ref: Numerical Recipes 2nd Ed. p.179

      real*8 on, off, ot, sq3o2, sig
      complex*16 a2, q, r, s, ba, bb, cz, ao3, t1, t2, t3

c     1/9, 1/54, 1/3, sqrt(3)/2:
      data on/0.1111111111111111d0/,   off/1.851851851851852d-2/
     ;   , ot/0.3333333333333333d0/, sq3o2/8.660254037844387d-1/
     ;   , cz/(0.d0,0.d0)/

      a2 = a * a
      q = (a2 - 3.d0 * b) * on
      r = ((2.d0*a2-9.d0*b) * a + 27.d0 * c) * off

      s = cdsqrt(r * r - q**3)
      sig = dsign(1.d0, dreal(dconjg(r)*s))
      s = sig * s

      ba = - (r + s) ** ot
      bb = cz
      if(ba .ne. cz)bb = q / ba

      ao3 = a * ot
      t1 = ba + bb
      t2 = (ba - bb) * dcmplx(0.d0, sq3o2)
      t3 = -(0.5d0 * t1 + ao3)

      roots(1) = t1 - ao3
      roots(2) = t3 + t2
      roots(3) = t3 - t2

      return
      end