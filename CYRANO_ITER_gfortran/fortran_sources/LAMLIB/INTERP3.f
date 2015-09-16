      subroutine interp3(xtab, incx, reftab, incref, nref, xabs, tab, inctab, ntab)

      implicit none

      integer nref, ntab, incx, incref, inctab

      double precision xtab(1+(nref-1)*incx), reftab(1+(nref-1)*incref), xabs(ntab), tab(1+(ntab-1)*inctab)

c     3-point (quadratic) interpolation in table:
c     Non-equidistant case.
C     Reference table of function: REFTAB(*), NREF distinct abscissae stored
C     in increasing order in XTAB. Increments incx, INCREF between values in XTAB 
C     and REFTAB, respectively.
C     XABS(NTAB): list of increasing abscissae between 1st and last element of
C     XTAB, where the interpolation has to be done.
C     Output: TAB(1:1+(ntab-1)*inctab),
C     interpolation of the list in REFTAB at the NTAB points XABS,
c     stored with increment inctab allowing output to the rows of a 2D array.
 
      integer i, itab, itab1, itab2, j, ix, ix1, ix2, isrchfge

      double precision x, x1, x2, x3, a1, a2, a3, x12, x23, x13, x1x, x2x, x3x

      external isrchfge

      if(nref.lt.3)write(6,*)'error interp2: the reference table only has', nref, ' entries' 
      if(xabs(1).lt.xtab(1) .or. xabs(ntab).gt.xtab(1+(nref-1)*incx))
     ;write(6,*)'interp2: warning: some values lie outside interpolating range'
     
      do 1 i = 1, ntab
      x = xabs(i)
      j = min0(isrchfge(nref, xtab, incx, x), nref-1)
      j = max0(2, j)
      ix = 1 + (j - 1) * incx
      ix1 = ix - incx
      ix2 = ix + incx
      itab = 1 + (j - 1) * incref
      itab1 = itab - incref
      itab2 = itab + incref
      x1 = xtab(ix1)
      x2 = xtab(ix)
      x3 = xtab(ix2)
      x12 = x1 - x2
      x23 = x2 - x3
      x13 = x1 - x3
      x1x = x1 - x
      x2x = x2 - x
      x3x = x3 - x
      a1 = x2x * x3x / (x12 * x13)
      a2 = x3x * x1x / (- x23 * x12)
      a3 = x1x * x2x / (x13 * x23)
c      a1 = (x2-x)*(x3-x)/((x2-x1)*(x3-x1))
c      a2 = (x3-x)*(x1-x)/((x3-x2)*(x1-x2))
c      a3 = (x1-x)*(x2-x)/((x1-x3)*(x2-x3))
      tab(1+(i-1)*inctab) = a1 * reftab(itab1) + a2 * reftab(itab) + a3 * reftab(itab2)
  1   continue

      return
      end
