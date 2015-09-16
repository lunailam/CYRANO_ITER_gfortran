      subroutine interp1(x1, x2, reftab, incref, nref, xabs, tab, ntab, type)

      implicit none

      character*1 type 
      
      integer nref, ntab, incref

      double precision x1, x2, reftab(1+(nref-1)*incref), xabs(ntab), tab(ntab)

C     Interpolate in table: 3-point (quadratic) Lagrange
C     Equidistant case.
C     Reference table of function: REFTAB(*), NREF equidistant abscissae
C     first X1, last X2. Increment INCREF between values in REFTAB.
C     XABS(NTAB): list of abscissae between X1 and X2.
C     type: if 'R' or 'r', central node taken to the right of abscissa whenever
C     possible; otherwise taken to its left.
C     Output: TAB(NTAB),
C     interpolation of the list in REFTAB at the NTAB points XABS.
      
 
      integer i, ita, itab, itab1, itab2, is
      
      double precision tastei, fi, p, pp1, pm1, a1, a2, a3

        if(nref.lt.3)then
        write(6,*)
     ;  'error interp1: cannot use quadratic interpolation, nref =', nref
        return
        end if
      
      tastei = dfloat(nref - 1) / (x2 - x1)
      is = 0
      if(type.ne.'r' .and. type.ne.'R')is = -1

         do 1 i = 1, ntab
         fi = (xabs(i) - x1) * tastei
         ita = min0(2 + is + int(fi), nref-1)
         ita = max0(ita, 2)
         p = fi - dfloat(ita - 1)
         itab = (ita - 1) * incref + 1
         itab1 = itab - incref
         itab2 = itab + incref
         pp1 = p + 1.
         pm1 = p - 1.
         a1 = p * pm1 * 0.5
         a2 = - pm1 * pp1
         a3 = p * pp1 * 0.5
         tab(i) =
     ;   a1 * reftab(itab1) + a2 * reftab(itab) + a3 * reftab(itab2)
  1      continue

      return
      end
