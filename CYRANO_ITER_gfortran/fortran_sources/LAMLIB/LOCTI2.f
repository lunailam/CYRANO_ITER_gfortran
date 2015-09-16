      subroutine locti2

      implicit none
      
c     Computes matrix 'bc2i' at local abscissa y
c     Trigono. data of magnetic field angle provided
c     by TABLES. Index in table is input through COMSWE.

      include 'pardim.copy'
      include 'comma2.copy'
      include 'commag.copy'
      include 'comro2.copy'
      include 'comfin.copy'
      include 'comin2.copy'
      include 'comswe.copy'
      include 'comphy.copy'

      integer i, j
      double precision t
      complex*16 ict,ist,ic1t,is1t
c     RESULTS ARE SENT THROUGH COMRO2

      t = sqrt2i
      si = sitab(intab)
      co = cotab(intab)
      si1 = si1tab(intab)
      co1 = co1tab(intab)

      ist = dcmplx(0.d0,si*t)
      ict = dcmplx(0.d0,co*t)
      is1t = dcmplx(0.d0,si1*t)
      ic1t = dcmplx(0.d0,co1*t)

      call zset( 6*5, czero, bc2i, 1 )

      bc2i(1,1) = t
      bc2i(1,2) = ict
      bc2i(1,4) = - ist

cPL 21/1/2005 what about contribution from Erho' to E+' and E-'? 
c     Ignored here, must be found from quadratic basis function derivatives for M23 elements, wherever needed!
      bc2i(2,2) = ic1t
      bc2i(2,3) = ict
      bc2i(2,4) = - is1t
      bc2i(2,5) = - ist

      bc2i(3,1) = t
      bc2i(3,2) = - ict
      bc2i(3,4) = ist

      bc2i(4,2) = - ic1t
      bc2i(4,3) = - ict
      bc2i(4,4) = is1t
      bc2i(4,5) = ist

      bc2i(5,2) = dcmplx(si,0.d0)
      bc2i(5,4) = dcmplx(co,0.d0)

      bc2i(6,2) = dcmplx(si1,0.d0)
      bc2i(6,3) = dcmplx(si,0.d0)
      bc2i(6,4) = dcmplx(co1,0.d0)
      bc2i(6,5) = dcmplx(co,0.d0)

      do 2 i = 1, 6
      do 2 j = 1, 5
      bc2ih(j,i) = dconjg( bc2i(i,j) )
    2 continue

      return
      end
