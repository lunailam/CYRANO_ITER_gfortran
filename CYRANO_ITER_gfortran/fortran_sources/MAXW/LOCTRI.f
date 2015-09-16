      subroutine loctri(iele, elem)

      implicit none
      logical elem
      integer iele

c     Computes matrix 'bci' of plasma term at Gauss points of element. 
c     6 by 6 block(s): rows are +,+',-,-',//,//'; columns are rho, rho', theta, theta', phi, phi'.
c     Trigono. data of magnetic field angle provided by tables. Index in table computed from INTAB,
C     INTAB input in COMSWE set to first Gauss point.
c     Circular cross-section
c     Results are sent through comrot
c     If ELEM is .false., calculation at current table index INTAB only.

      include 'pardim.copy'
      include 'comma2.copy'
      include 'commag.copy'
      include 'comrot.copy'
      include 'comfin.copy'
      include 'comin2.copy'
      include 'comswe.copy'
      include 'comphy.copy'

      integer ig1, ig2, i, j, intabl
      double precision t
      complex*16 ict,ist,ic1t,is1t

      t = sqrt2i

      call zset(6*6*(maxgau+1)*2, czero, bci1, 1)
c      call zset(360, czero, bci1, 1)

      if(elem)then
      ig1 = 0
      ig2 = ngauss
      else
      ig1 = 1
      ig2 = 1
      end if

      do 100 ig = ig1, ig2
      intabl = intab + ig - 1
      si = sitab(intabl)
      co = cotab(intabl)
      si1 = si1tab(intabl)
      co1 = co1tab(intabl)

      ist = dcmplx(0.d0,si*t)
      ict = dcmplx(0.d0,co*t)
      is1t = dcmplx(0.d0,si1*t)
      ic1t = dcmplx(0.d0,co1*t)

      bci(1,1,ig) = t
      bci(1,3,ig) = ict
      bci(1,5,ig) = - ist

      bci(2,2,ig) = t
      bci(2,3,ig) = ic1t
      bci(2,4,ig) = ict
      bci(2,5,ig) = - is1t
cPL 21/1/2005 bug:
c      bci(2,5,ig) = - ist
      bci(2,6,ig) = - ist

      bci(3,1,ig) = t
      bci(3,3,ig) = - ict
      bci(3,5,ig) = ist

      bci(4,2,ig) = t
      bci(4,3,ig) = - ic1t
      bci(4,4,ig) = - ict
      bci(4,5,ig) = is1t
      bci(4,6,ig) = ist

      bci(5,3,ig) = dcmplx(si,0.d0)
      bci(5,5,ig) = dcmplx(co,0.d0)

      bci(6,3,ig) = dcmplx(si1,0.d0)
      bci(6,4,ig) = dcmplx(si,0.d0)
      bci(6,5,ig) = dcmplx(co1,0.d0)
      bci(6,6,ig) = dcmplx(co,0.d0)

      do 1 i=1,6
      do 1 j=1,6
      bcih(i,j,ig) = dconjg(bci(j,i,ig))
  1   continue

 100  continue

      return
      end
