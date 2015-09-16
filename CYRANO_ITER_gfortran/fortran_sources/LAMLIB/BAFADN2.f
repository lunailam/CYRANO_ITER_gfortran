      subroutine bafadn(it, x, baf, bafsin)

      implicit none

      integer it

      double precision x
c
c     For basis function type # it (=1 or 2),
c     Computes basis functions and derivatives in local reduced coordinate X varying on (0,1) -> BAF;
c     Computes (BAF(X)-BAF(0))/X, (BAF'(X)-BAF'(0))/X -> BAFSIN
c
      include 'pardim.copy'

      double precision baf(maxbaf,0:2), bafsin(maxbaf,0:1)
c     The second index of BAF and BAFSIN labels the derivatives.

      include 'comfin.copy'
c      
c     Basis function library: ntyp=2 available types are
c     (1) Hermite cubic, (2) Lagrange quadratic
c      
      integer i, j, ini, iu

      double precision
     ;  coe(maxdeg+1,maxbaf,ntyp), b(maxbaf), b1(maxbaf), b2(maxbaf), t
      
      data ini/0/
      save ini, coe
      
C     Coefficients of basis functions, for each type:
C     maximum of MAXBAF basis functions for each basis function set type, with maximum
C     degree MAXDEG
C     NTYP types of basis functions, NBAF basis functions for each type;
C     NELTYP types of finite element
      data
     ;  coe/ 1.d0, 0.d0,-3.d0, 2.d0
     ;     , 0.d0, 1.d0,-2.d0, 1.d0
     ;     , 0.d0, 0.d0, 3.d0,-2.d0
     ;     , 0.d0, 0.d0,-1.d0, 1.d0
     ;
     ;     , 1.d0,-3.d0, 2.d0, 0.d0
     ;     , 0.d0, 4.d0,-4.d0, 0.d0
     ;     , 0.d0,-1.d0, 2.d0, 0.d0
     ;     , 4*0.d0
     ;  /

      if(ini.eq.0)then
c     Determine the leading power of x in each basis function, for all types:
      do iu = 1, ntyp
        do j = 1, nbaf(iu)
        i = 1
  1     continue
          if(coe(i,j,iu).eq.0)then
            if(i.lt.deg(iu)+1)then
            i = i + 1
            go to 1
            else
            write(6,*)'Error BAFADN: a basis function is identically zero'
            stop
            end if
          end if
        leadeg(j,iu) = i - 1
        end do      
      end do
      ini = 1
      end if
      
c     Compute coefficients of basis functions and 2 derivatives for all b.f.
c     of current type IT: 
c     (NB: update the size of BAF if additional elt. types are added which
c     require higher derivatives.)

        do j = 1, nbaf(it)
        t = coe(deg(it)+1,j,it)
        b(j) = t
        b1(j) = deg(it) * t
        b2(j) = deg(it) * (deg(it)-1) * t
          do i = 1, deg(it)
          b(j) = b(j) * x + coe(deg(it)+1-i, j, it)
          end do
          do i = 1, deg(it)-1
          b1(j) = b1(j) * x + (deg(it)-i) * coe(deg(it)+1-i, j, it)
          end do
          do i = 1, deg(it)-2
          b2(j) = b2(j) * x + 
     ;           (deg(it)-i) * (deg(it)-i-1) * coe(deg(it)+1-i, j, it)
          end do

        baf(j,0) = b(j)
        baf(j,1) = b1(j)
        baf(j,2) = b2(j)
        end do
      
c     Compute coefficients of (basis function  minus value at 0) / x and
c             coefficients of (basis function' minus value at 0) / x.
c     (NB: update the size of BAFSIN if additional types are added which
c     require higher derivatives.)

        do j = 1, nbaf(it)
        t = coe(deg(it)+1,j,it)
        b(j) = t
        b1(j) = deg(it) * t
          do i = 1, deg(it)-1
          b(j) = b(j) * x + coe(deg(it)+1-i, j, it)
          end do
          do i = 1, deg(it)-2
          b1(j) = b1(j) * x + (deg(it)-i) * coe(deg(it)+1-i, j, it)
          end do
        bafsin(j,0) = b(j)
        bafsin(j,1) = b1(j)
        end do
      
      return
      end
