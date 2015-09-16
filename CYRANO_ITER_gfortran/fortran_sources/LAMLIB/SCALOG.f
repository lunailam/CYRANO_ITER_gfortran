      subroutine scalog(n, a, ipath)

      implicit none
      integer n, ipath
      double precision a(n)

C     Transforms the input real vector A(N) for approximate log plot
C     when A(I)>=<0.

      integer i
      double precision absa
      
      if(ipath.eq.1)then
      do 3 i = 1, n
      absa = abs(a(i))
      	if(absa.gt.1.)then
c      		if(a(i).gt.0.d0)then
c      		a(i) = dlog10(a(i))
c      		else
c      		a(i) = - dlog10(-a(i))
c      		end if
	    a(i) = dsign(dlog10(absa), a(i))
      	else
      	a(i) = 0.d0
      	end if
  3   continue

      else if(ipath.eq.2)then
C     Transfo. is B:= ( A + sqrt(A**2+4) ) / 2
      do 1 i = 1, n
      	if( a(i).lt.-1.d2 )then
C	     ...asymptotic limit
      	a(i) = - dlog10(-a(i))
      	else if( a(i).gt. 1.d2 )then
      	a(i) = dlog10(a(i))
      	else
      	a(i) = dlog10( (a(i)+dsqrt(a(i)**2+4.d0))*0.5d0 )
      	end if
   1  continue
      end if

      return
      end
