      subroutine polint(xa, ya, n, x, y, dy)

      implicit none
	integer n
	double precision dy, x, y, xa(n), ya(n)
c
c     n-point (i.e. n-1 th degree) polynomial interpolation.
c     Numerical Recipes 3.1 page 103
c     Input: arrays xa, ya of length n, and a value x
c     Output: value y, interpolation of ya by (n-1)th degree polynomial at x,
c             and error estimate dy
c
	integer nmax, i, m, ns
c     Largest anticipated value of n:
	parameter (nmax = 10)
      double precision den, dif, dift, ho, hp, w, c(nmax), d(nmax)

	ns = 1
	dif = dabs(x - xa(1))
c     Find index ns of closest table entry, and initialize c, d:
	  do i = 1, n
	  dift = dabs(x - xa(i))
	    if(dift .lt. dif)then
	    ns = i
	    dif = dift
	    end if
        c(i) = ya(i)
	  d(i) = ya(i)
	  end do
	y = ya(ns)
	ns = ns - 1
	  do m = 1, n-1
	    do i = 1, n-m
	    ho = xa(i) - x
	    hp = xa(i+m) - x
	    w = c(i+1) - d(i)
	    den = ho - hp
c     This error can only occur if two input xa's are identical within roundoff:
	    if(den.eq.0.)write(6,*)'failure in polint'
	    den = w / den
	    d(i) = hp * den
	    c(i) = ho * den
	    end do
c     After each column in the tableau is completed, we decide which correction,
c     c or d, we want to add to our accumulating value of y, i.e. which path
c     to take through the tableau - forking up or down. Take the most 'straight
c     line'route through the tableau to its apex, updating ns to keep track of
c     where we are. This route keeps the partial approximation centered (insofar
c     as possible) on the target x. The last dy added is thus the error indication.
	    if(2*ns .lt. n-m)then
	    dy = c(ns+1)
	    else
	    dy = d(ns)
	    ns = ns - 1
	    end if
	  y = y + dy
	  end do

	return
	end
