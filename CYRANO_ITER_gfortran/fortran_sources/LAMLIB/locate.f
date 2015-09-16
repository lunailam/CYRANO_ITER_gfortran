      subroutine locate(xx, n, x, j)

c     Numerical Recipes 3.4 p.111
c     Search an ordered table by bisection.

      implicit none

	integer j, n
	double precision x, xx(n)
c
c     Input: xx(1:n), value x
c     Output: index j such that x is between xx(j) and xx(j+1)
c     xx must be monotonic increasing or decreasing
c     j=0 or j=n is returned to indicate x out of range
c
      integer jl, jm, ju

	jl = 0
	ju = n + 1

  10    if(ju-jl.gt.1)then
        jm = (ju + jl) / 2
	   if((xx(n).ge.xx(1)) .eqv. (x.ge.xx(jm)))then
	   jl = jm
	   else
	   ju = jm
	   end if
	  goto 10
	  end if

	  if(x.eq.xx(1))then
	  j = 1
	  else if(x.eq.xx(n))then
	  j = n - 1
	  else
	  j = jl
	  end if

	return
	end