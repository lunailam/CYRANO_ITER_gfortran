      subroutine hunt(xx, n, x, jlo)

c     Numerical Recipes 3.4 p.111
c     Search an ordered table, starting from guess position jlo

      implicit none

	integer jlo, n
	double precision x, xx(n)
c
c     Input: xx(1:n), value x, jlo initial guess of index
c     Output: index jlo such that x is between xx(jlo) and xx(jlo+1)
c     xx must be monotonic increasing or decreasing
c     jlo=0 or jlo=n is returned to indicate x out of range
c
      integer inc, jhi, jm
	logical ascnd

	ascnd = xx(n) .ge. xx(1)
	  if(jlo.le.0 .or. jlo.gt.n)then
	  jlo = 0
	  jhi = n + 1
	  goto 3
	  end if
	inc = 1

        if(x.ge.xx(jlo) .eqv. ascnd)then
  1     jhi = jlo + inc
          if(jhi.gt.n)then
	    jhi = n + 1
	    else if(x.ge.xx(jhi) .eqv. ascnd)then
          jlo = jhi
	    inc = inc + inc
	    goto 1
	    end if
	  else
	  jhi = jlo
  2     jlo = jhi - inc
          if(jlo.lt.1)then
	    jlo = 0
	    else if(x.lt.xx(jlo) .eqv. ascnd)then
	    jhi  = jlo
	    inc = inc + inc
	    goto 2
	    end if
	  end if
  3     if(jhi - jlo .eq. 1)then
        if(x.eq.xx(n))jlo = n - 1
	  if(x.eq.xx(1))jlo = 1
	  return
	  end if
	jm = (jhi + jlo) / 2
	  if(x.ge.xx(jm) .eqv. ascnd)then
	  jlo = jm
	  else
	  jhi = jm
	  end if
	goto 3

	end


c NB after hunt or locate:
c k = min(max(j-(m-1)/2,1),n+1-m) where m is number of points in intepolation 
c then call polint(xx(k), yy(k), m, ...)