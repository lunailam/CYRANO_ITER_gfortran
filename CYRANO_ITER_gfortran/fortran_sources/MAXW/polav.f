      double precision function polav(vec, incvec, n, polsym)

	double precision vec(1+(n-1)*incvec)

	integer n

	logical polsym

c     Poloidal average of n elements of array vec, stored with increment incvec assumed > 0
c     When polsym = .true., computes (vec(1) + 2*sum(vec(2:n-1)) + vec(n)) / (2*(n-1))
c     When polsym = .false., computes sum(vec(1:n-1)) / (n-1)

      integer ilast, ibutlast, i, iv

        if(incvec .eq. 1)then
	  ilast = n
	  else
	  ilast = 1+(n-1)*incvec
	  end if
      ibutlast = ilast - incvec

	polav = 0.d0
        if(polsym)then
	  iv = 1
	    do i = 2, n-1
	    iv = iv + incvec
	    polav = polav + vec(iv)
	    end do
	  polav = (2.d0 * polav + vec(1) + vec(ilast)) / dfloat(2*(n-1))
	  else
	  iv = 1
	    do i = 1, n-1
	    polav = polav + vec(iv)
	    iv = iv + incvec
	    end do
	  polav = polav / dfloat(n-1)
	  end if

	return

	end