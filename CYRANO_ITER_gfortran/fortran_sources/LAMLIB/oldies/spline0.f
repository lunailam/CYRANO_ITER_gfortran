	subroutine spline0(Nf, x, F, Np, xp, Fp)

c	Subroutine to intepolate a given function F at xp points 
c	using the cubic splines method with IMSL routines
c	method : 'not a knot'
c
c		Nf (IN) : Number of points of the given function F(x) 
c		x  (IN) : Abcissa of given function					 
c		F  (IN) : Given function 
c		Np (IN) : Number of points to evaluate F
c		xp (IN) : Points (abcissa) to evaluate F
c		Fp (OUT): F(xp) - Interpolated values of F at the xp points

	IMPLICIT NONE

!	Variables 
	integer, intent(in) :: Nf, Np
	real*8,  intent(in) :: x(Nf), F(Nf), xp(Np)
	real*8, intent(out) :: Fp(Np)
	real*8 :: break(Nf), coef(4, Nf), DCSVAL
	integer k

	call DCSDEC(Nf, x, F, 0, 0, 0, 0, break, coef)
	
	do k = 1,Np
	   Fp(k) = DCSVAL(xp(k), Nf-1, break, coef)
	end do


	end subroutine spline0