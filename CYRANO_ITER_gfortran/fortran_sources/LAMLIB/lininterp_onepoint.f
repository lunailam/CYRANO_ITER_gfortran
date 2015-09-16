     
      SUBROUTINE lininterp_onepoint(N, x, y, xvalue, yvalue)

      IMPLICIT NONE

c	Linear interpolation of function y(x) at point = xvalue
c	xvalue must lie between x(0) and x(N)
c	Result exported in yvalue(xvalue)

! Input

      integer, intent(in) :: N                                                            
	real*8, intent(in)  :: x(N), y(N)	         ! lower and upper limits 
	real*8, intent(in)  :: xvalue                ! value where f(x) has 
	                                              ! to be evaluated
! Output	                                                                 

	complex*16, intent(out)  :: yvalue			 ! f(x_value) 

! Polynomial coefficients (y = ax + b)
	integer :: j, indx_1, indx_2
c	complex*16 :: acoef,bcoef
	real*8 :: alpha, beta, x1, x2, y1, y2

! Interpolation

	if (xvalue .lt. x(1) .or. xvalue .gt. x(N)) then
		print *, 'interp_onepoint problem : value out of range!'
	    print *, 'xvalue=', xvalue 
		yvalue = 0.0d0
	else

c	First find (x1,x2,y1,y2) values of the given function y(x), defining
c	the interval of the linear interpolation (nearest to xvalue)
	indx_1 = 0
	do j = 1, N
	   if (x(j) .eq. xvalue) then  ! exact match (no interpolation)
		    yvalue = y(j)
			goto 1000
	   end if

	   if (x(j) .lt. xvalue) then
		   indx_1 = indx_1 + 1	   ! index of x1 (nearest value < xvalue)
	   end if
	   	
	end do 

	indx_2 = indx_1 + 1			   ! index of x2 (nearest value > xvalue)

	    x1 = x(indx_1)
		y1 = y(indx_1)
		x2 = x(indx_2)
		y2 = y(indx_2)

	   alpha = (xvalue-x2) / (x1-x2)
	   beta =  (xvalue-x1) / (x2-x1)
c	   beta = 1.0d0 - alpha
	   yvalue = alpha*y1 + beta*y2
	
	end if

1000	continue

	return

      END SUBROUTINE lininterp_onepoint
