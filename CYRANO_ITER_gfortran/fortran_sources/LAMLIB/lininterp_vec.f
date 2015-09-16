     
      SUBROUTINE lininterp_vec(x1, x2, y1_vec, y2_vec, N, xvalue, yvalues)

      IMPLICIT NONE

! Input

	integer :: N	                                                                 
	real*8, intent(in)  :: x1, x2	                 ! lower and upper limits 
	complex*16, intent(in)  :: y1_vec(N), y2_vec(N)  ! yi = f(xi)
	real*8, intent(in)  :: xvalue                    ! value where f(x) has 
	                                                 ! to be evaluated
! Output	                                                                 

	complex*16, intent(out)  :: yvalues(N)			 ! f(x_value) 

! Polynomial coefficients (y = ax + b)

	complex*16 :: acoef(N),bcoef(N)
	real*8 :: alpha, beta

! Interpolation

c	acoef=( y2_vec - y1_vec ) / (x2 - x1)
c	bcoef= y1_vec - acoef * x1
c	yvalues = acoef * xvalue + bcoef

	  alpha = (xvalue-x2) / (x1-x2)
	  beta =  (xvalue-x1) / (x2-x1)
c	beta = 1.0d0 - alpha
	yvalues = alpha*y1_vec + beta*y2_vec


      END SUBROUTINE lininterp_vec
