     
      SUBROUTINE lininterp_array(x1, x2, y1_vec, y2_vec, M, N, xvalue, yvalues)

      IMPLICIT NONE

! Input

	integer :: N,M	                                                                 
	real*8, intent(in)  :: x1, x2	                 ! lower and upper limits 
	complex*16, intent(in)  :: y1_vec(M,N), y2_vec(M,N)  ! yi = f(xi)
	real*8, intent(in)  :: xvalue                    ! value where f(x) has 
	                                                 ! to be evaluated
! Output	                                                                 

	complex*16, intent(out)  :: yvalues(M,N)			 ! f(x_value) 

! Polynomial coefficients (y = ax + b)

	complex*16 :: acoef(M,N),bcoef(M,N)

! Interpolation

	acoef=( y2_vec - y1_vec ) / (x2 - x1)

	bcoef= y1_vec - acoef * x1
	yvalues = acoef * xvalue + bcoef

      END SUBROUTINE lininterp_array
