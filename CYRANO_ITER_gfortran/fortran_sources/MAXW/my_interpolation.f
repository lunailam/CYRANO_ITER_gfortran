     
      SUBROUTINE my_interpolation(Ndat, xdata, Fdata, Nval, xvalue, Fvalue)

c		Given the (xdata,Fdata) data values, compute the Fvalue of the
c         function F at the intermediate xvalue point
c
c

      IMPLICIT NONE

! Input

	integer :: Ndat, Nval                ! Number of I/O points                                                               
	real*8, intent(in)  :: xdata(Ndat), Fdata(Ndat),  ! Given data points
     ;                       xvalue              ! x-values to evaluate the data 

! Output	                                                                 

	real*8, intent(out)  :: Fvalue			 ! F(x_values) 


! Polynomial coefficients (y = ax + b)

	complex*16 :: acoef,bcoef
	real*8 :: alpha, beta, x1, x2
	integer :: ix1, isrchfge


	if(xvalue.le.xdata(1))then  ! use first value
		
		Fvalue = Fdata(1)
		return

	elseif(xvalue.ge.xdata(Ndat))then  ! use last value

		Fvalue = Fdata(Ndat)	
		return

	else   ! xvalue inside xdata points
c	-----------------------------------


c		Find closest data points to given xvalue (ix1 >= xvalue > ix1-1)
		ix1 = isrchfge(Ndat, xdata, 1, xvalue)

	       if(xdata(ix1).eq.xvalue)then   ! EXACT match
	    
		       Fvalue = Fdata(ix1)   
			   return
			     
	       else
		   
		       x1 = xdata(ix1-1)
			   x2 = xdata(ix1)   
	           alpha = (xvalue-x2) / (x1-x2)
	           beta =  (xvalue-x1) / (x2-x1)

	           Fvalue = alpha*Fdata(ix1-1) + beta*Fdata(ix1)

		   end if

	end if
c	------------------------------------


	return

      END SUBROUTINE my_interpolation
