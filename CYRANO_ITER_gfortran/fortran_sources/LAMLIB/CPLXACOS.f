

! *************************************************************** !
!                                                                 !
! SUBROUTINE CPLXACOS ( entra(IN), sai(OUT) )						!
!                                                                 !
! External function to calculate the coeficients P_ell using      !
! the rectangles method to Approximate the integral               !
!                                                                 ! 
! *************************************************************** !
                                               
      SUBROUTINE CPLXACOS (entra, sai) 

! Input	                                                                 
	real*8, intent(in) :: entra    ! - Orbit parameters

! Output
	complex*16, intent(out) :: sai		! - T_ell coeficient
							   
! Constants
      real,    parameter :: pi = 3.14159265d0
	complex, parameter :: myi = (0.0d0,1.0d0)

! Aproximations for computing acos(X) --------------	
	
		if (entra > 1) then
			sai = myi * log(entra - myi*myi * dsqrt(-1.0d0+entra**2))
c			print *,'cosXo>1=',entra,' -> Xo=',sai
		end if

		if (entra < -1) then

			sai = pi+myi * log(-entra - myi*myi * dsqrt(-1.0d0+entra**2))
c			print *,'cosXo<-1=',entra,' -> Xo=',sai
		end if

		if (entra <= 1 .and. entra >= -1) then
		    sai = dacos (entra)
c			print *, '|cosXo|<1=,',entra,' -> Xo=', sai
		end if

	                                                              
      END SUBROUTINE CPLXACOS                                       
                                                                  
! *************************************************************** !

