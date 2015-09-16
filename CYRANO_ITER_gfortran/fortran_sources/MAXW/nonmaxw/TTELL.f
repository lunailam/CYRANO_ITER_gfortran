
! *************************************************************** !
!                                                                 !
! SUBROUTINE TT_ell ( AA, CCosXo, Ell, FLAG, z_sad, TTELL )       !
!                                                                 !
! External function to calculate the coeficients T_ell using      !
! the saddlepoint / steepest descent method                       !
!                                                                 ! 
! *************************************************************** !
                                               
      SUBROUTINE TT_ell (AA, CCosXo, Ell, FLAG, z_sad, TTELL) 
	
	IMPLICIT NONE

! Input	                                                                 
	real*8, intent(in)  :: AA, CCosXo   ! Orbit parameters
	integer, intent(in) :: Ell			! X - harmonic
	integer, intent(in) :: FLAG			! Flag for root search
! Output
	complex*16, intent(out) :: z_sad,   ! Saddlepoint
     &	                       TTELL    ! T_ell coeficient
							   
! Constants
      real*8,     parameter :: pi = 3.14159265d0
	complex*16, parameter :: myi = (0.0d0,1.0d0)

! Parameters
	
	real*8 :: x_sad, y_sad         ! real and imag. roots of saddlepoint eq.
	complex*16 ::  tf, tf2, xaux   ! - auxiliary functions    
	real*8 :: auxXo,step 	                                      
      integer :: j, rootstat, sinal


! Calculates the saddlepoints : 1st root branch (cosXo>1) ->SADDLE_ONE)	

	if (FLAG .eq. 1) then
			
		sinal = -sign(1.0d0,CCosXo)
		call SADDLE_ONE (AA, CCosXo, Ell, y_sad)
		x_sad = 0.0d0

	end if

! Calculates the saddlepoints : 2nd root branch (cosXo>1) -> SADDLE_TWO

	if (FLAG .eq. 2) then

		sinal = -sign(1.0d0,CCosXo)
		call SADDLE_TWO (AA, CCosXo, Ell, 'cub', y_sad)
c		x_sad = dacos (CCosXo*dcosh(y_sad)/(1+2*dsinh(y_sad)**2))
		auxXo = CCosXo*dcosh(y_sad)/(1.0d0+2.0d0*dsinh(y_sad)**2)
		call cplxacos(auxXo, xaux)
	    x_sad = dreal(xaux)
	
	end if

! Calculates the saddlepoints: 2nd root branch (cosXo<=1) -> SADDLE_THREE

	if (FLAG .eq. 3) then

		sinal = sign(1.0d0,CCosXo)
		call SADDLE_THREE (AA, CCosXo, Ell, 'cub', y_sad)
c	    print*, y_sad
c		x_sad = dacos (CCosXo*dcosh(y_sad)/(1+2*dsinh(y_sad)**2))
		auxXo = CCosXo*dcosh(y_sad)/(1.0d0+2.0d0*dsinh(y_sad)**2)
		call cplxacos(auxXo, xaux)
	    x_sad=dreal(xaux)
	
	end if


! Complex root	
	 z_sad = x_sad + myi * y_sad	! complex root (saddlepoints)

! Calculates the T_ell coeficients (steepest descent method)

      tf = AA**2 * (cdcos(z_sad) - CCosXo)**2 - myi * dfloat(Ell) * z_sad
c	if (abs(tf)>100) then
c	tf = 100.0d0 + myi*0.0d0
c	print *,'WARNING: T_ell underflow artificially avoided!' 
c	end if

	tf2 = -2.0d0 * ( cdcos(2.0d0*z_sad) - cdcos(z_sad) * CCosXo) 
       
	TTELL = sign(1.0d0,1.0d0-dabs(CCosXo)) *myi * pi * 
     ; cdsqrt(2.0d0/AA**2/tf2) * cdexp(-tf) * dfloat(sinal)  ! T_ell coeficient

c	z_sad=cdexp(-tf)
	END SUBROUTINE TT_ell                                       
                                                                  
! *************************************************************** !
