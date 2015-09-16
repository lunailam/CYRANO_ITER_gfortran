

! *************************************************************** !
!                                                                 !
! SUBROUTINE PP_ell ( AA(IN), CCosXo(IN), NNchi(IN), TELL(OUT) )  !
!                                                                 !
! External function to calculate the coeficients P_ell using      !
! the rectangles method to Approximate the integral               !
!                                                                 ! 
! *************************************************************** !
                                               
      SUBROUTINE PP_ell (MMbar, EEll, Eps, Bmean, PPELL) 

! Input	                                                                 
	integer, intent(in) :: MMbar, EEll  ! - Chi - harmonic
	real*8, intent(in) :: Eps, Bmean    ! - Orbit parameters

! Output
	complex*16, intent(out) :: PPELL  ! - T_ell coeficient
							   
! Constants
      real,    parameter :: pi = 3.14159265
	complex, parameter :: myi = (0,1)

	integer, parameter :: NNpi = 1000  ! - Number of points used to find
									   !   saddlepoints between (0 - pi)
	real*8, parameter :: step = pi / NNpi 

	real*8, dimension (NNpi + 1) :: theta, Chi, Bo  ! imaginary axis for FFsaddle
	complex*16, dimension (NNpi + 1) :: FFF      ! FFsaddle = d(expoent)/dz = 0 
								              ! is the saddlepoint equation                                  
      integer :: j
    

! Loop to calculate the expression FFsaddle(z) = d(expoent)/dz --------------	
	
	do j = 1, NNpi ! beginning j-loop 

       theta(j) = (real(j)-1) * step  

       Chi(j) = acos ( (cos(theta(j)) + Eps) / 
     &               ( 1 + Eps * cos(theta(j)) ) )
       Bo(j)  = Bmean / ( 1 + Eps * cos(theta(j)) )
       !Bo(j)  = Bmean
       FFF(j) = cos ( MMbar*theta(j) ) * cos ( EEll*Chi(j) ) / Bo(j) 

	if (EEll .eq. 0) then
	FFF(j) = 0.5* cos ( MMbar*theta(j) ) * cos ( EEll*Chi(j) ) / Bo(j)
	end if

	end do ! end of j-loop  	  

	     
      PPELL = 2/pi * (sum(FFF) - 0.5 * FFF(1)) * step
	
	!print *, PPELL, PPELL2	
	                                                              
      END SUBROUTINE PP_ell                                       
                                                                  
! *************************************************************** !

