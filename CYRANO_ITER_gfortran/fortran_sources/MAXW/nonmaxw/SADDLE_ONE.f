

! **************************************************************** !
!                                                                  !
! SUBROUTINE SADDLE_ONE ( A(IN), CCosXo(IN), Ell(IN), ROOT(OUT) )  !
!                                                                  !
!                                                                  ! 
! **************************************************************** !
                                               
      SUBROUTINE SADDLE_ONE (AA, CCosXo, Ell, ROOT) 
	
	IMPLICIT NONE

! Input	                                                                 
	real*8, intent(in) :: AA, CCosXo  ! Orbit parameters
	integer, intent(in) :: Ell		  ! X - harmonic

! Output
	real*8, intent(out) :: ROOT

! Parameters

	integer, parameter :: NN = 50000   ! Number of points used to find
 								         ! saddlepoints  

	real*8 :: y(NN),		  ! imaginary axis for Fsad
     ;          Fsad(NN)		  ! Fsad = Im ( d(expoent)/dz ) = 0 

	real*8 :: a, b, step
	complex*16 :: Xo		                                     
      integer :: j, rootstat 

! Initial values
      rootstat = 0
	ROOT = 0.0d0
	call cplxacos(CCosXo, Xo)

! Defining the grid to find the roots

	step = abs(Xo)/dfloat(NN)
	

! First element : Fsad(Xo)

	j = 1
	y(j)=abs(Xo)
	a=dcosh(y(j))
	b=dsinh(y(j))

	Fsad(j) = b * (a-CCosXo) + dfloat(Ell)/(2*AA*AA)

! Other elements : Fsad(2:NN)
! Note: The root search is now performed backwards (Xo -> 0)

	do j = 2, NN   ! beginning j-loop --------------------------

	   y(j) = abs(Xo)-(dfloat(j)-1) * step 
	   a=dcosh(y(j))
	   b=dsinh(y(j))

	   Fsad(j) = b * (a-CCosXo) + dfloat(Ell)/(2*AA*AA)
     	 
	   if (Fsad(j-1) * Fsad(j) < 0) then                             
	 	  ROOT = (y(j)+ y(j-1)) / 2.0d0  
		  rootstat = 1                                                    
	   end if                                                    
         if (rootstat .eq. 1)  EXIT   
	 
	end do ! end of j-loop ------------------------------------- 

cccc	  call FINDROOT(FFsaddle, y, NNpi, yy_saddle)

	if (rootstat <= 0) then
	print *,   "SAD1 : No complex root for A =", AA, ", ell =", Ell
	write(7,*) "SAD1 : No complex root for A =", AA, ", ell =", Ell
      ROOT = 0.0d0
	end if

	END SUBROUTINE SADDLE_ONE