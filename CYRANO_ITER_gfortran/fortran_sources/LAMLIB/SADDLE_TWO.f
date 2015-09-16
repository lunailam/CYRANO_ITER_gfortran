

! **************************************************************** !
!                                                                  !
! SUBROUTINE SADDLE_TWO ( A(IN), CCosXo(IN), Ell(IN), ROOT(OUT) )  !
!                                                                  !
!                                                                  ! 
! **************************************************************** !
                                               
      SUBROUTINE SADDLE_TWO (AA, CCosXo, Ell, PROC, ROOT) 

	IMPLICIT NONE
! Input	                                                                 
	real*8, intent(in) :: AA, CCosXo  ! Orbit parameters
	integer, intent(in) :: Ell		  ! X - harmonic
	character(3), intent(in) :: PROC  ! solver: 'ern' or 'cub'

! Output
	real*8, intent(out) :: ROOT

! Parameters

	integer, parameter :: NN = 5000   ! Number of points used to find
									   ! saddlepoints  
	real*8 :: y(NN),		  ! imaginary axis for Fsad
     ;          Fsad(NN)		  ! Fsad = Im( d(expoent)/dz ) = 0 
	real*8 :: a, b, step, mu, coef(4), cubroot(3)
	complex*16 :: Xo		                                     
      integer :: j, rootstat, number 

! Initial values
      rootstat = 0
	ROOT = 0.0d0
	call cplxacos(CCosXo, Xo)




	if(PROC .eq. 'ern') then ! ------------------------------------

! Defining the grid to find the roots

	step = 10.0d0/dfloat(NN)
	
c	if (AA>301 .and. abs(CCosXo)<1) then
c		step=0.01/real(NN)
c	print *, 'A=',AA,'-->Restricted root search (y_sad -> 0)'
c	end if

! First element : Fsad(0)

	j = 1
c	y(j) = 0.0d0
	y(j) = abs(Xo)/2.0d0 
	a=dcosh(y(j))
	b=dsinh(y(j))

      Fsad(j) = - a*b + 2*b*a * (CCosXo*a/(2*a**2-1))**2 
     ;          - CCosXo**2*b*a/(2*a**2-1) + dfloat(Ell)/(2*AA*AA)

! Other elements : Fsad(2:NN)

	do j = 2, NN   ! beginning j-loop -----------

         y(j) = (dfloat(j)-1) * step  + abs(Xo)/2.0d0 
	   a=dcosh(y(j))
	   b=dsinh(y(j))

         Fsad(j) = - a*b + 2*b*a * (CCosXo*a/(2*a**2-1))**2 
     ;             - CCosXo**2*b*a/(2*a**2-1) + dfloat(Ell)/(2*AA*AA)
     	 
	   if (Fsad(j-1) * Fsad(j) < 0) then                             
	 	  ROOT = (y(j)+ y(j-1)) / 2  
		  rootstat = 1                                                    
	   end if                                                    
         if (rootstat .eq. 1)  EXIT   
	 
	end do ! end of j-loop ----------------------

	end if ! PROC = 'ern' --------------------------------------------


	if(PROC .eq. 'cub') then ! ---------------------------------------
	
	! Cubic equation

	  mu = dfloat(Ell)/(AA*AA)
	  coef(1) = -mu
	  coef(2) = 1.0d0-CCosXo**2
	  coef(3) = -mu
	  coef(4) = 1.0d0

	  call PA03A(coef,cubroot,number)
	  ROOT = 0.5d0 * dlog(cubroot(1)+dsqrt(1.0d0+cubroot(1)**2))	
	  if (number>0) rootstat = 1 

	end if ! PROC = 'cub' --------------------------------------------


	if (rootstat <= 0) then
	  print *, "SAD2 : No complex root for A =", AA, ", ell =", Ell
	  write(7,*) "SAD2 : No complex root for A =", AA, ", ell =", Ell
        ROOT = 0.0d0
      endif

	END SUBROUTINE SADDLE_TWO