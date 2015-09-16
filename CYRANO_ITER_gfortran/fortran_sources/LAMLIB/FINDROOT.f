

! *************************************************************** !
!                                                                 !
! SUBROUTINE SUBROUTINE FINDROOT (FF, xx, NN, raiz)               !
!                                                                 !
! *************************************************************** !
                                               
      SUBROUTINE FINDROOT (FF, xx, NN, raiz) 

! Input	 
	integer, intent(in) :: NN		   ! Chi - harmonic                                                                
	real*8,  intent(in) :: FF(NN), xx(NN)  ! Orbit parameters

! Output
	real*8, intent(out) :: raiz

! Variables
	real*8  :: xmin, xmed, xmax
	integer :: Nmin, Nmed, Nmax
	integer :: j

! Loop to calculate the expression FFsaddle(z) = d(expoent)/dz --------------	
	
c	print *, xx

      Nmin=1;
      Nmax=NN;
      Nmed=int((Nmin+Nmax)/2);

	do j = 1,100

		if (FF(Nmin)*FF(Nmax)<0) then
			Nmed=int((Nmin+Nmax)/2)
			xmed=xx(Nmed)
   
   				if (FF(Nmin)*FF(Nmed)<0) then
					Nmax=Nmed;
					xmax=xx(Nmax)
					xmin=xx(Nmin)
				end if
				if (FF(Nmax)*FF(Nmed)<0) then
					Nmin=Nmed;
					xmin=xx(Nmin)
					xmax=xx(Nmax)
				end if

		end if

		raiz =  (xmax+xmin)/2 
c	print *, xmin, xmax, raiz
	
	end do

	


	END SUBROUTINE FINDROOT
	                                       
! *************************************************************** !                                                                
