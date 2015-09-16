! *************************************************************** !
!                                                                 !
! SUBROUTINE derivate (Nx, F, x, dFdx)										!
!                                                                 !
! *************************************************************** !
                                               
      SUBROUTINE derivate (Nx, F, x, dFdx)	 

      IMPLICIT NONE


! Input	  
	integer, intent(in) :: Nx
	real*8,  intent(in) :: F(*), x(*)

! Output
	real*8,  intent(out) :: dFdx(*)

! Parameters

	integer :: j


! Main program


	dFdx(1) = (F(2)-F(1)) / (x(2)-x(1))	! First point

	do j = 2 , Nx - 1
	   dFdx(j) = (F(j+1)-F(j-1)) / (x(j+1)-x(j-1))
	end do

	dFdx(Nx) = (F(Nx)-F(Nx-1)) / (x(Nx)-x(Nx-1)) ! Last point


      END SUBROUTINE derivate

! *********************************************************************** !
