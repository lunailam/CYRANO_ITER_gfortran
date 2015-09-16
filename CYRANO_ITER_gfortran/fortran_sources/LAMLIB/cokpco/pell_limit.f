     
      SUBROUTINE pell_limit(m_indx, L_min, L_max)

      IMPLICIT NONE

      include 'pardim.copy'
	include 'comgeo.copy' 

! Input
     	integer, intent(in) :: m_indx	                 

! Output	                                                                 
	integer, intent(out) :: L_min, L_max 

! Variables
	integer k, k2, MAXP
	real*8 :: biggest
! Main program
	MAXP = 128

	L_min=0
	L_max=MAXP
	biggest = maxval(real(gcdr(1:256,m_indx)))

	   do k = 0, MAXP
	      if (dreal(gcdr(k,m_indx)) > biggest*1.d-3) then
	         L_min = k
	         goto 666
		  end if
	   end do

666	   continue
		
	   do k = 0, MAXP
		  k2=MAXP-k
	      if (dreal(gcdr(k2,m_indx)) > biggest*1.d-3) then
	         L_max = k2
	         goto 777
	      end if
	   end do
			
777	   continue	


      END SUBROUTINE pell_limit