! *************************************************************** !
!                                                                 !
! SUBROUTINE read_prof (profile, specie, x, file_result, sortie )  !
!                                                                 !
! *************************************************************** !
                                               
      SUBROUTINE read_prof(profile, specie, x, file_result, sortie) 

      IMPLICIT NONE

! Input	  
	character(*), intent(in) :: profile
	character(4), intent(in) :: specie	                                                               
	real*8, intent(in)  :: x
   
! Output
	integer, intent(out) :: file_result		                                                                 
	real*8, intent(out) :: sortie


! Parameters

	integer :: j, OpenStat, filesize
	character(50) :: filename
	real*8 :: xall(100), nall(100)

! Main program

	filename = '../../Profiles/' // profile // '_' // specie // '.dat'

	file_result = 0
	j = 1
	
	open(UNIT = 99, FILE = filename, STATUS = "OLD",
     ;     IOSTAT = OpenStat, ACTION = "READ")
     	      
		  if (OpenStat > 0) then
c		      print *,'Error reading file:', filename
c	          print *,'-> Using standard value instead!'
			  file_result = 0
	          sortie = 0.0d0
	          close(99)
	          go to 333
	      else
	          read ( 99, *), xall(1), nall(1)
			  if ( x .eq. xall(1) ) then
			       sortie = nall(1)
				   goto 222
			  end if
		      
			  do 
	             j = j + 1
	      	     read ( 99, *, END = 222), xall(j), nall(j)
				 if ( x .eq. xall(j) ) then
				      sortie = nall(j)
					  goto 222
	             end if
				 if ( x < xall(j) .and. x > xall(j-1) ) then
c					 sortie = (nall(j)+nall(j-1)) / 2
                       call lininterp_vec( xall(j-1), xall(j), 
     ;                      nall(j-1), nall(j), 1, x, sortie)
				     goto 222
				 end if

	          end do
222			  continue
			  file_result = 1
			  filesize = j - 1
		  end if

	close(99)


333	continue


      END SUBROUTINE read_prof

! *********************************************************************** !
