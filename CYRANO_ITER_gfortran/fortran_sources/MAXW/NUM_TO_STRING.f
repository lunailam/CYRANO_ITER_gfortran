      SUBROUTINE NUM_TO_STRING (number, string)

      IMPLICIT NONE

	! Input	                                                                 
	integer, intent(in)  :: number

	! Output	                                                                 
	character(20), intent(out)  :: string

	print *, 'number=', number

	if (number < 10) then
		string = char(number + 48)
	end if
		
	if (number >=10 .and. number<100) then 
		string = char(number/10 + 48) // 
     ;             char(number-10 * int(number/10) + 48)
	end if

	if (number >=100 .and. number<1000) then 
		string = char(number/100 + 48) // 
     ;             char((number-100 * int(number/100))/10 + 48) //
     ;             char( number-10 * int(number/10) + 48)
	end if

	
	print *, 'string=', string

	end SUBROUTINE NUM_TO_STRING
