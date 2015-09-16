      SUBROUTINE INT_TO_STRING (number, string)

      IMPLICIT NONE

	! Input	                                                                 
	integer, intent(in)  :: number

	! Output	                                                                 
	character(3), intent(out)  :: string

c	print *, 'number=', number

	if (number < 10) then
c		string = char(48) // char(48) // char(number + 48)
          string = char(number + 48)
c	end if
		
	else if (number >=10 .and. number<100) then 
c		string = char(48) // char(number/10 + 48) // 
c     ;             char(number-10 * int(number/10) + 48)
		string = char(number/10 + 48) // 
     ;             char(number-10 * int(number/10) + 48)
c	end if

	else if (number >=100 .and. number<1000) then 
		string = char(number/100 + 48) // 
     ;             char((number-100 * int(number/100))/10 + 48) //
     ;             char( number-10 * int(number/10) + 48)
c	end if

c	if (number >= 1000) then
	else
	string = '***'
	end if

	
c	print *, 'string=', string

	end SUBROUTINE INT_TO_STRING
