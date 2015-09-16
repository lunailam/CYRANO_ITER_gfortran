      SUBROUTINE INT_TO_STRING4 (number, string)

      IMPLICIT NONE

	! Input	                                                                 
	integer, intent(in)  :: number

	! Output	                                                                 
	character(4), intent(out)  :: string

c	print *, 'number=', number

	if (number < 10) then
		string = char(48) // char(48) // char(48) // char(number + 48)
c          string = char(number + 48)
	end if
		
	if (number >=10 .and. number<100) then 
		string = char(48) // char(48) // char(number/10 + 48) // 
     ;             char(number-10 * int(number/10) + 48)
c		string = char(number/10 + 48) // 
c     ;             char(number-10 * int(number/10) + 48)
	end if

	if (number >=100 .and. number<1000) then 
		string = char(48) // char(number/100 + 48) // 
     ;             char((number-100 * int(number/100))/10 + 48) //
     ;             char( number-10 * int(number/10) + 48)
	end if

	if (number >=1000 .and. number<10000) then 
		string = char(number/1000 + 48) // 
     ;			 char((number-1000 * int(number/1000))/100 + 48) // 
     ;             char((number-100 * int(number/100))/10 + 48) //
     ;             char( number-10 * int(number/10) + 48)
	end if

	if (number >= 10000) then
	string = '****'
	end if

	
c	print *, 'string=', string

	end SUBROUTINE INT_TO_STRING4
