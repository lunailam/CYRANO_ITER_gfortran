     
      SUBROUTINE find_datafile(cosXoin, file1)

      IMPLICIT NONE

      include 'pardim.copy'
	include 'comswe.copy'  
	include 'comgeo.copy' 

! Input
     	real*8, intent(in)  :: cosXoin	                 

! Output	                                                                 
	character(100) :: file1 

! Variables
	character(3) :: rrr
	character(100) :: fileaux(500)    
      real*8 :: cosXo(500)
	integer i, OpenStat
		
! Main program

	call int_to_string2( nint(100*abscis(intab)) , rrr )

! Initial values
c	i = 1
c	cosXo(i) = 0.0d0
c	fileaux(i) = '../../M12tables/r' // folder // '/m12_cosXo0p00.dat'

! Opening index_cos.dat
	i=0

	open (UNIT = 12, FILE = '../../M12tables/r' // rrr // '/index_cos.dat',
     ;      IOSTAT = OpenStat, STATUS = "OLD", ACTION = "READ")
	      if (OpenStat > 0) then
		     print *, 'Error reading file: ../../M12tables/r' // rrr 
     ;                   // '/index_cos.dat'
	         stop
	      end if
	do
	  i=i+1
        read ( 12, "(F10.5, G80.10)", END = 777), cosXo(i), fileaux(i)
	  if ( abs(cosXo(i)-cosXoin)<1d-3 ) then
	    file1 = fileaux(i)
	    goto 777
	  end if
	
	end do

777	continue


	close(12)

      END SUBROUTINE find_datafile
