     
      SUBROUTINE set_filename(ROOTDIR, cosXoin, FILENAME)

      IMPLICIT NONE

c     include 'pardim.copy'
c	include 'comswe.copy'  
c	include 'comgeo.copy' 

! Input
     	real*8, intent(in)  :: cosXoin	                 
	character(20), intent(in) :: ROOTDIR 
! Output	                                                                 
	character(150) :: FILENAME 

! Variables
	character :: auxchar0, auxchar1, auxchar2, auxchar3(3)	
	real*8 :: cosXo

! Main program

	cosXo = dabs(cosXoin)

	call int_to_string(int(cosXo), auxchar3 )

	  auxchar2 = char( int(cosXo) + 48 )
		
	  auxchar1 = char( int(10*cosXo)
     ;                   - 10 * int(cosXo) 
     ;                   + 48 )
	    
	  auxchar0 = char( int(100*cosXo)
     ;	               - 10 * int(10*cosXo) 
     ;				   + 48 )
	
	if (cosXoin >= 0.0d0) then

	  if (cosXo < 10.0d0) then
	      FILENAME = ROOTDIR // '/m12_cosXo'  // 
     ;                 auxchar2 // 'p' // auxchar1 // auxchar0 // '.dat'
	  end if

	  if (cosXo >= 10.0d0 .and. cosXo < 100.0d0) then
	      FILENAME = ROOTDIR // '/m12_cosXo'  // 
     ;                 auxchar3(1) // auxchar3(2) // 'p' // 
     ;                 auxchar1 // auxchar0 // '.dat'
	  end if
		
	  if (cosXo >= 100.0d0) then
	      FILENAME = ROOTDIR // '/m12_cosXo' // 
     ;                 auxchar3(1) // auxchar3(2) // '0p' // 
     ;                 auxchar1 // auxchar0 // '.dat'
	  end if	

	end if

!	case cosXo < 0
	if (cosXoin < 0.0d0) then

	  if (cosXo < 10.0d0) then
	      FILENAME = ROOTDIR // '/m12_cosXo'  // 
     ;                 auxchar2 // 'p' // auxchar1 // auxchar0 // 'neg.dat'
	  end if

	  if (cosXo >= 10.0d0 .and. cosXo < 100.0d0) then
	      FILENAME = ROOTDIR // '/m12_cosXo'  // 
     ;                 auxchar3(1) // auxchar3(2) // 'p' // 
     ;                 auxchar1 // auxchar0 // 'neg.dat'
	  end if
		
	  if (cosXo >= 100.0d0) then
	      FILENAME = ROOTDIR // '/m12_cosXo' // 
     ;                 auxchar3(1) // auxchar3(2) // '0p' // 
     ;                 auxchar1 // auxchar0 // 'neg.dat'
	  end if	

	end if



c	print *, cosXo, FILENAME


      END SUBROUTINE set_filename
