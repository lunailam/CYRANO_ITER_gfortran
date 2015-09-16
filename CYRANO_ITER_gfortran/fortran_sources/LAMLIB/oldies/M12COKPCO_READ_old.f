! *************************************************************** !
!                                                                 !
! SUBROUTINE M12cokpco_read ( cosXo, mdiff, Nm, Avector, NA, sortie )   !
!                                                                 !
! *************************************************************** !
                                               
      SUBROUTINE M12cokpco_read (cosXoin, mdiff, Nm, Avector, NA, 
     ;           filename, sortie) 

	USE DFPORT		! necessary for 'etime' function

      IMPLICIT NONE

      include 'pardim.copy'
      include 'commod.copy'
      include 'cokpco.copy'
      include 'comswe.copy'      
      include 'comgeo.copy' 

! Input	                                                                 
	real*8, intent(in)  :: cosXoin
	integer, intent(in) :: Nm,NA
	integer, intent(in) :: mdiff(Nm)  
	real*8, intent(in)  :: Avector(NA)
	character(*), intent(in) :: filename		   
! Output	                                                                 
	complex*16, intent(out) :: sortie(NA,Nm)

! Constants
	complex*16, parameter :: myi = (0.0d0,1.0d0)

c	real, dimension (2) :: T1	! - Variables for evaluating
c	real :: tempo, tempo0               !   the elapsed CPU time
! Parameters

      real*8 :: A, sigma_A, cosXotest, cosXo, sigma_cos
	complex*16 :: data1(NA,Nm) , data2(NA,Nm) 	

	integer :: j, i	
	integer :: mdiffaux(2*Nm)
	character(100) :: GGs, GG2s, data_file1, data_file2 
	real*8 :: cosXo1, cosXo2
	character(2) :: charaux, charaux2
	integer :: filesize
	real*8 ::     Atab1(500)  !,Atab2(500)
	complex*16 :: m12tab1(500,Nm)  !,m12tab2(500,Nm)
	real*8 :: x1, x2
	complex*16 :: y1(Nm), y2(Nm)
	logical :: interpcos
	real*8 :: alpha, beta
! Main program

c	tempo0 = etime (T1)

		cosXo = abs(cosXoin)
	    sigma_cos   = sign(1.0d0, cosXoin)
c	    sigma1mcos = sign(1.0d0, 1-cosXo)
		
		if (cosXo > 100.0d0) cosXo = 100.0d0
		if (cosXo > 0.99 .AND. cosXo < 1.00) cosXo = 0.99d0
		if (cosXo > 1.00 .AND. cosXo < 1.01) cosXo = 1.01d0
c		call cplxacos(cosXo, Xo)


! Locating the correct data file(s) and verifying if interp. is needed
	
	cosXotest=0.20d0
	call find_datafile(cosXo, data_file1, data_file2)
c	print *, 'cosXo=', cosXotest
c	print *, data_file1
c	print *, data_file2
	
	if (data_file1 .eq. data_file2) then
	   interpcos = .FALSE.
	else
	   interpcos = .TRUE.
	end if

! Auxiliary vectors for reading the data files: FORMAT identifiers

	call INT_TO_STRING(2*Nm,  charaux)
	call INT_TO_STRING(2*Nm+1,charaux2)
	     GGs  = "(G11.5," // charaux // "G20.10)"
	     GG2s = "(" // charaux2 // "G20.10)"

! Reading the first data file (desired A-values will be sorted later)


	j=0
	open(UNIT = 14, FILE = data_file1, STATUS = "OLD", ACTION = "READ")
        
	  read ( 14, GGs), cosXo1, mdiffaux(1:2*Nm)
	  do 
	    j = j + 1
	    read ( 14, GG2s, END = 222), Atab1(j), m12tab1(j,1:Nm)
	  end do 
222	  continue
	  filesize = j - 1

	close(14)

! Sorting the desired A-values (input Avector)

	do i = 1, NA
	
	   A = dabs(Avector(i))
	   sigma_A = sign(1.0d0, Avector(i))

	   if (A > 1.0d+3) A = 5.0d+3
	   if (A < 6.0d-4) A = 6.0d-4


	   do j = 2, filesize

	     if (Atab1(j) .eq. A) then   ! Exact match : no A interpolation
	        data1(i,1:Nm) = m12tab1(j,1:Nm);
	     end if

	     if (Atab1(j) < A .and. Atab1(j-1) > A) then ! Linear interpolation 
	        x1 = Atab1(j-1)
	        x2 = Atab1(j)
		    y1 = m12tab1(j-1,1:Nm)
		    y2 = m12tab1(j,1:Nm)
	        call lininterp_vec( x1, x2, y1, y2, Nm, A, data1(i,1:Nm))
	     end if

	   end do  ! j

	end do  ! i



! Reading the second data file (if interp. in cosXo is needed)

	if (interpcos .eq. .TRUE.) then
	  
	  j=0
	  open(UNIT = 14, FILE = data_file2, STATUS = "OLD", ACTION = "READ")
          
		read ( 14, GGs), cosXo2, mdiffaux(1:2*Nm)
		do 
	      j = j + 1
	      read ( 14, GG2s, END = 333), Atab1(j), m12tab1(j,1:Nm)
	    end do 
333	    continue
	  
	  close(14)
	
	end if

! Sorting the desired A-values (input Avector) for data_file2

	do i = 1, NA
	
	   A = dabs(Avector(i))
	   sigma_A = sign(1.0d0, Avector(i))
	   if (A > 1.0d+3) A = 5.0d+3
	   if (A < 6.0d-4) A = 6.0d-4

	   do j = 2, filesize

	     if (Atab1(j) .eq. A) then   ! Exact match : no A interpolation
	        data2(i,1:Nm) = m12tab1(j,1:Nm);
c	        sortie(i,1:Nm) = sigma_A * data1(i,1:Nm) 
c			print *, sortie(i,1:Nm)
	     end if

	     if (Atab1(j) < A .and. Atab1(j-1) > A) then ! Linear interpolation 
	        x1 = Atab1(j-1)
	        x2 = Atab1(j)
		    y1 = m12tab1(j-1,1:Nm)
		    y2 = m12tab1(j,1:Nm)
	        call lininterp_vec( x1, x2, y1, y2, Nm, A, data2(i,1:Nm))
c			sortie(i,1:Nm) = sigma_A * data1(i,1:Nm) 
c			print *, 'Interpolation in A value'
c	        print *, y1
c	        print *, y2
c	        print *, (y1+y2)/2.0d0
c		    print *, sortie(i,1:Nm)
	     end if

	   end do  ! j

	end do  ! i


	if (interpcos .eq. .FALSE.) then
	  
	  sortie = sigma_A * data1
	
	else

	  alpha = (cosXo-cosXo2) / (cosXo1-cosXo2)
	  beta =  (cosXo-cosXo1) / (cosXo2-cosXo1)
c	  sortie = sigma_A * (alpha*data1+beta*data2)
	sortie = sigma_A * data1
	end if 

	if (cosXoin < 0) then
	   sortie = -myi * dconjg(-sortie / myi)
	end if


!	Writing data to output files
 
	if (filename .eq. 'nofile') then
	print *, ' --> No file output!'
	go to 2001
	end if

666	continue

! Auxiliary vector for file output: mdiffaux = [0 0 1 1 2 2 ....]

	do k = 1, Nm
		mdiffaux(2*k-1) = mdiff(k)
		mdiffaux(2*k)   = mdiff(k)
	end do 


c  Saves data to file 'm12.dat' -> format : [mbar  Re(M12)  Im(M12)]

	open (UNIT = 13, FILE = filename, STATUS = "REPLACE", 
     ;      ACTION = "WRITE")

	    write ( 13, GGs), cosXoin, mdiffaux(1:2*Nm)
		do j = 1, NA
	        write ( 13, GG2s), Avector(j), sortie(j,1:Nm)
	    end do 
   
	close (13)

c	print *, 'Data written to ', filename

2001	continue

c	tempo = etime (T1)-tempo0
c	print *, 'Time used', tempo, 'seconds'
c	write(7,*) 'Time used', tempo, 'seconds'
c	print *

      END SUBROUTINE M12cokpco_read

! *********************************************************************** !
