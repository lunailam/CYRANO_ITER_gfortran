c     ! *************************************************************** !
c     !                                                                 !
c	! SUBROUTINE M12cokpco_read ( data_file, : File with M12 tables	  !
c	!						      mdiff,     : m1-m2 vector			  !
c	!							  Nm,        : length of mdiff		  !
c	!						 	  Avector,   : A (k//) vector	      !
c	!							  NA,		 : length of A		      !
c	!							  sortie )   : M12 matrix (OUTPUT)    !
c	!                                                                 !
c	! Purpose: Read the M12 files generated in M12cokpco_tables       !
c	!          and compute the actual M12 matrix by interpolation	  !
c	!		   of the A (or k//) values required.					  !
c	!																  !
c	! *************************************************************** !
                                               
      subroutine M12cokpco_read (data_file, mdiff, Nm, Avector, NA, sortie) 

      implicit none

      include 'pardim.copy'
      include 'commod.copy'
      include 'cokpco.copy'
      include 'comswe.copy'      
      include 'comgeo.copy' 

!	Input	                                                                 
	character(*), intent(in) :: data_file
	integer, intent(in) :: Nm,NA
	integer, intent(in) :: mdiff(Nm)  
	real*8, intent(in)  :: Avector(NA)
	   
!	Output	                                                                 
	complex*16, intent(out) :: sortie(NA,Nm)

!	Variables
      real*8 :: A, sigma_A , cosXo
	complex*16 :: data1(NA,Nm)  	
	integer :: j, i	
	integer :: mdiffaux(2*Nm)
	integer :: filesize
	real*8 :: x1, x2
	complex*16 :: y1(Nm), y2(Nm)
	integer :: OpenStat, counter
	real*8 :: alpha, beta

ccccccccccccccccccccccc  Main program cccccccccccccccccccccc
c     Atab_size=214
	counter=0
	data1 = 0

c     1) Reading the (unformatted) data file (in ...\M12tables)
c	   (desired A-values will be sorted later)

	open(UNIT = 144, FILE = data_file, STATUS = "OLD", 
     ;     FORM = 'UNFORMATTED', IOSTAT = OpenStat, ACTION = "READ")
	     if (OpenStat > 0) then
		    print *, 'Error reading file: ', data_file
	        stop
	     end if
	     read ( 144 ), cosXo, mdiffaux(1:2*Nm)
	     do j = 1, 2000
	       read ( 144, END = 222), Atab(j), m12tab(j,1:Nm)
	     end do 
222	     continue
	     filesize = j - 1
	close(144)

c     2) Sorting the desired A-values (input Avector)

	do i = 1, NA
	
	   A = dabs(Avector(i))
	   sigma_A = sign(1.0d0, Avector(i))

	   if (A > 1.0d+3) A = 5.0d+3
	   if (A < 6.0d-4) A = 6.0d-4

	   do j = 2, filesize

	     if (Atab(j) .eq. A) then   ! Exact match : no A interpolation
			counter=counter+1
	        data1(i,1:Nm) = m12tab(j,1:Nm);
	        goto 3333
	     end if

	     if (Atab(j) < A .and. Atab(j-1) > A) then ! Linear interpolation 
	        x1 = Atab(j-1)
	        x2 = Atab(j)
		    y1 = m12tab(j-1,1:Nm)
		    y2 = m12tab(j,1:Nm)
	        call lininterp_vec( x1, x2, y1, y2, Nm, A, data1(i,1:Nm))
c	        alpha = (A-x2) / (x1-x2)
c	        beta =  (A-x1) / (x2-x1)
c	        data1(i,1:Nm) = alpha*y1 + beta*y2

	        goto 3333
	     end if

	   end do  ! j: Atab loop

3333	continue

	end do  ! i: Avector loop

	sortie = sigma_A * data1


      end subroutine M12cokpco_read

! *********************************************************************** !
