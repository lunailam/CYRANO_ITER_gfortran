! *************************************************************** !
!                                                                 !
! subroutine WRITE_M12_TOFILE	(mdiffvec, NM, Avec, NA)  	        !
!                                                                 !
!                                                                 !    
! *************************************************************** !

      SUBROUTINE WRITE_M12_TOFILE (mdiffvec, NM, Avec, NA)

	USE DFPORT		! necessary for 'etime' function

      IMPLICIT NONE

      include 'pardim.copy'
      include 'commod.copy'
      include 'comswe.copy'      
      include 'comgeo.copy' 
      include 'compla.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'cokpco.copy'

c	Input
	real*8 :: Avec(*)
	integer :: mdiffvec(*), NA, NM                                                                 

c	Variables
	integer :: j, NXo, ispec
	integer :: mdiffaux(2*Nm)
	complex*16 :: m12dummy(NA,NM)
	real*8 :: cosXo(100)
	character(150) :: FILE_NAME, FILE_NAME2
	character(4) :: panam4
	character(20) :: specandpol(2*maxspe)
	character(100) :: GGs, GG2s   
	character(2) :: charaux, charaux2
      real*8 :: wcbar(nspec)	
	integer :: OpenStat
	real :: T1(2)	            !   Variables for evaluating
	real :: tempo, tempo0       !   the elapsed CPU time

c ----------------------------- MAIN PROGRAM ----------------------------

	tempo0 = etime (T1)

c     1) Reference angle vector: cosXo
	
	if (ALLCOS_TABLES .eq. .TRUE.) then 
	   ! computing all values of cosXo [0..100]

	   do j = 1, 20
	      cosXo(j) = 0.05d0 * real(j-1)			
	   end do
	      cosXo(21) = 0.97d0
	      cosXo(22) = 0.99d0
	      cosXo(23) = 1.00d0
	      cosXo(24) = 1.01d0
	      cosXo(25) = 1.03d0
	      cosXo(26) = 1.05d0
	      cosXo(27) = 1.07d0
	      cosXo(28) = 1.09d0
	      cosXo(29) = 1.11d0
	      cosXo(30) = 1.15d0
	      cosXo(31) = 1.20d0
	      cosXo(32) = 1.25d0
	      cosXo(33) = 1.30d0
	      cosXo(34) = 1.35d0
	      cosXo(35) = 1.40d0
	      cosXo(36) = 1.45d0
	   do j = 1, 35
	      cosXo(36+j) = 1.25d0 + 0.25d0 * real(j)			
	   end do
	   do j = 1, 10
	      cosXo(71+j) = 10.0d0 + 1.0d0 * real(j)			
	   end do
	   do j = 1, 16
	      cosXo(81+j) = 20.0d0 + 5.0d0 * real(j)			
	   end do
	   NXo = 97

	else  ! ALLCOS_TABLES = .FALSE.
         ! -> computing only the (2*nspec) required values of cosXo

	   do ispec = 1, nspec

	      wcbar(ispec)        =  qom(ispec) * bbar(intab)  ! wc_bar
	      cosXo(ispec)        = (1 - omegag/wcbar(ispec))  ! cosXo (p = +1)
     ;                          * bbar(intab) / delb(intab)  ! left pol.
	      cosXo(ispec+nspec)  = (1 + omegag/wcbar(ispec))  ! cosXo (p = -1)
     ;                          * bbar(intab) / delb(intab)  ! right pol.
	      panam4 = paname(ispec)
		  specandpol(ispec)       = panam4 // '_lefty.dat'
		  specandpol(ispec+nspec) = panam4 // '_right.dat' 
		  	   
	   end do
	   NXo = 2*nspec

	end if   ! METHOD = allcos


!     3) Creating index file for cosXo values (index_cos.dat)

	open (UNIT = 1555, FILE = TABFOLDER // '/index_cos.dat', 
     ;      STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	      if (OpenStat > 0) then
		    print *, 'Error writing file: ', TABFOLDER // '/index_cos.dat'
	        print *, 'Folder ', TABFOLDER, ' may not exist!'
	      end if

c     4) Calculating M12 matrix for each cosXo value and storing 
c        the data in '../../M12tables/rXXXX/M12$$$.dat'

	do j = 1, NXo  ! cosXo loop  ! +++++++++++++++++++++++++++++++++++++

c        4.1) Create the correct file name 
	   if (ALLCOS_TABLES .eq. .TRUE.) then 
	      ! Creating filename with cosXo values 
	      call set_filename( TABFOLDER, cosXo(j) ,FILE_NAME)
	   else
	      ! Creating filename with specie + polariz. values
	      FILE_NAME = TABFOLDER // '/M12' // specandpol(j) 
	   end if
	   
	   write ( 1555, "(F10.5,'  ',A80)"), cosXo(j), FILE_NAME

c	   4.2) Calculate the M12 matrix	
	   m12dummy = 0
	   call M12cokpco_write (cosXo(j), 
     ;                         mdiffvec(1:NM), NM, 
     ;                         Avec(1:NA), NA,
     ;                         m12dummy(1:NA,1:NM) )

c     5) Writing data to output file
	
c        5.1) Auxiliary vector for file output: 
c           mdiffaux = [-klim -klim .... -1 -1 0 0 +1 +1 .... +klim +klim]

	   do k = 1, NM
		  mdiffaux(2*k-1) = mdiffvec(k)
		  mdiffaux(2*k)   = mdiffvec(k)
	   end do 

c        5.2) Auxiliary vector for file output: FORMAT identifier

	   call INT_TO_STRING(2*NM,   charaux)
	   call INT_TO_STRING(2*NM+1, charaux2)
	        GGs  = "(G11.5," // charaux // "G20.8)"
	        GG2s = "(" // charaux2 // "G20.8)"

c        5.3) Save FORMATTED M12 data to file to be read with MATLAB
c           Format: [cosXo     mm_i       mm_i    ...
c                      A    Re(M12)_i  Im(M12)_i  ...]

	   open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
     ;         IOSTAT = OpenStat, ACTION = "WRITE")
	         if (OpenStat > 0) then
		        print *, 'Error writing file: ', FILE_NAME
	            stop
		     end if
		     write ( 13, GGs), cosXo(j), mdiffaux(1:2*NM)
		     do k = 1, NA
	            write ( 13, GG2s), Avec(k), m12dummy(k,1:NM)
	         end do 
	   close (13)

c	   5.4) Save UNFORMATTED data (FASTER) to be read in M12cokpco_read.f

	   if(READ_TABLES) then
	      FILE_NAME2 = FILE_NAME(1:35) // ".unf"
	      open (UNIT = 13, FILE = FILE_NAME2, STATUS = "REPLACE", 
     ;            FORM = 'UNFORMATTED', IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		            print *, 'Error writing file: ', FILE_NAME2
	                stop
		        end if
 	            write ( 13 ), cosXo(j), mdiffaux(1:2*NM)
		       do k = 1, NA
	              write ( 13 ), Avec(k), m12dummy(k,1:NM)
	           end do 
 	      close (13)
	   end if


	end do ! j: end of cosXo loop ++++++++++++++++++++++++++++++++++++++++


	close (1555) ! close 'index_cos.dat'


cccccccccccc  Only for testing 1 cosXo value ccccccccccccccccc
c
c	call M12cokpco_write (2.0d0, mdiffvec(1:NM), NM, 
c     ;                            Avec(1:NA),NA,
c     ;                            '../../M12tables/m12_cosXo2p00.dat',
c     ;                             m12dummy)
c	stop
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	tempo = etime (T1)-tempo0
	print *, 'Time used :', tempo, 'seconds'

      END SUBROUTINE WRITE_M12_TOFILE
