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
	complex*16 :: m12dummy(NA,NM)
	real*8 :: cosXo(100)
	character(150) :: FILE_NAME
	character(4) :: panam4
	character(20) :: specandpol(2*maxspe)
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
		  specandpol(ispec)       = panam4 // '_left.dat'
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
	   call M12cokpco_write (cosXo(j), 
     ;                         mdiffvec(1:NM), NM, 
     ;                         Avec(1:NA), NA,
c     ;                         FILE_NAME,
     ;                         m12dummy)

	end do ! end of cosXo loop ++++++++++++++++++++++++++++++++++++++++


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
