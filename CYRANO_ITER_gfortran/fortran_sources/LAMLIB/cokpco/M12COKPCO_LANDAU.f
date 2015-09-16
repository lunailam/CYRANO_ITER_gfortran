! ******************************************************************* !
!																	!
! SUBROUTINE M12cokpco_landau ( mdiff,	: m1-m2 vector				!
!                               Nm,		: length of mdiff			!
!							  zvector,	: z=w/(k//.vT) vector		!
!							  kvector,	: k// vector				!
!							  Nz,		: length of z and k//		!
!							  filename,	: 'datafile.dat' or 'nofile'!
!							  sortie	: output vector				!
!						    )										!
!																	!
! ******************************************************************* !
                                               
      SUBROUTINE M12cokpco_landau (mdiff, Nm, zvector, kvector, Nz, 
     ;           filename, sortie) 

	USE DFPORT		! necessary for 'etime' function

      IMPLICIT NONE

      include 'pardim.copy'
      include 'comgeo.copy'
	include 'comequ.copy'
	include 'comphy.copy'

! Input	                                                                 
	integer, intent(in) :: Nm,Nz
	integer, intent(in) :: mdiff(Nm)  
	real*8, intent(in)  :: zvector(Nz), kvector(Nz)
	character(*), intent(in) :: filename		   
! Output	                                                                 
	complex*16, intent(out) :: sortie(Nz,Nm)

	real, dimension (2) :: T1	! - Variables for evaluating
	real :: tempo, tempo0               !   the elapsed CPU time
! Parameters

      real*8  :: kpar          ! Orbit parameters
	complex*16 :: Pell(Nm)   ! - P_ell coeficient
      complex*16 :: I2(Nz), ZFC(Nz)  


	real*8 :: aux1(Nz), aux2(Nz)
	integer :: intaux(Nz)
	integer :: j, k		! - loop index
	integer :: mdiffaux(2*Nm), mm
	character(100) :: GGs, GG2s
	character(2) :: charaux, charaux2


! Main program

	tempo0 = etime (T1)

!	Subroutine riftab calculates the Fried-Conte dispersion function ZFC
!	of argument zvector = omega / |k//.vT|
	call riftab(Nz, zvector, kvector, ZFC, aux1, aux2, intaux, onlyab_now)

	I2 = (zvector/omegag) * zvector * (1 + zvector * ZFC)
c	Obs: Extra multiplication by (zvector/omegag) to include 
c	     the  1 / |k//.vT| factor
	
	do k = 1, Nm ! beginning of m1-m2 loop ----------------------

	   mm = mdiff(k)
	   Pell(k) = gcdr(0,mm)
	   sortie(:,k) = - ci * Pell(k) * I2(:)

	end do ! i- end of m1-m2 loop ---------------------------

! Writing data to file

	if (filename .eq. 'nofile') goto 2001

! Auxiliary vector for file output: mdiffaux = [0 0 1 1 2 2 ....]

	do k = 1, Nm
		mdiffaux(2*k-1) = mdiff(k)
		mdiffaux(2*k)   = mdiff(k)
	end do 

! Auxiliary vector for file output: FORMAT identifier

	call INT_TO_STRING(2*Nm,  charaux)
	call INT_TO_STRING(2*Nm+1,charaux2)
	     GGs  = "(G11.5," // charaux // "G20.10)"
	     GG2s = "(" // charaux2 // "G20.10)"

c  Saves data to file 'm12.dat' -> format : [mbar  Re(M12)  Im(M12)]

	open (UNIT = 13, FILE = filename, STATUS = "REPLACE", 
     ;      ACTION = "WRITE")

	    write ( 13, GGs), 0, mdiffaux(1:2*Nm)
		do j = 1, Nz
	        write ( 13, GG2s), zvector(j), sortie(j,1:Nm)
	    end do 
   
	close (13)

c	print *, 'Data written to file:', filename

2001	continue

c	tempo = etime (T1)-tempo0
c	print *, 'Time used', tempo, 'seconds'
c	print *
  
      END SUBROUTINE M12cokpco_landau

! *********************************************************************** !
