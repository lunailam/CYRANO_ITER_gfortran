! ********************************************************************* !
!																	  !
! SUBROUTINE M12cokpco_tables											  !
!																	  !																	!
! ********************************************************************* !
   
                                              
      SUBROUTINE M12cokpco_tables 

	USE DFPORT	  ! necessary for 'etime' function

      IMPLICIT NONE

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'comfin.copy'
      include 'compla.copy'
	include 'commod.copy'
      include 'comswe.copy'
      include 'cokpco.copy'
      include 'complp.copy'
      include 'coequi.copy'

c	Variables
	real*8  :: Avec(2000)
	integer :: OpenStat
      integer :: i, j, mdiffvec(201), Nmdiff, NAvec
	logical :: resp
	real :: T1(2)	            !   Variables for evaluating
	real :: tempo, tempo0       !   the elapsed CPU time
	character(4) :: rrr  
	character(33) :: shell_comm

! ------------------------------- MAIN PROGRAM -------------------------------
  
	tempo0 = etime (T1)

c     1) Poloidal mode couplings (mdiff = m_i - m_j)
c	   (only calculate for positive mdiff values)

	Nmdiff = klim + 1  ! Nmdiff = 2*klim + 1  
	do j = 1, Nmdiff
         mdiffvec(j) = (j-1) ! mdiffvec(j) = (j-1) - klim 
	end do

c     2) Create vector of A = wc_delta / (k//.vT) values (A>0)
c        (Inverse thermal Doppler shift)  

	Avec(1) = 10000.0d0
	Avec(2) = 5000.0d0
	Avec(3) = 2500.0d0
	Avec(4) = 2000.0d0
	Avec(5) = 1500.0d0
	Avec(6) = 1250.0d0
	
	if (FINEAA_TABLES .eq. .TRUE.) then ! Fine grid (~ 1500 points, slower)

 	   do j = 1, 10
	      Avec(6+j) = 1000.0d0 / dfloat(j)
	   end do
	   do j = 1, 1000
	      Avec(16+j) = 100.0d0 / dfloat(2*j)
	   end do
	   do j = 1, 490
	      Avec(1016+j) = 0.05d0 - dfloat(j)*0.0001d0
	   end do
	   NAvec = 1506

	else  ! FINEAA_TABLES .eq. .FALSE.  ! Rough grid (~ 200 points, faster)

	   do j = 1, 40
	      Avec(6+j) = 1000.0d0 / dfloat(j)
	   end do
	   do j = 4, 50
	      Avec(43+j) = 99.0d0 / dfloat(j)
	   end do
         do j = 10, 100
	      Avec(84+j) = 19.0d0 / dfloat(j)
	   end do
	   do j = 1, 10
	      Avec(184+j) = 0.15d0 / dfloat(j)  
	   end do	
	   do j = 1, 20
	      Avec(194+j) = 0.10d0 / dfloat(10*j)  
	   end do
	   NAvec = 214
	
	end if  ! FINEAA_TABLES .eq. .TRUE.
	Atab_size = NAvec


c     3) Create file with index of radial elements
 
	open (UNIT = 2222, FILE = '../../M12tables/index_radii.dat', 
     ;      STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	      if (OpenStat > 0) then
	        print *,'Error writing file: ../../M12tables/index_radii.dat'
	      end if
		  write(2222,'(i8)'), i_rp

c     4) First radial element (magnetic axis) + + + + + + + + + + + + + + + + + + +
	
	intab = 2  ! cERN WARNING: delb(intab=1) = 0 -> divide by zero

	print *, '               ===================== '				
	print *, '               Generating M12 Tables '				
	print *, '               ===================== '				
      print *															
      print *, '------------------------------------------------'	    
      print *, 'Radial index = ', intab								
      print *, 'rho(m) = ', abscis(intab)	

c     4.0) Create folder for data output at radial point: '../../M12tables/rXXXX/'    
	call int_to_string4 (nint(1000*abscis(intab)), rrr)	
	TABFOLDER = '../../M12tables/r' // rrr
	shell_comm = 'mkdir ..\..\M12tables\r' // rrr
	resp = SYSTEM(shell_comm)
	if (resp .eq. -1) then
	    print *, 'Error creating folder: ', TABFOLDER
	    stop
	end if

	write(2222,'(i8, g15.6, A23)'), intab, abscis(intab), TABFOLDER
c	print *, intab, abscis(intab), COKFOLDER

c     4.1) Geom. coefficients (Plm)
	call fougdr_ern(1)
						
c     4.2) Calculate the M12 matrix for several A and cosXo values.
c          - A range : A = 10000 - 0.005 (automatic)
c          - cosXo range : if ALLCOS_TABLES = .T. -> cosXo = 0 - 100.0 (automatic)
c                          if ALLCOS_TABLES = .F. -> use ONLY 2*nspec needed values
c          - Files are written in folder '../../M12tables'

	call WRITE_M12_TOFILE(mdiffvec, Nmdiff, Avec, NAvec)                                  
      print *, '------------------------------------------------'	    


c     5) Other radial elements (inside plasma) + + + + + + + + + + + + + + + + + +
c        NB: - if ALLRAD_TABLES = .T. -> consider all radial points 
c            - if ALLRAD_TABLES = .F. -> consider only element nodes 

	if (.not. ALLRAD_TABLES) then ! Consider only elem. boundaries -> ifiabs(j)

          do j = 2, nele 	! Loop over nodes						
		   
		   i = ifiabs(j)												
		   ireg = iregoa(i)											
		   intab = i												
		   if(vacuum(ireg)) goto 3333

c		   5.0) Create folder for data output at radial point: '../../M12tables/rXXXX/'    
	       call int_to_string4(nint(1000*abscis(intab)), rrr)	
	       TABFOLDER = '../../M12tables/r' // rrr
	       shell_comm = 'mkdir ..\..\M12tables\r' // rrr
	       resp = SYSTEM(shell_comm)
		   write(2222,'(i8, g15.6, A23)'), intab, abscis(intab), TABFOLDER

c            5.1) Geom. coefficients		   									
		   call fougdr_ern(1)											
	       print *														
	       print *, '------------------------------------------------'  
	       print *, 'Radial block = ', intab							
	       print *, 'rho(m) = ', abscis(intab), abscno(intab)			

c            5.2) Calculate the M12 matrix for several A and cosXo values.	       
		   call WRITE_M12_TOFILE(mdiffvec, Nmdiff, Avec, NAvec)                               
	       print *, '------------------------------------------------'	
	
		end do ! loop over nodes (j)	

	else  ! ALLRAD_TABLES = .TRUE. -> Considering all points

          do j = 3, nabsci  ! Loop over all elements						

		   intab = j													
		   if( vacuum(iregoa(j)) ) goto 3333
		   
		   call int_to_string4(nint(1000*abscis(intab)), rrr)	
	       TABFOLDER = '../../M12tables/r' // rrr
	       shell_comm = 'mkdir ..\..\M12tables\r' // rrr
	       resp = SYSTEM(shell_comm)							
		   write(2222,'(i8, g15.6, A23)'), intab, abscis(intab), TABFOLDER		   
		   call fougdr_ern(1)												
	       
		   print *														
	       print *, '------------------------------------------------'  
	       print *, 'Radial index = ', intab						
	       print *, 'rho(m) = ', abscis(intab), abscno(intab)		    
	       call WRITE_M12_TOFILE(mdiffvec, Nmdiff, Avec, NAvec)                              
	       print *, '------------------------------------------------'	
		
		end do ! all points (j)	

	end if

3333	  continue

	  close (2222) ! close 'index_radii.dat'

	  print *
	  print *
	  print *,'Radial block = ',intab, '-> rho = ',abscis(intab)
	  print *,'Plasma-Vacuum Interface (vT = 0 -> A = +/-Inf)'
	  print *
	  print *, '---> M12 tables generated succesfully!!'
	  print *
	  tempo = etime (T1)-tempo0
	  print *, 'Total time used :', tempo, 'seconds'
	  print *
	  print *

      END SUBROUTINE M12cokpco_tables

! *********************************************************************** !
