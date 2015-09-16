! ************************************************************* !
!															  !
!	SUBROUTINE M12COKPCO(FLAG, i_spec)						  !
!                                                               !
!	Purpose: Calculates the plasma global response in const.  !
!			 k// coordinates in a (A or k//) vs.(m1-m2) mesh. !
!															  !
!	Input:	 FLAG = .T. - response for only one species		  !
!			 FLAG = .F. - response summed over all species	  !
!			 i_spec -> species index for FLAG = .T.			  !
!															  !
!	Output:  m12left(p=+1), m12right(p=-1), m12landau(p=0)	  !
!															  !
!	Calls:   int_to_string2, M12cokpco_write, M12cokpco_read, !
!			 M12cokpco_landau								  !
!															  !
! ************************************************************* !
                                               
      subroutine M12COKPCO(FLAG, i_spec)

c	USE DFPORT		! necessary for 'SYSTEM' function

	implicit none

	include 'pardim.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'compla.copy'
      include 'commod.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'comphy.copy'
	include 'cokpco.copy'

c	Input
	logical, intent(in) :: FLAG
	integer, intent(in) :: i_spec

c	Variables

	integer NA, NM, ispec, kp, NZ,i,j
	integer :: p(3), mdiffvector(2*maxcou+1), ilow, iup
	real*8 :: wcdelta(nspec), wcbar(nspec), vtherm(nspec), 
     ;          wplasma(nspec), cosXo(nspec), Acoef(2*maxpom-1,nspec) 
	real*8 :: zvector(2*maxpom-1,nspec), kvector(2*maxpom-1), rm
	character(100) :: FILE_NAME, DATA_FILE
	character(4) :: rrr
	character(4) :: panam4
	character(10) :: polname(2), polname_u(2)
	integer :: OpenStat
	integer :: mdiffaux(4*maxcou+2)
	character(100) :: GGs, GG2s   
	character(2) :: charaux, charaux2
c	real :: T1(2)	            !   Variables for evaluating
c	real :: tempo, tempo0       !   the elapsed CPU time
	real :: factor

c ------------------------- MAIN PROGRAM ------------------------------

c	tempo0 = etime (T1)
	m12left = czero
	m12right = czero
	m12landau = czero
      m12matrix = czero

c     0) Allow evaluation of M12 with only one species (FLAG = T)  

	if(FLAG) then
	   ilow = i_spec
	   iup  = i_spec
	else
	   ilow = 1
	   iup  = nspec
	end if

c     1) Poloidal mode differences (mdiff = m_i-m_j)
c	   (only use positive mdiff values, symmetrie applied later)

	NM = klim + 1  ! NM = 2*klim + 1
	do j = 1, NM
	   mdiffvector(j) = (j - 1)  ! mdiffvector(j) = (j - 1) - klim
	end do
c	print *, NM
c	print *, mdiffvector
c	stop

	if (WRITE_SCREEN) then
	  print *
	  print *, '---------------------------------------------------'
	  print *, 'Radial point =', intab, ' -> rho(m) = ', abscis(intab)
	  print *, '---------------------------------------------------'
	end if
	if (WRITE_REPORT) then
	  write(7777,*)
	  write(7777,*) '---------------------------------------------------'
	  write(7777,*) 'Radial point =', intab, ' -> rho(m) = ', abscis(intab)
	  write(7777,*) '---------------------------------------------------'
	end if


c     2) Create folder for the radial point: 
c		'../../M12run/rXXXX/', where XXXX is rho(mm)     

	if (WRITE_OUTPUT) then
	    open (UNIT = 15, FILE = COKFOLDER // '/index_cos.dat', 
     ;          STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	       if (OpenStat > 0) then
		      print *, 'Error writing file: ', COKFOLDER // '/index_cos.dat'
	          print *, 'Folder ', COKFOLDER, ' may not exist!'
	          stop
	       end if
	end if 

c     3) Wave polarizations (left, right)
	
      p = (/+1, -1, 0/) 
	polname(1) = "_lefty.dat"
      polname(2) = "_right.dat"
	polname_u(1) = "_lefty.unf"
      polname_u(2) = "_right.unf"	

	do kp = 1, 2 ! Polarization loop ++++++++++++++++++++++++++++++++++++++
	  
	   if (WRITE_SCREEN) then
	     print *
	     write(*,'(G24.1, G3.2)') '---> Polarization : p =', p(kp)
	   end if
	   if (WRITE_REPORT) then
	     write(7777,*)
	     write(7777,'(G24.1, G3.2)') '---> Polarization : p =', p(kp)
	   end if

c     3) Calculate A = p.wcdelta/(k//.vT) and cosXo=(1-p*wgen/wcbar)*Bbar/deltaB 
c        for all plasma species (ispec)

c     3.1) Species parameters:

	  do ispec = ilow, iup
	     vtherm(ispec) = vttab(intab,ispec)				! Thermal velocity
	     wplasma(ispec)= omptab(intab,ispec)			! Plasma frequency
	     wcdelta(ispec)= qom(ispec) * delb(intab)		! wc_delta
	     wcbar(ispec)  = qom(ispec) * bbar(intab)		! wc_bar
	     cosXo(ispec)  = (1 - p(kp)*omegag/wcbar(ispec))! Ref.angle Xo
     ;                        * bbar(intab) / delb(intab) 
	  end do

c     3.2) List of k// and A based on average pol.mode (rm):

        i = 0
        do rm = dfloat(minf(iel)), dfloat(msup(iel)), 0.5d0
		 i = i + 1
		 allkpa(i) = hachi(intab) * (rm + n * qfactor(intab))  ! k// values
		 do ispec = ilow, iup
			if ( allkpa(i)*vtherm(ispec) .eq. 0.0d0) then
				Acoef(i,ispec) = sign(1.0d5, p(kp)*wcdelta(ispec) )
			else
		 		Acoef(i,ispec) = p(kp) * wcdelta(ispec) / 
     ;                           ( dabs(allkpa(i))*vtherm(ispec) )
			end if
	     end do ! ispec
	  end do ! rm(i)
	  NA = 2*(modva2-modva1)+1

c     4) Calculating M12 data  ---------------------------

c	The M12 matrix is calculated in subroutine M12cokpco_write 
c	for p=+1 and p=-1 (p=0 is calculated later)

	  do ispec = ilow, iup ! Loop over species -----------------------

		panam4 = paname(ispec)
	    if (WRITE_SCREEN) then
	       print *,'Specie: ', paname(ispec)
	       print *,'CosXo = ', cosXo(ispec) 
	       print *,'Acoef(min-max):',Acoef(1,ispec),Acoef(NA,ispec)
	    end if
	    if (WRITE_REPORT) then
	       write(7777,*)'Specie: ', paname(ispec)
	       write(7777,*)'CosXo = ', cosXo(ispec) 
	       write(7777,*)'Acoef(min-max):',Acoef(1,ispec),Acoef(NA,ispec)
	    end if

          if (READ_TABLES .eq. .FALSE.) then 
		   ! Calculate M12 matrix explicitly
		   call M12cokpco_write ( cosXo(ispec), 
     ;                              mdiffvector(1:NM), NM, 
     ;                              Acoef(1:NA,ispec), NA,
     ;                              m12matrix(ispec,1:NA,1:NM) )

	    else  

		   ! Reading data stored in table files
 		   call int_to_string4(nint(1000*abscis(intab)), rrr)
	       DATA_FILE = "../../M12tables/r" // rrr //
     ;                   "/M12" // panam4 // polname_u(kp) 

		   call M12cokpco_read ( DATA_FILE, 
     ;                             mdiffvector(1:NM), NM, 
     ;                             Acoef(1:NA,ispec), NA,
     ;                             m12matrix(ispec,1:NA,1:NM) )
	
		end if 

c     5) Writing data to output file

		if (WRITE_OUTPUT) then

              FILE_NAME = COKFOLDER  // "/M12" // panam4 // polname(kp)
		    write ( 15, "(F10.5,'  ',A80)"), cosXo(ispec), FILE_NAME

c             5.1) Auxiliary vector for file output: 
c                  mdiffaux = [0 0 +1 +1 .... +klim +klim]
	        do k = 1, NM
		       mdiffaux(2*k-1) = mdiffvector(k)
		       mdiffaux(2*k)   = mdiffvector(k)
	        end do 

c             5.2) Auxiliary vector for file output: FORMAT identifier

	        call INT_TO_STRING(2*NM,  charaux)
	        call INT_TO_STRING(2*NM+1,charaux2)
	             GGs  = "(G11.5," // charaux // "G20.8)"
	             GG2s = "(" // charaux2 // "G20.8)"

c             5.3) Save FORMATTED M12 data to file:

			open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
				  write (13, GGs), cosXo(ispec), mdiffaux(1:2*NM)
				  do k = 1, NA
					 write (13, GG2s), Acoef(k,ispec), m12matrix(ispec,k,1:NM)
				  end do 
			close (13)

		end if ! (WRITE_OUTPUT)

c	    Species (and radius) dependent factor:	     
	    factor = eps0 * dPsidr_n(intab)*abscis(intab) / hachi(intab) *
     ;             wplasma(ispec)**2 /  wcdelta(ispec) 

		m12matrix(ispec,1:NA,1:NM) = m12matrix(ispec,1:NA,1:NM) * factor 
     ;		                         * (-ci) * 0.1216186d-6 * fregag/b0

	  end do ! ispec : loop over species ----------------------------

	  if (p(kp) .eq. 1) then
		  m12left(1:NA,1:NM) = sum(m12matrix(ilow:iup,1:NA,1:NM), DIM=1)
	  end if
        if (p(kp) .eq. -1) then
		  m12right(1:NA,1:NM) = -sum(m12matrix(ilow:iup,1:NA,1:NM), DIM=1)
	  end if

	end do ! kp : end of Polarization loop +++++++++++++++++++++++++++++++++


c     6) Parallel term (p = 0) ++++++++++++++++++++++++++++++++++++++++++++++++++++
      m12matrix = czero
	
	if (WRITE_SCREEN) then
	  print *
	  print *, '---> Landau term : p = 0'
	end if
	if (WRITE_REPORT) then
	  write(7777,*)
	  write(7777,*) '---> Landau term : p = 0'
	end if

c       List of k// and z = w / (k//.vT) based on average pol.mode (rm):

        i = 0
        do rm = dfloat(minf(iel)), dfloat(msup(iel)), 0.5d0
           i = i + 1
	     do ispec = ilow, iup
		   if ( allkpa(i)*vtherm(ispec) .eq. 0.0d0) then
			  zvector(i,ispec) = sign( 1.0d5, p(kp)*wcdelta(ispec) )
		   else
			  zvector(i,ispec) = omegag / (dabs(allkpa(i))*vtherm(ispec)) 
		   end if
		   kvector(i) = 1
		 end do ! ispec
	  end do ! rm(i)
       NZ = 2*(modva2-modva1)+1	  
	  	
        do ispec = ilow, iup ! Loop over species -----------------------

	    if (WRITE_SCREEN) then
	      print *, 'Specie: ', paname(ispec)
	      print *, 'z(min-max):',zvector(1,ispec),zvector(NZ,ispec)
	    end if
	    if (WRITE_REPORT) then
	      write(7777,*) 'Specie: ', paname(ispec)
	      write(7777,*) 'z(min-max):',zvector(1,ispec),zvector(NZ,ispec)
	    end if

	    if (WRITE_OUTPUT) then
             panam4 = paname(ispec)
		   FILE_NAME = COKFOLDER  // "/M12" // panam4 // "_landau.dat" 
		else
		   FILE_NAME = 'nofile'
	    end if

          call M12cokpco_landau (mdiffvector(1:NM), NM, 
     ;                           zvector(1:NZ,ispec), kvector(1:NZ), NZ, 
     ;                           FILE_NAME, 
     ;                           m12matrix(ispec,1:NZ,1:NM)) 

cERN	CHECK CHECK CHECK CHECK ------------------
c	    Species (and radius) dependent factor (p=0):	     
	    factor = eps0 * dPsidr_n(intab)*abscis(intab) / hachi(intab) *
     ;             wplasma(ispec)**2 /  omegag 

		m12matrix(ispec,1:NZ,1:NM) = m12matrix(ispec,1:NZ,1:NM) * factor 
     ;								* (-2.0d0) * 0.1216186d-6 * fregag/b0

	  end do ! ispec : loop over species ----------------------------
	  
	  m12landau(1:NZ,1:NM) = sum(m12matrix(ilow:iup,1:NZ,1:NM), DIM=1)
		
c	End of Landau term ++++++++++++++++++++++++++++++++++++++++++++++++++++
	

c     7) Apply (m1-m2) symmetrie to all matrices 
c	  (1) Shift calculated values to correct storage place (all right)
c	  (2) Fill the missing (negative) m1-m2 terms

		m12left(1:NA, NM:NM+klim)   = m12left(1:NA, 1:NM)
	    m12left(1:NA, 1:NM-1)       = m12left(1:NA, NM:2)
		m12right(1:NA, NM:NM+klim)  = m12right(1:NA, 1:NM)
	    m12right(1:NA, 1:NM-1)      = m12right(1:NA, NM:2)
		m12landau(1:NZ, NM:NM+klim) = m12landau(1:NZ, 1:NM)
	    m12landau(1:NZ, 1:NM-1)     = m12landau(1:NZ, NM:2)

	
	close(15)  ! close 'index_cos.dat' file
	
      END SUBROUTINE M12COKPCO

! *********************************************************************** !

