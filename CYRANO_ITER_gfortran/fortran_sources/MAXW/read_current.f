! *************************************************************** !
!                                                                 !
! SUBROUTINE read_current 										!
!             
!  reads current density profile from file (...\Profiles\Jtoroidal.dat) !
!  Store in COMMON jtor_tab(1:nabsci)
!  Use quadratic interpolation INTERP2 to CYRANO abscissa   !

!  *************************************************************** !
                                               
      SUBROUTINE read_current

      IMPLICIT NONE

      include 'pardim.copy'
	include 'comgeo.copy'
	include 'comfin.copy'
	include 'comfic.copy'
	include 'coequi.copy'
	include 'compla.copy'
	include 'commag.copy'
	include 'complp.copy'


! Parameters

	integer :: j, OpenStat, NP
	real*8 :: rhoJ_flu(1000), Jtor_flu(1000)
cccccccccccccccccccccccccccccc Main program cccccccccccccccccccccccccccccc
	
	print *
	print *, '--> Importing current density!'
	print *
	write(nofile,*)	
	write(nofile,*), '----> Importing current density(READ_CURRENT = T)'
	write(nofile,*)



cERN	 07/03/05 : Only import current density (jtot) from FILES 
c	 (All rest will be computed automatically in TABLES and genshap2)
c	 New switch : READ_JTOR

c	Read Toroidal current density from FILE
	open(UNIT = 99, FILE = '../../Profiles/Jtoroidal.dat', STATUS = "OLD",
     ;     IOSTAT = OpenStat, ACTION = "READ")
		 if(OpenStat>0) then
			 print*, 'Jtoroidal file not found!'
			 goto 666
		 else
		     do j = 1,1000
		        read ( 99, *, END = 666), rhoJ_flu(j), Jtor_flu(j)
		     end do
		 end if
	
	close (99)
666	continue
	NP = j - 1

	if(OpenStat>0) then	
c	File ot found: Use standard analytical expression
		print *, ' Using analytical expression for current density.'
		jtor_tab(1:i_rw) = (1-abscno(1:i_rw)*abscno(1:i_rw))**alpha
	else
c		Interpolate to CYRANO abscissa
		call interp2(rhoJ_flu, 1, Jtor_flu, 1, NP, 
     ;                 abscno(1:i_rw), jtor_tab(1:i_rw), i_rw )
c		do j = 1, nabsci
c		   call lininterp_onepoint(NP, rhoJ_flu, Jtor_flu, 
c    ;                               abscno(j), jtor_tab(j)) 
c		end do


	end if
	print *, '... done!'


	open (UNIT = 777, FILE = paplda // 'jtor.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		write(777,*), 'Toroidal current density'
		write(777,*), ' rad ', ' pla ' 		
		write(777,"(i5, i5)"), nabsci, i_rp
		write(777,*), 'rho(m)', ' jtor (A/m2)'
		do j = 1, nabsci
			write(777,"(g16.6, g16.6)"), abscis(j) , jtor_tab(j)
		end do
	close(777)


1111	continue

      END SUBROUTINE read_current

