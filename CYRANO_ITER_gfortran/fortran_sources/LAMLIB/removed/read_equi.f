! *************************************************************** !
!                                                                 !
! SUBROUTINE read_equi (Nflu)										!
!             
!	Fixed number of poloidal points Npol = 257    
                                                !
! *************************************************************** !
                                               
      SUBROUTINE read_equi(Nflu) 

      IMPLICIT NONE

      include 'pardim.copy'
	include 'comgeo.copy'
	include 'comfin.copy'
	include 'comfic.copy'
	include 'coequi.copy'
	include 'compla.copy'

! Input	  
	integer, intent(in) :: Nflu

! Parameters

	integer :: j, k, OpenStat, Npol_flu, Npol_wall
	real*8 :: results(8), rhotest(500), R_wall(1000), Z_wall(1000),
     ;          R0_wall, Z0_wall, R1_wall, KAP_wall, R2_wall,
     ;		  R3_wall, R4_wall, rho_wall

cccccccccccccccccccccccccccccc Main program cccccccccccccccccccccccccccccc
	
	print *
	print *, '--> Importing EXTERNAL(FLUSH) equilibrium!'
	print*
	write(nofile,*)	
	write(nofile,*), '----> Importing EXTERNAL(FLUSH) equilibrium (READ_FLUX = T)'
	write(nofile,*)

c     1) Read FLUSH data --------------------------------------

c     1.1) R coordinate (closed surfaces)
	open(UNIT = 99, FILE = 'Profiles/Rflush.dat', STATUS = "OLD",
     ;     IOSTAT = OpenStat, ACTION = "READ")
		 if (Openstat > 0) then
			print *, '(R,Z) profiles not found. Using D-shape!'
	        geneq = .FALSE.
			dshape = .TRUE.
			goto 1111 
		 end if
		 read ( 99, *), rhoR_flu(1:Nflu)
	     do j = 1,1000
		     read ( 99, *, END = 222), R_flu(1:Nflu,j)
		 end do
	close (99)
222	continue

c     1.2) Z coordinate (closed surfaces)
	open(UNIT = 99, FILE = 'Zflush.dat', STATUS = "OLD",
     ;     IOSTAT = OpenStat, ACTION = "READ")
	     read ( 99, *), rhoZ_flu(1:Nflu)
		 do j = 1,1000
		     read ( 99, *, END = 333), Z_flu(1:Nflu,j)
		 end do
	close (99)
333	continue
	Npol_flu = j-1	! FLUSH poloidal points

c     1.3) R coordinate (Wall)
	open(UNIT = 99, FILE = 'Profiles/Rwall.dat', STATUS = "OLD",
     ;     IOSTAT = OpenStat, ACTION = "READ")
	     do j = 1,1000
		     read ( 99, *, END = 444), R_wall(j)
		 end do
	close (99)
444	continue

c     1.4) Z coordinate (Wall)
	open(UNIT = 99, FILE = 'Profiles/Zwall.dat', STATUS = "OLD",
     ;     IOSTAT = OpenStat, ACTION = "READ")
	     do j = 1,1000
		     read ( 99, *, END = 555), Z_wall(j)
		 end do
	close (99)
555	continue
	Npol_wall = j-1			! WALL poloidal points
	rho_wall = rhowal/ap	! WALL minor radius (normalized)




c     2) Lao-Hirshmann-Wieland (LHW) fit ------------------------------------

c     2.1) Closed surfaces (except axis)
	do k = 2, Nflu
		call LHWFIT(Npol_flu, R_flu(k,1:Npol_flu), Z_flu(k,1:Npol_flu),
     ;                results)
		R0_flu(k) = results(1)	! R magnetic (Rgeo-Shaf(r)-0.5*delta/a*rho^2)
		Z0_flu(k) = results(2)  ! Z magnetic (approx. constant = Zaxis)
		R1_flu(k) = results(3)  ! minor radius (rho)
		KAP_flu(k)= results(4)	! Elongation (kappa)	
		R2_flu(k) = results(5)	! Triangularity (0.5*delta/a*rho^2)
		R3_flu(k) = results(6)	! Quadrangularity 
		R4_flu(k) = results(7)	! ????
	end do

c     2.2) Axis
		R0_flu(1) = Raxis		! R magnetic 
		Z0_flu(1) = Zaxis	    ! Z magnetic 
		R1_flu(1) = 0.0d0		! minor radius 
		KAP_flu(1)= KAP_flu(2)	! Elongation 	
		R2_flu(1) = 0.0d0		! Triangularity 
		R3_flu(1) = 0.0d0		! Quadrangularity 
		R4_flu(1) = 0.0d0		! ???

c     2.3) WALL		
		call LHWFIT(Npol_wall, R_wall, Z_wall, results)
		R0_wall = results(1)	! R magnetic 
		Z0_wall = results(2)	! Z magnetic 
		R1_wall = results(3)	! minor radius 
		KAP_wall= results(4)	! Elongation 	
		R2_wall = results(5)	! Triangularity 
		R3_wall = results(6)	! Quadrangularity 
		R4_wall = results(7)	! ????


c     3) Saving FLUSH coefficients to files ------------------------------

	open(UNIT = 990, FILE = 'Profiles/check/R0flu.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 991, FILE = 'Profiles/check/R1flu.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 992, FILE = 'Profiles/check/R2flu.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 993, FILE = 'Profiles/check/R3flu.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")		 
	open(UNIT = 994, FILE = 'Profiles/check/R4flu.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 995, FILE = 'Profiles/check/Z0flu.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")     		 
	open(UNIT = 996, FILE = 'Profiles/check/KAPflu.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
		 do j = 1,Nflu
		     write ( 990, *), rhoR_flu(j), R0_flu(j)	! R magnetic
		     write ( 991, *), rhoR_flu(j), R1_flu(j)	! minor radius 		 
		     write ( 992, *), rhoR_flu(j), R2_flu(j)	! Triangularity
		     write ( 993, *), rhoR_flu(j), R3_flu(j)	! Quadrangularity 		 
		     write ( 994, *), rhoR_flu(j), R4_flu(j)	! ??? 		 
		     write ( 995, *), rhoR_flu(j), Z0_flu(j)	! Z magnetic		 
		     write ( 996, *), rhoR_flu(j), KAP_flu(j)	! Elongation 		 		 
		 end do
		 write ( 990, *), rho_wall, R0_wall
		 write ( 991, *), rho_wall, R1_wall
		 write ( 992, *), rho_wall, R2_wall		 
		 write ( 993, *), rho_wall, R3_wall		 
		 write ( 994, *), rho_wall, R4_wall
		 write ( 995, *), rho_wall, Z0_wall		 
		 write ( 996, *), rho_wall, KAP_wall		 		 	
	close(990)
	close(991)
	close(992)
	close(993)
	close(994)
	close(995)
	close(996)

c     4) Additional points for PLASMA-WALL interpolation ----------------------
c	   N.B: Last point of FLUSH (separatrix) is neglected. It is 
c		    substituted by (R_sep+R_wall)/2 for smoother fitting.

c     4.1) Last FLUSH point Rsep <- (Rsep + Rwall)/2
	   R0_flu(Nflu)   = (R0_flu(Nflu)   + R0_wall) / 2	
	   Z0_flu(Nflu)   = (Z0_flu(Nflu)   + Z0_wall) / 2	
	   R1_flu(Nflu)   = (R1_flu(Nflu)   + R1_wall) / 2	
	   R2_flu(Nflu)   = (R2_flu(Nflu)   + R2_wall) / 2	
	   R3_flu(Nflu)   = (R3_flu(Nflu)   + R3_wall) / 2
	   R4_flu(Nflu)   = (R4_flu(Nflu)   + R4_wall) / 2
	   KAP_flu(Nflu)  = (KAP_flu(Nflu)  + KAP_wall) / 2	
	   rhoR_flu(Nflu) = (rhoR_flu(Nflu) + rho_wall) / 2

c     4.2) Artificial point (Rsep + 3*Rwall)/4
	   R0_flu(Nflu+1)   = (R0_flu(Nflu)   + R0_wall) / 2	
	   Z0_flu(Nflu+1)   = (Z0_flu(Nflu)   + Z0_wall) / 2  
	   R1_flu(Nflu+1)   = (R1_flu(Nflu)   + R1_wall) / 2  
	   R2_flu(Nflu+1)   = (R2_flu(Nflu)   + R2_wall) / 2	
	   R3_flu(Nflu+1)   = (R3_flu(Nflu)   + R3_wall) / 2	
	   R4_flu(Nflu+1)   = (R4_flu(Nflu)   + R4_wall) / 2	
	   KAP_flu(Nflu+1)  = (KAP_flu(Nflu)  + KAP_wall) / 2
	   rhoR_flu(Nflu+1) = (rhoR_flu(Nflu) + rho_wall) / 2

c     4.3) Wall 
	   R0_flu(Nflu+2)   = R0_wall
	   Z0_flu(Nflu+2)   = Z0_wall
	   R1_flu(Nflu+2)   = R1_wall 
	   R2_flu(Nflu+2)   = R2_wall
	   R3_flu(Nflu+2)   = R3_wall	
	   R4_flu(Nflu+2)   = R4_wall	
	   KAP_flu(Nflu+2)  = KAP_wall
	   rhoR_flu(Nflu+2) = rho_wall

	
c     5) Interpolate FLUSH coefs. with Cyrano abcissa -------------------------

	call spline0(Nflu+2,rhoR_flu,R0_flu, i_rw,abscno(1:i_rw),R0_cyr(1:i_rw))
	call spline0(Nflu+2,rhoR_flu,Z0_flu, i_rw,abscno(1:i_rw),Z0_cyr(1:i_rw))
	call spline0(Nflu+2,rhoR_flu,R1_flu, i_rw,abscno(1:i_rw),R1_cyr(1:i_rw))
	call spline0(Nflu+2,rhoR_flu,R2_flu, i_rw,abscno(1:i_rw),R2_cyr(1:i_rw))
	call spline0(Nflu+2,rhoR_flu,R3_flu, i_rw,abscno(1:i_rw),R3_cyr(1:i_rw))
	call spline0(Nflu+2,rhoR_flu,R4_flu, i_rw,abscno(1:i_rw),R4_cyr(1:i_rw))
	call spline0(Nflu+2,rhoR_flu,KAP_flu,
     ;             i_rw,abscno(1:i_rw),KAP_cyr(1:i_rw))
c	call spline0(Nflu+2,rhoR_flu,rhoR_flu, i_rw,abscno(1:i_rw),rhotest(1:i_rw))


c     6) Derivatives of LHW coefficients -------------------------------------		

      call derivate(nabsci, R0_cyr, abscis, dR0_cyr)
      call derivate(nabsci, Z0_cyr, abscis, dZ0_cyr)
	call derivate(nabsci, R1_cyr, abscis, dR1_cyr)
	call derivate(nabsci, R2_cyr, abscis, dR2_cyr)
	call derivate(nabsci, R3_cyr, abscis, dR3_cyr)
	call derivate(nabsci, R4_cyr, abscis, dR4_cyr)
	call derivate(nabsci, KAP_cyr, abscis, dKAP_cyr)

c     7) Saving CYRANO (interpolated) coefficients ---------------------------

	open(UNIT = 990, FILE = 'Profiles/check/R0cyr.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 991, FILE = 'Profiles/check/R1cyr.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 992, FILE = 'Profiles/check/R2cyr.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 993, FILE = 'Profiles/check/R3cyr.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 994, FILE = 'Profiles/check/R4cyr.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 995, FILE = 'Profiles/check/Z0cyr.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
	open(UNIT = 996, FILE = 'Profiles/check/KAPcyr.dat', 
     ;     STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
		 do j = 1,nabsci
		     write ( 990, *), abscno(j), R0_cyr(j)
		     write ( 991, *), abscno(j), R1_cyr(j)
		     write ( 992, *), abscno(j), R2_cyr(j)
		     write ( 993, *), abscno(j), R3_cyr(j)
		     write ( 994, *), abscno(j), R4_cyr(j)
		     write ( 995, *), abscno(j), Z0_cyr(j)
		     write ( 996, *), abscno(j), KAP_cyr(j)
		 end do
	close(990)
	close(991)
	close(992)
	close(993)
	close(994)
	close(995)
	close(996)
c	-------------------------------------------------------------------


cERN	 07/03/05 : 
c	 Only import current density (jtot) from FILES in read_current.f
c	 All rest will be computed automatically in TABLES and genshap2



c	Magnetic field and PITCH angle
c	Poloidal points are fixed (npp=257) but radial need interpolation

c     1.1) B toroidal (closed surfaces)
c	open(UNIT = 99, FILE = '../../Profiles/Btorflush.dat', STATUS = "OLD",
c     ;     IOSTAT = OpenStat, ACTION = "READ")
c	     read ( 99, *), rhoB_flu(1:Nflu)
c		 do j = 1,1000
c		     read ( 99, *, END = 666), Btor_flu(1:Nflu,j)
c		 end do
c	close (99)
c666	continue

c     1.2) Pitch angle (closed surfaces)
c	open(UNIT = 99, FILE = '../../Profiles/Pitchflush.dat', STATUS = "OLD",
c     ;     IOSTAT = OpenStat, ACTION = "READ")
c	     read ( 99, *), rhoP_flu(1:Nflu)
c		 do j = 1,1000
c		     read ( 99, *, END = 777), Pitch_flu(1:Nflu,j)
c		 end do
c	close (99)
c777	continue

c     1.3) Safety factor (closed surfaces)
c	open(UNIT = 99, FILE = '../../Profiles/qprofflush.dat', STATUS = "OLD",
c     ;     IOSTAT = OpenStat, ACTION = "READ")
c	     read ( 99, *), rhoQ_flu(1:Nflu)
c	     read ( 99, *), Q_flu(1:Nflu)
c	close (99)


c	 1) Magnetic field Btor

c     1.4) B toroidal (closed surfaces)
c	rhoB_flu(Nflu+1) = (rhoR_flu(Nflu) + rho_wall) / 2
c	rhoB_flu(Nflu+2) = rho_wall

c	Btor_flu(Nflu+1,:) = Btor_flu(Nflu,:)
c	Btor_flu(Nflu+2,:) = Btor_flu(Nflu,:)

c	Pitch_flu(Nflu+1,:) = Pitch_flu(Nflu,:)
c	Pitch_flu(Nflu+2,:) = Pitch_flu(Nflu,:)

c	Radial interpolation to Cyrano abscissa 
c	do j = 1, 257
c		call spline0(Nflu+2, rhoB_flu,       Btor_flu(1:Nflu+2,j), 
c    ;                 i_rw,   abscno(1:i_rw), Btor_cyr(1:i_rw,j))
c		call spline0(Nflu+2, rhoB_flu,       Pitch_flu(1:Nflu+2,j), 
c     ;                 i_rw,   abscno(1:i_rw), Pitch_cyr(1:i_rw,j))
c	end do



1111	continue

      END SUBROUTINE read_equi

