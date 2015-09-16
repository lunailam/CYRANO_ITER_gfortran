! **************************************************************** !
!                                                                  !
! SUBROUTINE M12cokpco_write ( cosXo,   : value of cos(Xo)		 !
!                              mdiff,   : vector with m1-m2 values !
!                              Nm,      : length of mdiff vector	 !
!                              Avector, : vector with A values	 !
!                              NA,      : length of Avector		 !
!                              filename : name of output file		 !
!                              sortie ) : M12(A,m1-m2) matrix      !
!                                                                  !
! **************************************************************** !
                                               
      SUBROUTINE M12cokpco_write (cosXoin, mdiff, Nm, Avector, NA, 
     ;           filename, sortie) 

      IMPLICIT NONE

      include 'pardim.copy'
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
      real*8,     parameter :: pi = 3.14159265d0
	complex*16, parameter :: myi = (0.0d0,1.0d0)

! Parameters

      real*8  :: A, cosXo					! Orbit parameters
	integer, parameter :: Nchi = 200	! Number of chi-harmonics 
	integer :: Lmax,mm					! Max. value of Ell for each A value
	integer    :: ell(Nchi+51)			! Chi-harmonic index
	complex*16 :: Tell(Nchi+51)			! T_ell coeficient
	complex*16 :: Pell(Nchi+51,Nm)	    ! P_ell coeficient
	complex*16 :: PRODU(Nm)				! Product (P_ell x T_ell)
	complex*16 :: z_sad(Nchi+51)		! Saddlepoint location
	complex*16 :: yy(Nchi+51)			! y values for Clenshaw's formula

	real*8 :: B, C, D, btilde, mutilde, aux, mu, deltamu, 
     ;          aux2, sigmacos, sigma1mcos, sigmaA, alpha, pi2  
	complex*16 :: Fclen, CNN, Xo
	integer :: N, Laux(NA)	
	integer :: root_branch
	integer :: j, k, i, kk, k2	! - loop index
	integer :: mdiffaux(2*Nm)
	character(100) :: GGs, GG2s   
	character(2) :: charaux, charaux2
	integer :: Ell_min(Nm), Ell_max(Nm)
	character(8) :: METHOD

! Main program

!     1) Method to be used:

c	METHOD='clenshaw'  ! Use CLENSHAW's backward recursion
      METHOD='directsu'  ! Calculate directly the sum (P_ell x T_ell)
					   ! using ONLY the recursion for T_ell


!     2) Definitions and limits for cosXo: 
!        (Use always |cosXo| and apply symmetrie later)

	cosXo = abs(cosXoin)
      sigmacos   = sign(1.0d0, cosXoin)
      sigma1mcos = sign(1.0d0, 1-cosXo)
	if (cosXo > 100.0d0) cosXo = 100.0d0
	call cplxacos(cosXo, Xo)


!	3) NEW: Determine Ell_min and Ell_max for truncating the sum
!	  based on the values of the geom.coef. P_ell(ell,m1-m2). 
!       NB: This is only used if METHOD is NOT Clenshaw.
!       Default values are: Ell_min=0 and Ell_max=npfft

	if (METHOD .ne. 'clenshaw') then
	   Ell_min=0		
	   Ell_max=npfft
	   do i=1,Nm
	      mm = mdiff(i)
	      call pell_limit(mm, Ell_min(i), Ell_max(i) )
	   end do
	end if


!     4) Loop over A-values 
	
	do j = 1, NA  ! beginning of A-loop +++++++++++++++++++++++++++++

!	   Definitions and limits for A 
!        (Use |A| and apply symmetrie later)
	   A = dabs(Avector(j))
	   sigmaA = sign(1.0d0, Avector(j))
	   if (A > 1.0d+4) A = 1.0d+4
	   if (A < 5.0d-4) A = 5.0d-4


!	   Approximation for small A (smallA.doc) -----------------
!	   (Use only the first 2 terms in the sum, 
!         since the others vanish)

	   if ( abs(A*cosXo) < 0.01) then
	 	  
		  pi2 = pi**1.5d0
		  sortie(j,1:Nm) = (A*sigmaA*pi2 - myi*A**2*cosXo*pi2) 
     ;                        * gcdr(0,mdiff(1:Nm)) * sigma1mcos 
     ;                        + myi*A**2/2.0d0*pi2
     ;                        * gcdr(1,mdiff(1:Nm)) * sigma1mcos
		  goto 3333 ! Go to next A value

	   end if

!        --------------------------------------------------------


!	   Approximation for very large A (verylargeA.doc) --------
!	   ( NB: NOT valid for cosXo = 1)
	   
	   if (cosXo .ne. 1.0d0) then
		   
		   if (A > 1000.0d0) then
			   
			   A = 10000.0d0
			   Lmax = 10
	           do k = Lmax-2, Lmax +1 ! ell-loop - - - - -  
	              ell(k) = k-1				 
			      Tell(k) = myi*pi * exp(myi * dfloat(ell(k)) * Xo) 
     ;                      / (A * sin(Xo)) * sigma1mcos
			   end do !  - - - - - - - - - - - - - - - - - 

		       if (METHOD .ne. 'clenshaw') then
	              ! Calculate the other elements
		          do k = 1, Lmax-3  ! ell-loop - - - - -  
	              kk = Lmax-2-k
	              ell(kk) = kk-1
	              ! Backward recurrence for Tell
	              Tell(kk) = Tell(kk+4) - 2*(ell(kk)+2)/(A**2)*Tell(kk+2) 
     ;                         -2*CosXo*(Tell(kk+3)-Tell(kk+1)) 
			      end do  !  - - - - - - - - - - - - - - - - - 
		       end if

			   go to 1001
		   
		   end if ! (A > 1000)
		
	   end if ! (cosXo .ne. 1)

	   if (cosXo .eq. 1.0d0 .and. A > 1000.0d0) then
	       ! No approx. for this case. 
             ! Only avoid too high A values.
		   A = 1000.0d0
	   end if

!        -----------------------------------------------------------


!	   Approximation for large A (largeA.doc) --------------------
!        (Use only the first 3 terms in the sum, 
!         since the T_ell spectrum is rather narrow.
!	    NB: ONLY valid for cosXo > 1)

	   if ( cosXo > 1.99 .AND. abs(A*cosXo) > 8) then
	 	  
	     aux2 = sqrt(cosXo**2-1)
	     alpha = A**2 * cosXo**2

c          sortie(j,1:Nm) = +myi*pi/aux2* gcdr(0,mdiff(1:Nm)) 
c     ;                     +myi*pi/aux2* (cosXo-aux2)* gcdr(1,mdiff(1:Nm))
c     ;                     +myi*pi/aux2* (cosXo-aux2)*gcdr(2,mdiff(1:Nm))
c     ;                   2*cosXo*(cosXo-aux2-1/2/cosXo) * gcdr(2,mdiff(1:Nm))

	     sortie(j,1:Nm) = + myi*pi/cosXo * gcdr(0,mdiff(1:Nm)) 
     ;                      * (cosXo/aux2 + 
     ;        1/(2*alpha)*1/2*cosXo**3*(2*cosXo**2+1)/aux2/(cosXo**2-1)**2) 
     ;                      + myi*pi/cosXo * gcdr(1,mdiff(1:Nm)) 
     ;                      * (cosXo/aux2 *(cosXo-aux2) +
     ;        1/(2*alpha) * 3/2*cosXo**4 / aux2 / (cosXo**2-1)**2) 
     ;                      + myi*pi/cosXo * gcdr(2,mdiff(1:Nm) )
     ;                      * (0.5/aux2*(cosXo-aux2-1/2/cosXo))
		
		 goto 3333 ! Go to next A value
	
	   end if

!	----------------------------------------------------


! -------------------- NOT USED ------------------------------------

c       Aproximation for large |A.cosXo| and cosXo>2

c	  if (abs(cosXo) > 2.0d0 .and. abs(A*(1.0d0-cosXo))>50) then 

c	    Lmax = 5
c	    do k = Lmax-2, Lmax +1 ! beginning of ell-loop 
c	       ell(k) = k-1
ccc		   Tell(k) = -pi/A * exp( -abs(Xo)*dfloat(ell(k)) ) / 
ccc     ;                  sinh(abs(Xo)) * sigmacos**(ell(k)+1)
c			   Tell(k) = myi*pi * exp(myi * dfloat(ell(k)) * Xo) 
c     ;                   / (A*sigmaA * sin(Xo)) * sigma1mcos
c	    end do ! k: end of ell-loop 
ccc		print *, '--> Approx large |A.cosXo| = ', abs(A*(1.0d0-cosXo))
	    
c		if (METHOD .ne. 'clenshaw') then
c		  do k = 1, Lmax-3  ! beginning of (second) ell-loop
c	      kk = Lmax-2-k
c	      ell(kk) = kk-1
	      !print *, kk, ":", ell2(kk)
	      ! Backward recurrence for Tell
c	      Tell(kk) = Tell(kk+4) - 2*(ell(kk)+2)/(A**2)*Tell(kk+2) 
c     &                  -2*CosXo*(Tell(kk+3)-Tell(kk+1)) 
	
c		  end do  ! end of (second) ell-loop
c		end if

c		go to 1001
c	  end if ! A.cosXo>50 

c	---------------------------------------------------------------


!     5) Beginning of numerical calculations

!     5.1) Value of maximum Chi-harmonic (Lmax) for starting 
!          the Clenshaw's recursion  

	   Lmax = int (19.17d0*dabs(A) + 60.0d0)  ! empirical expression
	
	   if (dfloat(Lmax)/2.0 .ne. Lmax/2) then ! only even Lmax values
	       Lmax = Lmax+1
	   end if
	   if (Lmax > Nchi) then		! upper limit
	       Lmax = Nchi
	   end if
	   if (Lmax < 10) then			! lower limit
	       Lmax = 10
	   end if

	   mu = dfloat(Lmax+1) / (2.0d0*A*A)  ! mu parameter


!	5.2) Choose the correct root branch for saddlepoints method
!	    (see saddle_roots.doc)

	   if (cosXo > 1) then
		
		  aux = abs(cosXo)/4.0 + sqrt(0.5d0+(cosXo/4.0)**2)	
		  btilde = log( aux + sqrt(aux**2-1.0d0) )
		  mutilde = sinh(btilde)**3 / cosh(btilde)
		  deltamu=0.5*mutilde
		
		  if (mu<mutilde) then
			 root_branch = 1
			 ! Avoid points near mu=mutilde
			 if (mu > mutilde - 0.5*deltamu) then
			    Lmax = Lmax - 50
			    if (Lmax < 10) Lmax = 10
			 end if
		  else !(mu>=mutilde)
			 root_branch = 2
	         ! Avoid points near mu=mutilde
			 if (mu <= mutilde +5*deltamu) then
			    Lmax = Lmax + 50
			    if (Lmax > Nchi) Lmax = Nchi
			 end if
		  end if
	
	   end if ! cosXo > 1
	
	   if (cosXo <= 1) then
		   root_branch=3   ! (similar to root_branch 2)
	   end if


c	if (METHOD .eq. 'clenshaw') then

!       5.3) Loop to calculate the 4 last T_ell coefficients 

		do k = Lmax-2, Lmax +1 ! beginning of ell-loop .-.-.-.-.-.
 		   ell(k) = k-1
		   if (abs(A) <= 1000.0d0) then
		     call TT_ell (A,cosXo,ell(k),root_branch,z_sad(k),Tell(k))

		   end if
		end do ! k: end of ell-loop .-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	
c	else ! Using only Recursion in Tell

cc		Lmax=30
cc		mu = dfloat(Lmax+1) / (2.0d0*A*A)
		

cc		if (mu<mutilde) then
cc			root_branch=1
cc			if (mu > mutilde - 5*deltamu) then
cc			Lmax=Lmax-20
cc			end if
cc		else !(mu>=mutilde)
cc			root_branch=2
cc			if (mu <= mutilde + 5*deltamu) then
cc			Lmax=Lmax+50
cc			end if
cc		end if

          
c	   do k = Lmax - 2, Lmax + 1 ! beginning of (first) ell-loop
c	      ell(k) = k -1
c           call TT_ell (A, CosXo, ell(k), root_branch, z_sad(k), Tell(k))
c           print *, ell(k)," -> ", Tell(k)!,"/",imag(z_sad(k))
   
c	   end do ! end of (first) ell-loop  
	
	if (METHOD .ne. 'clenshaw') then

	   do k = 1, Lmax-3  ! beginning of (second) ell-loop
	      kk = Lmax-2-k
	      ell(kk) = kk-1
	      !print *, kk, ":", ell2(kk)
	      ! Backward recurrence for Tell
c           call TT_ell (A, CosXo, ell(kk), root_branch, z_sad(kk), Tell(kk))


	      Tell(kk) = Tell(kk+4) - 2*(ell(kk)+2)/(A**2)*Tell(kk+2) 
     &                  -2*CosXo*(Tell(kk+3)-Tell(kk+1)) 


		end do  ! end of (second) ell-loop
	      
	end if  ! METHOD .ne. 'clenshaw'



1001	continue

!	Symmetries ------------------------ WARNING -----------------

	if (cosXoin < 0) then
	Tell=dconjg(Tell)
	end if

	if (sigmaA < 0) then
	Tell= -dconjg(Tell)
	end if


	do i = 1, Nm ! beginning of m1-m2 loop ----------------------

	mm = mdiff(i)

c     Pell coeficients calculated in FOUGDR2 and stored in gcdr

	Pell(1:10,i) = gcdr(0:9,mm)
	Pell(11:Nchi+51,i) = 0.0d0 + myi*0.0d0
 

	if (METHOD .eq. 'clenshaw') then
		! Clenshaw's backward recursion 

		B = 2 * CosXo
		C = 1/(A*A)
		
		yy(1) = Pell(1,i)
		yy(2) = Pell(2,i) + B * yy(1)
		yy(3) = Pell(3,i) + B * yy(2) - 4*C * yy(1)
		yy(4) = Pell(4,i) + B * yy(3) - 6*C * yy(2) - B * yy(1)

		do k = 5, Lmax +1 ! beginning of ell-loop .-.-.-.-.-.-.-.-
	    ell(k) = k-1
		yy(k) = Pell(k,i) + B * yy(k-1) - (2*ell(k))*C * yy(k-2) 
     ;                      - B * yy(k-3) + yy(k-4)
      	end do	! k: end of ell-loop .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	
		N = Lmax+1
		CNN = Pell(N,i)
	Fclen = yy(N-3) * (Tell(N-3) - B*Tell(N-2)+2*(N-1)*C*Tell(N-1))
     ;      + yy(N-2) * (Tell(N-2) - B*Tell(N-1))
     ;	  + yy(N-1) *  Tell(N-1) + (CNN + yy(N-4)) * Tell(N)   
  	
		sortie(j,i) = -myi * A*sigmaA * Fclen
		
	else
		   
		   PRODU(i) = sum(Tell(:)*Pell(:,i))		   
		   sortie(j,i) = -myi * A*sigmaA * PRODU(i) 		
	
		end if ! METHOD = clenshaw
		
		
	end do ! i- end of m1-m2 loop ---------------------------

3333	continue
	Laux(j)=Lmax

      end do ! j- end of A-loop ++++++++++++++++++++++++++++++++++++++++++


!	Writing data to file
 
	if (filename .eq. 'nofile') then
	print *, ' --> No file output!'
	go to 2001
	end if


! Auxiliary vector for file output: mdiffaux = [0 0 1 1 2 2 ....]

	do k = 1, Nm
		mdiffaux(2*k-1) = mdiff(k)
		mdiffaux(2*k)   = mdiff(k)
	end do 

! Auxiliary vector for file output: FORMAT identifier

	call INT_TO_STRING(2*Nm,  charaux)
	call INT_TO_STRING(2*Nm+1,charaux2)
	     GGs  = "(G11.5," // charaux // "G20.8)"
	     GG2s = "(" // charaux2 // "G20.8)"

c  Saves data to file 'filename' -> format : [A  Re(M12)  Im(M12)]

	open (UNIT = 13, FILE = filename, STATUS = "REPLACE", 
     ;      ACTION = "WRITE")

	    write ( 13, GGs), cosXoin, mdiffaux(1:2*Nm)
		do j = 1, NA
	        write ( 13, GG2s), Avector(j), sortie(j,1:Nm)
	    end do 
   
	close (13)


	open (UNIT = 11, FILE = 'Lmax.dat', STATUS = "REPLACE", 
     ;      ACTION = "WRITE")

		do j = 1, NA
	        write ( 11, *), Avector(j), Laux(j)
	    end do 
   
	close (11)


2001	continue

      END SUBROUTINE M12cokpco_write

! *********************************************************************** !
