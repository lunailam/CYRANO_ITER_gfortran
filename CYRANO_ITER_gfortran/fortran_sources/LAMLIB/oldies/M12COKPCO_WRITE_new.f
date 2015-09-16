! **************************************************************** !
!                                                                  !
! SUBROUTINE M12cokpco_write ( cosXoin, : value of cos(Xo)		 !
!                              mdiff,   : vector with m1-m2 values !
!                              Nm,      : length of mdiff vector	 !
!                              Avector, : vector with A values	 !
!                              NA,      : length of Avector		 !
!                              sortie ) : M12(A,m1-m2) matrix      !
!                                                                  !
! **************************************************************** !
                                               
      SUBROUTINE M12cokpco_write (cosXoin, mdiff, Nm, Avector, NA, sortie)

      IMPLICIT NONE

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'comphy.copy'
      include 'cokpco.copy'

!	Input	                                                                 
	real*8, intent(in)  :: cosXoin
	integer, intent(in) :: Nm,NA
	integer, intent(in) :: mdiff(Nm)  
	real*8, intent(in)  :: Avector(NA)
  
!	Output	                                                                 
	complex*16, intent(out) :: sortie(NA,Nm)

!     Variables

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
     ;          aux2, sigmacos, sigma1mcos, sigmaA, alpha, 
     ;          pi15, aux3, aux4, aux5, aux6  
	complex*16 :: Fclen, CNN, Xo
	integer :: N, Laux(NA)	
	integer :: root_branch
	integer :: j, k, i, kk, k2	! - loop index
c	integer :: mdiffaux(2*Nm)
c	character(100) :: GGs, GG2s   
c	character(2) :: charaux, charaux2
c	character(120) :: filename2
	integer :: Ell_min(Nm), Ell_max(Nm)
	character(8) :: METHOD
c	integer :: OpenStat
c	real*8 :: dbinom
c	external dbinom  

! --------------------------- MAIN PROGRAM -------------------------------

!     1) Method to be used:

c	 METHOD='clenshaw'  ! Use CLENSHAW's backward recursion (SLOWER)
c						! (not working very well yet -> oscillations in M12)
       METHOD='directsu'  ! Calculate directly the sum (P_ell x T_ell)
					    ! using ONLY the recursion for T_ell (working OK)

!     2) Definitions and limits for cosXo: 
!        (Use always |cosXo| and apply symmetries later)

	 cosXo = abs(cosXoin)
       sigmacos   = sign(1.0d0, cosXoin)
       sigma1mcos = sign(1.0d0, 1-cosXo)  ! WARNING 1-|cosXo| ???
	 if (cosXo > 100.0d0) cosXo = 100.0d0
	 call cplxacos(cosXo, Xo)

!     2.1) Additional definitions (needed for cosXo>1): 
	 if (cosXo > 1) then		
	    aux = abs(cosXo)/4.0 + sqrt(0.5d0+(cosXo/4.0)**2)	
	    btilde = log( aux + sqrt(aux**2-1.0d0) )
	    mutilde = sinh(btilde)**3 / cosh(btilde)
	    deltamu = 0.5*mutilde
	 end if


!     3) Geom. coefficients : calculated in FOUGDR2 and 
!        stored in gcdr (-npfft/2 ... +npfft/2, -klim ... +klim)
!	   (NB: Serie is temporarely truncated at ell=100)
	 Pell(1:101,1:Nm)    = gcdr(0:100,mdiff)
	 Pell(102:Nchi,1:Nm) = dcmplx(0.0d0, 0.0d0)	

!     3b) NEW: Determine Ell_min and Ell_max for truncating the sum
!	    based on the values of the geom.coef. P_ell(ell,m1-m2). 
!         NB: This is only used if METHOD is NOT Clenshaw.
!         Default values are: Ell_min=0 and Ell_max=npfft

c	 if (METHOD .ne. 'clenshaw') then
	    Ell_min = 0		
	    Ell_max = npfft/2
	    do i = 1,Nm
	       call pell_limit(mdiff(i), Ell_min(i), Ell_max(i) )
	       Ell_max(i)=min(Ell_max(i),100)  ! Pell's are truncated at ell=100
		   Ell_max(i)=max(Ell_max(i),2)    ! Minimum value ell=2
		end do
c	 end if

c	print *
c	do k = 1,Nm
c	print *, mdiff(k), Ell_min(k), Ell_max(k)
c	end do
c	stop

!     4) Loop over A-values 
	
	do j = 1, NA  ! beginning of A-loop +++++++++++++++++++++++++++++

!	   Definitions and limits for A 
!        (Use |A| and apply symmetries later)
	   A = dabs(Avector(j))
	   sigmaA = sign(1.0d0, Avector(j))
	   if (A > 1.0d+4) A = 1.0d+4
	   if (A < 5.0d-4) A = 5.0d-4


!        --------------------------------------------------------
!	   Approximation for small A (see smallA.doc) 
!	   (Use only the first 2 terms in the sum, 
!         since the intergrals for ell>2 vanish)

	   if ( abs(A*cosXo) < 0.01) then
	 	  
		  pi15 = pi**1.5d0
		  sortie(j,1:Nm) = (A*sigmaA*pi15 - ci*A*A*cosXo*pi15) 
     ;                        * gcdr(0,mdiff(1:Nm)) !* sigma1mcos 
     ;                        + ci*A*A/2.0d0*pi15
     ;                        * gcdr(1,mdiff(1:Nm)) !* sigma1mcos

		  goto 3333 ! Go to next A value

	   end if
!        --------------------------------------------------------

!        --------------------------------------------------------
!	   Approximation for very large A (see verylargeA.doc)
!	   ( NB: Only used for cosXo < 1)
	   
	   if (cosXo < 1.0d0 .and. A > 1000.0d0) then
			   A = 10000.0d0
			   Lmax = maxval(Ell_max)+1      !cccccccc WARNING
	           
			   do k = Lmax-2, Lmax +1 ! ell-loop - - - - -  
	              ell(k) = k-1				 
			      Tell(k) = ci*pi * exp(ci * dfloat(ell(k)) * Xo) 
     ;                      / (A * sin(Xo)) * sigma1mcos
			   end do !  - - - - - - - - - - - - - - - - - 

c		       if (METHOD .ne. 'clenshaw') then
	              ! Calculate the other elements
		          do k = 1, Lmax-3  ! ell-loop - - - - -  
	              kk = Lmax-2-k
	              ell(kk) = kk-1
	              ! Backward recurrence for Tell
	              Tell(kk) = Tell(kk+4) - 2*(ell(kk)+2)/(A**2)*Tell(kk+2) 
     ;                         -2*CosXo*(Tell(kk+3)-Tell(kk+1)) 
			      end do  !  - - - - - - - - - - - - - - - - - 
c		       end if
			   go to 1001
		   
		   end if ! (cosXo<1 and A > 1000)
		
	   if (cosXo .eq. 1.0d0 .and. A > 1000.0d0) then
	       ! No approx. for this case. 
		   A = 1000.0d0
	   end if
!        -----------------------------------------------------------

!        -----------------------------------------------------------
!	   Approximation for large A (see largeA.doc) 
!	   (NB: ONLY valid for cosXo > 1)

	   if ( cosXo > 1.0 .and. abs(A*(1-cosXo)) > 8) then ! cccccccccccc ERNLAST
	 	  
	     aux2 = sqrt(cosXo**2-1)
	     aux3 = cosXo + aux2
	     aux4 = aux3**2-1


              do k = minval(Ell_min)+1, maxval(Ell_max)+1 ! ell-loop - - - - -  
	           ell(k) = k-1				 
			   aux5 = (cosXo-aux2)**ell(k)
			   aux6 = aux3**(3-ell(k))

cERN      Float Underflow for too large values of ell	
		       Tell(k) = 1.0d0/A * (-pi)/aux2 * aux5
     ;                   + 0.5d0/A**3 * (-8*pi) * aux6 / (aux4)**5
cERN		Time consumption of dbinom is VERY LARGE
c		(replaced by equivalent expressions -> 2x FASTER)
c     ;				   * ( dbinom(2+ell(k),0)*dbinom(4,2)*1.0d0
c     ;                      +dbinom(2+ell(k),1)*dbinom(3,2)*aux4
c     ;                      +dbinom(2+ell(k),2)*dbinom(2,2)*aux4*aux4
c     ;                      )
     ;				    * ( 6.0d0
     ;                      +(2+ell(k))*3.0d0*aux4
     ;                      +(2+ell(k))*(1+ell(k))*0.5d0*aux4*aux4
     ;                      )
cERN	24 May 2004
			Tell(k)=-Tell(k);
	        end do !  - - - - - - - - - - - - - - - - - 
		 
		 goto 1001
	
	   end if
        
!        -----------------------------------------------------------


!      5) Numerical calculations

!      5.1) Value of maximum Chi-harmonic (Lmax) for starting 
!           the Clenshaw's recursion  

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


!	 5.2) Choose the correct root branch for saddlepoints method
!	    (see saddle_roots.doc)

	   if (cosXo > 1) then
			
		  if (mu<mutilde) then
			 root_branch = 1
			 ! Avoid points near mu=mutilde
			 if (mu > mutilde - 5*deltamu) then ! cccccccccccc ERNLAST
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


!      5.3) Loop to calculate the 4 last T_ell coefficients 

	   do k = Lmax-2, Lmax +1 ! ell-loop - - - - - - - - 
 		   ell(k) = k-1
		   if ( A <= 1000.0d0 ) then
		     call TT_ell (A,cosXo,ell(k),root_branch,z_sad(k),Tell(k))
		   end if
	   end do ! end of ell-loop - - - - - - - - - - - -
	
	
         ! Calculate all other T_ell elements for direct sum
c		maybe here Lmax could be redefined as Ell_max, 
c		if the backward recurrence holds
	      do k = 1, Lmax-3  ! ell-loop - - - - - - - - - - 
	         kk = Lmax-2-k
	         ell(kk) = kk-1
	         ! Backward recurrence for Tell
	         Tell(kk) = Tell(kk+4) - 2*(ell(kk)+2)/(A**2)*Tell(kk+2) 
     &                    -2*CosXo*(Tell(kk+3)-Tell(kk+1)) 
	      end do  ! end of ell-loop - - - - - - - - - - - -
	      

1001	   continue

!	 5.4) Aplly symmetries in A and cosXo (Symmetries.doc)

	   if (cosXoin < 0) then
	      Tell=-(-1)**ell * dconjg(Tell)
	   end if

	   if (sigmaA < 0) then
	      Tell= -dconjg(Tell)
	   end if
	
cc	goto 3333

!      6) Loop over m1-m2 values

	   do i = 1, Nm ! m1-m2 loop * * * * * * * * * * * * * * * * * * * * 
	     
			 ! Perform direct sum (truncated in P_ell)
		   
		     PRODU(i) = sum(Tell(Ell_min(i)+1:Ell_max(i)+1)
     ;                       *Pell(Ell_min(i)+1:Ell_max(i)+1,i) )	
		     sortie(j,i) = -ci * A * sigmaA * PRODU(i) 		
				
	   end do ! (i) end of m1-m2 loop * * * * * * * * * * * * * * * * * * *

3333	   continue
cc	   sortie(j,1:Nm) = -ci * A*sigmaA * MATMUL(Tell(1:65),Pell(1:65,1:Nm)) 
c	   Laux(j)=Lmax  ! just checking Lmax values

      end do ! (j) end of A-loop ++++++++++++++++++++++++++++++++++++++++++++



! ######## All data OUTPUT tranferred to M12COKPCO and WRITE_M12_TOFILE, 
!          since here we (now) compute M12 only for positive m1-m2 values,
!		 and we want to store the M12 matrix for all m1-m2 values.


2001	continue

      END SUBROUTINE M12cokpco_write

! *********************************************************************** !
