      subroutine OUTGRID

c     Builds radial and poloidal grid for all 2D OUTPUT routines
c	according to the logicals OUTGAU and OUTNPFT:
c     1) Radial
c	   OUTGAU  = TRUE : radial grid at Gauss pts. and el. boundaries (as in TABLES.f)
c	   OUTGAU  = FALSE: radial grid only at element boundaries
c     2) Poloidal
c	   OUTNPFT = TRUE : poloidal grid at NPFFT pts. (as in TABLES.f)
c	   OUTNPFT = FALSE: poloidal grid defined by variable NPLOTH (needs interpolation)
c
c	Grid variables:
c     1) for calculations: abscis(1:nabscis), polang(1:npfft)
c     2) for plotting:     absout(1:ist11),   polplo(1:nploth)
c						 absdis(1:???) - for 1D dispersion

      implicit none
	   
      include 'pardim.copy'
      include 'comgeo.copy'
      include 'comreg.copy'
      include 'comswe.copy'
      include 'complp.copy'
      include 'comfin.copy'
      include 'comin2.copy'
      include 'comphy.copy'
      include 'complo.copy'    
      include 'coequi.copy'	   
      include 'compla.copy'


	integer :: j, i, ipo, iplom1, ian, elem_indx(nele+1)
	real*8  :: polst, Rout(nabplo,npfft+1), Zout(nabplo,npfft+1)
     ;				, Rcok(nabplo,npfft+1), Zcok(nabplo,npfft+1)

C**********************************************************************

c	WARNING: Poloidal interpolation DISABLED!
c	For the moment, impose OUTNPFT = TRUE in ALL OUTPUT routines
	OUTNPFT = .true.
 
c     1) Radial GRID definition

	if(OUTGAU)then	! use previous radial grid (TABLES.f)
c					  -> Gauss pts and elem. boundaries (region bound. twice)

	   ist11 = nabsci						! Total number of points							
	   do ireg = 1, nreg					! Points per region (istp)						
	      istp(ireg) = 1 + (ngauss+1) * (ilael(ireg)-ifiel(ireg)+1)
	   end do	
	   absout(1:ist11) = abscis(1:nabsci)	! Minor radius (rho)
	   absono(1:ist11) = abscno(1:nabsci)	! y = rho/a
	   absr(1:ist11)   = eqt(1:nabsci,1,1)	! Outer major radius at Z=Zax			
	   intaplot(1:ist11) = (/ (i,i=1,nabsci) /)	! Index of radial points
	   ipla = istp(1)						! Index for plasma radius
c	   print *, nabsci, ist11
c	   print *, absr(1:ist11)
c	   print *, istp(1:nreg)
c	   print *, ipla	
c	   print *

	else	! create NEW radial grid  
c			  -> Axis + Element boundaries + Wall (plasma radius is in reg.1)

         ist11 = nele + 1						! Total number of points
	   do ireg = 1, nreg					! Points per region (istp)
	      istp(ireg) = ilael(ireg) - ifiel(ireg) + 1
	   end do
	   istp(1) = istp(1) + 1				! include r_p in 1st region
c	   Element index vector (incl. wall point)
	   elem_indx(1:nele+1) = (/ ifiabs(1:nele), nabsci /) 
	   absout(1:ist11) = abscis(elem_indx)	  ! Minor radius (rho)		
	   absono(1:ist11) = abscno(elem_indx)	  ! y = rho/a	
	   absr(1:ist11) = eqt(elem_indx,1,1)	  ! Outer major radius at Z=Zax
	   intaplot(1:ist11) = elem_indx(1:ist11) ! Index of radial points	
	   ipla = istp(1)						  ! Index of plasma radius 
	   intaplot(ipla) = intaplot(ipla) - 1    ! Take 'internal' plasma boundary

c	   print *, nabsci, ist11
c	   print *, absr(1:ist11)
c	   print *, istp(1:nreg)
c	   print *, ipla

	end if  ! (OUTGAU = true)


c	print *, intaplot(istp(2))
c	do j = 1, ist11
c	print *, intaplot(j), eqt(intaplot(j),1,1), vttab(intaplot(j),1)
c	end do
c	stop


c     2) Poloidal grid definition (NOT TESTED YET)

      if(OUTNPFT) then  ! use previous poloidal grid (TABLES.f)

		nploth = npfft			! already in Cyrano.f line 335 
		polplo(1:nploth+1) = polang(1:npfft+1)	! poloidal angle  
          iplom1 = nploth + 1		! Total number of pol. points (0-2pi)

	else  !define new poloidal GRID (needs further pol. interpolation)

	    polst = twopi / float(nploth)
		do ipo = 1, nploth
		   polplo(ipo) = (ipo - 1) * polst
		end do
		iplom1 = nploth + 1
		polplo(iplom1) = twopi

	end if


c     3) Generate new (R,Z) grid for OUTPUT (if necessary) 

c     3.1) CASE 1: Both OUTGAU and OUTNPFT are TRUE (no changes in grid)

	if (OUTGAU .and. OUTNPFT) then
	   Rout(1:ist11, 1:nploth+1) = eqt(1:nabsci, 1:npfft+1, 1) 	! R coord.
	   Zout(1:ist11, 1:nploth+1) = eqt(1:nabsci, 1:npfft+1, 2)  ! Z coord.
	   goto 1001				 ! Jump to write (R,Z) files
	end if  ! (OUTGAU & OUTNPFT = TRUE)

c     3.2) CASE 2: OUTGAU = FALSE and OUTNPFT = TRUE 
c		 (only change the RADIAL grid -> element boundaries)

	if (OUTGAU .eqv. .false. .and. OUTNPFT .eqv. .true.) then
	   Rout(1:ist11, 1:nploth+1) = eqt(intaplot(1:ist11), 1:npfft+1, 1) ! R coord.
	   Zout(1:ist11, 1:nploth+1) = eqt(intaplot(1:ist11), 1:npfft+1, 2) ! Z coord.
	   goto 1001				 ! Jump to write (R,Z) files
	end if ! (OUTGAU = F & OUTNPFT = T)

c     3.3) CASE 3: OUTGAU = TRUE and OUTNPFT = FALSE  (#### NOT TESTED)
c		 (only change the POLOIDAL grid -> interpolation)

	if(OUTGAU .eqv. .true. .and. OUTNPFT .eqv. .false.) then   

	   if(.not.updsym) then
            do i = 1, ist11	! radial loop - - - - - - 
		     intab = intaplot(i)
               call interp1(polang, polang(npfft+1), eqt(intab,1,1), nabplo,  
     ;                      npfft+1, polplo, Rout(i,1:iplom1), iplom1, 'r') ! R coord.
               call interp1(polang, polang(npfft+1), eqt(intab,1,2), nabplo, 
     ;		   		      npfft+1, polplo, Zout(i,1:iplom1), iplom1, 'r') ! Z coord.
            end do ! - - - - - - - - - - - - - - - - -
         else	! Use up-down symmetry:
            ian = nploth / 2 + 1
c            nan = iplom1 - ian + 1
            do i = 1, ist11 ! radial loop - - - - - - - - -
               intab = intaplot(i)
               call interp1(polang, polang(npfft/2+1), eqt(intab,1,1), nabplo, 
     ;                      npfft/2+1, polplo, Rout(i,1:ian), ian, 'r') ! R coord.
               do j = 1, ian
                  Rout(i,iplom1-j+1) = Rout(i,j)   ! R coord. (lower half)
               end do
               call interp1(polang, polang(npfft/2+1), eqt(intab,1,2), nabplo,  
     ;                      npfft/2+1, polplo, Zout(i,1:ian), ian, 'r') ! Z coord.
               do j = 1, ian
                  Zout(i,iplom1-j+1) = - Zout(i,j) ! Z coord. (lower half)
               end do
            end do ! (i) radial loop  - - - - - - - - - - - -
         end if  ! (.not.updsym)
	   goto 1001
	end if  ! ! (OUTGAU = T & OUTNPFT = F)

c     3.4) CASE 4: Both OUTGAU and OUTNPFT are false  (#### NOT TESTED)
c		 (change RADIAL and POLOIDAL grids )

	if(OUTGAU .eqv. .false. .and. OUTNPFT .eqv. .false.) then   

c	   First aplly new RADIAL grid (WARNING: nploth maybe not npfft!!)
	   Rout(1:ist11, 1:npfft+1) = eqt(intaplot(1:ist11), 1:npfft+1, 1) ! R coord.
	   Zout(1:ist11, 1:npfft+1) = eqt(intaplot(1:ist11), 1:npfft+1, 2) ! Z coord.

c	   Now proceed with poloidal interpolation of the new vectors (Rout,Zout)
c	   NB: The interp. solutions will overwrite the input vectors (Rout,Zout)
	   if(.not.updsym) then
            do i = 1, ist11	! radial loop - - - - - - 
		     intab = intaplot(i)
               call interp1(polang, polang(npfft+1), Rout(i,1:npfft+1), nabplo,  
     ;                      npfft+1, polplo, Rout(i,1:iplom1), iplom1, 'r') ! R coord.
               call interp1(polang, polang(npfft+1), Zout(i,1:npfft+1), nabplo, 
     ;		   		      npfft+1, polplo, Zout(i,1:iplom1), iplom1, 'r') ! Z coord.
            end do ! - - - - - - - - - - - - - - - - -
         else	! Use up-down symmetry:
            ian = nploth / 2 + 1
c            nan = iplom1 - ian + 1
            do i = 1, ist11 ! radial loop - - - - - - - - -
               intab = intaplot(i)
               call interp1(polang,polang(npfft/2+1),Rout(i,1:npfft/2+1),nabplo, 
     ;                      npfft/2+1, polplo, Rout(i,1:ian), ian, 'r') ! R coord.
               do j = 1, ian
                  Rout(i,iplom1-j+1) = Rout(i,j)   ! R coord. (lower half)
               end do
               call interp1(polang,polang(npfft/2+1),Zout(i,1:npfft/2+1),nabplo,  
     ;                      npfft/2+1, polplo, Zout(i,1:ian), ian, 'r') ! Z coord.
               do j = 1, ian
                  Zout(i,iplom1-j+1) = - Zout(i,j) ! Z coord. (lower half)
               end do
            end do ! (i) - - - - - - - - - - - - - - - - - - - 
         end if  ! (.not.updsym)

	end if  ! (OUTGAU & OUTNPFT = FALSE)


1001	continue

c	----------------------------- OUTPUT --------------------------------

c     4) Write (R,Z) GRID for plot routines to files			

c     4.1) R(rho,theta) coordinate
	open (UNIT = 777, FILE = paplda // 'Rcyr.dat', 
     ;	  STATUS = "REPLACE", ACTION = "WRITE")
		  write(777,*), 'R coordinate'
		  write(777,*), ' rad ' , ' pol ', ' pla ' 
		  write(777,"(i5, i5, i5)"), ist11, nploth+1, ipla 
		  do i = 1, ist11	
		     write(777,2222)(Rout(i,ipo), ipo = 1, nploth+1)
		  end do
	close(777)

c    	4.2) Z(rho,theta) coordinate
	open (UNIT = 777, FILE = paplda // 'Zcyr.dat', 
     ;	  STATUS = "REPLACE", ACTION = "WRITE")
		  write(777,*), 'Z coordinate'
		  write(777,*), ' rad ' , ' pol ', ' pla ' 
		  write(777,"(i5, i5, i5)"), ist11, nploth+1, ipla 
		  do i = 1, ist11	
			 write(777,2222)(Zout(i,ipo), ipo = 1, nploth+1)
		  end do
	close(777)


cERN_24JAN05: 
c	Special case for COKPCO: (R,Z) in thetabar mesh	

	if(cokpco)then ! Interpolate (R,Z) to thetabar mesh (overwrite old values)

		do i = 1, ist11	! radial loop - - - - - - 
		   intab = intaplot(i)
		   call interp2(ckt(i,1:nploth+1,3), 1, Rout(i,1:nploth+1), 1, nploth+1, 
     ;		            polang, Rcok(i,1:nploth+1), nploth+1) ! R coord.
		   call interp2(ckt(i,1:nploth+1,3), 1, Zout(i,1:nploth+1), 1, nploth+1, 
     ;		            polang, Zcok(i,1:nploth+1), nploth+1) ! Z coord.             
          end do ! - - - - - - - - - - - - - - - - -

c     4.3) R(rho,theta) coordinate in thetabar mesh
	open (UNIT = 777, FILE = paplda // 'Rcok.dat', 
     ;	  STATUS = "REPLACE", ACTION = "WRITE")
		  write(777,*), 'R coordinate'
		  write(777,*), ' rad ' , ' pol ', ' pla ' 
		  write(777,"(i5, i5, i5)"), ist11, nploth+1, ipla 
		  do i = 1, ist11	
		     write(777,2222)(Rcok(i,ipo), ipo = 1, nploth+1)
		  end do
	close(777)

c    	4.4) Z(rho,theta) coordinate in thetabar mesh
	open (UNIT = 777, FILE = paplda // 'Zcok.dat', 
     ;	  STATUS = "REPLACE", ACTION = "WRITE")
		  write(777,*), 'Z coordinate'
		  write(777,*), ' rad ' , ' pol ', ' pla ' 
		  write(777,"(i5, i5, i5)"), ist11, nploth+1, ipla 
		  do i = 1, ist11	
			 write(777,2222)(Zcok(i,ipo), ipo = 1, nploth+1)
		  end do
	close(777)

	end if	! (cokpco = TRUE)

      return

 2222 format(200(g14.6)) 

      end subroutine OUTGRID

