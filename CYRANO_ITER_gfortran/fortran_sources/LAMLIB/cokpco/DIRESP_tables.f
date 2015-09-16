      SUBROUTINE DIRESP_tables(LX,xx)

      IMPLICIT NONE

c     Build tables of several x-dependent quantities for 
c     given x-mesh : xx = [0...Xmax]

c     INPUT  : LX (length of xx), xx (vector with x values)

c     OUPUT (COMMONS):
c       - aatab(x)    : small 'a' from paper
c       - a2mb2(x)    : (a^2-b^2) from paper
c       - kefftab(x)  : kp or kt (passing / trapped)
c       - cosXMtab(x) : cosXM from paper 
c       - Ktab(x) and Etab(x)   : Elliptic integrals (Using IMSL routines)
c       - N0tab(x) and N1tab(x) : N0 and N1 coefficients

c     OBS: Using IMSL elliptic integrals with argument m = k^2 (k: modulus) 																					
      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'cokpco.copy'
      include 'comgeo.copy'	


	 
!	Input	                                                                 
	integer, intent(in) :: LX
	real*8, intent(in)  :: xx(LX)

     
      integer j, k, OpenStat, Lmax
	integer, parameter :: Nchi = 25	! Number of chi-harmonics (ell = 0... Nchi-1) 
      character(100) :: FILE_NAME
	real*8 :: xmax, xsep, k2, DELK, DELE
        real*8 ellfc, ellec
	external DELK, DELE


c	Maximum chi-harmonic index (ell=Lmax)
	Lmax = Nchi-1

c	Radius dependent quantities
      xnsep = - 2.d0 * delb(intab) / bmax(intab) 
      xmax = B0 / bmin(intab)
      xsep = (1.d0-dabs(xnsep))*xmax


c     x - loop ===============================================

	do j = 1, LX   

c        xtab(j) = dfloat(j-1)/dfloat(LX-1)*xmax;
         xtab(j) = xx(j)         

         if(xtab(j) .eq. 0.0d0)then     ! Special case: x = 0
c        ----------------------------------------------------        
            aatab(j)    = 1.0d0
            a2mb2(j)    = 1.0d-10       ! actually 0.0d0   
		  cosXMtab(j) = -1.0d10       ! actually -Inf 
            kefftab(j)  = 0.0d0         ! actually 0.0d0
            Ktab(j)     = pi/2
	      Etab(j)     = pi/2
	      N0tab(j)    = pi
		  N1tab(j)    = 0        


         elseif(xtab(j) .eq. xmax)then  ! Special case: x = Xmax
c        -------------------------------------------------------   
            aatab(j)    = 1.0d-10       ! actually 0.0d0
            a2mb2(j)    = 2*xmax*delb(intab)/B0
		  cosXMtab(j) = 1.0d0
            kefftab(j)  = 0.0d0         ! actually 0.0d0
            Ktab(j)     = pi/2
	      Etab(j)     = pi/2
		  N0tab(j)    = pi/sqrt(a2mb2(j)) 
		  N1tab(j)    = pi/sqrt(a2mb2(j)) 
		  ! Recursion for Nell series (trapped -> forwards)
		  NL(0,j) = N0tab(j)
		  NL(1,j) = N1tab(j)
		  do k = 2, Lmax	! (Lmax = Nchi-1)
		     NL(k,j) =  dfloat(2*k-2)    / (dfloat(k)-0.5d0) * cosXMtab(j) * NL(k-1,j) 
     ;	  		     - (dfloat(k)-1.5d0) / (dfloat(k)-0.5d0) * NL(k-2,j)
		  end do  


         elseif(xtab(j) .eq. xsep)then  ! Special case: x = xsep
c        -------------------------------------------------------   
            aatab(j)    = sqrt(1.0d0-xtab(j)*bmin(intab)/B0)
            a2mb2(j)    = 1.0d0-xtab(j)*bmin(intab)/B0 
		  cosXMtab(j) = -1.0d0
            kefftab(j)  = 1.0d0-1.d-8  ! actually 1.d0, but then problems with PI3
            Ktab(j)     = 10.6d0       ! actually +Inf, but K(1-1d-8)=10.6
	      Etab(j)     = 1.0d0
		  N0tab(j)    = 2/aatab(j) * Ktab(j)
	      N1tab(j)    = 2/aatab(j) * ( 2*Etab(j)-Ktab(j) )
	  

         else   ! General case : 0< x <xsep(pass) and xsep <x <Xmax(trap)
c        ----------------------------------------------------------------               
            
		  aatab(j) = sqrt(1.0d0-xtab(j)*bmin(intab)/B0)
            a2mb2(j) = 2*xtab(j)*delb(intab)/B0
		  cosXMtab(j) = 1/delb(intab)*(bbar(intab)-B0/xtab(j))

            if(xtab(j) < xsep)then ! >>>>>>>>>>>>>>>>>>>  passing orbits
               
               kefftab(j) = sqrt(a2mb2(j)) / aatab(j)
               k2 = kefftab(j)**2 
	         ! IMSL elliptic integrals (m=k^2)   
cJAC               Ktab(j) = DELK(k2)
cJAC	         Etab(j) = DELE(k2) 
c     Carlson.f routines (check argument is corect k2 or k  xxxxxxxxxxxxxxxxxxxxxxxx)
               Ktab(j) = ellfc(k2)
               Etab(j) = ellec(k2)

		     N0tab(j) = 2/aatab(j) * Ktab(j)
		     N1tab(j) = 2/aatab(j) * ( Ktab(j) + 2/k2*(Etab(j)-Ktab(j)) )

cERN             Avoid float overflow for large cosXM 			 
			 if(abs(N1tab(j))<1.d-2)then
						NL(2:Lmax,j) = 0.0d0
						goto 666
				  end if  

			 ! Recursion for Nell series (passing = backwards)
			 ! (Lmax = Nchi-1, ell = 0 ... Lmax)
cERN                  TO BE REVISED : Lmax = Lmax(cosXM) !!!!!!!!!!!!!!!!
			   NL(Lmax, j) = 1.0d0
			   NL(Lmax-1, j) = 1.d0
			   do k = Lmax-2, 0, -1
				  
				  NL(k,j) =  dfloat(2*k+2)    / (dfloat(k)+0.5d0) * cosXMtab(j) * NL(k+1,j) 
     ;					  - (dfloat(k)+1.5d0) / (dfloat(k)+0.5d0) * NL(k+2,j)
			   end do
			   ! Normalize to No	
				do k = 2, Nchi-1
				   NL(k,j) = NL(k,j) * (N0tab(j) / NL(0,j))
				end do
666				continue
				NL(0,j) = N0tab(j)
				NL(1,j) = N1tab(j)


                                                 
            else                 ! >>>>>>>>>>>>>>>>>>>>>>>  trapped orbits    
               
               kefftab(j) = aatab(j) / sqrt(a2mb2(j))
		     k2 = kefftab(j)**2
	         ! IMSL elliptic integrals (m=k^2)   
cJAC               Ktab(j) = DELK(k2)
cJAC	         Etab(j) = DELE(k2)    
c     Carlson.f routines (check argument is corect k2 or k  xxxxxxxxxxxxxxxxxxxxxxxx)
               Ktab(j) = ellfc(k2)
               Etab(j) = ellec(k2)

		     N0tab(j) = 2/aatab(j) * kefftab(j) * Ktab(j)
		     N1tab(j) = 2/aatab(j) * kefftab(j) * ( 2*Etab(j)-Ktab(j) )
			 ! Recursion for Nell series (trapped = forwards)
			   NL(0,j) = N0tab(j)
			   NL(1,j) = N1tab(j)
			   do k = 2, Lmax	! (Lmax = Nchi-1)
				  
				  NL(k,j) =  dfloat(2*k-2)    / (dfloat(k)-0.5d0) * cosXMtab(j) * NL(k-1,j) 
     ;					  - (dfloat(k)-1.5d0) / (dfloat(k)-0.5d0) * NL(k-2,j)
			   
			   end do      
				           
            end if


         end if ! Special cases (x=0, x=Xmax or x=xsep) ----------------------
             
		      
	end do    ! x-loop


c     =========================================================================


c     Write data to file (only for debugging) ----------------------

	if (WRITE_OUTPUT) then
          FILE_NAME = COKFOLDER  // "/diresp_table.dat" 
	    open (UNIT = 99, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;  	    BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;          IOSTAT = OpenStat, ACTION = "WRITE")
	          if (OpenStat > 0) then
		        print *, 'Error writing file: ', FILE_NAME
	              stop
		    end if	  
		    write(99,*),"Elliptic integ. 1st and 2nd kind (IMSL)"
		    write(99,*),"xsep  ", "  xmax" 
                write(99,*), xsep, xmax
                write(99,*),lxg
                do j = 1,lxg 
                   write(99,"(9G16.8)"), xtab(j), aatab(j), a2mb2(j), kefftab(j)
     ;                                 , Ktab(j), Etab(j)
     ;					   , cosXMtab(j), N0tab(j), N1tab(j) 
	          end do
           close(99)


          FILE_NAME = COKFOLDER  // "/Nell_table.dat" 
	    open (UNIT = 99, FILE = FILE_NAME, STATUS = "REPLACE",
     ;          IOSTAT = OpenStat, ACTION = "WRITE")
	          if (OpenStat > 0) then
		        print *, 'Error writing file: ', FILE_NAME
	            stop
		      end if	  
		      write(99,*),"Table of Nell coefficients, lines: x, cols: ell = [0 ... Nchi-1]"
		      write(99,*),"rho=", abscis(intab) 
			  write(99,*),"xsep=", xsep 
                write(99,*),"xmax=", xmax		      
			  write(99,*),"Nx ", "  Nchi" 
                write(99,*), lxg, Nchi
                do j = 1,lxg 
                   write(99,"(200G16.8)"), xtab(j), NL(0:Lmax,j) 
	          end do
           close(99)


	end if


      return
 
      END SUBROUTINE DIRESP_tables
