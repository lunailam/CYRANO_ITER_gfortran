      subroutine outpow_m12std  !(isp, pic10, plan, picm10)

      implicit none
      integer isp

C      Power computation based on stored M12 matrices (../M12std)

c	 OBS: Here, the COMMON variable m12matrix(1:nspec,1:NA,1:NM) is re-used,
c		  but in slightly different context than in the M12COKPCO routine.
c		  The first index, originally planned as a species index, is used here 
c		  as the polarization index: i=1 for p=+1, i=2 for p=-1 and i=3 for p=0, 
c		  since the mode convolution is done inside the loop over species.
c		  This is OK because the minimum number of species is 3
c		  The power density of the individual species is stored in powerdens(ispec)


        include 'pardim.copy'
        include 'comgdr.copy'
	include 'complp.copy'
	include 'comgeo.copy'
	include 'comswe.copy'
	include 'comant.copy'      
	include 'commod.copy'        
        include 'dynou2.copy'
	include 'compla.copy'       
	include 'commag.copy' 
	include 'comphy.copy' 
	include 'cokpco.copy'  
	include 'compow.copy'  
	include 'comequ.copy'  		     
                              
      integer :: kr, OpenStat, ispe, im, jv, ix, km, dummyi, Nxk, irad,
     ;           kstop1, kstop2, m1, m2, shiftm, Nx, s

      real*8 :: rhoaux, rhocheck, kp(2*nmoant-1), dummy, maver, t1, t2, t3
	character(100) :: FILE_NAME
      character(4) :: auxstr


c	For power with M12 matrices
	character(100) :: GGs, GG2s
	character(2) :: charaux, charaux2
	integer :: mdiffaux(4*maxcou+2)
c	complex*16 :: m12aux(2*maxcou+1,0:2*maxcou+1)
	integer NA, NMM
	real*8 powerdens(1:nspec), sumaux, power_modes(modva2-modva1+1,1:nspec)

	integer:: ir(200), Nrad
	real*8:: raio(200), powerfrac(1:nspec), deltarho, fac
	character*20 :: FOLDER(200)
	character(4) :: panam4
ccccccccccccccccccccccccccccccccccccc
c     Non-Maxwellian Species index
c	isp = 1				! Among non-Maxwellian species
c      ispe = ispgdr(isp)  ! Among all species	
cccccccccccccccccccccccccccccccccccc


c	NA = 9
c	NMM = 9

	close(6666)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Global factor applied to VMAT in ASPLAS (L.168, 626, 1210)
        fac = glofac * twopi ** 2
        if(.not.cyl)fac = fac * ra

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c	 1) Read radii that contain M12 information (Gauss points)

	     open (UNIT = 6666, FILE = './M12std/index_radii.dat', 
     ;           STATUS = "UNKNOWN", IOSTAT = OpenStat, ACTION = "READ")
	       if (OpenStat > 0) then
		      print *, 'Error reading file: ./M12std/index_radii.dat'
	          stop
	       end if	   	
		   do k=1,200
	          read(6666,"(i8, g15.6, A20)",END=1234), ir(k), raio(k), FOLDER(k)
             end do

1234	continue

		close (6666)

c		Number of radial points (Gauss points)
		Nrad = k-1



ccccccccccccccccccccccccccccccccc Power using M12 matrices ccccccccccccccccccccccccccccccccc


c	Open file for power output
            FILE_NAME = trim(paplda)  // "/RFpower_M12std.dat" 
		open (UNIT = 9999, FILE = FILE_NAME, STATUS = "NEW",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		      end if
		  write(9999,*),"Power absorption [W/m] per species (FLR0) using M12 matrices - OUTPOW_M12std"
		  write(9999,"(3A5)"),"nspec","rad","pla"
		  write(9999,"(3i5)"), nspec, Nrad, Nrad
		  write(9999,*), "rho(m)  ", paname(1:nspec), "total"

c	Open file for power output (power per mode)
            FILE_NAME = trim(paplda)  // "/RFpower_M12std_permode.dat" 
		open (UNIT = 9966, FILE = FILE_NAME, STATUS = "NEW",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		      end if
		  write(9966,*),"Total power absorption [W/m] per mode (FLR0) using M12 matrices - OUTPOW_M12std"
		  write(9966,"(3A5)"),"Nmod","rad","pla"
		  write(9966,"(3i5)"), modva2-modva1+1, Nrad, Nrad
		  write(9966,"(g15.6,200i5)"), "rho(m)  ", (k, k = modva1, modva2)


	print*, '... Power using stored M12 matrices (standard coord.) (/M12std)'

	powerfrac(1:nspec)=0.d0

      do kr = 1, Nrad  ! Loop over Gauss points
c     =================


c	   Need to update iel !!!
         minf(iel) = modva1
         msup(iel) = modva2
         shiftm = -minf(iel)+1
	   NA = 2*(modva2-modva1)+1  ! Number of k// values 
	   NMM = klim + 1		      ! Number of m1-m2 values
      
         irad = ir(kr)
         rhoaux = raio(kr)

         call int_to_string4(nint(1000*abscis(irad)),auxstr)
	   STDFOLDER = './M12std/r' // trim(auxstr)



c        5.2) Auxiliary vector for file input: FORMAT identifier
	   call INT_TO_STRING(2*NMM+1, charaux)
	   call INT_TO_STRING(2*NMM+2, charaux2)
	        GGs  = "(G11.5," // charaux // "G20.8)"
	        GG2s = "(" // charaux2 // "G20.8)"


c	   Loop over species -------------------------------------	

	    do s = 1, nspec
		panam4 = paname(s)

		   if(my_nspgdr>0 .and. s.eq. my_ispgdr(1) .and. abscis(irad)>0.03 ) then  ! + + + + + + + + + + + 

c		   Non-Maxwellian species

c			p=+1:
			FILE_NAME = trim(STDFOLDER)  // "/NonMaxM12" // 'SPEC_lefty.dat'
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "OLD",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "READ")
	              if (OpenStat > 0) then
		              print *, 'Error reading file: ', FILE_NAME
	                  stop
				  end if
				  read (13, GGs), dummy, dummy, mdiffaux(1:2*NMM)
				  do im = 1, NA
					 read (13, GG2s),  dummyi, kp(im), m12matrix(1,im,1:NMM) ! ! i=1 (p=+1)
				  end do 
			close (13)	

ccccccccccccccccccccccc Use Maxwellian p=-1 and p=0 matrices for now cccccccccccccccccccc
c			p=-1:
			FILE_NAME = trim(STDFOLDER)  // '/M12' //  panam4 // '_right.dat'
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "OLD",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "READ")
	              if (OpenStat > 0) then
		              print *, 'Error reading file: ', FILE_NAME
	                  stop
				  end if
				  read (13, GGs), dummy, dummy, mdiffaux(1:2*NMM)
				  do im = 1, NA
					 read (13, GG2s), kp(im), dummy, m12matrix(2,im,1:NMM) ! i=2 (p=-1)
				  end do 
			close (13)	

c			p=0:
			FILE_NAME = trim(STDFOLDER)  // '/M12' //  panam4 // '_landau.dat'
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "OLD",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "READ")
	              if (OpenStat > 0) then
		              print *, 'Error reading file: ', FILE_NAME
	                  stop
				  end if
				  read (13, GGs), dummy, dummy, mdiffaux(1:2*NMM)
				  do im = 1, NA
					 read (13, GG2s), kp(im), dummy, m12matrix(3,im,1:NMM) ! i=3 (p=0)
				  end do 
			close (13)	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


		   else  ! Maxwellian species + + + + + + + + + + + + + + + 

c			p=+1:
			FILE_NAME = trim(STDFOLDER)  // '/M12' //  panam4 // '_lefty.dat'
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "OLD",
     ;              IOSTAT = OpenStat, ACTION = "READ")
	              if (OpenStat > 0) then
		              print *, 'Error reading file: ', FILE_NAME
	                  stop
				  end if
				  read (13, GGs), dummy, dummy, mdiffaux(1:2*NMM)
				  do im = 1, NA
					 read (13, GG2s), dummy, dummy, m12matrix(1,im,1:NMM) ! i=1 (p=+1)
				  end do 
			close (13)	

c			p=-1:
			FILE_NAME = trim(STDFOLDER)  // '/M12' //  panam4 // '_right.dat'
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "OLD",
     ;              IOSTAT = OpenStat, ACTION = "READ")
	              if (OpenStat > 0) then
		              print *, 'Error reading file: ', FILE_NAME
	                  stop
				  end if
				  read (13, GGs), dummy, dummy, mdiffaux(1:2*NMM)
				  do im = 1, NA
					 read (13, GG2s), kp(im), dummy, m12matrix(2,im,1:NMM) ! i=2 (p=-1)
				  end do 
			close (13)	

c			p=0:
			FILE_NAME = trim(STDFOLDER)  // '/M12' //  panam4 // '_landau.dat'
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "OLD",
     ;              IOSTAT = OpenStat, ACTION = "READ")
	              if (OpenStat > 0) then
		              print *, 'Error reading file: ', FILE_NAME
	                  stop
				  end if
				  read (13, GGs), dummy, dummy, mdiffaux(1:2*NMM)
				  do im = 1, NA
					 read (13, GG2s), kp(im), dummy, m12matrix(3,im,1:NMM) ! i=3 (p=0)
				  end do 
			close (13)	


		end if  ! Maxwellian / non-Maxwellian + + + + + + + + + + + + + + + + + + + 





		  powerdens(s)=0.d0
	      power_modes(:,s)=0.d0
	      
		  do m1 = minf(iel), msup(iel)
	      do m2 = minf(iel), msup(iel)

            maver = dfloat(m1+m2)/2.d0
            im = (m1+m2) - 2*minf(iel) +1  ! m_average index (1 .... NA)
	      km = iabs(m1-m2) +1            ! mdiff index (m12matrix starts at k=1 for m1-m2=0 to k=NMM)

            t1 = dreal(dconjg(xpmp(1,irad,m2+shiftm)) * xpmp(1,irad,m1+shiftm))  ! E+(m2).E+(m1)
		  t2 = dreal(dconjg(xpmp(3,irad,m2+shiftm)) * xpmp(3,irad,m1+shiftm))  ! E-(m2).E-(m1)
		  t3 = dreal(dconjg(xpmp(5,irad,m2+shiftm)) * xpmp(5,irad,m1+shiftm))  ! E//(m2).E//(m1)

            powerdens(s) = powerdens(s) 
     ;					 - dimag(m12matrix(1,im,km)) * t1  ! p = +1
     ;					 - dimag(m12matrix(2,im,km)) * t2  ! p = -1
     ;					 - dimag(m12matrix(3,im,km)) * t3  ! p = 0
c		  Minus sign comes from (+i) factor used in CYRANO 	


c		  Power per mode	      
            power_modes(m1+shiftm,s)= power_modes(m1+shiftm,s) 
     ;					 - dimag(m12matrix(1,im,km)) * t1  ! p = +1
     ;					 - dimag(m12matrix(2,im,km)) * t2  ! p = -1
     ;					 - dimag(m12matrix(3,im,km)) * t3  ! p = 0

		  

c		    if(jv.eq.1 .and.ix.eq.1)then
c				print*, 'E(m1=',m1, ')*E(m2=',m2,')'
c				print*, 'x M12(im=',im,',k=',km,')'
c			end if

            end do

	      end do


	 end do  ! s: species loop


c 		  WARNING: Power density is in W/m, but we want W/m2 for int(p.r.dr)=Ptot
c		           --> ADD factor 1/rho
c				   --> ADD global factor fac
c				   --> ADD factor 1/ap ??????????????????????????

					powerdens(1:nspec)=fac*powerdens(1:nspec)/rhoaux/rnorm

					power_modes(:,1:nspec) = fac*power_modes(:,1:nspec)/rhoaux/rnorm


		  write(9999,"(25g16.6)"),  rhoaux, powerdens(1:nspec), sum(powerdens(1:nspec))

		  write(9966,"(200g16.6)"), rhoaux, (sum(power_modes(k+shiftm,1:nspec)),k=modva1,modva2)
c     ;                                        , sum( sum( power_modes(:,1:nspec) ) )  


ccc	      Power fractions (rectangles integration)
	      if(kr>2)then
	         deltarho = raio(kr) - raio(kr-1)
		     powerfrac(1:nspec) = powerfrac(1:nspec) + powerdens(1:nspec)*rhoaux*deltarho
		  end if

		  if(kr.eq.Nrad)then
		  write(9999,*),  'Power fractions : int[p(r).r.dr]'
		  write(9999,"(25g16.6)"),  powerfrac(1:nspec), sum(powerfrac(1:nspec))
	      end if



      end do ! kr: Loop over Gauss points
c     ======

	close(9999)
	close(9966)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




      return
      end
