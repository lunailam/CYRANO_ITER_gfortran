      subroutine read_QLFPheader
      
      implicit none

      include 'pardim.copy'
      include 'comfic.copy'
      include 'comfin.copy'
      include 'compla.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'comin2.copy'
      include 'comphy.copy'
      include 'comreg.copy'	! ERN
	include 'cokpco.copy'	! ERN
      include 'complp.copy'
      include 'coequi.copy'

c      Read QLFP data and interpolate to Cyrano v-grid
c	 Interpolation to Cyrano radial grid is done afterwards 
c	 in interp_QLFP (ASPLAS) at each radial eleme

c	 LAst modified : 11/02/06     

c	 WARNING: Array sizes defined locally, should UPDATE later
c	 WARNING: Maximum number of QLFPfiles = maxrdr (pardim.copy) ~ 50

      integer kr, jv, ix, isp, i, j, igx, igv, nradqlfp, 
     ; nradqlfp_cyr, nabsv_cyr, N1, N2
	real*8 :: dummy, vcyr(maxnev,maxrdr), v_min, v_max , v_12 
	character(100) :: QLFPDATA(maxrdr),  RUNID, WHEN
	integer :: OpenStat
	real*8 :: F0cyr(maxnex,maxnev,maxrdr,4), vel, aux1(maxrdr), aux2(maxrdr)
        integer find_nearest

	isp = 1
c	ispa = ispgdr(isp)



c	 1) Read QLFP header file (Nsurfaces, data file NAMES) -------------------

      
	print*, 'Reading QLFP data in ', trim(QLFPDIR) // 'QLFP_header.dat', ' ...'
	write(nofile,*)
	write(nofile,*)'Reading QLFP data in ', trim(QLFPDIR) // 'QLFP_header.dat', ' ...'


	open (UNIT = 223, FILE = trim(QLFPDIR) // 'QLFP_header.dat', 
     ;      STATUS = "OLD", IOSTAT = OpenStat, ACTION = "READ")
	            if (OpenStat > 0) then
		            print *, 'Error reading file: ', trim(QLFPDIR) // 'QLFP_header.dat'
	                stop
		        end if

		read(223,"(A45)"), RUNID
		read(223,"(A30)"), WHEN
		read(223,*), nradqlfp

		do kr = 1, nradqlfp   ! Radial loop =============================================	
c		   read(223,"(f16.6, A21)") raqlfp(kr), QLFPDATA(kr)
		   read(223,*) raqlfp(kr), QLFPDATA(kr)	
                   write(nofile,*)raqlfp(kr), trim(QLFPDATA(kr)), N1+N2+2 	
                end do
	close(223)


	
c	----------------------------------------------------------------------------


c	 2) Loop over QLFP mag surfaces --------------------------------------------

	do kr = 1, nradqlfp	

      print *, kr
c	    2.1) Read QLFP data files and store results in vgta, xgta and F0 

		open (UNIT = 224, FILE = trim(QLFPDIR) // trim(QLFPDATA(kr)), 
     ;          STATUS = "OLD", IOSTAT = OpenStat, ACTION = "READ")
	          if (OpenStat > 0) then
		          print *, 'Error reading file: ', QLFPDATA(kr)
	              stop
		      end if
			  read(224,*), dummy, aux1(kr), aux2(kr)  ! rho, n, vT
			  read(224,*), dummy				  ! xnsep
			  read(224,*), nabsx, nabsv
			  ! v-grid
			  read(224,"(A18, 100(g18.6))"), dummy, vgta(1:nabsv,kr)
			  ! fo
			  do ix = 1, nabsx
			     read(224,"(100(g18.6))"), xgta(ix,kr), F0(ix,1:nabsv,kr,1,isp)
			  end do
			  read(224,*) 
			  ! dfo/dv
			  do ix = 1, nabsx
			     read(224,"(100(g18.6))"), xgta(ix,kr), F0(ix,1:nabsv,kr,2,isp)
			  end do
			  read(224,*) 
			  ! dfo/dx
			  do ix = 1, nabsx
			     read(224,"(100(g18.6))"), xgta(ix,kr), F0(ix,1:nabsv,kr,3,isp)
			  end do

		close(224)




c	    2.2) NEW(Jan06): Interpolate the QLFP velocity grid to CYRANO v-grid
c		     The CYRANO v-grid is defined HERE (same for ALL QLFP radii)	
c			 (Repeat v-values at region boundaries for consist. with old runs)

		v_min = 1.d3			! minimum v
		v_12  = 2.5d6			! v-region interface
		v_max = 8.0d6			! maximum v
		N1 = 15					! N elements in region 1
		N2 = 65					! N elements in region 2
		nabsv_cyr = N1 + N2 + 2 ! Total number of v-points
	
		! Build v-grid
		call cuteq(v_min, v_12, N1, vcyr(1,kr))
		call cuteq(v_12, v_max, N2, vcyr(2+N1,kr))


c		2.3) Linear interp. of fo and derivatives to CYRANO v-grid

		do jv = 1, nabsv_cyr 

		    vel = vcyr(jv,kr)

			if(vel.le.vgta(nabsv,kr))then  ! interpolate

				do ix = 1, nabsx

c				call my_interpolation(nabsv, vgta(1:nabsv,kr), F0(ix,1:nabsv,kr,1,isp), 
c     ;				                  1,     vel,              F0cyr(ix,jv,kr,1)   )
c				call my_interpolation(nabsv, vgta(1:nabsv,kr), F0(ix,1:nabsv,kr,2,isp), 
c     ;				                  1,     vel,              F0cyr(ix,jv,kr,2)   )
c				call my_interpolation(nabsv, vgta(1:nabsv,kr), F0(ix,1:nabsv,kr,3,isp), 
c     ;				                  1,     vel,              F0cyr(ix,jv,kr,3)   )

c				   Quadratic interp.
		           call interp2(vgta(1:nabsv,kr), 1, F0(ix,1:nabsv,kr,1,isp), 1, nabsv, 
     ;                            vel,                 F0cyr(ix,jv,kr,1), 1 )
		           call interp2(vgta(1:nabsv,kr), 1, F0(ix,1:nabsv,kr,2,isp), 1, nabsv, 
     ;                            vel,                 F0cyr(ix,jv,kr,2), 1 )
		           call interp2(vgta(1:nabsv,kr), 1, F0(ix,1:nabsv,kr,3,isp), 1, nabsv, 
     ;                            vel,                 F0cyr(ix,jv,kr,3), 1 )


				end do  ! ix

			else  ! velocity point outside QLFP data range

				F0cyr(1:nabsx,jv,kr,1:3) = 0.d0

			end if

		end do  ! jv


c	Update global variables with new v-GRID (for other routines)

	F0(1:nabsx,1:nabsv_cyr,kr,1:3,isp) = F0cyr(1:nabsx,1:nabsv_cyr,kr,1:3)
	vgta(1:nabsv_cyr,kr) = vcyr(1:nabsv_cyr,kr)
	nabsv = nabsv_cyr


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	end do	! kr: radial QLFP loop ------------------------------------------

      print*, '.... OK!'

cERN	Output added 23/06/05

			open (UNIT = 77, FILE = paplda // 'QLFP_distribution.dat', STATUS = "REPLACE",
cc     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', "QLFP_distribution.dat"
	                  stop
				  end if
	    
	  				  write(77,*),"Non-Maxwellian dist. interpolated to CYRANO v-grid"
	  				  write(77,*),"(read_QLFPheader.f)"
					  write(77,*), nradqlfp, nabsv, nabsx

				  do kr = 1, nradqlfp	! ---------------------------------

					  write(77,*)
					  write(77,*) raqlfp(kr), aux1(kr), aux2(kr)

						write(77,*)'F0(xn:lines,v:cols)'
						write(77,2222)'*************',vgta(1:nabsv,kr)										
						do ix = 1,nabsx
						write(77,2222) xgta(ix,kr),F0(ix,1:nabsv,kr,1,isp)
						end do

						write(77,*)'dF0/dv(xn:lines,v:cols)'				
						do ix = 1,nabsx
						write(77,2222) xgta(ix,kr),F0(ix,1:nabsv,kr,2,isp)
						end do

						write(77,*)'dF0/dx(xn:lines,v:cols)'				
						do ix = 1,nabsx
						write(77,2222) xgta(ix,kr),F0(ix,1:nabsv,kr,3,isp)
						end do

					end do	! kr = 1, nradqlfp  -----------------------

	close (77)
c
c




	nraddr = nradqlfp


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cNEW  Determine the index of Cyrano radial points that 'match' the
c     surfaces given by BATCH (needed for generating the OUTPUT coefs. for Dirk)
cNEW  New COMMON variable intaqlfp
c     Careful: consider only Gauss points! (because DIRESP is only called at Gauss pts)
        do j = 1, nradqlfp
         intaqlfp(j) = find_nearest(nabsci, abscis, raqlfp(j))
         ! Exclude element boundaries
           if(any(intaqlfp(j).eq.ifiabs(1:nele)))then
            intaqlfp(j) = intaqlfp(j) - 1
           end if
           ! Correct last element (plasma edge)
           if(intaqlfp(j).eq.i_rp)then
            intaqlfp(j) = intaqlfp(j) - 1
           end if
c         print*, raqlfp(j), abscis(intaqlfp(j))
      end do
c      stop

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      return

        

c 1000 format(1h , 4(g14.6,2x))
 2222 format(1000(g18.6))        

      end
