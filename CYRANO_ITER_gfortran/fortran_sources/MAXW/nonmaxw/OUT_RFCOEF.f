      subroutine out_rfcoef(isp)  !(isp, pic10, plan, picm10)

      implicit none
      integer isp
C      double precision pic10, plan, picm10

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Still have problems with COMMON activex(1:NA,jv,ix,-NM:NM) definition
c	Size is to big > 2GB !!!!!!!!!
c	For the moment, activex(40-1,1:maxvg,1:501,-20:20) (COMGDR.copy)
c	There are 2 solutions (to be implemented):
c		1) Try DIRECT access to read (k//) data inside the (x,v) loop
c		2) Re-do the DIRESP computation at the QLFP radii
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


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
     ;           kstop1, kstop2, m1, m2, shiftm, Nx

      real*8 :: rhoaux, rhocheck, kp(2*nmoant-1), v(nabsv), 
     ;          dummy, xxn(maxnex+150), t
	real*8 :: activ(1:maxnex+150, 1:klim+1), RFcoef(nabsv,501)

	character(100) :: FILE_NAME
      character(4) :: auxstr
      real*8 :: omca, xtir1, xmax, dx, xnaux(501)
	real*8 :: maver
	real*8 :: powerdens, dfdv, dfdx, dfdv_prev, dfdx_prev, x, xprev, DELV,v_2,v_3
	

c	For power with M12 matrices
	character(100) :: GGs, GG2s
	character(2) :: charaux, charaux2
	integer :: mdiffaux(4*maxcou+2)
	complex*16 :: m12aux(2*maxcou+1,0:2*maxcou+1)
	integer NA, NMM, kmaux
	real*8 powerdens2, sumaux, m12test(2*maxcou+1,1:klim+1), vel
	real*8 FAC1, FAC2, FACTOR, powerdens3, fac, sumaux2(klim+1), NOVO


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Global factor applied to VMAT in ASPLAS (L.168, 626, 1210)
	fac = glofac * twopi ** 2 * ra
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Non-Maxwellian Species index (only one species for the moment)
	isp = 1				! Among non-Maxwellian species 
      ispe = ispgdr(isp)  ! Among all species	(usually = 2)

c	Number of modes
	NA = 2*(modva2-modva1)+1  ! Number of m_average values 
	NMM = klim + 1		      ! Number of m1-m2 values


c     Define the xn-grid for output (max. 501)
       Nx = 501
	 dx = 2.d0/(Nx-1)
       do k = 1, Nx
          xnaux(k) = -1.d0 + dfloat(k-1) * dx   
       end do


c	Open file for power output
            FILE_NAME = trim(paplda)  // "/RFpower_coefs.dat" 
		open (UNIT = 9999, FILE = FILE_NAME, STATUS = "NEW",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		      end if


cc	Open file for M12check
            FILE_NAME = trim(paplda)  // "M12check.dat" 
		open (UNIT = 1231, FILE = FILE_NAME, STATUS = "NEW",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		      end if


cc	Open file for power output based on local M2 matrix
            FILE_NAME = trim(paplda)  // "/RFpower_m12check.dat" 
		  open (UNIT = 8888, FILE = FILE_NAME, STATUS = "NEW",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		      end if


c	Open file RFcoef HEADER
            FILE_NAME = trim(paplda)  // "RFcoef_header.dat" 
		open (UNIT = 8001, FILE = FILE_NAME, STATUS = "NEW",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		      end if
            write(8001,"(I3)")nraddr


c     Loop over BATCH surfaces

      do kr = 1, nraddr
c     =================
      
      irad = intaqlfp(kr)

      xnsep = - 2.d0 * delb(irad) / bmax(irad) 
	xmax = B0 / bmin(irad)
      omca = qom(2) * B0
      xtir1 = omca / omegag
c       xaux(1:Nx) = (1.d0 - dabs(xaux(1:Nx))) * xmax 

      rhoaux = abscis(irad)

	FAC1 = dPsidr_n(irad) * abscis(irad) / (4.d0*pi*B0) 
	FAC2 = (2.d0*pi)**2 * qom(2)**2 * MH*AMASS(2)/2.d0

cccccccc Same as DIRESP_standard
c	Normalization factor
	NOVO = PI * (EEL*ZCH(ispe))**2 * dPsidr_n(intab) * abscis(intab)
     ;       / (2.d0*MH*AMASS(ispe) * B0)
c	Divide by ASPLAS factor
c      NOVO = NOVO/fac
cccccccccccccccccccccccccccccc


      call int_to_string4(nint(1000*abscis(irad)),auxstr)
      
            FILE_NAME = paplda  // "/RFtemp_" // auxstr // ".dat" 
		open (UNIT = 9090, FILE = FILE_NAME, STATUS = "OLD",
     ;            IOSTAT = OpenStat, ACTION = "READ")
	            if (OpenStat > 0) then
		          print *, 'Error reading file: ', FILE_NAME
	                stop
		      end if	  
c		read(9090,*),dummychar
            read(9090,*)rhocheck
	      read(9090,*)
	      read(9090,*), dummyi, dummyi
	      read(9090,*)

            do im = 1, 2*nmoant-1  ! loop over k// ------------------
               read(9090,*)kp(im) 
               
               do jv = 1,nabsv  ! v-loop - - - - - - - - 
                  read(9090,*)
                  read(9090,*) v(jv)
                  read(9090,*)
	            read(9090,*) Nxk
               	
                  do ix = 1, Nxk
	               read(9090,"(2G16.8)") xxn(ix)	! (k//,v) dependence 'ommitted'
	               do km = 1, klim+1
	               read(9090,"(I10,G20.10)") dummyi, activ(ix,km)
                     end do
	            end do

c     --------------------------------------------------------------------------------               
               ! Interpolate from xxn(k//,v) to constant xn-grid
			 ! (pad with zeros when needed _w0) 
		      do km = -klim, klim
	             kmaux = iabs(km)+1
                  do k = 1, Nx
                     call my_interpolation(Nxk, xxn(1:Nxk), activ(1:Nxk,kmaux),
     ;                                     1,     xnaux(k)          , activex(im,jv,k,km) )

                  end do
	            end do
c     ----------------------------------------------------------------------------               

               end do ! jv: v-loop - - - - - - - - 
                 
            end do  ! im : loop over k// ----------------------------
            
                   
      close(9090)



c            FILE_NAME = paplda  // "/RFcheck_" // auxstr // ".dat" 
c		open (UNIT = 9090, FILE = FILE_NAME, STATUS = "NEW",
c     ;            IOSTAT = OpenStat, ACTION = "WRITE")
c	            if (OpenStat > 0) then
c		          print *, 'Error writing file: ', FILE_NAME
c	                stop
c		      end if	  
c
c            write(9090,*),rhocheck
c	      write(9090,*)
c	      write(9090,*), 2*nmoant-1, nabsv
c	      write(9090,*)
c
c           do im = 1, 2*nmoant-1  ! loop over k// ------------------
c               write(9090,*), kp(im) 
c               
c               do jv = 1,nabsv
c                  write(9090,*)
c                  write(9090,*), v(jv)
c                  write(9090,*)
c	            write(9090,*), Nx
               	
c                  do ix = 1, Nx !Nxk
c	               write(9090,"(2G16.8)") xnaux(ix) !xx(im,jv,ix), xxn(im,jv,ix)
c	               do km = -klim, klim
c	               write(9090,"(I10,G16.8)"), km, activex(im,jv,ix,km)
c                     end do
c	            end do
               
c               end do
                 
c            end do  ! im : loop over k// ----------------------------
            
                   
c      close(9090)








ccccccccccccccccccccccccccccccccc	Compute RF coefs  cccccccccccccccccccccccccccccccccc

c	Need to update iel !!!
      minf(iel) = modva1
      msup(iel) = modva2
      shiftm = -minf(iel)+1

      do jv = 1,nabsv
         
	   FACTOR = FAC1*FAC2*v(jv)  ! v(jv)**5 / v(jv)**4 
	         
         do ix = 1, Nx

            RFcoef(jv,ix) = 0.d0

c		  Double loop over poloidal modes - - - - - - 
            do m1 = minf(iel), msup(iel)
	      do m2 = minf(iel), msup(iel)

            maver = dfloat(m1+m2)/2.d0
            im = (m1+m2) - 2*minf(iel) +1  ! m_average index
	      km = m1-m2                      ! mdiff index

c		  Only p=+1 here [ E+(m2).D21.E+(m1) ]
            t = dreal(dconjg(xpmp(1,irad,m2+shiftm)) * xpmp(1,irad,m1+shiftm))

            RFcoef(jv,ix) = RFcoef(jv,ix) + activex(im,jv,ix,km) * t / FACTOR 
c		  Add global factor 'fac' from VMAT		     

c		    if(jv.eq.1 .and.ix.eq.1)then
c				print*, 'E(m2=',m2, ')*E(m1=',m1,')'
c				print*, 'x M12(im=',im,',k=',km,')'
c			end if

            end do
	      end do
c		  - - - - - - - - - - - - - - - - - - - - - - -	

         end do
	end do


c	Save RF coefs. to file

            FILE_NAME = trim(paplda)  // "RFgood_" // auxstr // ".dat" 
	      write(8001,"(F12.6,A20)"),rhocheck, "'RFgood_" // auxstr // ".dat'" 
		  open (UNIT = 9090, FILE = FILE_NAME, STATUS = "NEW",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		        end if
		      write(9090,*), "'Hydrogen'", MH*AMASS(2), EEL*ZCH(2)
			  write(9090,*), rhocheck, dentab(irad,2)*1.d20, vttab(irad,2)  ! rho, n, vT
			  write(9090,*), xnsep, xtir1,xmax			    ! xnsep
			  write(9090,*), Nx, nabsv


c	   RF coefs:

         write(9090,"(200g18.6)"), '************',v(1:nabsv)      

	   write(9090,*)
	   write(9090,*), 'J=dPsi/dr.v3/(4.pi.Bo)'
	   write(9090,"(200(g18.6))"), (FAC1*v(jv)**3, jv=1,nabsv)
         

	   write(9090,*)
	   write(9090,*), '-1/v'
	   write(9090,"(200(g18.6))"), (-1.d0/v(jv), jv=1,nabsv)

	   write(9090,*)
	   write(9090,*), '+2/v2.(x-wc/w)'
	   do ix = 1, Nx
	      x = (1.d0 - dabs(xnaux(ix)) ) * xmax
	      write(9090,"(200(g18.6))"), xnaux(ix), (2.d0/v(jv)**2*(x-xtir1), jv=1,nabsv)
	   end do

	   write(9090,*)
	   write(9090,*), 'RFoef=a.sum[Em1.D21.Em2]' 
	   do ix = 1, Nx
	      write(9090,"(200(g18.6))"), xnaux(ix), (RFcoef(jv,ix), jv=1,nabsv)
	   end do

	close(9090)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccc CHECK M12 cccccccccccccccccccccccccccccccccccccccc

cc	Build F0i from F0(radqlfp)
c      call interp_QLFP_noout(irad)

c	write(1231,*), abscis(irad)


cc	Re-Compute M12(mbar,m1-m2) matrix, to CHECK interp, etc... !!!!!!!!!!!!!!
   

c      do im = 1, 2*nmoant-1 
c      do km = -klim,klim

c	   sumaux = 0.d0

c	   do jv = 1,nabsv
        
c            vel = v(jv)

c            do ix = 1, Nx
c
cc	      Interpolate F0i to the given xnaux(ix) point 
c	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 2, isp), 1, xnaux(ix), dfdv)
c	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 3, isp), 1, xnaux(ix), dfdx)
c		       x = (1.d0 - dabs(xnaux(ix)) ) * xmax          
c                dfdx = dfdx *xmax  ! BATCH gives df/dxn, not df/dx
c
c		  Interpolate F0i to the previous xnaux(ix-1) point
c            if(ix .gt. 1)then
c	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 2, isp), 1, xnaux(ix-1), dfdv_prev)
c	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 3, isp), 1, xnaux(ix-1), dfdx_prev)
c                 xprev = (1.d0 - dabs(xnaux(ix-1)) ) * xmax
c		       dfdx_prev = dfdx_prev *xmax
c           end if

c	      TRAPEZE integration in x, RECTANGLES in v

c		  v_3 = v(jv)*v(jv)*v(jv)
c	      v_2 = v(jv)*v(jv)
c		  if( jv < lvg)then
c			  DELV = v(jv+1) - v(jv)
c		  else
c			  DELV = v(jv) - v(jv-1)
c		  end if
c
c       if(ix .gt. 1)then
c            sumaux = sumaux + 0.5d0
c     ;      * ( activex(im,jv,ix,km)  * (vel**3*dfdv -2.d0*vel**2*(x-xtir1)*dfdx) 
c     ;        + activex(im,jv,ix-1,km)* (vel**3*dfdv_prev-2.d0*vel**2*(xprev-xtir1)*dfdx_prev) ) 
c     ;      * (xnaux(ix)-xnaux(ix-1))*xmax * DELV 

c	  end if

c            end do
c	   end do

c	m12test(im,iabs(km)+1) = sumaux

      
c	end do
c	end do

c      do im = 1, 2*nmoant-1  ! loop over k// ------------------
c        write(1231,"(I10,20G18.6)"), im, kp(im), m12test(im,1:klim+1) 
c      end do
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cc	Build F0i from F0(radqlfp)
      call interp_QLFP_noout(irad)

	write(1231,*), rhoaux
     

      call int_to_string4(nint(1000*abscis(irad)),auxstr)
      
            FILE_NAME = paplda  // "/RFtemp_" // auxstr // ".dat" 
		open (UNIT = 9090, FILE = FILE_NAME, STATUS = "OLD",
     ;            IOSTAT = OpenStat, ACTION = "READ")
	            if (OpenStat > 0) then
		          print *, 'Error reading file: ', FILE_NAME
	                stop
		      end if	  
c		read(9090,*),dummychar
            read(9090,*)rhocheck
	      read(9090,*)
	      read(9090,*), dummyi, dummyi
	      read(9090,*)

            do im = 1, 2*nmoant-1  ! loop over k// ------------------
               read(9090,*)kp(im) 

			sumaux2(1:klim+1) = 0.d0
			               
               do jv = 1,nabsv  ! v-loop - - - - - - - - 
                  read(9090,*)
                  read(9090,*) v(jv)
                  read(9090,*)
	            read(9090,*) Nxk

               	vel = v(jv)
				if( jv < lvg)then
					DELV = v(jv+1) - v(jv)
				else
					DELV = v(jv) - v(jv-1)
				end if


                  do ix = 1, Nxk

	               read(9090,"(2G16.8)") xxn(ix)	! (k//,v) dependence 'ommitted'
	               do km = 1, klim+1
	               read(9090,"(I10,G20.10)") dummyi, activ(ix,km)
                     end do

c			   activ(ix,km) = activ(ix,km) * NOVO * vel

cc	      Interpolate F0i to the given xxn(ix) point 
	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 2, isp), 1, xxn(ix), dfdv)
	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 3, isp), 1, xxn(ix), dfdx)
		       x = (1.d0 - dabs(xxn(ix)) ) * xmax          
                dfdx = dfdx *xmax  ! BATCH gives df/dxn, not df/dx

c		  Interpolate F0i to the previous xxn(ix-1) point
            if(ix .gt. 1)then
	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 2, isp), 1, xxn(ix-1), dfdv_prev)
	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 3, isp), 1, xxn(ix-1), dfdx_prev)
                 xprev = (1.d0 - dabs(xxn(ix-1)) ) * xmax
		       dfdx_prev = dfdx_prev *xmax
            end if

        if(ix .gt. 1)then
            sumaux2(1:klim+1) = sumaux2(1:klim+1) + 0.5d0
     ;      * ( activ(ix,  1:klim+1)  * (vel**3*dfdv -2.d0*vel**2*(x-xtir1)*dfdx) 
     ;        + activ(ix-1,1:klim+1)  * (vel**3*dfdv_prev-2.d0*vel**2*(xprev-xtir1)*dfdx_prev) ) 
     ;      * (xxn(ix)-xxn(ix-1))*xmax * DELV 
	  end if



	            end do  ! ix -loop

               end do ! jv: v-loop - - - - - - - - 
                 
	      m12test(im,1:klim+1) = sumaux2(1:klim+1) / fac ! ASPLAS factor 
		
		  write(1231,"(I10,20G18.6)"), im, kp(im), m12test(im,1:klim+1)

            end do  ! im : loop over k// ----------------------------
            
      
      close(9090)



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccc Power using M12 matrices computed above cccccccccccccccccccccccccccccc



		  powerdens2=0.d0
	      do m1 = minf(iel), msup(iel)
	      do m2 = minf(iel), msup(iel)

            maver = dfloat(m1+m2)/2.d0
            im = (m1+m2) - 2*minf(iel) +1  ! m_average index
	      km = iabs(m1-m2)+1                      ! mdiff index

            t = dreal(dconjg(xpmp(1,irad,m2+shiftm)) * xpmp(1,irad,m1+shiftm))
c            powerdens2 = powerdens2 + dimag(m12aux(im,km)) * t
		   powerdens2 = powerdens2 + m12test(im,km) * t

c 		  WARNING: Power density is in W/m, but we want W/m2 for int(p.r.dr)=Ptot
c		           --> ADD factor 1/rho
c				   --> ADD global factor fac (we saved M12/fac for compatib. with VMAT)
c				   --> ADD factor 1/ap ??????????????????????????	      

		   
            end do
	      end do

			powerdens2 = fac*powerdens2/rhoaux/rnorm

				print*, 'M12check:', rhoaux, powerdens2
				write(8888,*),  rhoaux, powerdens2


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccc Power using RFcoefs ccccccccccccccccccccccccccccccccccccccc

c	print*, '... Power using RF coefs (/Plot_data)'

c	Build F0i from F0(radqlfp)
      call interp_QLFP_noout(irad)


c	Integration over x and v (power density)	

	powerdens = 0.d0
	powerdens3 = 0.d0

      do jv = 1,nabsv
         
	   FACTOR = FAC1*FAC2*v(jv)

         do ix = 1, Nx

c	      Interpolate F0i to the given xnaux(ix) point 
	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 2, isp), 1, xnaux(ix), dfdv)
	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 3, isp), 1, xnaux(ix), dfdx)
		       x = (1.d0 - dabs(xnaux(ix)) ) * xmax          
                 dfdx = dfdx *xmax  ! BATCH gives df/dxn, not df/dx

c		  Interpolate F0i to the previous xnaux(ix-1) point
            if(ix .gt. 1)then
	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 2, isp), 1, xnaux(ix-1), dfdv_prev)
	      call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 3, isp), 1, xnaux(ix-1), dfdx_prev)
                 xprev = (1.d0 - dabs(xnaux(ix-1)) ) * xmax
		       dfdx_prev = dfdx_prev *xmax
            end if

c	      TRAPEZE integration in x, RECTANGLES in v

		  v_3 = v(jv)*v(jv)*v(jv)
	      v_2 = v(jv)*v(jv)
		  if( jv < lvg)then
			  DELV = v(jv+1) - v(jv)
		  else
			  DELV = v(jv) - v(jv-1)
		  end if

        if(ix .gt. 1)then

           powerdens3 = powerdens3 
     ;      + 0.5d0*(   RFcoef(jv,ix)  * (-1.d0/v(jv)*dfdv     +2.d0/v_2*(x-xtir1)*dfdx)   
     ;                + RFcoef(jv,ix-1)* (-1.d0/v(jv)*dfdv_prev+2.d0/v_2*(xprev-xtir1)*dfdx_prev)  ) 
     ;      * (xnaux(ix)-xnaux(ix-1))*xmax * DELV * FACTOR * v(jv)**4

	  end if


	   end do

	end do

c		write output to RFpow file	
				print*, 'RFcoef  :',rhocheck, -powerdens3/rhocheck/rnorm
				write(9999,*), rhocheck, -powerdens3/rhocheck/rnorm



      end do ! kr: Loop over BATCH surfaces
c     ======

		close(9999)
	      close(8001)
	close(8888)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc













cccccccccccccccc Power using M12 matrices computed above cccccccccccccccccccccccccccccc






c	Yet another way to compute the power (based on M12 matrices stored in M12run)

c	print*, '... Power using interpolated M12 matrices (OUT_RFCOEF)'


c      do kr = 1, nraddr
c     =================
      
c      irad = intaqlfp(kr)  ! intab index of QLFP radii
c      rhoaux = abscis(irad)

c      call int_to_string4(nint(1000*abscis(irad)),auxstr)
	
c	COKFOLDER = '../../M12run/r' // auxstr
	   
c	   NA = 2*(modva2-modva1)+1  ! Number of k// values 
c	   NMM = klim + 1		      ! Number of m1-m2 values

c         5.2) Auxiliary vector for file output: FORMAT identifier
c	    call INT_TO_STRING(2*NMM+1, charaux)
c	    call INT_TO_STRING(2*NMM+2, charaux2)
c	         GGs  = "(G11.5," // charaux // "G20.8)"
c	         GG2s = "(" // charaux2 // "G20.8)"

cc		p=+1:
c          FILE_NAME = trim(STDFOLDER)  // "/NonMaxM12" // 'SPEC_lefty.dat'
c			open (UNIT = 13, FILE = FILE_NAME, STATUS = "OLD",
cc     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
c     ;              IOSTAT = OpenStat, ACTION = "READ")
c	              if (OpenStat > 0) then
c		              print *, 'Error reading file: ', FILE_NAME
c	                  stop
c				  end if
c				  read (13, GGs), dummy, dummy, mdiffaux(1:2*NMM)
c				  do im = 1, NA
c					 read (13, GG2s),  dummyi, kp(im), m12aux(im,0:NMM-1)
c				  end do 
c			close (13)	
	
c		print*, 'rho=',abscis(irad)
c	do im = 1, NA

c	print*, kp(im), dimag(m12aux(im,0:4))

c	end do


c		  powerdens2=0.d0
c	      do m1 = minf(iel), msup(iel)
c	      do m2 = minf(iel), msup(iel)

c            maver = dfloat(m1+m2)/2.d0
c            im = (m1+m2) - 2*minf(iel) +1  ! m_average index
c	      km = iabs(m1-m2)+1                      ! mdiff index

c           t = dreal(dconjg(xpmp(1,irad,m2+shiftm)) * xpmp(1,irad,m1+shiftm))
cc            powerdens2 = powerdens2 + dimag(m12aux(im,km)) * t
c		   powerdens2 = powerdens2 + m12test(im,km) * t

c 		  WARNING: Power density is in W/m, but we want W/m2 for int(p.r.dr)=Ptot
c		           --> ADD factor 1/rho
c				   --> ADD global factor fac
c				   --> ADD factor 1/ap ??????????????????????????	      

c		   powerdens2 = fac*powerdens2/rhoaux/rnorm


c            end do
c	      end do

c				print*,  rhoaux, -powerdens2
c				write(8888,*),  rhoaux, -powerdens2

c      end do ! kr: Loop over BATCH surfaces
c     ======

c	close(8888)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      return
      end
