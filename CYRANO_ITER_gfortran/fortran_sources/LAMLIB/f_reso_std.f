
c	-------------------------------------------------

c     Function to compute 1/(w-wc-k//v//)*0.5*(vperp/v)^2 dt
c     where dt = r / (v//.sin(PITCH)) * dtheta
c     for a given theta value (poloidal angle) 

c	This version for standard coordinates k//=k//(theta)
c	k// is stored in COMMON kptab(1:npfft+1)

c     Other parameters are passed through COMMONS:
c	   VNOW XNOW sigNOW qomNOW 
c	   <kpNOW> -> mavNOW 


      FUNCTION f_reso_std(theta)

	implicit none

          
      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgeo.copy'
      include 'comgdr.copy' 


      real*8, intent(in) :: theta
	real*8 :: f_reso_std

	integer :: ipo
	real*8 :: vpar2, vperp2, vpar, wc, aux, v2, btt, sintt, kpp

	integer isrchfge
	external isrchfge

		 ! Find closest polang values
           ipo = isrchfge(npfft/2+1, polang, 1, theta)

	     if(polang(ipo).eq.theta)then  ! exact match on polang table

              btt = bmotab(intab,ipo)
              sintt = eqt(intab,ipo,12)
			kpp = kptab(ipo)

cERN  March 07:   include ipo==1 possibility -----

	     elseif(ipo.eq.1)then

              btt = bmotab(intab,ipo)
              sintt = eqt(intab,ipo,12)
			kpp = kptab(ipo)

c     -------------------------------------------

           else                          ! interpolate on 4 closest points

	        call spline0(4,polang(ipo-1:ipo+2),bmotab(intab,ipo-1:ipo+2),
     ;                     1, theta, btt)
	        call spline0(4,polang(ipo-1:ipo+2),eqt(intab,ipo-1:ipo+2,12),
     ;                     1, theta, sintt)
	        call spline0(4,polang(ipo-1:ipo+2),kptab(ipo-1:ipo+2),
     ;                     1, theta, kpp)
           end if


		  v2 = VNOW*VNOW
		  vpar2  = v2 * (1.d0 - XNOW*btt/B0)
cccccccccccccccccccccccccccccccccccccccccc
c		vpar2  = dmax1(vpar2,1.d0)
ccccccccccccccccccccccccccccccccccccccccc
		  vperp2 = v2 - vpar2
		  vpar   = sigNOW * dsqrt(abs(vpar2))
		  wc     = qomNOW * btt
	
	      aux = 0.5d0 * vperp2 / v2 / abs(vpar) / sintt
     ;                 / (omegag - wc - kpp*vpar)

	      f_reso_std = aux

      END FUNCTION f_reso_std
c	-------------------------------------------------


