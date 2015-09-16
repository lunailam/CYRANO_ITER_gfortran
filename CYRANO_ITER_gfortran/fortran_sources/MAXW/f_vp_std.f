
c	-------------------------------------------------

c     Function to compute v// at given pol. angle theta

c     Other parameters are passed through COMMONS:
c	   VNOW XNOW sigNOW qomNOW kpNOW 


      FUNCTION f_vp_std(theta)

	implicit none

      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'comgeo.copy'
     


      real*8, intent(in) :: theta
	real*8 :: f_vp_std

	integer :: ipo
	real*8 :: vpar2, vperp2, vpar, v2, btt

	integer isrchfge
	external isrchfge

		 ! Find closest polang values
           ipo = isrchfge(npfft/2+1, polang, 1, theta)

	     if(polang(ipo).eq.theta)then  ! exact match on polang table

              btt = bmotab(intab,ipo)

           else                          ! interpolate on 4 closest points

			if(ipo>1)then

	        call spline0(4,polang(ipo-1:ipo+2),bmotab(intab,ipo-1:ipo+2),
     ;                     1, theta, btt)

			else

c			btt = bmotab(intab,1)
	        call spline0(3,polang(ipo:ipo+2),bmotab(intab,ipo:ipo+2),
     ;                     1, theta, btt)
			end if
			

           end if


ccccccccccccccccccccccccccccccccc



		  v2 = VNOW*VNOW
		  vpar2  = v2 * (1.d0 - XNOW*btt/B0)
c		  vperp2 = v2 - vpar2
		  vpar   = sigNOW * dsqrt(abs(vpar2))
c		  wc     = qomNOW * btt
	
c	      aux = omegag - wc - kpp*vpar

	      f_vp_std = vpar
c		  f_ic_std = vpar2
c		f_ic_std = omegag			

      END FUNCTION f_vp_std
c	-------------------------------------------------


