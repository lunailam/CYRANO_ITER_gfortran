
c	-------------------------------------------------

c     Function to compute (w-wc-k//v//) at given pol. angle theta

c	This version for standard coordinates k//=k//(theta)
c	k// is stored in COMMON kptab(1:npfft+1)

c     Other parameters are passed through COMMONS:
c	   VNOW XNOW sigNOW qomNOW kpNOW 


      FUNCTION f_ic_std(theta)

	implicit none

     
      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgeo.copy'
      include 'comgdr.copy' 


      real*8, intent(in) :: theta
	real*8 :: f_ic_std

	integer :: ipo
	real*8 :: vpar2, vperp2, vpar, wc, aux, v2, btt, sintt, kpp, Raux

	integer isrchfge
	external isrchfge

		 ! Find closest polang values
           ipo = isrchfge(npfft/2+1, polang, 1, theta)

	     if(polang(ipo).eq.theta)then  ! exact match on polang table

              btt = bmotab(intab,ipo)
			kpp = kptab(ipo)
c              sintt = eqt(intab,ipo,12)

           else                          ! interpolate on 4 closest points

			if(ipo>1)then

	        call spline0(4,polang(ipo-1:ipo+2),bmotab(intab,ipo-1:ipo+2),
     ;                     1, theta, btt)
	        call spline0(4,polang(ipo-1:ipo+2),kptab(ipo-1:ipo+2),
     ;                     1, theta, kpp)

			else

	        call spline0(3,polang(ipo:ipo+2),bmotab(intab,ipo:ipo+2),
     ;                     1, theta, btt)
	        call spline0(3,polang(ipo:ipo+2),kptab(ipo:ipo+2),
     ;                     1, theta, kpp)
			end if
			

           end if


ccccccccccccccccccccccccccccccccc

c	if(circ)then

c		Raux = r0+abscis(intab)*cos(theta) 
c		btt = r0*b0/ Raux /  eqt(intab,1,14)
c		kpp = kptab(1)
c	end if	


		  v2 = VNOW*VNOW
		  vpar2  = v2 * (1.d0 - XNOW*btt/B0)
c		  vperp2 = v2 - vpar2
		  vpar   = sigNOW * dsqrt(abs(vpar2))
		  wc     = qomNOW * btt
	
	      aux = omegag - wc - kpp*vpar

	      f_ic_std = aux
c		  f_ic_std = vpar2
c		f_ic_std = omegag			

      END FUNCTION f_ic_std
c	-------------------------------------------------


