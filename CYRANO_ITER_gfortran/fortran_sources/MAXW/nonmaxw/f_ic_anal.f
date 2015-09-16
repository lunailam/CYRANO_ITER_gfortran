
c	-------------------------------------------------

c     Function to compute (w-wc-k//v//) at given pol. angle theta

c     Other parameters are passed through COMMONS:
c	   VNOW XNOW sigNOW qomNOW mavNOW 


      FUNCTION f_ic_anal(theta)

	implicit none

      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'comgeo.copy'
      include 'comant.copy'     


      real*8, intent(in) :: theta
	real*8 :: f_ic_anal

	real*8 :: Raux, btt, kpp
	real*8 :: vpar2, vpar, wc,  v2, sintt, costt

c	Geometry

		sintt = eqt(intab,1,12)
	    costt = eqt(intab,1,14)
		rho = abscis(intab)

		Raux = r0 + rho*dcos(theta) 
		btt = r0*b0/Raux / costt
		kpp = mavNOW*sintt + motoan(1)/Raux*costt

		v2 = VNOW*VNOW
		vpar2  = v2 * (1.d0 - XNOW*btt/b0)
		vpar   = sigNOW * dsqrt(abs(vpar2))
		wc     = qomNOW * btt
	
	    f_ic_anal = omegag - wc - kpp*vpar


		

      END FUNCTION f_ic_anal
c	-------------------------------------------------


