
c	-------------------------------------------------

c     Subroutine to expand the denominator B(theta) = w-wc-k//v//
c	at given pol. angle theta

c	Returns the vector yout = [a0 a1 ...] of the Taylor expansion
c	of the function B = w-wc-k//v// = b0 + b1(th-theta) + ...
c	around input angle 'theta'
c	Only first order for the moment !!!!!

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Other parameters are passed through COMMONS:
c	   VNOW XNOW sigNOW qomNOW mavNOW 
ccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE expand_B(theta, yout)

	implicit none

      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'comgeo.copy'
      include 'comant.copy'     


      real*8, intent(in) :: theta
      real*8, intent(out) :: yout(1:2) !!! First order

	real*8 :: Raux, btt, kpp
	real*8 :: vpar2, vpar, wc,  v2, sintt, costt, vperp2

c	Geometry

		sintt = eqt(intab,1,12)
	    costt = eqt(intab,1,14)
		rho = abscis(intab)

		Raux = r0 + rho*dcos(theta) 
		btt = r0*b0/Raux / costt
		kpp = mavNOW*sintt + motoan(1)/Raux*costt
		wc     = qomNOW * btt

		v2 = VNOW*VNOW
		vpar2  = v2 * (1.d0 - XNOW*btt/b0)
		vpar   = sigNOW * dsqrt(dabs(vpar2))
		vperp2 = v2 - vpar2


c         Value at theta
          yout(1) = dabs(vpar) * (omegag - wc - kpp*vpar);

c         First derivative
          yout(2) = dabs(vpar) * ( - qomNOW - vpar*motoan(1)*costt**2/(r0*b0) + 
     ;                  kpp*v2*XNOW/(2.d0*vpar*b0) )
     ;                * b0*r0/costt*rho*dsin(theta)/Raux**2
     ;             - (omegag-wc-kpp*vpar) 
     ;                * v2/(2.d0*vpar)*XNOW*r0/costt *rho*dsin(theta)/Raux**2


      END SUBROUTINE expand_B
c	-------------------------------------------------


