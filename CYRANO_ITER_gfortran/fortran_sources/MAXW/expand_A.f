
c	-------------------------------------------------

c     Subroutine to expand the numerator 
c	A(theta) = vperp2/v^2 * cos(mdiff*theta)
c	at given pol. angle theta

c	Returns the vector yout = [a0 a1 ...] of the Taylor expansion
c	of the function A = vperp2/v^2*cos() = a0 + a1(th-theta) + ...
c	around input angle 'theta'
c	Only first order for the moment !!!!!

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Other parameters are passed through COMMONS:
c	   VNOW XNOW sigNOW qomNOW mavNOW 
ccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE expand_A(theta, md, yout)

	implicit none

      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'comgeo.copy'
      include 'comant.copy'     

c	integer, intent(in) :: md
      real*8, intent(in)  :: theta, md
      real*8, intent(out) :: yout(1:2) !!! First order

	real*8 :: Raux, btt, kpp
	real*8 :: vpar2, vpar, wc,  v2, sintt, costt, vperp2

c	Geometry
c		sintt = eqt(intab,1,12)
	    costt = eqt(intab,1,14)
		rho = abscis(intab)

c	Field, etc...
		Raux = r0 + rho*dcos(theta) 
		btt = r0*b0/Raux / costt
c		kpp = mavNOW*sintt + motoan(1)/Raux*costt

c	Orbits
		v2 = VNOW*VNOW
		vpar2  = v2 * (1.d0 - XNOW*btt/b0)
		vpar   = sigNOW * dsqrt(dabs(vpar2))
		vperp2 = v2 - vpar2

c	Output
c         Value at theta
          yout(1) = vperp2/v2*dcos(md*theta);

c         First derivative
          yout(2) = XNOW*r0*rho/costt * dcos(md*theta) * dsin(theta)/Raux**2
     ;            - md*vperp2/v2*dsin(md*theta);


      END SUBROUTINE expand_A
c	-------------------------------------------------


