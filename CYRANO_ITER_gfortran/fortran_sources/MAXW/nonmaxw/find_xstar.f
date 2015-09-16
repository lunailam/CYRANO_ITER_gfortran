
c	-------------------------------------------------

c     Function to compute the double root resonance parameter xstar,
c	which is the tangent resonance near the banana tip (xmax>xstar>xtir). 

c	This version is based on the cubic equation obtained for 
c	the major radius R from the IC resonance equation
c	w-wc-k//v//=0 in circular geometry, where ktheta=m/r an B = BoRo/R/cosO

c     The cubic equation is a*R^3 + b*R^2 + c*R + d = 0
c	where
c	a = w^2/v^2 - kth^2*sinO^2
c     b = - 2*w*qom/v^2*alpha*Ro*Bo - 2*kth*N*sinO*cosO + x*kth^2*sinO^2*alpha*Ro
c	c = qom^2/v^2*alpha^2*ro^2*Bo^2 - N^2*cosO^2 + 2*x*kth*N*sinO*cosO*alpha*Ro
c	d = x*N^2*cosO^2*alpha*Ro

c     with kth=m/rho, cosO and sinO (pitch angle), N = tor. number, 
c     alpha = 1/cosO (so that B = BoRo/R*alpha)



      FUNCTION find_xstar(qoms, vin, mav, xo, xm)

	implicit none

 
      include 'pardim.copy'
      include 'commag.copy'
      include 'comant.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgeo.copy'
      include 'comma2.copy' 

      real*8, intent(in) :: qoms, vin, mav, xo, xm
	real*8 :: find_xstar

	integer :: j
	real*8 :: a,b,c,d, NN, kth, x, pp,qq,rr, bigA, bigB
	real*8 :: ini, fim, teste, limit, alp

	NN = motoan(1) 			
	kth = mav / abscis(intab)   ! ko   = m/r
      co = eqt(intab,1,14)
	si = eqt(intab,1,15)
      alp = 1/ co

	a = omegag**2/vin**2 - kth**2 * si**2

c      baux = 

c	Search for discriminant D = bigA^3/27 + bigB^2/4 = 0

	limit = -1.d-8
	ini = xo  ! (xtir = w/wc0)
	fim = xm  ! (Xmax)

	do j=1,100
         x = 0.5d0*(ini+fim)

	   b = -2*omegag*qoms/vin**2*alp*r0*B0 
     ;       -2*kth*NN*si*co + x*kth**2*si*si*alp*r0
       
	   c =  qoms**2/vin**2*alp**2*r0**2*B0**2 - NN**2*co**2
     ;      + 2*x*kth*NN*si*co*alp*r0
     
         d = x*NN**2*co**2*alp*r0     

	   pp = b / a
	   qq = c / a
	   rr = d / a
	  
	   bigA	= (3.d0*qq-pp*pp) / 3.d0
	   bigB = (2.d0*pp*pp*pp - 9.d0*pp*qq + 27.d0*rr) / 27.d0 

c	   DISCRIMINANT:
	   teste = bigA*bigA*bigA/27.d0 + bigB*bigB/4.d0

	   if(teste.gt.limit)then
            fim = x
	   else
	      ini = x
	   end if
         if(dabs(teste-limit)/dabs(limit)<1.d-2)goto 1233

	end do

1233	   continue
	           
	   find_xstar = x



      END FUNCTION find_xstar
c	-------------------------------------------------


