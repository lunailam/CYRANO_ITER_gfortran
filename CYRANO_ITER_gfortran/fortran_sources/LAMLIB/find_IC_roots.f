
c	-------------------------------------------------

c     Subroutine to compute the IC resonance roots for circular geometry

c	This version is based on the cubic equation obtained for 
c	the major radius R from the IC resonance equation
c	(w-wc-k//v//)^2=0 in circular geometry, where ktheta=m/r an B = BoRo/R/cosO

c     The cubic equation is a*R^3 + b*R^2 + c*R + d = 0
c	where
c	a = w^2/v^2 - kth^2*sinO^2
c     b = - 2*w*qom/v^2*alpha*Ro*Bo - 2*kth*N*sinO*cosO + x*kth^2*sinO^2*alpha*Ro
c	c = qom^2/v^2*alpha^2*Ro^2*Bo^2 - N^2*cosO^2 + 2*x*kth*N*sinO*Ro
c	d = x*N^2*cosO*Ro

c     with kth=m/rho, cosO and sinO (pitch angle), N = tor. number, 
c     alpha = 1/cosO (so that B = BoRo/R*alpha)



      SUBROUTINE find_IC_roots(qoms, vin, xin, mav, roots, Nr)

	implicit none

      include 'pardim.copy'
      include 'commag.copy'
      include 'comant.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgeo.copy'
      include 'comma2.copy'   


      real*8, intent(in)   :: qoms, vin, xin, mav
	real*8, intent(out)  :: roots(3)
	integer, intent(out) :: Nr

	integer :: j, number
	real*8 :: a,b,c,d, NN, kth, x
	real*8 :: alp, coef(4), cubroot(3)

	NN = motoan(1) 			
	kth = mav / abscis(intab)   ! ko   = m/r
      co = eqt(intab,1,14)  ! circular
	si = eqt(intab,1,15)  ! circular
      alp = 1/ co


c	Cubic equation

		x = xin

	   a = omegag**2/vin**2 - kth**2 * si**2

	   b = -2*omegag*qoms/vin**2*alp*r0*B0 
     ;       -2*kth*NN*si*co + x*kth**2*si*si*alp*r0
       
	   c =  qoms**2/vin**2*alp**2*r0**2*B0**2 - NN**2*co**2
     ;      + 2*x*kth*NN*si*r0
     
         d = x*NN**2*co*r0     

	  coef(1) = d
	  coef(2) = c
	  coef(3) = b
	  coef(4) = a

	  call PA03A(coef,cubroot,number)

c		Nr = number
            Nr = 0
	roots = cubroot

c	Roots are for R, we want for theta
	
	do j = 1, 3
	   if( dabs(cubroot(j)-r0)/abscis(intab) .le. 1.d0)then
	      Nr=Nr+1
            roots(j) = dacos( (cubroot(j)-r0)/abscis(intab) )
	   else
            roots(j) = 0.d0
	   end if

	end do

c	  if(number.eq.1)then
c	     Nr = 1
c	     roots(1) = cubroot(1)
c		 roots(2) = 0.d0
c       elseif(number>1)then
c	     Nr = 2
c	     roots(1) = cubroot(1)
c		 roots(2) = cubroot(2)
c	  else
c	     print*, 'Problems in cubic root finder PA03A!'
c	     print*, 'Stop at find_IC_roots for (xin,vin)', xin, vin
c	  end if

	  return

      END SUBROUTINE find_IC_roots
c	-------------------------------------------------


