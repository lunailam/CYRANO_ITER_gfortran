
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

c     -------------------------------------------------
c     NEW: Only return root(s) for given sigma = sig_in
c     -------------------------------------------------

      SUBROUTINE find_IC_roots_onesig(qoms, vin, xin, sig_in, mav, root_good, Nr)

	implicit none

      include 'pardim.copy'
      include 'commag.copy'
      include 'comant.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgeo.copy'
      include 'comma2.copy'   


      real*8, intent(in)   :: qoms, vin, xin, sig_in, mav
	real*8, intent(out)  :: root_good(2)
	integer, intent(out) :: Nr

	integer :: j, number, Naux
	real*8 :: a,b,c,d, NN, kth, x
	real*8 :: alp, coef(4), cubroot(3), roots(3)
      real*8 :: Baux, Raux, kp, vp

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
c	      roots = cubroot

c	Roots are for R, we want for theta
	Naux = 0
	do j = 1, 3
	   if( dabs(cubroot(j)-r0)/abscis(intab) .le. 1.d0)then
            roots(j) = dacos( (cubroot(j)-r0)/abscis(intab) )
	   else
            roots(j) = 0.d0
	   end if

	end do

c     Check which root corresponds to INPUT sig_in:

      Nr = 0
	root_good(1:2)=0.d0

      do j=2,3

         Raux = r0 + abscis(intab)*dcos(roots(j))
         Baux = r0*b0/Raux/co
	   kp = kth*si + NN/Raux*co 
	   vp = vin*dsqrt(dabs(1.d0-x*Baux/b0))

         if(dsign(1.d0,(omegag-qoms*Baux)/(kp*vp)) .eq. dsign(1.d0,sig_in))then
            root_good(Nr+1) = roots(j)
	      Nr=Nr+1
 	   endif

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

      END SUBROUTINE find_IC_roots_onesig
c	-------------------------------------------------


