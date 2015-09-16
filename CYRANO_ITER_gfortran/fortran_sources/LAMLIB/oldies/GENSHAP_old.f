      subroutine genshap(k)

      implicit none


c     Arbitrary geometry imported from JET FLUSH; 

c     THETA is the pitch angle
c     RHOCUR: plasma current profile limit radius


      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comswe.copy'
      include 'comphy.copy'
      include 'coequi.copy'

      integer, intent(in) :: k 
      real*8 ::  sinp, cosp, sin2p, cos2p, x

        cosp = dcos(phi)
        sinp = dsin(phi)
        cos2p = dcos(2.0d0*phi)
	  sin2p = dsin(2.0d0*phi)
	  rho = dmax1(rho,1e-6)

c	  R coordinate 
	  r = R0_cyr(k) +  R1_cyr(k)*cosp + R2_cyr(k)*cos2p +
     ;			       R3_cyr(k)*sinp + R4_cyr(k)*sin2p		
c	  Z coordinate 
	  z = Z0_cyr(k) + (R1_cyr(k)*sinp - R2_cyr(k)*sin2p +
     ;			       R3_cyr(k)*cosp + R4_cyr(k)*cos2p) * KAP_cyr(k)
c     ;                  -Zaxis	

c	----------- R derivatives ----------------

c	  1) dR/dr 	  
	  drrho = dR0_cyr(k) + dR1_cyr(k)*cosp + dR2_cyr(k)*cos2p +
     ;			           dR3_cyr(k)*sinp + dR4_cyr(k)*sin2p	

c	  2) 1/r.dR/dtheta (R1_cyr(k) == rho)
	  drthn = 1.0d0        *(-sinp) + R2_cyr(k)/rho*(-2*sin2p) +
     ;	      R3_cyr(k)/rho * cosp  + R4_cyr(k)/rho*( 2*cos2p)

c	  3) dR/dtheta 	  
c	  drth = R1_cyr(k)*(-sinp) + R2_cyr(k)*(-2*sin2p) +
c     ;	     R3_cyr(k)*  cosp  + R4_cyr(k)*( 2*cos2p)	
        drth = rho * drthn

c	  4)d2R/dtheta2 	  
	  drth2 = R1_cyr(k)*(-cosp) + R2_cyr(k)*(-4*cos2p) +
     ;	      R3_cyr(k)*(-sinp) + R4_cyr(k)*(-4*sin2p)	

c	  5)d2R/drdtheta
	  drrhoth = dR1_cyr(k)*(-sinp) + dR2_cyr(k)*(-2*sinp) +
     ;		  dR3_cyr(k)*  cosp  + dR4_cyr(k)*( 2*cosp)

c	----------------- Z derivatives --------------------

c	  1) dZ/dr 
	  dzrho = dZ0_cyr(k) + 
     ;          (dR1_cyr(k)*sinp - dR2_cyr(k)*sin2p +
     ;		   dR3_cyr(k)*cosp + dR4_cyr(k)*cos2p) * KAP_cyr(k) +
     ;		  ( R1_cyr(k)*sinp -  R2_cyr(k)*sin2p +
     ;		    R3_cyr(k)*cosp +  R4_cyr(k)*cos2p) * dKAP_cyr(k)

c	  2) 1/r.dZ/dtheta (R1_cyr(k) == rho) 	  
	  dzthn = (1.0d0         *  cosp  - R2_cyr(k)/rho * ( 2*cos2p) +
     ;	       R3_cyr(k)/rho *(-sinp) + R4_cyr(k)/rho * (-2*sin2p) ) 
     ;           * KAP_cyr(k)

c	  3) dZ/dtheta 	
c	  dzth = (R1_cyr(k)*  cosp  - R2_cyr(k)* ( 2*cos2p) +
c     ;	      R3_cyr(k)*(-sinp) + R4_cyr(k)* (-2*sin2p)) * KAP_cyr(k)
        dzth = rho * dzthn

c	  4) d2Z/dtheta2 	
	  dzth2 = (R1_cyr(k)*(-sinp) - R2_cyr(k)*(-4*sin2p) +
     ;	       R3_cyr(k)*(-cosp) + R4_cyr(k)*(-4*cos2p)) * KAP_cyr(k)

c       5) d2Z/drdtheta
	  dzrhoth = (dR1_cyr(k)*  cosp  - dR2_cyr(k)*( 2*cos2p) +
     ;		     dR3_cyr(k)*(-sinp) + dR4_cyr(k)*(-2*sin2p)) * KAP_cyr(k) +
     ;		    ( R1_cyr(k)*  cosp  -  R2_cyr(k)*( 2*cos2p) +
     ;		      R3_cyr(k)*(-sinp) +  R4_cyr(k)*(-2*sin2p)) * dKAP_cyr(k)

c	  -- General quantities --

	  nrho2 = drrho*drrho + dzrho*dzrho
        nt2n = drthn*drthn + dzthn*dzthn
        nt2 = rho * rho * nt2n
        ntn = dsqrt(nt2n)
        nt = rho * ntn
        jacn = drrho*dzthn - dzrho*drthn
        jac = rho * jacn
cccccccccccccccccccccccccccccccc
c	  This variable is not used anywhere!!!
ccc        jacav = kappa * rho	
cccccccccccccccccccccccccccccccccccc
        g12n = drrho*drthn + dzrho*dzthn
        g12 = rho * g12n
c       newmu is mu of the notes, mu = (Ztt Rt - Rtt Zt) / NT2 
        newmu =  (dzth2*drth - drth2*dzth) / dmax1(nt2n,1e-10) 
c	  C = Zrt.Rt - Rrt.Zt
        cn = dzrhoth*drth - drrhoth*dzth
c       lambda is as in the notes, lambda = (ZrtRt - RrtZt)/NT2 = C/NT2
        lambda = cn / dmax1(nt2n,1e-10)

c	this is Raxis/R
      r0or = Ra/r
       
      return
      end subroutine genshap
