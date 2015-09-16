      subroutine dshap1(ind)

      implicit none
      integer ind

c     D-shaped geometry; compute the "purely geometrical" terms

c     THETA is the angle between total and toroidal fields
c     RHOCUR: plasma current profile limit radius
c     Conversion R,Z to RHO,PHI (ind=1)
c                RHO,PHI to R,Z (ind=2)

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comswe.copy'
      include 'comphy.copy'

      double precision 
     ;  sinp, cosp
     ;, sp2, x, sht, omx2
c     ;, f, totcur, gfac
         
c      external f, totcur, gfac

        cosp = dcos(phi)
        sinp = dsin(phi)
        sp2 = sinp * sinp
        x = rho * api
        omx2 = 1.d0 - x * x
c       Now Shafranov shift inside region 1 only:
          if(ireg.eq.1)then
            if(ishsht.eq.0)then
c           Shafranov shift: discontinuous deriv. at x=1!
            shsh = shsh0 * omx2
            sht = shsh0 * (-2.d0 * x * api)
            else if(ishsht.eq.1)then
c           To try: continuous deriv. at x=1; changed below as well: Jacobian...
            shsh = shsh0 * omx2 ** 2
            sht = shsh0 * 2.d0 * omx2 * (-2.d0 * x * api)
            else
            shsh = 0.d0
            sht = 0.d0
            end if
          else
          shsh = 0.d0
          sht = 0.d0
          end if

        rmr0 = rho * cosp + shsh - rho * sp2 * delta * x
        rmra = rmr0 - shsh0      
cERN       z = z0*(1-rho**2/rhowal**2) + kappa * rho * sinp
	z =  z0 + kappa * rho * sinp

      drrho = cosp + (sht - 2. * x * delta * sp2)
      drthn = - sinp * (1.d0 + delta * x * 2.d0 * cosp)
      drth = rho * drthn
cERN      dzrho = -z0*2*rho/rhowal**2+kappa * sinp
	dzrho = kappa * sinp
      dzthn = kappa * cosp
      dzth = rho * dzthn

      nrho2 = drrho * drrho + dzrho * dzrho
      nt2n = drthn * drthn + dzthn * dzthn
      nt2 = rho * rho * nt2n
      ntn = dsqrt(nt2n)
      nt = rho * ntn
c     g12 = drrho * drth + dzrho * dzth
c     jac = drrho * dzth - dzrho * drth
      jacn = drrho * dzthn - dzrho * drthn
c?    jacn = kappa * (1.d0 + sht * cosp)
      jac = rho * jacn
      jacav = kappa * rho
      g12n = drrho * drthn + dzrho * dzthn
      g12 = rho * g12n
cERN  cn = dzthn * delta * api * 2.d0 * sinp * cosp
c	cn = C / rho^2
	cn = kappa * delta * api * 2.d0 * sinp * cosp * cosp
c     newmu is mu of the notes, NT2 * mu = Ztt Rt - Rtt Zt
      newmu = kappa * (1.d0 + 2.d0*delta*x*cosp**3) / nt2n
c     lambda is as in the notes on geometry, NT2 * lambda = ZrtRt - RrtZt = C
c     (Not to be confused with lambda factor of thesis!)
      lambda = cn / nt2n

        if(cyl)then
        ror0 = 1.d0
        r0or = 1.d0
        r = rmr0
        else
        r = r0 + rmr0
c       these are now R/Ra and Ra/R:
        ror0 = 1.d0 + rmra * r0i  
        r0or = 1.d0 / ror0 ! r0or = Ra/R
        end if

	  
      return
      end
