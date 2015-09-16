      subroutine comag1(ind)

      implicit none
      integer ind

c     magnetic configuration - circular concentric surfaces

c     THETA: angle between total and toroidal fields (don't confuse with poloidal angle, noted phi in the code!)
c     ind: switch-  =1: cartesian--> polar
c                   =2: polar --> cartesian
c     
      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comphy.copy'

      double precision 
     ;  mu, geom, dtgth, fpol, gdia, totcur, gfac, f, g
     ;, cop, sip
c     ;, curr

      external totcur, gfac, f, g

        if(ind .eq. 1)then
c       RMR0 and Z are input

        rho = dsqrt(rmr0*rmr0+z*z)
          if(rho.eq.0.d0)then
          phi = 0.d0
          else
          mu = dasin(z/rho)
      	if(rmr0.lt.0.d0)then
            phi = mu
      	else
      	phi = pi - mu
      	end if
          end if
        cop = dcos(phi)
        sip = dsin(phi)

        else if(ind .eq. 2)then
c       RHO and PHI are given
        cop = dcos(phi)
        sip = dsin(phi)
        rmr0 = rho * cop
        z    = rho * sip
        end if
      
      r = rmr0 + r0
      drrho = cop
      drthn = - sip
      dzrho = sip
      dzthn = cop
      epsi = rho * r0i
      y = rho * rnori
      gdia = g(rho)

        if(cyl)then
        ror0 = 1.d0
        r0or = 1.d0
        btor = b0 * gdia
        else
        r = r0 + rmr0
        ror0 = 1.d0 + rmr0 * r0i
        r0or = 1.d0 / ror0
        btor = b0 * r0or * gdia
        end if

c      curr = totcur(rho)
      geom = gfac(rho)    ! This is parameter LAMBDA (thesis after eq. 1.25)
      fpol = f(rho)

      tgth = fpol / gdia
      theta = datan(tgth)
      if(b0.lt.0.d0)theta = theta - pi*dsign(1.d0,theta)
c     NB: THETA is angle between magnetic field and toroidal unit vector.

      si = sin(theta)
      co = cos(theta)

      bpol = btor * tgth

      if(j0.ne.0.d0)then
        if(rho .lt. 1.d-8)then
          if(cyl)then
          q = 0.
c         tgthor undefined
          else
          q = b0 * gdia / (pi*2.d-7*j0*r0)
          tgthor = r0i / q
          end if
        else
        q = geom / tgth
        tgthor = tgth / rho
        end if
      else
      tgthor = 0.d0
      end if

      bmodul = btor / co

      if(rho .lt. 1.d-8)then
      dtgth=0.d0
      else
      dtgth = (f(rho*1.001d0)/g(rho*1.001d0)-f(rho*0.999d0)/g(rho*0.999d0)) * 5.d2 / y
      end if

      co1 =   co**2 * dtgth
      si1 =   co1 * co
      co1 = - co1 * si

      return
      end
