      subroutine comag2

      implicit none

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comphy.copy'

      double precision f, g, dtgth, fpol, gdia
      external f, g

      rho = y * rnorm
      epsi = rho * r0i
      fpol = f(rho)
      gdia = g(rho)

      theta = datan(fpol/gdia)
      if(b0.lt.0.d0)theta = theta - dsign(pi,theta)
c     NB: THETA is angle between magnetic field and toroidal unit vector.

      si = sin(theta)
      co = cos(theta)

      if(rho .lt. 1.d-6)then
      dtgth=0.d0
      else
      dtgth = (f(rho*1.001d0)/g(rho*1.001d0)-f(rho*0.999d0)/g(rho*0.999d0)) * 5.d2 / y
      end if

      co1 =   co * co * dtgth
      si1 =   co1 * co
      co1 = - co1 * si

      return
      end
