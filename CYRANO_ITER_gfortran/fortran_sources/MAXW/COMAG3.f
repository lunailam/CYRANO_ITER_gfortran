      subroutine comag3

      implicit none

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comphy.copy'

      double precision f, g, fpol
      external f, g

      rho = y * rnorm
      epsi = rho * r0i
      rmr0 = rho * dcos(phi)
      z = rho * dsin(phi)

      	if(cyl)then
      	ror0 = 1.
      	r0or = 1.
      	btor = b0
      	else
      	r = r0 + rmr0
      	ror0 = 1. + rmr0 * r0i
      	r0or = 1. / ror0
      	btor = b0 * r0or
      	end if

      fpol = f(rho)
      tgth = fpol / g(rho)
      bpol = btor * tgth
      bmodul = abs(btor) * sqrt(1.+fpol*fpol)

      return
      end
