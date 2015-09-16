      double precision function f(x)

      implicit none
      double precision x

c     Specific to circular concentric equilibrium.

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comphy.copy'

      double precision gfac, totcur
      external gfac, totcur
      
      f = 0.d0

      if(x.eq.0.d0)return

      ap = rhocur

c Problem with GFAC in general geom.:
c      if(geneq)then
c      else if(dshape)then
c      	if(x.lt.0.001*ap) then
c      	f = j0*mu0/b0 * (kappa**2/(kappa**2+1.)) * x
c      	else
c      	f = 2.d-7*totcur(x)/(b0*ra*gfac(x))
c      	end if
c      else
      
      	if(cyl)then
      		if(x .lt. 0.001*ap) then
      		f = j0 * mu0 / b0 * 0.5d0 * x
      		else
      		f = 2.d-7 * totcur(x) / (b0 * x)
      		end if
      	else
      		if(x .lt. 0.001*ap) then
      		f = j0 * mu0 / b0 * 0.5d0 * x *(1.d0 - (x*r0i)**2 * 0.5d0)
      		else
c      		f = 2.d-7 * totcur(x) / (b0 * ra * gfac(x))
      		f = 2.d-7 * totcur(x) / (b0 * r0 * gfac(x))
      		end if
      	end if
      
c      end if

      return
      end
