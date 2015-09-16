      double precision function fno(x)

      implicit none
      double precision x

c     This is function f(x) (also in lamlib) divided by rho.
c     Specific to circular concentric equilibrium.

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comphy.copy'

      double precision gfacn, totcun
      external gfacn, totcun
      
      fno = 0.d0

      if(x.eq.0.d0)return

      ap = rhocur

c Problem with GFAC in general geom.:
c      if(geneq)then
c      else if(dshape)then
c      	if(x.lt.0.001*ap) then
c      	fno = j0*mu0/b0 * (kappa**2/(kappa**2+1.))
c      	else
c      	fno = 2.d-7*totcun(x)/(b0*ra*gfacn(x))
c      	end if
c      else
      
      	if(cyl)then
      		if(x .lt. 0.001*ap) then
      		fno = j0 * mu0 / b0 * 0.5d0
      		else
      		fno = 2.d-7 * totcun(x) / b0
      		end if
      	else
      		if(x .lt. 0.001*ap) then
      		fno = j0 * mu0 / b0 * 0.5d0 * (1.d0 - (x*r0i)**2 * 0.5d0)
      		else
c      		fno = 2.d-7 * totcun(x) / (b0 * ra * gfacn(x))
      		fno = 2.d-7 * totcun(x) / (b0 * r0 * gfacn(x))
      		end if
      	end if
      
c      end if

      return
      end
