      double precision function circ_shear(x)

      implicit none
      double precision x

c     Specific to circular concentric equilibrium.
c     Evaluates the shear d ln q / d rho at radial abscissa rho = x
c     (Should work in the cylindrical limit as well)
cPL30/1/2005

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comphy.copy'

      double precision roa, roa2, tra
      

      if(dabs(x).le.1.d-6)then
      circ_shear = 0.d0      

      else
      ap = rhocur
      circ_shear = 2.d0 / (x * (1 - (x*r0i)**2))    
        
	  if(x .lt. ap)then
	  roa = x / rhocur
	  roa2 = roa ** 2
	  tra = 1 - roa2
        circ_shear = circ_shear - 2.d0 * (alpha+1.d0) * roa * tra**alpha / (rhocur * (1.d0-tra**(alpha+1.d0)))
	  end if

        if(cyl)then
        else
        end if
      end if

      return
      end
