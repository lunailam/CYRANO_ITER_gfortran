      double precision function gfac(x)

      implicit none
      double precision x

c     GFAC is the function called capital LAMBDA in my thesis, cf.near eq.1.5
C     X is the radial coordinate. This routine only applies to circular
C     equilibrium with concentric surfaces.
C     See routine TABLES for the general case.

      include 'pardim.copy'
      include 'comgeo.copy'

      double precision xx

      if(cyl)then
      gfac=0.d0

      else
      xx = x * r0i
      gfac = xx / sqrt(1.d0-xx*xx)

      end if

      return
      end
