      subroutine dshap2

      implicit none

c     D-shaped or any geometry: 'magnetic' terms at a grid point. 
c     Must be preceded by tabulation of 'geometrical' terms!
c     Grid point indices INTAB, INTABP input through COMSWE.

c     THETA is the angle between total and toroidal fields
c     RHOCUR: plasma current profile limit radius

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comswe.copy'
      include 'comphy.copy'

c      double precision f, g, totcur, gfac
c      external f, g, totcur, gfac

C     NB: F (tabulated) includes sign of B0 and sign of plasma current.
C         Jacobian JAC assumed > 0; NT is > 0.

        if(cyl)then
        btor = b0  * eqta1d(intab,2)
        else
        btor = b0 * eqta1d(intab,2) * r0orta(intab,intabp)
        end if

      tgthor = eqt(intab,intabp,7) * eqta1d(intab,1) 
     ;      / (eqt(intab,intabp,9) * eqta1d(intab,2))	

      tgth = tgthor * abscis(intab)

      theta = datan(tgth)
      if(b0.lt.0.d0)theta = theta -  dsign(pi, theta)
c     NB: THETA is angle between magnetic field and toroidal unit vector.

      si = dsin(theta)
      co = dcos(theta)

      bpol = btor * tgth
      bmodul = btor / co
c      bmodul = dsqrt(btor * btor + bpol * bpol)

c     Calculation of q and psi is done in routine tables (1d table)

      return
      end
