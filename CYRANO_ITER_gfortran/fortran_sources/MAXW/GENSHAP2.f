      subroutine genshap2

      implicit none

c     General shaped geometry: 'magnetic' terms at a grid point. 
c     Grid point indices INTAB, INTABP input through COMSWE.

c     THETA:  angle between total and toroidal fields - pitch angle
c	BTOR :  Toroidal field component



ccccccccccc  IMPORT Btor and PITCH from FLUSH
cERN	 2D Profiles of Btoroidal and PITCH angle read from
c	 files in 'read_equi.f'.
c	 Poloidal grid is fixed but radial needs interpolation

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comswe.copy'
      include 'comphy.copy'




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
