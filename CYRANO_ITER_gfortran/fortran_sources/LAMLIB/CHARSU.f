      subroutine charsu(charv)

      implicit none
      
      double precision charv(4)
         
C     Left hand cutoff  (charv(1) = L - n//^2), 
c	right hand cutoff (charv(2) = R - n//^2), 
c	FW perp.resonance (charv(3) = S - n//^2), 
c	SW perp. resonance(charv(4) = S).

c     REAL PART OF DISPERSION USED, OR COLD TENSOR.
c     Position passed through comswe: indices intab, intabp, ireg
c     LOCAL VALUE OF k// FOR M=0 MODE IS USED.
c     Species with a general equilibrium distribution are considered Maxwellian
c     here!

      include 'pardim.copy'
      include 'comphy.copy'
      include 'compla.copy'
      include 'comreg.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comgeo.copy'
      include 'commod.copy'
      include 'comswe.copy'

cERN      logical inscr
cERN      external inscr
	integer idum
      double precision kpar, keta, npar2, neta
     ;, lco, rco, fwpr, swpr
      complex*16 sigsd(3,3), sigod(3,3), splten(3,3,5), lstix, rstix, pstix
     ;, sstix, dstix

      lco = 1.
      rco = 1.
      fwpr = 1.
      swpr = 1.

      if(.not.vacuum(ireg))then
      kpar = kphi * r0orta(intab, intabp)
c      keta = - kpar * eqt(intab, intabp,15)
      kpar = kpar * eqt(intab, intabp,14)
      npar2 = (kpar / k0)**2
c      neta = (keta / k0)**2
        if(coldpl(ireg))then
        call colten(sigsd, .false., idum, .true., .false.)
        else
        call hotten(kpar, sigsd, sigod, .false., idum, .true., .false.)
        end if
      lstix = 1. + sigsd(1,1)
      rstix = 1. + sigsd(2,2)
c      pstix = 1. + sigsd(3,3)
      sstix = 0.5*(rstix+lstix)
c      dstix = 0.5*(rstix-lstix)

      lco = dreal(lstix) - npar2
      rco = dreal(rstix) - npar2
      fwpr = dreal(sstix) - npar2
      swpr = dreal(sstix)
      charv(1) = lco		! LH cutoff
      charv(2) = rco		! RH cutoff
      charv(3) = fwpr		! FW perp. resonance
      charv(4) = swpr		! SW perp. resonance

      end if
      return
      end
