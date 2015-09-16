      DOUBLE PRECISION FUNCTION FCO()

      implicit none
      
C     I.C. CUT-OFF CHARACTERISTIC FUNCTION at current table point.
C     REAL PART OF DISPERSION USED, OR COLD TENSOR.
c     Position passed through comswe: indices intab, intabp, ireg
C     LOCAL VALUE OF k// FOR M=0 MODE IS USED.
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

      LOGICAL INSCR
      EXTERNAL INSCR
      integer idum
      DOUBLE PRECISION KPAR,KETA,NPAR,NETA,X
      COMPLEX*16 SIGSD(3,3), SIGOD(3,3), SPLTEN(3,3,5), LSTIX, RSTIX, PSTIX
     ;, SSTIX, DSTIX

      FCO = 1.

      IF(.NOT.VACUUM(IREG))THEN
      KPAR = KPHI * r0orta(intab, intabp)
      KETA = - KPAR * eqt(intab, intabp,15)
      KPAR = KPAR * eqt(intab, intabp,14)
      NPAR = (KPAR / K0)**2
      NETA = (KETA / K0)**2
        if(coldpl(ireg))then
        CALL COLTEN(SIGSD, .FALSE., IDUM, .true., .false.)
        else
        CALL HOTTEN(KPAR, SIGSD, SIGOD, .FALSE., IDUM, .true., .false.)
        end if
      LSTIX = 1. + SIGSD(1,1)
      RSTIX = 1. + SIGSD(2,2)
c      PSTIX = 1. + SIGSD(3,3)
c      SSTIX = 0.5*(RSTIX+LSTIX)
c      DSTIX = 0.5*(RSTIX-LSTIX)

      fco = (npar - dreal(lstix)) * (npar - dreal(rstix))

C     FCO = ( PSTIX*(RSTIX*LSTIX-2*SSTIX*NPAR+NPAR**2-NETA*SSTIX
C    ; +NETA*NPAR) - NETA*(RSTIX*LSTIX-(NPAR+NETA)*SSTIX) )
C    ;     * (1.-OMC/(OMEGAG+CI*DAMP(2,IREG)))
c      FCO =  (RSTIX * LSTIX - 2. * SSTIX * NPAR + NPAR ** 2 - NETA * SSTIX
c     ;        + NETA * NPAR)
c     ;     * (1. - OMC / (OMEGAG + CI * DAMP(2,IREG)))

      END IF
      RETURN
      END
