      DOUBLE PRECISION FUNCTION FLH()

      implicit none
      
c  (Artisanal...)
C     Perpendicular resonance characteristic function at current table point.
C     REAL PART OF DISPERSION USED, OR COLD TENSOR
C     Position passed through comswe: indices intab, intabp, ireg
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
      DOUBLE PRECISION KPAR, X
      COMPLEX*16
     ;  SIGSD(3,3), LSTIX, RSTIX
     ;, SIGOD(3,3), SPLTEN(3,3,5)

      FLH = 1.d0
      IF(.NOT.VACUUM(IREG))THEN
      KPAR = KPHI * eqt(intab, intabp,14) * r0orta(intab, intabp)
        if(coldpl(ireg))then
        CALL COLTEN(SIGSD, .FALSE., IDUM, .true., .false.)
        else
        CALL HOTTEN(KPAR, SIGSD, SIGOD, .FALSE., IDUM, .true., .false.)
        end if
      LSTIX = 1.d0 + SIGSD(1,1)
      RSTIX = 1.d0 + SIGSD(2,2)

C BAD FUNCTION (POLE) FOR ROOT SEARCH IN COLD PLASMA:
      FLH = 0.5 * DREAL(LSTIX + RSTIX) - KPAR ** 2 / K02

C     TESTING...:
c      FLH = (0.5*DREAL(LSTIX+RSTIX) - KPAR**2/K02)
c     ;          * (1.d0 - OMC/(OMEGAG+CI*DAMP(2,IREG)))
      END IF
      RETURN
      END
