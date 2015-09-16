      SUBROUTINE COLTEN(SIGSD, ONESPE, ISP, gdrtoo, onlyab)

      IMPLICIT NONE
      LOGICAL ONESPE, gdrtoo, onlyab
      INTEGER ISP
      COMPLEX*16 SIGSD(3,3)
      
C     COMPUTES +-// COMPONENTS OF SUSCEPTIBILITY TENSOR
C     IN THE COLD PLASMA LIMIT at current table point.
C     ARGUMENTS:
c     Input: onespe: switch (true: one species, else all species)
c            isp:    if onespe is true, index of the required species.
c            gdrtoo:
c            onlyab: if .true., only compute (collisional) absorption
c     Through common comswe: radial and poloidal table indices intab, intabp.
c     Output:  3*3 complex diag.matrix SIGSD


      include 'pardim.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'compla.copy'
      include 'comswe.copy'
      include 'compri.copy'
      include 'comfic.copy'
      include 'comphy.copy'

      LOGICAL  INSCR, INSC
      INTEGER I, J, ISP1, ISP2, ISPE, IWR
     ;, spli(maxspe), is, ncl
      DOUBLE PRECISION OMC, COEFF
      COMPLEX*16 OMCOL

      IWR = NOFILE
      call zset(3*3, czero, sigsd, 1)

      if(gdrtoo)then
        IF(ONESPE)THEN
        NCL = 1
        spli(1) = isp
        ELSE
        NCL = NSPEC
          do i = 1, ncl
          spli(i) = i
          end do
        END IF
      else
        IF(ONESPE)THEN
        if(spegdr(isp))return
        NCL = 1
        spli(1) = isp
        ELSE
        NCL = NSPEC - NSPGDR
        call icopy(ncl, ispnod, 1, spli, 1)
        END IF
      end if

      bmodul = BMOTAB(intab, intabp)
c     phi must be defined to find whether we are inside a screen:
      phi = polang(intabp)
      INSC = INSCR(IREG)

      DO 2 IS = 1, ncl
C     ~~~~~~~~~~~~~~~~
      ispe = spli(is)
        IF(INSC)THEN
      	OMP(ISPE) = OMPTAB2(INTAB,ISPE)
      	ELSE
      	OMP(ISPE) = OMPTAB(INTAB,ISPE)
      	END IF
      OMC = QOM(ISPE) * BMODUL
      OMCOL = OMEGAG + CI * DAMP(ISPE,IREG)

      COEFF = OMP(ISPE) ** 2 / OMEGAG
        if(onlyab)then
      SIGSD(1,1) = SIGSD(1,1) - dimag(COEFF / (OMCOL - OMC))
      SIGSD(2,2) = SIGSD(2,2) - dimag(COEFF / (OMCOL + OMC))
      SIGSD(3,3) = SIGSD(3,3) - dimag(COEFF /  OMCOL)
        else
      SIGSD(1,1) = SIGSD(1,1) - COEFF / (OMCOL - OMC)
      SIGSD(2,2) = SIGSD(2,2) - COEFF / (OMCOL + OMC)
      SIGSD(3,3) = SIGSD(3,3) - COEFF /  OMCOL
        end if
   2  CONTINUE
C     ~~~~~~~~

      RETURN
      END
