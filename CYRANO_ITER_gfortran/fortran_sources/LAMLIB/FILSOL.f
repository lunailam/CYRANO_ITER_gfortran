      SUBROUTINE FILSOL

      IMPLICIT NONE

C     FILLS ALL DOF IN SOLUTION VECTOR USING BOUNDARY CONDITIONS
C     AT MAGNETIC AXIS.

      include 'pardim.copy'
      include 'comusr.f'
      include 'comswe.copy'
      include 'comequ.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'commod.copy'
      include 'comant.copy' 
      include 'comfic.copy'
      include 'comfin.copy'
      include 'comfou.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comphy.copy'

      CHARACTER*3 ELETYP
      INTEGER LBL, NRHSS, ICOLO, IBULO, M1, M2, INDEXT, JA, ID, I2, I3, IBU
      DOUBLE PRECISION AXL

      CALL RDSOL(IMOTO, 1, BEL, LDBEL, LBL, NRHSS)

      ELETYP=STYP(1)
      ICOLO = ICONN(1)
      IBULO = IBUB(1)
      IEL=1
C     -----
      M1=MINF(IEL)
      M2=MSUP(IEL)
      INDEXT = NMODE(IEL-1)*ICOLO + MAX0(NMODE(IEL-1), NMODE(IEL))*IBULO

      IF(CROWN)THEN
      IF(ELETYP.NE.'HEC')RETURN

      DO 2 JA = 1, neffan + NSCREE*NMOSCR
C     ===================================

                DO 41 MR = 1, NMODE(IEL)
C               -------------
      M = MR - 1 + MINF(IEL)
c      ID = ISHIBC + (MR-1)*ICOLO
      ID = (MR-1) * ICOLO
      BEL(3+ID, JA) = BEL(1+ID, JA)
  41  CONTINUE

    2 CONTINUE
C     ========
      ELSE
      DO 1 JA = 1, neffan + NSCREE*NMOSCR
C     ===================================

                DO 42 MR = 1, NMODE(IEL)
C               -------------
      M = MR - 1 + MINF(IEL)
c      ID = ISHIBC + (MR-1)*ICOLO
      ID = (MR-1)*ICOLO
                   IF(ELETYP.EQ.'HEC')THEN
      IF( M.EQ.0 .OR. (M/2)*2.NE.M ) GOTO 42
      IF(M.LT.-1)THEN
      AXL=-(M+2.D0)/(M-2.D0)
      I2=2
      I3=4
      BEL( I3+ID, JA ) = AXL * BEL( I2+ID, JA )
      END IF
      IF(M.GT.1)THEN
      AXL=-(M-2.D0)/(M+2.D0)
      I2=4
      I3=2
      BEL( I3+ID, JA ) = AXL * BEL( I2+ID, JA )
      END IF
                   END IF
                   IF( ELETYP.EQ.'M23' .AND. IABS(M).EQ.1 )THEN
      BEL( 2+ID, JA ) = CI * M * BEL( 1+ID, JA )
                   END IF
                   IF( ELETYP.EQ.'M23' .AND. (M/2)*2.NE.M )THEN
C BUBBLE WAS ELIMINATED FOR ODD M
      IBU = ISHIBC + NMODE(IEL-1)*ICOLO + (MR-1)*IBULO + 1
      BEL(IBU,JA) = 0.75D0*BEL(1+ID,JA)
     ;            + 0.25D0*BEL(1+ID+INDEXT,JA)
                   END IF
  42  CONTINUE

    1 CONTINUE
C     ========
      END IF

C     VOIR:
      CALL WTSOL(IMOTO,1,BEL,LDBEL,LBL,NRHSS)
      RETURN
      END
