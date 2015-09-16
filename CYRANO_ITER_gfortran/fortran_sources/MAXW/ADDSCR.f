      SUBROUTINE ADDSCR

      IMPLICIT NONE

c     REVOIR EN 3D!

      include 'pardim.copy'
      include 'dynou2.copy'
      include 'comreg.copy'
      include 'comfin.copy'
      include 'comant.copy'
      include 'commod.copy'
      include 'comswe.copy'
      include 'comequ.copy'

      INTEGER I, IS, IM, IA, MP, MML, JS, JM, JA, JB, J, J2, LBL
     ;, IPVT(MAXSCR*MAXSCM), NRHSS, INFO

      DOUBLE PRECISION XTRA, Y
C     ;, RCOND 

      COMPLEX*16 ACST(MAXSCR*MAXSCM,MAXSCR*MAXSCM)
     ;          ,BCST(MAXSCR*MAXSCM,MAXANT)
     ;          ,TEMP
     ;          ,POLFAC

      include 'comphy.copy'

C     REDUCES FUNDAMENTAL SOLUTION USING CONSTRAINTS ON FARADAY SCREENS

      WRITE(6,*)'ENTER ADDSCR'
      if(maxscr.eq.0 .or. nscree.eq.0)return
      CALL ZSET( (MAXSCR*MAXSCM)**2, CZERO, ACST, 1 )
      CALL ZSET( (MAXSCR*MAXSCM)*MAXANT, CZERO, BCST, 1 )

      DO IS = 1, NSCREE
      IREG = IRBSCR(IS)
      Y = RX0(IREG)
      DO IM = 1, NMOSCR
      IA = (IS-1)*NMOSCR + IM

        DO MP = 1, NMOANT
        MML = MOSCR(IM) - MOANT(MP)
        POLFAC = LAMWID(IS) / RX0M(IREG) * CDEXP( -CI*MML*THEAS(IS) )
        XTRA = MML * LAMWID(IS) / RX0M(IREG) * 0.5D0
        IF( XTRA .NE. 0.D0 ) POLFAC = POLFAC * SIN(XTRA) / XTRA

          IF( NSCLAM(IS) .EQ. 1 )THEN
          XTRA = 0.D0
          ELSE
          XTRA = MML * DTHES(IS) / ( NSCLAM(IS) - 1 ) * 0.5D0
          END IF

          IF( dabs(DMOD(XTRA,PI)) .lt. 1.e-15 )THEN
c          IF( DMOD(XTRA,PI) .EQ. 0.D0 )THEN
          POLFAC = POLFAC * NSCLAM(IS)
          ELSE
          POLFAC = POLFAC
     ;    * DSIN( NSCLAM(IS)*XTRA ) / DSIN( XTRA )
          END IF

          DO JS = 1, NSCREE
            DO JM = 1, NMOSCR
            JA = (JS-1)*NMOSCR + JM
c           Shield modes between themselves:
            ACST(IA,JA) = ACST(IA,JA) + POLFAC * XRTPB(4,IREG,MP,NEFFAN+JA,1)
            END DO
          END DO

          DO JB = 1, NEFFAN
c         Shield modes with antenna modes: each matrix element receives from
c         a single antenna. Toroidal positions not relevant.
          BCST(IA,JB) = BCST(IA,JB) - POLFAC * XRTPB(4,IREG,MP,JB,1)
          END DO

        END DO

      END DO
      END DO

c     FACTORIZE ACST:
      call zgetrf(NSCREE*NMOSCR, NSCREE*NMOSCR, ACST, MAXSCR*MAXSCM, IPVT
     ;, info)
      if(info.ne.0)write(6,*)'addscr called zgetrf: info=',info
      
c     Compute A**-1 * B and store in B:
      call zgetrs('No transpose', NSCREE*NMOSCR, NEFFAN, ACST
     ;, MAXSCR*MAXSCM, IPVT, BCST, MAXSCR*MAXSCM, info)

c     Compute reduced fundamental solution and update XRTPB:
      DO 7 ITHOMA = 1, NTHOMA
      CALL RDSOL(IMOTO, ITHOMA, BEL, MAXBLL, LBL, NRHSS)
        DO J = 1, NEFFAN
          DO J2 = 1, NSCREE*NMOSCR
          TEMP = BCST(J2,J)
            DO I = 1, LBL
            BEL(I,J) = BEL(I,J) + TEMP * BEL(I,J2+NEFFAN)
            END DO
          END DO
        END DO 

      CALL WTSOL(IMOTO,ITHOMA,BEL,MAXBLL,LBL,NRHSS)
    7 CONTINUE

      DO 8 J = 1, NEFFAN
      DO 8 IREG = 1, NREG
      DO 8 I = 1, NDOF-1
      DO 8 J2 = 1, NSCREE*NMOSCR
      DO 8 MR = 1, NMODE(IMAX)
      XRTPB(I,IREG,MR,J,1) = XRTPB(I,IREG,MR,J,1)
     ; + BCST(J2,j) * XRTPB(I,IREG,MR,J2+NEFFAN,1)
    8 CONTINUE

      WRITE(6,*)'EXIT ADDSCR'
      RETURN
      END
