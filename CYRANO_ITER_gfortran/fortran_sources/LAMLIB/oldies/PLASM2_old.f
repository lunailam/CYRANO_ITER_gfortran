      SUBROUTINE PLASMB(PLA, PLAS, IPOINT, K, V, IV1, LDV, SWITCH, nfft, add)

      IMPLICIT NONE
      logical add
      INTEGER IPOINT, K, IV1, LDV, SWITCH, nfft
      COMPLEX*16 PLA(6,6), PLAS(6,6), V(IV1:IV1+LDV-1,19)

c     Builds volume plasma term matrix PMA
c     from matrix of poloidal harmonics V.
c     Harmonics of current order K are stored in Kth row of V.
c     Circular plasma cross section (magnetic field angle independent
c     of poloidal angle).

      include 'PARDIM.COPY'
      include 'COMSWE.COPY'
      include 'COMGEO.COPY'
      include 'COMSUB.COPY'
      include 'COMMA2.COPY'
      include 'COMPLA.COPY'
      include 'COMPHY.COPY'

      CHARACTER*3 ELETYP
      INTEGER KP, i, j
      DOUBLE PRECISION TEM
      COMPLEX*16 TEM2
     ;, PMA(6,6), PMAS(6,6)

      CALL ZSET(36, CZERO, PMA, 1)

      ELETYP = STYP(ISUBR)

c      KP = IABS(IPOINT)
      kp = ipoint
      if(ipoint.lt.0)kp = ipoint + nfft

      IF(SWITCH.EQ.1)THEN
C     ~~~~~~~~~~~~~~~~~~~
c     For volume terms

      IF(GLOMAX)THEN
C     --------------
      PMA(1,1) = - Y * V(KP,1)
      PMA(3,3) = - Y * V(KP,2)
      PMA(5,5) = - Y * V(KP,3)
        if(.not.circ .and. eletyp.eq.'M23')then
        PMA(1,3) = - Y * V(KP,4)
        PMA(1,5) = - Y * V(KP,5)
        PMA(3,5) = - Y * V(KP,6)
        PMA(3,1) = - PMA(1,3)
        PMA(5,1) = - PMA(1,5)
        PMA(5,3) = PMA(3,5)
        end if
c to see: this is for circular!
        IF( .NOT.COLDPL(IREG) .AND. FLROPS(IREG) )THEN
        CALL ZSET(36, CZERO, PMAS, 1)
        PMAS(1,1) = - Y * V(KP,17)
        PMAS(2,2) = - Y * V(KP,18)
        PMAS(3,3) = - Y * V(KP,19)
        PMAS(1,2) = - Y * V(KP,14)
        PMAS(1,3) = - Y * V(KP,7)
        PMAS(1,4) = - Y * V(KP,10)
        PMAS(1,5) = - Y * V(KP,8)
        PMA(1,6)  = - Y * V(KP,11)
        PMA(2,2)  = - Y * V(KP,4)
        PMA(2,4)  = - Y * V(KP,9)
        PMAS(3,4) = - Y * V(KP,15)
        PMAS(3,5) = - Y * V(KP,12)
        PMA(3,6)  = - Y * V(KP,13)
        PMA(4,4)  = - Y * V(KP,5)
        PMA(6,6)  = - Y * V(KP,6)

        PMAS(2,3) = - PMAS(1,4)
        PMAS(2,1) =   PMAS(1,2)
        PMAS(3,1) =   PMAS(1,3)
        PMAS(4,1) =   PMAS(1,4)
        PMAS(5,1) = - PMAS(1,5)
        PMA(6,1)  = - PMA(1,6)
        PMA(3,2)  =   PMA(2,3)
        PMA(4,2)  =   PMA(2,4)
        PMAS(4,3) =   PMAS(3,4)
        PMAS(5,3) = - PMAS(3,5)
        PMA(6,3)  = - PMA(3,6)

          IF(K.NE.0)THEN
          TEM = - K * CO
          PMAS(1,1) = PMAS(1,1) + TEM*(-V(KP,14))
          PMA(1,2)  = PMA(1,2) +  TEM*(-V(KP,4) )
          PMAS(1,3) = PMAS(1,3) + TEM*(-V(KP,10))
          PMA(1,4)  = PMA(1,4) +  TEM*  V(KP,9)
          PMAS(3,1) = PMAS(3,1) + TEM*(-V(KP,10))
          PMA(3,2)  = PMA(3,2) +  TEM*(-V(KP,9) )
          PMAS(3,3) = PMAS(3,3) + TEM*  V(KP,15)
          PMA(3,4)  = PMA(3,4) +  TEM*  V(KP,5)
          PMA(5,1)  = PMA(5,1) +  TEM*( V(KP,11))
          PMA(5,3)  = PMA(5,3) +  TEM*(-V(KP,13))
          PMAS(5,5) = PMAS(5,5) + TEM*  V(KP,16)
          END IF
C       (Real scalar times complex vector:)
        CALL zdscal(36, K0RN2, PMAS, 1)
        END IF

      ELSE
C     ----
c to see in non circ!
      PMA(1,1) = - Y * V(KP,1)
      PMA(3,3) = - Y * V(KP,2)
      PMA(5,5) = - Y * V(KP,3)
      PMA(1,3) = - Y * V(KP,4)
      PMA(1,5) = - Y * V(KP,5)
      PMA(3,5) = - Y * V(KP,6)
      PMA(3,1) = PMA(1,3)
      PMA(5,1) = PMA(1,5)
      PMA(5,3) = PMA(3,5)
      END IF
C     ------
      CALL zdscal(36, K0RN2, PMA, 1)

      ELSE IF(SWITCH.EQ.2)THEN
C     ~~~~~~~~~~~~~~~~~~~~~~~~
c     For boundary conditions
        IF(COLJUM)THEN

          IF(GLOMAX)THEN
          PMA(1,1) = V(KP,1)
          PMA(3,3) = V(KP,2)
          PMA(5,5) = V(KP,3)
            if(.not.circ .and. eletyp.eq.'M23')then
            PMA(1,3) = V(KP,4)
            PMA(1,5) = V(KP,5)
            PMA(3,5) = V(KP,6)
            PMA(3,1) = - PMA(1,3)
            PMA(5,1) = - PMA(1,5)
            PMA(5,3) = PMA(3,5)
            end if
          ELSE
          PMA(1,1) = V(KP,1)
          PMA(3,3) = V(KP,2)
          PMA(5,5) = V(KP,3)
          PMA(1,3) = V(KP,4)
          PMA(1,5) = V(KP,5)
          PMA(3,5) = V(KP,6)
          PMA(3,1) = PMA(1,3)
          PMA(5,1) = PMA(1,5)
          PMA(5,3) = PMA(3,5)
          END IF

        ELSE
c       Jump condition with flr contributions
        TEM2 = CI * OMEGAG * EPS0 * RNORM * 0.5D0
        PMA(1,1) = TEM2 * V(KP,1)
        PMA(1,2) = TEM2 * V(KP,2)
        PMA(1,3) = TEM2 * V(KP,3)
        PMA(1,4) = TEM2 * V(KP,4)
        PMA(3,1) = TEM2 * (-V(KP,3))
        PMA(3,2) = TEM2 * V(KP,4)
        PMA(3,3) = TEM2 * V(KP,5)
        PMA(3,4) = TEM2 * V(KP,6)
        PMA(5,1) = TEM2 * V(KP,7)
        PMA(5,3) = TEM2 * V(KP,8)
        PMA(5,6) = TEM2 * V(KP,9)
        END IF

      ELSE IF(SWITCH.EQ.3)THEN
C     ~~~~~~~~~~~~~~~~~~~~~~~~
C     Case SWITCH=1 *(-1/Y)

        IF(GLOMAX)THEN
        PMA(1,1) = V(KP,1)
        PMA(3,3) = V(KP,2)
        PMA(5,5) = V(KP,3)
          if(.not.circ .and. eletyp.eq.'M23')then
          PMA(1,3) = V(KP,4)
          PMA(1,5) = V(KP,5)
          PMA(3,5) = V(KP,6)
          PMA(3,1) = - PMA(1,3)
          PMA(5,1) = - PMA(1,5)
          PMA(5,3) = PMA(3,5)
          end if

          IF( .NOT.COLDPL(IREG) .AND. FLROPS(IREG) )THEN
          CALL ZSET( 36, CZERO, PMAS, 1 )
          PMAS(1,1) = V(KP,17)
          PMAS(2,2) = V(KP,18)
          PMAS(3,3) = V(KP,19)
          PMAS(1,2) = V(KP,14)
          PMAS(1,3) = V(KP,7)
          PMAS(1,4) = V(KP,10)
          PMAS(1,5) = V(KP,8)
          PMAS(1,5) = V(KP,8)
          PMA(1,6)  = V(KP,11)
          PMA(2,2)  = V(KP,4)
          PMA(2,4)  = V(KP,9)
          PMAS(3,4) = V(KP,15)
          PMAS(3,5) = V(KP,12)
          PMA(3,6)  = V(KP,13)
          PMA(4,4)  = V(KP,5)
          PMA(6,6)  = V(KP,6)
    
          PMAS(2,3) = - PMAS(1,4)
          PMAS(2,1) =   PMAS(1,2)
          PMAS(3,1) =   PMAS(1,3)
          PMAS(4,1) =   PMAS(1,4)
          PMAS(5,1) = - PMAS(1,5)
          PMA(6,1)  = - PMA(1,6)
          PMAS(3,2) =   PMAS(2,3)
          PMA(4,2)  =   PMA(2,4)
          PMAS(4,3) =   PMAS(3,4)
          PMAS(5,3) = - PMAS(3,5)
          PMA(6,3)  = - PMA(3,6)
    
            IF(K.NE.0)THEN
            YINV = ABSCNI(INTAB)
            TEM =   K * CO * YINV
            PMAS(1,1) = PMAS(1,1) + TEM*(-V(KP,14))
            PMAS(1,2) = PMAS(1,2) + TEM*(-V(KP,4) )
            PMAS(1,3) = PMAS(1,3) + TEM*(-V(KP,10))
            PMAS(1,4) = PMAS(1,4) + TEM*  V(KP,9)
            PMAS(3,1) = PMAS(3,1) + TEM*(-V(KP,10))
            PMAS(3,2) = PMAS(3,2) + TEM*(-V(KP,9) )
            PMAS(3,3) = PMAS(3,3) + TEM*  V(KP,15)
            PMAS(3,4) = PMAS(3,4) + TEM*  V(KP,5)
            PMAS(5,1) = PMAS(5,1) + TEM*( V(KP,11))
            PMAS(5,3) = PMAS(5,3) + TEM*(-V(KP,13))
            PMAS(5,5) = PMAS(5,5) + TEM*  V(KP,16)
            END IF
          CALL zdscal(36, K0RN2, PMAS, 1)
          END IF

        ELSE
        PMA(1,1) = V(KP,1)
        PMA(3,3) = V(KP,2)
        PMA(5,5) = V(KP,3)
        PMA(1,3) = V(KP,4)
        PMA(1,5) = V(KP,5)
        PMA(3,5) = V(KP,6)
        PMA(3,1) = PMA(1,3)
        PMA(5,1) = PMA(1,5)
        PMA(5,3) = PMA(3,5)
        END IF

      CALL zdscal(36, K0RN2, PMA, 1)

      END IF
C     ~~~~~~

      if(.not.add)then
      CALL ZCOPY(36, PMA, 1, PLA, 1)
      else 
        do i = 1, 6
          do j = 1, 6
          pla(i,j) = pla(i,j) + pma(i,j)
          end do
        end do
      end if

      IF(.NOT.COLDPL(IREG) .AND. FLROPS(IREG))THEN
        if(.not.add)then
        CALL ZCOPY(36, PMAS, 1, PLAS, 1)
        else 
          do i = 1, 6
            do j = 1, 6
            plas(i,j) = plas(i,j) + pmas(i,j)
            end do
          end do
        end if
      END IF

      RETURN
      END

