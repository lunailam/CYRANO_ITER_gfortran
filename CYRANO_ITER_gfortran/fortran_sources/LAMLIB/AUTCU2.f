      SUBROUTINE AUTCUT

      IMPLICIT NONE

C     This routine generates a minimal mesh according to the local
C     dispersion relation.
C     It partly relies on physical knowledge about the dispersion,
C     so use in different scenarii should always be checked!
C     14/09/89 without Bernstein mode.

      include 'pardim.copy'
      include 'comswe.copy'
      include 'comfic.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comfin.copy'
      include 'comin2.copy'
      include 'comnud.copy'
      include 'compla.copy' 
      include 'comber.copy'
      include 'commod.copy'
      include 'comant.copy'
      include 'comphy.copy'

      LOGICAL RISC(MAXREG), LOCUT, KIMOLO, INSC, INSCR

      EXTERNAL INSCR

      INTEGER MSUPL(0:MAXNEL), TMSUPL(0:MAXNEL,MAXREG)
     ;, MINFL(0:MAXNEL), TMINFL(0:MAXNEL,MAXREG)
     ;, NMODL(0:MAXNEL), TNMODL(0:MAXNEL,MAXREG)
     ;, IREGB, IR1, IR2, IR3, ISIGN, NROOTP
     ;, I1, I2, IM, IM2, JA, MMAXL, MMINL, MMAXLR, MMININ, NMMAIN
     ;, IAHA, ILO, IELELO, KT, IED, IELG
     ;, IDAMAX, I1ST, NSUM

      DOUBLE PRECISION SECOND, REMAIN

      DOUBLE PRECISION TCUT(0:MAXNEL,MAXREG), LEL, LELIN1(MAXPOM),
     ; LELIN2(MAXPOM), LELIN, LELPRE,
     ; TDAM(MAXPOM),TDAMMA, WKB, WKBI,
     ; KRLOR(MAXPOM),KRLOI(MAXPOM),KPRLO,KRLOMA,KPRE(MAXPOM),
     ; CUT(0:MAXNEL),YBOU,
     ; KILOMA,KIPRE(MAXPOM),KPILO

      COMPLEX*16 KRHOM, KRLOC(MAXPOM)
      
      EXTERNAL LONG

      WRITE(NOFILE,*)'ENTER AUTCUT ; time=',SECOND()-TIMIN

      NSREG = 0
      SX0(0) = RX0(0)
      SX0M(0) = RX0M(0)

      DO 1 IREG = 1, NREG
      NSREG = NSREG + 1
      RISC(IREG)=.FALSE.
      NS(IREG) = 1
      ISUBR = IREG
      SX0(ISUBR) = RX0(IREG)
      SX0M(ISUBR) = RX0M(IREG)

        IF( ISTYP(ISUBR) .EQ. 1 )THEN
        STYP(ISUBR) = 'HEC'
        IDDL(ISUBR) = 2 * NDOF
        ICONN(ISUBR) = NDOF
        ELSE IF( ISTYP(ISUBR) .EQ. 2 )THEN
        STYP(ISUBR) = 'M23'
        IDDL(ISUBR) = 2 * NDOF - 1
        ICONN(ISUBR) = NDOF - 1
        END IF

      IBUB(ISUBR) = IDDL(ISUBR) - 2 * ICONN(ISUBR)
   1  CONTINUE

      TDAMMA = - DLOG(DAMCUT)

CC    REGIONS ARE CUT IN ELTS.STARTING FROM ANTENNAE AND SCREENS:
C     DO 2 IA = 1, NTOANT + NSCREE
C     IF( IA.LE.NTOANT )THEN
C     IR1 = IRBANT(IA)
C     ELSE
C     IR1 = IRBSCR(IA-NTOANT)
C     END IF
C     IREGB = IR1
C     IR2 = IR1 + 1

CC    Left and right of antenna or screen:
C     Sweep region boundaries inwards:

      DO 2 IREGB = NREG-1, 1, -1
      IR1 = IREGB
      IR2 = IREGB + 1

C     Look at both sides of each boundary:
      DO 2 IREG = IR2, IR1, -1
      ISIGN = 2 * ( IREG - IREGB ) - 1
      LOCUT=.TRUE.
      KIMOLO = KILLMO
C     Find whether other region boundary is active :
C     In that case, use full mode set in region.

        DO JA = 1, NTOANT + NSCREE
          IF( JA.LE.NTOANT )THEN
          IR3 = IRBANT(JA)
          ELSE
          IR3 = IRBSCR(JA-NTOANT)
          END IF
        IF(IR3.EQ.IREGB+ISIGN)KIMOLO=.FALSE.
        END DO

      IF(.NOT.RISC(IREG) .AND. LOCUT)THEN
C     -----------------------------------
        IF( VACUUM(IREG))THEN
        NROOTP = 1
        ELSE IF(COLDPL(IREG) .OR. .NOT.FLROPS(IREG))THEN
        NROOTP=2
      	ELSE
      	NROOTP=3
      	END IF
      Y = RX0(IREGB)
      YBOU = RX0(IREGB+ISIGN)
      ISUBR = IREG

      IF( PASSIV(IREGB) )THEN
C     -----------------------
      MMAXL = TMSUPL(0,IR2)
      MMINL = TMINFL(0,IR2)
      MMAXLR = TNMODL(0,IR2)
      ELSE
C     ----
      	IF(MONOMO)THEN
      	MMAXL = MOANT(1)
      	MMINL = MMAXL
      	ELSE
C	     IF( (IPL.EQ.0. .OR. KPHI.EQ.0.) .AND. MODVA1*MODVA2.LT.0 )THEN
CC       Modes m and -m require same mesh when no plasma current or axisym:
C	     MMINL = 0
C	     MMAXL = MAX0( IABS(MODVA1), IABS(MODVA2) )
C	     ELSE
      	MMINL = MODVA1
      	MMAXL = MODVA2
C	     END IF
      	END IF
      MMAXLR=MMAXL+1-MMINL
      MMININ = MMINL
      NMMAIN = MMAXLR
      END IF
C     ------

      WRITE(NOFILE,*)
     ;'Cutting region ',IREG,' from boundary ',IREGB,
     ;' with starting mode set ',MMINL,MMAXL
      CUT(0) = Y
      MSUPL(0) = MMAXL
      MINFL(0) = MMINL
      NMODL(0) = MMAXLR
      IF(.NOT.PASSIV(IREGB))CALL DSET( NMMAIN, 0., TDAM, 1 )
      IELELO = 0
C     WRITE(6,*)IELELO,Y,MMAXL
      REMAIN = RL(IREG)
      KPRLO = 0.

  10  CONTINUE
C     ========

      MSTUD1 = MMINL
      MSTUD2 = MMAXL
      MSTUST = 1
      CALL MAXROO(MMINL,MMAXL,MONOMO,KRLOC,NROOTP)

      DO 16 IM = 1, MMAXLR
      IM2 = IM + MMINL - MMININ
      KRLOC(IM) = KRLOC(IM) * RNORM
      KRLOR(IM) = DREAL(KRLOC(IM))
      KRLOI(IM) = DIMAG(KRLOC(IM))
      KRLOMA = KRLOR(IM)
      KILOMA = KRLOI(IM)
      	IF(IELELO.GT.0)THEN
      	KPRLO = ABS(KRLOMA - KPRE(IM2)) / LELPRE
      	KPILO = ABS(KILOMA - KIPRE(IM2)) / LELPRE
      	ELSE
      	KPRLO = 10. * ELPWL * KRLOMA**2 / PI
      	KPILO = 20. * ELPDL * KILOMA**2
      	END IF

      	IF(KRLOMA.NE.0.D0)THEN
      	WKB = KPRLO / ( KRLOMA**2 * ELPWL ) * 2.*PI
      	ELSE IF(KPRLO.NE.0.D0)THEN
      	WKB = 1.D40
      	ELSE
      	WKB = 0.
      	END IF

      	IF(WKB.LE. 0.1)THEN
      	LELIN1(IM) = ELPWL * KRLOMA * 0.5 / PI
      	ELSE IF(WKB .GE. 10.)THEN
      	LELIN1(IM) = SQRT(KPRLO*ELPWL*0.25/PI)
      	ELSE
      	LELIN1(IM) = KPRLO / (( SQRT(1.D0+2.*WKB) - 1.D0 )*KRLOMA)
      	END IF

      	IF(KILOMA.NE.0.D0)THEN
      	WKBI = KPILO / ( KILOMA**2 * ELPDL )
      	ELSE IF(KPILO.NE.0.D0)THEN
      	WKBI = 1.D40
      	ELSE
      	WKBI = 0.
      	END IF

      	IF(WKBI.LE. 0.1)THEN
      	LELIN2(IM) = ELPDL * KILOMA
      	ELSE IF(WKBI .GE. 10)THEN
      	LELIN2(IM) = SQRT(KPILO*ELPDL*0.5)
      	ELSE
      	LELIN2(IM) = KPILO / (( SQRT(1.D0+2.*WKBI) - 1.D0 )*KILOMA)
      	END IF

  16  CONTINUE

      I1 = idamax(MMAXLR,LELIN1,1)
      I2 = idamax(MMAXLR,LELIN2,1)
      LELIN = DMAX1( LELIN1(I1), LELIN2(I2) )
      LEL = DMIN1( REMAIN, 1./LELIN )
      IF(LEL.LT.LELMIN)
     ;WRITE(NOFILE,*)'Caution: element density is lower than'
     ;//' dispersion based value in local elt.',IELELO+1
C     Puts a limit on element smallness: (input LELMIN)
      LEL = DMAX1(LEL,LELMIN)
C     Trick for sharp interface: to generalize !
C     IF(IELELO.LE.20 .AND. .NOT.VACUUM(IREG) .AND. PLSTEP(IREGB))THEN
C     LEL = AMIN1(LEL,LELEDG)
C     END IF
      Y = CUT(IELELO) + ISIGN * LEL
      REMAIN = ISIGN * ( YBOU - Y )
      IF(REMAIN.LT.0.D0)THEN
      LEL = ISIGN * ( YBOU - CUT(IELELO) )
      REMAIN = 0.D0
      Y = YBOU
      END IF
      IELELO = IELELO + 1
      WRITE(6,*)IELELO,': Re on ',I1-1+MMINL,',Im on ',I2-1+MMINL,' LEL=
     ; ',LEL
      CUT(IELELO) = Y

      DO 17 M = MMINL, MMAXL
      IM = M + 1 - MMINL
      IM2 = M + 1 - MMININ
      TDAM(IM2) = TDAM(IM2) + KRLOI(IM) * LEL
      KPRE(IM2) = KRLOR(IM)
      KIPRE(IM2) = KRLOI(IM)
  17  CONTINUE

C     Drop the extreme modes when they are damped to DAMCUT*
C     their value at antenna:
      IF(KIMOLO)THEN

      	IF( TDAM(MMAXL+1-MMININ).GE.TDAMMA )THEN
      	WRITE(NOFILE,*)'Mode ',MMAXL,' dropped from local element ',
     ;	IELELO,' at rho=',Y*RNORM
      	MMAXL = MMAXL - 1
      	MMAXLR = MMAXLR - 1
      	END IF

      	IF( .NOT.MONOMO
     ;	.AND. TDAM(MMINL+1-MMININ).GE.TDAMMA )THEN
      	WRITE(NOFILE,*)'Mode ',MMINL,' dropped from local element ',
     ;	IELELO,' at rho=',Y*RNORM
      	MMINL = MMINL + 1
      	MMAXLR = MMAXLR - 1
      	END IF

      END IF

C     WRITE(6,*)IELELO,CUT(IELELO),MMAXL
      MSUPL(IELELO) = MMAXL
      MINFL(IELELO) = MMINL
      NMODL(IELELO) = MMAXLR

      IF( REMAIN .GT. 0. )THEN
      LELPRE = LEL
      IF( MMAXLR .GT. 0 )GOTO 10
      IELELO = IELELO + 1
      MSUPL(IELELO) = 0
      MINFL(IELELO) = 0
      NMODL(IELELO) = 0
      END IF

      CUT(IELELO) = RX0(IREGB+ISIGN)
      IELE(ISUBR) = IELELO
C     WRITE(6,*)ISUBR,IELE(ISUBR)
      IAHA = IELE(ISUBR) * ((1-ISIGN)/2)

      DO 11 IEL = 0, IELE(ISUBR)
      ILO = IAHA + ISIGN * IEL
      TCUT(IEL,IREG) = CUT( ILO )
      TMINFL(IEL,IREG) = MINFL( ILO )
      TMSUPL(IEL,IREG) = MSUPL( ILO )
      TNMODL(IEL,IREG) = NMODL( ILO )
  11  CONTINUE

      IF( IELELO.GT.1 )THEN
      RISC(IREG) = .TRUE.
C     ...the job is done.
      ELSE
C     Do not allow less than 2 elts. per region:
      IELE(ISUBR) = 2
      KT = IELE(ISUBR)
      CALL CUTEQ( SX0(ISUBR-1), SX0(ISUBR), KT, TCUT(0,IREG) )
      DO 14 IEL = 0, IELE(ISUBR)
      TMSUPL(IEL,IREG) = MSUPL( IAHA )
      TMINFL(IEL,IREG) = MINFL( IAHA )
      TNMODL(IEL,IREG) = NMODL( IAHA )
  14  CONTINUE
      RISC(IREG) = .TRUE.
      END IF

      END IF
C     ------

   2  CONTINUE

C     CUT REMAINING REGIONS:
      DO 3 IREG = 1, NREG
      	IF( .NOT. RISC(IREG) )THEN
      	ISUBR = IREG
      	IED = I1ST(ISUBR)
      	KT = IELE(ISUBR)
      	CALL CUTEQ( SX0(ISUBR-1), SX0(ISUBR), KT, TCUT(0,IREG) )

      	DO 7 IEL = 1, KT
      		IF(MONOMO)THEN
      		TMINFL(IEL,IREG)= MOANT(1)
      		TMSUPL(IEL,IREG)= MOANT(1)
      		TNMODL(IEL,IREG)= 1
      		ELSE
      		TMINFL(IEL,IREG)= MODVA1
      		TMSUPL(IEL,IREG)= MODVA2
      		TNMODL(IEL,IREG)= MODVA2-MODVA1+1
      		END IF
   7  CONTINUE

      RISC(IREG) = .TRUE.
      END IF
   3  CONTINUE

      NELE = NSUM( NSREG, IELE )
      WRITE(6,*)'NELE=',NELE
      IF( NELE .GT. MAXNEL )STOP 'MAXNEL too small'

C     STORE THE CUT IN GLOBAL VECTOR FX0:
      FX0(0) = RX0(0)
      DO 4 IREG = 1, NREG
      ISUBR = IREG
      IED = I1ST(ISUBR)
      KT = IELE(ISUBR)
      	DO 5 IEL = 1, KT
      	IELG = IEL + IED
      	FX0(IELG) = TCUT(IEL,IREG)
      	MSUP(IELG) = TMSUPL(IEL,IREG)
      	MINF(IELG) = TMINFL(IEL,IREG)
      	NMODE(IELG) = TNMODL(IEL,IREG)
   5  	CONTINUE
      CALL LONGEURS( KT, TCUT(0,IREG), FL(IED+1), 1 )
      WRITE(NOFILE,2001) ISUBR, KT, STYP(ISUBR)
      WRITE(NOFILE,*)(FX0(IELG),IELG=IED,KT+IED)
   4  CONTINUE

      MSUP(0)=MSUP(1)
      MINF(0)=MINF(1)
      NMODE(0)=NMODE(1)
      CALL LONGEURS( NSREG, SX0, SL, 1 )

      WRITE(NOFILE,*)'EXIT AUTCUT ; time=',SECOND()-TIMIN
      RETURN

 2001 FORMAT(1H ,'SUBREGION ',I1,': ',I3,' ELEMENTS OF TYPE ',A3)

      END
C
C********************************************************
C
      SUBROUTINE MAXROO(MMINL,MMAXL,MONOMO,KRHOM,NROOTP)

      IMPLICIT NONE

C     VALUE OF RADIAL WAVE VECTOR RELEVANT FOR MESH SIZE

      include 'pardim.copy'

      LOGICAL MONOMO
      INTEGER MMINL, MMAXL, NROOTP
      COMPLEX*16 KRHOM(MAXPOM)

      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comswe.copy'
      include 'commod.copy'
      include 'compla.copy'
      include 'comphy.copy'

      LOGICAL INSCR
      INTEGER IPHI, NPHI, IRO1, IRO
      DOUBLE PRECISION A, REM, IMM, KP1, KP2
      COMPLEX*16 AC, CRI1, CRI2, CRI3, CRI4, KRHO2(3,MAXPOM)
      
      EXTERNAL INSCR

C     Watch dispersion on chords PHI=0, PHI=PI :
      NPHI = 1
      IF(.NOT.POLSYM)NPHI = 2
      DO 3 IPHI = 1, NPHI
      PHI = (IPHI-1)*PI
      CALL COMAG2
      CALL COMAG3
      CALL YPRO(INSCR(IREG))
      CALL DISPER(1, KRHO2, NROOTP, .false.)
C     Choose k**2 root with largest real part
C     slow wave doesn't influence mesh in standard situations.
      IF(IPHI.NE.1)THEN
      IRO1 = 1
      ELSE
      IRO1 = 2
      CALL ZCOPY(MMAXL-MMINL+1, KRHO2, 3, KRHOM, 1)
      END IF
      DO 1 IRO = IRO1, NROOTP
      DO 1 M = MMINL, MMAXL
      MR = M + 1 - MMINL
      IF( DREAL(KRHOM(MR)) .LT. DREAL(KRHO2(IRO,MR)) )
     ; KRHOM(MR) = KRHO2(IRO,MR)
   1  CONTINUE
   3  CONTINUE
      DO 2 M = MMINL, MMAXL
      MR = M + 1 - MMINL
      KRHOM(MR) = SQRT(KRHOM(MR))
      KRHOM(MR) = DCMPLX( DABS(DREAL(KRHOM(MR))), DABS(DIMAG(KRHOM(MR))) )
   2  CONTINUE
C     IF(.NOT.CYL)THEN
C     KP1 = (N/(R0 + RHO))**2
C     KP2 = (N/(R0 - RHO))**2
C     END IF
C
C     MMAXLR=MMAXL+1
C     DO 1 M = 0, MMAXL
C     MR = M + 1
C     IF( Y .EQ. 0. .AND. M .NE. 0 )THEN
C     KRHOM(MR) = (0.D0,1.D40)
C    ;       + SQRT( CUN*(K02 - KPHI**2 ) )
C     ELSE
C
C     A = K02
C     IF( RHO.NE.0. ) A = A - (M/RHO)**2
C     AC = A
C
C	     IF( CYL )THEN
C	     CRI3 = SQRT( AC - KPHI**2 )
C	     REM = DABS(DREAL(CRI3))
C	     IMM = DABS(DIMAG(CRI3))
C	
C	     ELSE
C	     CRI1 = SQRT( AC - KP1 )
C	     CRI2 = SQRT( AC - KP2 )
C	     REM = AMAX1( DABS(DREAL(CRI1)),
C	    ;             DABS(DREAL(CRI2)) )
C	     IMM = AMAX1( DABS(DIMAG(CRI1)),
C	    ;             DABS(DIMAG(CRI2)) )
C		     IF( .NOT.MONOMO )THEN
C		     CRI3 = SQRT( CUN*( K02 - (N/(R0 + RHO))**2 ) )
C		     CRI4 = SQRT( CUN*( K02 - (N/(R0 - RHO))**2 ) )
C		     REM = DAMAX1( REM, DABS(DREAL(CRI3)), DABS(DREAL(CRI4)) )
C		     IMM = AMAX1( IMM, DABS(DIMAG(CRI3)), DABS(DIMAG(CRI4)) )
C		     END IF
C	
C	     END IF
C     KRHOM(MR) = DCMPLX( REM, IMM )
C     END IF
C  1  CONTINUE

      RETURN
      END

C     FUTUR: MAILLAGE "LARMOR"
C AJOUTER SWITCH FLROPS
C     Y=FX0(IEL)
C     PHI=PI
C     CALL YPRO(INSCR(IREG))
C     CALL COMAG2
C     CALL COMAG3
C     IF( .NOT. ( MONOMO .OR. VACUUM(IREG) .OR. COLDPL(IREG) ) )THEN
C     RLAN=MH/EEL*VT(ISPRES,1)/BMODUL/RNORM
C          IF(RLAN.EQ.0.D0)THEN
C          MRL=MSTUD2
C          ELSE
C          MRL=IDINT(2.D0*PI*Y/RLAN)
C          END IF
C     MSUP(IEL) = MAX0( MIN0( MRL, MSTUD2 ), 1 )
C     MINF(IEL) = - MSUP(IEL)
C     NMODE(IEL) = MSUP(IEL) - MINF(IEL) + 1
C
C     CHECKS VALIDITY OF THEORY : LOCAL DISPERSION OF BERNSTEIN MODE
C     DO 10 J=MINF(IEL),MSUP(IEL)
C     KPAR=DFLOAT(J)*SI/RHO-DFLOAT(N)*CO/R
C     CALL BERDIS(KPAR)
C 10  CONTINUE
C     IF( .NOT. VACUUM(IREG) )THEN
C     TRA = BMODUL
C     IF(.NOT.CYL) TRA = TRA * ROR0
C     T = TRA * FLOAT(NHA) * ZCH(ISPRES) / AMASS(ISPRES) * EEL/MH
C     TR0 = ( T / OMEGAG - 1.D0 ) / Y
C
C     IF( DABS(TR0).LE.1.D0 .AND. .NOT.COLDPL(IREG)
C    ;   .AND. IPL.NE.0.D0 .AND. VT(ISPRES,1).NE.0.D0 )THEN
C            ...a resonant layer of finite thickness is present.
C     THE0=DACOS(TR0)
C     PHI=THE0
C
C     CALL COMAG2
C     CALL COMAG3
C     CASE OF PLASMA CURRENT =  ZERO :   SEE LATER RESONANCE LOOK
C     KPHILO=KPHI
C     IF(.NOT.CYL)KPHILO=KPHI*R0/R
C     CORES(IEL)= KPHILO*Y*RNORM*CO/SI
C     T33=T*PI*RNORM*Y**2*DSIN(PHI)/XBOU/SI/VT(ISPRES,1)
C     SL1(IEL)=CORES(IEL)-T33
C     SL2(IEL)=CORES(IEL)+T33
C     SL1(IEL)=MAX0( NINT(SL1(IEL)), MINF(IEL)-1 )
C     SL2(IEL)=MIN0( NINT(SL2(IEL)), MSUP(IEL)+1 )
C     ELSE
C     CORES(IEL)=0.
C     SL1(IEL)=0
C     SL2(IEL)=0
C     END IF
C     WRITE(NOFILE,100) IEL,MRL,MINF(IEL),MSUP(IEL)
C    ;                 ,INT(CORES(IEL)),INT(SL1(IEL)),INT(SL2(IEL))
C     END IF
