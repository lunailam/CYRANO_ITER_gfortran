C
C++++++TO BE INCLUDED IN CYRANO++++++++++++++++++++++++++++++++++++++
C
       SUBROUTINE NEWTEN(CASE, 
     +                   OMC, OMP, OM, 
     +                   KPERP, RKPAR, 
     +                   VELEC, VPERP, VPAR, 
     +                   V0, UDRIFT, 
     +                   KNEW)

      IMPLICIT NONE
      CHARACTER*5 CASE
      DOUBLE PRECISION OMC, OMP, OM, RKPAR, VELEC, VPERP, VPAR, V0, UDRIFT
      COMPLEX*16 KPERP, KNEW(3,3)
      
C
C      THIS ROUTINE IS AN INTERFACE BETWEEN CYRANO AND
C      RAYMOND'S GENERAL DIELECTRIC TENSOR ROUTINES
C
       COMPLEX*16     KPERPD, KS(3,3)
       COMPLEX*16     L(3,3), LINVER(3,3), PROD(3,3)
C
C      COMMONS PASSING INFO TO RAYMOND'S ROUTINES
C
       DOUBLE PRECISION OMCD, OMPD, RKPARD, OMD, VC, VPARR, VPERPR
       COMMON/CSPECI/OMCD, OMPD
       COMMON/WAVEVE/KPERPD, RKPARD, OMD
C
C      CASE =  'CASE1'  : MAXWELLIAN DISTRIBUTION,   ALPHAS IN D
C              'CASE2'  : SLOWING DOWN DISTRIBUTION, ALPHAS IN D
C              'CASE3'  : MAXWELLIAN DISTRIBUTION,   D BEAM IN D
C              'CASE4'  : SLOWING DOWN DISTRIBUTION, D BEAM IN D
C      OMC/OMP : CYCLOTRON AND PLASMA FREQ. OF SPECIES
C      OM      : ANTENNA FREQ.
C      KPERP   : PERPENDICULAR WAVE VECTOR
C      RKPAR   : PARALLEL WAVE VECTOR
C      VELEC   : THERMAL VELOCITY OF THE ELECTRONS
C      VPERP   : PERPENDICULAR THERMAL VELOCITY
C      VPAR    : PARALLEL THERMAL VELOCITY
C      V0      : SLOWING DOWN BIRTH VELOCITY
C      VC      : SLOWING DOWN CRITICAL VELOCITY (COMPUTED AT EACH STEP)
C

       DATA  L/     (0.707106781,0.),(+0.707106781,0.),(0.,0.),
     +              (0.,0.707106781),(0.,-0.707106781),(0.,0.),
     +              (0.000000000,0.),(0.,+0.000000000),(1.,0.)/,
     +       LINVER/(+0.707106781,0.),(0.,-0.707106781),(0.,0.),
     +              (+0.707106781,0.),(0.,+0.707106781),(0.,0.),
     +              (+0.000000000,0.),(0.,+0.000000000),(1.,0.)/
C
C      TO PASS VARIABLES, THE FOLLOWING IDENTIFICATIONS ARE NEEDED:
C      (PLACING A COMMON IN CYRANO WOULD AVOID TIME LOSS HERE)
C
       OMD   =OM
       OMCD  =OMC
       OMPD  =OMP
       KPERPD=KPERP
       RKPARD=RKPAR
C
C      RAYMOND USES KT/M AND NOT 2KT/M FOR THERMAL VELOCITIES**2 WHILE
C      IN CYRANO & BRAYCOH, THE DEFINITION 2KT/M IS USED.
C
       VPERPR=VPERP/1.414213562
       VPARR =VPAR /1.414213562
C
C      CALL OF RAYMOND'S ROUTINES
C
       IF(CASE.EQ.'CASE1'.OR.CASE.EQ.'CASE3')   THEN
               CALL NMKSMA(VPERPR,VPARR,UDRIFT,KS)
         ELSE
               IF(CASE.EQ.'CASE2') THEN
                        VC=0.09*VELEC
                  ELSE
                        VC=0.07*VELEC
               ENDIF
               CALL NMKSSD(V0,VC,UDRIFT,KS)
       ENDIF
C      WRITE(6,*)' UNROTATED TENSOR:'
C      WRITE(6,20)KS
 20    FORMAT(3(1P,6E12.2/))
C
C      WE COMPUTE THE ROTATED TENSOR
C
       CALL MCRCR(3,3,KS,3,3,3,LINVER,3,3,3,PROD,3)
       CALL MCRCR(3,3,L,3,3,3,PROD,3,3,3,KNEW,3)
       RETURN
       END
C===================================================================
C
C
      SUBROUTINE NMKSSD(V0,VC,UDRIFT,KS)
C
      IMPLICIT NONE
C
      INTEGER    I,J
      REAL*8     AUX,AUX1,UDRIFT,V0,VC
      REAL*8     FJM1,FJ,FJP1,FNORM,UJ,UJP1
      COMPLEX*16 KS(3,3)
C
C
      INTEGER    NFUN,NFUNMX
      PARAMETER (NFUNMX = 200)
      REAL*8     ARI(NFUNMX),FCI(NFUNMX),FUI(NFUNMX),FVI(NFUNMX)
      COMMON/CDIFUN/ARI,FCI,FUI,FVI,NFUN
C
C**********************************************************************
C*Definition of the representation parameter  NFUN                *****
C*                                                                *****
      NFUN=20
C* This definition corresponds to a relative error generally better ***
C* than 5e-3, except for small Vc/V0 ratios for which the error   *****
C* can become larger near the origin (v below vc<<v0)             *****
C* Use NMKST for more checking.                                   *****
C*                                                                *****
C* The constructed function is normalized to 1.                   *****
C**********************************************************************
C
      AUX=V0/NFUN
      AUX1=AUX/VC
      FNORM=1.D0/DLOG(1.D0+(V0/VC)**3)
      FJM1=1.D0
      FJ=1.D0
      DO 400 I=1,NFUN
      ARI(I)=1.D0
      FUI(I)=UDRIFT
      FVI(I) = I*AUX
      J=NFUN-I+1
      FJP1 = FJ
      FJ   = FJM1
      FJM1 = FNORM * DLOG(1.D0+( (J-1)*AUX1 )**3)
      UJP1 = (DFLOAT(J+1)/DFLOAT(J))**3-1.D0
      UJ   = 1.D0 - (DFLOAT(J-1)/DFLOAT(J))**3
      FCI(J) =  (FJ-FJM1)/UJ - (FJP1-FJ)/UJP1
C     IF(J.EQ.1.OR.J.EQ.20)WRITE(6,*)J,FCI(J)
400   CONTINUE
C
      CALL NMKS(KS)
C
      RETURN
      END
C NMKSMA
C
C Computation of the species dielectric tensor contribution Ks for a
C Maxwellian distribution using the NMKS representation
C The number of elements used is NFUN and the maximum range is given by
C XIMAX. These values are FIXED here in a DATA statement: see below.
C For checking accuracy of the representation, etc., use NMKST.
C
C                                                        2
C Input variables: VTHP Perpendicular thermal velocity (V =KT/m)
C                  VTHL Parallel      thermal velocity    "
C                  UDRIFT Parallal drift velocity
C
C Output variables: KS=Ks
C
      SUBROUTINE NMKSMA(VTHP,VTHL,UDRIFT,KS)
C
      IMPLICIT NONE
C
      INTEGER    I,J
      REAL*8    AR,AUX,NMKSTF,UDRIFT,VTHP,VTHL,XIMAX
      REAL*8    FJM1,FJ,FJP1,FJNORM,UJ,UJP1
      COMPLEX*16 KS(3,3)
C
C
      INTEGER    NFUN,NFUNMX
      PARAMETER (NFUNMX = 200)
      REAL*8    ARI(NFUNMX),FCI(NFUNMX),FUI(NFUNMX),FVI(NFUNMX)
      COMMON/CDIFUN/ARI,FCI,FUI,FVI,NFUN
C
C**********************************************************************
C*Definition of the representation parameters NFUN                *****
C*                                        and XIMAX= max(v/vthp)  *****
      DATA XIMAX/4./
      NFUN=20
C* This definition corresponds to a max. relative error on function ***
C* values of 0.018 for the largest v/vth (near 4.)                *****
C* and is <1.D-2 for nearly all the range.                        *****
C*                               3                                *****
C* For this XIMAX=4., Integral(dv Maxw.) = 0.9988660...           *****
C* but the reconstructed maxw. is exactly normalized to 1.        *****
C**********************************************************************
C
      AR=VTHP/VTHL
C
      AUX=XIMAX/NFUN
      FJM1 = NMKSTF(XIMAX)
      FJNORM=FJM1
      FJM1=1.D0
      FJ=1.D0
      DO 400 I=1,NFUN
      ARI(I)=AR
      FUI(I)=UDRIFT
      FVI(I) = I*AUX*VTHP
      J=NFUN-I+1
      FJP1 = FJ
      FJ   = FJM1
      FJM1 = NMKSTF((J-1)*AUX)/FJNORM
      UJP1 = (DFLOAT(J+1)/DFLOAT(J))**3-1.D0
      UJ   = 1.D0 - (DFLOAT(J-1)/DFLOAT(J))**3
      FCI(J) =  (FJ-FJM1)/UJ - (FJP1-FJ)/UJP1
C     IF(J.EQ.1.OR.J.EQ.20)WRITE(6,*)J,FCI(J)
400   CONTINUE
C
      CALL NMKS(KS)
C
      RETURN
      END
C
C
C NMKSTF                        2
C               _   _         -z /2
C F(z) = erf(z/V2)-V(2/pi) z e
C
C
      DOUBLE PRECISION FUNCTION NMKSTF    (Z)
C
      IMPLICIT NONE
C
      EXTERNAL DERF
      REAL*8  derf,S2,S2PI,Z
C
      DATA S2/ 1.4142 13562 37309 505D0/,
     .     S2PI/ 7.9788 45608 02865 356D-1/
C
      NMKSTF = derf(Z/S2) - S2PI*Z*DEXP(-0.5D0*Z**2)
      RETURN
      END
C NMKS
C
C Computation of the contribution Ks of species s to the
C                    dielectric tensor K = I + Sum Ks
C                                               s
C
C    Output is Ks=KSUM ; all elements are defined
C
C Summation is made over the distribution function representation
C            f = Sum Ci fi        i=1,...,NFUN
C                 i
C                                         2    2   2       2
C            fi = (3ri/4pi) vi**-3 gamma( vi - v - ri (v//-Ui)  )
C                                              p
C  where ri is the anisotropy factor r=vperp/vpar
C        Ui is the parallel drift velocity
C
C  ri is stored in ARI
C  Ui is stored in FUI
C  Ci is stored in FCI
C  vi is stored in FVI      max dimension is NFUNMX
C  NFUN, Ci and Vi are input through the common /CDIFUN/
C
C
C
      SUBROUTINE NMKS(KSUM)
C
      IMPLICIT NONE
C
      INTEGER    I,IERR,ISTOP,J,N,NMAX
      REAL*8    A,AC,AR,AUX,BN,ERR,MOPOM2,OMEGA,USV,V0
      COMPLEX*16 DKS(3),KPERP,KPRL,KPRL2,KS(3,3),KSUM(3,3),MN(3,3)
C
      REAL*8    OMC,OMP,RKLL,RELER2
      COMMON/CSPECI/OMC,OMP
      COMMON/WAVEVE/KPERP,RKLL,OMEGA
      PARAMETER (RELER2=1.D-10)
C RELER2 is the square of the relative error used for truncation of the
C KS series over n.
C
      INTEGER    IFUN,NFUN,NFUNMX
      PARAMETER (NFUNMX=200)
      REAL*8    ARI,FCI,FUI,FVI
      COMMON/CDIFUN/ARI(NFUNMX),FCI(NFUNMX),FUI(NFUNMX),FVI(NFUNMX),
     .              NFUN
C
      DATA NMAX/100/
C          NMAX   is the maximum abs(n) in the summation
C     DATA PI3/9.4247 77960 76937 97D0/
C
C
C
      IF(NFUN.GT.NFUNMX) GOTO 900
      MOPOM2=-(OMP/OMEGA)**2
C
      DO 180 I=1,3
      DO 180 J=1,3
180   KSUM(I,J)=0.
C
C loop 400 over distribution functions *********************************
      DO 400 IFUN=1,NFUN
      V0=FVI(IFUN)
      AR=ARI(IFUN)
      USV=AR*FUI(IFUN)/V0
C
C
      AC=OMC/(RKLL*V0)
      A =(OMEGA-RKLL*FUI(IFUN)) /(RKLL*V0)
      KPRL=KPERP*V0/OMC
      KPRL2=KPRL**2
      ISTOP=0
C
C n=0
C
      CALL NMMN(0,A,KPRL2,KPRL,AC,MN,AR,USV)
        DO 200 I=1,3
          DO 200 J=I,3
200       KS(I,J)=MN(I,J)
C
        DO 300 N=1,NMAX
C N>0
        BN=A-N*AC
        CALL NMMN(N,BN,KPRL2,KPRL,AC,MN,AR,USV)
          DO 220 I=1,3
          DKS(I)=MN(I,I)
            DO 220 J=I,3
220         KS(I,J)=KS(I,J)+MN(I,J)
C N<0
        BN=A+N*AC
        CALL NMMN(-N,BN,KPRL2,KPRL,AC,MN,AR,USV)
          DO 240 I=1,3
          DKS(I)=MN(I,I)+DKS(I)
            DO 240 J=I,3
240         KS(I,J)=KS(I,J)+MN(I,J)
C
C       Convergence tests
        ERR=0.D0
        IERR=0
          DO 260 I=1,3
          AUX=KS(I,I)*CONJG(KS(I,I))
          IF(AUX.NE.0.D0) THEN
          ERR = ERR + DKS(I)*CONJG(DKS(I))/AUX
          IERR=IERR+1
          ENDIF
260       CONTINUE
C         WRITE(6,*)'ERR=',ERR
C
          IF(ERR.LT.RELER2) THEN
          IF(ISTOP.GT.2) GOTO 320
          ISTOP=ISTOP+1
          ELSE
          ISTOP=0
          ENDIF
C
300     CONTINUE
      WRITE(6,*) '  NMKS --- series has not converged at order n=',
     .       NMAX
      WRITE(6,*) '  ****** --- relative error= ', ERR
320   CONTINUE
C     write(6,*) ' n=',N
C
        DO 340 I=1,3
        DO 340 J=I,3
340     KS(I,J)=-1.5D0*KS(I,J)
        DO 350 I=1,3
350     KS(I,I)=(1.D0,0.D0)+KS(I,I)
C
        DO 370 I=1,3
        DO 370 J=I,3
370     KSUM(I,J)=KSUM(I,J)+FCI(IFUN)*MOPOM2*KS(I,J)
C
400   CONTINUE
C loop 400 end *********************************************************
C
      KSUM(2,1)=-KSUM(1,2)
      KSUM(3,1)= KSUM(1,3)
      KSUM(3,2)=-KSUM(2,3)
C     WRITE(6,*)'MOPOM2=',MOPOM2
C     WRITE(6,*)'KSUM'
C     WRITE(6,*)KSUM(1,1),KSUM(1,2),KSUM(1,3)
C     WRITE(6,*)KSUM(2,1),KSUM(2,2),KSUM(2,3)
C     WRITE(6,*)KSUM(3,1),KSUM(3,2),KSUM(3,3)
      RETURN
C
900   WRITE(6,*) 'NMKS --- Number of component distribution functions ',
     .           'NFUN=',NFUN,' exceeds NFUNMX=',NFUNMX
999   WRITE(6,*) 'NMKS --- Execution is stopped'
      STOP
      END
C NMMN
C
C Computation of the Mn matrix, output in MN
C
C In new notations Mn = Ksn
C
C Only the upper triangular is defined
C
C AR is the anisotropy ratio AR=r=Vthp/Vthl
C USV is the ratio of drift to "thermal" velocity  USV=Ui/Vi
C
C
      SUBROUTINE NMMN(N,BN,KPRL2,KPRL,AC,MN,AR,USV)

C                                     AC = OMEGAC/(K// V0)
C
C
      IMPLICIT NONE
C
      INTEGER    N
      REAL*8    AC,AR,BN,NA,USV,YN
      COMPLEX*16 CI,GNL(2),HNLK(4,3),KPRL,KPRL2,MN(3,3),Z,Z2
C
      DATA CI/(0.D0,1.)/
C
C
C
        IF(CDABS(KPRL2).EQ.0.D0) THEN
        Z=1.D-20
        Z2=1.D-40
        ELSE
        Z=KPRL
        Z2=KPRL2
        ENDIF
      YN=BN*AR
      NA=N*AC/AR
C
      CALL NMHNLK(N,YN,Z2,HNLK,GNL)
C
      Z=AR/Z
      Z2=Z*Z
C
      MN(1,1)=          (HNLK(2,1) + NA * HNLK(1,1))
      MN(1,2)=          (HNLK(2,2) + NA * HNLK(1,2))
      MN(2,2)= (      2*(HNLK(2,3) + NA * HNLK(1,3))
     .         -N**2 *          MN(1,1)
     .         +          GNL(2)   + NA *  GNL(1)    )*Z2
      MN(1,3)=          (HNLK(3,1) + NA * HNLK(2,1))
      MN(2,3)=          (HNLK(3,2) + NA * HNLK(2,2))
      MN(3,3)=           HNLK(4,1) + NA * HNLK(3,1)
C
      IF(USV.NE.0.D0) THEN
        MN(3,3)=MN(3,3) + 2*USV*MN(1,3) + USV**2 * MN(1,1)
        MN(1,3)=MN(1,3) +   USV*MN(1,1)
        MN(2,3)=MN(2,3) +   USV*MN(1,2)
      ENDIF
C
      MN(1,1)=   N**2 * MN(1,1) *Z2
      MN(1,2)=   N*CI * MN(1,2) *Z2
      MN(1,3)=   N    * MN(1,3) *Z
      MN(2,3)=  -CI   * MN(2,3) *Z
C     WRITE(6,*)'MN(1,1),MN(1,2)=',MN(1,1),MN(1,2)
C
      RETURN
      END
C NMSJ0
C
C Computation of S sub j0 of formula Eq.(28) of EUR-FU/XII-80/87-78
C                         but with l=0
C
C
C The computation is made recursively when Y is kept fixed and
C                          J is incremented by 1 at each call.
C
C
      COMPLEX*16 FUNCTION NMSJ0  (J,Y)
C
      IMPLICIT NONE
C
      INTEGER    I,J,JO,M,IMAX,NBOUND,NTERM(12)
      REAL*8    ABSY,ABSYO,AUX,AUX1,AUXS,NMCPJ,RKLL,RKLLO,SJM(80),
     .           Y,YVEC(80),YBOUND(11),YO,Y2
      COMPLEX*16 CSCALP,IPI,SJ0,S00
C
      COMPLEX*16 KPERP
      REAL*8 OMEGA
      COMMON/WAVEVE/KPERP,RKLL,OMEGA
C
      DATA IPI/(0.,3.1415 92653 58979 D0)/,SJ0/(0.D0,0.D0)/
      DATA ABSYO,YO/2*1.D17/,RKLLO/0.D0/,
     .     YBOUND/1.D8,    1.D4,    4.D1,    6.31D0, 3.415D0,
     .            2.512D0, 1.848D0, 1.585D0, 1.446D0, 1.360D0, 1.300D0/
      DATA JO/0/,NBOUND/11/,
     .     NTERM/  1,  2,  5, 10, 15,
     .            20, 30, 40, 50, 60, 70, 80/
C
C
      ABSY=DABS(Y)
        IF(ABSY.GT.1.26D0) THEN
C       The expansion for large Y is used
          IF(ABSY.NE.ABSYO) THEN
C         new value of Y
C           determine IMAX= (max value of m) +1
            DO 200 I=1,NBOUND
              IF(ABSY.GT.YBOUND(I)) THEN
              M=I
              GOTO 250
              ENDIF
200         CONTINUE
          M=NBOUND+1
250       IMAX=NTERM(M)
C         Build the two sequences YVEC= y**-2m and SJM=sjm
          YVEC(1)=1.D0
          SJM(1)=NMCPJ(J,0)
            IF(IMAX.GT.1) THEN
            Y2=1.D0/Y**2
              DO 300 I=2,IMAX
              M=I-1
              YVEC(I)=Y2*YVEC(I-1)
300           SJM(I)=NMCPJ(J,M)
            ENDIF
            YO=Y
            ABSYO=ABSY
          ELSE
C         Old value of y
            IF(J.EQ.JO) GOTO 550
            IF(J.EQ.JO+1) THEN
C           compute sjm recursively
            AUX=JO+1
            AUX1=JO+0.5
            DO 350 I=1,IMAX
350         SJM(I)=AUX*SJM(I)/(AUX1+I)
            ELSE
            SJM(1)=NMCPJ(J,0)
              IF(IMAX.GT.1) THEN
                DO 400 I=2,IMAX
400             SJM(I)=NMCPJ(J,I-1)
              ENDIF
            ENDIF
          ENDIF
          SJ0=-CSCALP(SJM,YVEC,IMAX)/Y
          JO=J
C
        ELSE
C       abs(y) is lower than 1.26
          IF(Y.NE.YO.OR.RKLL.NE.RKLLO) THEN
C         new value of y (r k// has changed)
C         Construct S00
          AUX=1-Y**2
          AUX1=Y+1
            IF(AUX1.EQ.0.D0) THEN
            AUX1=1.D-30
            AUX=0.D0
            IF(J.EQ.0)
     .       WRITE(6,*) ' NMSJ0 warning ****** Y=-1 !! S00 is infinite!'
            ENDIF
          AUX1=(Y-1)/AUX1
            IF(AUX1.EQ.0.D0) THEN
            AUX1=1.D-30
            AUX=0.D0
            IF(J.EQ.0)
     .       WRITE(6,*) ' NMSJ0 warning ****** Y=+1 !! S00 is infinite!'
            ENDIF
          S00=DLOG(ABS(AUX1))
            IF(ABSY.LT.1) THEN
              IF(RKLL.GT.0.D0) THEN
              S00=S00+IPI
              ELSE
              S00=S00-IPI
              ENDIF
            ENDIF
            IF(J.GT.0) THEN
C           Computation from the expansion formula in (1-y2)
            AUX1=1.D0
            AUXS=0.D0
              DO 450 M=1,J
              AUXS=AUXS+AUX1*NMCPJ(J-M,0)
450           AUX1=AUX1*AUX
            SJ0=-Y*AUXS+AUX1*S00
            ELSE
            SJ0=S00
            ENDIF
          ABSYO=ABSY
          YO=Y
          RKLLO=RKLL
          ELSE
C         previous y value
            IF(J.EQ.JO) GOTO 550
            IF(J.EQ.JO+1) THEN
            SJ0=AUX*SJ0-Y*NMCPJ(JO,0)
            ELSE
C           Computation from the expansion formula in (1-y2)
            AUX1=1.D0
            AUXS=0.D0
              DO 500 M=1,J
              AUXS=AUXS+AUX1*NMCPJ(J-M,0)
500           AUX1=AUX1*AUX
            SJ0=-Y*AUXS+AUX1*S00
            ENDIF
          ENDIF
        ENDIF
C
      JO=J
550   NMSJ0=SJ0
C     WRITE(6,*)'NMSJ0=',NMSJ0
      RETURN
      END
C
C SCALAR PRODUCT OF TWO COMPLEX VECTORS
C
      COMPLEX*16 FUNCTION CSCALP  (V1,V2,NDIM)
C
      IMPLICIT NONE
C
      REAL*8 V1(1),V2(1)
      INTEGER I,NDIM
C
      CSCALP=0.
        IF(NDIM.GT.0) THEN
          DO 200 I=1,NDIM
200       CSCALP=CSCALP+V1(I)*V2(I)
        ENDIF
      RETURN
      END
C NMBN0
C
C Computation of b sub n0 of formula Eq.(7) of EUR-FU/XII-80/87-78
C
C bn0= 1 / { 2**2n  (n!)**2 }
C
C
C
      DOUBLE PRECISION FUNCTION NMBN0  (N)
C
      IMPLICIT NONE
C
      INTEGER I,N,NA,NNO,NNA,NO
      REAL*8 BN0I,X,dgamma
C
      DATA BN0I/1.D0/,NO/0/
C
      NA=IABS(N)
      NNO=NA-NO
        IF(NNO.NE.0) THEN
        NNA=IABS(NNO)
          IF(NNA.GT.10) THEN
          X=NA+1
          BN0I=2.**(-2*N)/dgamma(X)**2
          ELSE
            IF(NNO.GT.0) THEN
C           Forward recursion
            X=2*NO
              DO 200 I=1,NNA
              X=X+2.D0
200           BN0I=BN0I/X**2
            ELSE
C           Backward recursion
            X=2*NO
              DO 250 I=1,NNA
              BN0I=BN0I*X**2
250           X=X-2.D0
            ENDIF
          ENDIF
        ENDIF
C
      NMBN0=BN0I
      NO=NA
      RETURN
      END
C NMCPJ
C
C COMPUTATION OF C(P,J)= GAMMA(P+1)*GAMMA(J+1/2)/GAMMA(P+J+3/2)
C
C THE COMPUTATION IS NORMALLY MADE RECURSIVELY STARTING FROM
C THE EARLIER COMPUTED VALUE CPJO UNLESS THE PATH FROM THE OLD
C (IPO,IJO) VALUES TO THE PRESENT ONES (IP,IJ) IS TOO LONG
C (MORE THAN 100 STEPS IN TOTAL). IN THIS LAST CASE CPJ IS COMPUTED
C DIRECTLY FROM THE GAMMA FUNCTION.
C
C
      DOUBLE PRECISION FUNCTION NMCPJ  (IP,IJ)
C
      IMPLICIT REAL *8 (A-H,O-Z)
C
      DATA IPO,IJO,CPJO/0,0,2.D0/,ICNT/0/
C
      IF(ICNT.GT.1000.AND.IP+IJ .LT. 56) GOTO 500
      IPD=IP-IPO
      IJD=IJ-IJO
      IF(IABS(IPD)+IABS(IJD).GT.40.AND.IP+IJ.LT.56) GOTO 500
C
150   IF(IPD) 200,250,160
160   IPO=IPO+1
      DO 180 I=IPO,IP
180   CPJO=CPJO*I/(IJO+I+0.5D0)
      GOTO 250
200   IPD=-IPD
      DO 220 I=1,IPD
      CPJO=CPJO*(IJO+IPO+0.5D0)/IPO
220   IPO=IPO-1
250   CONTINUE
C
      IF(IJD) 400,330,300
300   IJO=IJO+1
      DO 310 I=IJO,IJ
310   CPJO=CPJO*(I-0.5D0)/(I+IP+0.5D0)
330   ICNT=ICNT+1
340   NMCPJ=CPJO
        IF(NMCPJ.NE.0.D0) THEN
        IJO=IJ
        IPO=IP
        ELSE
        IJO=0
        IPO=0
        CPJO=2.D0
        ENDIF
      RETURN
C
400   IJD=-IJD
      DO 410 I=1,IJD
      CPJO=CPJO*(IJO+IP+0.5D0)/(IJO-0.5D0)
410   IJO=IJO-1
      GOTO 330
500   ICNT=0
      CPJO=dgamma(IP+1.D0)*dgamma(IJ+0.5D0)/dgamma(IJ+IP+1.5D0)
      GOTO 340
      END
C NMHNLK
C
C Computation of the functions
C            k
C           Hnl = Sum bnm (kprl2)**(n+m) (n+m)**k S(n+m)l
C
C    and    Gnl = Sum bnm (kprl2)**(n+m+1) S(n+m+1)l
C
C with k=0,1,2   and l=0,1,2,3 ; n stands for iabs(n)
C                            k
C The outputs are stored in Hnl = HNLK(l+1,k+1)
C                       and Gnl = GNL(l+1)
C
C
      SUBROUTINE NMHNLK(N,BN,KPRL2,HNLK,GNL)
C
      IMPLICIT NONE
C
      INTEGER    I,IERR,ISTOP,I1,J,MM,MMAX,N,NABS
      REAL*8    AUX,BN,ERR,NMBN0,NMCPJ,RELER2
      COMPLEX*16 CAUX,DGNL(2),DHNLK(4,3),GNL(2),HNLK(4,3),
     .           KPRL2,KP2J,NMSJ0,SJ(4),SJ1(2)
C
      PARAMETER (RELER2=1.D-16)
C
      DATA MMAX/101/
C          MMAX-1 is the maximum order in kprl2 in the summation over m.
C
C
C
      NABS=IABS(N)
      J=NABS
        IF(CDABS(KPRL2).EQ.0.D0.AND.J.EQ.0) THEN
        KP2J=1.D0
        ELSE
        KP2J=(KPRL2**J)*NMBN0(NABS)
        ENDIF
      DO 200 I=1,4
      DO 200 I1=1,3
200   HNLK(I,I1)=0.
      DO 220 I=1,2
220   GNL (I)   =0.
      ISTOP=0
C
C Loop: summation over m=MM-1
C
        DO 500 MM=1,MMAX
        SJ(1)=NMSJ0(J,BN)
        SJ(2)=NMCPJ(J,0) + BN*SJ(1)
        SJ(3)=             BN*SJ(2)
        SJ(4)=NMCPJ(J,1) + BN*SJ(3)
        SJ1(1)=NMSJ0(J+1,BN)
        SJ1(2)=NMCPJ(J+1,0) + BN*SJ1(1)
C
          DO 300 I=1,4
300       DHNLK(I,1)= KP2J * SJ(I)
        CAUX=J
          DO 310 I=1,3
310       DHNLK(I,2)= CAUX * DHNLK(I,1)
          DO 320 I=1,2
320       DHNLK(I,3)= CAUX * DHNLK(I,2)
C
        KP2J=KP2J*KPRL2
          DO 330 I=1,2
330       DGNL(I)= KP2J * SJ1(I)
C
          DO 340 I=1,4
            DO 340 I1=1,3
340         HNLK(I,I1)= HNLK(I,I1) + DHNLK(I,I1)
          DO 350 I=1,2
350       GNL (I)   = GNL(I) + DGNL(I)
C       Convergence tests
        ERR=0.D0
        IERR=0
        AUX=HNLK(2,3)*CONJG(HNLK(2,3))
          IF(AUX.NE.0.D0) THEN
          ERR = ERR + DHNLK(2,3)*CONJG(DHNLK(2,3))/AUX
          IERR=IERR+1
          ENDIF
        AUX=HNLK(3,2)*CONJG(HNLK(3,2))
          IF(AUX.NE.0.D0) THEN
          ERR = ERR + DHNLK(3,2)*CONJG(DHNLK(3,2))/AUX
          IERR=IERR+1
          ENDIF
        AUX=GNL(1)*CONJG(GNL(1))
          IF(AUX.NE.0.D0) THEN
          ERR = ERR + DGNL(1)*CONJG(DGNL(1))/AUX
          IERR=IERR+1
          ENDIF
C
          IF(ERR.LT.RELER2) THEN
          IF(ISTOP.GT.2) GOTO 520
          ISTOP=ISTOP+1
          ELSE
          ISTOP=0
          ENDIF
C
        AUX = ((NABS+MM-0.5D0)/(2*NABS+MM))/(MM*(NABS+MM))
        KP2J = -KP2J * AUX
500     J=J+1
      I=MMAX-1
      WRITE(6,*) '  NMHNLK --- series are not converged at order ',
     .           I,' in kprl2'
      WRITE(6,*) '  ****** --- relative error= ', ERR
520   CONTINUE
      RETURN
      END
