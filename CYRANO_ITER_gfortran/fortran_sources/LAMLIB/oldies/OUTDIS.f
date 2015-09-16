      SUBROUTINE OUTDIS

      IMPLICIT NONE

C     BUILDS OUTPUT TABLES, to be plotted elsewhere.
 
      include 'PARDIM.COPY'
      include 'DYNOU2.COPY'
      include 'COMPOW.COPY'
      include 'COMGEO.COPY'
      include 'COMANT.COPY'
      include 'COMMOD.COPY'
      include 'COMREG.COPY'
      include 'COMSUB.COPY'
      include 'COMSWE.COPY'
      include 'COM3DI.COPY'
      include 'COMFOU.COPY'
      include 'COMRO2.COPY'
      include 'COMPLO.COPY'
      include 'complp.copy'
      include 'COMMAG.COPY'
      include 'COMMA2.COPY'
      include 'COMPLA.COPY'
      include 'COMFIN.COPY'
      include 'COMFIC.COPY'
      include 'COMEQU.COPY'
      include 'COMIN2.COPY'
      include 'COMPHY.COPY'
       
      CHARACTER*9  CHARFI(2), CHARCO(6), CHARPA(3)
 
      INTEGER 
     ;  I, J, L
     ;, IPLOM1, IPOL, NJF
     ;, JSH, IREL
     ;, IDECA, NPLDIS
     ;, NROM
     ;, ian, nan, nppp
     ;, NSUM, I1ST, CEILQ
     ;, IDAMIN, IDMIN, IDMAX, isrchfge
C     ;, K, IPLO, NKUP, NIPHI, INCPLO, IS1, IS2, IPL1
C     ;, IPL2, IPL3, IPL4, NKULO, IKUD, IKU
C     ;, INCRHO, IEQU, IPHI, IEQ2, NROOTP

      DOUBLE PRECISION 
     ;  RLE, THFAC
     ;, FNIN1I, FNPLO
     ;, TABMIN, TABMAX, TABMI2, TABMA2, ratp(maxthp+1), yatp(maxthp+1)
     ;, FCO, FLH, CABS2, RVENO
     ;, polst
     ;, X1(4,npfft+1)
c     ;, X2(npfft+1)

      DOUBLE PRECISION 
     ;  ST2T2I
     ;, TRAPEZE, GAUSS
 
      COMPLEX*16
     ;  ZDOTC, zdotu
     ;, TRA1
 
      EXTERNAL 
     ;  NSUM, I1ST, CEILQ, TRAPEZE, GAUSS, FCO, FLH, RVENO
     ;, ZDOTC, ZDOTU, ZSET, DSET, DSUM, mucrvz
     ;, IDAMIN, IDMIN, IDMAX, isrchfge
 
      PARAMETER(NROM=3)

C      CHARACTER*7 JOBNAM
      CHARACTER*48 TEXT(14), title, yla(maxspe), yladis(2*nrom*maxpom)
 
      DOUBLE PRECISION 
     ;  SECOND

      EXTERNAL SECOND
C
      DATA 
     ;  CHARFI/'Electric ','Magnetic '/
     ;, CHARCO/'radial ','poloidal ','toroidal '
     ;, 'L.H.(+) ','R.H.(-) ','parallel '/
     ;, CHARPA/'real','imag','modulus'/
 
      CABS2(TRA1) = DREAL(TRA1 * DCONJG(TRA1))
 
c	----------------------- MAIN ------------------------------

cERN	NB : The (R,Z) grid for the 2D OUTPUT was already built and
c	     written to Rcyr.dat and Zcyr.dat in subroutine OUTGRID.f.
c		 All relevant variables are passed from OUTGRID through 
c		 COMMONS (nploth, polplo, istp(ireg), ist11, etc... )
c		 Respective computations were commented with cERN
 
      WRITE(NOFILE,*)'ENTER OUTDIS ; time=',SECOND()-TIMIN
C                     ------------
 
      ST2T2I = 0.5D0 * SQRT2I
c      CALL GRCHRC(0.2,0.,IDUM)
      
      do i = 1, nspec
      yla(i) = '                                                '
      end do
      do i = 1, 2*nrom*maxpom
      yladis(i) = '                                                '
      end do
 
      IMOTO = 1
      IF(CYL)THEN
      KPHI = KTOAN(IMOTO)
      ELSE
      N = MOTOAN(IMOTO)
      KPHI = N * R0I
      END IF
  
C     Index of current picture:
      IGLPLO = 0
      IPLOM1 = nploth + 1
 
C     GENERATE PLOT ABSCISSAE:
C     ------------------------(EXTRA POINTS=REGION BOUNDARIES)
c     new:
c     impose output at Gauss points to avoid radial interpolations:
cERN      OUTGAU = .true.

C     total number of points:
cERN       	IF(OUTGAU)THEN
c       NEW: absout becomes identical to abscis, absono to abscno:
cERN       	IST11 = nabsci
c      	IST11 = NELE * NGAUSS + 2 * NREG
cERN       	ELSE
cERN       	IST11 = NELE * ( NINTER + 1 ) + NREG
cERN       	END IF
cERN      FNIN1I = 1. / FLOAT(NINTER + 1)

      FNPLO = FLOAT(NPLOTH)
      THFAC = twopi / FNPLO
 
c      I = 0
cERN       DO 20 IREG = 1, NREG
C     --------------------
C       Number of points in each region:
cERN       	IF(OUTGAU)THEN
cC	    Output at Gauss points for each element + region boundaries:
c      	ISTP(IREG) = 2 + NGAUSS * (ILAEL(IREG)-IFIEL(IREG)+1)
c       Output at Gauss points and element boundaries; 
c      	internal region boundaries twice:
cERN       	ISTP(IREG) = (NGAUSS + 1) * (ILAEL(IREG) - IFIEL(IREG) + 1) + 1
c      	I = I + 1
c      	intaplot(i) = IREG + (ifiel(ireg)-1)*(NGAUSS+1)
c      	ABSR(I) = eqt(intaplot(i),1,1)
c      	ABSOUT(I) = RX0M(IREG-1)
c      	ABSONO(I) = RX0(IREG-1)
c      	  DO IEL = IFIEL(IREG), ILAEL(IREG)
c      	  JSH = IREG + (IEL-1)*(NGAUSS+1)
c      		DO J = 1, NGAUSS
c      		I = I + 1
c      	    intaplot(i) = j + jsh
c          	ABSR(I) = eqt(intaplot(i),1,1)
c      		ABSONO(I) = ABSCNO(JSH+J)
c      		ABSOUT(I) = ABSCIS(JSH+J)
c            END DO
c          END DO
        
c      	ELSE
cC	    NINTER equally spaced internal points + boundaries in each element:
c      	obsolete!
c      	ISTP(IREG) = 1 + (NINTER + 1) * (ILAEL(IREG) - IFIEL(IREG) + 1)
c      	  DO IEL = IFIEL(IREG), ILAEL(IREG)
c      	  RLE = FL(IEL)
c      		DO J = 0, NINTER
c      		I = I + 1
c      		ABSONO(I) = FX0(IEL-1) + J * FNIN1I * RLE
c      		ABSOUT(I) = ABSONO(I) * RNORM
c            END DO
c          END DO
cERN       	END IF
c      I = I + 1
c      intaplot(i) = intaplot(i-1) + 1
c      ABSR(I) = eqt(intaplot(i),1,1)
c      ABSOUT(I) = RX0M(IREG)
c      ABSONO(I) = RX0(IREG)
cERN   20  CONTINUE
C     --------
c  New:
cERN         if(outgau)then
cERN         call dcopy(ist11, abscis, 1, absout, 1)
cERN         call dcopy(ist11, abscno, 1, absono, 1)
cERN         call dcopy(ist11, eqt(1,1,1), 1, absr, 1)
cERN          do i = 1, ist11
cERN           intaplot(i) = i
cERN           end do
cERN         end if
      
c     Poloidal angles for 2D plots:
cERN       POLST = twopi / fnplo
cERN       POLPLO(IPLOM1) = twopi
cERN         do IPOL = 1, NPLOTH
cERN         POLPLO(IPOL) = (ipol - 1) * polst
cERN         end do

c     Write 2D table of R and Y at plot points:
      write(76,*)ist11
      write(76,*)iplom1
c     Elongation (to open graphic window):
        if(circ)then
        write(76,*)1.d0
        else
        write(76,*)kappa
       end if
        if(.not.updsym)then
C       Interpolate equilibrium at plot points, no up-down symmetry.
C       This is only correct if OUTGAU=.T.! Otherwise 2D interpol.
          do i = 1, ist11
          INTAB = intaplot(i)
          CALL INTERP1(POLANG, POLANG(NPFFT+1), eqt(INTAB,1,1), NABPLO,
     ;    NPFFT+1, POLPLO, ratp, iplom1, 'R')
          CALL INTERP1(POLANG, POLANG(NPFFT+1), eqt(INTAB,1,2), NABPLO,
     ;    NPFFT+1, POLPLO, yatp, iplom1, 'R')
          write(76,2000)(ratp(j), yatp(j), j = 1, iplom1)
          end do
        else
C       Use up-down symmetry:
C       This is only correct if OUTGAU=.T.! Otherwise 2D interpol.
        IAN = NPLOTH / 2 + 1
        NAN = IPLOM1 - IAN + 1
          do i = 1, ist11
          INTAB = intaplot(i)
          CALL INTERP1(POLANG, POLANG(NPFFT/2+1), eqt(INTAB,1,1), NABPLO,
     ;    NPFFT/2+1, POLPLO, ratp, IAN, 'R')
c          CALL DCOPY(ian, ratp, 1, ratp(nan), -1)
            do j = 1, ian
            ratp(iplom1-j+1) = ratp(j)
            end do
          CALL INTERP1(POLANG, POLANG(NPFFT/2+1), eqt(INTAB,1,2), NABPLO,
     ;    NPFFT/2+1, POLPLO, yatp, IAN, 'R')
            do j = 1, ian
            yatp(iplom1-j+1) = - yatp(j)
            end do
          write(76,2000)(ratp(j), yatp(j), j = 1, iplom1)
          end do
        end if
 
      IF(.NOT.VACUUM(1))THEN
c      CALL DSET((MAXTHP+1)*NABPLO*2, -75.75D20, TAB, 1)
      CALL DSET((MAXTHP+1)*ISTP(1), 0.d0, TAB, 1)
      CALL DSET((MAXTHP+1)*ISTP(1), 0.d0, TAB(1,1,2), 1)
C     Table to locate cut-offs and resonances:   IN CENTRAL REGION ONLY
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C                                         (PROBLEMS WITH JUMPS!)
c     Resp. L-n//2, R-n//2, S-n//2 using real parts and
c     n//=nphi (refractive index of m=0 mode)
      i = 0
      npp = npfft + 1
      nppp = iplom1
        if(updsym)then
        npp = npft2 + 1
        nppp = nploth / 2 + 1
        end if
         
      do ireg = 1, 1
        do irel = 1, istp(ireg)
        i = i + 1
cERN        intab = i
        intab = intaplot(i)
          do ipol = 1, npp
          intabp = ipol
          call charsu(x1(1,ipol))
          end do
          do l = 1, 3
          call interp1(polang, polang(npp), x1(l,1), 4,
     ;    npp, polplo, tab(1,i,l), nppp, 'R')
            if(updsym)then
              do j = 1, ian
              tab(iplom1-j+1,i,l) = tab(j,i,l)
              end do
            end if
          end do
          
        end do
      end do
      end if
 
      fnplo = dfloat(nploth)
 
      if(.not.replay)then
c     plot equilibrium profiles and local dispersion:
c     absdis, abspro built in groots
      npldis = nele + nreg
      if(.not.polsym)npldis = 2 * (nele + nreg)
 
cERN      if(.not.glovac)then
cERNc     -------------------      

cERN      	IGLPLO = IGLPLO + 1
cERN      	title = 'Density profiles                                '
cERN      	WRITE(nofile,1000) IGLPLO, title
cERN      	PLOLEG(IGLPLO) = title
cERN        call wrp1d(70, IGLPLO, NPLDIS, ABSDIS, PROFIL(1,1)
cERN     ;, 2*(MAXNEL+MAXREG), 1, NSPEC
cERN     ;, title
cERN     ;, 'm         ', yla
cERN     ;, 'm**-3     ', PANAME
cERN     ;, 0, TEXT,.false.)

cERN        if(.not.glocol)then
cERN      	IGLPLO = IGLPLO + 1
cERN      	title = 'Temperature profiles                            '
cERN      	WRITE(NOFILE,1000) IGLPLO, title
cERN      	PLOLEG(IGLPLO) = title
cERN        call wrp1d(71, IGLPLO, NPLDIS, ABSDIS, PROFIL(1,NSPEC+1)
cERN     ;, 2*(MAXNEL+MAXREG), 1, NSPEC
cERN     ;, title
cERN     ;, 'm         ', yla
cERN     ;, 'eV        ', PANAME
cERN     ;, 0, TEXT,.false.)
cERN        end if
        
cERN      end if
cERNc     ------   

cERN      	IGLPLO = IGLPLO + 1
cERN      	title = 'Magnetic induction profile                      '
cERN      	WRITE(nofile,1000) IGLPLO, title
cERN      	PLOLEG(IGLPLO) = title
cERN        call wrp1d(72, IGLPLO, NPLDIS, ABSDIS, PROFIL(1,2*NSPEC+1)
cERN     ;, 2*(MAXNEL+MAXREG), 1, 1
cERN     ;, title
cERN     ;, 'R, m      ', yla
cERN     ;, 'Tesla     ', yla
cERN     ;, 0, TEXT,.false.)

cERN        IF(IPL.NE.0.D0)then
cERN        IGLPLO = IGLPLO + 1
cERN        title = 'q profile                                       '
cERN        WRITE(nofile,1000) IGLPLO, title
cERN        PLOLEG(IGLPLO) = title
cERN        write(73,*)title
cERN        write(73,*)nabsci, 2
cERN          DO i = 1, nabsci
cERN          write(73,2000)abscis(i), eqta1d(i,4)
cERN			write(73,2000)abscis(i), qfactor(i)
cERN          END DO
c        call wrp1d(73, IGLPLO, NPLDIS, ABSPRO, PROFIL(1,2*NSPEC+2)
c     ;, 2*(MAXNEL+MAXREG), 1, 1
c     ;, title
c     ;, 'm         ', yla
c     ;, '          ', yla
c     ;, 0, TEXT,.false.)

cERN        IGLPLO = IGLPLO + 1
cERN        title = 'Flux function                                   '
cERN        WRITE(nofile,1000) IGLPLO, title
cERN        PLOLEG(IGLPLO) = title

cERN        write(73,*)title
cERN        write(73,*)nabsci, 2
cERN          DO i = 1, nabsci
cERN			write(73,2000)abscis(i), EQTA1D(i,7)
cERN			write(73,2000)abscis(i), Psir(i)
cERN          END DO
cERN        END IF

cERN      if(cokpco)then
cERN      IGLPLO = IGLPLO + 1
cERN      title = 'H length                                     '
cERN      WRITE(nofile,1000) IGLPLO, title
cERN      PLOLEG(IGLPLO) = title
cERN      write(73,*)title
cERN      write(73,*)nabsci, 2
cERN        DO i = 1, nabsci
cERN        write(73,2000)abscis(i), hachi(i)
cERN        END DO
cERN      end if
 
c      NPLDIS = NELE + NREG - 1
      NPLDIS = NELE + NREG
      IF(.NOT.POLSYM)NPLDIS = 2 * NPLDIS
      NJF = 8
      IF(GLO0OR)NJF = 6

      DO J = 1, NJF * ((MSTUD2 - MSTUD1) / MSTUST + 1)
      CALL SCALOG(NPLDIS, ROOTS(1,J), 2)
      END DO

cERN  TO BE REMOVED ################################	
c     All root real and imag parts:
      IGLPLO = IGLPLO + 1
      title = 'Dispersion roots                               '
      PLOLEG(IGLPLO) = title
      WRITE(NOFILE,1000) IGLPLO, title
      i = 2*NROM*((mstud2-mstud1)/mstust+1)
      call wrp1d(74, IGLPLO, NPLDIS, ABSDIS, ROOTS(1,1)
     ;, 2*(MAXNEL+MAXREG), 1, i
     ;, title
     ;, 'm         ', yladis
     ;, '          ', 'Log scale, linear near 0                        '
     ;, 0, TEXT, .false.)

c      DO 1116 M = MSTUD1, MSTUD2, MSTUST
c      MR = (M - MSTUD1) / MSTUST + 1
c      IKUD = NROM*(MR-1)
c      X2 = RX0M(NREG)
c      	IF(POLSYM)THEN
c        NIPHI = 1
c        IPLO = 0
c        X1 = 0.
c     	ELSE
c      	NIPHI = 2
c      	IPLO = NPLDIS/2
c      	X1 = - RX0M(NREG)
c      	END IF
c      NKUP = 0
c
c      IGLPLO = IGLPLO + 1
c      title = 'Dispersion: roots real part                     '
c      PLOLEG(IGLPLO) = title
c      WRITE(NOFILE,1000) IGLPLO, title

c      DO 112 IPHI = 1, NIPHI
c      INCPLO = 3 - 2 * IPHI
c      	DO IREG = 1, NREG
c      		IF(VACUUM(IREG))THEN
c      		NROOTP = 1
c      		ELSE IF(COLDPL(IREG) .OR. .NOT.FLROPS(IREG))THEN
c      		NROOTP=2
c      		ELSE
c      		NROOTP=3
c      		END IF
c      	IS1 = NSUM(IREG-1,NS)+1
c      	IS2 = NSUM(IREG,NS)
c      	IPL1 = IPLO+INCPLO
c      	IPLO = IPLO + NSUM(NS(IREG), IELE(IS1))*INCPLO
c      	IF(IREG.GT.1)IPLO = IPLO + INCPLO
c      	IPL2 = IPLO
c      	IPL3 = MIN0(IPL1,IPL2)
c      	IPL4 = MAX0(IPL1,IPL2)
c      	NKULO = NROOTP
c      	NKUP = NKUP + NKULO
c TEMPORAIRE: ECRIT PAR REGION (NOMBRE DE RACINES EST VARIABLE)
c        call wrp1d(NOFILE,IGLPLO,IPL4-IPL3+1,ABSDIS(IPL3),ROOTS(IPL3,2*IKUD+1)
c     ; ,2*(MAXNEL+MAXREG),2,NKULO,TEXT(13),' ',TEXT(2),TEXT(9),' ',0,TEXT
c     ; ,.false.)
c      	END DO
c     	IPLO = NPLDIS/2 + 1
c 112  CONTINUE
C TEMPORAIRE: ECRIT POUR 1 RACINE
C      call wrp1d(NOFILE, IGLPLO, NPLDIS, ABSDIS, ROOTS(1,2*IKUD+1)
C     ;, 2*(MAXNEL+MAXREG), 1, 1
C     ;, title
C     ;, 'm         ', yla
C     ;, '          ', 'Log scale, linear near 0                        '
C     ;, 0, TEXT,.false.)
C
c 1116 CONTINUE

c      DO 1115 M = MSTUD1, MSTUD2, MSTUST
c      MR = (M - MSTUD1)/MSTUST + 1
c      IKUD = NROM*(MR-1)
c      IPLO = 0
c      NIPHI = 1
c      	IF(.NOT.POLSYM)THEN
c      	NIPHI = 2
c      	IPLO = NPLDIS/2
c      	END IF
c      NKUP = 0
c
c      IGLPLO = IGLPLO + 1
c      title = 'Dispersion: roots imag part                     '
c      PLOLEG(IGLPLO) = title
c      WRITE(NOFILE,1000) IGLPLO, title

c      DO 113 IPHI=1,NIPHI
c      INCPLO = 3-2*IPHI
c      DO 13 IREG = 1, NREG
c      	IF( VACUUM(IREG))THEN
c      	NROOTP = 1
c      	ELSE IF(COLDPL(IREG) .OR. .NOT.FLROPS(IREG))THEN
c      	NROOTP=2
c      	ELSE
c      	NROOTP=3
c      	END IF
c      IS1 = NSUM(IREG-1,NS)+1
c      IS2 = NSUM(IREG,NS)
c      IPL1 = IPLO+INCPLO
c      IPLO = IPLO + NSUM(NS(IREG), IELE(IS1))*INCPLO
c      IF(IREG.GT.1)IPLO = IPLO + INCPLO
c      IPL2 = IPLO
c      IPL3 = MIN0(IPL1,IPL2)
c      IPL4 = MAX0(IPL1,IPL2)
cC     NKULO = NROOTP*(MSTUD2-MSTUD1+1)
c      NKULO = NROOTP
c      NKUP = NKUP + NKULO
c      IKUD = NROOTP*(MR-1)
cC TEMPORAIRE: ECRIT PER REGION (NOMBRE DE RACINES EST VARIABLE)
c        call wrp1d(NOFILE,IGLPLO,IPL4-IPL3+1,ABSDIS(IPL3),ROOTS(IPL3,2*IKUD+2)
c     ; ,2*(MAXNEL+MAXREG),2,NKULO,TEXT(13),' ',TEXT(2),TEXT(9),' ',0,TEXT
c     ;, .false.)
c  13  	CONTINUE
c      IPLO = NPLDIS/2 + 1
c 113  CONTINUE
C TEMPORAIRE: ECRIT POUR 1 RACINE
C      call wrp1d(NOFILE, IGLPLO, NPLDIS, ABSDIS, ROOTS(1,2*IKUD+2)
C     ;, 2*(MAXNEL+MAXREG), 1, 1
C     ;, title
C     ;, 'm         ', yla
C     ;, '          ', 'Log scale, linear near 0                        '
C     ;, 0, TEXT, .false.)

c 1115 CONTINUE
 
C      IF(.NOT.GLO0OR)THEN
C
C      DO 117 M = MSTUD1, MSTUD2, MSTUST
C      MR = (M - MSTUD1) / MSTUST + 1
C      IKU = 2*NROM*( (MSTUD2-MSTUD1)/MSTUST+1 ) + 2 * MR - 1

c      IPLO = 0
c      NIPHI = 1
c      X1 = 0.
c      X2 = RX0M(NREG)
c      	IF(.NOT.POLSYM)THEN
c      	NIPHI = 2
c      	IPLO = NPLDIS/2
c      	X1 = - RX0M(NREG)
c      	END IF
c      NKUP = 0
c      DO 119 IPHI = 1, NIPHI
c      INCPLO = 3 - 2 * IPHI
c      	DO 120 IREG = 1, NREG
c      	IS1 = NSUM(IREG-1,NS)+1
c      	IS2 = NSUM(IREG,NS)
c      	IPL1 = IPLO+INCPLO
c      	IPLO = IPLO + NSUM(NS(IREG), IELE(IS1))*INCPLO
c      	IF(IREG.GT.1)IPLO = IPLO + INCPLO
c      	IPL2 = IPLO
c      	IPL3 = MIN0(IPL1,IPL2)
c      	IPL4 = MAX0(IPL1,IPL2)
c      	NKULO = 2
c      		IF( .NOT.VACUUM(IREG) .AND. .NOT.COLDPL(IREG) .AND. FLROPS(IREG)
c     ;		  )THEN
c      		NKUP = NKUP + NKULO
c        call wrp1d(NOFILE,IGLPLO,IPL4-IPL3+1,ABSDIS(IPL3),ROOTS(IPL3,iku)
c     ; ,2*(MAXNEL+MAXREG),2,NKULO,TEXT(13),' ',TEXT(2),TEXT(9),' ',0,TEXT)
c        call wrp1d(NOFILE,IGLPLO,IPL4-IPL3+1,ABSDIS(IPL3),ROOTS(IPL3,iku+1)
c     ; ,2*(MAXNEL+MAXREG),2,NKULO,TEXT(13),' ',TEXT(2),TEXT(9),' ',0,TEXT
c     ;, .false.)
c      		END IF
c 120  	CONTINUE
c      IPLO = NPLDIS/2 + 1
c 119  CONTINUE

C      IGLPLO = IGLPLO + 1
C      title = 'Krho max. * Larmor                              '
C      PLOLEG(IGLPLO) = title
C      WRITE(NOFILE,1000) IGLPLO, title
C      call wrp1d(NOFILE, IGLPLO, NPLDIS, ABSDIS, ROOTS(1,iku)
C     ;, 2*(MAXNEL+MAXREG), 1, 1
C     ;, title
C     ;, 'm         ', yla
C     ;, '          ', yla, 0, TEXT, .false.)
C
C      IGLPLO = IGLPLO + 1
c23456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
C      title = 'Keta * Larmor                                   '
C      PLOLEG(IGLPLO) = title
C      WRITE(NOFILE,1000) IGLPLO, title
C      call wrp1d(NOFILE, IGLPLO, NPLDIS, ABSDIS, ROOTS(1,iku+1)
C     ;, 2*(MAXNEL+MAXREG), 1, 1
C     ;, title
C     ;, 'm         ', yla
C     ;, '          ', yla, 0, TEXT, .false.)
C
C 117  CONTINUE
C      END IF
 
      end if  ! .not.replay
 
      IF(.NOT.VACUUM(1))THEN
C     Fast wave cut-offs and perpendicular resonance:
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      TABMIN = TAB2(idmin((MAXTHP+1)*ISTP(1), TAB, 1))
      TABMAX = TAB2(idmax((MAXTHP+1)*ISTP(1), TAB, 1))
      IDECA = (MAXTHP+1) * NABPLO
      TABMI2 = TAB2(IDECA + idmin((MAXTHP+1)*ISTP(1), TAB(1,1,2), 1))
      TABMA2 = TAB2(IDECA + idmax((MAXTHP+1)*ISTP(1), TAB(1,1,2), 1))
CERN        IF(TABMIN*TABMAX.LT.0. .OR. TABMI2*TABMA2.LT.0.)THEN
        IGLPLO = IGLPLO + 1
        WRITE(NOFILE,1000) IGLPLO, 'Characteristic surfaces             '
        PLOLEG(IGLPLO) = 'Characteristic surfaces             '

        write(77,*)istp(1)
        write(77,*)iplom1
c       Elongation (to open graphic window):
          if(circ)then
          write(77,*)1.d0
          else if(dshape)then
          write(77,*)kappa
	    else
c         Must provide a valid elongation value when geneq is true!!!
          write(77,*)kappa
          end if
          do i = 1, istp(1)
            do j = 1, iplom1
            write(77,2000)tab(j,i,1), tab(j,i,2), tab(j,i,3)
            end do
          end do

cERN       end if
      end if
      
      return
 1000 format(1H ,'Plot #',I4,'  ',A50)
 2000 format(1h ,2(2x,g14.6))
      end
