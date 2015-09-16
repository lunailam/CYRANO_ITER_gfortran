      SUBROUTINE POLFFT(V, LDV, ONESPE, ISP, gdrtoo, SWITCH, TRAFO, IPATH, nfft
     ;, onlyab)

      IMPLICIT NONE

      LOGICAL ONESPE, TRAFO, gdrtoo, onlyab

      INTEGER LDV, ISP, SWITCH, IPATH, nfft

      COMPLEX*16 V(LDV,19)

      include 'PARDIM.COPY'
      include 'COMSUB.COPY'
      include 'COMGEO.COPY'
      include 'COMMAG.COPY'
      include 'COMMA2.COPY'
      include 'COMANT.COPY'
      include 'COMMOD.COPY'
      include 'COMPLA.COPY'
      include 'COMIN2.COPY'
      include 'COMFIC.COPY'
      include 'COMSWE.COPY'
      include 'COMPHY.COPY'

c     29/10/2000: added argument ONLYAB: if .true., only use damping terms
c     in dielectric tensor (imag.part of plasma dispersion function) and
c     ignore reactive terms.
c
c     This version with loop over average M
c
c     Susceptibility tensor: - for species isp (onespe=.t.)
c                            - all species     (onespe=.f.)
c     Fast poloidal fourier transform : call to cray routine cfft2
c     This version uses vectorized plasma dispersion function.
c     Output:
c                 v(ldv,19) (with spectra valid up to nfft/2 th harmonic)

      LOGICAL MUL, TRAFOL, USUAL, vssh

      CHARACTER*3 ELETYP

      INTEGER
     ;  NFC, INDEX(19,2), NTR(2), ISP1, ISP2
     ;, IPO, IKIN, NFFT1, MFFT, NPOI, IR, IS, I, J, I2, J2, LRES, NTRANS
     ;, NSETZ, NFFTME

c     valid for up-down asymmetric equilibria:
      PARAMETER(NFC = NPFFT+1)

      DOUBLE PRECISION 
     ;  TRA1, TRA6
     ;, FFSTEP, ANGLES(NFC)
     ;, KPER, KPER2, KPAR
     ;, KPARV(NFC), KPERV(NFC)
     ;, TEM, TEM3, TEMK, TEMK2, co2, si2, s7

      COMPLEX*16 SIGSD(3,3,NFC), SIGOD(3,3,NFC), V1(19), WORK(5*NPFFT/2)
     ;, a, b, c

      DATA NFFTME/0/
      
      save nfftme, work
 
CC    WRITE(NOFILE,*)'ENTER POLFFT'
      eletyp = styp(isubr)
	  
      MAVE = DFLOAT(MAV2) * 0.5
 
c     IPATH to short-circuit index list building:
      IF(IPATH .NE. 1)GOTO 8

c     Don't FFT at magnetic axis:
      TRAFOL = TRAFO .AND. Y.GT.1.D-40
c     Following is true if a torus, not at magn. axis and TRAFO=.t.:
      MUL = TRAFOL .AND. .NOT.CYL
 
      RHO = abscis(intab)
      YINV = ABSCNI(INTAB)
      RHOINV = YINV * RNORI
        IF(ONESPE)THEN
        ISP1 = ISP
        ISP2 = ISP
        ELSE
        ISP1 = 1
        ISP2 = NSPEC
        END IF

c     Define lists of coefficients to poloidal-Fourier-transform:

      IF(SWITCH.EQ.1)THEN
C     +++++++++++++++++++
c     This case deals with volume terms

c       Label resonant and nonresonant dielectric tensor elements:      
        IF(GLOMAX)THEN
C       --------------
c       case Raymond Koch's nonmaxwellian tensor not in use
        
          IF(COLDPL(IREG) .OR. .NOT.FLROPS(IREG))THEN
c         Zeroth order in Larmor: only ±1 cyclotron resonances.

            if(circ .or. eletyp.eq.'HEC')then
c           3 terms in tensor:
c           - diagonal(L,R,P) in +-// coordinates
c           - elements S,±iD,P
C           Number of Fourier transforms:
            NTRANS = 3
            NTR(1) = 0
            NTR(2) = 0
              DO LRES = -1, 1
C               'RES' was computed in cyrano.f (LRES=-2 to 2). (It does not
C               detect Landau resonance, but the corresponding terms are smooth
C               poloidally.)
                IF(RES(LRES))THEN
c               Resonant terms:
                NTR(1) = NTR(1) + 1
c               Cyclotron resonance: index: value 2, 3, 1 for resp. -, //, +
c               associated with LRES=-1,0,1
                INDEX(NTR(1),1) = 3 - LRES*(3*LRES+1) / 2
                ELSE
c               Nonresonant terms:
                NTR(2) = NTR(2) + 1
                INDEX(NTR(2),2) = 3 - LRES*(3*LRES+1) / 2
                END IF
              END DO
            else
c           6 terms in tensor in rho, theta, phi coordinates:
            ntrans = 6
            NTR(1) = 0
            NTR(2) = 0
c           Quick, not optimized: terms assumed all resonant or nonresonant
              if(res(-1) .or. res(0) .or. res(1))then
              NTR(1) = 6
                do i = 1, 6
                INDEX(i,1) = i
                end do
              else
              NTR(2) = 6
                do i = 1, 6
                INDEX(i,2) = i
                end do
              end if
            end if

          ELSE
c         Dielectric tensor to second order in kperp*Larmor inclusive:
c         NB: here resonances are treated as mutually exclusive
          NTRANS = 19
c         The 19 coefficients are:
c@
            IF(RES(1))THEN
c           Fundamental ion cyclotron resonance is present
            NTR(1) = 12
            INDEX(1,1) = 1
            INDEX(2,1) = 3
            INDEX(3,1) = 4
            INDEX(4,1) = 6
            INDEX(5,1) = 7
            INDEX(6,1) = 8
            INDEX(7,1) = 9
            INDEX(8,1) = 10
            INDEX(9,1) = 11
            INDEX(10,1) = 14
            INDEX(11,1) = 16
            INDEX(12,1) = 17
            NTR(2) = 7
            INDEX(1,2) = 2
            INDEX(2,2) = 5
            INDEX(3,2) = 12
            INDEX(4,2) = 13
            INDEX(5,2) = 15
            INDEX(6,2) = 18
            INDEX(7,2) = 19
          
            ELSE IF(RES(2))THEN
c           Second harmonic ion cyclotron resonance is present
            NTR(1) = 3
            INDEX(1,1) = 4
            INDEX(2,1) = 14
            INDEX(3,1) = 17
            NTR(2) = 16
            INDEX(1,2) = 1
            INDEX(2,2) = 2
            INDEX(3,2) = 3
            INDEX(4,2) = 5
            INDEX(5,2) = 6
            INDEX(6,2) = 7
            INDEX(7,2) = 8
            INDEX(8,2) = 9
            INDEX(9,2) = 10
            INDEX(10,2) = 11
            INDEX(11,2) = 12
            INDEX(12,2) = 13
            INDEX(13,2) = 15
            INDEX(14,2) = 16
            INDEX(15,2) = 18
            INDEX(16,2) = 19

            ELSE
c           No fundamental or second harmonic ion cyclotron resonance is present
c           All coefficients are assumed poloidally smooth
            NTR(1) = 0
            NTR(2) = 19
              DO I = 1, NTR(2)
              INDEX(I,2) = I
              END DO
            END IF
          END IF
       
        ELSE
C       ----
c       Case GLOMAX=.f.: some species use R.Koch's nonmaxwellian tensor
        NTRANS = 6
          IF(RES(1))THEN
c         Fundamental ion cyclotron resonance is present
          NTR(1) = 1
          INDEX(1,1) = 1
          NTR(2) = 5
          INDEX(1,2) = 2
          INDEX(2,2) = 3
          INDEX(3,2) = 4
          INDEX(4,2) = 5
          INDEX(5,2) = 6
           
          ELSE
c         No fundamental ion cyclotron resonance is present
c         All coefficients are assumed poloidally smooth - check why!
          NTR(1) = 0
          NTR(2) = 6
          INDEX(1,2) = 1
          INDEX(2,2) = 2
          INDEX(3,2) = 3
          INDEX(4,2) = 4
          INDEX(5,2) = 5
          INDEX(6,2) = 6
          END IF
        END IF
C       ------

      ELSE IF(SWITCH.EQ.2)THEN
C     ++++++++++++++++++++++++
c     This case deals with boundary terms

      IF(COLJUM)THEN
c     'Cold-like' jump condition:
 
        IF(GLOMAX)THEN
C       --------------
c       No call to R.K.'s nonmaxwellian routines
          if(circ .or. eletyp.eq.'HEC')then
          NTRANS = 3
            IF(RES(1))THEN
            NTR(1) = 1
c           ++ term is resonant
            INDEX(1,1) = 1
            NTR(2) = 2
c           and -- and // terms are smooth
            INDEX(1,2) = 2
            INDEX(2,2) = 3
            ELSE
c           all terms are smooth
            NTR(1) = 0
            NTR(2) = 3
            INDEX(1,2) = 1
            INDEX(2,2) = 2
            INDEX(3,2) = 3
            END IF
          else
c         6 terms in tensor in rho, theta, phi coordinates:
          ntrans = 6
          NTR(1) = 0
          NTR(2) = 0
            if(res(-1) .or. res(0) .or. res(1))then
            NTR(1) = 6
              do i = 1, 6
              INDEX(i,1) = i
              end do
            else
            NTR(2) = 6
              do i = 1, 6
              INDEX(i,2) = i
              end do
            end if
          end if

        ELSE
C       ----
c       At least a species requires R.K.'s nonmaxwellian routines
        NTRANS = 6
          IF(RES(1))THEN
          NTR(1) = 1
          INDEX(1,1) = 1
          NTR(2) = 5
          INDEX(1,2) = 2
          INDEX(2,2) = 3
          INDEX(3,2) = 4
          INDEX(4,2) = 5
          INDEX(5,2) = 6

          ELSE
          NTR(1) = 0
          NTR(2) = 6
          INDEX(1,2) = 1
          INDEX(2,2) = 2
          INDEX(3,2) = 3
          INDEX(4,2) = 4
          INDEX(5,2) = 5
          INDEX(6,2) = 6
          END IF
        END IF
C       ------

      ELSE
c     There are flr terms in the jump conditions
      NTRANS = 9
        IF(RES(1))THEN
c       Fundamental ion cyclotron resonance is present
        NTR(1) = 6
        INDEX(1,1) = 1
        INDEX(2,1) = 2
        INDEX(3,1) = 3
        INDEX(4,1) = 4
        INDEX(5,1) = 7
        INDEX(6,1) = 9
        NTR(2) = 3
        INDEX(1,2) = 5
        INDEX(2,2) = 6
        INDEX(3,2) = 8

        ELSE IF(RES(2))THEN
c       Second harmonic ion cyclotron resonance is present, but no fundamental
        NTR(1) = 2
        INDEX(1,1) = 1
        INDEX(2,1) = 2
        NTR(2) = 7
        INDEX(1,2) = 3
        INDEX(2,2) = 4
        INDEX(3,2) = 5
        INDEX(4,2) = 6
        INDEX(5,2) = 7
        INDEX(6,2) = 8
        INDEX(7,2) = 9

        ELSE
c       No fundamental or second harmonic ion cyclotron resonance is present
c       All coefficients are assumed poloidally smooth
        NTR(1) = 0
        NTR(2) = 9
          DO I = 1, 9
          INDEX(I,2) = I
          END DO
        END IF
      END IF
      END IF
 
      IPATH  =  2
  8   CONTINUE

c     Initialize output table:
      CALL ZSET(LDV * NTRANS, CZERO, V, 1)
 
c     Loop over resonant and non resonant terms:
      DO 1000 IKIN = 1, 2

c     Choose number of poloidal points: 
      IF(TRAFOL)THEN
C     --------------

c     23/10/2000: simplified the logic:
        IF(IKIN.EQ.1 .OR. (NSCREE.GT.0 .AND. vssh))THEN
c       Use allowed maximum of Fourier points NPFFT for the resonant terms or if
c       a Faraday shield generates poloidal inhomogeneities:
c       (NPFFT must be a power of 2 and >=8)
        NFFT = NPFFT
        ELSE
c       Number of poloidal points for fft of smooth terms, NFTNR, is user-supplied:
        NFFT = min0(NFTNR, NPFFT)
        END IF
c--------------------------------------------------------------------------
c      Old mess: 
c
cc     At least 2 points per wavelength: 4 points per wl. are taken
c      NFFT1 = 4 * KLIM
cc     Case of a different plasma in screens, second source of poloidal
cc     inhomogeneity:
c      vssh = .false.
c
c        IF(NSCREE.GT.0)THEN
cC       Faraday shields with vacuum inside: take into account
cC       to find number of poloidal Fourier points.
cC       VSSH=.t. if there is vacuum in one of the screens.
c          DO IS = 1, NSCREE
c	      if(vscrsh(is))vssh = .true.
cc         10 points over poloidal extent of shield:
c          IF(vscrsh(is) .AND. IREG.GT.IRBSCR(IS))
c     ;    NFFT1 = DMAX1(DFLOAT(NFFT1), 10*(2.D0*PI/DTHES(IS)))
c          END DO
c        END IF
c
c        IF(IKIN.EQ.1 .OR. (NSCREE.GT.0 .AND. vssh))THEN
c        NFFT = NFFT1
cc       At least 8 Fourier points for correct call of library routine:
c        NFFT = MAX0(8, NFFT)
c        ELSE
cc       Number of pol.coeffs. for smooth terms, NFTNR, is user-supplied:
c        NFFT = NFTNR
c        END IF
c
c        if(NFFT.GT.NPFFT)then
c        WRITE(NOFILE,*)
c     ;  'WARNING: WORK SPACE TOO SHORT FOR FFT ;'
c     ;  //' SPECTRUM MAY BE POORLY TRUNCATED; NFFT = ', NFFT
c        end if
c      NFFT = MIN0(NFFT, NPFFT)
cC      WRITE(NOFILE,*)'Will use nfft = ', nfft
c 
cc     NFFT then taken a power of 2:
c      TRA6 = DLOG(DFLOAT(NFFT)) / DLOG(2.D0)
c      TRA1 = DINT(TRA6)
c      IF(TRA1.NE.TRA6)TRA6 = TRA6 + 1
c      MFFT = DINT(TRA6)
cc     At least 8 Fourier points for correct call of library routine:
c      MFFT = MAX0(3, MFFT)
c      NFFT = 2 ** MFFT
cC      WRITE(NOFILE,*)'Will use closest power of 2, nfft = ', nfft
c
c--------------------------------------------------------------------------

      TEM3 = 1.D0 / DFLOAT(NFFT)
     
c     Only re-initialize fft routine if number of points changed:
c     to check: second argument of cfft2: ±1?
        IF(NFFTME .NE. NFFT)THEN
        CALL CFFT2(1, 1, NFFT, V, WORK, V)
c       Store number of Fourier points for comparison at next routine call:
        NFFTME = NFFT
        END IF

        IF(updsym)THEN
c       Use up-down plasma symmetry:
        NPOI = NFFT / 2 + 1
        ELSE
        NPOI = NFFT + 1
        END IF
       
      ELSE
C     ----
c     Case no fft is required (at magnetic axis, or if poloidally symmetric problem):
      NFFT = 1
      TEM3 = 1.d0
      NPOI = 1

      END IF
C     ------

      FFSTEP = twopi / DFLOAT(NFFT)
c     Poloidal increment in equilibrium tables:
      ITAINC = NPFFT / NFFT
CC    WRITE(NOFILE,*)'FFT: M=',MFFT,' N=',NFFT,' 2*PI/N=',FFSTEP
CC    write(nofile,*)'polfft: itainc=',itainc,' npfft=',npfft,' nfft=',nfft
 
c     -----------------------------------------------------------------------------
c     Build lists of relevant k// on current surface and compute dielectric tensor:
c     -----------------------------------------------------------------------------

      IF(circ)THEN
c     Magnetic angle is constant on flux surfaces
      SI = SITAB(INTAB)
      CO = COTAB(INTAB)
        IF(CYL)THEN
c       k// and keta are constant on flux surfaces
c       kpar in m**-1:
        KPAR = MAVE * SI * RHOINV + KPHI * CO
c       kper * rnorm:
        KPER = MAVE * YINV * CO - KPRN * SI
        CALL DSET(NPOI, KPAR, KPARV, 1)
        CALL DSET(NPOI, KPER, KPERV, 1)
        END IF
      END IF
 
      CALL DCOPY(NPOI, POLANG, ITAINC, ANGLES, 1)
      IPO = 1
      IF(GLOMAX)THEN
c     R.K.'s nonmaxwellian tensor is NOT in use

        IF(.NOT.POLSYM)THEN
          IF(circ)THEN
            DO IR = 1, NPOI
            r0or = r0orta(INTAB,IPO)
            KPARV(IR) = MAVE * SI * RHOINV + KPHI * CO * R0OR
            KPERV(IR) = MAVE * YINV * CO - KPRN * R0OR * SI
            IPO = IPO + ITAINC
            END DO
          ELSE
            DO IR = 1, NPOI
            r0or = r0orta(INTAB,IPO)
            s7 = 1.d0 / eqt(INTAB,IPO,7)
c           k// = mave*sinThn/Ntn + n*cosTh/R:
            KPARV(IR) = MAVE * eqt(INTAB,IPO,12) * s7
     ;                + KPHI * eqt(INTAB,IPO,14) * R0OR
c           kprn:
c           keta*rnorm = mave*cosTh/(y*Ntn) - rnorm * n*sinTh/R:
            KPERV(IR) = MAVE * eqt(INTAB,IPO,14) * s7 * yinv
     ;                - KPRN * eqt(INTAB,IPO,15) * R0OR
            IPO = IPO + ITAINC
            END DO
          END IF
        END IF
c     Compute dielectric tensor on set of poloidal points:
      CALL HOTTEV(NPOI, kparv, SIGSD, SIGOD, ONESPE, ISP, gdrtoo, onlyab)
 
      ELSE
c     Case GLOMAX=.f., i.e. R.K.'s nonmaxwellian tensor is in use
        IF(circ)THEN
          DO IR = 1, NPOI
          intabp = ipo
c          r0or = r0orta(INTAB,IPO)
          KPARV(IR) = MAVE * SI * RHOINV + KPHI * CO * R0OR
          KPERV(IR) = MAVE * YINV * CO - KPRN * R0OR * SI
c          phi = polang(ipo)
c          BMODUL = BMOTAB(INTAB,IPO)
          CALL HOTTEN(KPARV(IR), SIGSD(1,1,IR), SIGOD(1,1,IR), ONESPE, ISP
     ;  , gdrtoo, onlyab)
          IPO = IPO + ITAINC
          END DO
        ELSE
          DO IR = 1, NPOI
          intabp = ipo
c          r0or = r0orta(INTAB,IPO)
          s7 = 1.d0  / eqt(INTAB,IPO,7)
          KPARV(IR) = MAVE * eqt(INTAB,IPO,12) * s7
     ;              + KPHI * eqt(INTAB,IPO,14) * R0OR
          KPERV(IR) = MAVE * eqt(INTAB,IPO,14) * s7 * yinv
     ;              - KPRN * eqt(INTAB,IPO,15) * R0OR
c          phi = polang(ipo)
c          BMODUL = BMOTAB(INTAB,IPO)
          CALL HOTTEN(KPARV(IR), SIGSD(1,1,IR), SIGOD(1,1,IR), ONESPE, ISP
     ;  , gdrtoo, onlyab)
          IPO = IPO + ITAINC
          END DO
        END IF
      END IF

c     General geom.: must perform coord. transfo. 
c     FLR terms: to do for M23!
        IF(.not.circ .and. STYP(ISUBR).EQ.'M23')THEN
        IPO = 1
          DO IR = 1, NPOI
c         Stix's S:
          a = 0.5d0 * (sigsd(1,1,ir) + sigsd(2,2,ir))
c         Stix's D:
          b = 0.5d0 * (sigsd(1,1,ir) - sigsd(2,2,ir))
c         Stix's P:
          c = sigsd(3,3,ir)
          co = eqt(intab,ipo,14)
          si = eqt(intab,ipo,15)
          co2 = co * co
          si2 = si * si
c         sigsd becomes diel.tensor matrix in r,t,p coordinates:
          sigsd(1,1,ir) = a
          sigsd(1,2,ir) = ci * b * co
          sigsd(1,3,ir) = - ci * b * si
          sigsd(2,2,ir) = a * co2 + c * si2
          sigsd(2,3,ir) = (c - a) * si * co
          sigsd(3,3,ir) = a * si2 + c * co2
          sigsd(2,1,ir) = - sigsd(1,2,ir)
          sigsd(3,1,ir) = - sigsd(1,3,ir)
          sigsd(3,2,ir) =   sigsd(2,3,ir)
          IPO = IPO + ITAINC
          END DO
c         if(flrops(ireg))then
c         end if
        END IF

c     ---------------------------------------------------------
c     Compact storage of distinct tensor elements in vector V1:
c     ---------------------------------------------------------

      IPO = 1
        if(circ)then
c       ~~~~~~~~~~~~
        DO IR = 1, NPOI
        KPER = KPERV(IR)
        KPER2 = KPER ** 2
        
        CALL ZSET(NTRANS, CZERO, V1, 1)
        TEM = TEM3
        
        IF(SWITCH.EQ.1)THEN
c       Volume terms

        IF(MUL)TEM = TEM / r0orta(INTAB,IPO)
        TEMK = TEM * KPER
        TEMK2 = TEMK * KPER
        
          IF(GLOMAX)THEN
c         L,R,P:
          V1(1) = SIGSD(1,1,IR) * TEM
          V1(2) = SIGSD(2,2,IR) * TEM
          V1(3) = SIGSD(3,3,IR) * TEM
            IF(.NOT.COLDPL(IREG) .AND. FLROPS(IREG))THEN
c           FLR terms:
            V1(17) = - SIGOD(1,1,IR) * TEMK2
            V1(18) = - SIGOD(2,2,IR) * TEMK2
            V1(19) = - SIGOD(3,3,IR) * TEMK2
            V1(4) = - SIGOD(1,1,IR) * TEM
            V1(5) = - SIGOD(2,2,IR) * TEM
            V1(6) = - SIGOD(3,3,IR) * TEM
            V1(7) =   SIGOD(1,2,IR) * TEMK2
            V1(8) = - SIGOD(1,3,IR) * TEMK
            V1(9) = - SIGOD(1,2,IR) * TEM
            V1(10) = - SIGOD(1,2,IR) * TEMK
            V1(11) =   SIGOD(1,3,IR) * TEM
            V1(12) =   SIGOD(2,3,IR) * TEMK
            V1(13) =   SIGOD(2,3,IR) * TEM
            V1(14) =   SIGOD(1,1,IR) * TEMK
            V1(15) = - SIGOD(2,2,IR) * TEMK
            V1(16) = - SIGOD(3,3,IR) * TEMK
            END IF
          ELSE
c         Some RK-nonmaxwellian species
          V1(1) = SIGSD(1,1,IR) * TEM
          V1(2) = SIGSD(2,2,IR) * TEM
          V1(3) = SIGSD(3,3,IR) * TEM
          V1(4) = SIGSD(1,2,IR) * TEM
          V1(5) = SIGSD(1,3,IR) * TEM
          V1(6) = SIGSD(2,3,IR) * TEM
          END IF
 
        ELSE IF(SWITCH.EQ.2)THEN
c       Boundary terms
        TEMK = TEM * KPER

          IF(COLJUM)THEN
c         No FLR terms in jump conditions
            IF(GLOMAX)THEN
c           L,R,P:
            V1(1) = SIGSD(1,1,IR) * TEM
            V1(2) = SIGSD(2,2,IR) * TEM
            V1(3) = SIGSD(3,3,IR) * TEM
            ELSE
c           Some RK-nonmaxwellian species
            V1(1) = SIGSD(1,1,IR) * TEM
            V1(2) = SIGSD(2,2,IR) * TEM
            V1(3) = SIGSD(3,3,IR) * TEM
            V1(4) = SIGSD(1,2,IR) * TEM
            V1(5) = SIGSD(1,3,IR) * TEM
            V1(6) = SIGSD(2,3,IR) * TEM
            END IF
          ELSE
c         Some FLR terms in jump conditions
          V1(1) = - SIGOD(1,1,IR) * TEMK
          V1(2) =   SIGOD(1,1,IR) * TEM
          V1(3) = - SIGOD(1,2,IR) * TEMK
          V1(4) =   SIGOD(1,2,IR) * TEM
          V1(5) =   SIGOD(2,2,IR) * TEMK
          V1(6) =   SIGOD(2,2,IR) * TEM
          V1(7) =   SIGOD(1,3,IR) * TEM
          V1(8) =   SIGOD(2,3,IR) * TEM
          V1(9) =   SIGOD(3,3,IR) * TEM
          END IF
        END IF
 
c         Contrib. of current poloidal point stored in table V:
          DO J = 1, NTR(IKIN)
          J2 = INDEX(J,IKIN)
          V(IR,J2) = V1(J2)
          END DO

        IPO = IPO + ITAINC
        END DO

        else
c       ~~~~
c       D-shaped and general equilibria
        DO IR = 1, NPOI
        KPER = KPERV(IR)
        KPER2 = KPER ** 2
        
        CALL ZSET(NTRANS, CZERO, V1, 1)
        
        IF(SWITCH.EQ.1)THEN
c       Volume terms

c       Normalized Jacobian / nfft:
        TEM = TEM3 * eqt(intab,ipo,9) / r0orta(INTAB,IPO)
        TEMK = TEM * KPER
        TEMK2 = TEMK * KPER
        
          IF(GLOMAX)THEN
c         L,R,P for HEC element, tensor in r,t,p coord.for M23 element:
          V1(1) = SIGSD(1,1,IR) * TEM
          V1(2) = SIGSD(2,2,IR) * TEM
          V1(3) = SIGSD(3,3,IR) * TEM
            if(eletyp .eq. 'M23')then
            V1(4) = SIGSD(1,2,IR) * TEM
            V1(5) = SIGSD(1,3,IR) * TEM
            V1(6) = SIGSD(2,3,IR) * TEM
            end if
c            IF(.NOT.COLDPL(IREG) .AND. FLROPS(IREG))THEN
c            to do!
c            END IF
          ELSE
c         Some RK-nonmaxwellian contributions
          V1(1) = SIGSD(1,1,IR) * TEM
          V1(2) = SIGSD(2,2,IR) * TEM
          V1(3) = SIGSD(3,3,IR) * TEM
          V1(4) = SIGSD(1,2,IR) * TEM
          V1(5) = SIGSD(1,3,IR) * TEM
          V1(6) = SIGSD(2,3,IR) * TEM
          END IF
 
        ELSE IF(SWITCH.EQ.2)THEN
c       Boundary terms
        TEM = TEM3

          IF(COLJUM)THEN
            IF(GLOMAX)THEN
c           L,R,P for HEC element, tensor in r,t,p coord.for M23 element:
            V1(1) = SIGSD(1,1,IR) * TEM
            V1(2) = SIGSD(2,2,IR) * TEM
            V1(3) = SIGSD(3,3,IR) * TEM
              if(eletyp .eq. 'M23')then
              V1(4) = SIGSD(1,2,IR) * TEM
              V1(5) = SIGSD(1,3,IR) * TEM
              V1(6) = SIGSD(2,3,IR) * TEM
              end if
            ELSE
c           Some RK-nonmaxwellian contributions
            V1(1) = SIGSD(1,1,IR) * TEM
            V1(2) = SIGSD(2,2,IR) * TEM
            V1(3) = SIGSD(3,3,IR) * TEM
            V1(4) = SIGSD(1,2,IR) * TEM
            V1(5) = SIGSD(1,3,IR) * TEM
            V1(6) = SIGSD(2,3,IR) * TEM
            END IF
          ELSE
c         to do!
          END IF
        END IF
 
          DO J = 1, NTR(IKIN)
          J2 = INDEX(J,IKIN)
          V(IR,J2) = V1(J2)
          END DO

        IPO = IPO + ITAINC
        END DO
        end if
c       ~~~~~~

      IF(TRAFOL)THEN
c     CALL CRAY FAST FOURIER TRANSFORM, N=2**M CASE:
C     ----------------------------------------------
      IF(updsym .and. NFFT.GE.4)THEN
c     Uses up-down symmetry:
        DO I = 1, NTR(IKIN)
        I2 = INDEX(I,IKIN)
          DO IR = 2, NFFT / 2
          V(NFFT-IR+2,I2) = V(IR,I2)
          END DO
        END DO
      END IF
 
c to check: second argument of cfft2: ±1?
        DO I = 1, NTR(IKIN)
        I2 = INDEX(I,IKIN)
        CALL CFFT2(0, 1, NFFT, V(1,I2), WORK, V(1,I2))
        END DO
      END IF

 1000 CONTINUE
 
CC    WRITE(NOFILE,*)   'EXIT POLFFT'
 
      RETURN
      END
