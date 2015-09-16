      subroutine polfft(v, ldv, onespe, isp, gdrtoo, switch, trafo, ipath
     ;, onlyab)

c 3/4/2004: suppressed one argument before last (nfft), obsolete.

ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)
      implicit none

      logical onespe, trafo, gdrtoo, onlyab

      integer ldv, isp, switch, ipath, nfft

      complex*16 v(ldv,19)

      include 'pardim.copy'
      include 'comsub.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comant.copy'
      include 'commod.copy'
      include 'compla.copy'
      include 'comin2.copy'
      include 'comfic.copy'
      include 'comswe.copy'
      include 'comphy.copy'
      include 'comfft.copy'

c     29/10/2000: added argument ONLYAB: if .true., only use damping terms
c     in dielectric tensor (imag.part of plasma dispersion function) and
c     ignore reactive terms.
c
c     This version with loop over average M
c
c     Susceptibility tensor: - for species isp (onespe=.t.)
c                            - all species     (onespe=.f., dummy isp)
c     Fast poloidal Fourier transform : call to IMSL double complex routine
c     (formerly call to cray routine cfft2)
c
c     This version uses vectorized plasma dispersion function.
c     Output: v(ldv,19) 
c     (with spectra valid - but not accurate - up to nfft/2 th harmonic,
c      where nfft=npfft for resonant terms and nftnr for nonresonant ones.)
c
c     switch: = 1: computes volume terms
c             = 2: computes boundary terms
c     trafo:
c
c     ipath: input and output (=2 on output). 
c            When =1: prepares index lists; when not 1, bypasses index list building.

      logical mul, trafol, vssh

      character*3 eletyp

      integer
     ;  nfc, index(19,2), ntr(2), isp1, isp2
     ;, ipo, ikin, npoi, ir, is, i, j, i2, j2, lres, ifc, ntrans
c     ;, nfftme

c     valid for up-down asymmetric equilibria:
      parameter(nfc = npfft+1)

      double precision 
     ;  ffstep, angles(nfc)
     ;, kper, kper2, kpar
     ;, kparv(nfc), kperv(nfc)
     ;, tem, tem3, temk, temk2, co2, si2, s7

      complex*16 sigsd(3,3,nfc), sigod(3,3,nfc), v1(19)
     ;, a, b, c

c      data nfftme/0/
c      save nfftme
 
cc    write(nofile,*)'enter polfft'
      eletyp = styp(isubr)
	  
c     mav2 is passed through commod:
      mave = dfloat(mav2) * 0.5
 
c     IPATH not 1 to short-circuit index list building:
      if(ipath .ne. 1)goto 8

c     Don't FFT at magnetic axis:
      trafol = trafo .and. y.gt.1.d-40
c     Following is true if a torus, not at magn. axis and TRAFO=.t.:
      mul = trafol .and. .not.cyl
 
      rho = abscis(intab)
      yinv = abscni(intab)
      rhoinv = yinv * rnori
        if(onespe)then
        isp1 = isp
        isp2 = isp
        else
        isp1 = 1
        isp2 = nspec
        end if

c     Define lists of coefficients to poloidal-Fourier-transform:

      if(switch.eq.1)then
c     +++++++++++++++++++
c     This case deals with volume terms

c       Label resonant and nonresonant dielectric tensor elements:      
        if(glomax)then
c       --------------
c       case Raymond Koch's nonmaxwellian tensor not in use
        
          if(coldpl(ireg) .or. .not.flrops(ireg))then
c         Zeroth order in Larmor: only ±1 cyclotron resonances.

            if(circ .or. eletyp.eq.'HEC')then
c           3 terms in tensor:
c           - diagonal(L,R,P) in +-// coordinates
c           - or elements S,±iD,P
c           Number of Fourier transforms:
            ntrans = 3
            ntr(1) = 0
            ntr(2) = 0
              do lres = -1, 1
c             'RES' was computed in cyrano.f (lres = -2 to 2). 
c             (It does not detect Landau resonance, but the corresponding terms are smooth poloidally.)
c               Cyclotron resonance: index: value 2, 3, 1 for resp. -, //, +
c               associated with lres=-1,0,1
	          if(lres.eq.-1)then
	          ifc = 2
	          else if(lres.eq.0)then
	          ifc = 3
	          else
	          ifc = 1
	          end if
                if(res(lres))then
c               Resonant terms:
                ntr(1) = ntr(1) + 1
c               Cyclotron resonance: index: value 2, 3, 1 for resp. -, //, +
c               associated with lres=-1,0,1
                index(ntr(1),1) = ifc
cPL  too clever to be readable:    index(ntr(1),1) = 3 - lres*(3*lres+1) / 2
                else
c               Nonresonant terms:
                ntr(2) = ntr(2) + 1
                index(ntr(2),2) = ifc
cPL                index(ntr(2),2) = 3 - lres*(3*lres+1) / 2
                end if
              end do
            else  ! Noncircular or M23 element:
c           6 independent terms in tensor in rho, theta, phi coordinates:
            ntrans = 6
            ntr(1) = 0
            ntr(2) = 0
c           Quick, not optimized: terms assumed all resonant or nonresonant
              if(res(-1) .or. res(0) .or. res(1))then
              ntr(1) = 6
                do i = 1, 6
                index(i,1) = i
                end do
              else
              ntr(2) = 6
                do i = 1, 6
                index(i,2) = i
                end do
              end if
            end if

          else
c         Dielectric tensor to second order in kperp*Larmor inclusive:
c         NB: here resonances are treated as mutually exclusive
cPL2June04: this limitation must be removed for e.g. fundamental H together with 2nd harmonic D and other important scenarios!
c         NB following excludes electron cyclotron resonances and anomalous cyclotron damping on p < 0
          ntrans = 19
c         The 19 coefficients are:
c@
            if(res(1))then
c           Fundamental ion cyclotron resonance is present
            ntr(1) = 12
            index(1,1) = 1
            index(2,1) = 3
            index(3,1) = 4
            index(4,1) = 6
            index(5,1) = 7
            index(6,1) = 8
            index(7,1) = 9
            index(8,1) = 10
            index(9,1) = 11
            index(10,1) = 14
            index(11,1) = 16
            index(12,1) = 17
            ntr(2) = 7
            index(1,2) = 2
            index(2,2) = 5
            index(3,2) = 12
            index(4,2) = 13
            index(5,2) = 15
            index(6,2) = 18
            index(7,2) = 19
          
            else if(res(2))then
c           Second harmonic ion cyclotron resonance is present
            ntr(1) = 3
            index(1,1) = 4
            index(2,1) = 14
            index(3,1) = 17
            ntr(2) = 16
            index(1,2) = 1
            index(2,2) = 2
            index(3,2) = 3
            index(4,2) = 5
            index(5,2) = 6
            index(6,2) = 7
            index(7,2) = 8
            index(8,2) = 9
            index(9,2) = 10
            index(10,2) = 11
            index(11,2) = 12
            index(12,2) = 13
            index(13,2) = 15
            index(14,2) = 16
            index(15,2) = 18
            index(16,2) = 19

            else
c           No fundamental or second harmonic ion cyclotron resonance is present
c           All coefficients are assumed poloidally smooth
c           NB excludes electron cyclotron resonances and anomalous cyclotron damping on p <0
            ntr(1) = 0
            ntr(2) = 19
              do i = 1, ntr(2)
              index(i,2) = i
              end do
            end if
          end if
       
        else
C       ----
c       Case GLOMAX=.f.: some species use R.Koch's nonmaxwellian tensor
        ntrans = 6
          if(res(1))then
c         Fundamental ion cyclotron resonance is present
          ntr(1) = 1
          index(1,1) = 1
          ntr(2) = 5
          index(1,2) = 2
          index(2,2) = 3
          index(3,2) = 4
          index(4,2) = 5
          index(5,2) = 6
           
          else
c         No fundamental ion cyclotron resonance is present
c         All coefficients are assumed poloidally smooth - check why!
          ntr(1) = 0
          ntr(2) = 6
          index(1,2) = 1
          index(2,2) = 2
          index(3,2) = 3
          index(4,2) = 4
          index(5,2) = 5
          index(6,2) = 6
          end if
        end if
c       ------

      else if(switch.eq.2)then
c     ++++++++++++++++++++++++
c     This case deals with boundary terms

      if(coljum)then
c     'Cold-like' jump condition:
 
        if(glomax)then
c       --------------
c       No call to R.K.'s nonmaxwellian routines
cPL          if(circ .or. eletyp.eq.'HEC')then
cPL     Assume between plasma and vacuum and use rtp coord. for jump equations
          ntrans = 3
            if(res(1))then
            ntr(1) = 1
c           ++ term is resonant
            index(1,1) = 1
            ntr(2) = 2
c           ... and -- and // terms are smooth
            index(1,2) = 2
            index(2,2) = 3
            else if(res(-1))then
            ntr(1) = 1
c           -- term is resonant (anomalous case or ECRH)
            index(1,1) = 2
            ntr(2) = 2
c           ... and ++ and // terms are smooth
            index(1,2) = 1
            index(2,2) = 3
            else
c           all terms are smooth
            ntr(1) = 0
            ntr(2) = 3
            index(1,2) = 1
            index(2,2) = 2
            index(3,2) = 3
            end if
c          else
cc         6 independent terms in tensor in rho, theta, phi coordinates:
c          ntrans = 6
c          ntr(1) = 0
c          ntr(2) = 0
c            if(res(-1) .or. res(0) .or. res(1))then
c            ntr(1) = 6
c              do i = 1, 6
c              index(i,1) = i
c              end do
c            else
c            ntr(2) = 6
c              do i = 1, 6
c              index(i,2) = i
c              end do
c            end if
cPL          end if

        else
c       ----
c       At least a species requires R.K.'s nonmaxwellian routines
        ntrans = 6
          if(res(1))then
          ntr(1) = 1
          index(1,1) = 1
          ntr(2) = 5
          index(1,2) = 2
          index(2,2) = 3
          index(3,2) = 4
          index(4,2) = 5
          index(5,2) = 6

          else
          ntr(1) = 0
          ntr(2) = 6
          index(1,2) = 1
          index(2,2) = 2
          index(3,2) = 3
          index(4,2) = 4
          index(5,2) = 5
          index(6,2) = 6
          end if
        end if
c       ------

      else
c     There are flr terms in the jump conditions
      ntrans = 9
        if(res(1))then
c       Fundamental ion cyclotron resonance is present
        ntr(1) = 6
        index(1,1) = 1
        index(2,1) = 2
        index(3,1) = 3
        index(4,1) = 4
        index(5,1) = 7
        index(6,1) = 9
        ntr(2) = 3
        index(1,2) = 5
        index(2,2) = 6
        index(3,2) = 8

        else if(res(2))then
c       Second harmonic ion cyclotron resonance is present, but no fundamental
        ntr(1) = 2
        index(1,1) = 1
        index(2,1) = 2
        ntr(2) = 7
        index(1,2) = 3
        index(2,2) = 4
        index(3,2) = 5
        index(4,2) = 6
        index(5,2) = 7
        index(6,2) = 8
        index(7,2) = 9

        else
c       No fundamental or second harmonic ion cyclotron resonance is present
c       All coefficients are assumed poloidally smooth
c       NB excludes electron cyclotron resonances and anomalous cyclotron damping on p <0
        ntr(1) = 0
        ntr(2) = 9
          do i = 1, 9
          index(i,2) = i
          end do
        end if
      end if
      end if
 
      ipath  =  2
  8   continue

c     Initialize output table:
      call zset(ldv * ntrans, czero, v, 1)
 
c     Loop over resonant and non resonant terms:
      do 1000 ikin = 1, 2
ccc	if(ntr(ikin) .gt. 0)then

c     Choose number of poloidal points: 
      if(trafol)then
c     --------------

c     23/10/2000: simplified the logic:
	vssh = .false.
        if(nscree.gt.0)then
c       Faraday shields with vacuum inside: take into account
c       to find number of poloidal Fourier points.
c       VSSH=.t. if there is vacuum in one of the screens.
          do is = 1, nscree
	      if(vscrsh(is) .and. ireg.gt.irbscr(is))vssh = .true.
          end do
        end if

        if(ikin.eq.1 .or. (nscree.gt.0 .and. vssh))then
c       Use allowed maximum of Fourier points NPFFT for the resonant terms or if
c       the volume inside a Faraday shield generates poloidal inhomogeneities:
c       (NPFFT must be a power of 2 and >=8)
        nfft = npfft
        else
c       Number of poloidal points for fft of smooth terms, NFTNR, is user-supplied:
c       Must be a power of 2!
        nfft = nftnr
        end if
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

      tem3 = 1.d0 / dfloat(nfft)
     
c        if(nfftme .ne. nfft)then
cc       Store number of Fourier points for comparison at next routine call:
c        nfftme = nfft
c        end if

        if(updsym)then
c       Use up-down plasma symmetry:
        npoi = nfft / 2 + 1
        else
        npoi = nfft + 1
        end if
       
      else
c     ----
c     Case no fft is required (at magnetic axis, or if poloidally symmetric problem):
      nfft = 1
      tem3 = 1.d0
      npoi = 1

      end if
c     ------

      ffstep = twopi / dfloat(nfft)
c     Poloidal increment in equilibrium tables:
      itainc = npfft / nfft
CC    WRITE(NOFILE,*)'FFT: M=',MFFT,' N=',NFFT,' 2*PI/N=',FFSTEP
CC    write(nofile,*)'polfft: itainc=',itainc,' npfft=',npfft,' nfft=',nfft
 
c     -----------------------------------------------------------------------------
c     Build lists of relevant k// on current surface and compute dielectric tensor:
c     -----------------------------------------------------------------------------

      if(circ)then
c     Magnetic angle is constant on flux surfaces
      si = sitab(intab)
      co = cotab(intab)
        if(cyl)then
c       k// and keta are constant on flux surfaces
c       kpar in m**-1:
        kpar = mave * si * rhoinv + kphi * co
c       kper * rnorm:
        kper = mave * yinv * co - kprn * si
        call dset(npoi, kpar, kparv, 1)
        call dset(npoi, kper, kperv, 1)
        end if
      end if
 
      call dcopy(npoi, polang, itainc, angles, 1)
      ipo = 1
      if(glomax)then
c     R.K.'s nonmaxwellian tensor is NOT in use

        if(.not.polsym)then
          if(circ)then
            do ir = 1, npoi
            r0or = r0orta(intab,ipo)
            kparv(ir) = mave * si * rhoinv + kphi * co * r0or
            kperv(ir) = mave * yinv * co - kprn * r0or * si
            ipo = ipo + itainc
            end do
          else
            do ir = 1, npoi
            r0or = r0orta(intab,ipo)
            s7 = 1.d0 / eqt(intab,ipo,7)
c           k// = mave*sinthn/ntn + n*costh/r:
            kparv(ir) = mave * eqt(intab,ipo,12) * s7
     ;                + kphi * eqt(intab,ipo,14) * r0or
c           kprn:
c           keta*rnorm = mave*costh/(y*ntn) - rnorm * n*sinth/r:
            kperv(ir) = mave * eqt(intab,ipo,14) * s7 * yinv
     ;                - kprn * eqt(intab,ipo,15) * r0or
            ipo = ipo + itainc
            end do
          end if
        end if
c     Compute dielectric tensor on set of poloidal points:
      call hottev(npoi, kparv, sigsd, sigod, onespe, isp, gdrtoo, onlyab)
 
      else
c     Case GLOMAX=.f., i.e. R.K.'s nonmaxwellian tensor is in use
        if(circ)then
          do ir = 1, npoi
          intabp = ipo
c          r0or = r0orta(intab,ipo)
          kparv(ir) = mave * si * rhoinv + kphi * co * r0or
          kperv(ir) = mave * yinv * co - kprn * r0or * si
c          phi = polang(ipo)
c          bmodul = bmotab(intab,ipo)
          call hotten(kparv(ir), sigsd(1,1,ir), sigod(1,1,ir), onespe, isp
     ;  , gdrtoo, onlyab)
          ipo = ipo + itainc
          end do
        else
          do ir = 1, npoi
          intabp = ipo
c          r0or = r0orta(intab,ipo)
          s7 = 1.d0  / eqt(intab,ipo,7)
          kparv(ir) = mave * eqt(intab,ipo,12) * s7
     ;              + kphi * eqt(intab,ipo,14) * r0or
          kperv(ir) = mave * eqt(intab,ipo,14) * s7 * yinv
     ;              - kprn * eqt(intab,ipo,15) * r0or
c          phi = polang(ipo)
c          bmodul = bmotab(intab,ipo)
          call hotten(kparv(ir), sigsd(1,1,ir), sigod(1,1,ir), onespe, isp
     ;  , gdrtoo, onlyab)
          ipo = ipo + itainc
          end do
        end if
      end if

        if(.not.circ .and. styp(isubr).eq.'M23')then
c       General geom.: must perform coord. transfo. for M23
c       FLR terms: to do for M23!
        ipo = 1
          do ir = 1, npoi
c         Stix's S-1=(R+L)/2-1:
          a = 0.5d0 * (sigsd(1,1,ir) + sigsd(2,2,ir))
c         Stix's D=(R-L)/2:
          b = 0.5d0 * (sigsd(2,2,ir) - sigsd(1,1,ir))
cPL          b = 0.5d0 * (sigsd(1,1,ir) - sigsd(2,2,ir))
cPL       no net change: just made D consistent with Stix!
c         Stix's P-1:
          c = sigsd(3,3,ir)
          co = eqt(intab,ipo,14)
          si = eqt(intab,ipo,15)
          co2 = co * co
          si2 = si * si
c         Caution: sigsd becomes diel.tensor matrix in r,t,p coordinates:
c           (    S-1,         -iD co,                    iD  si )
c           (  iD co,  co2(S-1) + si2(P-1),         si co (P-S) )
c           ( -iD si,          si co (P-S),  si2(S-1) + co2(P-1))

          sigsd(1,1,ir) = a
          sigsd(1,2,ir) = - ci * b * co
          sigsd(1,3,ir) =   ci * b * si
cPL          sigsd(1,2,ir) =  ci * b * co
cPL          sigsd(1,3,ir) = - ci * b * si
cPL       no net change: just made D consistent with Stix!
          sigsd(2,2,ir) = a * co2 + c * si2
          sigsd(2,3,ir) = (c - a) * si * co
          sigsd(3,3,ir) = a * si2 + c * co2
          sigsd(2,1,ir) = - sigsd(1,2,ir)
          sigsd(3,1,ir) = - sigsd(1,3,ir)
          sigsd(3,2,ir) =   sigsd(2,3,ir)
          ipo = ipo + itainc
          end do
c         if(flrops(ireg))then
c         end if
        end if

c     ---------------------------------------------------------
c     Compact storage of distinct tensor elements in vector V1:
c     ---------------------------------------------------------

      ipo = 1
        if(circ)then
c       ~~~~~~~~~~~~
        do ir = 1, npoi
        kper = kperv(ir)
        kper2 = kper ** 2
        
        call zset(ntrans, czero, v1, 1)
        tem = tem3
        
        if(switch.eq.1)then
c       Volume terms

        if(mul)tem = tem / r0orta(intab,ipo)
        temk = tem * kper
        temk2 = temk * kper
        
          if(glomax)then
c         l,r,p:
          v1(1) = sigsd(1,1,ir) * tem
          v1(2) = sigsd(2,2,ir) * tem
          v1(3) = sigsd(3,3,ir) * tem
            if(.not.coldpl(ireg) .and. flrops(ireg))then
c           FLR terms:
            v1(17) = - sigod(1,1,ir) * temk2
            v1(18) = - sigod(2,2,ir) * temk2
            v1(19) = - sigod(3,3,ir) * temk2
            v1(4) = - sigod(1,1,ir) * tem
            v1(5) = - sigod(2,2,ir) * tem
            v1(6) = - sigod(3,3,ir) * tem
            v1(7) =   sigod(1,2,ir) * temk2
            v1(8) = - sigod(1,3,ir) * temk
            v1(9) = - sigod(1,2,ir) * tem
            v1(10) = - sigod(1,2,ir) * temk
            v1(11) =   sigod(1,3,ir) * tem
            v1(12) =   sigod(2,3,ir) * temk
            v1(13) =   sigod(2,3,ir) * tem
            v1(14) =   sigod(1,1,ir) * temk
            v1(15) = - sigod(2,2,ir) * temk
            v1(16) = - sigod(3,3,ir) * temk
            end if
          else
c         Some RK-nonmaxwellian species
          v1(1) = sigsd(1,1,ir) * tem
          v1(2) = sigsd(2,2,ir) * tem
          v1(3) = sigsd(3,3,ir) * tem
          v1(4) = sigsd(1,2,ir) * tem
          v1(5) = sigsd(1,3,ir) * tem
          v1(6) = sigsd(2,3,ir) * tem
          end if
 
        else if(switch.eq.2)then
c       Boundary terms
        temk = tem * kper

          if(coljum)then
c         No FLR terms in jump conditions
            if(glomax)then
c           L-1,R-1,P-1:
            v1(1) = sigsd(1,1,ir) * tem
            v1(2) = sigsd(2,2,ir) * tem
            v1(3) = sigsd(3,3,ir) * tem
            else
c           Some RK-nonmaxwellian species
            v1(1) = sigsd(1,1,ir) * tem
            v1(2) = sigsd(2,2,ir) * tem
            v1(3) = sigsd(3,3,ir) * tem
            v1(4) = sigsd(1,2,ir) * tem
            v1(5) = sigsd(1,3,ir) * tem
            v1(6) = sigsd(2,3,ir) * tem
            end if
          else
c         Some FLR terms in jump conditions
          v1(1) = - sigod(1,1,ir) * temk
          v1(2) =   sigod(1,1,ir) * tem
          v1(3) = - sigod(1,2,ir) * temk
          v1(4) =   sigod(1,2,ir) * tem
          v1(5) =   sigod(2,2,ir) * temk
          v1(6) =   sigod(2,2,ir) * tem
          v1(7) =   sigod(1,3,ir) * tem
          v1(8) =   sigod(2,3,ir) * tem
          v1(9) =   sigod(3,3,ir) * tem
          end if
        end if
 
c         Contrib. of current poloidal point stored in table V:
          do j = 1, ntr(ikin)
          j2 = index(j,ikin)
          v(ir,j2) = v1(j2)
          end do

        ipo = ipo + itainc
        end do

        else
c       ~~~~
c       D-shaped and general equilibria
        do ir = 1, npoi
        kper = kperv(ir)
        kper2 = kper * kper
        
        call zset(ntrans, czero, v1, 1)
        
        if(switch.eq.1)then
c       Volume terms

c       Normalized Jacobian / nfft:
        tem = tem3 * eqt(intab,ipo,9) / r0orta(intab,ipo)
        temk = tem * kper
        temk2 = temk * kper
        
          if(glomax)then
c         L-1,R-1,P-1 for HEC element
c@ CHECK: tensor in r,t,p coord.for M23 element: must have proper ldif in genera
          v1(1) = sigsd(1,1,ir) * tem
          v1(2) = sigsd(2,2,ir) * tem
          v1(3) = sigsd(3,3,ir) * tem
            if(eletyp .eq. 'M23')then
            v1(4) = sigsd(1,2,ir) * tem
            v1(5) = sigsd(1,3,ir) * tem
            v1(6) = sigsd(2,3,ir) * tem
            end if
c            IF(.NOT.COLDPL(IREG) .AND. FLROPS(IREG))THEN
c            to do!
c            END IF
          else
c         Some RK-nonmaxwellian contributions
          v1(1) = sigsd(1,1,ir) * tem
          v1(2) = sigsd(2,2,ir) * tem
          v1(3) = sigsd(3,3,ir) * tem
          v1(4) = sigsd(1,2,ir) * tem
          v1(5) = sigsd(1,3,ir) * tem
          v1(6) = sigsd(2,3,ir) * tem
          end if
 
        else if(switch.eq.2)then
c       Boundary terms
        tem = tem3
c @ any jacobian here? NO
          if(coljum)then
            if(glomax)then
c              if(eletyp .eq. 'HEC')then
cc             L-1,R-1,P-1 for HEC element:
              v1(1) = sigsd(1,1,ir) * tem
              v1(2) = sigsd(2,2,ir) * tem
              v1(3) = sigsd(3,3,ir) * tem
c              else if(eletyp .eq. 'M23')then
c See later    Tensor with rows in r,t,p coord.for M23 element
c              v1(4) = sigsd(1,2,ir) * tem
c              v1(5) = sigsd(1,3,ir) * tem
c              v1(6) = sigsd(2,3,ir) * tem
c              end if
            else
c           Some RK-nonmaxwellian contributions
            v1(1) = sigsd(1,1,ir) * tem
            v1(2) = sigsd(2,2,ir) * tem
            v1(3) = sigsd(3,3,ir) * tem
            v1(4) = sigsd(1,2,ir) * tem
            v1(5) = sigsd(1,3,ir) * tem
            v1(6) = sigsd(2,3,ir) * tem
            end if
          else
c         FLR jumps to do in noncircular!
          end if
        end if
 
          do j = 1, ntr(ikin)
          j2 = index(j,ikin)
          v(ir,j2) = v1(j2)
          end do

        ipo = ipo + itainc
        end do
        end if
c       ~~~~~~

      if(trafol)then
c     Call fft, n=2**m case:
c     ---------------------
        if(updsym .and. nfft.ge.4)then
c       Uses up-down symmetry:
          do i = 1, ntr(ikin)
          i2 = index(i,ikin)
            do ir = 2, nfft / 2
            v(nfft-ir+2,i2) = v(ir,i2)
            end do
          end do
        end if
 
c to check: second argument of cfft2: ±1?
	  if(nfft.eq.nftnr)then
c       Nonresonant terms:
          do i = 1, ntr(ikin)
          i2 = index(i,ikin)
cccccccccccccccc ERNESTO cccccccccccccccccccccccccccc
cERN      call cfft2(0, 1, nftnr, v(1,i2), worknr, v(1,i2))
	    call df2tcb(nftnr, v(1,i2), v(1,i2), worknr2, cpy) ! -- IMSL --
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          end do
c	  2003 fix: store harmonics as if always NPFFT, and pad with zeros:
	    do i = 1, ntr(ikin)
	    i2 = index(i, ikin)
	      do j = 1, nfft/2
	      v(npfft-j+1,i2) = v(nfft-j+1,i2)
	      end do
	      do j = nfft/2+1, npfft-nfft/2
	      v(j,i2) = czero
	      end do
	    end do
	  elseif(nfft.eq.npfft)then
          do i = 1, ntr(ikin)
          i2 = index(i,ikin)
ccccccccccccccccccc ERNESTO cccccccccccccccccccccccccccc
cERN          call cfft2(0, 1, npfft, v(1,i2), cworkp, v(1,i2))
	        call df2tcb(npfft, v(1,i2), v(1,i2), cwork2, cpy) ! -- IMSL --
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          end do
	  else
	  write(6,*)'Error in POLFFT: nfft is neither npfft nor nftnr!'
	  write(6,*)'nfft=',nfft
	  stop
	  end if
ccccccccccccccccccc ERNESTO cccccccccccccccccccccccccccc
c	print *, v(1:10,1)
c	print *, 'STOP at POLFFT'
c	stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
      end if

ccc      end if  ! ntr(ikin) > 0
 1000 continue  ! ikin = 1, 2
 
cc    write(nofile,*)   'exit polfft'
 
      return
      end
