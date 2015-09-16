      subroutine outrff

ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)
      implicit none

c     Builds output tables, to be plotted elsewhere.
c     New, shortened version (1998) of OUTPUT, free of graphics software.
 
      include 'pardim.copy'
      include 'dynou2.copy'
      include 'compow.copy'
      include 'comgeo.copy'
      include 'comant.copy'
      include 'commod.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comswe.copy'
      include 'com3di.copy'
      include 'comfou.copy'
      include 'comro2.copy'
      include 'complo.copy'
      include 'complp.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comfin.copy'
      include 'comfic.copy'
      include 'comequ.copy'
      include 'comin2.copy'
      include 'comphy.copy'
      include 'comfft.copy'
       
      logical inl, inr
 
      character*3 eletyp
      character*9 charfi(2), charco(6), charpa(3)
	character(100) GGs 
      integer 
     ;  i, i2, ip, it, j, aqp(3), aux
     ;, ipol, ii, jj, ii2
     ;, is1, is2
     ;, idllo, icolo, ibulo, icpblo, ia, ja, ja1, ja2, jc, im, imd, idecax
     ;, idecar, idecam, jae, jae1, jae2
     ;, ireg1, ifirse, ie1, nmodc, m1, m2, mrg, mabs, mrd
     ;, ima, iblk, ireq, j1, j2
     ;, lbl, ic, istatu, ie2, nintma, jd, j22, j22i, img
     ;, j22m1i, j22p1, j22p1i, minfc, msupc, k1, k2, kp, ki, impk, nblk, ist
     ;, istg, mra
     ;, j22m1, nrhss
     ;, nwr
     ;, ibl
     ;, jm
     ;
     ;, nsum, i1st, ceilq
     ;, idamin, idmin, idmax, isrchfge

      double precision 
     ;  rle, rlei, absstp, absmax, wmin
     ;, fnin1i
     ;, cabs2, rveno, dsinc
     ;, bafplo(21,0:maxinp+1), bafpln(maxbaf,0:2,ntyp,0:maxinp+1), tr
     ;, bafpls(maxbaf,0:1,ntyp,0:maxinp+1)
     ;, xtra, fac4, st2t2i
     ;, kmax, kstep
     ;, trapeze, gauss
     ;, second
     ;, tdelta, talpha, tbeta

      complex*16
     ;  toi, toj
     ;, poyjum(0:maxreg)
     ;, c1, c2, ctra2(maxant)
     ;, zdotc, zdotu, cdeid
     ;, tra1, sum1, sum2, sum3, torfac, polfa2, fac, cbl, cbli, tinpcp
     ;, tofacv(maxant)
     ;, cbl2, cbl2i, ztra
     ;, tofac, polfac(maxpom)
     ;, zimp2(maxant,maxant,maxtom), fac3
     ;, eeta
     ;, corr
c     ;, erho

      complex*16 foutra(-npft2:npft2,4)
	double precision t, ta(npfft,4)
 
      external 
     ;  nsum, i1st, ceilq, trapeze, gauss, rveno
     ;, zdotc, zdotu, zset, dset, dsum, mucrvz, cdeid, tofac, dsinc
     ;, idamin, idmin, idmax
 
      integer 
     ;  trncp1, trncp2, trncp5
 
      external second
c
      data charfi/'Electric ','Magnetic '/
     ;    ,charco/'radial ','poloidal ','toroidal ',
     ;            'L.H.(+) ','R.H.(-) ','parallel '/
     ;    ,charpa/'real','imag','modulus'/
 
      cabs2(tra1) = dreal(tra1 * dconjg(tra1))
 
c**********************************************************************
 
      write(nofile,*)'Enter OUTRFF'
      
      call zset((ndof-1)*maxpom*maxreg*maxexc*maxtom, czero, xrtpb, 1)
      st2t2i = 0.5d0 * sqrt2i

c     Field tor. modes RTP components at region boundaries:
c     NB: with 'HEC', only works in circular!
c     XRTPB used below in impedance calc.; must be written if 'HEC' subregion
c     touches an antenna!
      do imoto = 1, ntotor
      do ireg = 1, nreg
c     ====================
      y = rx0(ireg)
      isubr = nsum(ireg, ns)
      eletyp = styp(isubr)
      iel = ilael(ireg)
      intab = iel * (ngauss+1) + ireg
      if(circ .and. eletyp .eq. 'HEC')call loctr2
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo

      call rdsol(imoto, iel+1+2*(ireg-1), bel, maxbll, lbl, nrhss)

      rle = fl(iel)
c     Fill xrtpb for 'effective' number of antennae used in solving system:
      do ja = 1, neffan + nscree * nmoscr
        if(eletyp.eq.'M23')then
          do j2 = 1, ndof-1
            do im = 1, nmode(iel)
            idecax = (im - 1) * icolo
            xrtpb(j2,ireg,im,ja,imoto) = bel(j2+idecax, ja)
            end do
          end do
        else if(eletyp.eq.'HEC')then
          if(circ)then
            do jc = 1, ndof
              do ic = 1, ndof - 1
                do im = 1, nmode(iel)
                idecax = (im - 1) * icolo
                xrtpb(ic,ireg,im,ja,imoto) = xrtpb(ic,ireg,im,ja,imoto) + bc2(ic,jc) * bel(jc+idecax, ja)
                end do
              end do
            end do
          else
c@ 16/12/03 at work...
c         Compute poloidal FT of si, co and their radial derivative:
c         RCFFT2 normalisation:
c	    t = 1.d0 / (2.d0 * npfft)
c         IMSL normalisation:
	    t = 1.d0 / npfft
	      do i = 1, npfft
	      co = eqt(intab,i,14)
	      si = eqt(intab,i,15)
	      si1 = eqt(intab,i,16) * rnorm
	      co1 = - si1 * si / co
	      ta(i,1) = t * co
	      ta(i,2) = t * si
	      ta(i,3) = t * co1
	      ta(i,4) = t * si1
            end do
	      do i = 1, 4
cERN	      call rcfft2(0, -1, npfft, ta(1,i), workm, foutra(0,i))
	      call df2trf(npfft, ta(1,i), ta(1,i), work2) ! -- IMSL --
c           Here we want a Fourier transform with integral of exp(-i k theta)!
	      foutra(0,i) = dcmplx(ta(1,i), 0.d0)
	        do j = 1, npft2-1
	        foutra(j,i) = dcmplx(ta(2*j,i), ta(2*j+1,i))
              end do
	      foutra(npft2,i) = dcmplx(ta(npfft,i), 0.d0)
	      end do
cPL 7/4/2004: this loop was missing:
		do i = 1, 4
	        do j = 1, npft2
	        foutra(-j,i) = dconjg(foutra(j,i))
              end do
	      end do
cPL
c           Convolution with +, -, // components to get rho, theta, phi:
            do im = 1, nmode(iel)
	        do jm = 1, nmode(iel)
cPL 7/4/2004: this seems the correct sign:
	        k = jm - im
cPL	        k = im - jm
              idecax = (jm - 1) * icolo
              call bc2spectrum(foutra(k,1), foutra(k,2), foutra(k,3), foutra(k,4), k) 
                do jc = 1, ndof
                  do ic = 1, ndof - 1
                  xrtpb(ic,ireg,im,ja,imoto) = xrtpb(ic,ireg,im,ja,imoto) + bc2(ic,jc) * bel(jc+idecax, ja)
                  end do
                end do
              end do
            end do
          end if
        end if
      end do
      end do
      end do
c     ======

      if(wrisol)then
      write(nofile,*)' '
      write(nofile,*)'xrtpb:'
      write(nofile,*)'-----'

      imoto = 1
      do ireg = 1, nreg
      isubr = nsum(ireg,ns)
      iel = i1st(isubr) + iele(isubr)
        do im = 1, nmode(iel)
        nwr = ceilq(5,4)
          do i = 1, nwr
          write(nofile,182)(j+i-1,xrtpb(j+i-1,ireg,im,1,1),j=1,5+1-i,nwr)
          end do
        end do
      end do
      end if
 
      if(nscree .gt. 0) then
      write(nofile,*)'Number of Faraday shields =',NSCREE
        if(monomo)then
        write(nofile,*)'Screen current is modulated by mode M=',moscr(1)
        else
        write(nofile,*)'Screen current is modulated by modes M=', moscr(1), ' to m=', moscr(nmoscr)
        end if

c     Enforce screen constraints if not an old solution: check again!
c     This modifies XRTPB and solution vector:
      if(.not.ploold)call addscr
      else
      write(nofile,*)'No Faraday shield.'
      end if

C@    VOIR:
C     SOL. DEJA ECRITE!
C     IF(KEEPSO)WRITE(NOLOFI,1003)((X(J,K),J=1,NGELIM),K=1,neffan)
 
c-------------------------------------c
c     Antenna input impedance matrix  c
c-------------------------------------c
c !To update if HEC touches antenna! + feeders in general geom
      call zset(maxtom*maxant**2, czero, zimp, 1)
      call zset(maxtom*maxant**2, czero, zimpfe, 1)
      do 572 imoto = 1, ntotor
c     ************************
        if(cyl)then
        kphi = ktoan(imoto)
        else
        n = motoan(imoto)
        kphi = n * r0i
        end if
 
c       Contributions of antenna straps:
        do ja2 = 1, neffan
        ireg = irbant(ja2)
        na = ilael(ireg)
        rhoant = rx0m(ireg)
          if(circ)then
c         Circular concentric case:
          fac = - twopi ** 2 * rhoant
            do ja1 = 1, neffan
              do i = 1, nmoant
              im = moant(i) + 1 - minf(na)
              zimp(ja1,ja2,imoto) = zimp(ja1,ja2,imoto)
     ;        + fac * rcura(i,ja2) * dconjg(xrtpb(2,ireg,im,ja1,imoto))
c             (For poloidal antennae)
              end do
            end do
	    else
c         D shape and general geometry:
c to do!!!
          fac = - twopi ** 2 * rhoant
            do ja1 = 1, neffan
              do i = 1, nmoant
              im = moant(i) + 1 - minf(na)
              zimp(ja1,ja2,imoto) = zimp(ja1,ja2,imoto)
     ;        + fac * rcura(i,ja2) * dconjg(xrtpb(2,ireg,im,ja1,imoto))
c             (For poloidal antennae)
              end do
            end do
	    end if
        end do
        if(samant)then
c       Deals with array of identical antennae:
          do ja2 = 1, ntoant
            do ja1 = 1, ntoant
            zimp(ja1,ja2,imoto) = zimp(1,1,imoto)
            end do
          end do
        end if

c       Feeder contributions:
        do 533 ia = 1, neffan
c       =====================
          if(feeder(ia))then
c         NB: geom.current variation along feeders is canceled by volume element
            do m = minf(imax), msup(imax)
            mr = m + 1 - minf(imax)
              if(falen(ia))then
                if(anetyp(ia).eq.'opc')then
c               Inlet contribution only:
                polfac(mr) = - cdeid(-m*thea1(ia))
                else if(anetyp(ia).eq.'shc')then
c               Outlet and inlet contributions:
                polfac(mr) = cdeid(-m*thea2(ia)) - cdcos(betal(ia)) * cdeid(-m*thea1(ia))
                end if
              else
              polfac(mr) =
     ;        cdeid(-m*theaa(ia)) * (-2.d0*ci * dsin(m*dthea(ia)*0.5d0))
              end if
            end do
 
          ireg = irbant(ia)
          ireg1 = ireg + 1
          rhoant = rx0m(ireg)
c         torfac = tofac(ia)
C         Toroidal delta-function antenna at phi=0:
          torfac = cun
 
        do 34 ireg = ireg1, nreg
c       ------------------------
        ifirse = ifiel(ireg)
        is1 = nsum(ireg-1, ns) + 1
        is2 = is1 + ns(ireg) - 1
        ithoma = ifiel(ireg)+2*(ireg-1)
        call rdsol(imoto, ithoma, bel, maxbll, lbl, nrhss)

        do 34 isubr = is1, is2
c       ----------------------
        eletyp = styp(isubr)
        icolo = iconn(isubr)
        ibulo = ibub(isubr)
        ie1 = i1st(isubr)+1
 
        do 34 iel = ie1, i1st(isubr) + iele(isubr)
c       ------------------------------------------
        rle = fl(iel)
        ielm = iel
        if(nmode(iel-1) .gt. nmode(iel))ielm = iel - 1
        nmodc = nmode(ielm)

          if( nmodc.gt.0 )then
c         ....................
          ithoma = iel + 1 + 2 * (ireg - 1)
          call rdsol(imoto, ithoma, ber, maxbll, lbl, nrhss)

            do 44 ja = 1, neffan
            sum1 = czero
            m1 = minf(ielm)
            m2 = msup(ielm)
            sum2 = czero

              do 35 m = m1, m2
              inl = nmode(iel-1).gt.0 .and. m.ge.minf(iel-1) .and. m.le.msup(iel-1)
              mrg = m + 1 - minf(iel-1)
              mr = m + 1 - m1
              mabs = m + 1 - minf(imax)
              inr = nmode(iel).gt.0 .and. m.ge.minf(iel) .and. m.le.msup(iel)
              mrd = m + 1 - minf(iel)
              idecax = (mrg - 1) * icolo
              idecam = nmode(iel-1) * icolo + (mr - 1) * ibulo + 1
              idecar = (mrd - 1) * icolo
c general: geom. factors missing; following assumes rf field has a constant
c multiplier in the integrand.
                if(eletyp .eq. 'M23')then
                sum3 = 4.d0 * bel(idecam,ja)
                if(inl)sum3 = sum3 + bel(1+idecax, ja)
                if(inr)sum3 = sum3 + ber(1+idecar, ja)
                polfa2 = polfac(mabs) * rle / 6.d0
                else if(eletyp .eq. 'HEC')then
                sum3 = czero
                if(inl)sum3 = sum3
     ;           +               bel(1+idecax,ja) + bel(3+idecax,ja)
     ;           + rle / 6.d0 * (bel(2+idecax,ja) + bel(4+idecax,ja))
                if(inr)sum3 = sum3
     ;           +               ber(1+idecar,ja) + ber(3+idecar,ja)
     ;           - rle / 6.d0 * (ber(2+idecar,ja) + ber(4+idecar,ja))
                polfa2 = polfac(mabs) * rle * st2t2i
                end if
              sum2 = sum2 + polfa2 * conjg(sum3)
 35           continue

c           Sum of feeder contributions over relevant elements:
            zimpfe(ja,ia,imoto) = zimpfe(ja,ia,imoto) - torfac * sum2 * rnorm
 44         continue
c
c           Shift right node to left for next element:
            do j = 1, nrhs
              do i = 1, lblock(ithoma)
              bel(i,j) = ber(i,j)
              end do
            end do
          end if
c         ......
  
 34     continue
c       --------

            end if
 533      continue
c         ========
          do ja = 1, neffan
            do ia = 1, neffan
            zimp(ja,ia,imoto) = zimp(ja,ia,imoto) + zimpfe(ja,ia,imoto)
            end do
          end do
          if(samant)then
c         Deals with array of identical antennae:
            do ja = 1, ntoant
              do ia = 1, ntoant
              zimpfe(ja,ia,imoto) = zimpfe(1,1,imoto)
              zimp(ja,ia,imoto) = zimp(1,1,imoto)
              end do
            end do
          end if

        do ia = 1, ntoant
          do ja = 1, ntoant
          zimp2(ia,ja,imoto) = zimp(ia,ja,imoto)
          end do
        end do
        
        do 57 ia = 1, ntoant
          if(falen(ia))then
          cbl = cdcos(betal(ia))
            if(cbl.ne.(0.d0,0.d0))then
            cbli = 1. / cbl
            else
            cbli = 1.d99
            end if
          end if
        do 57 ja = 1, ntoant
        if(falen(ia) .and. anetyp(ia).eq.'SHC')
     ;  zimp(ja,ia,imoto) = zimp(ja,ia,imoto) * cbli
          if(falen(ja) .and. anetyp(ja).eq.'SHC')then
          cbl2 = dconjg(cdcos(betal(ja)))
            if(cbl2.ne.czero)then
            cbl2i = 1.d0 / cbl2
            else
            cbl2i = 1.d99
            end if
          zimp(ja,ia,imoto) = zimp(ja,ia,imoto) * cbl2i
          end if
c       zimpto(ja,ia) = zimpto(ja,ia) + zimp(ja,ia) * coeimp
        zimp3(ja,ia,imoto) = zimp(ja,ia,imoto)
 57     continue

        do ia = 1, ntoant
        toi = dconjg(tofac(ia))
          do ja = 1, ntoant
          toj = toi * tofac(ja)
          zimp(ia,ja,imoto) = zimp(ia,ja,imoto) * toj
          zimp2(ia,ja,imoto) = zimp2(ia,ja,imoto) * toj
          zimpfe(ia,ja,imoto) = zimpfe(ia,ja,imoto) * toj
          end do
        end do
  572 continue
c     ********

c     HERE, BUILD WEIGHT FACTORS FOR INTERPOLATION
c     AND TOTAL IMPEDANCE ZIMPTO in 3D runs
 
      write(nofile,*)'A priori antenna parameters: (used when FALEN=.T. only) '
      write(nofile,*)'----------------------------'
      write(nofile,*)'Antenna ; a priori impedance (eng.convention);'
     ; // ' char.impedance ; electrical length'
        do i = 1, ntoant
        if(falen(i))write(nofile,*)i, zap(i), zcapri(i), betal(i)
        end do
 
      write(nofile,*)'Computed antenna parameters:'
      write(nofile,*)'============================'
      write(nofile,*)'a) Delta-function antennae at PHI=0:'
      write(nofile,*)'   --------------------------------'
      write(nofile,*)'Array partial input impedance matrices:'
      do imoto = 1, ntotor
        do i = 1, ntoant
        write(nofile,500)(zimp3(i,j,imoto), j = 1, ntoant)
        end do
      end do
      
      write(nofile,*)'b) Given array 2D impedance matrices:'
      write(nofile,*)'   ---------------------------------'
      write(nofile,*)'Array partial impedance matrix in "large current'
     ; // ' basis":'
      do imoto = 1, ntotor
        do i = 1, ntoant
        write(nofile,500)(zimp2(i,j,imoto), j = 1, ntoant)
        end do
      end do

      write(nofile,*)'Array partial impedance matrix in "large current'
     ; // ' basis": feeder contribution'
      do imoto = 1, ntotor
        do i = 1, ntoant
        write(nofile,500)(zimpfe(i,j,imoto), j = 1, ntoant)
        end do
      end do

      write(nofile,*)'Input impedance matrices:'
      do imoto = 1, ntotor
        do i = 1, ntoant
        write(nofile,500)(zimp(i,j,imoto), j = 1, ntoant)
        end do
      end do

c     to write: zimpto for 3d runs
 
      entry repou2
C     ************
c     This entry point to repeat an output with different antenna
c     currents.
      write(nofile,*)'enter repout ; time=',second()-timin
      write(nofile,*)'Antenna basis currents:'
      write(nofile,*)(tancur(i),i=1,ntoant)
 
      call dset(maxpom, 0.d0, monorm, 1)
      call dset(2*maxpom, 0.d0, ancplo, 1)

        if(plot2d.or.plopow.or.transp)then
c       Table of exp(i*m*theta), assuming uniform poloidal mesh on (0,2pi):
          do m = minf(imax), msup(imax)
          ima = m - minf(imax) + 1
          tabexp(ima,1) = cun
          tabexp(ima,nploth+1) = cun
cPL 3/05/04            tabexp(ima,ipol) = cdeid(m * polplo(ipol))
	    tdelta = dfloat(m) * twopi / dfloat(nploth)
	    talpha = 2.d0 * dsin(0.5d0*tdelta)**2
	    tbeta = dsin(tdelta)
            do ipol = 2, nploth
	      corr = dcmplx(talpha, -tbeta) * tabexp(ima,ipol-1)
            tabexp(ima,ipol) = tabexp(ima,ipol-1) - corr
cPL 3/05/04            tabexp(ima,ipol) = cdeid(m * polplo(ipol))
            end do
          end do
c       Table of exp(i*m*theta), assuming poloidal mesh equidistant on (0,2pi):
c        tr = CDEID(POLPLO(2))
c        TABEXP(1,1) = cun
c        TABEXP(1,nploth+1) = cun
c        tabexp(1,2) = tr ** minf(imax)
c          do IMA = 2, NMODE(IMAX)
c          TABEXP(ima,1) = cun
c          TABEXP(ima,nploth+1) = cun
c          tabexp(ima,2) = tr * tabexp(ima-1,2)
c          end do
c          do IPOL = 3, nploth
c            do ima = 1, nmode(imax)
c            TABEXP(IMA,IPOL) = tabexp(ima,2) * tabexp(ima,ipol-1)
c            end do
c          end do
        END IF
 
      trncp1 = 3 * ncomp + 1
      trncp2 = 3 * ncomp + 2
      trncp5 = 3 * ncomp + 5
      call dset(nabplo, 0.d0, dpoydr, 1)
cPL not used:      call dset(nabplo*ndof, 0.d0, yana, 1)

c     For antenna current spectrum plot: 
      wmin = 0.5 * dza(idamin(ntoant,dza,1))
c     wmin = dmin1(wmin, zaa(idamin(ntoant,zaa,1)))
      if(wmin.eq.0.d0)wmin = 0.1
c     Ad hoc limit:
      kmax = twopi / wmin
      kstep = 2. * kmax / dfloat(nptosp-1)
        if(cyl)then
c       Abscissa: m**-1
        absmax = kmax
        absstp = kstep
        else
c       Abscissa: mode index; using outer wall equatorial major radius:
c        absmax = kmax * eqt(nabsci,1,1)
c        absstp = max(1.d0, kstep * eqt(nabsci,1,1))
        absmax = (nptosp - 1) / 2
        absstp = 1
        end if
        do i = 1, nptosp
        abstom(i) = - absmax + (i - 1) * absstp
        end do
        do i = 1, nptosp
        ztra = czero
          do ia = 1, ntoant
            if(cyl)then
            torfac = cdeid(- abstom(i) * zaa(ia))
            xtra   = abstom(i) * dza(ia) * 0.5
            else
            torfac = cdeid(- abstom(i) * phiaa(ia))
            xtra   = abstom(i) * dphia(ia) * 0.5
            end if
          torfac = torfac * dsinc(xtra)
          ztra = ztra + tancur(ia) * torfac
          end do
        antosp(i) = cdabs(ztra)
        end do
 
c     Compute total input power: for shorted antenna, this uses the short-circuit current.
      tinpcp = czero
        do imoto = 1, ntotor
        call mucrvz(ntoant, ntoant, zimp2(1,1,imoto), maxant,
     ;              ntoant, tancur, 1, ntoant, ctra2)
        tinpcp = tinpcp + 0.5d0 * zdotc(ntoant, tancur, 1, ctra2, 1)
        end do
      write(nofile,*)'Total input complex power=', TINPCP, ' VA'
      tinpow = dreal(tinpcp)
      losles = dabs(tinpow) .le. 1.d-10 * dabs(dimag(tinpcp))
        if(losles)then
        write(nofile,*)'No significant input active power for this run'
        else
        if(.not.glovac .and. tinpow.lt.0.d0)write(nofile,*)'caution: plasma is feeding antennas'
        end if

c--------------------------------------------------------------------------c
c     Compute field modal components at plot abscissae:                    c
c     Electric ones in all cases, magnetic ones in circular.               c
c     For type HEC, second derivatives are also computed, stored in XPMP2  c
c     Local magnetic field in DSHAPE is computed in a following section.   c 
c--------------------------------------------------------------------------c

C     Once for all elements, basis functions at plot abscissae:
cERN      if(outgau)then  
        do i = 0, ngauss + 1
        call bafad(aga(i), bafplo(1,i))
          do it = 1, ntyp
          call bafadn(it, aga(i), bafpln(1,0,it,i), bafpls(1,0,it,i))
          end do
        end do
cERN      else
cERN        fnin1i = 1.d0 / dfloat(ninter+1)
cERN        do i = 0, ninter + 1
cERN        tr = i * fnin1i
cERN        call bafad(tr, bafplo(1,i))
cERN        end do
cERN      end if
 
      call dset(nabplo*(3*nficom+5), 0.d0, yout, 1)
      call dset(nabplo*3*nficom, 0.d0, youth, 1)

c     Open file for solutions storage: unitgl was closed in main:
      if(.not.monoto)then
	aqp(1) = unitgl
	call aqopen(aqp, 3, unitgl, istatu)
      end if
	iblk = 1
      ireq = 0

      do 304 imoto = 1, ntotor
c     ========================
        if(cyl)then
        kphi = ktoan(imoto)
        else
        n = motoan(imoto)
        kphi = n * r0i
        end if
      kprn = kphi * rnorm
      kprn2 = kprn * kprn
c       Antenna toroidal factors for current value of kphi or n, times
c       antenna basis current:
        do ia = 1, ntoant
        tofacv(ia) = tofac(ia) * tancur(ia)
        end do

      if(imoto.gt.1)call aswa(1)
      call zset((ndof-1)*maxpom*nabplo, czero, xrtp, 1)
      call zset(3*maxpom*nabplo, czero, hrtp, 1)
      call zset(ndof*maxpom*nabplo, czero, xpmp, 1)
      call zset(nficom*maxpom*nabplo, czero, xpmp2, 1)
      call zset(maxbll*maxrhs, czero, bel, 1)
      call zset(maxbll*maxrhs, czero, ber, 1)

cccccccccccccc
      i = 0
cccccccccccccc

      do 11 ireg = 1, nreg
c     ====================
      is1 = nsum(ireg-1, ns) + 1
      is2 = is1 + ns(ireg) - 1
c     Starting diag. block index:
      ithoma = ifiel(ireg) + 2 * (ireg - 1)
c     Read corresp. solution vector block ('left part' of first element of this
c     region):
      call rdsol(imoto, ithoma, bel, maxbll, lbl, nrhss)

      do 11 isubr = is1, is2
c     ======================
      it = istyp(isubr)
      eletyp = styp(isubr)
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo
      
      ie1 = i1st(isubr) + 1
      ie2 = ie1 + iele(isubr) - 1

      do 11 iel = ie1, ie2
c     ====================
c     Current diag. block index:
      ithoma = iel + 1 + 2 * (ireg - 1)
c     Read corresp. solution vector block ('right part' of current element):
      call rdsol(imoto, ithoma, ber, maxbll, lbl, nrhss)
      rle = fl(iel)
        do ip = -2, 1
        rlpow(ip) = rle ** ip
        end do
      rlei = rlpow(-1)
cERN        if(outgau)then  
          if(iel.eq.ifiel(ireg))then
c         Include left node for first element in this region:
          j1 = 0
          else
c         Don't for others:
          j1 = 1
          end if
c         Always include right node:
          j2 = ngauss+1
cERN        else
cERN          if(iel.eq.ilael(ireg))then
cERN          nintma = ninter + 1
cERN          else
cERN          nintma = ninter
cERN          end if
cERN        j1 = 0
cERN        j2 = nintma
cERN        end if
 
        if(nmode(iel-1).gt.nmode(iel))then
        ielm = iel - 1
        else
        ielm = iel
        end if
      nmodc = nmode(ielm)
c
        do 305 j = j1, j2
c       -----------------
        i = i + 1

        do 305 ja = 1, ntoant
c       ---------------------
        toj = tofacv(ja)
        jae = ja
        if(samant)jae = 1
c       Dealing with eventual radial variations of number of pol. modes:
c       Copy solution vector into reference vector with nmodc pol. modes.
        call zset(2*maxbll, czero, solvec, 1)
        idecax = 1 + (nmodc - nmode(iel-1)) * icolo
        call zcopy(nmode(iel-1)*icolo, bel(1,jae), 1, solvec(idecax), 1)
          if(ibulo .gt. 0)then
          idecam = 1 + nmodc * icolo
          call zcopy(nmodc*ibulo, bel(nmode(iel-1)*icolo+1,jae), 1, solvec(idecam), 1)
          end if
        idecar = 1 + nmodc * icpblo + (nmodc - nmode(iel)) * icolo
        call zcopy(nmode(iel)*icolo, ber(1,jae), 1, solvec(idecar), 1)
        
          if(eletyp.eq.'HEC')then
            do jd = 1, nficom
            j22 = 2 * jd
            j22m1 = j22 - 1
c             Contributions from left node:
              do img = 1, nmode(iel-1)
              m = img - 1 + minf(iel-1)
              im = m + 1 - minf(ielm)
              idecax = (img - 1) * icolo
              j22i   =   j22 + idecax
              j22m1i = j22m1 + idecax
c             E+,-,//:
              xpmp(j22m1,i,im) = xpmp(j22m1,i,im) + toj * (
     ;        bel(j22m1i,jae) * bafplo(1,j)
     ;      + bel(j22i,jae)   * bafplo(2,j) * rlpow(1)    )
c             d/dy E+,-,//:
              xpmp(j22,i,im) = xpmp(j22,i,im) + toj * (
     ;        bel(j22m1i,jae) * bafplo(5,j) * rlpow(-1)
     ;      + bel(j22i,jae)   * bafplo(6,j)           )
c             d2/dy2 E+,-,//:
              xpmp2(jd,i,im) = xpmp2(jd,i,im) + toj * rlei * (
     ;        bel(j22m1i,jae) * bafplo(15,j) * rlpow(-1)
     ;      + bel(j22i,jae)   * bafplo(16,j)                 )
              end do

c             Contributions from right node:
              do imd = 1, nmode(iel)
              m = imd - 1 + minf(iel)
              im =  m + 1 - minf(ielm)
              idecar = (imd - 1) * icolo
              j22i = j22 + idecar
              j22m1i = j22m1 + idecar
c             E+,-,//:
              xpmp(j22m1,i,im) = xpmp(j22m1,i,im) + toj * (
     ;        ber(j22m1i,jae) * bafplo(3,j)
     ;      + ber(j22i,jae)  *  bafplo(4,j) * rlpow(1)    )
c             d/dy E+,-,//:
              xpmp(j22,i,im) = xpmp(j22,i,im) + toj * (
     ;      + ber(j22m1i,jae) * bafplo(7,j) * rlpow(-1)
     ;      + ber(j22i,jae)  *  bafplo(8,j)            )
c             d2/dy2 E+,-,//:
              xpmp2(jd,i,im) = xpmp2(jd,i,im) + toj * rlei * (
     ;        ber(j22m1i,jae) * bafplo(17,j) * rlpow(-1)
     ;      + ber(j22i,jae)   * bafplo(18,j)                 )
              end do
            end do

          else if(eletyp.eq.'M23')then
c           Contributions from left, central and right nodes to Erho:
            do img = 1, nmode(iel-1)
            m = img - 1 + minf(iel-1)
            im = m + 1 - minf(ielm)
            idecax = (img - 1) * icolo + 1
            xrtp(1,i,im) = xrtp(1,i,im) + toj * bel(idecax,jae) * bafplo(9,j)
            end do

            do im = 1, nmodc
            idecam = nmode(iel - 1) * icolo + (im - 1) * ibulo + 1
            xrtp(1,i,im) = xrtp(1,i,im) + toj * bel(idecam,jae) * bafplo(10,j)
            end do

            do imd = 1, nmode(iel)
            m = imd - 1 + minf(iel)
            im = m + 1 - minf(ielm)
            idecar = (imd - 1) * icolo + 1
            xrtp(1,i,im) = xrtp(1,i,im) + toj * ber(idecar,jae) * bafplo(11,j)
            end do
c           Magnetic axis: contributions from left, central and right nodes to dErho/dy:
            if(iel.eq.1 .and. j.eq.0)then
              do im = 1, nmodc
              erpa(im) = czero
              end do

              do img = 1, nmode(iel-1)
              m = img - 1 + minf(iel-1)
              im = m + 1 - minf(ielm)
              idecax = (img - 1) * icolo + 1
              erpa(im) = erpa(im) + toj * bel(idecax,jae) * bafpln(1,1,2,j) * rlei
              end do

              do im = 1, nmodc
              idecam = nmode(iel - 1) * icolo + (im - 1) * ibulo + 1
              erpa(im) = erpa(im) + toj * bel(idecam,jae) * bafpln(2,1,2,j) * rlei
              end do
              
              do imd = 1, nmode(iel)
              m = imd - 1 + minf(iel)
              im = m + 1 - minf(ielm)
              idecar = (imd - 1) * icolo + 1
              erpa(im) = erpa(im) + toj * ber(idecar,jae) * bafpln(3,1,2,j) * rlei
              end do
            end if
 
            do jd = 1, (ndof - 2) / 2
            j22 = 2 * jd
            j22p1 = j22 + 1
 
c             Contributions from left node to other components:
              do img = 1, nmode(iel-1)
              m = img - 1 + minf(iel-1)
              im = m + 1 - minf(ielm)
              idecax = (img - 1) * icolo
              j22i = j22 + idecax
              j22p1i = j22p1 + idecax
c             Etheta, Ephi:
              xrtp(j22,i,im) = xrtp(j22,i,im) + toj * (
     ;        bel(j22i,jae)   * bafplo(1,j)
     ;      + bel(j22p1i,jae) * bafplo(2,j) * rle     )
c             d/dy Etheta, d/dy Ephi:
              xrtp(j22p1,i,im) = xrtp(j22p1,i,im) + toj * (
     ;        bel(j22i,jae)   * bafplo(5,j) * rlei
     ;      + bel(j22p1i,jae) * bafplo(6,j)               )
              end do
 
c             Contributions from right node:
              do imd = 1, nmode(iel)
              m = imd - 1 + minf(iel)
              im = m + 1 - minf(ielm)
              idecar = (imd - 1) * icolo
              j22i = j22 + idecar
              j22p1i = j22p1 + idecar
c             Etheta, Ephi:
              xrtp(j22,i,im) = xrtp(j22,i,im) + toj * (
     ;        ber(j22i,jae)  * bafplo(3,j)
     ;      + ber(j22p1i,jae)* bafplo(4,j) * rle      )
c             d/dy Etheta, d/dy Ephi:
              xrtp(j22p1,i,im) = xrtp(j22p1,i,im) + toj * (
     ;        ber(j22i,jae)   * bafplo(7,j) * rlei
     ;      + ber(j22p1i,jae) * bafplo(8,j)               )
              end do

            end do
          end if
 305    continue
c       --------

c       Shift right to left for next element:
        do jj = 1, nrhs
          do ii = 1, lblock(ithoma)
          bel(ii,jj) = ber(ii,jj)
          end do
        end do

  11  continue
c     ========

      if(wrisol)then
c     Rho, theta, phi components of solution in 'M23' subregions;
c     zero in 'HEC' subregions:  
      write(nofile,*)' '
      write(nofile,*)'xrtp:'
      write(nofile,*)'-----'

      do ic = 1, 5
        do im = 1, nmode(2)
        write(nofile,*)'ic=', ic, '; im=', im
          do i = 1, ist11
          write(nofile,2000)
     ;    absr(i), dreal(xrtp(ic,i,im)), dimag(xrtp(ic,i,im))
          end do
        end do
      end do

c      nwr = ceilq(ist11,4)
c      do 183 ic = 1, 5
c      do 183 im = 1, nmode(2)
c      write(nofile,*)'ic=', ic, '; im=', im
c      do 184 i = 1, nwr
c      write(nofile,182)(j+i-1, xrtp(ic,j+i-1,im), j=1, ist11+1-i,nwr)
c184   continue
c183   continue
      end if
 
c------------------------------------------------------------------c
c     Circular case:                                               c
c     Coord. transformation (rtp to +-// or vice versa) for modes  c
c     and calculation of magnetic field modes;                     c
c                                                                  c
c     CIRC  and DSHAPE: elec. and magn. field at theta=0, phi=0    c
c------------------------------------------------------------------c
      fac3 = dcmplx(0.d0, -1.d0 / (omegag * mu0 * rnorm))
      i = 0
      do ireg = 1, nreg
c     =================
      is1 = nsum(ireg-1, ns) + 1
      is2 = is1 + ns(ireg) - 1

      do isubr = is1, is2
c     ====================
      eletyp = styp(isubr)
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo
      ie1 = i1st(isubr) + 1
      ie2 = ie1 + iele(isubr) - 1
 
      do iel = ie1, ie2
c     ==================
      rle = fl(iel)
cERN      if(outgau)then  
		if(iel.eq.ifiel(ireg))then
c         Include left node for first element in this region:
          j1 = 0
          else
c         Don't for others:
          j1 = 1
          end if
c         Always include right node:
          j2 = ngauss+1
cERN      else
cERN        if(iel.eq.ilael(ireg))then
cERN           nintma = ninter + 1
cERN        else
cERN           nintma = ninter
cERN        end if
cERN        j1 = 0
cERN        j2 = nintma
cERN      end if

      ielm = iel
      if(nmode(iel-1) .gt. nmode(iel))ielm = iel - 1
      nmodc = nmode(ielm)
      minfc = minf(ielm)
      msupc = msup(ielm)
      
        do j = j1, j2
c       -------------
        i = i + 1
c        intab = intaplot(i)
        intab = i
      
        y = absono(i)
        yinv = abscni(intab)
        call zset(ndof-1, czero, eloc, 1)
        call zset(ndof, czero, epmloc, 1)
        call zset(nficom, czero, hloc, 1)
        call zset(nficom, czero, elocth, 1)
 
      if(circ)then
c     ~~~~~~~~~~~~
      call loctr2
      call locti2
      call fouco
 
        if (eletyp .eq. 'HEC') then
c       Transform modes from +-// to r,t,p:
          do jc = 1, ndof
            do ic = 1, ndof - 1
            c1 = bc2(ic,jc)
            call zaxpy(nmodc, c1, xpmp(jc,i,1), ndof*nabplo, xrtp(ic,i,1), (ndof-1)*nabplo)
c             do im = 1, nmodc
c             xrtp(ic,i,im) = xrtp(ic,i,im) + c1 * xpmp(jc,i,im)
c             end do
            end do
          end do
        else if (eletyp .eq. 'M23') then
c       Transform modes from r,t,p to +-//:
cPL21/01/05: NB be aware that the contributions from Erho' to E+' and E-' are missing for the M23 type.
          do jc = 1, ndof - 1
            do ic = 1, ndof
            c1 = bc2i(ic,jc)
            call zaxpy(nmodc, c1, xrtp(jc,i,1), (ndof-1)*nabplo, xpmp(ic,i,1), ndof*nabplo)
c             do im = 1, nmodc
c             xpmp(ic,i,im) = xpmp(ic,i,im) + c1 * xrtp(jc,i,im)
c             end do
            end do
          end do
        end if
 
c     Compute magnetic field modes, and local electric and magnetic components at theta=0, phi=0:
      if(.not.cokpco)then  ! standard coordinates
c     ...................
      do im = 1, nmodc
      m = im - 1 + minfc
c  a) Terms which contribute to k=0 only:
      hrtp(2,i,im) = hrtp(2,i,im) - xrtp(5,i,im)
      hrtp(3,i,im) = hrtp(3,i,im) + xrtp(3,i,im)

        if(y .eq. 0.) then
        hrtp(1,i,im) = hrtp(1,i,im) + ci * m * xrtp(5,i,im)
cPL20/01/05: implement the limit (i omega mu0) * Hphi(0) = d/dy(2*Etheta-im*Erho)
          if(styp(1) .eq. 'HEC')then
          hrtp(3,i,im) = (-ci*sqrt2i) * ((m+2.d0) * xpmp(2,i,im) + (m-2.d0) * xpmp(4,i,im))
	    else if(styp(1) .eq. 'M23')then
c         dErho/dy(0) previously stored in erpa:
          hrtp(3,i,im) = 2.d0 * xrtp(3,i,im) - ci * m * erpa(im)
	    end if
        else
        hrtp(1,i,im) = hrtp(1,i,im) + ci*m*yinv * xrtp(4,i,im)
        hrtp(3,i,im) = hrtp(3,i,im) + (xrtp(2,i,im) - ci * m * xrtp(1,i,im)) * yinv
        end if

c  b) Terms which couple different harmonics:
      k1 = - min0( ncrot, m-minf(ielm) )
      k2 =   min0( ncrot, msup(ielm)-m )
        do k = k1, k2
        kp = iabs(k)
        impk = im + k
        hrtp(1,i,im) = hrtp(1,i,im) - ci * kprn * bfou(kp) * xrtp(2,i,impk)
        hrtp(2,i,im) = hrtp(2,i,im) + ci * kprn * bfou(kp) * xrtp(1,i,impk) - cfou(kp) * xrtp(4,i,impk)
cPL20/01/05
        if(y .ne. 0.d0)then
cPL20/01/05 dfou is odd!
cPL20/01/05    hrtp(1,i,im) = hrtp(1,i,im) + ci * yinv * dfou(kp) * xrtp(4,i,impk)
c        hrtp(1,i,im) = hrtp(1,i,im) + ci * yinv * dfou(kp) * dsign(1.d0,dfloat(k)) * xrtp(4,i,impk)
cPL28/01/05 using negative harmonics:
        hrtp(1,i,im) = hrtp(1,i,im) + ci * yinv * dfou(k) * xrtp(4,i,impk)

cPL20/01/05: added the correct limit at y=0:
	  else if(kp.eq.1)then
        hrtp(1,i,im) = hrtp(1,i,im) - ci * 0.5 * bfou(0) * rnorm * r0i* dsign(1.d0,dfloat(k)) * xrtp(4,i,impk)         
	  end if
        end do 

      hrtp(1,i,im) = hrtp(1,i,im) * fac3
      hrtp(2,i,im) = hrtp(2,i,im) * fac3
      hrtp(3,i,im) = hrtp(3,i,im) * fac3
 
      end do

      else      ! cokpco case
c     ....
      do im = 1, nmodc
      m = im - 1 + minfc

c  a) Terms which contribute to k=0 only:
      hrtp(2,i,im) = hrtp(2,i,im) - xrtp(5,i,im)
      hrtp(3,i,im) = hrtp(3,i,im) + xrtp(3,i,im)

        if(y .eq. 0.) then
        hrtp(1,i,im) = hrtp(1,i,im) + ci * m * xrtp(5,i,im)
cPL20/01/05: implement the limit (i omega mu0) * Hphi(0) = d/dy(2*Etheta-i*m*Erho) -in*q0/R0(Erho*cos(theta)-Etheta*sin(theta))
c@CHECK again
          if(styp(1) .eq. 'HEC')then
          hrtp(3,i,im) = (-ci*sqrt2i) * ((m+2.d0) * xpmp(2,i,im) + (m-2.d0) * xpmp(4,i,im))
	    else if(styp(1) .eq. 'M23')then
c         dErho/dy(0) previously stored in erpa:
          hrtp(3,i,im) = 2.d0 * xrtp(3,i,im) - ci * m * erpa(im)
	    end if
        else
        hrtp(1,i,im) = hrtp(1,i,im) + ci*m*yinv * xrtp(4,i,im)
        hrtp(3,i,im) = hrtp(3,i,im) + (xrtp(2,i,im) - ci * m * xrtp(1,i,im)) * yinv
        end if

c  b) Terms which couple different harmonics:
      k1 = - min0( ncrot, m-minf(ielm) )
      k2 =   min0( ncrot, msup(ielm)-m )
        do k = k1, k2
        kp = iabs(k)
        impk = im + k
        hrtp(1,i,im) = hrtp(1,i,im) - ci * kprn * bfou(kp) * xrtp(2,i,impk)
        hrtp(2,i,im) = hrtp(2,i,im) + ci * kprn * bfou(kp) * xrtp(1,i,impk) - (cfou(kp) + ci * n * dcokfou(k)) * xrtp(4,i,impk)
	  hrtp(3,i,im) = hrtp(3,i,im) + ci * n * dcokfou(k) * xrtp(2,i,impk)
cPL20/01/05:
        if(y .ne. 0.d0)then
cPL20/01/05 dfou is odd!
cPL20/01/05    hrtp(1,i,im) = hrtp(1,i,im) + ci * yinv * dfou(kp) * xrtp(4,i,impk)
c        hrtp(1,i,im) = hrtp(1,i,im) + ci * yinv * (dfou(kp) * dsign(1.d0,dfloat(k)) + n * acokfou(k)) * xrtp(4,i,impk)
cPL28/01/05 using negative harmonics:
        hrtp(1,i,im) = hrtp(1,i,im) + ci * yinv * (dfou(k) + n * acokfou(k)) * xrtp(4,i,impk)
	  hrtp(3,i,im) = hrtp(3,i,im) - ci * n * yinv * acokfou(k) * xrtp(1,i,impk)
cPL20/01/05: added the limit at y=0:
	  else if(kp.eq.1)then
        hrtp(1,i,im) = hrtp(1,i,im) - ci * 0.5 * r0i * xrtp(4,i,impk) * (bfou(0) * rnorm * dsign(1.d0,dfloat(k)) - n * qtab(i))
	  hrtp(3,i,im) = hrtp(3,i,im) - ci * (0.5 * r0i * n * qtab(i)) * xrtp(1,i,impk)
	  end if
        end do 

      hrtp(1,i,im) = hrtp(1,i,im) * fac3
      hrtp(2,i,im) = hrtp(2,i,im) * fac3
      hrtp(3,i,im) = hrtp(3,i,im) * fac3
 
      end do
      end if
c     ......
            
c     Electric and magnetic components at theta=0, phi=0:                 (also valid for cokpco=.true.)
      do im = 1, nmodc
        do ii = 1, nficom
        ii2 = 2 ** (ii - 1)
c       eloc has 5 rows, field stored in #1, 2, 4:
        eloc(ii2,1) = eloc(ii2,1) + xrtp(ii2,i,im)
c       hloc has 5 rows, field stored in #1, 2, 3:
        hloc(ii,1) = hloc(ii,1) + hrtp(ii,i,im)
        end do
        do ii = 1, ndof
c       epmloc has 6 rows, field / radial derivative stored resp. in #(1,3,5) / #(2,4,6):
        epmloc(ii,1) = epmloc(ii,1) + xpmp(ii,i,im)
        end do         
      end do

      else
c     ~~~~ (DSHAPE or GENEQ cases)

      call foucog(2)
      
c     Compute magnetic field modes and Poynting flux:
      call magmod(i, poynt(1:nmodc,i))

c       Electric and magnetic components at theta=0, phi=0:
        if(eletyp .eq. 'HEC')then
          do ii = 1, ndof
            do im = 1, nmodc
            epmloc(ii,1) = epmloc(ii,1) + xpmp(ii,i,im)
            end do
          end do
        eloc(1,1) = sqrt2i * (epmloc(1,1) + epmloc(3,1))
        eeta = - ci * sqrt2i * (epmloc(1,1) - epmloc(3,1))
        eloc(2,1) =   eeta * eqt(intab,1,14) + epmloc(5,1) * eqt(intab,1,15)
        eloc(4,1) = - eeta * eqt(intab,1,15) + epmloc(5,1) * eqt(intab,1,14)
c         Hloc in rho, eta, //:
          do ii = 1, nficom
            do im = 1, nmodc
            hloc(ii,1) = hloc(ii,1) + hrtp(ii,i,im)
            end do
          end do
c       (radial derivatives not evaluated)

        else if(eletyp .eq. 'M23')then
          do ii = 1, ndof - 1
            do im = 1, nmodc
c2004       m = im - 1 + minfc
            eloc(ii,1) = eloc(ii,1) + xrtp(ii,i,im)
            end do
          end do
c         Hloc in rho, theta, phi:
          do ii = 1, nficom
            do im = 1, nmodc
            hloc(ii,1) = hloc(ii,1) + hrtp(ii,i,im)
            end do
          end do
        eeta = eloc(2,1) * eqt(intab,1,14) - eloc(4,1) * eqt(intab,1,15)
        epmloc(1,1) = sqrt2i * (eloc(1,1) + ci * eeta)
        epmloc(3,1) = sqrt2i * (eloc(1,1) - ci * eeta)
        epmloc(5,1) = eloc(2,1) * eqt(intab,1,15) + eloc(4,1) * eqt(intab,1,14)
        end if        
         
      end if
C     ~~~~~~
 
      do ii = 1, ncomp
      i2 = 3 * (ii - 1)
cPL2005 NB: check that the following are always increments from 0, can be replaced by value affectation:
        if(ii.le.3)then
C       R,T,P components of electric and magnetic fields:
        ii2 = 2 ** (ii - 1)
        yout(i,i2+1) = yout(i,i2+1) + dreal(eloc(ii2,1))
        yout(i,i2+2) = yout(i,i2+2) + dimag(eloc(ii2,1))
        yout(i,i2+3) = yout(i,i2+3) + cdabs(eloc(ii2,1))
        youth(i,i2+1) = youth(i,i2+1) + dreal(hloc(ii,1))
        youth(i,i2+2) = youth(i,i2+2) + dimag(hloc(ii,1))
        youth(i,i2+3) = youth(i,i2+3) + cdabs(hloc(ii,1))
        else if(ii.eq.4)then
c       E+ component:
        yout(i,i2+1) = yout(i,i2+1) + dreal(epmloc(1,1))
        yout(i,i2+2) = yout(i,i2+2) + dimag(epmloc(1,1))
        yout(i,i2+3) = yout(i,i2+3) + cdabs(epmloc(1,1))
        end if
c     Numerical norm of solution:
      enorm(i,ii) = yout(i,i2+3) ** 2
      end do
       
c  16    continue
            end do  ! j
          end do  ! iel
        end do  ! isubr
      end do  ! ireg
c     ======

      if(.not.monoto)then
      nblk = ceilq(nwn*5*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqwrite(aqp,xrtp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*3*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqwrite(aqp,hrtp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*ndof*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqwrite(aqp,xpmp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*nficom*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqwrite(aqp,xpmp2,iblk,nblk,ireq,0,istatu)
      iblk = iblk + nblk
      end if


	
cERN  ccccccccccc Writing 1D fields to individual files cccccccccccc

	GGs = "(g16.6, g16.6, g16.6)"

	open (UNIT = 701, FILE = paplda // 'Erad_1D.dat', STATUS = "REPLACE", ACTION = "WRITE")
	write(701,*),"Radial electric field"
	open (UNIT = 702, FILE = paplda // 'Epol_1D.dat', STATUS = "REPLACE", ACTION = "WRITE")
	write(702,*),"Poloidal electric field"
	open (UNIT = 703, FILE = paplda // 'Etor_1D.dat', STATUS = "REPLACE", ACTION = "WRITE")
	write(703,*),"Toroidal electric field"
	open (UNIT = 704, FILE = paplda // 'Eplus_1D.dat', STATUS = "REPLACE", ACTION = "WRITE")
	write(704,*),"L.H. electric field"
	open (UNIT = 705, FILE = paplda // 'Hrad_1D.dat', STATUS = "REPLACE", ACTION = "WRITE")
	write(705,*),"Radial magnetic field"
	open (UNIT = 706, FILE = paplda // 'Hpol_1D.dat', STATUS = "REPLACE", ACTION = "WRITE")
	write(706,*),"Poloidal magnetic field"
	open (UNIT = 707, FILE = paplda // 'Htor_1D.dat', STATUS = "REPLACE", ACTION = "WRITE")
	write(707,*),"Toroidal magnetic field"

	  do j = 1, nabsci
	  aux = j
	  write(701,GGs), eqt(aux,1,1) , yout(aux,1), yout(aux,2)
	  write(702,GGs), eqt(aux,1,1) , yout(aux,4), yout(aux,5)		
	  write(703,GGs), eqt(aux,1,1) , yout(aux,7), yout(aux,8)
	  write(704,GGs), eqt(aux,1,1) , yout(aux,10), yout(aux,11)		 
	  write(705,GGs), eqt(aux,1,1) , youth(aux,1), youth(aux,2)
	  write(706,GGs), eqt(aux,1,1) , youth(aux,4), youth(aux,5)
	  write(707,GGs), eqt(aux,1,1) , youth(aux,7), youth(aux,8)		
	  end do

	close(701)
	close(702)
	close(703)
	close(704)
	close(705)
	close(706)	
	close(707)	

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c-----------------------------------------------c
c     Field modulus squared,                    c
c     Inward Poynting flux at plot abscissae:   c
c-----------------------------------------------c
ccc      FAC = CI * twopi ** 2 / (2.D0 * OMEGAG * MU0)
c     Calculation using magn. field:
      fac = - 0.5 * twopi ** 2 * rnorm
      if(.not.cyl)fac = fac * ra
      totpoy = czero
 
      i = 0
      call dset(nabplo*nmode(imax), 0.d0, tablo, 1)
      ist = 1
      do 45 ireg = 1, nreg
c     ====================
      is1 = nsum(ireg-1, ns) + 1
      is2 = is1 + ns(ireg) - 1

      do 36 isubr = is1, is2
c     ======================
      eletyp = styp(isubr)
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo
      ie1 = i1st(isubr) + 1
      ie2 = ie1 + iele(isubr) - 1

      do 36 iel = ie1, ie2
c     ====================
      rle = fl(iel)
      ielm = iel
      if(nmode(iel-1) .gt. nmode(iel)) ielm = iel - 1

cERN        if(outgau)then
          if(iel.eq.ifiel(ireg))then
c         include left node for first element in this region:
          j1 = 0
          else
c         don't for others:
          j1 = 1
          end if
c         always include right node:
          j2 = ngauss+1
cERN        else
cERN          if(iel.eq.ilael(ireg))then
cERN          nintma = ninter + 1
cERN          else
cERN          nintma = ninter
cERN          end if
cERN        j1 = 0
cERN        j2 = nintma
cERN        end if

        do 36 j = j1, j2
c       ----------------
        i = i + 1
        y = absono(i)   ! ERN absono -> abscno

c         Field modulus squared:
          if(eletyp.eq.'M23')then
            do mr = 1, nmode(ielm)
            m = mr - 1 + minf(ielm)
            mra = m + 1 - minf(imax)
            tablo(i,mra) = cabs2(xrtp(1,i,mr)) + cabs2(xrtp(2,i,mr)) + cabs2(xrtp(4,i,mr))
            end do
          else if(eletyp.eq.'HEC')then
            do mr = 1, nmode(ielm)
            m = mr - 1 + minf(ielm)
            mra = m + 1 - minf(imax)
            tablo(i,mra) = cabs2(xpmp(1,i,mr)) + cabs2(xpmp(3,i,mr)) + cabs2(xpmp(5,i,mr))
            end do
          end if

c         Poynting inward flux:
          if(circ)then           ! Circular case: same expression in standard and const. k// coordinates
          call fouco
            do mr = 1, nmode(ielm)
            m = mr - 1 + minf(ielm)
            mra = m + 1 - minf(imax)
            poynt(mr,i) = czero
c           conj(etheta):
            c1 = dconjg(xrtp(2,i,mr))
c           conj(Ephi):
            c2 = dconjg(xrtp(4,i,mr))
c            poynt(mr,i) = poynt(mr,i) - ci*y*kprn * c2 * xrtp(1,i,mr)
 
              if(cyl) then
              k1 = 0
              k2 = 0
              else
              k1 = - min0(1, m-minf(ielm))
              k2 =   min0(1, msup(ielm)-m)
              end if

              do k = k1, k2
              mpkr = mr + k
              ki = iabs(k)
c             Poynting using HRTP:
              poynt(mr,i) = poynt(mr,i) + y * afou(ki) * (c1*hrtp(3,i,mpkr) - c2*hrtp(2,i,mpkr))

c             Poynting using XRTP only:
c                if(ki .eq. 1) then
c                poynt(mr,i) = poynt(mr,i) + y * 0.5d0 * rnor0 * c2 * xrtp(4,i,mpkr)
c                end if
c              poynt(mr,i) = poynt(mr,i) + afou(ki) * (
c     ;        c1 * (xrtp(2,i,mpkr) + y * xrtp(3,i,mpkr) - ci * (m+k) * xrtp(1,i,mpkr)) + c2 * y * xrtp(5,i,mpkr) )
              end do

            poynt(mr,i) = poynt(mr,i) * fac
            end do
            
c         else
c         Was done above in call magmod
          end if

c         Sum over modes, gives total inward Poynting flux:
          do mr = 1, nmode(ielm)
          yout(i,trncp1) = yout(i,trncp1) + dreal(poynt(mr,i))
          yout(i,trncp2) = yout(i,trncp2) + dimag(poynt(mr,i))
          end do 
          
  36    continue
c       ========

cERN  ccccccccccccccccc Writing Poynting Flux ccccccccccccccccccccccc
	open (UNIT = 708, FILE = paplda // 'Poynting_flux.dat', STATUS = "REPLACE", ACTION = "WRITE")
      write(708,*),"Inward Poynting flux"
	  do j = 1,nabsci
        aux = j
	  write(708,GGs), eqt(aux,1,1) , yout(aux,trncp1), yout(aux,trncp2)
	  end do
	close(708)	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN      if(outgau)then
      istg = ist - ngauss
c      istg = ist + 1 - ngauss
        do iel = ifiel(ireg), ilael(ireg)
        rle = fl(iel) * rnorm
        istg = istg + ngauss + 1
c        istg = istg + ngauss
          do im = 1, nmode(imax)
c         check weight rho in integrand!
          monorm(im) = monorm(im) + gauss(rle, absout(istg), tablo(istg,im), 1, 1)   ! ERN absout->abscis
          end do
        end do
cERN      else
cERN        do im = 1, nmode(imax)
cERN        monorm(im) = monorm(im)
cERN     ;  + trapeze(istp(ireg), absout(ist), tablo(ist,im), 1, 1)
cERN        end do
cERN      end if

      ist = ist + istp(ireg)  ! cERN ####### maybe here
c     Difference of inward Poynting flux over each region: real part should
c     equal nett absorbed power in region (needs accounting for eventual feeders
      poydir(ireg) = dcmplx(yout(i,trncp1) - yout(i-istp(ireg)+1,trncp1),
     ;                      yout(i,trncp2) - yout(i-istp(ireg)+1,trncp2))
c     Sum differences of inward Poynting flux over regions: real part should
c     equal total absorbed power and total input power, minus feeder power.
      totpoy = totpoy + poydir(ireg)
  45  continue
c     ========

c     Inward Poynting flux differences across each element:
      open(unit=80, file = paplda // 'Inward Poynting flux differences across each element (OUTRFF)', status='unknown')
      write(80,*)'Inward Poynting flux differences across each element and boundary:'
      write(80,*)'Element, block, rho at elt. centre, real, imag'
      write(80,*)nele + nreg, 5
c     ...and inward Poynting flux at element nodes:
      open(unit=82, file = paplda // 'Inward Poynting flux at element nodes (OUTRFF)', status='unknown')
      write(82,*)'Inward Poynting flux at element nodes'
      write(82,*)'rho at elt. node, real, imag'
      write(82,*)nele + nreg, 3
      i = 0
      call dset(maxnbl, 0.d0, poydie, 1)
        do ireg = 1, nreg
c       Block index in solution vector:
        ibl = ifiel(ireg) + 2 * (ireg - 1)
          do iel = ifiel(ireg), ilael(ireg)
cERN            if(outgau)then
              if(iel.eq.ifiel(ireg))then
c             Include left node for first element in this region:
              j1 = 0
              else
c             Don't for others:
              j1 = 1
              end if
c             Always include right node:
              j2 = ngauss+1
cERN            else
cERN              if(iel.eq.ilael(ireg))then
cERN              nintma = ninter + 1
cERN              else
cERN              nintma = ninter
cERN              end if
cERN            j1 = 0
cERN            j2 = nintma
cERN            end if

c         i points on right node of current element, i-j2 on left node:
          i = i + j2 - j1 + 1
          poydie(ibl) = dcmplx(yout(i,trncp1) - yout(i-j2,trncp1),
     ;                         yout(i,trncp2) - yout(i-j2,trncp2))
          write(80,*)iel, ibl, rnorm*(fx0(iel-1)+0.5*fl(iel)), dreal(poydie(ibl)), dimag(poydie(ibl))
          write(82,*)rnorm*fx0(iel), yout(i,trncp1), yout(i,trncp2)
          ibl = ibl + 1
          end do
c       This is the difference across region boundary:
        if(ireg .lt. nreg)then
        poydie(ibl) = dcmplx(yout(i+1,trncp1) - yout(i,trncp1),
     ;                       yout(i+1,trncp2) - yout(i,trncp2))
        write(80,*)iel, ibl, rnorm*(fx0(iel-1)+0.5*fl(iel)), dreal(poydie(ibl)), dimag(poydie(ibl))
        write(82,*)rnorm*fx0(iel), yout(i+1,trncp1), yout(i+1,trncp2)
	  end if
        end do        
	close(unit=80, status='KEEP')
	close(unit=82, status='KEEP')
c
c      write(75,*)' '
c      write(75,*)'abscissae and table yout:'
c      write(75,*)'------------------------'
c     Number of radial points and number of columns written, incl. abscissae:
c      write(75,*)ist11, 3*ncomp+3
c      do i = 1, ist11
c      write(75,2000)absr(i), (yout(i,j),j=1,3*ncomp+2)
c      end do
      write(75,*)ist11, 3*ncomp+3+9+1
c2004      write(75,*)ist11, 3*ncomp+3+9
      do i = 1, ist11
c2004      write(75,2000)absr(i), (yout(i,j),j=1,3*ncomp+2), (youth(i,j),j=1,9)
      write(75,2000)absr(i), (yout(i,j),j=1,3*ncomp+2), (youth(i,j),j=1,9), abscis(i)
      end do
 
c      fac4 = 1.d0 / twopi ** 2
c      if(.not.cyl)fac4 = fac4 * r0i

      poyjum(0) = dcmplx(yout(1,trncp1), yout(1,trncp2))
      ist = 1
        do ireg = 1, nreg
        i = ist
        dpoydr(i) = (yout(i+1,trncp1) - yout(i,trncp1)) / (absout(i+1) * (absout(i+1) - absout(i)))
c       (First absout(i+1) in above line to avoid trouble at axis)
          do i = ist + 1, ist + istp(ireg) - 2
CCC         dFlux/(rho*drho *(2pi)**2 *RA): compare to POWDEN in OUTPOW
c  now    dFlux/(rho*drho): compare to POWDEN in OUTPOW (new scaling)
          dpoydr(i) = (yout(i+1,trncp1) - yout(i-1,trncp1)) / (absout(i) * (absout(i+1) - absout(i-1)))
          end do
C       Set IST to 1st point of next region:
        ist = ist + istp(ireg)
        i = ist - 1
        dpoydr(i) = (yout(i,trncp1) - yout(i-1,trncp1)) / (absout(i) * (absout(i) - absout(i-1)))
c       Jumps of complex inward Poynting flux at region boundaries:
        poyjum(ireg) = dcmplx(yout(i,trncp1), yout(i,trncp2))
        if(ireg.lt.nreg)poyjum(ireg) = poyjum(ireg) - dcmplx(yout(ist,trncp1), yout(ist,trncp2))
        end do

 304  continue
c     ========

      write(nofile,*)'Input complex power based on internal boundaries (no account of feeders!): ', totpoy
      write(nofile,*)'Inward Poynting jumps at region boundaries:'
        do i = 0, nreg
        write(nofile,*) i, poyjum(i)
        end do
      WRITE(NOFILE,*)'Electric field modes L2 norm:'
      WRITE(75,*)'Electric field modes L2 norm'
      WRITE(75,*)nmode(imax), 2
      WRITE(NOFILE,*)'  M       NORM'
        DO IM = 1, NMODE(IMAX)
        MONORM(IM) = DSQRT(MONORM(IM) / (RX0M(NREG)-RX0M(0)))
        WRITE(NOFILE,3000) IM-1+MINF(IMAX), MONORM(IM)
 3000   format(1h , i4, 2x, g13.5)
        WRITE(75,*) float(IM-1+MINF(IMAX)), MONORM(IM)
          DO I = 1, IST11
          TABLO(I,IM) = DSQRT(TABLO(I,IM))
          END DO
        END DO
      WRITE(79,*)'Electric field modes norm'
      WRITE(79,*)ist11, nmode(imax)
      write(79,2001)(absout(i), i = 1, ist11)
      write(79,2001)(minf(imax)+im-1, im = 1, nmode(imax))
        do i = 1, ist11
        write(79,2001)(tablo(i,im), im = 1, nmode(imax))
        end do

cERN  cccccccccccc Write E-field poloidal harmonics ccccccccccccccccc
	open (UNIT = 701, FILE = paplda // 'Efield_pol_harmonics.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
	write(701,*),"Electric field poloidal harmonics (L2 norm)"
        do i = 1, nmode(imax)
        write(701,"(i5,g16.6)") i-1+minf(imax), monorm(i)
        end do
	close(701)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Detect zero components: 
      write(nofile,*)'outez:'
        if(plosol)then
          do i = 1, trncp2
          outez(i) = rveno(yout(1,i), ist11, 2) .eq. 0.d0
          end do
        outez(3*nficom+3) = rveno( monorm, nmode(imax), 2 ) .eq. 0.d0
        write(nofile,*)(outez(i), i = 1, trncp2)
        end if
 
      return 

 182  format(1h ,4(i4,', (',g11.3,';',g11.3,') '))
 500  format(1h ,4('(',g23.15,',',g23.15,') '))
 2000 format(1h ,25(2x,g13.5))
c2004 2000 format(1h ,21(2x,g13.5))
 2001 format(1h ,10(1x,g12.4))
 1000 format(1h ,'plot #',i4,'  ',a50)
      end

