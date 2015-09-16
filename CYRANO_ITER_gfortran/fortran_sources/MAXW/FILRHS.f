      subroutine filrhs(atabou)

      implicit none

      logical atabou

c     Builds righ-hand side for current finite element
c     atabou: switch, .true. at a boundary, .false. in volume

c     N.B.: it is assumed that the first element (axis, crown) is free
c     from sources.

      include 'pardim.copy'
      include 'comusr.f'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comant.copy'
c include 'comfou.copy'
      include 'comfic.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comfin.copy'
      include 'commod.copy'
      include 'comequ.copy'
      include 'comswe.copy'
      include 'comin2.copy'
      include 'comphy.copy'
c
      logical inl, inr

      integer 
     ;  ja, irea, imod, isc, irea1, isuba, icolo, iela
     ;, nmodc, ncstrt, i, im, j, ideca, ifirse, ilaste, ibulo
     ;, ielma, mrm, mr2, idecam, ideca2
     ;, nsum
c     ;, it, index(12)

      double precision 
     ;  rle
     ;, jaint1(maxgau), jaint2(maxgau), ntint1(maxgau), ntint2(maxgau)
     ;, gint1(maxgau), gint2(maxgau)
     ;, jon1(maxgau), jon2(maxgau), gon1(maxgau), gon2(maxgau)
     ;, coi1(maxgau), coi2(maxgau), sii1(maxgau), sii2(maxgau) 
     ;, phibarae(maxgau,2), tth1, tth2
c     ;, tr1(maxgau), tr2(maxgau)
c     ;, thba, theab1, theab2, dtheab

      complex*16 
     ;  fac, fac0, facv, facv2, fac2, facju, far, fat
     ;, vec(6), tra, trac, trho, tthe, tpl, tmi, tpa, fapl, fami, fapa, eminpbmp(maxgau,2)
     ;, cdeid
c     ;, t

      external nsum, cdeid

c     Gives power units to system quadratic form:
      fac2 = ci * glofac * omegag * mu0 * rnorm * twopi ** 2

c     Loop over independent sources: antennae and Faraday shield modes
      do 1 ja = 1, neffan + nscree * nmoscr

c       Identify the corresponding region boundary: 
        if(ja .le. neffan)then
        irea = irbant(ja)
        else
        imod = mod(ja - neffan, nmoscr)
        isc = (ja - neffan - imod) / nmoscr
        if(imod.ne.0)isc = isc + 1
        irea = irbscr(isc)
        end if
      irea1 = irea + 1
      isuba = nsum(irea,ns)
      icolo = iconn(isuba)
      iela = ilael(irea)
      nmodc = nmode(iela)

c     
      if(atabou .and. irea.eq.ireg .and. nmodc.gt.0)then
c     ..................................................
c     We are at a region boundary containing the current antenna...
      ncstrt = nmodc * ncstr(irea)

      y = rx0(irea)
      fac = fac2 * y
c      intab = iela * (ngauss + 1) + irea
      intab = ifiabs(ifiel(irea1)) - 1

        if(circ)then
c       ~~~~~~~~~~~~
        si = sitab(intab)
        co = cotab(intab)

          do i = 1, nmoant
          mr = moant(i) + 1 - minf(iela)

c         Right part of last element in region, = left part of current (boundary) block:
          ideca = (mr - 1) * icolo
          write(nofile,1000)moant(i), cura(i,ja), rcura(i,ja), rcurat(i,ja), chaa(i,ja) / eps0

c         Source terms from surface term in weak form:
            if(styp(isuba).eq.'HEC')then
c           (R/R0)*j+,m:
            vec(1) = fac * ci * sqrt2i * (co * rcura(i,ja) - si * rcurat(i,ja))
            vec(2) = czero
c           (R/R0)*j-,m:
            vec(3) = - vec(1)
            vec(4) = czero
c           (R/R0)*j//,m:
            vec(5) = fac * (si * rcura(i,ja) + co * rcurat(i,ja))
            vec(6) = czero

            else if(styp(isuba).eq.'M23')then
            vec(1) = czero
c           (R/R0)*jtheta,m:
            vec(2) = fac * rcura(i,ja)
            vec(3) = czero
c           (R/R0)*jphi,m:
            vec(4) = fac * rcurat(i,ja)
            vec(5) = czero

            else if(styp(isuba).eq.'CAR')then
            print *, 'CAR element at antenna: RHS contribution to write!'
	      stop
            end if

            do j = 1, icolo
            bel(ideca+j,ja) = bel(ideca+j,ja) + vec(j)
            end do
          end do

        else  ! Non-circular case
c       ~~~~
          do i = 1, nmoant
          mr = moant(i) + 1 - minf(iela)
          ideca = (mr - 1) * icolo
          write(nofile,1000)moant(i), cura(i,ja), rcura(i,ja), rcurat(i,ja), chaa(i,ja)/eps0

c if toroidal antenna current: to do
          
            if(styp(isuba).eq.'HEC')then
            vec(1) = fac * ci * sqrt2i * rcurac(i,ja)
            vec(2) = czero
            vec(3) = - vec(1)
            vec(4) = czero
            vec(5) = fac * rcuras(i,ja)
            vec(6) = czero
            else if(styp(isuba).eq.'M23')then
            vec(1) = czero
            vec(2) = fac * rcura(i,ja)
            vec(3) = czero
            vec(4) = fac * rcurat(i,ja)
            vec(5) = czero
            end if

            do j = 1, icolo
            bel(ideca+j,ja) = bel(ideca+j,ja) + vec(j)
            end do
          end do
        end if
c       ~~~~~~

c       Source terms in jump conditions between regions:
c       -----------------------------------------------
        facju = dcmplx(0.d0, - omegag * mu0 * rnorm)
c       (facju is -1/(2*fac) of old notes)
c       JUANT: charge/eps0, 0, jtheta, 0, Rjphi
c       @check: shouldn't it be Rjphi/R0?
c     Antenna: 5 constraints, normalized as follows in INTBOU:
c     [Erho + ...(plasma)] = ...                    (charge / epsilon0)
c     [Etheta] = ...                                (0)
c     [dEtheta/dy -dErho/(y dtheta)] = ...          (-i omega mu0 rnorm jtheta)
c     [Ephi] = ...                                  (0)
c     (R/R0)*[dEphi/dy -rnorm*dErho/(R dphi)] = ... (-i omega mu0 rnorm Rjphi/R0)

        do mr = 1, nmodc
        ideca = nmode(iela) * icolo + (mr - 1) * ncstr(irea)
c       Nonzero contributions; note different factor on radial component:
        bel(ideca+1,ja) = bel(ideca+1,ja) + juant(1,mr,ja)
        bel(ideca+3,ja) = bel(ideca+3,ja) + facju * juant(3,mr,ja)
        bel(ideca+5,ja) = bel(ideca+5,ja) + facju * juant(5,mr,ja)
c          do j2 = 1, ncstr(irea)
c          bel(ideca+j2,ja) = bel(ideca+j2,ja) + facju * juant(j2,mr,ja)
c          end do
        end do

      end if
c     ......
  1   continue


c     Feeders:
c     ========
c     For all antennas and eventual Faraday shield modes,...
      do 2 ja = 1, neffan + nscree * nmoscr
        if(ja .le. neffan)then
c       ja is antenna index
        irea = irbant(ja)
        else
c       isc is faraday shield index
        imod = mod(ja - neffan, nmoscr)
        isc = (ja - neffan - imod) / nmoscr
        if(imod.ne.0)isc = isc + 1
        irea = irbscr(isc)
        end if
      irea1 = irea + 1
      iela = ilael(irea)
      nmodc = nmode(iela)

c     ... find whether there is a feeder (volume) contribution to the current element:
      if(.not.atabou .and. ja.le.neffan .and. feeder(ja) .and. iel.ge.ifiel(irea1))then
c     (n.b. when FALEN is true, remember the different current normaliz. according to antenna end type.)
      irea = ireg
      isuba = isubr
      iela = iel

      ifirse = ifiel(irea)
      ilaste = ilael(irea)
      icolo = iconn(isuba)
      ibulo = ibub(isuba)
      rle = fl(iela)
      ielma = iela
      if(nmode(iela-1).gt.nmode(iela))ielma = iela - 1
        if(nmode(ielma).gt.0)then
c 5/12/03: following replaced by construction of 
c          theab1(), theab2(), thba(), dtheab() in filant;
c          (idem for Faraday shields).
c          if(cokpco)then
c         approx. Thebars:
c          intab = irea + (ngauss + 1) * nsum(isuba,iele)
c          theab1 = ckt(intab, ipolin(ja), 3)
c          theab2 = ckt(intab, ipolou(ja), 3)
c          thba = 0.5d0 * (theab1 + theab2)
c          dtheab = theab2 - theab1
c          else
c          theab1 = thea1(ja)
c          theab2 = thea2(ja)
c          thba = theaa(ja)
c          dtheab = dthea(ja)
c          end if
c        if(theab1(ja) .lt. 0.d0)theab1(ja) = theab1(ja) + twopi
c        if(theab2(ja) .lt. 0.d0)theab2(ja) = theab2(ja) + twopi

c          if(.not.circ)then
c approx: these factors vary inside each element:
c          tr1 = eqt(intab,ipolin(ja),9) / eqt(intab,ipolin(ja),7)
c          tr2 = eqt(intab,ipolou(ja),9) / eqt(intab,ipolou(ja),7)
c          end if
          
        fac0 = ci * glofac * omegag * mu0 * rnorm * rle
c     (toroidal delta function antenna, flat spectrum tofac(ja) = 1)
c       ; * tofac(ja)
        
          if(circ)then
c         ~~~~~~~~~~~~
          if(.not.cokpco)then  ! Standard coordinates:
c         ...................
          do im = 1, nmoant
          m = moant(im)
          mr = m + 1 - minf(iela-1)
          mrm = m + 1 - minf(ielma)
          mr2 = m + 1 - minf(iela)
          ideca = (mr - 1) * icolo
          idecam = nmode(iela-1) * icolo + (mrm - 1) * ibulo
          ideca2 = nmode(iela-1) * icolo + nmode(ielma)*ibulo + (mr2-1)*icolo
    
          inl = nmode(iela-1).gt.0 .and. m.ge.minf(iela-1) .and. m.le.msup(iela-1)
          inr = nmode(iela).gt.0 .and. m.ge.minf(iela) .and. m.le.msup(iela)
    
          facv = fac0 * cdeid(- m * thba(ja))

            if(.not.falen(ja))then  ! Case no propagation effects:
            facv = facv * (- 2.d0 * ci * dsin(m * dtheab(ja) * 0.5d0))

            else  ! With propagation effects:
c           N.B. current density variation along feeders is canceled by volume element.
            tra = cdeid(m * dtheab(ja) * 0.5d0)
              if(anetyp(ja).eq.'SHC')then
c             outlet and inlet contributions:
              facv = facv * (dconjg(tra) - cdcos(betal(ja)) * tra)
              else if(anetyp(ja).eq.'OPC')then
c             inlet contribution only:
              facv = facv * (- tra)
              end if
            end if

            if(styp(isuba).eq.'M23')then
            if(inl)
     ;      bel(ideca+1,ja) = bel(ideca+1,ja) + facv / 6.d0
            bel(idecam+1,ja) = bel(idecam+1,ja) + facv * 2.d0 / 3.d0
            if(inr)
     ;      bel(ideca2+1,ja) = bel(ideca2+1,ja) + facv / 6.d0

            else if(styp(isuba).eq.'HEC')then
            facv2 = facv * sqrt2i
              if(inl)then
              bel(ideca+1,ja) = bel(ideca+1,ja) + facv2 / 2.d0
              bel(ideca+2,ja) = bel(ideca+2,ja) + facv2 / 12.d0 * rle
              bel(ideca+3,ja) = bel(ideca+3,ja) + facv2 / 2.d0
              bel(ideca+4,ja) = bel(ideca+4,ja) + facv2 / 12.d0 * rle
              end if
              if(inr)then
              bel(ideca2+1,ja) = bel(ideca2+1,ja) + facv2 / 2.d0
              bel(ideca2+2,ja) = bel(ideca2+2,ja) - facv2 / 12.d0 * rle
              bel(ideca2+3,ja) = bel(ideca2+3,ja) + facv2 / 2.d0
              bel(ideca2+4,ja) = bel(ideca2+4,ja) - facv2 / 12.d0 * rle
              end if
            end if
          end do

          else     ! Constant k// coordinates:
c         ....
c         exp[-in(phibar-phi)] at Gauss points on each feeder:
          intab = ifiabs(iel)
	    tth1 = thea1(ja)
          if(tth1.lt.0.d0)tth1 = tth1 + twopi
          tth2 = thea2(ja)
          if(tth2.lt.0.d0)tth2 = tth2 + twopi
            do ig = 1, ngauss
            intab = intab + 1
cPL25/1/05 more accurate:
            call interp1(polang, polang(npfft+1), ckt(intab,1,4), nabplo, npfft+1, tth1, phibarae(ig,1), 1, 'R')
	      eminpbmp(ig,1) = cdeid(-n*(phibarae(ig,1) - phibarac(ja)))
            call interp1(polang, polang(npfft+1), ckt(intab,1,4), nabplo, npfft+1, tth2, phibarae(ig,2), 1, 'R')
	      eminpbmp(ig,2) = cdeid(-n*(phibarae(ig,2) - phibarac(ja)))
c	      eminpbmp(ig,1) = cdeid(-n*(ckt(intab,ipolin(ja),4) - phibarac(ja)))
c	      eminpbmp(ig,2) = cdeid(-n*(ckt(intab,ipolou(ja),4) - phibarac(ja)))
	      end do

          do im = 1, nmoant
          m = moant(im)
          mr = m + 1 - minf(iela-1)
          mrm = m + 1 - minf(ielma)
          mr2 = m + 1 - minf(iela)
          ideca = (mr - 1) * icolo
          idecam = nmode(iela-1) * icolo + (mrm - 1) * ibulo
          ideca2 = nmode(iela-1) * icolo + nmode(ielma)*ibulo + (mr2-1)*icolo
    
          inl = nmode(iela-1).gt.0 .and. m.ge.minf(iela-1) .and. m.le.msup(iela-1)
          inr = nmode(iela).gt.0 .and. m.ge.minf(iela) .and. m.le.msup(iela)
    
          facv = fac0 * cdeid(- m * thba(ja))
          tra = cdeid(m * dtheab(ja) * 0.5d0)
          trac = dconjg(tra)

           if(falen(ja))then
           write(nofile,*)'FILRHS circ: propagation effects along antenna not implemented yet with cokpco, are ignored'
           end if

          intab = ifiabs(iel)
            do ig = 1, ngauss
            intab = intab + 1
              if(styp(isuba).eq.'M23')then
c             Weight for radial feeder current component:
              trho = wga(ig) * (trac * eminpbmp(ig,2) - tra * eminpbmp(ig,1))
              far = facv * trho
c             M23 basis function type is 2 for radial component:
                if(inl)then
                bel(ideca+1,ja) = bel(ideca+1,ja) + far * bafgn(1,0,2,ig)
	          end if
              bel(idecam+1,ja) = bel(idecam+1,ja) + far * bafgn(2,0,2,ig)
                if(inr)then
                bel(ideca2+1,ja) = bel(ideca2+1,ja) + far * bafgn(3,0,2,ig)
                end if

c              Test integral of basis functions same as analytical value used above in std coordinates:
c              write(nofile,*)'cokpco test integral basis functions: must print three 1'
c              write(nofile,*)sum(wga(1:ngauss)*bafgn(1,0,2,1:ngauss))*6.d0
c              write(nofile,*)sum(wga(1:ngauss)*bafgn(2,0,2,1:ngauss))*1.5d0
c              write(nofile,*)sum(wga(1:ngauss)*bafgn(3,0,2,1:ngauss))*6.d0
c              The test was successfully passed!

              else if(styp(isuba).eq.'HEC')then	        
	        tpl =  wga(ig) * sqrt2i * (trac * eminpbmp(ig,2) - tra * eminpbmp(ig,1))
              fapl = facv * tpl
                if(inl)then
                bel(ideca+1,ja) = bel(ideca+1,ja) + fapl * bafgn(1,0,1,ig)
                bel(ideca+2,ja) = bel(ideca+2,ja) + fapl * bafgn(2,0,1,ig) * rle
                bel(ideca+3,ja) = bel(ideca+3,ja) + fapl * bafgn(1,0,1,ig)
                bel(ideca+4,ja) = bel(ideca+4,ja) + fapl * bafgn(2,0,1,ig) * rle
                end if
                if(inr)then
                bel(ideca2+1,ja) = bel(ideca2+1,ja) + fapl * bafgn(3,0,1,ig)
                bel(ideca2+2,ja) = bel(ideca2+2,ja) + fapl * bafgn(4,0,1,ig) * rle
                bel(ideca2+3,ja) = bel(ideca2+3,ja) + fapl * bafgn(3,0,1,ig)
                bel(ideca2+4,ja) = bel(ideca2+4,ja) + fapl * bafgn(4,0,1,ig) * rle
                end if
              end if
            end do
          end do
          end if
c         ......

          else  ! Non-circular case
c         ~~~~
          if(cokpco)write(6,*)'FILRHS: cokpco case not dealt with yet in non-circular... Expect anything!'

c         Interpolate equilibrium to 2 feeders poloidal positions:
c          intab = ireg + (ngauss + 1) * (iel - 1)
          intab = ifiabs(iel)
            do ig = 1, ngauss
            intab = intab + 1
c           Jn at inlet:
            call interp1(polang, polang(npfft+1), eqt(intab,1,9), nabplo, npfft+1, thea1(ja), jaint1(ig), 1, 'R')
c     ;    , npfft+1, theab1, jaint1(ig), 1, 'R')
c           Ntn at inlet:
            call interp1(polang, polang(npfft+1), eqt(intab,1,7), nabplo, npfft+1, thea1(ja), ntint1(ig), 1, 'R')
c     ;    , npfft+1, theab1, ntint1(ig), 1, 'R')
c           Gn at inlet:
            call interp1(polang, polang(npfft+1), eqt(intab,1,8), nabplo, npfft+1, thea1(ja), gint1(ig), 1, 'R')
c           cos Theta at inlet:
            call interp1(polang, polang(npfft+1), eqt(intab,1,14), nabplo, npfft+1, thea1(ja), coi1(ig), 1, 'R')
c           sin Theta at inlet:
            call interp1(polang, polang(npfft+1), eqt(intab,1,15), nabplo, npfft+1, thea1(ja), sii1(ig), 1, 'R')
            jon1(ig) = jaint1(ig) / ntint1(ig)
            gon1(ig) = gint1(ig) / ntint1(ig)
            end do
            if(udsant(ja))then
            call dcopy(ngauss, jaint1, 1, jaint2, 1)
            call dcopy(ngauss, ntint1, 1, ntint2, 1)
            call dcopy(ngauss, jon1, 1, jon2, 1)
            call dcopy(ngauss, coi1, 1, coi2, 1)
            call dcopy(ngauss, sii1, 1, sii2, 1)
c             Gn is antisymmetric:
	        do ig = 1, ngauss
	        gint2(ig) = - gint1(ig)
	        gon2(ig) = - gon1(ig)
	        end do
            else
c            intab = ireg + (ngauss + 1) * (iel - 1)
            intab = ifiabs(iel)
              do ig = 1, ngauss
              intab = intab + 1
c             Jn at outlet:
              call interp1(polang, polang(npfft+1), eqt(intab,1,9), nabplo, npfft+1, thea2(ja), jaint2(ig), 1, 'L')
c     ;      , npfft+1, theab2, jaint2(ig), 1, 'L')
c             Ntn at outlet:
              call interp1(polang, polang(npfft+1), eqt(intab,1,7), nabplo, npfft+1, thea2(ja), ntint2(ig), 1, 'L')
c     ;      , npfft+1, theab2, ntint2(ig), 1, 'L')
c             Gn at outlet:
              call interp1(polang, polang(npfft+1), eqt(intab,1,8), nabplo, npfft+1, thea2(ja), gint2(ig), 1, 'L')
c             cos Theta at outlet:
              call interp1(polang, polang(npfft+1), eqt(intab,1,14), nabplo, npfft+1, thea2(ja), coi2(ig), 1, 'R')
c             sin Theta at outlet:
              call interp1(polang, polang(npfft+1), eqt(intab,1,15), nabplo, npfft+1, thea2(ja), sii2(ig), 1, 'R')
              jon2(ig) = jaint2(ig) / ntint2(ig)
              gon2(ig) = gint2(ig) / ntint2(ig)
              end do
            end if
          do im = 1, nmoant
          m = moant(im)
          mr = m + 1 - minf(iela-1)
          mrm = m + 1 - minf(ielma)
          mr2 = m + 1 - minf(iela)
          ideca = (mr - 1) * icolo
          idecam = nmode(iela-1) * icolo + (mrm - 1) * ibulo
          ideca2 = nmode(iela-1) * icolo + nmode(ielma)*ibulo + (mr2-1)*icolo
    
          inl = nmode(iela-1).gt.0 .and. m.ge.minf(iela-1) .and. m.le.msup(iela-1)
          inr = nmode(iela).gt.0 .and. m.ge.minf(iela) .and. m.le.msup(iela)
    
          facv = fac0 * cdeid(- m * thba(ja))
          tra = cdeid(m * dtheab(ja) * 0.5d0)
          trac = dconjg(tra)

           if(falen(ja))then
           write(nofile,*)'Propagation effects along antenna not implemented yet in general geometry, are ignored'
           end if

            do ig = 1, ngauss

              if(styp(isuba).eq.'M23')then
c             Weight for radial feeder current component:
              trho = wga(ig) * (trac * jon2(ig) - tra * jon1(ig))
c Old experiment; why is this better?:
ccc            t = wga(ig) * (- trac * jon2(ig) + tra * jon1(ig))
c             Weight for poloidal feeder current component
c             NB present in noncircular geometry!
              tthe = wga(ig) * (trac * gon2(ig) - tra * gon1(ig))
              far = facv * trho
              fat = facv * tthe
c             M23 basis function type is 2 for radial and 1 for poloidal component!
                if(inl)then
                bel(ideca+1,ja) = bel(ideca+1,ja) + far * bafgn(1,0,2,ig)
                bel(ideca+2,ja) = bel(ideca+2,ja) + fat * bafgn(1,0,1,ig)
                bel(ideca+3,ja) = bel(ideca+3,ja) + fat * bafgn(2,0,1,ig) * rle
	          end if
              bel(idecam+1,ja) = bel(idecam+1,ja) + far * bafgn(2,0,2,ig)
                if(inr)then
                bel(ideca2+1,ja) = bel(ideca2+1,ja) + far * bafgn(3,0,2,ig)
                bel(ideca2+2,ja) = bel(ideca2+2,ja) + fat * bafgn(3,0,1,ig)
                bel(ideca2+3,ja) = bel(ideca2+3,ja) + fat * bafgn(4,0,1,ig) * rle
                end if
              else if(styp(isuba).eq.'HEC')then	        
	        tpl =  wga(ig) * sqrt2i * (trac * (jon2(ig) + ci * coi2(ig)*gon2(ig)) - tra * (jon1(ig) + ci * coi1(ig)*gon1(ig)))
	        tmi =  wga(ig) * sqrt2i * (trac * (jon2(ig) - ci * coi2(ig)*gon2(ig)) - tra * (jon1(ig) - ci * coi1(ig)*gon1(ig)))
	        tpa =  wga(ig) * (trac * sii2(ig) * gon2(ig) - tra * sii1(ig) * gon1(ig))
              fapl = facv * tpl
              fami = facv * tmi
              fapa = facv * tpa
                if(inl)then
                bel(ideca+1,ja) = bel(ideca+1,ja) + fapl * bafgn(1,0,1,ig)
                bel(ideca+2,ja) = bel(ideca+2,ja) + fapl * bafgn(2,0,1,ig) * rle
                bel(ideca+3,ja) = bel(ideca+3,ja) + fami * bafgn(1,0,1,ig)
                bel(ideca+4,ja) = bel(ideca+4,ja) + fami * bafgn(2,0,1,ig) * rle
                bel(ideca+5,ja) = bel(ideca+5,ja) + fapa * bafgn(1,0,1,ig)
                bel(ideca+6,ja) = bel(ideca+6,ja) + fapa * bafgn(2,0,1,ig) * rle
                end if
                if(inr)then
                bel(ideca2+1,ja) = bel(ideca2+1,ja) + fapl * bafgn(3,0,1,ig)
                bel(ideca2+2,ja) = bel(ideca2+2,ja) + fapl * bafgn(4,0,1,ig) * rle
                bel(ideca2+3,ja) = bel(ideca2+3,ja) + fami * bafgn(3,0,1,ig)
                bel(ideca2+4,ja) = bel(ideca2+4,ja) + fami * bafgn(4,0,1,ig) * rle
                bel(ideca2+5,ja) = bel(ideca2+5,ja) + fapa * bafgn(3,0,1,ig)
                bel(ideca2+6,ja) = bel(ideca2+6,ja) + fapa * bafgn(4,0,1,ig) * rle
                end if
              end if
            end do
          end do
          end if
c         ~~~~~~
    
        end if
      end if
   2  continue

      return

 1000 format(1h , i4, 4(1h(, g15.7, 1h,, g15.7, 1h), 2x))
      end
