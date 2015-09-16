      subroutine filant

ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)
      implicit none

c     Fills modal excitation vector according to antenna poloidal configuration.
c     (During the main calculation, the antennae are assumed infinitely thin
c      toroidally. The toroidal spectral factors are applied during processing
c      of the solution.)

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'comant.copy'
      include 'compla.copy'
      include 'commag.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comfin.copy'
      include 'comswe.copy'
      include 'comfou.copy'
      include 'commod.copy'
      include 'comfic.copy'
      include 'comin2.copy'
      include 'comphy.copy'
      include 'comfft.copy'

      character*3 eletyp

      integer 
     ;  i, im, impk, is, iexc, i1, j, ja, jant, jm, k1, k2, km, kp, mml, mmk
     ;, isrchfge, isrchfgt, nsum
     ;, ntf, itf
     ;, iregb
      
      double precision 
     ;  xtra, pi2i, dnfi
     ;, tth1, tth2, jonin
     ;, jaint1, jaint2, ntint1, ntint2, gint1, gint2, tr1, tr2, fa
cPL 25/8/04:
     ;, vsthbar(npfft), store(npfft+1), rfo(-npft2:npfft), ifo(-npft2:npfft), epsa, keps
     ;, rcs(maxpom), ics(maxpom)
     ;, dsinc
c 5/12/03: now in common comant:
c     ;, thba, theab1, theab2, dtheab
      
      complex*16 
     ;  torfac, polfac, polfa2, polfa3, polfa4, polfa5, polcur, polcha
     ;, ccsc
     ;, cxtra, csxox, xtra2, tra
     ;, foutra(-npft2:npfft,4), ft(maxpom), einpbac
     ;, cdeid, cdsinc
     
      external isrchfge, isrchfgt, nsum, dsinc, cdeid, cdsinc

      write(nofile,*)'enter filant'
      pi2i = 0.5d0 * oopi

c       Loop over the antennae:
        do i = 1, ntoant
        iregb = irbant(i)
        rhoant = rx0m(iregb)
c NB: correct use of inputs thea1, thea2 requires no 2*pi angle jump on antenna!!!
c Idem on eventual screens
          if(udsant(i))then
          theaa(i) = 0.d0
          else
          theaa(i) = (thea1(i) + thea2(i)) * 0.5d0
          end if
        dthea(i) = thea2(i) - thea1(i)
          if(cyl)then
          dphia(i) = 0.d0
          else
c         Antenna toroidal position is given by distance in theta=0 plane,
c         not by angle:
c         Radial index at antenna:
c 5/12/03 Using ifiabs table: (-1 to get same point as before. 
c         Remember there are 2 table points at each region boundary)
          intab = ifiabs(ifiel(iregb+1)) - 1
c          intab = (ngauss + 1) * ilael(iregb) + iregb
          phiaa(i) = zaa(i) / eqt(intab,1,1)
          dphia(i) = dza(i) / eqt(intab,1,1)
          end if

          if(falen(i))then
c         A priori values for antenna propagation parameters:

C         Characteristic impedance:
          zcapri(i) = sqrt(laprio(i) / caprio(i))

C         Electrical length per unit poloidal angle:
C         NB: assumes poloidal antenna!
c         For circular case: (else must find antenna poloidal length first)
c         Small loss approximation to propagation constant (times radius):
          deltaa(i) = rhoant * (omegag * dsqrt(laprio(i)*caprio(i)) - ci * 0.5 * raprio(i) / zcapri(i) )
c         Electrical length (complex, radian):
          betal(i) = deltaa(i) * dthea(i)
          tra = cdsin(betal(i)) / cdcos(betal(i))
c           Single-strap input impedance: (engineers'convention)
      	if(anetyp(i).eq.'SHC')then
      	zap(i) = ci * zcapri(i) * tra
      	else if(anetyp(i).eq.'OPC')then
      	zap(i) = - ci * zcapri(i) / tra
      	end if
          if(.not.circ)write(nofile,*)'general geometry and cokpco: case FALEN=.T. remains to deal with!'
          end if

        end do  ! antenna loop

        if(monomo)then
        nmoant = 1
        modva1 = moant(1)
        modva2 = moant(1)
        write(nofile,*)'mode m=',moant(1),' ,kphi=',kphi

        else
        nmoant = modva2 - modva1 + 1
	    if (nmoant .gt. maxpom) then
	    write(nofile,*)'FILANT: too many poloidal modes! I stop. nmoant=', nmoant, 'Recompile with MAXPOM >=', nmoant
	    end if
          do i = 1, nmoant
          moant(i) = modva1 + i - 1
          end do
        write(nofile,*)'Mode set: m=',modva1,',',modva2,' kphi=',kphi
        end if

      call zset(7*maxpom*maxexc, czero, curare, 1)

      do 2 j = 1, neffan
c     ==================
      iregb = irbant(j)
      y = rx0(iregb)
      rhoant = rx0m(iregb)
      isubr = nsum(iregb, ns)
      eletyp = styp(isubr)
c      intab = iregb + (ngauss + 1) * ilael(iregb)
c 5/12/03 
      intab = ifiabs(ifiel(iregb+1)) - 1
 
c     Poloidal location of strap inlet and outlet expressed on [0, 2*pi[ for interpolation in poloidal tables
c     ... or for numerical approximation of Morse pulse.
      write(nofile,*)'FILANT: assuming low field side antenna'
      tth1 = thea1(j)
      if(tth1.lt.0.d0)tth1 = tth1 + twopi
c     (tth1 strictly < 2 pi)
      ipolin(j) = isrchfgt(npfft+1, polang, 1, tth1) - 1
c      ipolin(j) = isrchfgt(npfft+1, polang, 1, tth1)
c     (ipolin > 1 and <= npfft+1)
      tth2 = thea2(j)
      if(tth2.lt.0.d0)tth2 = tth2 + twopi
      ipolou(j) = isrchfgt(npfft+1, polang, 1, tth2)
c      ipolou(j) = isrchfgt(npfft+1, polang, 1, tth2) - 1
     
cPL 10/12/04: circular geometry with std or const.k// coord
        if(circ)then
c       ------------
        call fouco
c        write(nofile,*)'In filant; Fourier coeffs. a,b,c:'
c        write(nofile,*)(afou(i),i=0,ncrot)
c        write(nofile,*)(bfou(i),i=0,ncrot)
c        write(nofile,*)(cfou(i),i=0,ncrot)

          if(cokpco)then
          dnfi = 1.d0 / dfloat(npfft)
c         Phibar at antenna strap centre: switched off 19/1/05, result must be independent of phibarac!
c          call interp1(polang, polang(npfft+1), ckt(intab,1,4), nabplo, npfft+1, theaa(j), phibarac(j), 1, 'R')
          phibarac(j) = 0.d0
          einpbac = cdeid(n * phibarac(j))
c         Toroidal factor exp(-i*n*(phibar-phibar of strap centre)):
cPL27/1/05simplified using new table einphb:
          foutra(0:npfft,1) = dconjg(einphb(intab,1:npfft+1)) * einpbac * dnfi
c         This transform used in construction of poloidal current spectrum, used in jump conditions:
          foutra(0:npfft,2) = foutra(0:npfft,1) * r0orta(intab,1:npfft+1)

c          call cset(npfft, czero, foutra(0,1), 1)
cc         Incorporating antenna poloidal width, assuming only low field side antennas:
c            do i = 1, min(ipolin(j), ipolou(j))
c            foutra(i-1,1) = dconjg(einphb(intab,i)) * einpbac * dnfi
c            end do
c            do i = max(ipolin(j),ipolou(j)), npfft+1
c            foutra(i-1,1) = dconjg(einphb(intab,i)) * einpbac * dnfi
c            end do
c           To plot with array viewer:
	      rfo(0:npfft-1) = dreal(foutra(0:npfft-1,1))  
 	      ifo(0:npfft-1) = dimag(foutra(0:npfft-1,1))  
	      rfo(0:npfft-1) = dreal(foutra(0:npfft-1,2))  
 	      ifo(0:npfft-1) = dimag(foutra(0:npfft-1,2))  

	    call df2tcf(npfft, foutra(0,1), foutra(0,1), cwork2, cpy) ! -- IMSL -- exp[-im theta]
	    call df2tcf(npfft, foutra(0,2), foutra(0,2), cwork2, cpy) ! -- IMSL -- exp[-im theta]
            do i = 1, npft2
            foutra(-i,1) = foutra(npfft-i,1)
            foutra(-i,2) = foutra(npfft-i,2)
            end do
c         To plot with array viewer:
	    rfo(-npft2:npft2) = dreal(foutra(-npft2:npft2,1))  
	    ifo(-npft2:npft2) = dimag(foutra(-npft2:npft2,1))
	    rfo(-npft2:npft2) = dreal(foutra(-npft2:npft2,2))  
	    ifo(-npft2:npft2) = dimag(foutra(-npft2:npft2,2))

c         Estimate poloidal spectrum range for exp[-in*(phibar-phi)]: (just for information)
          epsa = rhoant * r0i
          keps = dsqrt((1.d0-epsa)/(1.d0+epsa))
	    q = qtab(intab)
	    write(nofile,*)'Estimated m limits of the exp[-in*(phibar-phi)] spectrum: ', -(1-keps)*n*q, ' to ', (1-keps)*n*q/keps
c  	    To complete when nq <0...
          write(6,*)'Done with circ cokpco filant'
          end if

        else
c       ----
c       Non-circular case.
c       Interpolate geometrical coeff. at strap inlet and outlet:
c 5/12/03: fixed a bug here (pol.angle at outlet was used for inlet).
c       Jn:
        call interp1(polang, polang(npfft+1), eqt(intab,1,9), nabplo, npfft+1, tth1, jaint1, 1, 'R')
c       Ntn:
        call interp1(polang, polang(npfft+1), eqt(intab,1,7), nabplo, npfft+1, tth1, ntint1, 1, 'R')
c       Gn:
        call interp1(polang, polang(npfft+1), eqt(intab,1,8), nabplo, npfft+1, tth1, gint1, 1, 'R')

        tr1 = jaint1 / ntint1

          if(udsant(j))then
          jaint2 = jaint1
          ntint2 = ntint1
          gint2 = - gint1
          tr2 = tr1
          else
c         Jn:
          call interp1(polang, polang(npfft+1), eqt(intab,1,9), nabplo, npfft+1, tth2, jaint2, 1, 'L')
c         Ntn:
          call interp1(polang, polang(npfft+1), eqt(intab,1,7), nabplo, npfft+1, tth2, ntint2, 1, 'L')
c         Gn:
          call interp1(polang, polang(npfft+1), eqt(intab,1,8), nabplo, npfft+1, tth2, gint2, 1, 'L')
          tr2 = jaint2 / ntint2
          end if

        call zset(4*(3*npft2+1), czero, foutra(-npft2,1), 1)
      
        dnfi = 1.d0 / dfloat(npfft)
c       Jn / Ntn at inlet:
        jonin = tr1
c       Various Fourier transforms, to be convoluted with Morse pulse later on:
          if(cyl)then
          ntf = 3
            do i = 1, npp
	      i1 = i - 1
c           For RCURA:
c           Ntn:
            foutra(i1,1) = jonin * eqt(intab,i,7) * dnfi
c           Ntn times cos and sin of magn. angle:
            foutra(i1,2) = foutra(i1,1) * eqt(intab,i,14)
            foutra(i1,3) = foutra(i1,1) * eqt(intab,i,15)
            end do
          else
          ntf = 4
            do i = 1, npp
	      i1 = i - 1
c           For RCURA:
c           checked this line 11/12/03
            foutra(i1,1) = jonin * eqt(intab,i,7) * dnfi
            foutra(i1,2) = foutra(i1,1) * eqt(intab,i,14)
            foutra(i1,3) = foutra(i1,1) * eqt(intab,i,15)
c           For CURA:
c           1/R:
            foutra(i1,4) = jonin / eqt(intab,i,1) * dnfi
            end do
          end if

c         NB program execution doesn't pass here any more in circular! THIS SECTION TO BE UPDATED
          if(cokpco)then
c         Phibar at antenna strap centre:
cPL disabled
c          call interp1(polang, polang(npfft+1), ckt(intab,1,4), nabplo, npfft+1, theaa(j), phibarac(j), 1, 'R')
          phibarac(j) = 0.d0
c         Interpolate tables to uniform thetabar mesh:
            do itf = 1, ntf
	      store(1:npp) = dreal(foutra(0:npp-1,itf))
            call interp2(ckt(intab,1,3), nabplo, store(1:npp), 1, npp, polang, vsthbar(1:npp), npp)
	      foutra(0:npp-1,itf) = dcmplx(store(1:npp),0.d0)
	      end do
c         Additional toroidal factor exp(-i*n*(phibar-phibar of strap centre)):
          store(1:npp) = ckt(intab,1:npp,4) - phibarac(j)
          call interp2(ckt(intab,1,3), nabplo, store(1:npp), 1, npp, polang, vsthbar(1:npp), npp)
            do itf = 1, ntf
cWRONG        foutra(0:npp-1,itf) = foutra(0:npp-1,itf) * cdeid(-n*vsthbar(1:npp))
              do i = 1, npp
	        i1 = i - 1
              foutra(i1,itf) = foutra(i1,itf) * cdeid(-n*(ckt(intab,i,4)-phibarac(j)))
              end do
            end do
          end if
c@
c     In khothb, store khi on equidistant Thetabar mesh: use contents of existing vector polang:
c     polang, khi and Thbar all range from 0 to 2pi; 
c     NB: Thbar was tabulated in ckt(,,3), khi in ckt(,,5).
c      call interp2(ckt(intab,1,3), nabplo, ckt(intab,1,5), nabplo, npp, polang, khothb, npp)

          if(updsym)then
c         Use up-down symmetry for the above quantities:
            do itf = 1, ntf
              do i = 1, npp-1
              foutra(npp+i-1,itf) = dconjg(foutra(npp-i-1,itf))
              end do
            end do
          end if
ccccccccccccccccccccc ERNESTO ccccccccccccccccccccccccccc
          do itf = 1, ntf
c          call cfft2(0, 1, npfft, foutra(0,itf), cworkp, foutra(0,itf))
	    call df2tcb(npfft, foutra(0,itf), foutra(0,itf), cwork2, cpy) ! -- IMSL --
            do i = 1, npft2
            foutra(-i,itf) = foutra(npfft-i,itf)
            end do
          end do
        end if
c       ------
c      print *, foutra(1:10,1)
c	print *, 'Stop at FILANT'
c	stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if(cokpco .and. .not.circ)then
c 5/12/03 Interpolate Thetabar at antenna strap inlet and outlet:
	  call interp1(polang, polang(npfft+1), ckt(intab,1,3), nabplo, npfft+1, tth1, theab1(j), 1, 'R')
	  call interp1(polang, polang(npfft+1), ckt(intab,1,3), nabplo, npfft+1, tth2, theab2(j), 1, 'R')
cERN	  BUG fixed here (theab1 -> theab1 - 2*pi)        
	  if(theab1(j) .gt. pi) theab1(j) = theab1(j) - 2*pi
cPL
	  if(theab2(j) .gt. pi) theab2(j) = theab2(j) - 2*pi
	  thba(j) = 0.5d0*(theab1(j) + theab2(j))
        dtheab(j) = theab2(j) - theab1(j)

        else
c       In circular geometry, thetabar = theta:
        theab1(j) = thea1(j)
        theab2(j) = thea2(j)
        thba(j) = theaa(j)
        dtheab(j) = dthea(j)
        end if

c     Delta-function antenna at phi=0: (actual antenna width and toroidal location are accounted for in output routines)
      torfac = pi2i

        do 6 i = 1, nmoant
c       ==================
        m = moant(i)

          if(.not.falen(j))then  ! Case of no propagation effects along straps
          polcha = czero
c         Fourier transform of 'Morse pulse'
          polcur = cdeid(- m * thba(j)) * dtheab(j) * pi2i * dsinc(m * dtheab(j) * 0.5d0)

          else                   ! Case of propagation effects along straps
c           general geometry: to do!
            if(.not.circ)then
            write(nofile,*)'FILANT: propagation effects along antenna not implemented in noncircular!'  
	      stop
            end if              	
c         falen, circular case:
          polfac = cdeid(- m * thba(j)) * dtheab(j) * pi2i
          polcha = polfac
          polcur = polfac
          cxtra = 0.5d0 * dtheab(j) * (m + deltaa(j))
	    csxox = cdsinc(cxtra)
c     	Remember BETAL is complex!
          xtra2 = cdexp(ci*0.5d0*betal(j))
          polfa2 = xtra2 * csxox
          cxtra = 0.5d0 * dthea(j) * (m - deltaa(j))
	    csxox = cdsinc(cxtra)
          polfa3 = csxox / xtra2

            if(anetyp(j).eq.'SHC')then
c           antenna in short-circuit; unit current is at short circuit:
            polcur = polcur * 0.5d0 * (polfa2 + polfa3)
            polcha = polcha * 0.5d0 * (polfa3 - polfa2) * ci
            else if(anetyp(j).eq.'OPC')then
            ccsc = cun / cdsin(betal(j))
c           insulated antenna end; unit current is at input:
            polcur =   polcur * (- ci * 0.5d0 * (polfa2 - polfa3) * ccsc)
            polcha = - polcha * 0.5d0 * (polfa2 + polfa3) * ccsc
            end if
          end if
          ft(i) = torfac * polcur

c 18/1/05: circular with constant k// coord. now treated inside circular option
          if(circ)then
c         ~~~~~~~~~~~~
          if(.not.cokpco)then            
          rcura(i,j) = ft(i)

            if(cyl)then
c           """""""""""
            cura(i,j) = ft(i)
            if(falen(j))chaa(i,j) = chaa(i,j) + torfac * polcha * (deltaa(j)/(ci * omegag * rhoant))
            else
c           """"
            polfac = czero
            polcha = czero
            k1 = max0(-ncrot, moant(1)-m )
            k2 = min0( ncrot, moant(nmoant)-m)
              do k = k1, k2
              mpk = m + k
              k2 = iabs(k)
              polfa2 = dthea(j) * cdeid(- mpk * theaa(j)) * bfou(k2) * pi2i
              polfa5 = polfa2
      	    if(.not.falen(j))then
      	    xtra = mpk * dthea(j) * 0.5d0
      	    polfa2 = polfa2 * dsinc(xtra)

      	    else
      	    cxtra = dthea(j) * 0.5d0 * (mpk + deltaa(j))
c     	    Remember BETAL is complex!
      	    xtra2 = cdexp(ci * betal(j) * 0.5d0)
      	    polfa4 = polfa4 * cdsinc(cxtra)
      	    polfa3 = cun / xtra2
      	    cxtra = dthea(j) * 0.5d0 * (mpk - deltaa(j))
      	    polfa3 = polfa3 * cdsinc(cxtra)	
      		if(anetyp(j).eq.'SHC')then
      		polfa5 = 0.5d0 * polfa5 * (polfa3 - polfa4) * ci
      		polfa2 = 0.5d0 * polfa2 * (polfa3 + polfa4)
      		else if(anetyp(j).eq.'OPC')then
                  ccsc = cun / cdsin(betal(j))
      		polfa5 = - polfa5 * (polfa4 + polfa3) * ccsc * 0.5d0
      		polfa2 = - polfa2 * (polfa4 - polfa3) * ccsc * 0.5d0 * ci
      		end if
      	    polcha = polcha + polfa5	
      	    end if
              polfac = polfac + polfa2
              end do
c           Resulting polfac is Fourier transform of R0/R times poloidal Morse pulse
            cura(i,j) = cura(i,j) + torfac * polfac * r0i
            if(falen(j))chaa(i,j) = chaa(i,j) + torfac * polcha * r0i * (deltaa(j)/(ci * omegag * rhoant))
            end if
c           """"""
c
          else ! cokpco circular
cPL25/1/05 replace following line by following convolution, more accurate (keeps 'exact' antenna in and out points):
c          rcura(i,j) = torfac * foutra(m,1)
          polfac = czero
          polfa2 = czero
            do k = -npft2, npft2
            mmk = m - k
	      cxtra = cdeid(- mmk * thba(j)) * dsinc(mmk * dtheab(j) * 0.5d0)
            polfac = polfac + foutra(k,1) * cxtra
            polfa2 = polfa2 + foutra(k,2) * cxtra
            end do
          rcura(i,j) = torfac * dtheab(j) * pi2i * polfac
          cura(i,j) = torfac * dtheab(j) * pi2i * polfa2 * r0i  ! To check!

cPL27/1/05 doesn't seem right; replaced by above.
cc         Get current density spectrum:
c          polfac = czero
c          k1 = max0(-ncrot, moant(1)-m )
c          k2 = min0( ncrot, moant(nmoant)-m)
c            do k = k1, k2
cc            mpk = m + k
c	      impk = i + k
ccPL26/1/05 replace following line, obsolete:
cc            polfac = polfac + foutra(mpk,1) * bfou(iabs(k))
c            polfac = polfac + rcura(impk,j) * bfou(iabs(k))
c            end do
c          cura(i,j) = cura(i,j) + polfac * r0i
	    end if

c19/1/05:
          else  ! Noncircular case, both std and constant k//
c         ~~~~
            if(cyl)then
            cura(i,j) = ft(i)

            else
            polfac = czero
            k1 = max0(-ncrot, moant(1)-m )
            k2 = min0( ncrot, moant(nmoant)-m)

              do k = k1, k2
              mpk = m + k
              k2 = iabs(k)
              polfa2 = cdeid(-mpk*thba(j)) * foutra(k,4) * dsinc(mpk * dtheab(j) * 0.5d0)
              polfac = polfac + polfa2
              end do
	      polfac = polfac * dtheab(j) * pi2i 
            cura(i,j) = torfac * polfac
            end if

          polfac = czero
          polfa3 = czero
          polfa4 = czero

            do k = k1, k2
            mpk = m + k
            k2 = iabs(k)
            polfa2 = cdeid(-mpk*thba(j)) * dsinc(mpk * dtheab(j) * 0.5d0)
c            polfa2 = dtheab(j) * cdeid(-mpk*thba(j)) * pi2i * dsinc(mpk * dtheab(j) * 0.5d0)
            polfac = polfac + polfa2 * foutra(k,1)
cPL20/8/04 2 bugs fixed here:
            polfa3 = polfa3 + polfa2 * foutra(k,2)
            polfa4 = polfa4 + polfa2 * foutra(k,3)
            end do
	    fa = dtheab(j) * pi2i
          polfac = polfac * fa 
          polfa3 = polfa3 * fa
          polfa4 = polfa4 * fa
          rcura(i,j) = torfac * polfac
          rcurac(i,j) = torfac * polfa3
          rcuras(i,j) = torfac * polfa4          
          end if
c         ~~~~~~

  6     continue ! Poloidal mode loop
c       ========      
c     For display in array viewer:
      rcs(1:nmoant) = dreal(rcura(1:nmoant,j))
      ics(1:nmoant) = dimag(rcura(1:nmoant,j))
      rcs(1:nmoant) = dreal(cura(1:nmoant,j))
      ics(1:nmoant) = dimag(cura(1:nmoant,j))
  2   continue   ! Antenna loop
c     ========

c     Case of no feeders: contribution to the charge.
c     This is for demonstration purposes only: feeders should always be there
c     to conserve current!
      if(circ)then
c     Not done in general geometry. No need for that.
        do jant = 1, neffan
        if(.not.feeder(jant))then
        iregb = irbant(jant)
        rhoant = rx0m(iregb)
        isubr = nsum(iregb,ns)
c        intab = iregb + (ngauss + 1) * nsum(isubr,iele)
        intab = ifiabs(ifiel(iregb+1)) - 1

          do i = 1, nmoant
          m = moant(i)
      	if(cyl)then
      	chaa(i,jant) = chaa(i,jant) + dfloat(m) / (omegag*rhoant) * cura(i,jant)
      	else
      	k1 = max0(-ncrot, moant(1)-m )
      	k2 = min0( ncrot, moant(nmoant)-m)
      	  do k = k1, k2
      	  mpk = m + k
      	  kp = iabs(k)
      	  if(kp.le.ncrot)chaa(i,jant) = chaa(i,jant) + dfloat(mpk) / (omegag*r0*rhoant) * rcura(j,jant) * bfou(kp)
      	  end do
      	end if
          end do
        end if
        end do
      end if
        
c       Axisymmetric Faraday shields:
        if(nscree .gt. 0)then
	    if(.not. circ .or. cokpco)then
	    write(nofile,*)'Faraday screen remains to implement in general geometry,'
	    write(nofile,*)'and for const. k// coordinates. I stop.'
	    call exit
	    end if
          if(monomo)then
          nmoscr = 1
          moscr(1) = 0
          else
          im = 1 + (nmoscr - 1) / 2
      	do i = 1,nmoscr
      	moscr(i) = i - im
            end do
          end if
        irscrm = nreg + 1

        do 10 is = 1, nscree
          irscrm = min(irscrm, irbscr(is)+1)
          ireg = irbscr(is)

          rhoant = rx0m(ireg)
          y = rx0(ireg)
          isubr = nsum(ireg,ns)
c          intab = ireg + (ngauss + 1) * nsum(isubr,iele)
          intab = ifiabs(ifiel(ireg+1)) - 1

	      if(circ)then
c  NB: correct use of inputs thes1, thsa2 requires no 2*pi angle jump on antenna!!!
      	theas(is) = (thes1(is) + thes2(is)) * 0.5d0
      	dthes(is) =  thes2(is) - thes1(is)
            thesb1(is) = thes1(is)
            thesb2(is) = thes2(is)
            thbs(is) = theas(is)
            dthesb(is) = dthes(is)
	      else
            tth1 = thes1(j)
            if(tth1.lt.0.d0)tth1 = tth1 + twopi
c           (tth1 strictly < 2 pi)
            tth2 = thes2(j)
            if(tth2.lt.0.d0)tth2 = tth2 + twopi
	      call interp1(polang, polang(npfft+1), ckt(intab,1,3), nabplo, npfft+1, tth1, thesb1(is), 1, 'R')
	      call interp1(polang, polang(npfft+1), ckt(intab,1,3), nabplo, npfft+1, tth2, thesb2(is), 1, 'R')
            thbs(is) = 0.5d0*(thesb1(is) + thesb2(is))
            dthesb(is) = thesb2(is) - thesb1(is)
	      end if

            if(circ)then
      	  call fouco
      	  else
	        write(nofile,*)'Faraday screen remains to implement in general geometry!'
c       	  call foucog(1)
c       	  to do in noncircular!
      	  end if

      	  do 10 im = 1, nmoscr
      	  iexc = neffan + (is - 1)*nmoscr + im

      		do jm = 1, nmoant
      		mml = moant(jm) - moscr(im)
c to do in general geom:
c      		polfac = lamwid(is) / rx0m(ireg) * cdeid(-mml*theas(is))
c 5/12/03: for const. k// coord. 
c Treatment of exp(i*n*phibar(theta)) remains to implement.
      		polfac = lamwid(is) / rx0m(ireg) * cdeid(-mml*thbs(is))
      		  if( nsclam(is) .eq. 1 )then
      		  xtra = 0.d0
      		  else
c      		  xtra = mml * dthes(is) / (nsclam(is) - 1) * 0.5d0
c 5/12/03: for const. k// coord.
      		  xtra = mml * dthesb(is) / (nsclam(is) - 1) * 0.5d0
      		  end if

      		  if( dmod(xtra,pi) .eq. 0.d0 )then
      		  polfac = polfac * nsclam(is)
      		  else
      		  polfac = polfac
     ;		   * dsin(nsclam(is)*xtra) / dsin(xtra)
c     ;		   * sin(nsclam(is)*xtra) / sin(xtra)
      		  end if
      		xtra = mml * lamwid(is) / rx0m(ireg) * 0.5d0
c      		if(xtra .ne. 0.d0)polfac = polfac * sin(xtra) / xtra
      		polfac = polfac * dsinc(xtra)
      		curat(jm,iexc) = polfac
      		if(cyl)rcurat(jm,iexc) = curat(jm,iexc)
            end do

          do jm = 1, nmoant
      	do km = 1, nmoant
      	kp = iabs(moant(km)-moant(jm))
            if(kp .le. ncrot)then
c see: toroidal strap
              if(circ)then
c Here rcurat is modes of R * jphi:
      	  if(.not.cyl)rcurat(jm,iexc) = rcurat(jm,iexc) + r0 * afou(kp) * curat(km,iexc)
              chaa(jm,iexc) = chaa(jm,iexc) + kphi / omegag * bfou(kp) * curat(km,iexc)
      	  end if
            else
c to write...        
            end if
            end do
          end do
   10 	continue
c     End of code for Faraday shields
      end if

      do ja = 1, neffan
        do i = 1, nmoant
        write(nofile,1000)moant(i), cura(i,ja), rcura(i,ja), rcurat(i,ja), chaa(i,ja) / eps0
c        write(12345,*) moant(i), dreal(cura(i,ja)), dreal(rcura(i,ja)), 
c     ; dreal(rcurat(i,ja)), dreal(chaa(i,ja) / eps0)
c	        write(12346,*) moant(i), dimag(cura(i,ja)), dimag(rcura(i,ja)), 
c     ; dimag(rcurat(i,ja)), dimag(chaa(i,ja) / eps0)
	  end do
      end do
 1000 format(1h , i4, 4(1h(, g15.7, 1h,, g15.7, 1h), 2x))

      write(nofile,*)'Exit filant'
      return
      end
