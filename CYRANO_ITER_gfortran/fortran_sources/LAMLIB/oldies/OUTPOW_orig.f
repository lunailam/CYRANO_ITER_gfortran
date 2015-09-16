      subroutine outpow

      implicit none

c     Builds output tables, to be plotted elsewhere.
c     New, shortened version (1998) of sub. OUTPUT, free of graphics software.
 
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
      include 'comgdr.copy'
      include 'compla.copy'
      include 'comfin.copy'
      include 'comfic.copy'
      include 'comequ.copy'
      include 'comin2.copy'
      include 'comphy.copy'
      include 'coequi.copy'       
      logical onespe
 
      character*3  eletyp
      character*9  charfi(2), charco(6), charpa(3)
      character(100) :: GGs	!ERN
	character(50) filename	!ERN
	character(4) pana4	!ERN
	 
      integer 
     ;  i, i1, j, ii, aqp(3), aux
     ;, iplom1, ipol
     ;, is1, is2
     ;, ideca
     ;, idllo, icolo, ibulo, icpblo
     ;, ie1
     ;, ima, iblk, ireq, j1, j2, ispe, ib
     ;, istatu, ie2, nintma
     ;, minfc, nblk, ist
     ;, istg, nan, ian, ipir, isp
     ;, nsum, i1st, ceilq
     ;, idamin, idmin, idmax

      double precision 
     ;  rle
     ;, flipi, tabmin, tabmax, tabpoi, flh
     ;, dsum, cabs2, rveno
     ;, scpow, powquo(nabplo), pcheck
     ;, fac5, trav(maxthp+1,2), t(2), tt
     ;, picm10, rm
     ;, trapeze, gauss

      complex*16
     ;  zdotc, zdotu
     ;, tra1
 
      external 
     ;  nsum, i1st, ceilq, trapeze, gauss, flh, rveno
     ;, zdotc, zdotu, zset, dset, dsum, mucrvz
     ;, idamin, idmin, idmax
 
      character*48 text(14)
      character*75 title, title2

	 
      double precision 
     ;  second

      external second
c
      data charfi/'electric ','magnetic '/
     ;    ,charco/'radial ','poloidal ','toroidal ',
     ;            'L.H.(+) ','R.H.(-) ','parallel '/
     ;    ,charpa/'real','imag','modulus'/
 
      cabs2(tra1) = dreal(tra1 * dconjg(tra1))
 
c**********************************************************************
      if(.not.(plopow.or.transp) .or. glovac)return
	write(nofile,*)
     ;'*******************************  ENTER OUTPOW  *******************************'
      
c     Power absorption density on magnetic surfaces for all species
c     -------------------------------------------------------------
      iplom1 = nploth + 1
      ig = 1
      onespe = .true.
      call dset((nspec+1)*3*nreg, 0.d0, regpow, 1)
      call dset(nele, 0.d0, elpow, 1)
      call dset((nspec+1)*2, 0.d0, parpow, 1)
      scpow = 0.d0
      call dset((nspec+2)*nabplo*2, 0.d0, powden, 1)
      flipi = 1.d0 / dfloat(nploth)


      do 100 ispe = 1, nspec
c     ======================
      call dset((maxthp+1)*nabplo*2, 0.d0, tab, 1)
	print *, paname(ispe),"...."
 
      iblk = 1
      if(.not.monoto)call aqwait(aqp,istatu)

      do 307 imoto = 1, ntotor
c     Read modal field components:
      if(.not.monoto)then
      nblk = ceilq(nwn*5*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqread(aqp,xrtp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*3*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqread(aqp,hrtp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*ndof*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqread(aqp,xpmp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*nficom*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqread(aqp,xpmp2,iblk,nblk,ireq,0,istatu)
      iblk = iblk + nblk
      call aqwait(aqp,istatu)
c  3d: define kprn!
      end if
 
      i = 0
      ist = 1
      do 1045 ireg = 1, nreg
c     ----------------------
      i = ist - 1
        if(.not.vacuum(ireg))then
        is1 = nsum(ireg-1, ns) + 1
        is2 = is1 + ns(ireg) - 1
          do 1036 isubr = is1, is2
          eletyp = styp(isubr)
          idllo = iddl(isubr)
          icolo = iconn(isubr)
          ibulo = ibub(isubr)
          icpblo = idllo-icolo
          ie1 = i1st(isubr) + 1
          ie2 = ie1 + iele(isubr) - 1
            do 1036 iel = ie1, ie2
            rle = fl(iel)
            ielm = iel
            if(nmode(iel-1) .gt. nmode(iel))ielm = iel - 1
            minfc = minf(ielm)
            ima = minfc + 1 - minf(imax)
cERN              if(outgau)then
                 j1 = 1
                 if(iel.eq.ifiel(ireg))j1 = 0
                 j2 = ngauss + 1
cERN              else
cERN                 nintma = ninter
cERN                 if(iel .eq. ilael(ireg))nintma = ninter + 1
cERN                 j1 = 0
cERN                 j2 = nintma
cERN              end if
             
            do j = j1, j2
            intab = j + (iel - 1) * (ngauss + 1) + ireg
            i = intab
            y = absono(i)
            rho = absout(i)
c           flr corrections care:
            yinv = abscni(intab)
            rhoinv = yinv * rnori
	print*,i
cERN  : Speed test - remove all poloidal interp. (later implement OUTNPFT)
 
cERN              if(.not.updsym)then
cERNc             interpolate equilibrium at plot points, no up-down symmetry.
cERNc             this is only correct if outgau=.t.! otherwise 2d interpol.
cERN              call interp1(polang, polang(npfft+1), bmotab(intab,1), nabplo,
cERN     ;        npfft+1, polplo, bmlvec, iplom1, 'R')
cERN              call interp1(polang, polang(npfft+1), r0orta(intab,1), nabplo,
cERN     ;        npfft+1, polplo, r0rvec, iplom1, 'R')
cERN                if(.not.circ)then
cERN                call interp1(polang, polang(npfft+1), eqt(intab,1,14), nabplo,
cERN     ;          npfft+1, polplo, cov, iplom1, 'R')
cERN                call interp1(polang, polang(npfft+1), eqt(intab,1,15), nabplo,
cERN     ;          npfft+1, polplo, siv, iplom1, 'R')
cERN                call interp1(polang, polang(npfft+1), eqt(intab,1,7), nabplo,
cERN     ;          npfft+1, polplo, ntnt, iplom1, 'R')
cERN                end if

				bmlvec(1:nploth+1) = bmotab(intab,1:npfft+1)
				r0rvec(1:nploth+1) = r0orta(intab,1:npfft+1)
				cov(1:nploth+1)    = eqt(intab,1:npfft+1,14)
				siv(1:nploth+1)    = eqt(intab,1:npfft+1,15)
				ntnt(1:nploth+1)   = eqt(intab,1:npfft+1,7)

cERN              else

c             Use up-down symmetry:
c             this is only correct if outgau=.t.! otherwise use 2d interpol.
cERN              ian = nploth / 2 + mod(nploth,2)
cERN              nan = nploth - ian + 1
cERN              call interp1(polang, polang(npfft/2+1), bmotab(intab,1), nabplo,
cERN     ;        npfft/2+1, polplo, bmlvec, ian, 'R')
cERN              call dcopy(ian, bmlvec, 1, bmlvec(nan), -1)
cERN              call interp1(polang, polang(npfft/2+1), r0orta(intab,1), nabplo,
cERN     ;        npfft/2+1, polplo, r0rvec, ian, 'R')
cERNc             now r0rvec contains Ra over R, not the inverse!
cERN              call dcopy(ian, r0rvec, 1, r0rvec(nan), -1)
cERN                if(.not.circ)then
cERN                call interp1(polang, polang(npfft/2+1), eqt(intab,1,14),
cERN     ;          nabplo, npfft/2+1, polplo, cov, ian, 'R')
cERN                call dcopy(ian, cov, 1, cov(nan), -1)
cERN                call interp1(polang, polang(npfft/2+1), eqt(intab,1,15),
cERN     ;          nabplo, npfft/2+1, polplo, siv, ian, 'R')
cERN                call dcopy(ian, siv, 1, siv(nan), -1)
cERN                call interp1(polang, polang(npfft/2+1), eqt(intab,1,7),
cERN     ;          nabplo, npfft/2+1, polplo, ntnt, ian, 'R')
cERN                call dcopy(ian, ntnt, 1, ntnt(nan), -1)
cERN                else
cERN                si = sitab(intab)
cERN                co = cotab(intab)
cERN                call dset(nploth, co, cov, 1)
cERN                call dset(nploth, si, siv, 1)
cERN                call dset(nploth, eqt(intab,1,7), ntnt, 1)
cERN                end if
cERN              end if



      if(nspgdr.gt.0 .and. spegdr(ispe))then
c     ======================================
c     Species with general dielectric response
c     compute x, v meshes at current radius and interpolate f0:
      call intf0
c     geometrical couplings for dielectric response and k// list:
cERN      call fougdr(2)
		call fougdr_ern(2)
      i1 = 0
        do rm = dfloat(minf(iel)), dfloat(msup(iel)), 0.5d0
          i1 = i1 + 1
cERN      allkpa(i1) = hachi(intab) * (rm + n * eqta1d(intab,4))
	    allkpa(i1) = hachi(intab) * (rm + n * qfactor(intab))
        end do
c     assuming 1 gdr species: total abs. at current radius
      call powgdr(1, pic10, plan, picm10)
c     check normaliz.:
      powden(i,ispe,1) = powden(i,ispe,1) + plan + pic10 + picm10
      powden(i,ispe,2) = powden(i,ispe,2) + plan + pic10 + picm10

      else if(cokpco)then
c	===================
c     Maxwellian in const. k// coordinates
c     geometrical couplings for dielectric response and k// list:
cERN      call fougdr(1)
		call fougdr_ern(1)
      i1 = 0
        do rm = dfloat(minf(iel)), dfloat(msup(iel)), 0.5d0
           i1 = i1 + 1
cERN       allkpa(i1) = hachi(intab) * (rm + n * eqta1d(intab,4))
		 allkpa(i1) = hachi(intab) * (rm + n * qfactor(intab))
        end do
      write(nofile,*)'OUTPOW line 244: Ernesto, this is not written yet!'
c     The rest of this 'else' needs to be written 
c     or done during call to powabs inside next else

        do ipol = 1, nploth
c       for dirty passing to powabs; intabp is passed as a plot point index:
        intabp = ipol
        bmodul = bmlvec(ipol)
        r0or = r0rvec(ipol)
        ror0 = 1.d0 / r0or
        call zcopy(nmode(ielm), tabexp(ima,ipol), 1, copisi, 1)
         if(.not.circ)then
         si = siv(ipol)
         co = cov(ipol)
          end if
        phi = polplo(ipol)
c write contrib for cokpco case here:
ccc      call powabs(i, onespe, ispe, .false.)
c     tab(,,1) to give absorption to order 0 in lr;
c     tab(,,2) to give total with absorption up to 2nd order in lr
        tab(ipol,i,1) = plan + pic10
c       tab(ipol,i,1) = plan + pic10 + psd + pcoda
        tab(ipol,i,2) = tab(ipol,i,1) + pic2 + pic11 + pic12 + pel2 + pttmp

c       times normalized volume element R/Ra * Jn: 
        trav(ipol,1) = tab(ipol,i,1) * ror0 * eqt(intab,ipol,9)
        trav(ipol,2) = tab(ipol,i,2) * ror0 * eqt(intab,ipol,9)
 
        end do

c     "flux surface averaged" absorption [but using volume element!]: 
      powden(i,ispe,1) = powden(i,ispe,1) + dsum(nploth, trav, 1) * flipi
      powden(i,ispe,2) = powden(i,ispe,2) + dsum(nploth, trav(1,2), 1) * flipi

      else
c     ====
      do ipol = 1, nploth
c     for dirty passing to powabs; intabp is passed as a plot point index:
      intabp = ipol
      bmodul = bmlvec(ipol)
      r0or = r0rvec(ipol)
      ror0 = 1.d0 / r0or
      call zcopy(nmode(ielm), tabexp(ima,ipol), 1, copisi, 1)
        if(.not.circ)then
        si = siv(ipol)
        co = cov(ipol)
        end if
      phi = polplo(ipol)
c write contrib of general diel. response and output for next step of qlfp
c elsewhere!!!
      call powabs(i, onespe, ispe, .false.)
c     tab(,,1) to give absorption to order 0 in lr;
c     tab(,,2) to give total with absorption up to 2nd order in lr
      tab(ipol,i,1) = plan + pic10
c      tab(ipol,i,1) = plan + pic10 + psd + pcoda
      tab(ipol,i,2) = tab(ipol,i,1) + pic2 + pic11 + pic12 + pel2 + pttmp

c     times normalized volume element R/Ra * Jn: 
      trav(ipol,1) = tab(ipol,i,1) * ror0 * eqt(intab,ipol,9)
      trav(ipol,2) = tab(ipol,i,2) * ror0 * eqt(intab,ipol,9)
 
      end do

c     "flux surface averaged" absorption [but using volume element!]: 
      powden(i,ispe,1) = powden(i,ispe,1) + dsum(nploth, trav, 1) * flipi
      powden(i,ispe,2) = powden(i,ispe,2) + dsum(nploth, trav(1,2), 1) * flipi

      end if
c     ======

      end do
 1036 continue
 
      end if
      ist = ist + istp(ireg)
 1045 continue
c     --------
 307  continue

	print*,'... end of radial loop'
c       against trouble with flr terms at axis when fields obtained from flr0:
        if(.not.flrops(1))then
        powden(1,ispe,2) = powden(2,ispe,2)
        call dcopy(iplom1, tab(1,2,2), 1, tab(1,1,2), 1)
        end if
c     Against trouble in general: 
c     Set 1st radial abscissa values to 2nd radial abscissa values:
      powden(1,ispe,1) = powden(2,ispe,1)
      call dcopy(iplom1, tab(1,2,1), 1, tab(1,1,1), 1)
 
      ist = 1
      do ireg = 1, nreg
	   if(.not.vacuum(ireg))then
		  if(flrops(ireg))then
			 ii = 2
		  else
			 ii = 1
		  end if
		  tt = 0.d0
cERN		  if(outgau)then  
		     istg = ist - ngauss
		     do iel = ifiel(ireg), ilael(ireg)
c			    Element length in m:
			    rle = rnorm * fl(iel)
			    istg = istg + ngauss + 1
c			    Integrate powden * rho over elements: 
			    t(1) = gauss(rle, absout(istg), powden(istg,ispe,1), 1, 1)
			    t(2) = gauss(rle, absout(istg), powden(istg,ispe,2), 1, 1)
c			    Powers per species per region:
			    regpow(ispe,1,ireg) = regpow(ispe,1,ireg) + t(1)
			    regpow(ispe,2,ireg) = regpow(ispe,2,ireg) + t(2)
c			    Sums for calculation of power abs. consistent with wave equation:
			    regpow(ispe,3,ireg) = regpow(ispe,3,ireg) + t(ii)
c			    Sum for total active power per element:
			    elpow(iel) = elpow(iel) + t(ii)
			    tt = tt + t(ii)
		     end do ! (iel)	
cERN            else
cERN		     t(1) = trapeze(istp(ireg), absout(ist), powden(ist,ispe,1), 1, 1)
cERN		     t(2) = trapeze(istp(ireg), absout(ist), powden(ist,ispe,2), 1, 1)
cERN		     parpow(ispe,1) = parpow(ispe,1) + t(1)
cERN		     parpow(ispe,2) = parpow(ispe,2) + t(2)
cERN		     tt = tt + t(ii)
cERN	      end if

c		  Powers per species:
		  parpow(ispe,1) = parpow(ispe,1) + regpow(ispe,1,ireg)
		  parpow(ispe,2) = parpow(ispe,2) + regpow(ispe,2,ireg)
c		  Total powers per region:
		  regpow(nspec+1,1,ireg) = regpow(nspec+1,1,ireg) + regpow(ispe,1,ireg)
		  regpow(nspec+1,2,ireg) = regpow(nspec+1,2,ireg) + regpow(ispe,2,ireg)
c		  Total power per region, consistent with wave equation::
		  regpow(nspec+1,3,ireg) = regpow(nspec+1,3,ireg) + regpow(ispe,3,ireg)
		  scpow = scpow + tt
         end if
	   ist = ist + istp(ireg)
      end do !(ireg)
 
c     Periodicity: copy values at theta=0 to theta=2*pi:
      call dcopy(ist11, tab, maxthp+1, tab(iplom1,1,1), maxthp+1)
      call dcopy(ist11, tab(1,1,2), maxthp+1, tab(iplom1,1,2), maxthp+1)
      
      parpow(nspec+1,1) = parpow(nspec+1,1) + parpow(ispe,1)
      parpow(nspec+1,2) = parpow(nspec+1,2) + parpow(ispe,2)
c     Energy balance check: only take mechanisms included in wave
c     pattern calc. into account -> scpow:
cc     temporary, for deuterium 1st harmonic estimate:
c        if(usenma .or. (equidf(ispe).eq.'maxwe' .and. amass(ispe).ne.2.) )
c     ;  then
c        scpow = scpow + parpow(ispe)
c          do i = 1, ist11
c          powden(i,nspec+2) = powden(i,nspec+2) + powden(i,ispe)
c          end do
c        end if

c	Total power
	powden(1:nabsci,nspec+1,1) = powden(1:nabsci,nspec+1,1)
     ;                          + powden(1:nabsci,ispe,1)
	powden(1:nabsci,nspec+1,2) = powden(1:nabsci,nspec+1,2)
     ;                          + powden(1:nabsci,ispe,2)

 
cERN	ccccc Writing 2D power densities to separate files cccccccccccccc

      if(plot2d)then
	   
	   pana4 = paname(ispe)
	   do ipir = 1, 2   ! FLR0 and FLR2
		 ideca = (ipir - 1) * (maxthp+1) * nabplo
		 tabmin = tab2(ideca + idmin((maxthp+1)*nabplo, tab(1,1,ipir), 1))
		 tabmax = tab2(ideca + idmax((maxthp+1)*nabplo, tab(1,1,ipir), 1))
		 if(tabmin.ne.tabmax  .and.
     ;     (abs(tabmin).gt.1.d-60 .or. abs(tabmax).gt.1.d-60))then
			iglplo = iglplo + 1
			if(ipir.eq.1)then
			   title = '2D power density - ' // paname(ispe) // ' - FLR0'
			   filename = 'power2D_' //	pana4 // '_FLR0.dat'
			else
			   title = '2D power density - ' // paname(ispe) // ' - FLR2'
			   filename = 'power2D_' //	pana4 // '_FLR2.dat'
			end if
			write(nofile,1000) iglplo, '2d ' // title
			ploleg(iglplo) = '2d ' // title
			
			print *, filename
			open (UNIT = 777, FILE = paplda // filename, 
     ;			  STATUS = "REPLACE", ACTION = "WRITE")			
			   write(777,*), title
			   write(777,*), ' rad ', ' pol ', ' pla '				  
			   write(777,"(i5,i5,i5)"), ist11, iplom1, ipla
			   do i = 1, ist11
	              aux = intaplot(i)
			      write(777,2222)(tab(j,aux,ipir), j = 1, iplom1)
			   end do
			close(777)

		 end if ! (tabmin .ne. tabmax & ...)
	   end do ! (ipir = 1,2) 

      end if ! (plot2d)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  100 continue  ! End of species loop (ispe = 1,nspec)
c     ========

cERN	Aplly geom. factor to 1D power densities

	if(cyl) then
	   fac5 = twopi ** 2
      else
	   fac5 = twopi ** 2 * ra
      end if

	parpow(1:nspec+1,:) = parpow(1:nspec+1,:) * fac5
	regpow(1:nspec+1,:,nreg) = regpow(1:nspec+1,:,nreg) * fac5
	elpow(1:nele) = elpow(1:nele) * fac5
	powden(1:nabsci,1:nspec+1,:) = powden(1:nabsci,1:nspec+1,:) * fac5

cERN	call dscal((maxspe+1)*2, fac5, parpow, 1)
cERN      call dscal((maxspe+1)*3*nreg, fac5, regpow, 1)
cERN      call dscal(nele, fac5, elpow, 1)
cERN        do isp = 1, nspec + 1
cERN        call dscal(ist11, fac5, powden(1,isp,1), 1)
cERN        call dscal(ist11, fac5, powden(1,isp,2), 1)
cERN        end do

c     Integral of (rho*powden drho) equals total power:
      scpow = scpow * fac5

c     Cumulative power at element right nodes:
	cumpow(0) = 0.d0
      cumpow(1) = elpow(1)
	  do iel = 2, nele
	  cumpow(iel) = cumpow(iel-1) + elpow(iel)
	  end do

c       Assumes never other flr region than the first one:
        if(flrops(1))then
        ii = 2
        else
        ii = 1
        end if

        if(losles)then
           pcheck = scpow
           write(nofile,*) 'split power balance check: ', scpow,
     ;     ' Watt found in plasma, with antenna input ', tinpow, ' watt'
c          write(nofile,*) 'input active power from poynting flux = ',
c     ;    dreal(totpoy), '- this disregards feeders!'
        else
cERN          do i = 1, istp(1)  ! don't use OUTGAU option here 
	 	    do i = 1, 1 + (ngauss+1)*(ilael(1)-ifiel(1)+1) ! istp(1) for OUTGAU=T
              if(dpoydr(i).ne.0.d0)then
                powquo(i) = powden(i,nspec+1,ii) / dpoydr(i)
              else
                powquo(i) = 1.d0
              end if
           end do

           pcheck = scpow / tinpow
           write(nofile,*) 
           write(nofile,*) 'Power balance check (OUTPOW): '
     ;     , 100.*pcheck, '% of input power found in plasma'
c          if(dreal(totpoy) .ne. 0.d0)
c     ;    write(nofile,*) 100. * scpow / dreal(totpoy)
c     ;,   ' % of poynting flux is found in plasma - nb: this disregards feeders!'
           write(nofile,*)
     ;     'Percentage of poynting flux difference found in each region:'
           write(nofile,*) 
     ;     '------------------------------------------------------------'
           write(nofile,*)'n.b.: will differ from 100% in a region with feeders,'
     ;     // 'and will be unreliable in a lossless region.'
           do ireg = 1, nreg
              write(nofile,*)ireg, 100. * regpow(nspec+1,3,ireg) / dreal(poydir(ireg))
           end do
           write(nofile,*) 'Power fraction to each species:'
           write(nofile,*) '------------------------------'
           write(nofile,*) 'scpow : ', scpow
           write(nofile,*) 'To order 0 in Larmor radius:'
           write(nofile,*) '---------------------------'
           write(nofile,*) 'Raw total power: ', parpow(nspec+1,1)

cPL        6/5/2004: switched some mess off from here!
           do isp = 1, nspec
              write(nofile,*)paname(isp),' : ',parpow(isp,1)
           end do      
           write(nofile,*) 'To order 2 in Larmor radius:'
           write(nofile,*) '---------------------------'
           write(nofile,*) 'Raw total power: ', parpow(nspec+1,2)
           do isp = 1, nspec
              write(nofile,*)paname(isp),' : ',parpow(isp,2)
           end do      
cPL        6/5/2004: switched some mess off above here!

        end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

cERN  cccccccccc Write 1D power absorption to individual files ccccccccccc

	GGs  = "(g16.6," // char(nspec+1+48) // "g16.6)"
c	NB: powden(:,nspec+1,:) is the total power summed over species
c		powden(:,:,1) is FLR0 and powden(:,:,2) is FLR2

c     1) FLR0 power density
	open (UNIT = 701, FILE = paplda // 'powerdens1D_FLR0.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		  write(701,*),"Average power absorption density per species (FLR0)"
		  write(701,"(3A5)"),"nspec","rad","pla"
		  write(701,"(3i5)"), nspec, nabsci, i_rp
		  write(701,*), "rho(m)  ", paname(1:nspec), "total"
            do i = 1, nabsci
	         aux = i
		     write(701,GGs), abscis(aux), powden(aux,1:nspec+1,1)
		  end do
	close(701)

c     2) FLR2 power density
	open (UNIT = 702, FILE = paplda // 'powerdens1D_FLR2.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		  write(702,*),"Average power absorption density per species (FLR2)"
		  write(702,"(3A5)"),"nspec","rad","pla"
		  write(702,"(3i5)"), nspec, nabsci, i_rp
		  write(702,*), "rho(m)  ", paname(1:nspec), "total"
            do i = 1, nabsci
			 aux = i
		     write(702,GGs), abscis(aux), powden(aux,1:nspec+1,2)
		  end do
	close(702)

c     3) dPoynt/dr
	open (UNIT = 703, FILE = paplda // 'dPoynting_dr.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		  write(703,*),"Radial derivative of inward Poynting flux (dPoynt/dr)"
            write(703,"(2A5)"),"rad","pla"
		  write(703,"(2i5)"), nabsci, i_rp
		  do i = 1, nabsci
	         aux = i
		     write(703,GGs), abscis(aux),  dpoydr(aux)
		  end do
	close(703)

c     4) Power quotient 
	open (UNIT = 701, FILE = paplda // 'power_quotient.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		  write(701,*),"sp.absorbing in wave calc./dPoydr"
	      write(701,"(2A5)"),"rad","pla"
		  write(701,"(2i5)"), nabsci, i_rp
            do i = 1, nabsci
			 aux = i
			 write(701,GGs), abscis(aux), powquo(aux)
		  end do
	close(701)

cERN  5) Total power per species (NEW) + central densities
	open (UNIT = 701, FILE = paplda // 'totalpower.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		  write(701,*),"Raw total power per species"
		  write(701,*), nspec
		  write(701,*), paname(1:nspec), "total"
            write(701,*), parpow(1:nspec+1,1)
		  write(701,*), parpow(1:nspec+1,2)
		  write(701,*), "Species central densities"
		  write(701,*), dentab(1,1:nspec)
	close(701)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Active power per element, written with block index in global matrix and
c     radial coordinate of element centre:
      open(unit=80
     ;, file = paplda // 'Element active powers (OUTPOW)', status='unknown')
      write(80,*)
     ;'Total absorbed power per element (from direct calculation in outpow):'
      write(80,*)nele+nreg, 4
        do iel = 1, nele
        write(80,2001)float(iel), rnorm*(fx0(iel-1)+0.5*fl(iel))
     ;              , elpow(iel), eltvol(iel)
        end do
	close(unit=80, status='KEEP')

c     Cumulative power at element nodes, written with radial abscissa:
      open(unit=80
     ;, file = paplda // 'Cumulative power at element nodes (OUTPOW)'
     ;, status='unknown')
      write(80,*)
     ;'Cumulative power vs rho (from direct calculation in outpow):'
      write(80,*)nele+1, 2
	  do iel = 0, nele
	  write(80,*)fx0(iel)*rnorm, cumpow(iel)
	  end do
	close(unit=80, status='KEEP')

	write(nofile,*)
     ;'*******************************  EXIT OUTPOW  *******************************'
      return

 1000 format(1h ,'plot #',i4,'  ',a50)
 2000 format(1h ,14(1x,g12.4))
 2001 format(1h ,4(2x,g16.8))
 2222 format(200(g14.6))  
      end
