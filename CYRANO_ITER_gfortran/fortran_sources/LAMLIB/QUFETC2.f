      subroutine qufetc2(onebl, iblock, toread, clean, ini, icurspe)

c      subroutine qufetc2(onebl, iblock, toread, clean, ini)

      implicit none

      logical onebl, toread, clean, ini
      integer iblock, icurspe
c
c     Computes residuals and quadratic forms for each element using solution
c     of linear system.
c
c     onebl=.true.: do it only for input block index iblock
c           .false.: for all blocks at once
c     iblock: only used if onebl=.true., dummy otherwise.
c             index of block to use
c     toread: 
c     toread=.true.: read element matrix and rhs with call to rdelt
c            .false.: the blocks are already stored in AEL and BEL
c     clean=.true.: set matrices AEL and BEL to zero after use
c           .false.: don't
c     ini = .true.: initialize various arrays and variables.
c           .false.: don't
c
c     icurspe: current species index, only used if computing element power species by species,
c              i.e. when recalc_by_species_now is true.

      include 'pardim.copy'

      include 'comusr.f'
      include 'comequ.copy'
      include 'comfin.copy'
      include 'comreg.copy'
      include 'comant.copy'
      include 'complp.copy'
      include 'comgeo.copy'
      include 'compla.copy'

      logical donewr, writenow

      integer i, is, j, js, k, ib, iel, lbl(2), lblt, ish, ia, ja, ireg
     ;, imres, imrres, imires, ibmres
     ;, ibmrr, ibmir
     ;, iel1, iel2, isp
      double precision cabs1, resman, rrman, irman
     ;, ar, ai, a1
      complex*16 
     ;  xsbmax(maxnbl), xsax(maxnbl), sysquf, soupow, matpow_raw, matpow
     ;, czero, zdum
     ;, residu(maxbll)
     ;, tofacv(maxant), tofac, zed, zed2, zb
     ;, xsb(maxnbl), solvec(2*maxbll)
     ;, rawzij(maxexc,maxexc), zij_quf(maxant,maxant)
     ;, source(2*maxbll)
     ;, xsax_sp(maxnbl,maxspe)
     ;, cbl, cbli, cbl2, cbl2i

      external tofac

      data czero/(0.d0,0.d0)/
c     ;, donewr/.false./
      save 
     ;  xsbmax, sysquf, soupow, residu, xsb, resman, rrman, irman
     ;, imres, imrres, imires, ibmres
     ;, ibmrr, ibmir
     ;, tofacv
     ;, xsax_sp
c     ;, donewr

      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))

        if(ini)then
        call zset(nthoma, czero, xsbmax, 1)
        call zset(nthoma, czero, xsax, 1)
        call zset(nthoma, czero, xsb, 1)
        call zset(maxnbl*nspec, czero, xsax_sp, 1)
c       For total quadratic form:
        sysquf = czero
c       For total source power:
        soupow = czero
c       For residual, maximum norm of residual, of its real and imag.parts:
        resman = 0.d0
        rrman = 0.d0
        irman = 0.d0
c       For contribution of current block to global residual:
        call zset(maxbll, czero, residu, 1)
c       Antenna toroidal factors for current value of kphi or n, times
c       antenna basis current:
c       (check / revise if any_falen!)
          do ia = 1, ntoant
          tofacv(ia) = tofac(ia) * tancur(ia)
          end do
c       For 'raw' impedance matrix:
          do ia = 1, nrhs
          do ja = 1, nrhs
          rawzij(ia,ja) = czero
          end do
          end do      
        end if

        if(onebl)then
        iel1 = iblock
        iel2 = iblock
        else
        iel1 = 1
        iel2 = nthoma - 1
        end if

c     Main loop: (over system blocks rather than elements)
      do iel = iel1, iel2
c     If required, read element and rhs matrices:
      if(toread .or. .not.onebl)call rdelt(1, iel, lbl(1), lbl(2), nrhs)
      lblt = lbl(1) + lbl(2)

      call zset(lblt, czero, solvec, 1)
      call zset(lblt, czero, source, 1)
      ish = 0
c     Loop over left and right element blocks:
      do ib = 1, 2
c     Use block ctom to store solution blocks:
      call rdsol(1, iel+ib-1, ctom, ldct, lbl(ib), nrhs)
c     'Raw' input impedance matrix:
        do ia = 1, nrhs
        do ja = 1, nrhs
          do i = 1, lbl(ib)
          rawzij(ia,ja) = rawzij(ia,ja) + dconjg(ctom(i,ia)) * bel(i+ish,ja)
          end do
        end do
        end do
c     Current excitation case: store solution in solvec:
        do j = 1, ntoant
	  js = j
	  if(samant)js = 1
        do i = 1, lbl(ib)
        solvec(i+ish) = solvec(i+ish) + tofacv(j) * ctom(i,js)
        end do
        end do
      ish = lbl(1)
      end do

c     Excitation case under study:
c     store corresponding source term in source:
cPL        do i = 1, lblt
cPL        source(i) = bel(i,1) * tofacv(1)
cPL        bel(i,1) = bel(i,1) * tofacv(1)
cPL        end do
        do j = 1, ntoant
	  js = j
	  if(samant)js = 1
cPL        do j = 2, nrhs
          do i = 1, lblt
          source(i) = source(i) + tofacv(j) * bel(i,js)
cPL        bel(i,1) = bel(i,1) + tofacv(j) * bel(i,j)
          end do
        end do
c     Initialize 2nd column of BEL to zero (used in elt quadratic form calc.):
      call zset(lblt, czero, bel(1,2), 1)

c     Work of sources on solution, for each element: transconj(x).b
        do k = 1, lblt
        xsb(iel) = xsb(iel) + dconjg(solvec(k)) * source(k)
cPL080404        xsb(iel) = xsb(iel) + dconjg(solvec(k)) * bel(k,1)
        end do

c Now assumed finished with contents of source

cPLc     Residual b-Ax for each element is computed in first col. of block BEL:
c     Residual b-Ax for each element is computed in source:
        do k = 1, lblt
          do i = 1, lblt
cPL          bel(i,1) = bel(i,1) - ael(i,k) * solvec(k)
          source(i) = source(i) - ael(i,k) * solvec(k)
          end do
        end do
c     Quadratic form transconj(x).(b-Ax) for each element:
        do k = 1, lblt
cPL        xsbmax(iel) = xsbmax(iel) + dconjg(solvec(k)) * bel(k,1)
        xsbmax(iel) = xsbmax(iel) + dconjg(solvec(k)) * source(k)
        end do
c     Product A.x for present element is computed in 2nd col. of block BEL:
        do k = 1, lblt
          do i = 1, lblt
          bel(i,2) = bel(i,2) + ael(i,k) * solvec(k)
          end do
        end do
c     Quadratic form transconj(x).A.x for present element:
        do k = 1, lblt
        xsax(iel) = xsax(iel) + dconjg(solvec(k)) * bel(k,2)
        end do
	    if(recalc_by_species_now)then
            do k = 1, lblt
            xsax_sp(iel,icurspe) = xsax_sp(iel,icurspe)
     ;                           + dconjg(solvec(k)) * bel(k,2)
            end do
	    end if
c     Maximum norm of global (assembled system) residual b-Ax:
      do i = 1, lbl(1)
cPL      residu(i) = residu(i) + bel(i,1)
      residu(i) = residu(i) + source(i)
      ar = dabs(dreal(residu(i)))
      ai = dabs(dimag(residu(i)))
      a1 = ar + ai
        if(ar.gt.rrman)then
        rrman = ar
        imrres = i
        ibmrr = iel
        end if
        if(ai.gt.irman)then
        irman = ai
        imires = i
        ibmir = iel
        end if
        if(a1.gt.resman)then
        resman = a1
        imres = i
        ibmres = iel
        end if
      end do

      do i = 1, lbl(2)
cPL      residu(i) = bel(i+lbl(1),1)
      residu(i) = source(i+lbl(1))
      end do
      
	  if(iel.eq.nthoma-1)then
          do i = 1, lbl(2)
          ar = dabs(dreal(residu(i)))
          ai = dabs(dimag(residu(i)))
          a1 = ar + ai
            if(ar.gt.rrman)then
            rrman = ar
            imrres = i
            ibmrr = iel
            end if
            if(ai.gt.irman)then
            irman = ai
            imires = i
            ibmir = iel
            end if
            if(a1.gt.resman)then
            resman = a1
            imres = i
            ibmres = iel
            end if
          end do
        end if

      sysquf = sysquf + xsbmax(iel)
      soupow = soupow + xsb(iel)
      end do

        if(clean)then
          do j = 1, lblt
          do i = 1, lblt
          ael(i,j) = czero
          end do
          end do
          do j = 1, nrhs
          do i = 1, lblt
          bel(i,j) = czero
          end do
          end do
        end if

	writenow = 
     ;      .not.onebl 
     ;.or. (.not.recalc_by_species_now .and. iblock.eq.nthoma-1)
     ;.or. (     recalc_by_species_now .and. iblock.eq.i_last_plasma_block 
     ;                                 .and. icurspe.eq.nspec)
      if(writenow)then
c      if(.not.donewr .and. (.not.onebl .or. iblock.eq.nthoma-1))then
c	donewr = .true.
c     Put active and reactive contributions as real and imaginary parts:
      soupow = dcmplx(0.d0,1.d0) * soupow
      sysquf = dcmplx(0.d0,1.d0) * sysquf

      write(6,*)' '
      write(6,*)'Element block quadratic forms:'
      write(6,*)'----------------------------- '

	  if(.not.onebl)then
c	  This case is called when stomat=.true.
        open(unit=80
     ;  , file = paplda // 'Element quadratic forms (stored blocks)'
     ;  , status='unknown')
        open(unit=82
     ;  , file = paplda // 'Block residuals (stored blocks)'
     ;  , status='unknown')
c	  else
	  else if(.not.recalc_by_species_now)then
cPL	  This case is called when recalc_total=.true.
        open(unit=80
     ;  , file = paplda // 'Element quadratic forms (re-computed blocks)'
     ;  , status='unknown')
        open(unit=82
     ;  , file = paplda // 'Block residuals (re-computed blocks)'
     ;  , status='unknown')

	    if(onlyab_now)then
	    write(80,*)'active power only'
	    write(82,*)'active power only'
	    else
	    write(80,*)'active and reactive powers'
	    write(82,*)'active and reactive powers'
	    end if

	  else  ! recalc_by_species_now:
c	  This case is called when recalc_by_species_now =.true.
cPL icurspe, paname(icurspe)
cPL     ;  , file = paplda // 'Element quadratic forms_' // paname(icurspe)
        open(unit=80
     ;  , file = paplda // 'Element quadratic forms per species'
     ;  , status='unknown')
        open(unit=82
     ;  , file = paplda // 'Block residuals (recomputed blocks by species)'
     ;  , status='unknown')

	    if(onlyab_now)then
	    write(80,*)'active power only'
	    write(82,*)'active power only'
	    else
	    write(80,*)'active and reactive powers'
	    write(82,*)'active and reactive powers'
	    end if

	  end if

      write(6,1000) sysquf
      write(80,1000) sysquf
 1000 format('Global system quadratic form (QUFETC2) = (',g24.16,1h,,g24.16,1h))

      write(6,1001) soupow
      write(80,1001) soupow
 1001 format('Total input power (work of sources on solution, QUFETC2) = ('
     ;,g24.16,1h,,g24.16,1h))

	if(.not.recalc_by_species_now)then
      write(80,*)'i x*.A.x'
c      write(80,*)'x*.A.x'
      write(80,*)nele + nreg, 5
      write(82,*)'i x*.(b-A.x)'
c      write(82,*)'x*.(b-A.x)'
      write(82,*)nele + nreg, 3
        ib = 0
          do ireg = 1, nreg
          do iel = ifiel(ireg), ilael(ireg)
          ib = ib + 1
c         Quadratic form i transconj(x).A.x  and work of sources on solution
c         i transconj(x).b for each block
c         The factor i ensures active+i*reactive:
          write(80,*)float(ib), -dimag(xsax(ib)), dreal(xsax(ib))
     ;                        , -dimag(xsb(ib)),  dreal(xsb(ib))
c          write(80,*)float(ib), dreal(xsax(ib)), dimag(xsax(ib))
c     ;                        , dreal(xsb(ib)), dimag(xsb(ib))
c         Weak form residual i transconj(x).(b-A.x) for each block:
          write(82,*)float(ib), -dimag(xsbmax(ib)), dreal(xsbmax(ib))
c          write(82,*)float(ib), dreal(xsbmax(ib)), dimag(xsbmax(ib))
          end do
c         Combine 2 boundary blocks, except for wall (which only has 1):
          ib = ib + 1
            if(ireg.lt.nreg)then
            zed = xsax(ib) + xsax(ib+1)
            zed2 = xsbmax(ib) + xsbmax(ib+1)
            zb = xsb(ib) + xsb(ib+1)
            else
            zed = xsax(ib)
            zed2 = xsbmax(ib)
            zb = xsb(ib)
            end if
          write(80,*)float(ib), -dimag(zed), dreal(zed), -dimag(zb), dreal(zb)
          write(82,*)float(ib), -dimag(zed2), dreal(zed2)
c          write(80,*)float(ib), dreal(z), dimag(z), dreal(zb), dimag(zb)
c          write(82,*)float(ib), dreal(z2), dimag(z2)
          ib = ib + 1
          end do

	else  ! the case recalc_by_species_now
      write(80,*)'i x*.A.x'
      write(80,*)nspec
	  do isp = 1, nspec
	  write(80,*)paname(isp)
	  end do
      write(80,*)nele, 3+2*(nspec+1)+2
      write(82,*)'i x*.(b-A.x)'
      write(82,*)nele, 3
        ib = 0
          do ireg = 1, nreg
          do iel = ifiel(ireg), ilael(ireg)
          ib = ib + 1
c         Quadratic form i transconj(x).A.x for each element and each species, total by element,
c         and work of sources on solution i transconj(x).b for each element:
c         The factor i ensures format active+i*reactive:
          write(80,*)float(iel), rnorm*(fx0(iel-1)+0.5*fl(iel)), eltvol(iel)
     ;  , (-dimag(xsax_sp(ib,isp)), dreal(xsax_sp(ib,isp)), isp=1,nspec)
     ;  , -dimag(xsax(ib)), dreal(xsax(ib))
     ;  , -dimag(xsb(ib)),  dreal(xsb(ib))
c         Weak form residual i transconj(x).(b-A.x) for each block:
          write(82,*)float(ib), -dimag(xsbmax(ib)), dreal(xsbmax(ib))
c          write(82,*)float(ib), dreal(xsbmax(ib)), dimag(xsbmax(ib))
          end do
c         Skip boundary blocks:
          ib = ib + 2
          end do

	end if  ! .not.recalc_by_species_now

	close(unit=80, status='KEEP')
	close(unit=82, status='KEEP')

      write(6,*)' '
      write(6,*)'Maximum norm of residual / of real part / of imag part:'
      write(6,*)'------------------------------------------------------'
c     Revised convention: factor i:
      write(6,100)'max. residual,  block and d.o.f. #:', resman, ibmres, imres
      write(6,100)'max. real part, block and d.o.f. #:', irman, ibmir, imires
      write(6,100)'max. imag part, block and d.o.f. #:', rrman, ibmrr, imrres
c      write(6,100)'max. residual,  block and d.o.f. #:', resman, ibmres, imres
c      write(6,100)'max. real part, block and d.o.f. #:', rrman, ibmrr, imrres
c      write(6,100)'max. imag part, block and d.o.f. #:', irman, ibmir, imires
  100 format(1h , a35, g13.5, 2x, i4, 2x, i4)
      write(6,*)' '
c      write(6,*)'Total input power (work of sources on solution) = ', soupow

        do ia = 1, nrhs
          do ja = 1, nrhs
          rawzij(ia,ja) = rawzij(ia,ja) * dcmplx(0.d0,2.d0)
          end do
        end do
      write(6,*)' '
      write(6,*)'Raw impedance matrix (effective wire antennas at phi=0):'
      write(6,*)'------------------------------------------------------- '
        do ia = 1, nrhs
        write(6,1002)(rawzij(ia,ja),ja=1,nrhs)
        end do
 1002 format(4(2h (,g24.16,1h,,g24.16,1h)))

c     Antenna toroidal location and shape factors to produce actual impedance matrix:
      write(6,*)' '
      write(6,*)'Array impedance matrix (QUFETC2):'
      write(6,*)'-------------------------------- '

        do ia = 1, ntoant
	  is = ia
	  if(samant)is = 1
          do ja = 1, ntoant
	    js = ja
	    if(samant)js = 1
          zij_quf(ia,ja) = rawzij(is,js) * dconjg(tofac(ia)) * tofac(ja)
          end do
        end do

        if(any_falen)then
        do ja = 1, ntoant
          if(falen(ja))then
          cbl = cdcos(betal(ja))
            if(cbl.ne.(0.d0,0.d0))then
            cbli = 1. / cbl
            else
            cbli = 1.d99
            end if
          end if
        do ia = 1, ntoant
        if(falen(ja) .and. anetyp(ja).eq.'SHC')
     ;  zij_quf(ia,ja) = zij_quf(ia,ja) * cbli
          if(falen(ia) .and. anetyp(ia).eq.'SHC')then
          cbl2 = dconjg(cdcos(betal(ia)))
            if(cbl2.ne.czero)then
            cbl2i = 1.d0 / cbl2
            else
            cbl2i = 1.d99
            end if
          zij_quf(ia,ja) = zij_quf(ia,ja) * cbl2i
          end if
        end do
        end do
        end if

        do ia = 1, ntoant
        write(6,1002)(zij_quf(ia,ja),ja=1,ntoant)
        end do

c temporary condition
	if(.not.recalc_by_species_now)then
      open(unit=80, file = paplda // 'Impedance matrix', status='unknown')
c	write impedance matrix real and imaginary parts for easy further manipulation:
	write(80,*) ntoant, 2*ntoant
        do ia = 1, ntoant
        write(80,1003)(dreal(zij_quf(ia,ja)), dimag(zij_quf(ia,ja)),ja=1,ntoant)
 1003 format(8(1h ,g23.15))
        end do
	close(unit=80, status='KEEP')
	end if

      matpow_raw = czero
        do ia = 1, ntoant
	  is = ia
	  if(samant)is = 1
          do ja = 1, ntoant
	    js = ja
	    if(samant)js = 1
          matpow_raw = matpow_raw + 0.5 * rawzij(is,js)
     ;                                  * dconjg(tofacv(ia)) * tofacv(ja)
          end do
        end do
      write(6,*)'Total input power (from raw impedance matrix, QUFETC2) = '
     ;        , matpow_raw
      write(6,*)' '

      matpow = czero
        do ia = 1, ntoant
          do ja = 1, ntoant
          matpow = matpow + 0.5 * zij_quf(ia,ja)
     ;                          * dconjg(tancur(ia)) * tancur(ja)
          end do
        end do
      write(6,*)'Total input power (from array impedance matrix, QUFETC2) = '
     ;        , matpow
      write(6,*)' '

      end if  ! writenow

      return
      end