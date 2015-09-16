      subroutine genera3

c	USE DFPORT		! ccccccccc ERNESTO cccccccccccccc

      implicit none

c     Generation of finite element problem
c     Assembly of curl.curl and plasma terms element by element

c     Includes loops over average pol. mode index.
C     Hot model with poloidal mode expansion
C     Computes element matrices and calls frontal solver.
C     27/10/2000: this routine has two modes of call:
C     gensys=.true.: assembles element matrices and calls frontal solver;
C            .false.: computes element matrices and calls qufetc2 (element
C            power calculation). There are two variants:
C            a) recalc_by_species_now=.false.: computes total powers
C            a) recalc_by_species_now=.true.:  computes powers per species
C
      include 'pardim.copy'

c     Special use of common COMSOL in following include: redefined larger!!
c     Assumes v, totop not used in that common elsewhere!!
C     NB: because continuation line is not allowed after cf77 include '...'
c     Careful here:
c     File COMUSR2 must be updated every time file COMUSR changes!
c      COMPLEX*16 V(6,6,2*MAXPOM**2), TOTOP(12,12,BLLEN)
      include 'comusr2.f'
c     ;, v, totop

      include 'comreg.copy'
      include 'comsub.copy'
      include 'comgeo.copy'
      include 'comequ.copy'
      include 'comant.copy'
      include 'comfin.copy'
      include 'comfic.copy'
      include 'comin2.copy'
      include 'comrot.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comfou.copy'
      include 'commod.copy'
      include 'commmk.copy'
      include 'compri.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'comphy.copy'
	include 'cokpco.copy'	!ERN

      character*3 eletyp 

      integer 
     ;  index(12), jndex(12)
     ;, i, j, idum, ig1, is1, is2, ifirse, ilaste
     ;, iel1, iel2, nmodc
     ;, idllo, icolo, ibulo, icpblo, nadd2, isca, lbl1, lbl2, lblt
     ;, ielet
     ;, nsum, i1st
     ;, ipvt(maxbll), info, kp1, l, ic, ib, jb, id
     ;, isp

      integer
     ;  lg, lm, ld, it, ll
c     integer ng, nm, nd
      
      double precision 
     ;  cabs1, rle, rlei, rownor(maxbll), scaf(1-maxbll:maxbll)
     ;, second, stiref

      double precision rtem, sc1(maxbll), ri, ri1

      complex*16 zdum

      complex*16 t
c     ;, to(12,12), toto(12*12), tos(12,12), totos(12*12)
c     ;, big(12,12), bigs(12,12), big1(12*12), big1s(12*12)

ccccccccccccccccccc ERNESTO (cokpco) cccccccccccccccccccccccccccccc
	real :: T1(2)	          !   Variables for evaluating
	real :: tempo, tempo0       !   the elapsed CPU time
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      equivalence (toto,to), (totos,tos), (big1,big), (big1s,bigs)

      external second, nsum, i1st

      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
 
      write(nofile,*)'Enter GENERA3 ; time=',second()-timin

      write(nofile,*)' '
        if(gensys)then
        write(nofile,*)'Construction of linear system'
        write(nofile,*)'-----------------------------'
        write(nofile,*)'Test dimensions:'
        write(nofile,*)'ldat, ldbt, ldct, ldat1, ldbt1, ldct1, ldael, lcael
     ;, ldbel, ldbe1, luszzz:'
        write(nofile,*)ldat, ldbt, ldct, ldat1, ldbt1, ldct1, ldael, lcael
     ;, ldbel, ldbe1, luszzz
        else
        write(nofile,*)
     ;  'Computation of wave absorption, recomputing element matrices'
        write(nofile,*)
     ;  '-------------------------------------------------------------'
	    if(onlyab_now)then
	    write(nofile,*)'Only recomputing active terms'
	    else
	    write(nofile,*)'Recomputing both active and reactive terms'
	    end if
	    if(recalc_by_species_now)then
	write(nofile,*)'Recomputing species by species (option recalc_by_species)'
	    else
	write(nofile,*)'Recomputing sum over species (option recalc_total)'
	    end if
        end if

 
c     Number of right-hand sides:
      nrhs = neffan + nscree * nmoscr
c 
c     Preparation of basis functions at Gauss points:
      call bafaga
 
      timsol = 0.
      timfac = 0.
      
        if(gensys)then
        call proini(nthoma, lblock, nrhs, 1, idum)
        else
c       Clean up element and rhs workspace:
        lblt = lblock(1) + lblock(2)
cERN          do j = 1, lcael
cPL          do i = 1, ldael
cPL          ael(i,j) = czero
cPL          end do
cERN          call zset(ldael, czero, ael(1,j), 1)
cERN          end do
          ael(1:ldael,1:lcael) = czero
cERN		do j = 1, nrhs
cPL          do i = 1, ldbel
cPL          bel(i,j) = czero
cPL          end do
cERN          call zset(ldbel, czero, bel(1,j), 1)
cERN          end do
          bel(1:ldbel,1:nrhs) = czero
        end if

      ithoma = 1
      lefnod = .false.       
      write(nofile,*)' '
      write(nofile,*)' '
 
      do 2000 ireg = 1, nreg
c     ======================

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN	14/03/05: If STDVAC = true, always use standard coordinates
c			 in vacuum (regions 2 and 3) (Restore at the end)
	if(ireg>1 .and. STDVAC)cokpco=.false.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write(nofile,7000)ireg
 7000 format('Starting assembly in region ', i2)
      is1 = nsum(ireg-1, ns) + 1
      is2 = is1 + ns(ireg) - 1
      ifirse = ifiel(ireg)
      ilaste = ilael(ireg)
      isubr = is1
      iel = ifiel(ireg)
      if(ireg.gt.1)call regjum(ireg-1, 'R')

      do 2001 isubr = is1, is2
c     ========================
      write(nofile,7001)isubr
 7001 format('Now in subregion ', i2)
      eletyp = styp(isubr)
      ielet = istyp(isubr)
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo
      
      lg = 0
      lm = icolo
      ld = icpblo
      do i = 1, nficom
      it = ibafty(i,ielet)
c       Localization of basis functions in element matrix of one mode pair
c       (i.e. inside 12*12 or 11*11 matrix, resp. for HEC and M23 types)
c       first index of JBFLOC: basis function; second index: field component.
        do ll = 1, ng(it)
        lg = lg + 1
        jbfloc(ll,i) = lg
        end do

        do ll = 1, nm(it)
        lm = lm + 1
        jbfloc(ll+ng(it),i) = lm
        end do

        do ll = 1, nd(it)
        ld = ld + 1
        jbfloc(ll+ng(it)+nm(it),i) = ld
        end do
      end do
      
 
C     Index list for assembly:
        do  i = 1, icolo
        index(i) = i
        jndex(i) = i
        end do
 
      iel1 = i1st(isubr) + 1
      iel2 = i1st(isubr) + iele(isubr)
      do 2002 iel = iel1, iel2
c     ========================

      write(nofile,7002)iel
 7002 format('Element', i4)
c     Element length:
      rle = fl(iel)
c     and a few powers thereof:
        do i = -2, 3
        rlpow(i) = rle ** i
        end do
      rlei = rlpow(-1)

      nmodc = max0(nmode(iel-1), nmode(iel))

      if(nmodc.gt.0)then
c     ..................
 
c     Special care for 1/y terms in first element, only if there is a magn. axis
c     boundary:
      if(crown .or. iel.gt.1)then
      ig1 = 1
      else
      ig1 = 0

c     To use in case sing terms had a special treatment for iel>1:
C      CALL BAFAD(-FX0(IEL-1)*RLEI, BAFG(1,0))
C      CALL BAFADN(it, -FX0(IEL-1)*RLEI, BAFGN(1,0,it,0), BAFSIN(1,0,it,0))
c     Now only for iel=1 at axis, dealt with in call to bafaga. 
c     Hence next 2 lines removed:
C      CALL BAFAD(0.d0, BAFG(1,0))
C      CALL BAFADN(it, 0.d0, BAFGN(1,0,it,0), BAFSIN(1,0,it,0))
      end if
 
      intab = 1 + (iel-1) * (ngauss+1) + ireg

cERN  NEW option for kappa profile (force kappa = 1 near magnetic axis) 
        y = abscno(intab)
	if(dshape .and. kapprof) kappa = kappa_orig - (kappa_orig-1) * exp(-50*y*y)
cERN  -----------------------------------------------------------------

        if(circ)then
c       Rotation to magnetic (+-//) axes at Gauss points:
c       OK for circular; depends on poloidal angle in general geometry; in this case, dealt with in FOUCOG, VACMBG.
        call loctra(iel, .true.)
        end if
 
c     Filling of LDIF:
      call dset(72*(ngauss+1), 0.d0, ldif1, 1)
        do ic = 1, nficom
        it = ibafty(ic,ielet)
          do ib = 1, nbafel(ic,ielet)
          jb = jbfloc(ib,ic)
          id = ideriv(ib,it)
          ri = rlpow(id)
          ri1 = rlpow(id - 1)
c           'Field' degrees of freedom:
            do ig = ig1, ngauss
            ldif(2*ic-1, jb, ig) = ri * bafgn(ib,0,it,ig)
            end do
c           'Radial derivative' degrees of freedom:
            if(eletyp.ne.'M23' .or. ic.ne.1 .or. flrops(ireg))then
c           Special case of curl.curl and no FLR plasma terms: 
c           Erho' never contributes
              do ig = ig1, ngauss
              ldif(2*ic, jb, ig) = ri1 * bafgn(ib,1,it,ig)
              end do
            end if
          end do
        end do


cccccccccccccccccccccc ERNESTO cccccccccccccccccccccccc

c     Assemble curl.curl and displacement current contribution:

      write(*,"(a13,i4,a17,g10.4,a2)") '     element#',iel, 
     ;            ' ( node radius = ', abscis(intab),' )'
      write(605,"(a13,i4,a17,g10.4,a2)") '     element#',iel, 
     ;            ' ( node radius = ', abscis(intab),' )'
c	  write(*,"a20,g10.4,a2") '    ( node radius = ', abscis(intab),' )'
c	  print *
	  tempo0 = second()-timin
      if(gensys .or. .not.onlyab_now)then
	  write(605,*), '    --> (GENERA3.f) Begin ASCURL:', tempo0, 'sec' 
c	  print *, kappa
	call ascurl
	  tempo  = second()-timin
          write(605,*), '    --> (GENERA3.f)   End ASCURL:', tempo, 'sec' 
	  write(605,*), '         Time used:', tempo-tempo0, 'sec'
cJAC	  print *
	end if

c     Assemble plasma contribution:
	  
      if(.not. vacuum(ireg))then  ! Only process plasma regions

        if(recalc_by_species_now)then  ! Postprocessing: element powers by species
cJAC 	  tempo0 = etime (T1)
cJAC	  print *, '    --> (GENERA3.f) Begin ASPLAS + QUFETC2' // 
cJAC     ;           ' with species loop:', tempo0, 'sec'
          do isp = 1, nspec
c         assemble for 1 species at a time:
cERN	    call asplasnew(.true., isp)
	    call asplas(.true., isp)

c         get plasma active and if required reactive power for each species and clean up ael, bel
	    call qufetc2(.true., ithoma, .false., .true., ithoma.eq.1 .and. isp.eq.1, isp)
  	    end do
cJAC 	  tempo = etime (T1)
cJAC        print *, '    --> (GENERA3.f)   End ASPLAS + QUFETC2' // 
cJAC     ;           ' with species loop:', tempo, 'sec'
cJAC	  print *, '         Time used:', tempo-tempo0, 'sec'
cJAC	  print *

        else  ! Usual treatment:
	  tempo0 = second()-timin
	  write(605,*) '    --> (GENERA3.f) Begin ASPLAS:', tempo0, 'sec' 
cERN	  call asplasnew(.false., idum)
	  call asplas(.false., idum)
	  tempo  = second()-timin
          write(605,*), '    --> (GENERA3.f)   End ASPLAS:', tempo, 'sec' 
	  write(605,*), '         Time used:', tempo-tempo0, 'sec'
cJAC	  print *
        end if
	end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

      if((gensys .or. .not.rawel1) .and. iel.eq.1)then
c     ================================================
c     First element at magnetic axis or crown inner boudary.
c     Note for hypothetic automatic meshing: nmode(0) = nmode(1) assumed!

c     Boundary conditions: 
        if(.not.bcinb1)then
c       Remove matrix singularity at magnetic axis by setting '1.' on diagonal
c       for non-assembled or eliminated dof:
        call filmao(ael, ldael)

        else
c       Build axis constraints in CBC, their conj.transposed in CBCTC:
        call axisbc(cbc, maxbll, lbl1, lbl2)
        write(nofile,*)'Axis constraints: ', lbl1, ' equations connecting max.', lbl2, ' columns'
          do i = 1, lbl1
            do j = 1, lbl2
            cbctc(j,i) = dconjg(cbc(i,j))
            end do
          end do
c        write(nofile,*)'cbc:'
c          do i = 1, lbl1
c	    write(nofile,*),i
c          write(nofile,1001)(cbc(i,j), j = 1, lbl2)
c          end do

C        write(nofile,*)'cbctc:'
C          do i = 1, lbl1
C          write(nofile,1001)(cbctc(j,i), j = 1, lbl2)
C          write(nofile,*)' '
C          end do
 1001     format(11(1h(,g11.3,1h,,g11.3,1h)))
c       Partial factoriz. of CBCTC with row pivoting restricted to first block;
ccc       This disables scaled pivoting (flaw in CGEFA3 for zero row...):
cc        call dset(LBLOCK(1), 1.d0, sc1, 1)
        call cgefa3(cbctc, 2*maxbll, lbl2, lbl1, lblock(1), ipvt, sc1, 0, info)
        if(info .ne. 0)write(nofile,*)'cgefa3 returned info =', info
C        write(nofile,*)'factored cbctc:'
C          do i = 1, lbl1
C          write(nofile,1001)(cbctc(j,i), j = 1, lbl2)
C          write(nofile,*)' '
C          end do
c 1002     format(12(1h(,g11.3,1h,,g11.3,1h)))
c      write(nofile,*)'ael block before elim'
c	do i = 1, lblock(1)+lblock(2)
c	write(nofile,1002)(ael(i,j),j=1,lblock(1)+lblock(2))
c	end do
c       Eliminate Lagrange multipliers, assuming no source terms in first elt.:
c       a) Gauss elimination on first elt. block row, with multipliers in CBCTC:
          do k = 1, lbl1
          kp1 = k + 1
          l = ipvt(k)
            do j = 1, lblock(1) + lblock(2)
            t = ael(l,j)
              if (l .ne. k) then
              ael(l,j) = ael(k,j)
              ael(k,j) = t
              end if
            call zaxpy(lbl2-k, t, cbctc(kp1,k), 1, ael(kp1,j), 1)
            end do
          end do
c      write(nofile,*)'ael block after elim'
c	do i = 1, lblock(1)+lblock(2)
c	write(nofile,1002)(ael(i,j),j=1,lblock(1)+lblock(2))
c	end do
c       b) Replace first rows of AEL by constraints; clean up:
          do i = 1, lbl1
            do j = lbl2 + 1, lblock(1) + lblock(2)
            ael(i,j) = czero
            end do
            do j = 1, lbl2
            ael(i,j) = cbc(i,j)
C           Important: these are equivalenced to solver arrays, must be clean:
            cbc(i,j) = czero
            cbctc(j,i) = czero
            end do
          end do
c     write(nofile,*)'ael block after replacement of first rows'
c	do i = 1, lblock(1)+lblock(2)
c	write(nofile,1002)(ael(i,j),j=1,lblock(1)+lblock(2))
c	end do
c       First element is now ready for FBT algorithm to start.
        end if

c     End of special code for first element.
      end if
C     ======
 
c     Current element is built. Now associated r.h.s. contribution:
      call filrhs(.false.)

C     IF(IREG.GE.NREG-1)THEN
C     WRITE(6,*)'ELT ',IEL,'; BLOCK ',ITHOMA
C     DO I = 1, LBLOCK(ITHOMA)+LBLOCK(ITHOMA+1)
C     WRITE(6,*)I,BEL(I,1)
C     END DO
C     END IF

C     [Here, coud add condensation of bubble unknowns (M23 elements).
C     Need add switch in data, open new AQIO file;
C     then final NADD2 = ICOLO * NMODE(IEL-1)]

C     Now call block factorization:
      nadd2 = ibulo * nmodc + icolo * nmode(iel-1)
      if(nadd2.gt.0)then

        if(gensys)then
c       --------------
          if(scale)then
          isca = 2
C   
C         Build scaling factors:
C   
          lbl1 = 0
          lbl2 = lblock(ithoma)
            if(iel.gt.1)then
            lbl1 = lblock(ithoma-1)
            call dcopy(lbl1,scaf(1),1,scaf(1-lbl1),1)
            end if
     
            do i = 1, lbl2
            rtem = cabs1(ael(i,i+lbl1*modele))
              if( rtem .eq. 0.d0 )then
              scaf(i) = 1.d0
              else
c              scaf(i) = 1.d0 / sqrt(rtem)
              scaf(i) = dsqrt(rtem)
              end if
            end do
c        do i = 1, lbl2
c          do j = 1-modele*lbl1, lbl2
c          rtra(j) = cabs1(ael(i,j+modele*lbl1))*scaf(j)
c          end do
c          if(modele.eq.0 .and. ithoma.gt.1)then
c            do j = 1-lbl1, 0
c            rtra(j) = cabs1(btom(i,j+lbl1))*scaf(j)
c            end do
c          end if
c        rownor(i) = scaf(i) * rtra( ismax(lbl1+lbl2,rtra(1-lbl1),1)
c       ;-lbl1)
c        end do

          else
          isca = 0
          end if
         

c       Perform one step in solver:
	  stiref = second()
ccccccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccc
  	  tempo0 = second()-timin
	  write(605,*), '    --> (GENERA3.f) Begin FBSTEP:', tempo0, 'sec' 
        call fbstep(isca, scaf(1), scaf(1))
	  tempo  = second()-timin
          write(605,*)'    --> (GENERA3.f)   End FBSTEP:', tempo, 'sec' 
	  write(605,*), '         Time used:', tempo-tempo0, 'sec'
cJAC	  print *
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	  
	  timfac = timfac + second() - stiref

        else if(.not.recalc_by_species_now)then  ! .not.gensys, case total power per element
c       ---------------------------------------
c       Case routine is called after solution of linear system:
c       Compute element quadratic form, residual, ...:
        call qufetc2(.true., ithoma, .false., .true., ithoma.eq.1, idum)
        end if
c       ------

      ithoma = ithoma + 1
      write(603,*)'                    ithoma=',ithoma

      end if

        if(.not.giplas .and. iel.gt.1 .and. .not.vacuum(ireg) .and. iel.ne.ilael(ireg))then
c       Under development (active only when data giplas=.false., i.e. analytical integration of products of basis functions):
c       Use results of right node to assemble left node in next element.
c       Not ready yet for the case recalc_by_species_now=T!
        lefnod = .true.
        call passem
        lefnod = .false.       
        end if 

      end if
c     ......

 2002 continue

 2001 continue

      isubr = is2
      iel = ilael(ireg)
ccccccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccc
      if(nmode(iel).gt.0)then
	  tempo0 = second()-timin
          print*,'Enter REGJUM'
	  write(605,*), '    --> (GENERA3.f) Begin REGJUM:', tempo0, 'sec' 
	call regjum(ireg,'L')
	  tempo  = second()-timin
          print*,'Exit REGJUM'
          write(605,*), '    --> (GENERA3.f)   End REGJUM:', tempo, 'sec' 
	  write(605,*), '         Time used:', tempo-tempo0, 'sec'
c	  print *
	end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

 2000 continue
c     ========
c    ==========
c     idcat = idcat + nmode(nele) * iconn(isubr)

      if(nmode(nele).gt.0)then
        if(gensys)then
	  stiref = second()
	  call fbstep(0, rownor, rownor)
	  timfac = timfac + second() - stiref
	  end if
      ithoma = ithoma + 1
      end if

      if(gensys)then
      write(nofile,*)'system solved; time= ',second()-timin
      write(nofile,*)'time spent in factorization routine= ',timfac
 
C     Fill all components of solution at inner boundary 
C     (magnetic axis or metal shell), when required:
      if(.not.bcinb1)call filsol
      
c      if(wrisol)then
c      write(nofile,*)' '
c      write(nofile,*)'Solution vectors:'
c      i2 = 0
c        do i = 1, nthoma
c        CALL RDSOL(IMOTO, i, BEL, LDBEL, LBL, NRHSS)
c        do i1 = 1, LBL
c        i2 = i2 + 1
c        write(nofile,1000)i2,(BEL(i1,j),j=1,nrhss)
c        end do
c        end do
c      end if
      end if  ! gensys

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN	14/03/05: Restore original cokpco value
	cokpco = cokpco_orig
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
c 1000 format(1h ,i4,(2x,g13.5))
      end
