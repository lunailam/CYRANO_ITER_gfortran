      SUBROUTINE DIRESP_active(isp)

	use dfport		! necessary for 'etime' function

      IMPLICIT NONE
      integer isp

C     Dielectric response of a species, using equilibrium of QLFP code
C     interpolated at local radial index INTAB of Cyrano by routine INTF0
C     isp is index of species among the ones requiring general diel. response
C     Absolute species index is ISPE																						

      include 'PARDIM.COPY'
      include 'COMREG.COPY'
      include 'COMSUB.COPY'
      include 'COMGEO.COPY'
      include 'COMEQU.COPY'
      include 'COMANT.COPY'
      include 'COMFIN.COPY'
      include 'COMFIC.COPY'
      include 'COMIN2.COPY'
      include 'COMROT.COPY'
      include 'COMMAG.COPY'
      include 'COMMA2.COPY'
      include 'COMPLA.COPY'
      include 'COMFOU.COPY'
      include 'COMMOD.COPY'
      include 'COMMMK.COPY'
      include 'COMPRI.COPY'
      include 'COMSWE.COPY'
      include 'COMGDR.COPY'
      include 'COMUSR2.F'
      include 'COMPHY.COPY'
      include 'COKPCO.COPY'
	 
      logical pass, cxr, toini

      save toini
      
      integer i, j, im, ix, jv, ip, iex, jev, nset0, isig, igv, igx
     ;, igv2, igx2, jev1, jev2, iex1, iex2
     ;, ispe, ii, irv, jvn, irx, ixn, icount, ic, ipov, ipox, ih
     ;, ji1, jm1, i1, im1, ib, jb
     ;, iplb, kstop1, kstop2, L, La
     ;, isrchfgt

C     act, react: im and re of iL      
c     See if they could be equivalenced to, say, VMAT in comusr2,
c     or other memory saving?
      double precision
     ;  act(2*maxpom-1,0:maxcou,-1:1)
     ;, react(2*maxpom-1,0:maxcou,-1:1)
     ;, ac(2*maxpom-1,0:maxcou,-1:1,3)
     ;, reac(2*maxpom-1,0:maxcou,-1:1,3)
      equivalence (act,ac), (react,reac)

      double precision 
     ;  v, v2, vi, xn, xna, x
      double precision 
     ;  u0, u0i, u1, u2, ru1, ru2, iu1, iu2, xtir1, khi
      double precision 
     ;  khi1, lamb, del2, del, te0, t1, t2, t3, t4, t5
     ;, cmr1, cmr2, smr1, smr12, smr2, smr22, cmr0
     ;, smr0, smr02, csc22
     ;, clmr0(0:maxcou), clmr1(0:maxcou), clmr2(0:maxcou)
     ;, er0, r00, r1, r2, r3, s0, s1, s1i, s2, s2i, omxnai, xmax
     ;, hp, hp2, hp3, hp4, hp5, hp6, omca, xdb, xdbi, sxdb, sxdbi
     ;, petita(5*maxnex), petitb(5*maxnex)
     ;, nze, N1, L0, L1, jze, J1, alp2, alt2
     ;, JL(-1:maxcou+3,5), LL(-1:maxcou+3,5) !, NL(5*maxnex,-1:maxcou+3)
     ;, tj1(maxcou+2), tj2(maxcou+2), tj3(maxcou+2), tj4(maxcou+2)
     ;, cmmv(5*maxnex)
     ;, kel, kelv(5*maxnex), s2k
     ;, maomi
     ;, ek(5*maxnex), ee(5*maxnex), cmm, w0, w1, w2
     ;, oo, o1, o2, o3, o4, o5, fis
     ;, bn(4,4)
     ;, ellfc, ellec, ellpic

      complex*16
     ;  JLc(-1:maxcou+1), LLc(-1:maxcou+1), cu1, cdellpic
     ;, ccmr1, csmr12
     
cERN	Some Variables added ccccccccccccccccccccccccccccccccccccccccccccc
	character(100) :: FILE_NAME
	character(4) :: panam4
	character(10) :: polname(2)
	integer :: OpenStat, NMM
	integer :: mdiffaux(4*maxcou+2), mdiffvector(2*maxcou+1)
	character(100) :: GGs, GG2s, GGaux  
	character(2) :: charaux, charaux2
	real*8 :: cosXo, Acoef(2*maxpom-1), aux1, aux2, aux3, omcdel, TESTF0, TESTF1
	real*8 :: actout1, actout2, sinout1, sinout2, X1, X2

	integer, parameter :: Nchi = 128	! Number of chi-harmonics 
c	complex*16 :: Pell(Nchi+51,klim+1)	! P_ell coeficient
	real*8 :: Pell(Nchi+51,klim+1)	    ! P_ell coeficient
	complex*16 :: PRODU(klim+1)			! Product (P_ell x T_ell)
	integer :: Ell_min(klim+1), Ell_max(klim+1), mdiff(klim+1)
	logical only_act 
	real :: TT(2)	            !   Variables for evaluating
	character(4) pana4
	real*8 :: AAleft(2*maxpom-1), x_tang, 
     ;          xn_tang_pos(2*maxpom-1), xn_tang_neg(2*maxpom-1)
	integer :: indx_k(2*maxpom-1)
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      equivalence (JLc, JL), (LLc, LL)

      external ellfc, ellec, ellpic, cdellpic

      data toini/.true./

cERN	Number of k// values (also given by NA = 2*nmoant-1)	
	NA = 2*(modva2-modva1)+1

cERN	 Reset Auxiliary variables
c	aux1 = 0.0d0
c	aux2 = 0.0d0
c	aux3 = 0.0d0
	m12left   = czero
	m12right  = czero
	m12landau = czero


c       Coefficients for NL recurrences:
        if(toini)then
          do j = 1, maxcou+2
          tj1(j) = dfloat(2*j) / (dfloat(j)-0.5d0)
          tj2(j) = (dfloat(j)+0.5d0) / (dfloat(j)-0.5d0)
          tj3(j) = dfloat(2*j) / (dfloat(j)+0.5d0)
          tj4(j) = (dfloat(j)-0.5d0) / (dfloat(j)+0.5d0)
          end do
        toini = .false.
        end if


      ispe = ispgdr(isp)
      
      nset0 = (2*maxpom-1)*(maxcou+1)*3
      call dset(3*nset0, 0.d0, ac, 1)      
      call dset(3*nset0, 0.d0, reac, 1)
      call dset(nset0, 0.d0, tact, 1)      
      call dset(nset0, 0.d0, treact, 1)
      clmr0(0) = 1.d0
      clmr1(0) = 1.d0
      clmr2(0) = 1.d0

c     Normalized x at separatrix (co-passing):
      xnsep = - 2.d0 * delb(intab) / bmax(intab) 
c     Indices of left and right limits of x intervals:
c     Resp. co-passing, 'co-trapped', 'counter-trapped', counter-passing
        if(gigdr)then
c       For Gauss points loop
        igl(1) = 1
        igl(2) = isrchfgt(lxg, xngaug, 1, xnsep)
        igl(3) = isrchfgt(lxg, xngaug, 1, 0.d0)
        igl(4) = isrchfgt(lxg, xngaug, 1, - xnsep)
        idl(4) = lxg
        idl(1) = igl(2) - 1
        idl(2) = igl(3) - 1
        idl(3) = igl(4) - 1
        igv2 = ngauv
        igx2 = ngaux
        else
c       For element nodes loop
        call icopy(nxgreg, ifiax, 1, igl, 1)
        call icopy(nxgreg, ilaax, 1, idl, 1)
        igv2 = 1
        igx2 = 1
        end if

c	Radius dependent quantities
      omca = qom(ispe) * B0
	omcdel = qom(ispe) * delb(intab)
      xtir1 = omca / omegag
      oo = 1.d0 / dabs(xtir1)
      te0 =  PI / (qom(ispe)*hachi(intab)*delb(intab))
      s0 = bbar(intab) / delb(intab)
      xmax = B0 / bmin(intab)
      maomi = bmax(intab) / bmin(intab)
      s1 = B0 / delb(intab)
      s1i = 1.d0 / s1
      s2 = bmin(intab) / delb(intab)
      s2i = delb(intab) / bmin(intab)
c	write(666,*), abscis(intab)
c	write(666,*), xmax
c	write(666,*), maomi

c     Normaliz. of radial response density:
      r00 = PI * (EEL*ZCH(ispe))**2
     ; * 4.d0 * dPsidr_n(intab) * abscis(intab)
     ; / (2.d0*MH*AMASS(ispe) * B0 * hachi(intab) * omegag)


cERN	10/05/05 OUTPUT for radial quantities ccccccccccccccccccccccccccccc
	if (WRITE_OUTPUT) then
         FILE_NAME = COKFOLDER  // "/radial_data.dat" 
			open (UNIT = 99, FILE = FILE_NAME, STATUS = "REPLACE",
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
		  
				  write(99,*),"Some rho-dependent quantities"
				  write(99,*) ' rho=   ', abscis(intab)  
				  write(99,*) ' wRF=   ', omegag
				  write(99,*) ' q/m=   ', qom(ispe)
				  write(99,*) ' Bo=    ', B0
				  write(99,*) ' delB=  ', delb(intab)
				  write(99,*) ' Bbar=  ', bbar(intab)
				  write(99,*) ' wc=    ', omca 
				  write(99,*) ' wcdel= ', omcdel
				  write(99,*) ' wcbar= ', qom(ispe) * bbar(intab)

				  write(99,*) ' Xmax=  ', xmax
				  write(99,*) ' xnsep= ', xnsep
				  write(99,*) ' xtir1= ', xtir1
				  write(99,*) ' r00=   ', r00
				  
				  write(99,*) ' Nk//=   ', NA
					do j = 1,NA
						write(99,*), allkpa(j)
					end do

			close(99)

	end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      write(6,*)'diresp: before generate NL passing'    
c	print*,'diresp: before generate NL:', etime (TT) 
	  
c     Generate the NORMALIZED NL on x mesh:
c     Passing regions:
      do ii = 1, 4, 3
        do i = igl(ii), idl(ii)
cERN
        xn = xngaug(i)

cERN	Avoid exact xn=xnsep
	if(ii.eq.1 .and. i.eq.idl(ii))xn = xngaug(i)-0.002d0
 	if(ii.eq.4 .and. i.eq.igl(ii))xn = xngaug(i)+0.002d0
cERN	Avoid exact xn=0
        xna = dmax1(dabs(xn),1e-6)
        petita(i) = dsqrt(xna)
        petitb(i) = dsqrt(maomi*(xna+xnsep))
        kel = dsqrt(1.d0-(petitb(i)/petita(i))**2)
        kelv(i) = kel
        call ellfec(kel, ek(i), ee(i))
c       Normalized NL:
        nze = 2.d0 / petita(i) * ek(i)
cERN        if(kel.ne.0.d0)then
	if(kel.gt.1.d-5)then
        cmm = - 0.5d0 * (kel + 1.d0 / kel)
        cmmv(i) = cmm
        s2k = dsqrt(2.d0 / kel)
        N1 = 2.d0 / petita(i) * (ek(i) +2.d0/(kel*kel)*(ee(i) - ek(i)))
        else
        N1 = 0.d0
c NB:   actual value is - infinity:
        cmmv(i) = 0.d0
        end if
c       backward recurrence is stable for passing: check starting index for
c       overall accuracy!
        NL(i,maxcou+3) = 0.d0
        NL(i,maxcou+2) = 1.d0
          do j = maxcou+2, 1, -1
          NL(i,j-1) = tj1(j)*cmm * NL(i,j) - tj2(j) * NL(i,j+1)
          end do
          do j = 1, maxcou+3
          NL(i,j) = NL(i,j) * (nze / NL(i,0))
          end do
        NL(i,0) = nze
        NL(i,-1) = NL(i,1)
        end do
      end do

      write(6,*)'diresp: before generate NL trapped'      
c     Trapped regions:
      do ii = 2, 3
        do i = igl(ii), idl(ii)
cERN
        xn = xngaug(i)

cERN	Avoid exact xn=xnsep
	if(ii.eq.2 .and. i.eq.igl(ii))xn = xngaug(i)+0.002d0
	if(ii.eq.3 .and. i.eq.idl(ii))xn = xngaug(i)-0.002d0
cERN	Avoid exact xn=0
	if(ii.eq.2 .and. i.eq.idl(ii))xn = xngaug(i)-0.002d0
 	if(ii.eq.3 .and. i.eq.igl(ii))xn = xngaug(i)+0.002d0

        xna = dabs(xn)
        petita(i) = dsqrt(xna)
c       Convention, actual b is complex for trapped:
        petitb(i) = 0.d0
        omxnai = 1.d0 / (1.d0 - xna)
        cmm = (1.d0 - xna * s0) * omxnai
        cmmv(i) = cmm
        kel = dsqrt(0.5d0 * (1.d0 - cmm))
        kelv(i) = kel
        call ellfec(kel, ek(i), ee(i))
        x = (1.d0 - dabs(xn)) * xmax
        sxdbi = dsqrt(s1 / x)
c       Normalized NL:
        NL(i,0) = 2.d0 * kel * ek(i) / petita(i)
        NL(i,1) = 2.d0 * kel * (2.d0*ek(i) - ee(i)) / petita(i)
        NL(i,-1) = NL(i,1)
          do j = 1, maxcou
          NL(i,j+1) = tj3(j)*cmm*NL(i,j) - tj4(j)*NL(i,j-1)
          end do
        end do
      end do


cERN	Write N_ell series to FILE in radial folder COKFOLDER ccccc
	if(WRITE_OUTPUT)then
		FILE_NAME = COKFOLDER  // "/Nell.dat"
		open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
     ;		  IOSTAT = OpenStat, ACTION = "WRITE")
			  if (OpenStat > 0) then
				  print *, 'Error writing file in DIRESP: ', FILE_NAME
				  stop
			  end if
			  write(13,*) ' NL series'
			  write(13,*) ' Nx ', ' Nchi '  
			  write(13,*) lxg, maxcou+2
			  write(13,*) ' radius ', ' xn_sep ', ' X(max) '  
			  write(13,*) abscis(intab), xnsep, xmax
			  do i = 1, lxg    
				 write(13,2222)xngaug(i), cmmv(i), (NL(i,j),j=0,maxcou+1)     
			  end do
		close(13)
	end if
c	ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cERN	NEW 06/05/05 
c	(v,x) grid OUTPUT

	if (WRITE_OUTPUT) then
         FILE_NAME = COKFOLDER  // "/VXloop.dat" 
			open (UNIT = 66, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if	  
				  write(66,*),"v,xn,k//,AC,... for Non-Maxwellians (DIRESP)"
				  write(66,*),nxgreg, nvgreg
				  write(66,*),(ielrv(i),i=1,nvgreg)
				  write(66,*),(ielrx(i),i=1,nxgreg)
				  write(66,*),NA
				  write(66,*) ' radius ', ' xn_sep ', ' X(max) '  
			      write(66,*) abscis(intab), xnsep, xmax
				  write(66,*), '    '
				  write(66,*) ' v ', ' xn ', ' F0i '  

cERN	Tangent resonance oputput
         FILE_NAME = COKFOLDER  // "/tang_res.dat" 
			open (UNIT = 123, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if  
				  write(123,*),"Tangent resonance location xn (DIRESP)"
					  write(123,*),lvg,NA
				  write(123,*) ' radius ', ' xn_sep ', ' X(max) '  
				  write(123,*) abscis(intab), xnsep, xmax
				  write(123,"(200G16.3)") 0.0d0, allkpa(1:NA)
c				  write(123,*) ' v ', ' xn ', ' x '  

	end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	TESTF0 = 0.d0
	TESTF1 = 0.d0


c	print*,'diresp: before 2D loop:', etime (TT) 
      write(6,*)'DIRESP: before 2D loop'  
c     V region loop:
      do irv = 1, nvgreg
        if(gigdr)then
C       jev is an element index
        jev1 = ifielv(irv)
        jev2 = ilaelv(irv)
        else
C       jev is a node index; jv = jev
        jev1 = ifiav(irv)
        jev2 = ilaav(irv)
        end if

c     X region loop:
      do irx = 1, nxgreg

C       This assumes 4 x regions!:
        if(irx.le.2)then
        isig = 1
        else 
        isig = -1
        end if
      fis = dfloat(isig)
      pass = irx.eq.1 .or. irx.eq.4
        if(gigdr)then
C       iex is an element index
        iex1 = ifielx(irx)
        iex2 = ilaelx(irx)
        else
C       iex is a node index; ix = iex
        iex1 = ifiax(irx)
        iex2 = ilaax(irx)
        end if


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN	10/05/05: BUG here: vlelr and xlelr are region lengths!
cERN			 We need the element lengths instead! (assumed constant in each region)
cERN  call normat(vlelr(irv), xlelr(irx))
	call normat(vlel(jev1), xlel(iex1))
      call intnod
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      
c     V element loop:
c     when gigdr=.F., loop over nodes instead
      do jev = jev1, jev2
c     ===================
C     right node index (gigdr=.t.):
      jvn = irv + jev

c     X element loop:
c     when gigdr=.F., loop over nodes instead
      do iex = iex1, iex2
c     -------------------
C     right node index (gigdr=.t.):
      ixn = irx + iex
C     f0 and first derivatives at Gauss points:
      icount = 0
      if(gigdr)call fagp(jev, iex, jvn, ixn, isp)

c     V Gauss point loop:
c     when gigdr=.F., igv2 = 1  and jv = jev    
      do igv = 1, igv2
C     ================
      jv = (jev - 1) * igv2 + igv
      v = vgaug(jv)
      v2 = v * v
      vi = 1.d0 / v
      er0 = r00 * v2
      r1 = er0 * v2 * PI

cERN	24/05/05: New procedure to treat the tangential resonances ccccccccccc

c	Reset variables
	indx_k(1:NA) = 0
	x_tang = 0.0d0
	xn_tang_pos(1:NA) =  0.0d0
	xn_tang_neg(1:NA) = 0.0d0

c	Big A from paper (p = +1)
	AAleft(1:NA) = (omegag-qom(ispe)*bmin(intab)) / (allkpa(1:NA)*v);

	do j = 1,2*nmoant-1	! k// - loop
		if(abs(AAleft(j)) < 1.0d0)then
			indx_k(j) = 1	! FLAG indicating tangential resonance for given (k//,v)
			x_tang = xtir1 * (1.d0-AAleft(j)**2)/(1.d0-allkpa(j)*v*AAleft(j)/omegag)
			xn_tang_pos(j) =  1-x_tang/xmax
			xn_tang_neg(j) = -1+x_tang/xmax
c			Find indices in xn abcissa
c			to do! ###################
		end if
	end do

c	xn_tang output
	if(iex.eq.1)then
		write(123,"(200G16.3)") v, xn_tang_neg(1:NA)    
	end if

cERN	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     X Gauss point loop:
c     when gigdr=.F., igx2 = 1  and ix = iex

        do igx = 1, igx2
C       ~~~~~~~~~~~~~~~~
        ix = (iex - 1) * igx2 + igx
        xn = xngaug(ix)


cERN	Avoid xn=xnsep
	if(xn.eq.xnsep .and. irx.eq.1) xn = xnsep-0.0002d0
	if(xn.eq.xnsep .and. irx.eq.2) xn = xnsep+0.0002d0
	if(xn.eq.-xnsep .and. irx.eq.3) xn = -xnsep-0.0002d0
	if(xn.eq.-xnsep .and. irx.eq.4) xn = -xnsep+0.0002d0

cERN	Avoid xn=0
	if(xn.eq.0.d0 .and. irx.eq.2) xn = -0.002d0
	if(xn.eq.0.d0 .and. irx.eq.3) xn = +0.002d0

        x = (1.d0 - dabs(xn)) * xmax
  
        khi1 = x / xtir1
C       x*delta/B0:
        xdb = x / s1
        xdbi = 1.d0 / xdb
        sxdb = dsqrt(xdb)
        sxdbi = 1.d0 / sxdb



cERN	24/05/05 deal with tangent resonance (to do!) 
c	if(irx .eq. 1 .or. irx .eq. 2)then	! xn_tang_neg 
c	if(xn > xn_tang_neg) then ! current element is singular
c	  indx_tg_neg(jv) = ix	! find index of current element
c	end if

c	else !(x-regs 3 and 4) ! xn_tang_pos 
	
c	if(xn>xn_tang_pos) then
	! special element
	! previous element is singular: special treatment
c	indx_tg_pos(jv) = ix-1 ! find index of previous element
c	end if

c	end if


cERN	VX output
	write(66,*),'     '
	write(66,"(4G16.8)")v,xn,F0i(ix, jv, 1, isp),F0i(ix, jv, 2, isp)



c	write(5656,*), "  v,  vlelr,  vgaug(j+1)-vgaug(j)"
c	write(5656,"(4G16.8)") v,vlel(jev1) ,vgaug(jv+1)-vgaug(jv)

c	write(4646,*), "  xn,  xlelr,  xngaug(i+1)-xngaug(i)"
c	write(4646,"(4G16.8)") xn,xlel(iex1) ,xngaug(ix+1)-xngaug(ix)



c	print *, vlelr(irv)* xlelr(irx)
c	print*, jv,ix        
c	print*, jv,ix


c     Find resonances:
c loop over ip=±1 to add below!?
cERN      ip = 1
cERN      khi = ip * khi1


c     Fragment stored in GIGDR.f: Gaussian integration in x and v, to update if
c     ever needed!
c     INCLUDE 'GIGDR.f'



      if(.not.gigdr)then
c     Semi analytical integration over elements
C     =========================================
c     Groupe selon coef. de f0, df0/dv, df0/dx (cf dernier indice ac et reac)



      r3 = er0 * PI * s1 * 0.25d0 * omegag / dabs(omca)


            
      do im = 1, 2*nmoant-1
C     >>>>>>>>>>>>>>>>>>>>>

c     p=0:
      u0i = (allkpa(im) * v) / omegag
      u0 = 1.d0 / u0i
      cmr0 = s0 - xdbi * (1.d0 - u0 * u0)
      smr02 = 1.d0 - cmr0 * cmr0
      csc22 = 2.d0 / (1.d0 - cmr0)
C       p=0 only absorbs if isig*k//>0, |u0i|>1, and x between...
c        if(isig .eq. dint(dsign(1.d0,u0i)) .and. dabs(u0i).gt.1.d0)then
        if(u0*fis.gt.petitb(ix) .and. u0*fis.lt.petita(ix))then
c	if(abs(u0)<1 .or. abs(petita(ix))>1)then
c	print *, petita(ix), petitb(ix), u0
c      end if
c	  if(u0.gt.petitb(ix) .and. u0.lt.petita(ix))then

        smr0 = dsqrt(smr02)
        clmr0(1) = cmr0
c      r2 = r1 * xdbi * hp
	print *, 'Active p=0'
c	Active p=0:
        t3 = er0 * PI * xdbi * dabs(u0)**3 / smr0
        ac(im,0,0,2) = - v * t3
        ac(im,0,0,3) =  2.d0 * x * t3
        ac(im,1,0,2) = cmr0 * ac(im,0,0,2)
        ac(im,1,0,3) = cmr0 * ac(im,0,0,3)
          do i = 2, maxcou
          clmr0(i) = 2.d0 * cmr0 * clmr0(i-1) - clmr0(i-2)
          ac(im,i,0,2) = ac(im,0,0,2) * clmr0(i)
          ac(im,i,0,3) = ac(im,0,0,3) * clmr0(i)
          end do
        end if


c     p=±1: 2 roots to resonance equations.
      do ip = -1, 1, 2
c     ````````````````
      khi = dfloat(ip) * khi1
      lamb = 0.5 * u0i * khi
      del2 = lamb * lamb + 1.d0 - khi
      cxr = del2 .lt. 0.d0		! complex root

        if(cxr)then

cERN		reactive calculation removed here!		

        else

        del = dsqrt(del2)
        u1 = lamb + del
        u2 = lamb - del
        t1 = dfloat(ip) / (qom(ispe) * delb(intab))
        t2 = omegag - dfloat(ip) * qom(ispe) * bbar(intab)
        cmr1 = t1 * (allkpa(im) * v * u1 - t2)
        cmr2 = t1 * (allkpa(im) * v * u2 - t2)

c       ac(im,0,-1:1)
c see: two roots can take place together

		actout1 = 0.0
		sinout1 = 0.0
	      X1 = 0.0
cERN          if(dabs(cmr1) .le. 1.d0 .and. x < 1.d0)then
		if(dabs(cmr1) .le. 1.d0)then
		X1 = dacos(cmr1)
            smr12 = 1.d0 - cmr1 * cmr1
          smr1 = dsqrt(smr12)
		sinout1 = smr1 
          clmr1(1) = cmr1
          t3 = r3 * (1.d0 - u1 * u1) / (2*del * smr1)
          t4 = v * t3
		actout1 = -t4
          t5 = t3 * 2.d0 * (x - ip * xtir1)
          ac(im,0,ip,2) = ac(im,0,ip,2) - t4
          ac(im,1,ip,2) = ac(im,1,ip,2) - t4 * cmr1
          ac(im,0,ip,3) = ac(im,0,ip,3) + t5
          ac(im,1,ip,3) = ac(im,1,ip,3) + t5 * cmr1
            do i = 2, maxcou
            clmr1(i) = 2.d0 * cmr1 * clmr1(i-1) - clmr1(i-2)
            ac(im,i,ip,2) = ac(im,i,ip,2) - t4 * clmr1(i)
            ac(im,i,ip,3) = ac(im,i,ip,3) + t5 * clmr1(i)
            end do
          end if

		actout2 = 0.0
		sinout2 = 0.0
	      X2=0.0
cERN          if(dabs(cmr2).le.1.d0 .and. x < 1.d0)then
          if(dabs(cmr2).le.1.d0)then
	    X2 = dacos(cmr2)
          smr22 = 1.d0-cmr2*cmr2
          smr2 = dsqrt(smr22)
		sinout2 = smr2
          clmr2(1) = cmr2
          t3 = r3 * (1-u2*u2) / (2*del*smr2)
          t4 = v * t3
		actout2 = -t4
          t5 = t3 * 2.d0 * (x - ip * xtir1)
          ac(im,0,ip,2) = ac(im,0,ip,2) - t4
          ac(im,1,ip,2) = ac(im,1,ip,2) - t4 * cmr2
          ac(im,0,ip,3) = ac(im,0,ip,3) + t5
          ac(im,1,ip,3) = ac(im,1,ip,3) + t5 * cmr2
            do i = 2, maxcou
            clmr2(i) = 2.d0 * cmr2 * clmr2(i-1) - clmr2(i-2)
            ac(im,i,ip,2) = ac(im,i,ip,2) - t4 * clmr2(i)
            ac(im,i,ip,3) = ac(im,i,ip,3) + t5 * clmr2(i)
            end do
          end if

        end if	! (cxr = TRUE)

c	VXloop output
          if(ip.eq.1)then
			write(66,"(7G16.8)"), allkpa(im), actout1, actout2, sinout1, sinout2, X1, X2
		end if

      end do	! polarization (ip = -1, +1)
c     ``````
      end do	! k// loop (im = 1 : 2*nmoant-1)
C     >>>>>>


c     Loops over 'top right' nodes of surrounding elements:
c     Implements analytical integration of products of basis functions times linear interpolation of remaining coefficients
c     The loops over j and i sum contributions of surrounding elements times dielectric coefficients at the current mesh node.


	

	aux1 = 0.0d0
	aux2 = 0.0d0
	aux3 = 0.0d0


cE	goto 8888

        do j = max(0, ifiav(irv)-jev+1), min(1, ilaav(irv)-jev)
        ipov = 1 - j
        ji1 = jev + j
        jm1 = ji1 - 1
        do i = max(0, ifiax(irx)-iex+1), min(1, ilaax(irx)-iex)
        ipox = 1 - i
        i1 = iex + i
        im1 = i1 - 1
      



        bn(1,1) = F0i(im1, jm1, 1, isp)	! F0
        bn(2,1) = F0i(im1, jm1, 2, isp)	! dF0/dv
        bn(1,2) = F0i(im1, jm1, 3, isp)	! dF0/dx
        bn(2,2) = F0i(im1, jm1, 4, isp)	! d2F0/dvdx
        bn(1,3) = F0i(i1, jm1, 1, isp)
        bn(2,3) = F0i(i1, jm1, 2, isp)
        bn(1,4) = F0i(i1, jm1, 3, isp)
        bn(2,4) = F0i(i1, jm1, 4, isp)
        bn(3,1) = F0i(im1, ji1, 1, isp)
        bn(4,1) = F0i(im1, ji1, 2, isp)
        bn(3,2) = F0i(im1, ji1, 3, isp)
        bn(4,2) = F0i(im1, ji1, 4, isp)
        bn(3,3) = F0i(i1, ji1, 1, isp)
        bn(4,3) = F0i(i1, ji1, 2, isp)
        bn(3,4) = F0i(i1, ji1, 3, isp)
        bn(4,4) = F0i(i1, ji1, 4, isp)

c---------------------------------------------------
cPL 13/5/05: item to check: when current node is at a region boundary, may find elements with a length different from 
c the value in current region!
c---------------------------------------------------      

c A few explanations (PL26/5/2005):

c im: index of average poloidal mode; 
c ih: corresponds to ell index of theory; 
c ip: cyclotron harmonic (-1,0,1)

c Last index of ac and reac: 1, 2 or 3, correspond to coefficient of f0, df0/dv, df0/dx respectively.

c bn array: has received coefficients of bicubic basis function representation of the distribution function on current element
c

c prods array:
c 1st index v basis function
c 2nd index v derivative index, 0 or 1
c 3rd index labels function of reduced variable ksi to integrate over v element: 0 for (1-ksi), 1 for ksi
c 4th index x basis function
c 5th index x derivative index, 0 or 1
c 6th index labels function of reduced variable ksi to integrate over x element: 0 for (1-ksi), 1 for ksi

c I confirm a factor * lv * lx is missing in these expressions / pull out of nested loops as much as possible!

c	------------------------------------------------------------
cERN	25/05/05: Reorganize next loops

cERNc       Loop over v basis functions
cERN        do ib = 1, 4
cERNc         Loop over x basis functions
cERN          do jb = 1, 4
cERN            do ip = -1, 1
cERN            do ih = 0, maxcou
cERN            do im = 1, 2*nmoant-1
cERN            Coeff. of f0, df0/dv, df0/dx:
cERN            tact(im, ih, ip) = tact(im, ih, ip) + bn(ib,jb) * (
cERN     ;        ac(im,ih,ip,2) * prods(ib,1,ipov,jb,0,ipox)
cERN     ;      + ac(im,ih,ip,3) * prods(ib,0,ipov,jb,1,ipox)    )
cERN            treact(im, ih, ip) = treact(im, ih, ip) + bn(ib,jb) * (
cERN     ;        reac(im,ih,ip,1) * prods(ib,0,ipov,jb,0,ipox)
cERN     ;      + reac(im,ih,ip,2) * prods(ib,1,ipov,jb,0,ipox)
cERN     ;      + reac(im,ih,ip,3) * prods(ib,0,ipov,jb,1,ipox)      )
cERN            end do
cERN            end do
cERN            end do
cERN          end do	! jb = 1:4
cERN        end do	! ib = 1:4

c       Loop over v basis functions
        do ib = 1, 4
c         Loop over x basis functions
          do jb = 1, 4
			aux1 = aux1 + bn(ib,jb) * prods(ib,0,ipov,jb,0,ipox)  * (vgaug(jv+1)-vgaug(jv)) * (xngaug(ix+1)-xngaug(ix))
			aux2 = aux2 + bn(ib,jb) * prods(ib,1,ipov,jb,0,ipox)  * (vgaug(jv+1)-vgaug(jv)) * (xngaug(ix+1)-xngaug(ix))
			aux3 = aux3 + bn(ib,jb) * prods(ib,0,ipov,jb,1,ipox)  * (vgaug(jv+1)-vgaug(jv)) * (xngaug(ix+1)-xngaug(ix))
          end do	! jb = 1:4
        end do	! ib = 1:4

        end do  ! i
        end do  ! j


cE    8888	continue
cE	aux2 = F0i(ix, jv, 2, isp)  * (vgaug(jv+1)-vgaug(jv)) * (xngaug(ix+1)-xngaug(ix))
cE	aux3 = F0i(ix, jv, 3, isp)  * (vgaug(jv+1)-vgaug(jv)) * (xngaug(ix+1)-xngaug(ix))




	tact(1:2*nmoant-1, 0:maxcou, -1:1) = tact(1:2*nmoant-1, 0:maxcou, -1:1)
     ;									+ac(1:2*nmoant-1, 0:maxcou, -1:1, 2)*aux2
c     ;                                    +ac(1:2*nmoant-1, 0:maxcou, -1:1, 3)*aux3

c	treact(1:2*nmoant-1, 0:maxcou, -1:1) = treact(1:2*nmoant-1, 0:maxcou, -1:1)
c     ;									+reac(1:2*nmoant-1, 0:maxcou, -1:1, 1)*aux1
c     ;									+reac(1:2*nmoant-1, 0:maxcou, -1:1, 2)*aux2
c     ;                                    +reac(1:2*nmoant-1, 0:maxcou, -1:1, 3)*aux3




c	--------------------------------------------------------------------------------
	if(ix<lxg .and. jv<lvg)then


	TESTF0 = TESTF0 + (F0i(ix+1, jv+1, 1, isp) + F0i(ix, jv, 1, isp) )/2 
     ;                  * (vgaug(jv+1)-vgaug(jv)) * (xngaug(ix+1)-xngaug(ix)) 
	TESTF1 = TESTF1 + (F0i(ix+1, jv+1, 1, isp) + F0i(ix, jv, 1, isp) )/2 
     ;                  * vlel(jv) * xlel(ix)
	end if



      end if



	ac(1:2*nmoant-1, 0:maxcou, -1:1, 1:3) = 0.d0


      end do	! igx = 1:igx2 (x Gauss point loop, igx2=1 for gigdr=F)
C     ~~~~~~
      end do	! igv = 1:igv2 (v Gauss point loop, igv2=1 for gigdr=F)
C     ======
      end do  ! iex = iex1:iex2 (x Element loop, nodes for gigdr=F)   
c     ------
      end do  ! jev = jev1:jev2 (v Element loop, nodes for gigdr=F)   
c     ======

      end do  ! irx = 1, nxgreg (x Region loop, nxreg=4)  
      end do  ! irv = 1, nvgreg (v Region loop, usually nxreg=2)  

cERN	Extra normalization for Tell series -----------------
	tact = tact / lvg / lxg
c	-----------------------------------------------------


	close (66)	! (v,x) grid output
	close (123)	! tang. resonance output

c	print*,'diresp: before convolution:', etime (TT) 

c     In the end: convolution with geom. coeffs...
c     Assembly of contributions: add to VMAT, same assembly
c     as Maxwellian contributions will follow in ASPLAS.

cERN  26Apr05  : Completely modified from here (same format as M12COKPCO.f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     1) Poloidal mode differences (mdiff = m_i-m_j)
c	   (only use positive mdiff values, symmetrie applied later)

	NMM = klim + 1  
	do j = 1, NMM
	   mdiff(j) = (j - 1)  
	end do


!     3) Geom. coefficients : calculated in FOUGDR2 and 
!        stored in gcdr (-npfft/2 ... +npfft/2, -klim ... +klim)
!	   (NB: Serie is temporarely truncated at ell=128 and considered REAL!!!)
c	 Pell(1:101,1:Nm)    = gcdr(0:100,mdiff)
c	 Pell(102:Nchi,1:Nm) = dcmplx(0.0d0, 0.0d0)	
	 Pell(1:129,1:NMM)    = dreal(gcdr(0:128,mdiff))
	 Pell(130:Nchi,1:NMM) = 0.0d0	


!     3b) NEW: Determine Ell_min and Ell_max for truncating the sum
!	    based on the values of the geom.coef. P_ell(ell,m1-m2). 
!         Default values are: Ell_min=0 and Ell_max=npfft/2

	    Ell_min = 0		
	    Ell_max = 128
	    do i = 1,NMM
	       call pell_limit(mdiff(i), Ell_min(i), Ell_max(i) )
	       Ell_max(i)=min(Ell_max(i),128)  ! Pell's are truncated at ell=128
		   Ell_max(i)=max(Ell_max(i),2)    ! Minimum upper value ell=2
		   if (Ell_min(i)>Ell_max(i)) Ell_min(i) = Ell_max(i)
		end do



c	( NA = 2*(modva2-modva1)+1 )

	do j = 1, NA  ! beginning of k//-loop +++++++++++++++++++++++++++++


!      6) Loop over m1-m2 values
	
	   do i = 1, NMM ! m1-m2 loop * * * * * * * * * * * * * * * * * * * * 
		     
c			Ell_min(i) = 0
c			Ell_max(i) = 9

			 ! Perform direct sum (truncated in P_ell)
		   
			 m12left(j,i)   = sum(tact(j,Ell_min(i)+0:Ell_max(i)+0,1)
     ;                              *Pell(Ell_min(i)+1:Ell_max(i)+1,i) )
     ;				   + ci * sum(treact(j,Ell_min(i)+0:Ell_max(i)+0,1)
     ;                                *Pell(Ell_min(i)+1:Ell_max(i)+1,i) )

			 m12right(j,i)  = sum(tact(j,Ell_min(i)+0:Ell_max(i)+0,-1)
     ;                              *Pell(Ell_min(i)+1:Ell_max(i)+1,i) )
     ;				   + ci * sum(treact(j,Ell_min(i)+0:Ell_max(i)+0,-1)
     ;                                *Pell(Ell_min(i)+1:Ell_max(i)+1,i) )

			 m12landau(j,i) = sum(tact(j,Ell_min(i)+0:Ell_max(i)+0,0)
     ;                              *Pell(Ell_min(i)+1:Ell_max(i)+1,i) )
     ;				   + ci * sum(treact(j,Ell_min(i)+0:Ell_max(i)+0,0)
     ;                                *Pell(Ell_min(i)+1:Ell_max(i)+1,i) )
	   end do ! (i) end of m1-m2 loop * * * * * * * * * * * * * * * * * * *

      end do ! (j) end of A-loop ++++++++++++++++++++++++++++++++++++++++++++



cERN	cccccccccccc VERY STRANGE NORMALIZATION ccccccccccc
	m12left = -ci * m12left ! / lvg / lxg * 5.44 
c	m12left =  ci * m12left / lvg / lxg * 49.00
c	m12left = -ci * m12left / lvg / lxg *	1.4467
c	ccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN		Removed 27Apr05
cERN       iplb = 0
cERN         DO MAV2 = 2*MINF(IEL), 2*MSUP(IEL)
cERN         im = mav2 - 2*MINF(IEL) + 1
cERN          IF(MAV2 .LE. MINF(IEL)+MSUP(IEL))THEN
cERN           KSTOP2 = MIN(KLIM, MAV2-2*MINF(IEL))
cERN           ELSE
cERN           KSTOP2 = MIN(KLIM, 2*MSUP(IEL)-MAV2)
cERN           END IF
cERN         IF(MOD(MAV2-KSTOP2,2).NE.0)KSTOP2 = KSTOP2 - 1
cERN         KSTOP1 = - KSTOP2
 
cERN           DO K = KSTOP1, KSTOP2, 2
cERN           IPLB = IPLB + 1
c@ wrong array bounds!
cERN            do L = k - ncresp, k + ncresp
cERN            La = iabs(L)
cERN            vmat(1,1,iplb) = vmat(1,1,iplb) + gcdr(L-k,k) * 
cERN     ;      dcmplx(act(im,La,1), react(im,La,1))
cERN            vmat(3,3,iplb) = vmat(3,3,iplb) + gcdr(L-k,k) * 
cERN     ;      dcmplx(act(im,La,-1), react(im,La,-1))
cERN            vmat(5,5,iplb) = vmat(5,5,iplb) + gcdr(L-k,k) * 
cERN     ;      dcmplx(act(im,La,0), react(im,La,0))
cERN            end do


cERN           END DO
cERN         END DO



cERN	26Apr05: Some output added

cc    ONLY for DEBUG: Writing tact series to file
c	(Should COMMENT later)
	if (WRITE_OUTPUT) then
         FILE_NAME = COKFOLDER  // "/ACTseries" // '_lefty.dat'
		GGaux  = "(i5, f10.5, 200f20.8)"
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
c				  write (13, "(200i20.0"), '0', '0', mdiffaux(1:2*NMM)			  
				  write(13,*)"ACTIVE series for all k//"
				  write(13,*)NA,100
				  write(13,*)abscis(intab)
			      write(13,*)vttab(intab,ispe)
				  write(13,*)omcdel
				  do k = 1, NA
					 write (13, GGaux),  k, allkpa(k), tact(k,0:99,1)
				  end do 
			close (13)	

         FILE_NAME = COKFOLDER  // "/TEST_F0.dat" 
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
c				  write (13, "(200i20.0"), '0', '0', mdiffaux(1:2*NMM)			  
				  write(13,*)"Integral F0 dv dx"
				  write(13,*),'TESTF0=',TESTF0				  
				  write(13,*),'TESTF1=',TESTF1
	write(13,*),'Theory'
	write(13,*),'n  =',1.0d20*dentab(intab,ispe)
	write(13,*),'vT =',vttab(intab,ispe)
	write(13,*),'INT=',1.0d20*dentab(intab,ispe)/vttab(intab,ispe)**2/3.1415
	write(13,*),'auxs'
	write(13,*),aux1
c	write(13,*),aux2
c	write(13,*),aux3

			close (13)	




	end if
			

cc     5) Writing M12 matrices to output file

	if (WRITE_OUTPUT) then

		pana4 = paname(ispe)		! Non-Maxwellian species

c         5.1) Auxiliary vector for file output: 
c              mdiffaux = [0 0 +1 +1 .... +klim +klim]
          do j = 1, NMM
		   mdiffaux(2*j-1) = mdiff(j)
		   mdiffaux(2*j)   = mdiff(j)
	    end do 

c         5.2) Auxiliary vector for file output: FORMAT identifier
	    call INT_TO_STRING(2*NMM+1, charaux)
	    call INT_TO_STRING(2*NMM+2, charaux2)
	         GGs  = "(G11.5," // charaux // "G20.8)"
	         GG2s = "(" // charaux2 // "G20.8)"

c		p=+1:
          FILE_NAME = COKFOLDER  // "/NonMaxM12" // 'SPEC_lefty.dat'
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
				  write (13, GGs), '0.0', '0.0', mdiffaux(1:2*NMM)
				  do k = 1, NA
					 write (13, GG2s),  k, allkpa(k), m12left(k,1:NMM)
				  end do 
			close (13)	    

c		p=-1:
          FILE_NAME = COKFOLDER  // "/NonMaxM12" // 'SPEC_right.dat'
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
				  write (13, GGs), '0.0', '0.0', mdiffaux(1:2*NMM)
				  do k = 1, NA
					 write (13, GG2s), k, allkpa(k), m12right(k,1:NMM)
				  end do 
			close (13)		

c		p=0:	
	    FILE_NAME = COKFOLDER  // "/NonMaxM12" // "SPEC_landau.dat"
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
				  write (13, GGs), '0.0','0.0', mdiffaux(1:2*NMM)
				  do k = 1, NA
					 write (13, GG2s),  k, allkpa(k), m12landau(k,1:NMM)
				  end do 
			close (13)
	
	end if



c     8) Apply (m1-m2) symmetrie to all matrices 
	  
c	  (1) Shift calculated values to correct storage place (all to the right)
		  m12left  (1:NA, NMM:NMM+klim) = m12left  (1:NA, 1:NMM)
		  m12right (1:NA, NMM:NMM+klim) = m12right (1:NA, 1:NMM)
		  m12landau(1:NA, NMM:NMM+klim) = m12landau(1:NA, 1:NMM)

c	  (2) Fill the missing (negative) m1-m2 terms according to symmetrie (even)
		  do j = 1, NMM-1
			 m12left  (1:NA, j) = m12left  (1:NA, NMM+klim-j+1)
			 m12right (1:NA, j) = m12right (1:NA, NMM+klim-j+1)
			 m12landau(1:NA, j) = m12landau(1:NA, NMM+klim-j+1)
		  end do

      return

2222  format(200(f14.6))   
      end