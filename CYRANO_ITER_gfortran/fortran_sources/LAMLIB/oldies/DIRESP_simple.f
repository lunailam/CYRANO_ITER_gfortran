      SUBROUTINE DIRESP_simple(isp)

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
     ;, JL(-1:maxcou+3,5), LL(-1:maxcou+3,5) ! NL(5*maxnex,-1:maxcou+3),
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
	real*8 :: DELV(nabsv), DELX(nabsx), zz1, zz2, tt1, tt2, sigma
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      equivalence (JLc, JL), (LLc, LL)

      external ellfc, ellec, ellpic, cdellpic

      data toini/.true./

cERN	Number of k// values (also given by NA = 2*nmoant-1)	
	NA = 2*(modva2-modva1)+1

cERN	 Reset Auxiliary variables

	m12left   = czero
	m12right  = czero
	m12landau = czero



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
                    write(99,*) ' cos(Pitch)= ', eqt(intab,1,14)	  ! at theta=0
				  write(99,*) ' Ip=    ' ,  eqta1d(intab,3)*abscis(intab)**2
				  write(99,*) ' Nk//=   ', NA
					do j = 1,NA
						write(99,*), allkpa(j)
					end do

			close(99)

	end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



 

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



cERN	Simple v-loop added -----------------------------
	do jv = 1, lvg
c	================
	jev = jv
	
	
		
      v = vgaug(jv)
      v2 = v * v
      vi = 1.d0 / v
      er0 = r00 * v2
      r1 = er0 * v2 * PI

	if( jv < lvg)then
		DELV(jv) = vgaug(jv+1)-vgaug(jv)
	else
		DELV(jv) = vgaug(jv)-vgaug(jv-1)
	end if


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


cERN	Simple x-loop added -----------------------------
	do ix = 1, lxg
c	~~~~~~~~~~~~~~~~
	iex = ix

        xn = xngaug(ix)
	  sigma = - dsign(1.d0, xn)	

cERN	Avoid xn=xnsep
	if(xn.eq.xnsep .and. irx.eq.1) xn = xnsep-0.0002d0
	if(xn.eq.xnsep .and. irx.eq.2) xn = xnsep+0.0002d0
	if(xn.eq.-xnsep .and. irx.eq.3) xn = -xnsep-0.0002d0
	if(xn.eq.-xnsep .and. irx.eq.4) xn = -xnsep+0.0002d0

cERN	Avoid xn=0
	if(xn.eq.0.d0 .and. irx.eq.2) xn = -0.002d0
	if(xn.eq.0.d0 .and. irx.eq.3) xn = +0.002d0


	if( ix < lxg)then
		DELX(ix) = xngaug(ix+1)-xngaug(ix)
	else
		DELX(ix) = xngaug(ix)-xngaug(ix-1)
	end if


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





c      if(.not.gigdr)then
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
cERN          if(dabs(cmr1) .le. 1.d0)then
		if(dabs(cmr1) .le. 1.d0 .and. dsign(1.d0,u1) .eq. sigma)then
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
c	    ac(im,0,ip,2) =  - t4
c        ac(im,1,ip,2) =  - t4 * cmr1
c         ac(im,0,ip,3) =  + t5
c          ac(im,1,ip,3) =  + t5 * cmr1
            do i = 2, maxcou
            clmr1(i) = 2.d0 * cmr1 * clmr1(i-1) - clmr1(i-2)
            ac(im,i,ip,2) = ac(im,i,ip,2) - t4 * clmr1(i)
            ac(im,i,ip,3) = ac(im,i,ip,3) + t5 * clmr1(i)
c            ac(im,i,ip,2) = - t4 * clmr1(i)
c            ac(im,i,ip,3) = + t5 * clmr1(i)
            end do
          end if

		actout2 = 0.0
		sinout2 = 0.0
	      X2=0.0
cERN          if(dabs(cmr2).le.1.d0)then
          if(dabs(cmr2).le.1.d0 .and. dsign(1.d0,u2) .eq. sigma)then
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
c          ac(im,0,ip,2) =  - t4
c         ac(im,1,ip,2) =  - t4 * cmr2
c          ac(im,0,ip,3) =  + t5
c          ac(im,1,ip,3) =  + t5 * cmr2


            do i = 2, maxcou
            clmr2(i) = 2.d0 * cmr2 * clmr2(i-1) - clmr2(i-2)
            ac(im,i,ip,2) = ac(im,i,ip,2) - t4 * clmr2(i)
            ac(im,i,ip,3) = ac(im,i,ip,3) + t5 * clmr2(i)
c            ac(im,i,ip,2) =  - t4 * clmr2(i)
c            ac(im,i,ip,3) =  + t5 * clmr2(i)
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



	   aux2 = F0i(ix, jv, 2, isp)  * DELX(ix) * DELV(jv)
	   aux3 = F0i(ix, jv, 3, isp)  * DELX(ix) * DELV(jv)

cE	Test 'real' Maxwellian

c	  zz1 = 1.d20*dentab(intab,2) / (sqrtpi * vttab(intab,2)) ** 3
c		  zz2 = - 2.d0 / vttab(intab,2)**2
c		  tt1 = zz1 * dexp(- (v / vttab(intab,2)) ** 2)
c		  tt2 = zz2 * v * tt1

c		aux2 = dmin1(tt2,-1.d-40)  * DELX(ix) * DELV(jv)
c	    aux3 = 0.0d0  * DELX(ix) * DELV(jv)

cE	----> Gives the same results as with numerical dFdv, since I 
c		inserted the Analytical expression in INTF0.f



	tact(1:2*nmoant-1, 0:maxcou, -1:1) = tact(1:2*nmoant-1, 0:maxcou, -1:1)
     ;									  +ac(1:2*nmoant-1, 0:maxcou, -1:1, 2)*aux2
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



c      end if	! gigdr

cERN	RESET variable ac after (x,v) contribution has been added to tact
c	%%%%%%%%%%%%% (ThIS TAKES A LOT OF TIME %%%%%%%%%%%%%%)

	ac(1:2*nmoant-1, 0:maxcou, -1:1, 1:3) = 0.0d0



      end do	! igx = 1:igx2 (simple x-loop)
C     ~~~~~~
      end do	! igv = 1:igv2 (simple v-loop)
C     ======


cERN	Extra normalization for Tell series -----------------
c	tact = tact ! / lvg / lxg
c	-----------------------------------------------------


	close (66)	! (v,x) grid output
	close (123)	! tang. resonance output


c	-------------------------------------------------------------
cc    ONLY for DEBUG: Writing tact series to file
c	(Should COMMENT later)

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
c	----------------------------------------------------------------











c	cccccccccccccc CONVOLUTION ccccccccccccccccccccccc

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



cERN	cccccccccccc CYRANO NORMALIZATION ccccccccccc
	m12left = -2*ci * m12left ! / lvg / lxg * 5.44 
c	m12left =  ci * m12left / lvg / lxg * 49.00
c	m12left = -ci * m12left / lvg / lxg *	1.4467
c	ccccccccccccccccccccccccccccccccccccccccccccccccccc



cERN	26Apr05: Some output added

	if (WRITE_OUTPUT) then

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
c				write(13,*),'auxs'
c				write(13,*),aux1
c				write(13,*),aux2
c				write(13,*),aux3

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