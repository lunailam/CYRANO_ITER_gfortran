      SUBROUTINE DIRESP_simple_chi(isp)

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
     ;, NL(5*maxnex,-1:maxcou+3), JL(-1:maxcou+3,5), LL(-1:maxcou+3,5)
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
c	real*8 :: AAleft(2*maxpom-1), x_tang, 
c     ;          xn_tang_pos(2*maxpom-1), xn_tang_neg(2*maxpom-1)
c	integer :: indx_k(2*maxpom-1)



	real*8 :: DELV(nabsv), DELX(101), zz1, zz2, tt1, tt2





c	For chi=constant mesh

	integer :: N_CHI
	real*8 :: wcbar, wcmin, wcmax, kpar, eps
	real*8 :: A, B, C, a1, b1, CHI_1, CHI_2, DELCHI, NEWFAC
      real*8 :: cosX(101), u(101), LAMBDAC(101) , chi(201)
	real*8 :: Tell_dv(2*maxpom-1,0:maxcou,-1:1,1:4), xsep
	complex*16 :: sos
c	real*8 :: xtab(1001), Etab(1001), Ktab(1001), kp, kt, keff(1001), 
c     ;          N0tab(1001), N1tab(1001), DELK, DELE
c	real*8 :: kp, kt, keff(1001), DELK, DELE
c	external DELK, DELE
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

cERN	added
	wcbar = qom(ispe) * bbar(intab)
	wcmin = (wcbar-omcdel)
	wcmax = (wcbar+omcdel)
      xsep = (1-dabs(xnsep))*xmax

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

	end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c     Build tables of elliptic integrals K(k) and E(k) for fixed x-mesh x=[0...Xmax]
c     Using IMSL elliptic integrals
c     Output through COMMONS : xtab(x), Ktab(x) and Etab(x) 
c     OBS: Using 1000 x-points to avoid further interpolation inside main loops

      call DIRESP_tables





cERN	External k//-loop added

      do im = 1, 2*nmoant-1 ! k//-loop >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     >>>>>>>>>>>>>>>>>>>>>

	kpar = allkpa(im)	


cERN	Simple v-loop added -----------------------------
	do jv = 1, lvg
c	================
	jev = jv
	
	
      v = vgaug(jv)
      v2 = v * v
      vi = 1.d0 / v
      er0 = r00 * v2
      r1 = er0 * v2 * PI

c	Element lengths	
	if( jv < lvg)then
		DELV(jv) = vgaug(jv+1)-vgaug(jv)
	else
		DELV(jv) = vgaug(jv)-vgaug(jv-1)
	end if

c	Auxiliary 'constants'
        a1 = omcdel/(kpar*v);
        b1 = (omegag-wcbar)/(kpar*v);
        eps = (kpar*v)/omegag;
        A = (omegag-wcmin)/(kpar*v);
        B = (omegag-wcmax)/(kpar*v);
        C = (omegag-wcbar)/omcdel;	! should go out of loop ....

	NEWFAC = v*v*v * r00 * pi * B0 /(4*delb(intab))


c	Build CHI grid according to CHI_max and CHI_min values  
 
	  call cplxacos(-1/a1-C,sos)
        CHI_2 = real(sos);		! max ????
	  call cplxacos(+1/a1-C,sos)
        CHI_1 = real(sos);		! min ????
        DELCHI = dabs(CHI_1-CHI_2)/50;	! CHI element length (constant)
        do k=1,51
           chi(k) = CHI_1 + (k-1)*DELCHI; ! Constant CHI mesh   
        end do
        N_CHI = 51;



c	Now at fixed (kpar,v)
	



cERN	Simple CHI-loop added -----------------------------

	do ix = 1, N_CHI
c	~~~~~~~~~~~~~~~~
	  iex = ix


            cosX(ix) = dcos(chi(ix));
            u(ix)    = a1 * cosX(ix) + b1;
            LAMBDAC(ix) = ((u(ix))**2 - 1) / (eps*u(ix)-1);

            x = LAMBDAC(ix) * xtir1
c            xn(ix) =  1-x(ix)/Xmax;


c       petita(ix) = (1-x*bmin(intab)/b0)
c       petitb(ix) = (1-x*bmax(intab)/b0)




c	old stuff -------------------
c        x = (1.d0 - dabs(xn)) * xmax
c        khi1 = x / xtir1
cC       x*delta/B0:
c        xdb = x / s1
c        xdbi = 1.d0 / xdb
c        sxdb = dsqrt(xdb)
c        sxdbi = 1.d0 / sxdb
c	-----------------------------




cERN	V - CHI output
	write(66,*),'     '
	write(66,"(5G16.8)"), v, chi(ix), F0i(ix, jv, 1, isp), F0i(ix, jv, 2, isp)
     ;                    , x



      r3 = er0 * PI * s1 * 0.25d0 * omegag / dabs(omca)


cc	TO MODIFY LATER -------------------------------------------------
cc     p=0:
c      u0i = (allkpa(im) * v) / omegag
c      u0 = 1.d0 / u0i
c      cmr0 = s0 - xdbi * (1.d0 - u0 * u0)
c      smr02 = 1.d0 - cmr0 * cmr0
c      csc22 = 2.d0 / (1.d0 - cmr0)
cC       p=0 only absorbs if isig*k//>0, |u0i|>1, and x between...
cc        if(isig .eq. dint(dsign(1.d0,u0i)) .and. dabs(u0i).gt.1.d0)then
c        if(u0*fis.gt.petitb(ix) .and. u0*fis.lt.petita(ix))then
cc	if(abs(u0)<1 .or. abs(petita(ix))>1)then
cc	print *, petita(ix), petitb(ix), u0
cc      end if
cc	  if(u0.gt.petitb(ix) .and. u0.lt.petita(ix))then
c
c        smr0 = dsqrt(smr02)
c        clmr0(1) = cmr0
cc      r2 = r1 * xdbi * hp
c	print *, 'Active p=0'
cc	Active p=0:
c        t3 = er0 * PI * xdbi * dabs(u0)**3 / smr0
c        ac(im,0,0,2) = - v * t3
c        ac(im,0,0,3) =  2.d0 * x * t3
c        ac(im,1,0,2) = cmr0 * ac(im,0,0,2)
c        ac(im,1,0,3) = cmr0 * ac(im,0,0,3)
c          do i = 2, maxcou
c          clmr0(i) = 2.d0 * cmr0 * clmr0(i-1) - clmr0(i-2)
c          ac(im,i,0,2) = ac(im,0,0,2) * clmr0(i)
c          ac(im,i,0,3) = ac(im,0,0,3) * clmr0(i)
c          end do
c        end if
c	------------------------------------------------------------------




c     p=±1: 2 roots to resonance equations.
      do ip = -1, 1, 2
c     ````````````````
c      khi = dfloat(ip) * khi1
c      lamb = 0.5 * u0i * khi
c      del2 = lamb * lamb + 1.d0 - khi
c      cxr = del2 .lt. 0.d0		! complex root



	    if(dabs(u(ix)) .le. 1.d0)then	! add contribution to Tell

             Tell_dv(im,0,ip,2) =  - NEWFAC * LAMBDAC(ix) * (A-B)/2
             Tell_dv(im,1,ip,2) =  - NEWFAC * LAMBDAC(ix) * (A-B)/2 * cosX(ix)

             ! ell series by cos(ell.X) recursion
             cmr1 = cosX(ix)
             clmr1(1) = cmr1
             do i = 2, maxcou
c                clmr1(i) = 2.d0 * cmr1 * clmr1(i-1) - clmr1(i-2)
ccccccccccccccccccccccccccccccccccccccccccccccccc
				clmr1(i) = dcos(i*chi(ix))
cccccccccccccccccccccccccccccccccccccccccccccccc
                Tell_dv(im,i,ip,2) = - NEWFAC * LAMBDAC(ix) * (A-B)/2 * clmr1(i)    
             end do

			
			actout1 = - NEWFAC * LAMBDAC(ix) * (A-B)/2

	    end if


c	VXloop output
c          if(ip.eq.1)then
c			write(66,"(7G16.8)"), allkpa(im), actout1
c		end if


      end do	! polarization (ip = -1, +1)
c     ``````
c      end do	! k// loop (im = 1 : 2*nmoant-1)
Cc     >>>>>>


	   aux2 = F0i(1, jv, 2, isp)    * DELV(jv) * DELCHI
	   aux3 = F0i(1, jv, 3, isp)    * DELV(jv) * DELCHI

c	Trapeze integration

c	if( jv < lvg)then
c	   aux2 = (F0i(1, jv+1, 2, isp) + F0i(1, jv, 2, isp))/2.0     * DELV(jv) * DELCHI
c	   aux3 = (F0i(1, jv+1, 3, isp) + F0i(1, jv, 2, isp))/2.0      * DELV(jv) * DELCHI
c	else
c	   aux2 = (F0i(1, jv, 2, isp) + F0i(1, jv-1, 2, isp))/2.0     * DELV(jv) * DELCHI
c	   aux3 = (F0i(1, jv, 3, isp) + F0i(1, jv-1, 2, isp))/2.0      * DELV(jv) * DELCHI
c	end if



cE	Test 'real' Maxwellian

c	  zz1 = 1.d20*dentab(intab,2) / (sqrtpi * vttab(intab,2)) ** 3
c		  zz2 = - 2.d0 / vttab(intab,2)**2
c		  tt1 = zz1 * dexp(- (v / vttab(intab,2)) ** 2)
c		  tt2 = zz2 * v * tt1

c		aux2 = dmin1(tt2,-1.d-40)  * DELCHI * DELV(jv)
c	    aux3 = 0.0d0  * DELCHI * DELV(jv)

cE	----> Gives the same results as with numerical dFdv, since I 
c		inserted the Analytical expression in INTF0.f



	tact(im, 0:maxcou, -1:1) = tact(im, 0:maxcou, -1:1)
     ;						 + 2 * Tell_dv(im, 0:maxcou, -1:1, 2) * aux2
c     ;                                    +ac(im, 0:maxcou, -1:1, 3)*aux3

c	treact(1:2*nmoant-1, 0:maxcou, -1:1) = treact(1:2*nmoant-1, 0:maxcou, -1:1)
c     ;									+reac(1:2*nmoant-1, 0:maxcou, -1:1, 1)*aux1
c     ;									+reac(1:2*nmoant-1, 0:maxcou, -1:1, 2)*aux2
c     ;                                    +reac(1:2*nmoant-1, 0:maxcou, -1:1, 3)*aux3



      end do	! (simple CHI-loop)
C     ~~~~~~


c     Reactive part (use x-mesh instead of CHI-mesh)

	do ix = 1, lxg
c	~~~~~~~~~~~~~~~~
	  iex = ix
        xn = xngaug(ix)
        x = (1.d0 - dabs(xn)) * xmax       

        khi1 = x/xtir1

        do ip = -1, 1, 2
c     ````````````````
        khi = dfloat(ip) * khi1
        lamb = 0.5 * u0i * khi
        del2 = lamb * lamb + 1.d0 - khi
        cxr = del2 .lt. 0.d0		! complex root

        if(cxr)then

cERN		complex roots! To be done!!!!!!!!!!!!!!!!!!!!!!!1		

        else

        del = dsqrt(del2)
        u1 = lamb + del
        u2 = lamb - del
        t1 = dfloat(ip) / (qom(ispe) * delb(intab))
        t2 = omegag - dfloat(ip) * qom(ispe) * bbar(intab)
        cmr1 = t1 * (allkpa(im) * v * u1 - t2)
        cmr2 = t1 * (allkpa(im) * v * u2 - t2)



        end if	! (cxr = TRUE)

c	VXloop output
          if(ip.eq.1)then
			write(66,"(7G16.8)"), allkpa(im), actout1, actout2, sinout1, sinout2, X1, X2
		end if


      end do	! polarization (ip = -1, +1)
c     ``````


      end do	! (simple x-loop) ----------------- REACTIVE part
C     ~~~~~~





   
c	print*, a1, b1, u(1)	  
      end do	!  (simple v-loop)
C     ======

      end do	! k// loop (im = 1 : 2*nmoant-1)
C     >>>>>>


	close (66)	! (v,x) grid output


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
		     
			Ell_min(i) = 0
			Ell_max(i) = 9

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
	m12left = -ci * m12left ! / lvg / lxg * 5.44 
c	m12left =  ci * m12left / lvg / lxg * 49.00
c	m12left = -ci * m12left / lvg / lxg *	1.4467
c	ccccccccccccccccccccccccccccccccccccccccccccccccccc


			

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