      SUBROUTINE DIRESP_standard_2_reac_test(isp)

c	use dfport		! necessary for 'etime' function

      IMPLICIT NONE
      integer isp

C     Dielectric response of a species, using equilibrium of QLFP code
C     interpolated at local radial index INTAB of Cyrano by routine INTF0
C     isp is index of species among the ones requiring general diel. response
C     Absolute species index is ISPE																						

      include 'pardim.copy'
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
      include 'comusr2.f'
      include 'comphy.copy'
      include 'cokpco.copy'
	include 'complp.copy'
	include 'compow.copy'		 

      logical pass, cxr, toini

      save toini
      
      integer i, j, im, ix, jv, ip, iex, jev, nset0, isig, igv, igx
     ;, igv2, igx2, jev1, jev2, iex1, iex2
     ;, ispe, ii, irv, jvn, irx, ixn, icount, ic, ipov, ipox, ih
     ;, ji1, jm1, i1, im1, ib, jb
     ;, iplb, kstop1, kstop2, L, La
     ;, isrchfgt


      double precision 
     ;  v, v2, vi, xn, xna, x
      double precision 
     ;  u0, u0i, u1, u2, u1_i, u2_i, ru1, ru2, iu1, iu2, xtir1, khi
      double precision 
     ;  lamb, del2, del, te0, t1, t2, t3, t4, t5
     ;, cmr2, cosX_i, cmr2_i, smr1, smr12, smr2, smr22, cmr0
     ;, smr0, smr02, csc22
     ;, clmr0(0:128), clmr1(0:128), clmr2(0:128)
     ;, er0, r00, r1, r2, r3, s0, s1, s1i, s2, s2i, omxnai, xmax
     ;, hp, hp2, hp3, hp4, hp5, hp6, omca, xdb, xdbi, sxdb, sxdbi
c     ;, petita(5*maxnex), petitb(5*maxnex)
c     ;, nze, N1, L0, L1, jze, J1, alp2, alt2
c     ;, LL(-1:maxcou+3,5) ! JL(-1:maxcou+3,5), NL(5*maxnex,-1:maxcou+3)
c     ;, tj1(maxcou+2), tj2(maxcou+2), tj3(maxcou+2), tj4(maxcou+2)
c     ;, cmmv(5*maxnex)
c     ;, kel, s2k
c     ;, maomi
c     ;, ek(5*maxnex), ee(5*maxnex), cmm, w0, w1, w2
c     ;, oo, o1, o2, o3, o4, o5, fis
c     ;, bn(4,4)
c     ;, ellfc, ellec, ellpic

c      complex*16
c     ;  JLc(-1:maxcou+1), LLc(-1:maxcou+1), cu1, cu2, cdellpic, cdelpc
c     ;, ccosX, csmr12, cosX_cx, cmr2_cx
     
cERN	Some Variables added ccccccccccccccccccccccccccccccccccccccccccccc
	character(100) :: FILE_NAME
	character(4) :: panam4
	character(10) :: polname(2)
	integer :: OpenStat, NMM, ind_k, ind_mdiff
	integer :: mdiffaux(4*maxcou+2), mdiffvector(2*maxcou+1)
	character(100) :: GGs, GG2s, GGaux  
	character(2) :: charaux, charaux2
	real*8 :: cosXo, Acoef(2*maxpom-1), aux1, aux2, aux3, omcdel, TESTF0, TESTF1
	real*8 :: actout1, actout2, sinout1, sinout2, X1, X2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer, parameter :: Nchi = 21	! Number of chi-harmonics (ell = 0... Nchi-1) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	character(4) pana4

      real*8 :: DELV(nabsv), DELX(nabsx), zz1, zz2, tt1, tt2



c	For CHI=constant mesh (ACTIVE)

	integer :: N_CHI, ind_x, indx
	real*8 :: wcbar, wcmin, wcmax, kpar, kpar2, eps, eps2, aux2reac, xnaux, xsep
	real*8 :: Apaper, Bpaper, A, B, C, a1, b1, CHI_min, CHI_max, DELCHI, NEWFAC
      real*8 :: cosX, u(202), LAMBDAX , chi(202), Aneg, Bneg, Apapern, Bpapern
c	real*8 :: Tell_act(302,0:Nchi,-1:1,1:4), xsep
c	complex*16 :: sos
	


      real*8 :: xvec(501), sigma, kaux(501), cx1, cx2, AAFAC, aux_u1, aux_u2
	
c	real*8 :: DELK, DELE, ellf, elle, ellpi
	

c	For reactive chi-grid
      integer :: N_tab,ind_1, ind_2, ss, Lmax
	real*8 :: dx, xvec_chi(1001)
	real*8 :: Lamax, aux, rechi, XX, YY, bb, cc, dd, cosX1max, cosX1max_im, 
     ;          X1max, X1max_im 	
	complex*16 :: auxim

c	real*8 :: aachi, abchi, keffchi, Nellchi(0:Nchi), ak, bk, Lsep, usep_1, usep_2

c	real*8 :: dxdchi

	real*8 :: xnfine(800), delaux  ! , auxfine(20), rtest(1:2*nmoant-1, 0:Nchi)
	integer :: Nfine

c	For special x-grid
	real*8 :: xstar, xnstar, wid, xtang, xntang, xtang2, xntang2
	integer :: i_xstar, i_xtang, i_xtang2, Nxx
	real*8 :: DELXF(500)
	integer :: Nsepp, Nsept, Ncplx, Ntan1, Ntan2, Nroot(lxg+500)

	real*8 :: treact2(1:2*nmoant-1,0:Nchi,-1:1), xngaugx(lxg+500), 
     ;          ICroots(lxg+500,3)

	logical :: tang_res, tang_res2
	character(4) :: IDPREV

      integer isrchfge
	external isrchfge

      integer find_nearest
	external find_nearest

      real*8 find_xstar
	external find_xstar

cNEW	Anisotropic case
	logical my_anisotropy
	real*8 :: dfdv, dfdx, dfdv_prev, dfdx_prev, xprev

cNEW  For exporting RF coefficients
      logical WRITE_COEFS
      character(4) :: auxstr
	real*8 :: dchidx, Tell_aux(702,0:Nchi), auxaux, sigmaprev, dummy

cNEW  For exporting RF coefs with 'velocity files' (optimized)
      character(3) :: auxstr2


c	Cauchy Principal value
	logical CAUCHY_PRINC, www
	real*8  testreac(501, klim+1), testact(501, klim+1), testaux(501,klim+1), 
     ;                               treact2test(33, klim+1), treact2aux(33,klim+1),
     ;                               NOVO, xntir, tact2test(33, klim+1)  
	real*8  xaux, auxtang, xteste, fac, delt
	complex*16   auxm12(100,100)
	integer ind2, indaux


	real*8  raiz, vperp2, vpar2, vpar, btt, sintt, deltwa, gammadd, auxres(600)
	integer ipo, mdiff
	real*8 f_reso_std
	external f_reso_std

	real*8 f_ic_std
	external f_ic_std


      integer Nraux(500)
	real*8 ICaux(500,2)


      logical ANAL_REAC


c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     external ellfc, ellec, ellpic, cdellpic, cdelpc
c	external DELK, DELE, ellf, elle, ellpi	

      data toini/.true./


c	my_anisotropy = .true.
      WRITE_COEFS = .false.
	only_active=.true.
	ANAL_REAC=.true.
	


c     CHECK if this is a BATCH surface for exporting the RF coefs.
c     (The BATCH surface indexes are in intabqlfp -> see read_QLFPheader.f)
cx      if(any(intab.eq.intaqlfp)) then
cx	   WRITE_COEFS = .true.
cx	end if


C	NEW variables for refining x-grid
c	Total number must be less than 75
	Nsepp = 12		! Number of pts near sep. - passing side
	Nsept = 8       ! Number of pts near sep. - trapped side
	Ncplx = 15		! Double root resonance
	Ntan1 = 20		! Tangent resonance u=A
	Ntan2 = 20		! Tangent resonance u=B

cccccccccccccccccccccc
      fac = glofac * twopi ** 2
      if(.not.cyl)fac = fac * ra
cccccccccccccccccccccc

c	Number of k// values (also given by NA = 2*nmoant-1)	
	NA = 2*(modva2-modva1)+1

c	Number of poloidal couplings (modva2-modva1+1 for all coupled)
	NMM = klim + 1


c     Index of k// for data output for DEBUGGING
      ind_k = 1
	ind_mdiff = 1 ! (mdiff=ind-1)

c	Reset some variables
	m12left   = czero
	m12right  = czero
	m12landau = czero

	treact2test(1:NA,1:NMM) = 0.d0
	tact2test(1:NA,1:NMM) = 0.d0
	treact2aux(1:NA,1:NMM) = 0.d0

c     Non-Maxwellian Species index
      ispe = ispgdr(isp)
      
c     Normalized x at separatrix (co-passing):
      xnsep = - 2.d0 * delb(intab) / bmax(intab) 

c	Radius dependent quantities
      omca = qom(ispe) * B0
	omcdel = qom(ispe) * delb(intab)
      xtir1 = omca / omegag
c      oo = 1.d0 / dabs(xtir1)
      te0 =  PI / (qom(ispe)*hachi(intab)*delb(intab))
      s0 = bbar(intab) / delb(intab)
      xmax = B0 / bmin(intab)
c      maomi = bmax(intab) / bmin(intab)
      s1 = B0 / delb(intab)
      s1i = 1.d0 / s1
      s2 = bmin(intab) / delb(intab)
      s2i = delb(intab) / bmin(intab)

cERN	added
	wcbar = qom(ispe) * bbar(intab)
	wcmin = (wcbar-omcdel)
	wcmax = (wcbar+omcdel)
      xsep  = (1-dabs(xnsep))*xmax
	Lamax = xmax/xtir1
c	Lsep  = xsep/xtir1
	xntir = -1.d0+xtir1/xmax 
		
c	Normalization factor
	NOVO = PI * (EEL*ZCH(ispe))**2
     ; * dPsidr_n(intab) * abscis(intab)
     ; / (2.d0*MH*AMASS(ispe) * B0)
c	Divide by ASPLAS factor
c      NOVO = NOVO/fac



c	New xn-grid refined at separatrix (xngaugx)
	xngaugx(1:lxg) = xngaug(1:lxg)
	Nxx = lxg
cccccccccccccccccccccccccccccccccccccccccccccc
c	  sigma = +1
c	  call insert_points2(Nxx, xngaugx(1:Nxx), xnsep, 'R', 2, 4, 2.d-4, xngaugx(1:Nxx+4))
c	  Nxx = Nxx + 4
c	  call insert_points2(Nxx, xngaugx(1:Nxx), xnsep, 'L', 2, 8, 2.d-4, xngaugx(1:Nxx+8))	
c	  Nxx = Nxx + 8
c	  sigma = -1
c	  call insert_points2(Nxx, xngaugx(1:Nxx), -xnsep, 'R', 2, 8, 2.d-4, xngaugx(1:Nxx+8))
c	  Nxx = Nxx + 8
c	  call insert_points2(Nxx, xngaugx(1:Nxx), -xnsep, 'L', 2, 4, 2.d-4, xngaugx(1:Nxx+4))
c	  Nxx = Nxx + 4
cccccccccccccccccccccccccccccccccccccccccccccc

c	UPDATE for REACTIVE:
cccccccccccccccccccccccccccccccccccccccccccccc
c	  sigma = +1
	  call insert_points2(Nxx, xngaugx(1:Nxx), xnsep, 'R', 2, 10, 2.d-4, xngaugx(1:Nxx+10))
	  Nxx = Nxx + 10
	  call insert_points2(Nxx, xngaugx(1:Nxx), xnsep, 'L', 2, 20, 2.d-4, xngaugx(1:Nxx+20))	
	  Nxx = Nxx + 20
c	  call insert_points2(Nxx, xngaugx(1:Nxx), xnsep, 'L', 2, 10, 2.d-4, xngaugx(1:Nxx+10))	
c	  Nxx = Nxx + 10
c	  sigma = -1
	  call insert_points2(Nxx, xngaugx(1:Nxx), -xnsep, 'R', 2, 20, 2.d-4, xngaugx(1:Nxx+20))
	  Nxx = Nxx + 20
c	  call insert_points2(Nxx, xngaugx(1:Nxx), -xnsep, 'R', 2, 10, 2.d-4, xngaugx(1:Nxx+10))
c	  Nxx = Nxx + 10
	  call insert_points2(Nxx, xngaugx(1:Nxx), -xnsep, 'L', 2, 10, 2.d-4, xngaugx(1:Nxx+10))
	  Nxx = Nxx + 10
cccccccccccccccccccccccccccccccccccccccccccccc





cERN  NEW: Temporary output for exporting the RF coefs. for BATCH
c     (Only at the mag. surfaces given by intaqlfp(nrad))

      if(WRITE_COEFS)then
            call int_to_string4(nint(1000*abscis(intab)),auxstr)
c      call int_to_string2(intab,auxstr)
          
            FILE_NAME = paplda  // "/RFtemp_" // auxstr // ".dat" 
		open (UNIT = 9090, FILE = FILE_NAME, STATUS = "REPLACE",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		      end if	  
c		write(9090,*),"RF coefs. for BATCH"
            write(9090,*),abscis(intab)
	      write(9090,*)
	      write(9090,*), NA, lvg
	      write(9090,*)

	end if

cERN	OUTPUT for radial quantities ccccccccccccccccccccccccccccc
	if (WRITE_OUTPUT) then
         FILE_NAME = trim(STDFOLDER)  // "/radial_data.dat" 
		open (UNIT = 99, FILE = FILE_NAME, STATUS = "REPLACE",
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                    stop
		        end if
		  	  write(99,*),"Some radius dependent quantities (DIRESP)"
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
			  write(99,*) ' NOVO=   ', NOVO 
                    write(99,*) ' cos(Pitch)= ', eqt(intab,1,14)	  ! at theta=0
			  write(99,*) ' Ip=    ' ,  eqta1d(intab,3)*abscis(intab)**2	
c			  write(99,*) ' H=    ' ,  1/hachi(intab)
			  write(99,*) ' Nk//=  ', NA
			  do j = 1,NA
			     write(99,*), allkpa(j)
			  end do
			  write(99,*) ' fac=  ', fac
		close(99)
	end if



ccccccccccccccccc Open files for OUTPUT ccccccccccccccccccccccccc
c     (Mainly for debugging, should remove later on!)

     
  

c	(v,x) loop OUTPUT (reactive fine grid) ----------------------------------
         FILE_NAME = trim(STDFOLDER)  // "/alpha2fine.dat" 
		open (UNIT = 76, FILE = FILE_NAME, STATUS = "REPLACE",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
		      end if	  
            write(76,*),"REACTIVE: xn,ur,cosX,ELLPI,Jo,J1... for Non-Maxwellians (DIRESP)"        
  


   
               
       
               
c     Total sums OUTPUT ---------------------------------------

c         FILE_NAME = trim(STDFOLDER)  // "/total_sums.dat" 
c		open (UNIT = 85, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;            IOSTAT = OpenStat, ACTION = "WRITE")
c            write(85,*),"Total sums over (x,v) and (chi,v) "        


         FILE_NAME = trim(STDFOLDER)  // "/testact.dat" 
		open (UNIT = 1550, FILE = FILE_NAME, STATUS = "REPLACE",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
            write(1550,*),"test direct integration"  
            write(1550,*),ind_mdiff-1
		   
         FILE_NAME = trim(STDFOLDER)  // "/testreac.dat" 
		open (UNIT = 1551, FILE = FILE_NAME, STATUS = "REPLACE",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
            write(1551,*),"test direct integration"  

         FILE_NAME = trim(STDFOLDER)  // "/Cauchy_details.dat" 
		open (UNIT = 1552, FILE = FILE_NAME, STATUS = "REPLACE",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
            write(1552,*),"Cauchydetails"  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         FILE_NAME = trim(STDFOLDER)  // "/Cauchy_details_anal.dat" 
		open (UNIT = 1559, FILE = FILE_NAME, STATUS = "REPLACE",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
            write(1559,*),"Cauchydetails (analytical integ.)"  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         FILE_NAME = trim(STDFOLDER)  // "/Cauchy_debugs.dat" 
		open (UNIT = 1553, FILE = FILE_NAME, STATUS = "REPLACE",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
            write(1553,*),"Cauchydebugs"  



c         FILE_NAME = trim(STDFOLDER)  // "/IC_roots.dat" 
c		open (UNIT = 1554, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;            IOSTAT = OpenStat, ACTION = "WRITE")
c            write(1554,*),"IC resonance roots"  


         FILE_NAME = trim(STDFOLDER)  // "/IC_roots_onesigma.dat" 
		open (UNIT = 1555, FILE = FILE_NAME, STATUS = "REPLACE",
     ;            IOSTAT = OpenStat, ACTION = "WRITE")
            write(1555,*),"IC resonance roots"  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c	External (m1+m2)/2-loop 

      do im = 1, 2*nmoant-1 ! k//-loop > > > > > > > > > > > > > > > > > > > > > > > > 
C     >>>>>>>>>>>>>>>>>>>>>

	mave = modva1 + 0.5d0*dfloat(im-1)  ! (m1+m2)/2	
	print*,im, mave
    
c     Is this used in the new routines?????
cccccccccccccc k// = k//(theta)  (COMMAG.copy)
	   if(circ)then
c           Magnetic angle is constant on flux surfaces
		  si = eqt(intab,1,15)
		  co = eqt(intab,1,14)
	      kptab(1:npfft+1) = mave * si * rhoinv + kphi * co * r0orta(intab,1:npfft+1)    
	   else
c           k// = mave*sinthn/ntn + n*costh/R:  (kphi=n/Ro))
            kptab(1:npfft+1) = mave * eqt(intab,1:npfft+1,12) / eqt(intab,1:npfft+1,7)
     ;           + kphi * eqt(intab,1:npfft+1,14) * r0orta(intab,1:npfft+1)
	   end if	
cccccccccccc



	!!!!!! k// at theta=0 !!!!!!
	kpar = kptab(1)
	!!!!!! k// at theta=pi !!!!!!
	kpar2 = kptab(npfft/2+1)

      if(im .eq. ind_k)then
	   write(76,*),mave, kpar, kpar2
c	   write(1554,*),mave, kpar, kpar2
	   write(1555,*),mave, kpar, kpar2
      end if

	write(9090,*), mave

c      write(85,*),'k// at theta=0=',kpar

c     For debugging orbit integrals	
      if(im.eq.ind_k)then
	   www=.true.
	else
	   www=.false.
	end if




cERN	Simple v-loop 

	do jv = 1, lvg   ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c	================
	jev = jv
      v = vgaug(jv)
      v2 = v * v

c	Element lengths	
	if( jv < lvg)then
	    DELV(jv) = vgaug(jv+1)-vgaug(jv)
	else
	    DELV(jv) = vgaug(jv)-vgaug(jv-1)
	end if

      if(im .eq. ind_k)then
	   write(76,*),v
c	   write(1554,*),v
	   write(1555,*),v
      end if



c	Reset fine grid to 'original' grid refined near separatrix
	xnfine(1:Nxx) = xngaugx(1:Nxx)
	Nfine = Nxx



ccccccccccccc New : veloc. files for RF coef ccccccccccccccccccccc

      if(WRITE_COEFS .and. im.eq.1)then
            call int_to_string2(jv,auxstr2)
         
            FILE_NAME = trim(trim(STDFOLDER)) // "/v_" // auxstr2 // ".dat" 

		open (UNIT = 2000+jv, FILE = FILE_NAME, STATUS = "REPLACE",
     ;	      IOSTAT = OpenStat, ACTION = "WRITE")
	            if (OpenStat > 0) then
		          print *, 'Error writing file: ', FILE_NAME
	                stop
			end if	  
            write(2000+jv,*)v 


      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc










cccccccccccccccccc TO BE REVISED cccccccccccccccccccccccccc

cccccccccccccccccccccccccccc	kpar = dabs(kpar)
	kpar = dabs(kpar)    ! LFS (Bmin)
	kpar2 = dabs(kpar2)  ! HFS (Bmax)

c	Auxiliary constants for fixed (kpar,v) 


        eps = (kpar*v)/omegag;
	  eps2 = (kpar2*v)/omegag;
        Apaper = (omegag-wcmin)/(kpar*v);
        Bpaper = (omegag-wcmax)/(kpar2*v);



c	Determine if there is a tangent resonance? (for u1, A<1)
c	UPDATE for u2 LATER !!!!!!!!!!!!
	tang_res=.false.
	if(dabs(Apaper)<1.0d0)tang_res=.true.

	tang_res2=.false.
	if(dabs(Bpaper)<1.0d0)tang_res2=.true.



cccccccccccccccccccccccccccc TO BE REVISED ccccccccccccccccccccccccccccc






1010	continue



cccccccccccccccccccccccccccccc ACTIVE + REACTIVE (fine x - mesh) ccccccccccccccccccccccccccccccccc

ccccccc Refined x-mesh, insert points near singularities


c     1) Delta=0 singularity -------------------------------------

c	   xstar = 2.d0*xtir1/(eps*eps)  * (1.d0-dsqrt(1-eps*eps)) 
c        (from k//=const. coordinates)      

c        Routine based on cubic equation for circular cross-section
         xstar = find_xstar(qom(ispe), v, mave, xtir1, xmax)

	   xnstar = -1.d0+xstar/xmax  
	   i_xstar = find_nearest(lxg, xnfine, xnstar)
c	   Width = Delta_x original
	   if(i_xstar.eq.1)i_xstar=2
	   wid = dabs(xnfine(i_xstar+1)-xnfine(i_xstar-1))

c	    sigma = +1
c	    call insert_points2(Nfine, xnfine(1:Nfine), xnstar, 'R', 2, Ncplx, 1.d-5, xnfine(1:Nfine+Ncplx))
c		Nfine = Nfine + Ncplx
	    call insert_points2(Nfine, xnfine(1:Nfine), xnstar, 'R', 2, 10, 1.d-5, xnfine(1:Nfine+10))
		Nfine = Nfine + 10
cccccccccccc FOR REACTIVE
	    call insert_points2(Nfine, xnfine(1:Nfine), xnstar, 'R', 2, 10, 1.d-5, xnfine(1:Nfine+10))
		Nfine = Nfine + 10
cccccccccccccccccccccccc

	    call insert_points2(Nfine, xnfine(1:Nfine), xnstar, 'L', 2, 10, 1.d-4, xnfine(1:Nfine+10))
		Nfine = Nfine + 10
	    call insert_points2(Nfine, xnfine(1:Nfine), xnstar, 'L', 2, 15, 1.d-5, xnfine(1:Nfine+15))
		Nfine = Nfine + 15

c	    sigma = -1
c	    call insert_points2(Nfine, xnfine(1:Nfine), -xnstar, 'L', 2, Ncplx, 1.d-5, xnfine(1:Nfine+Ncplx))
c		Nfine = Nfine + Ncplx
	    call insert_points2(Nfine, xnfine(1:Nfine), -xnstar, 'L', 2, 10, 1.d-5, xnfine(1:Nfine+10))
		Nfine = Nfine + 10
cccccccccccc FOR REACTIVE
	    call insert_points2(Nfine, xnfine(1:Nfine), -xnstar, 'L', 2, 10, 1.d-5, xnfine(1:Nfine+10))
		Nfine = Nfine + 10
cccccccccccccccccccccccccc

	    call insert_points2(Nfine, xnfine(1:Nfine), -xnstar, 'R', 2, 10, 1.d-4, xnfine(1:Nfine+10))
		Nfine = Nfine + 10
	    call insert_points2(Nfine, xnfine(1:Nfine), -xnstar, 'R', 2, 5, 1.d-5, xnfine(1:Nfine+5))
		Nfine = Nfine + 5	

c	goto 1441

c     3) Tangent resonance u = A
	if(tang_res)then
		
	   xtang = xtir1 * (1-Apaper**2)/(1-eps*Apaper)
	   xntang = -1.d0 + xtang/xmax
c	   i_xtang = isrchfge(Nfine, xnfine(1:Nfine), 1, xntang)

	   ! LEFT (ACTIVE: non-singular)
	   i_xtang = isrchfge(Nfine, xnfine(1:Nfine), 1, xntang)
	   if(i_xtang > 3)then
ccccccccccccccccccccccccccccccccccccccccc
c	      call insert_points2(Nfine, xnfine(1:Nfine), xntang, 'L', 2, 5, 1.d-5, xnfine(1:Nfine+5))
c	      Nfine = Nfine + 5
ccccccccccccccccccccccccccccccccccccccccccc
c	UPDATE for REACTIVE
ccccccccccccccccccccccccccccccccccccccccc
	      call insert_points2(Nfine, xnfine(1:Nfine), xntang, 'L', 3, 20, 1.d-5, xnfine(1:Nfine+20))
	      Nfine = Nfine + 20
	      call insert_points2(Nfine, xnfine(1:Nfine), xntang, 'L', 2, 10, 1.d-5, xnfine(1:Nfine+10))
	      Nfine = Nfine + 10
	      call insert_points2(Nfine, xnfine(1:Nfine), xntang, 'L', 2, 20, 1.d-5, xnfine(1:Nfine+10))
	      Nfine = Nfine + 20
ccccccccccccccccccccccccccccccccccccccccccc
 	   end if
         ! RIGHT (ACTIVE: singular)

	   call insert_points2(Nfine, xnfine(1:Nfine), xntang, 'R', 3, 20, 1.d-4, xnfine(1:Nfine+20))
	   Nfine = Nfine + 20
	   call insert_points2(Nfine, xnfine(1:Nfine), xntang, 'R', 2, 10, 1.d-4, xnfine(1:Nfine+10))
	   Nfine = Nfine + 10



	else
	   xtang  = 0.d0
	   xntang = 0.d0
	end if

c     4) Tangent resonance u = B
	if(tang_res2)then
		
	   xtang2 = xtir1 * (1-Bpaper**2)/(1-eps2*Bpaper)
	   xntang2 = 1.d0 - xtang2/xmax
c	   i_xtang2 = isrchfge(Nfine, xnfine(1:Nfine), 1, xntang2)

	   ! LEFT (ACTIVE: singular)
	   call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'L', 3, 20, 1.d-4, xnfine(1:Nfine+20))
	   Nfine = Nfine + 20
	   call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'L', 2, 10, 1.d-4, xnfine(1:Nfine+10))
	   Nfine = Nfine + 10


	   ! RIGHT (ACTIVE: non-singular)
	   i_xtang2 = isrchfge(Nfine, xnfine(1:Nfine), 1, xntang2)
	   if(i_xtang2<Nfine-0)then
cccccccccccccccccccccccccccccc
c	      call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'R', 2, 5, 1.d-5, xnfine(1:Nfine+5))
c	       Nfine = Nfine + 5
ccccccccccccccccccccccccccccccccccccccccccc
c	UPDATE for REACTIVE
ccccccccccccccccccccccccccccccccccccccccc
	   if(i_xtang2<Nfine-1)then
	      call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'R', 3, 20, 1.d-5, xnfine(1:Nfine+20))
	       Nfine = Nfine + 20
	      call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'R', 2, 10, 1.d-5, xnfine(1:Nfine+10))
	       Nfine = Nfine + 10
	      call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'R', 2, 10, 1.d-5, xnfine(1:Nfine+10))
	       Nfine = Nfine + 10
	   else
	      call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'R', 2, 20, 1.d-5, xnfine(1:Nfine+20))
	       Nfine = Nfine + 20
		call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'R', 2, 10, 1.d-5, xnfine(1:Nfine+10))
	       Nfine = Nfine + 10
		call insert_points2(Nfine, xnfine(1:Nfine), xntang2, 'R', 2, 20, 1.d-5, xnfine(1:Nfine+10))
	       Nfine = Nfine + 20
	   end if
cccccccccccccccccccccccccccccccccccc
	   end if

	else
	   xtang2  = 0.d0
	   xntang2 = 0.d0
	end if


cNEW	Insert points near xntir (begining of double root) ----------------------------
c	   sigma=+1	  
	   call insert_points2(Nfine, xnfine(1:Nfine), xntir, 'R', 2, 10, 1.d-5, xnfine(1:Nfine+10))
	   Nfine = Nfine + 10
	   call insert_points2(Nfine, xnfine(1:Nfine), xntir, 'L', 2, 5, 1.d-5, xnfine(1:Nfine+5))
	   Nfine = Nfine + 5

c	   sigma=-1	  
	   call insert_points2(Nfine, xnfine(1:Nfine), -xntir, 'L', 2, 5, 1.d-5, xnfine(1:Nfine+5))
	   Nfine = Nfine + 5
	   call insert_points2(Nfine, xnfine(1:Nfine), -xntir, 'R', 2, 10, 1.d-5, xnfine(1:Nfine+10))
	   Nfine = Nfine + 10
c	--------------------------------------------------------------------------









	if(im .eq. ind_k)then
		write(76,*) xntang, xntang2, Nfine
	    write(76,*) xtir1, xmax
		write(76,*) 'x*', xnstar, -1.d0+xteste/xmax
	    write(76,*) xstar, xteste
	end if



	if(im.eq.ind_k)then
		write(1550,"(I10)"), Nfine
		write(1551,"(I10)"), Nfine
c		write(1554,"(I10)"), Nfine
		write(1555,"(I10)"), Nfine
	    write(1553,*)
	    write(1553,*), mave, v
          write(1553,"(I10)"), Nfine
	endif

      if(im.eq.ind_k .and. v.gt.0.9d6 .and. v<1.1d6)then
         write(1552,*), kpar, v
         write(1552,"(I10)"), Nfine
         write(1559,*), kpar, v
         write(1559,"(I10)"), Nfine
c	   write(1553,*), kpar, v
c         write(1553,"(I10)"), Nfine
	endif

      if(WRITE_COEFS)then
         write(2000+jv,*)im,mave  
         write(2000+jv,*)Nfine
      end if


c	RF coefs (temporary)
      if(WRITE_COEFS)then
      write(9090,*)
      write(9090,*), v
      write(9090,*)
	write(9090,"(I10)"), Nfine
	end if



cccccccccccccccccccccccccccccc ACTIVE + REACTIVE (fine x - mesh) ccccccccccccccccccccccccccccccccc


c     Auxiliary variable for TRAPEZE integration
c	IDPREV = 'XXXX'

	testact(1:Nfine,1:NMM)  = 0.d0
	testreac(1:Nfine,1:NMM) = 0.d0
	testaux(1:Nfine,1:NMM) = 0.d0

	do ix = 1, Nfine  ! Fine x - loop ------------------------------------------------
c	~~~~~~~~~~~~~~~~


c	NEW: To avoid problems near xtir in standard coordinates ---------
	  delt = 2.d-5
cc	  sigma>0 (xn<0)
	  if( dabs(xnfine(ix) - xntir)  < delt)then
            xnfine(ix) = xntir + dsign(1.d0,xnfine(ix)-xntir) * delt
	      if(xnfine(ix)>xnstar-1.d-5)xnfine(ix)=xnstar+1.d-5     !! Look at insert points @ xstar
	  end if

cc	  sigma>0 (xn<0)
	  if( dabs(xnfine(ix) - (-xntir))  < delt)then
            xnfine(ix) = -xntir + dsign(1.d0,xnfine(ix)-(-xntir)) * delt
	      if(xnfine(ix)<-xnstar+1.d-5)xnfine(ix)=-xnstar-1.d-5   !! Look at insert points @ xstar
	  end if
c     ------------------------------------------------------------------

c	NEW: To avoid problems near xtang --------------------------------
	  delt = 2.d-5
cc	  sigma>0 (xn<0)
	  if( xntang.ne.0.d0 .and. dabs(xnfine(ix) - xntang)  < delt)then
            xnfine(ix) = xntang + dsign(1.d0,xnfine(ix)-xntang) * delt
c	      if(xnfine(ix)>xnstar-1.d-5)xnfine(ix)=xnstar+1.d-5     !! Look at insert points @ xstar
	  end if

cc	  sigma<0 (xn>0)
	  if( xntang2.ne.0.d0 .and. dabs(xnfine(ix) - (-xntang2))  < delt)then
            xnfine(ix) = -xntang2 + dsign(1.d0,xnfine(ix)-(-xntang2)) * delt
c	      if(xnfine(ix)<-xnstar+1.d-5)xnfine(ix)=-xnstar-1.d-5   !! Look at insert points @ xstar
	  end if
c     -----------------------------------------------------------------



ccccccccccccccccccccccccccccccccccccccccc

        xn = xnfine(ix)
	  sigma = -dsign(1.d0,xn)    
	  x = (1.d0 - dabs(xn) ) * xmax

        auxtang = 0.d0
	  if(sigma>0) auxtang = xtang
	  if(sigma<0) auxtang = xtang2



ccccccccccccccc TEST on differnet BRANCHES cccccccccccccccccccc


        Nth_oi=100
					
        if(x.le.xsep)then ! PASSING (1 root) ----------------------------------------------
          

c		 ===================================================================
c            ROOT:
             call find_IC_roots_onesig(qom(ispe), v, x, sigma, mave, ICaux(ix,1:2), Nraux(ix))
  		      if(im.eq.ind_k)then
                     write(1555,"(10G20.8)")v, xn, ICaux(ix,1:2), Nraux(ix)
	            end if  
                             
c	       ACTIVE:
             call oi_pass_act(x, v, sigma, mave, NMM, qom(ispe), auxtang, xtir1,
     ;                                         ICaux(ix,1), testact(ix,1:NMM), www)

         if(.not.only_active)then
c            REACTIVE:
            if(ANAL_REAC)then
		    call oi_pass_react_anal2(x, v, sigma, mave, NMM, qom(ispe), auxtang, xtir1,
     ;                                  ICaux(ix,1), testaux(ix,1:NMM), www)
	       else
		    call oi_pass_react     (x, v, sigma, mave, NMM, qom(ispe), auxtang, xtir1,
     ;                                  ICaux(ix,1), testreac(ix,1:NMM), www)
	       endif
c		 ====================================================================

                 endif

        elseif(x.gt.xsep .and. x.lt.xtir1)then  ! TRAPPED (1 root)  -------------------------------
               
                   
c		 ====================================================================
c            ROOT:
             call find_IC_roots_onesig(qom(ispe), v, x, sigma, mave, ICaux(ix,1:2), Nraux(ix))
		      if(im.eq.ind_k)then
                     write(1555,"(10G20.8)")v, xn, ICaux(ix,1:2), Nraux(ix)

	            end if  

c		 ACTIVE:
		 call oi_trap_act(x, v, sigma, mave, NMM, qom(ispe), auxtang, xtir1,
     ;                                         ICaux(ix,1), testact(ix,1:NMM), www)

          if(.not. only_active)then
c		 REACTIVE:	
             if(ANAL_REAC)then
		    call oi_trap_react_anal2(x, v, sigma, mave, NMM, qom(ispe), auxtang, xtir1,
     ;                                  ICaux(ix,1), testaux(ix,1:NMM), www)             
             else	    
		    call oi_trap_react     (x, v, sigma, mave, NMM, qom(ispe), auxtang, xtir1,
     ;                                  ICaux(ix,1), testreac(ix,1:NMM), www)
	       endif
c		 =====================================================================

              endif


        elseif(x.gt.xtir1 .and. x.lt.xstar)then  ! TRAPPED (double root)  -------------------------------
               

c		 =====================================================================
c            ROOTS:
	       call find_IC_roots_onesig(qom(ispe), v, x, sigma, mave, ICaux(ix,1:2), Nraux(ix))
		      if(im.eq.ind_k)then
	               write(1555,"(10G20.8)")v, xn, ICaux(ix,1:2), Nraux(ix)
	            end if   

c		 ACTIVE:
		 call oi_trap_act_double(x, v, sigma, mave, NMM, qom(ispe), xtang, xtang2, xtir1,
     ;                                         ICaux(ix,2), ICaux(ix,1), testact(ix,1:NMM), www)
                 if(.not. only_active)then
c		 REACTIVE:
             if(ANAL_REAC)then
		 call oi_trap_react_double_anal2(x, v, sigma, mave, NMM, qom(ispe), xtang, xtang2, xtir1,
     ;                                         ICaux(ix,2), ICaux(ix,1), testaux(ix,1:NMM), www)
             else
		 call oi_trap_react_double(x, v, sigma, mave, NMM, qom(ispe), xtang, xtang2, xtir1,
     ;                                         ICaux(ix,2), ICaux(ix,1), testreac(ix,1:NMM), www)
             endif

           end if
c		 =====================================================================


        elseif(x.ge.xstar .and. x.le.xmax)then  ! TRAPPED (complex roots)  -------------------------------

c		 =====================================================================

c		 ACTIVE:
             testact(ix,1:NMM) = 0.d0   ! Always non-resonant
                 if(.not.only_active)then
c		 REACTIVE:
             if(ANAL_REAC)then
		 call oi_trap_react_cplx_anal2(x, v, sigma, mave, NMM, qom(ispe), 
     ;                                    testaux(ix,1:NMM), www)
             else
		 call oi_trap_react_cplx     (x, v, sigma, mave, NMM, qom(ispe), 
     ;                                    testreac(ix,1:NMM), www)
             endif
                endif

			if(im.eq.ind_k)then 
c			      write(1554,"(10G20.8)")v, xn, 0.0d0, 0.0d0, 0
	                  write(1555,"(10G20.8)")v, xn, 0.d0, 0.d0, 0
					  co = eqt(intab,1,14)
                        write(1553,"(400G16.8)") x, xn, 0, 
     ;                  dacos(r0/abscis(intab)*(x/co-1.d0)), 0, 0
		    end if    
c			=====================================================================

        else    ! Should be nothing left

                print*, 'Problems: x is out of range'
	          print*,'(x,xmax)',x,xmax
				print*,'STOP at DIRESP_standard_2_test line 925.'
	          stop

        end if  ! (x < xsep) -------------------------------------------------------


cccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccc
        if(ANAL_REAC)then
           testreac(ix,1:NMM) = testaux(ix,1:NMM)
        endif
ccccccccccccccccccccccccccccccc


	if(im.eq.ind_k)then


c             print*, 2.d0*NOVO*v*testact(ix,ind_mdiff)

		write(1550,"(10G18.8)")v,xn,2.d0*NOVO*v*testact(ix,ind_mdiff)

		write(1551,"(10G18.8)")v,xn, 2.d0*NOVO*v*testreac(ix,ind_mdiff),
     ;                                   2.d0*NOVO*v*testaux(ix,ind_mdiff)

c		write(1554,"(10G20.8)")v, xn, ICroots(ix,1:3), Nroot(ix)
c		write(1555,"(10G20.8)")v, xn, ICaux(ix,1:2), Nraux(ix)

	endif

	testreac(ix,1:NMM) =  -2.d0*NOVO*v*testreac(ix,1:NMM)
	testaux(ix,1:NMM)  =  -2.d0*NOVO*v*testaux(ix,1:NMM)
	testact(ix,1:NMM)  =   2.d0*NOVO*v*testact(ix,1:NMM)


      if(WRITE_COEFS)then
         write(2000+jv,"(400g16.8)")xn, testact(ix,1:NMM)  
      end if


		if(WRITE_COEFS)then      

	    write(9090,"(G16.8)") xn
          
	    do j = 1, NMM
	       write(9090,"(I10, G20.10)"), j-1, testact(ix,j)
          end do

	    end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     ``````



c	Add (x,v) contribution to the REACTIVE (summing) variable 'treact'



cNEW	Anisotropic case -------------------------------------------------

c	   1) Interpolate F0i to the given x-fine point 

	   call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 2, isp), 1, xnfine(ix), dfdv)
	   call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 3, isp), 1, xnfine(ix), dfdx)
              dfdx = dfdx / xmax  ! QLFP gives df/dxn, not df/dx


         if(ix .gt. 1)then

	   call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 2, isp), 1, xnfine(ix-1), dfdv_prev)
	   call my_interpolation(lxg, xngaug(1:lxg), F0i(1:lxg, jv, 3, isp), 1, xnfine(ix-1), dfdx_prev)
              xprev = (1.d0 - dabs(xnfine(ix-1)) ) * xmax
	        dfdx_prev = dfdx_prev / xmax ! QLFP gives df/dxn, not df/dx

         end if



c     -> TRAPEZE integration

      if(.not.only_active)then
ccccccccccccccccccccccccccccccccccccc
        if(ix .gt. 1)then
           treact2test(im,1:NMM) = treact2test(im,1:NMM) 
     ;     + 0.5d0*(   testreac(ix,1:NMM) *   (v*v*v*dfdv-2.d0*v*v*(x-xtir1)*dfdx) 
     ;               + testreac(ix-1,1:NMM) * (v*v*v*dfdv_prev-2.d0*v*v*(xprev-xtir1)*dfdx_prev)) 
     ;     * (xnfine(ix)-xnfine(ix-1))*xmax * DELV(jv) 
	  end if
ccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccc
        if(ix .gt. 1)then
           treact2aux(im,1:NMM) = treact2aux(im,1:NMM) 
     ;     + 0.5d0*(   testaux(ix,1:NMM) *   (v*v*v*dfdv-2.d0*v*v*(x-xtir1)*dfdx) 
     ;               + testaux(ix-1,1:NMM) * (v*v*v*dfdv_prev-2.d0*v*v*(xprev-xtir1)*dfdx_prev)) 
     ;     * (xnfine(ix)-xnfine(ix-1))*xmax * DELV(jv) 
	  end if
ccccccccccccccccccccccccccccccc
        endif


ccccccccccccccccccccccccccccccccccccc
        if(ix .gt. 1)then
           tact2test(im,1:NMM) = tact2test(im,1:NMM) 
     ;     + 0.5d0*(   testact(ix,1:NMM) *   (v*v*v*dfdv-2.d0*v*v*(x-xtir1)*dfdx) 
     ;               + testact(ix-1,1:NMM) * (v*v*v*dfdv_prev-2.d0*v*v*(xprev-xtir1)*dfdx_prev)) 
     ;     * (xnfine(ix)-xnfine(ix-1))*xmax * DELV(jv) 
	  end if
ccccccccccccccccccccccccccccccc





      end do	! (fine x-loop) -------------------------------------------------
C     ~~~~~~

	if(WRITE_COEFS)then      
	   write(9090,*) 
      end if
            
c	Some v-dependent output
      if(im.eq.ind_k)then
         write(76,*), " "
      end if
 
  
      end do	!  (simple v-loop)  = = = = = = = = = = = = = = = = = = = = = = =
C     ======


      end do	! k//-loop (im = 1 : NA) > > > > > > > > > > > > > > > > > > > > >
C     >>>>>>




      close (76)
	close(1550)
	close(1551)
	close(1552)
	close(1553)
	close(1553)

	do jv = 1, lvg
	   close(2000+jv)
	end do

	close(9090)



c     1) Poloidal mode differences (mdiff = m_i-m_j)
c	   (only use positive mdiff values, symmetrie applied later)

c	NMM = klim + 1  
c	do j = 1, NMM
c	   mdiff(j) = (j - 1)  
c	end do


        FILE_NAME = trim(STDFOLDER)  // "/testmy_lefty.dat"
			open (UNIT = 13, FILE = FILE_NAME, STATUS = "REPLACE",
c     ;			  BUFFERED='YES', BLOCKSIZE = 32*512, BUFFERCOUNT = 50, 
     ;              IOSTAT = OpenStat, ACTION = "WRITE")
	              if (OpenStat > 0) then
		              print *, 'Error writing file: ', FILE_NAME
	                  stop
				  end if
				  do k = 1, NA
					 write (13, "(I8, 200G20.8)"),  k, allkpa(k), 
     ;                                                treact2aux(k,1:NMM)/fac, tact2test(k,1:NMM)/fac
				  end do 
			close (13)


c	Read REACTIVE from Maxwellian results
	if(only_active)then  

c       Auxiliary vector for file output: FORMAT identifier
	  call INT_TO_STRING(2*NMM+1,  charaux)
	  call INT_TO_STRING(2*NMM+2, charaux2)
	  GGs  = "(G11.5," // charaux // "G20.8)"
	  GG2s = "(" // charaux2 // "G20.8)"

           FILE_NAME = trim(STDFOLDER)  // "/M12Dbea_lefty.dat"
		 open (UNIT = 13, FILE = FILE_NAME, STATUS = "UNKNOWN",
     ;           IOSTAT = OpenStat, ACTION = "READ")
	           if (OpenStat > 0) then
		           print *, 'Error reading file: ', FILE_NAME
	               stop
			   end if
			   read(13, GGs), dummy, dummy, mdiffaux(1:2*NMM)
			   do j = 1, NA
				  read (13, GG2s), dummy, dummy, auxm12(j,1:NMM)
			   end do 
		 close (13)

		
c	 m12left(1:NA,1:NMM) = dimag(m12left(1:NA,1:NMM)) + dreal(auxm12(1:NA,1:NMM))
        treact2test(1:NA,1:NMM) = dreal(auxm12(1:NA,1:NMM)) * fac

	end if



      m12left(1:NA,1:NMM) = dcmplx(tact2test(1:NA,1:NMM), treact2test(1:NA,1:NMM))


cERN	cccccccccccc RADIAL and CYRANO (-i) NORMALIZATION ccccccccccc
	m12left   = -ci *  m12left  / fac   ! 1/fac for consistency with VMAT !!!!!
c	m12right  = -ci *  m12right
c	m12landau = -ci *  m12landau
c	ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


			

cc     5) Writing M12 matrices to output file

	if (WRITE_OUTPUT) then

		pana4 = paname(ispe)		! Non-Maxwellian species

c         5.1) Auxiliary vector for file output: 
c              mdiffaux = [0 0 +1 +1 .... +klim +klim]
          do j = 1, NMM
		   mdiffaux(2*j-1) = (j - 1) 
		   mdiffaux(2*j)   = (j - 1) 
	    end do 

c         5.2) Auxiliary vector for file output: FORMAT identifier
	    call INT_TO_STRING(2*NMM+1, charaux)
	    call INT_TO_STRING(2*NMM+2, charaux2)
	         GGs  = "(G11.5," // charaux // "G20.8)"
	         GG2s = "(" // charaux2 // "G20.8)"

c		p=+1:
          FILE_NAME = trim(STDFOLDER)  // "/NonMaxM12" // 'SPEC_lefty.dat'
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
          FILE_NAME = trim(STDFOLDER)  // "/NonMaxM12" // 'SPEC_right.dat'
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
	    FILE_NAME = trim(STDFOLDER)  // "/NonMaxM12" // "SPEC_landau.dat"
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

c	  (2) Fill the missing (negative) m1-m2 terms according to symmetry (even)
		  do j = 1, NMM-1
			 m12left  (1:NA, j) = m12left  (1:NA, NMM+klim-j+1)
			 m12right (1:NA, j) = m12right (1:NA, NMM+klim-j+1)
			 m12landau(1:NA, j) = m12landau(1:NA, NMM+klim-j+1)
		  end do

      return

2222  format(200(f14.6))   

      end
