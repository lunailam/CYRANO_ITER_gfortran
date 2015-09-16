      subroutine asplas(onespe, icurspe)

c      SUBROUTINE asplas

cJAC	use dfport		! necessary for 'etime' function

      implicit none

	logical onespe

	integer icurspe

c     Generation of finite element problem
c     Assembly of plasma terms for the current element 
      
c     Includes loops over average pol. mode index.
C     Hot model with poloidal mode expansion
C
C     07/09/89: worked out further optimization in small matrix products
c
c     8/4/2004: added arguments onespe, icurspe allowing assembly for a single species
c     Purpose: allow evaluation of power balance per species, integrated over each element.
c              This to be achieved in queftc2, called in genera3.
c     Inputs:
c	onespe = .false.: normal use, assemble contribution of all species
c              .true:   only assemble species # icurspe
c     icurspe: dummy if onespe is false; index of species to assemble if onespe is true.
c
      include 'pardim.copy'

c     NB: continuation line is not allowed after cf77 include '...' !
c     Careful here:
c     File COMUSR2 must be updated every time file COMUSR changes!
c      COMPLEX*16 V(6,6,2*MAXPOM**2), TOTOP(12,12,BLLEN)
      include 'comusr2.f'
c     ;, vmat, totop

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
	include 'cokpco.copy'

      character*3 eletyp

      logical inl, inr, jnl, jnr, trafo

      integer 
     ;  in(6), jn(6), index(12), jndex(12)
     ;, ipop1(3), jpop1(3), ipop2(3), jpop2(3), npopi, npopj
     ;, i, j, idum
     ;, nmodc, indext, m1, m2, ni, nj, ip, jp
     ;, ideca, jdeca, id, jd, nset0, kstop1, kstop2, iplb, iplbs, ipath
     ;, iplbmu, nwkbl, i1, i2, j1, j2, jndext
     ;, idllo, icolo, ibulo, icpblo, ielet
     ;, nsum, i1st
     ;, isp, nfft, ig2, nno, intast, ino, ipro, ic1, l1, ir1, it1, ib1
     ;, jb1, id1, ic2, l2, ir2, it2, ib2, jb2, id2
      
ccccccccccccccccccc VARIABLES (cokpco) cccccccccccccccccccccccccccccc

	integer :: k_index, m_index, OpenStat
	logical :: resp
	real :: T1(2)	            !   Variables for evaluating
	real :: tempo, tempo0       !   the elapsed CPU time
	character(4) :: rrr  
	character(30) :: shell_comm
	character(60) :: shell_comm2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision 
     ;  cabs1, rle
     ;, second
     ;, rm

      complex*16 zdum

      complex*16 axl(3), axc(3), tem3, fac
     ;, v2(0:npfft,28)
     ;, to(12,12), toto(12*12), tos(12,12), totos(12*12)
     ;, big(12,12), bigs(12,12), big1(12*12), big1s(12*12)
     ;, cri

      save nset0, m1, m2, jndex, index, ielet, eletyp, idllo, icolo, ibulo
     ;, icpblo

      equivalence (toto,to), (totos,tos), (big1,big), (big1s,bigs)

      external second, nsum, i1st

      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))

c	tempo0 = etime (T1)
 
cERN	Reset VMAT elements to zero (necessary for cokpco = T)
	vmat(1:ndof,1:ndof,1:2*maxpom**2) = czero

c	-----------------------------------------------------

c     FAC: common factor to all matrix elements
c     remember field spectrum has dimensions:
c       - V/m in toroidal geometry
c       - V   in a cylinder
      if(cokpco)then
	fac = rnorm
	else
      fac = glofac * twopi ** 2
      if(.not.cyl)fac = fac * ra
	end if

      eletyp = styp(isubr)
      ielet = istyp(isubr)
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo

      trafo = .not. polsym
      rle = fl(iel)
      nmodc = max0(nmode(iel-1), nmode(iel))

C     Index list for assembly:
      indext = nmodc * icpblo
      do i = 1, icolo
	   index(i) = i
	   jndex(i) = i
	   index(i+icpblo) = i + indext
	   jndex(i+icpblo) = i + indext
      end do 

cERN  ccccccccccccccc Opening files for cokpco OUTPUT cccccccccccccccccc

cJAC      if(cokpco) then

c	  if(WRITE_OUTPUT) then
c	     open (UNIT = 5555, FILE = '../../M12run/index_radii.dat', 
c     ;           STATUS = "REPLACE", IOSTAT = OpenStat, ACTION = "WRITE")
c	       if (OpenStat > 0) then
c		      print *, 'Error creating file: ../../M12run/index_radii.dat'
c	          stop
c	       end if	   	
c	  end if

c	  if(WRITE_REPORT) then
c           open (UNIT = 7777, FILE = '../../cokpco_report.txt', 
c     ;           STATUS = "REPLACE", ACTION = "WRITE")
c	     write(7777,*)'====================================================='
c		 write(7777,*)'Calculation of plasma response M12 for k// = constant'
c		 write(7777,*)'====================================================='
c		 write(7777,*)
c	     write(7777,*)'Number of poloidal mode diff.(m1-m2) =',2*klim+1
c		 write(7777,*)'Number of A (k//) per element =',2*(modva2-modva1)+1
c	     if (READ_TABLES .eq. .FALSE.) then
c		    write(7777,*)'Method : Calculating M12 values explicitly'
c	     else
c		    write(7777,*)'Method : Interpolating M12 values from tables'
c	     end if
c	   end if

c	   if (WRITE_SCREEN) then
c	     print *,'====================================================='
c	     print *,'Calculation of plasma response M12 for k// = constant'
c	     print *,'====================================================='
c	     print *
c	     print *,'Number of poloidal mode diff.(m1-m2) =',2*klim+1
c	     print *,'Number of A (or k//) per element =',2*(modva2-modva1)+1
c	     if (READ_TABLES .eq. .FALSE.) then
c		    print *, 'Method : Calculating M12 values explicitly'
c	     else
c		    print *, 'Method : Interpolating M12 values from table files'
c	     end if
c	   end if

c	end if ! cokpco = true

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     First element at magnetic axis or crown inner boudary
c     !nmode(0) = nmode(1) assumed!

c	========================================================== 
      if(iel.eq.1)then ! BEGINNING of first element computations
c     ==========================================================

      intab = 1 + (iel-1) * (ngauss+1) + ireg

      if(circ .and. eletyp.eq.'M23') call loctri(iel,.true.)
      if(klim.ge.nmodc-1)then
         nset0 = nmodc**2
      else
         nset0 = (klim + 1) * (2 * nmodc - klim) - nmodc
      end if
      write(nofile,*)'nbmax, bllen, 2*nset0:', nbmax, bllen, 2*nset0
      
c     Terms will be summed into totop:
      call zset(2*144*nset0, czero, totop, 1)
 
      m1 = minf(iel)
      m2 = msup(iel)
      kstop1 = 0
      kstop2 = 0

      do 717 ig = 1, ngauss  ! First element Gauss points loop
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        intab = ig + (iel-1) * (ngauss+1) + ireg
        y = abscno(intab)
        iplb = 0
        iplbs = nset0
        ipath = 1

c       1) Case of constant k// coordinates -------------------------------------

        if(cokpco) then  
	
           print*, '(ASPLAS) Not ready for cokpco in this version!'
cJAC	    tempo0 = etime (T1)
cJAC	    print *, '		 --> (ASPLAS.f) Begin FOUGDR + M12COKPCO:', tempo0

cJACc         1.1) Create folder for data output at radial point: '../../M12run/rXXXX/'    
cJAC		if (WRITE_OUTPUT) then
cJAC	       call int_to_string4 (nint(1000*abscis(intab)), rrr)	
cJAC	       COKFOLDER = '../../M12run/r' // rrr
cJAC	       shell_comm = 'mkdir ..\..\M12run\r' // rrr
cJAC	       resp = SYSTEM(shell_comm)
cJAC	       if (resp .eq. -1) then
cJAC	          print *, 'Error creating folder: ', COKFOLDER
cJAC		      stop
cJAC	       end if
cJAC		   write(5555,'(i8, g15.6, A20)'), intab, abscis(intab), COKFOLDER
cJACc	       print *, intab, abscis(intab), COKFOLDER
cJAC	    end if  ! WRITE_OUTPUT = true

cJACc		1.2) Geometrical couplings for plasma response in 'constant k//' coord.:
cJACc		     This call must be reviewed when non-Maxwellian behaviour will be 
cJACc		     implemented, because Max. and non-Max. need different geom. coef.
cJACc		     NEW: If reading M12 table files we don't need this call!!
cJAC		if (READ_TABLES .eq. .FALSE.) then
cJACcERN			call fougdr(1)
cJAC			call fougdr_ern(1)	! fougdr_ern: vectorized version (checked OK!)
cJAC		else
cJAC		    if (WRITE_OUTPUT) then
cJAC		        shell_comm2 = 'copy ..\..\M12tables\r' // rrr // '\Plm.dat ' 
cJAC    ;                            // '..\..\M12run\r' // rrr 
cJAC				resp = SYSTEM(shell_comm2)
cJAC			end if
cJAC		end if
c		     The geom. coef. are calculated using FFT with respect to variable Chi and 
c		     are stored in the matrix gcdr(L,mdiff), where L is the Chi-harmonic index
c		     and mdiff = m1-m2 is the poloidal mode difference.
c	           -> L range     : [-npfft/2 .... -1 0 +1 .... +npfft/2]
c	           -> mdiff range : [-klim .... -1 0 +1 .... +klim]	

c		1.3) Calculate M12 matrices (already summed for all species)
c              The results are stored in the common's (COKPCO.COPY):
c                m12left   (k_index, m_index) : for p = +1 polarization
c		       m12right  (k_index, m_index) : for p = -1 polarization
c                m12landau (k_index, m_index) : for p = 0  polarization
c	                   -> k_index : index for k// (or A)
c                        -> m_index : index for m1-m2 (poloidal couplings)

cJAC	    call M12COKPCO(onespe, icurspe)

c	         NB: If onespe = TRUE, M12COKPCO will return the response of a single
c		         species (identified as icurspe) in m12left, m12right and m12landau.
c		         Otherwise, it returns the standard response summed over all species.

c		1.4) Write M12 elements to VMAT matrix 
c		     To be careful because data in m12(i,j) matrices and in VMAT have different 
c              structures:
c	           -> k// (or A) index : 'i' in m12 is 'mav2-2*m1+1' in VMAT
c	           -> mdiff index : 'j' in m12 is 'k+klim+1' in VMAT

c		Loop over average pol. mode :
cJAC		do mav2 = 2*m1, 2*m2 ! mav2 = 2*(m_i+m_j) loop ++++++++++++++++
cJAC		   k_index = mav2-2*m1+1
cJAC      	   if(mav2 .le. m1+m2)then
cJAC		  	  kstop2 = min(klim, mav2-2*m1)
cJAC		   else
cJAC		      kstop2 = min(klim, 2*m2-mav2)
cJAC		   end if
cJAC		   if(mod(mav2-kstop2,2).ne.0) kstop2 = kstop2 - 1
cJAC		   kstop1 = - kstop2
cJAC		   do k = kstop1, kstop2, 2 ! k = m_i-m_j loop ------------
cJAC		      iplb = iplb + 1
cJACcERN		      k_index = mav2-2*m1+1
cJAC			  m_index = k+klim+1
cJAC	          vmat(1,1,iplb) = m12left  (k_index, m_index) ! p = +1
cJAC	          vmat(3,3,iplb) = m12right (k_index, m_index) ! p = -1
cJAC	          vmat(5,5,iplb) = m12landau(k_index, m_index) ! p = 0 
cJAC		   end do  ! end of k loop --------------------------------
cJAC          
cJAC		end do ! end of mav2 loop ++++++++++++++++++++++++++++++++++++

cJAC		iplbmu = iplb
	
cJAC		tempo = etime (T1)
cJAC		print *, '		 --> (ASPLAS.f)   End FOUGDR + M12COKPCO:', tempo
cJAC		print *, '                        Time used ', tempo-tempo0
	
	  end if  ! cokpco=true

c       End case of constant k// coordinates ------------------------------------


c	  2) Non-Maxwellians ???? -------------------------------------------------

        if(nspgdr.gt.0) then ! + + + + + + + 
          print*, '(ASPLAS)nspgrd not ready!'
cJAC		 if(onespe)then
cJACcPL@          to do: allow evaluation 1 species at a time:
cJAC              print *,'ASPLAS gdr: must implement calc. 1 species at a time'
cJAC	        stop
cJAC		 end if
cJACc		 Species with a general dielectric response:
cJACc		 Compute x, v meshes at current radius and interpolate equil. distr. f0:
cJAC		 call intf0
cJACc		 Add contribution to VMAT blocks. Currently 1 species allowed.
cJAC		 isp = 1
c		 do isp = 1, nspgdr
cJAC		 call diresp(isp)
cJACc		 end do
      
	  end if ! nspgdr > 0 + + + + + + + + 

c	  End Non-Maxwellians ???? ------------------------------------------------


c	  3) Case of standard angle variables -------------------------------------
c          Maxwellian species:

        if(.not.cokpco)then

c		 Loop over average pol. mode (mav2 = 2 maverage):
		 do mav2 = 2*m1, 2*m2 ! mav2 = 2*(m_i+m_j) loop + + + + +
			if(mav2 .le. m1+m2)then
				kstop2 = min(klim, mav2-2*m1)
			else
				kstop2 = min(klim, 2*m2-mav2)
			end if
			if(mod(mav2-kstop2,2).ne.0) kstop2 = kstop2 - 1
			kstop1 = - kstop2
c			3.1) Poloidal Fourier tf. of diel. coefficients:
c			call polfft(v2(0,1), npfft+1, .false., idum, .false., 1, trafo, ipath
c     ;		, onlyab_now)
			call polfft(v2(0,1), npfft+1, onespe, icurspe, .false., 1
     ;					,trafo, ipath, onlyab_now)
c			NB ig loop sets ipath to 1. ipath=2 on exit bypasses some inits in polfft at subsequent calls.
cPL			call polfft(v2(0,1), npfft+1, .false., idum, .false., 1, trafo, ipath
c     ;		, nfft, onlyab)
 
			do k = kstop1, kstop2, 2  ! k = m_i-m_j loop - - - - - - -
			   mpk = (mav2 + k) / 2
                 mpkr = mpk + 1 - m1
                 iplb = iplb + 1
c                Considers all blocks sing:
                 iplbs = iplbs + 1
c                kr = k - kstop1 + 1
                 m = mpk - k
                 mr = mpkr - k
c                3.2) Blocks are computed for current Gauss point and copied into VMAT:
                 call plasmb(vmat(1,1,iplb), vmat(1,1,iplbs), k, k, v2(0,1), 0
     ;                      , npfft+1, 1, .false.)
cPL     ;        , npfft+1, 1, nfft, .false.)
c                 write(603,*)mav2,k,vmat(1,1,1)
              end do ! end of k loop - - - - - - - - - - - - - - - - - -

		 end do  ! end of mav2 loop + + + + + + + + + + + + + + + + + + + + + +
 
		 if(.not.coldpl(ireg) .and. flrops(ireg))then
c			There are only 'singular' plasma blocks in this case (see PLASMB)
			iplbmu = iplbs
		 else
			iplbmu = iplb
		 end if

        end if ! cokpco = false

c	  END case of standard angle variables ------------------------------------


c	  4) ??? ------------------------------------------------------------------ 

c       To be done above (polfft, plasmb) in noncircular:
        if(circ .and. eletyp.eq.'M23')then
c          nwkbl = iplbmu
c          Workspace to mul3bl: amount of free 6*6 blocks in totop(12,12,...):
           nwkbl = 4 * (bllen - 2 * nset0)
           if(nwkbl .le. 0) write(6,*)'error asplas: nwkbl=', nwkbl
cPL25/8/04:
	     call mul3bl_fast(bcih(1,1,ig), vmat, bci(1,1,ig), iplbmu)
cPL           call mul3bl_fast(bcih(1,1,ig), vmat, bci(1,1,ig), iplbmu
cPL     ;               , totop(1,1,2*nset0+1), nwkbl)
        end if
c       no problem with cold plasma; quadratic elt. to extend for hot(erho')

c	  END of ??? --------------------------------------------------------------

c	  5) Store results in totop (mu3tf2bl)-------------------------------------
cERN
cERN	  tempo0 = etime (T1)
cERN	  print *, '		 --> (ASPLAS.f) Begin mu3tf2bl:', tempo0
	 
        tem3 = rle * wga(ig) * fac
cERN	  call mu3tf2bl(totop, tem3, ldif(1,1,ig), vmat, iplbmu, idllo)
c        print*, iplb, iplbs, iplbmu,idllo
        call mu3tf2bl_fast(totop, tem3, ldif(1,1,ig), vmat, iplbmu, idllo)
c		   subroutine mu3tf2bl(r, al, a, b, nblo, icola)
c		   r:= r + al * trans(a) * b * a
cERN
cERN	  tempo = etime (T1)
cERN	  print *, '		 --> (ASPLAS.f)   End mu3tf2bl:', tempo
cERN	  print *, '                        Time used ', tempo-tempo0

c	  END of Store results in totop (mu3tf2bl)---------------------------------


 717  continue  ! End of first element Gauss points loop
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      iplb = 0
      iplbs = nset0
      m1 = minf(iel)
      m2 = msup(iel)

        do mav2 = 2*m1, 2*m2
c       ----------------------------------
          if(mav2 .le. m1+m2)then
          kstop2 = min(klim, mav2-2*m1)
          else
          kstop2 = min(klim, 2*m2-mav2)
          end if
        if(mod(mav2-kstop2,2).ne.0)kstop2 = kstop2 - 1
        kstop1 = - kstop2
 
          do k = kstop1, kstop2, 2
          iplb = iplb + 1
          iplbs = iplbs + 1
          mpk = (mav2 + k) / 2
          mpkr = mpk - m1 + 1
          if(.not.bcinb1)call bouas1('c', nj, npopj, axc, jpop1, jpop2, jn)
 
c          jdeca = ishibc + icolo * ( mpkr - 1 )
          jdeca = icolo * ( mpkr - 1 )
 
c         kr = k - kstop1 + 1
          ideca = jdeca - k * icolo
          m = mpk - k
          mr = mpkr - k

            if(eletyp.eq.'M23')then
            index(icolo+1) = (m - m1) * ibulo + 1
     ;        + (msup(iel-1) - m + 1) * icolo
            jndex(icolo+1) = (mpk - m1) * ibulo + 1
     ;        + (msup(iel-1) - mpk + 1) * icolo
            end if
          call zcopy(144, totop(1,1,iplb), 1, to, 1)

            if(.not.coldpl(ireg) .and. flrops(ireg))then
c Voir facteur yinv??? 
cERN              do j = 1, 12
cERN                do i = 1, 12
cERN                to(i,j) = to(i,j) + totop(i,j,iplbs)
			to(1:12,1:12) = to(1:12,1:12) + totop(1:12,1:12,iplbs)
cERN                end do
CERN              end do
            end if

            if(.not.bcinb1)then
c           +++++++++++++++++++
            call bouas1('R', ni, npopi, axl, ipop1, ipop2, in)
              do ip = 1, npopi
                do j = 1, idllo
                to(ipop1(ip),j) = to(ipop1(ip),j) + to(ipop2(ip),j) * axl(ip)
                end do
              end do
              do jp = 1, npopj
                do i = 1, idllo
                to(i,jpop1(jp)) = to(i,jpop1(jp)) + to(i,jpop2(jp))*axc(jp)
                end do
              end do
 
              do j = 1, nj
              jd = jdeca + jndex(jn(j))
                do i = 1, ni
                id = ideca + index(in(i))
                ael(id,jd) = ael(id,jd) + to(in(i),jn(j))
                end do
                do i = icpblo + 1, idllo
                id = ideca + index(i)
                ael(id,jd) = ael(id,jd) + to(i,jn(j))
                end do
              end do
 
              do j = icpblo+1, idllo
              jd = jdeca + jndex(j)
                do i = 1, ni
                id = ideca + index(in(i))
                ael(id,jd) = ael(id,jd) + to(in(i),j)
                end do
                do i = icpblo+1, idllo
                id = ideca + index(i)
                ael(id,jd) = ael(id,jd) + to(i,j)
                end do
              end do
      
            else
c           ++++
C     New option at axis for general geometry; essential (+ nat. if required)
C     conds. enforced by call to AXISBC in GENERA3
c     Modif. by hot plasma to implement in AXISBC!
C     
              do j = 1, idllo
              jd = jdeca + jndex(j)
                do i = 1, idllo
                id = ideca + index(i)
                ael(id,jd) = ael(id,jd) + to(i,j)
                end do
              end do
        
            end if
c           ++++++
 
          end do
        end do
c       ------

      return
c	====================================================
      end if ! (iel = 1) END of first element computations
c     ====================================================


c     Other elements:

      mr = 1
      ielm = iel
      if(nmode(iel-1).gt.nmode(iel)) ielm = iel-1
      nmodc = nmode(ielm)

      m1 = minf(ielm)
      m2 = msup(ielm)
 
      if(klim.ge.nmodc-1)then
         nset0 = nmodc**2
      else
         nset0 = (klim+1)*(2*nmodc-klim) - nmodc
      end if

c     Terms will be summed into totop:

      if(giplas)then
c        Initialize the array:
         call zset(144*2*nset0, czero, totop, 1)
c        Gaussian quadrature over the element
         intab = 1 + (iel - 1) * (ngauss + 1) + ireg
c        Upper limit of radial loop over current element (ngauss internal Gauss points):
         ig2 = ngauss
c        Number of nodes considered and step size for radial index:
         nno = 1
         intast = 1
      else
c        Using analytical integrals of products of basis functions and linear
c        interpolation of plasma contribution over the element:
c        This prevents looping over Gauss points:
         ig2 = 1
         if(iel.eq.ifiel(ireg) .or. iel.eq.2) then
c           intab on left node; comput. for left and right nodes.
c           Assume abscissae in eqt etc unchanged.
            intab = (iel - 1) * (ngauss + 1) + ireg
c           Initialize array totop:
            call zset(144*2*nset0, czero, totop, 1)
c           Number of nodes considered (left and right) and step size for radial index (skips internal points):
            nno = 2
            intast = ngauss + 1
         else
c           intab on right node; assume abscissae in eqt etc unchanged:
            intab = iel * (ngauss + 1) + ireg
c           Deal with right node only, assuming contribution of left one already stored in totop 
c           at the end of previous element assembly. Hence do not reset totop to zero in this case.
            nno = 1
c           here a dummy value:
            intast = 1
         end if
      end if ! (giplas = true)

      if(circ .and. eletyp.eq.'M23') call loctri(iel, giplas)
      
	intab = intab - intast
c	=================================================
      do ino = 1, nno ! Loop over current element nodes
c	=================================================

c         Index ipro is used for accessing analytical integrals in array prodin:
          if(nno.eq.1)then
c         ipro only points on right node (corresp. integrand includes ksi^1): 
          ipro = 1
	    else
c         ipro successively points on left and right node (respective integrands include 1-ksi and ksi): 
          ipro = ino - 1
	    end if

        do 719 ig = 1, ig2  ! Other elements Gauss points loop
c       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        intab = intab + intast
        y = abscno(intab)
        iplb = 0
        iplbs = nset0
        ipath = 1


c	  1) Case of constant k// coordinates -----------------------------------
 
        if(cokpco) then        

           print*, '(ASPLAS) Not ready for cokpco in this version!'

cJACc	     tempo0 = etime (T1)
cJACc	     print *, '		---> (ASPLAS.f) Begin FOUGDR + M12COKPCO:', tempo0

c          1.1) Create folder for data output at radial point: '../../M12run/rXXXX/'    
cJAC	     if (WRITE_OUTPUT) then
cJAC	         call int_to_string4 (nint(1000*abscis(intab)), rrr)	
cJAC	         COKFOLDER = '../../M12run/r' // rrr
cJAC	         shell_comm = 'mkdir ..\..\M12run\r' // rrr
cJAC	         resp = SYSTEM(shell_comm)
cJAC	         if (resp .eq. -1) then
cJAC	            print *, 'Error creating folder: ', COKFOLDER
cJAC		        stop
cJAC	         end if
cJAC		     write(5555,'(i8, g15.6, A20)'), intab, abscis(intab), COKFOLDER
cJACc	         print *, intab, abscis(intab), COKFOLDER
cJAC	     end if  ! WRITE_OUTPUT = true

c		 1.2) Geometrical couplings for plasma response in 'constant k//' coord.:
c               This call must be reviewed when non-Maxwellian behaviour will be 
c               implemented, because Max. and non-Max. need different geom. coefficients.
c               NEW: If reading M12 table files avoid this call!
cJAC		 if (READ_TABLES .eq. .FALSE.) then
cJACcERN			 call fougdr(1)
cJAC			 call fougdr_ern(1)     ! fougdr_ern: vectorized version (checked OK!)
cJAC		 else
cJAC		     if (WRITE_OUTPUT) then
cJAC			     shell_comm2 = 'copy ..\..\M12tables\r' // rrr // '\Plm.dat ' 
cJAC    ;                             // '..\..\M12run\r' // rrr 
cJAC				 resp = SYSTEM(shell_comm2)
cJAC			 end if
cJAC		 end if
c			  The geom. coef. are calculated using FFT with respect to variable Chi and 
c			  are stored in the matrix gcdr(L,mdiff), where L is the Chi-harmonic index
c			  and mdiff = m1-m2 is the poloidal mode difference.
c				-> L range     : [-npfft/2 .... -1 0 +1 .... +npfft/2]
c				-> mdiff range : [-klim .... -1 0 +1 .... +klim]	

c		 1.3) Calculate M12 matrices (already summed for all species)
c			  The results are stored in the common's (COKPCO.COPY):
c				m12left   (k_index, m_index) : for p = +1 polarization
c				m12right  (k_index, m_index) : for p = -1 polarization
c				m12landau (k_index, m_index) : for p = 0  polarization
c						-> k_index : index for k// (or A)
c						-> m_index : index for m1-m2 (poloidal couplings)

cJAC		 call M12COKPCO(onespe, icurspe)

c			  NB: If onespe = TRUE, M12COKPCO will return the response of a single
c				  species (identified as icurspe) in m12left, m12right and m12landau.
c				  Otherwise, it returns the standard response summed over all species.

c		 1.4) Writing M12 elements to VMAT matrix 
c			  To be careful because data in m12(i,j) matrices and in VMAT have different 
c			  structures:
c				-> k// (Acoef) index : 'i' in m12 is 'mav2-2*m1+1' in VMAT
c				-> mdiff index : 'j' in m12 is 'k+klim+1' in VMAT

cERN				iplb=0  ! WARNING: CHECK this %%%%%%%%%%%%%%%%%%%%

c		 Loop over average pol. mode : (this loop is very fast)
cJAC		 do mav2 = 2*m1, 2*m2 ! mav2 = 2*(m_i+m_j) loop ++++++++++++++++
cJAC              k_index = mav2-2*m1+1
cJAC			if(mav2 .le. m1+m2)then
cJAC                 kstop2 = min(klim, mav2-2*m1)
cJAC			else
cJAC                 kstop2 = min(klim, 2*m2-mav2)
cJAC			end if
cJAC			if(mod(mav2-kstop2,2).ne.0) kstop2 = kstop2 - 1
cJAC			kstop1 = - kstop2
cJAC             do k = kstop1, kstop2, 2 ! k = m_i-m_j loop ------------
cJAC                 iplb = iplb + 1
cJACcERN		       k_index = mav2-2*m1+1
cJAC			   m_index = k+klim+1
cJAC	           vmat(1,1,iplb) = m12left  (k_index, m_index) ! p = +1
cJAC	           vmat(3,3,iplb) = m12right (k_index, m_index) ! p = -1
cJAC	           vmat(5,5,iplb) = m12landau(k_index, m_index) ! p = 0 

cJAC             end do  ! end of k loop --------------------------------
          
cJAC		 end do ! end of mav2 loop ++++++++++++++++++++++++++++++++++++

	     iplbmu = iplb ! WARNING: CHECK this
c	     tempo = etime (T1)
c	     print *, '		---> (ASPLAS.f)   End FOUGDR + M12COKPCO:', tempo
c	     print *, '                        Time used ', tempo-tempo0

	  end if ! cokpco = TRUE


c	  END Case of constant k// coordinates ------------------------------------


c	  2) Non-Maxwellians ???? -------------------------------------------------

        if(nspgdr.gt.0)then
           print*,'(asplas) NSPGR NOT READY!'
cJACc       Species with a general dielectric response:
cJAC           if(onespe)then
cJACcPL@          to do: allow evaluation 1 species at a time:
cJAC              print *,'ASPLAS gdr 2nd call: must implement calculation
cJAC     ;                 1 species at a time'
cJAC	        stop
cJAC           end if
cJACc          Compute x, v meshes at current radius and interpolate equil. distr. f0:
cJAC           call intf0
cJACc          Add contribution to VMAT blocks. Currently 1 species allowed.
cJAC           isp = 1
cJACc          do isp = 1, nspgdr
cJAC           call diresp(isp)
cJACc          end do
        end if

c	  END of Non-Maxwellians ???? ---------------------------------------------

c	  3) Case of standard angle variables -------------------------------------
c          Maxwellian species:

        if(.not.cokpco)then

c	     tempo0 = etime (T1)
c	     print *, '		---> (ASPLAS.f) Begin POLFFT + PLASMB:', tempo0

c         Loop over average pol. mode (mav2 = 2 maverage):
          do mav2 = 2*m1, 2*m2
c         --------------------
		   if(mav2 .le. m1+m2)then
			  kstop2 = min(klim, mav2-2*m1)
             else
			  kstop2 = min(klim, 2*m2-mav2)
		   end if
		   if(mod(mav2-kstop2,2).ne.0) kstop2 = kstop2 - 1
             kstop1 = - kstop2
 
c            Plasma terms poloidal integration (all species together):
             call polfft(v2(0,1), npfft+1, onespe, icurspe, .false., 1
     ;                  ,trafo, ipath, onlyab_now)
c            call polfft(v2(0,1), npfft+1, .false., idum, .false., 1, trafo, ipath
c     ;      , onlyab_now)
cPL     ;    , nfft, onlyab)
c             write(603,*)mav2,v2(0:10,1)
             do k = kstop1, kstop2, 2
			  iplb = iplb + 1
			  iplbs = iplbs + 1
			  mpk = (mav2 + k) / 2
			  m = mpk - k
			  mr = m + 1 - m1
			  mpkr = mr + k
c			  Blocks are computed for current Gauss point and copied into VMAT:
			  call plasmb(vmat(1,1,iplb), vmat(1,1,iplbs), k, k, v2(0,1), 0
     ;		            , npfft+1, 1, .false.)
cPL     ;					, npfft+1, 1, nfft, .false.)
c                          write(603,*),mav2,k,vmat(3,3,iplb)
             end do  ! k
          
		end do  ! mav2
c         ------

c	     tempo = etime (T1)
c	     print *, '		---> (ASPLAS.f)   End POLFFT + PLASMB:', tempo
c	     print *, '                        Time used ', tempo-tempo0

          if(.not.coldpl(ireg) .and. flrops(ireg))then
c         There are only 'singular' plasma blocks in this case (see PLASMB)
          iplbmu = iplbs
          else
          iplbmu = iplb
          end if
        
	  end if     ! cokpco = false

c	  END Case of standard angle variables ------------------------------------
c             write(603,*)iplb,iplbs
c             write(603,*)vmat(3,3,1:10)
cERN	  DEBUGGING
c	  if(abscis(intab)>0.3) then
c	     write(11111,"(f15.5)"), dreal(vmat(1,1,:))
c	     write(22222,"(f15.5)"), dreal(vmat(3,3,:))
c	     write(33333,"(f15.5)"), dreal(vmat(5,5,:))
c	     write(11110,"(f15.5)"), dimag(vmat(1,1,:))
c	     write(22220,"(f15.5)"), dimag(vmat(3,3,:))
c	     write(33330,"(f15.5)"), dimag(vmat(5,5,:))
c	     stop
c	  end if


c             stop

c	  4) ??? ------------------------------------------------------------------ 

c       To be done above (polfft, plasmb) in noncircular:
        if(circ .and. eletyp.eq.'M23')then
c       nwkbl = iplbmu
c       Workspace to mul3bl: amount of free 6*6 blocks in totop(12,12,...):
        nwkbl = 4 * (bllen - 2 * nset0)
        if(nwkbl .le. 0)write(6,*)'Error asplas: nwkbl=', nwkbl
cPL25/8/04:
        call mul3bl_fast(bcih(1,1,ig), vmat, bci(1,1,ig), iplbmu)
cPL        call mul3bl_fast(bcih(1,1,ig), vmat, bci(1,1,ig), iplbmu,
cPL     ;              totop(1,1,2*nset0+1), nwkbl)
        end if

c	  End of ??? -------------------------------------------------------------- 


c	  5) Store results in totop (mu3tf2bl) ------------------------------------

        if(giplas) then ! Gaussian integration over element + + + + + + + + + + + 
          
		tem3 = rle * wga(ig) * fac
cERN
cERN		tempo0 = etime (T1)
cERN		print *, '		---> (ASPLAS.f) Begin mu3tf2bl:', tempo0          

cERN 		call mu3tf2bl(totop, tem3, ldif(1,1,ig), vmat, iplbmu, idllo)
		call mu3tf2bl_fast(totop, tem3, ldif(1,1,ig), vmat, iplbmu, idllo)
c			 subroutine mu3tf2bl(r, al, a, b, nblo, icola)
c			 r:= r + al * trans(a) * b * a

cERN
cERN		tempo = etime (T1)
cERN		print *, '		---> (ASPLAS.f)   End mu3tf2bl:', tempo
cERN		print *, '                        Time used ', tempo-tempo0

        else ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  
             ! Analytical integration over element, using linear interpolation of equilibrium (under development)
            do ic1 = 1, nficom                  ! Loop over test function components
            do l1 = 0, 1                        ! Loop over test function derivatives (0th and 1st)
            ir1 = 2 * ic1 - 1 + l1              !   Corresponding row index in vmat
            it1 = ibafty(ic1,ielet)
              do ib1 = 1, nbafel(ic1,ielet)     ! Loop over basis functions of current test function component
              jb1 = jbfloc(ib1,ic1)
              id1 = ideriv(ib1,it1) - l1
                do ic2 = 1, nficom              ! Loop over field components
                do l2 = 0, 1                    ! Loop over field derivatives (0th and 1st)
                ir2 = 2 * ic2 - 1 + l2          !   Corresponding column index in vmat
                it2 = ibafty(ic2,ielet)
                  do ib2 = 1, nbafel(ic2,ielet) ! Loop over basis functions of current field component
                  jb2 = jbfloc(ib2,ic2)
                  id2 = ideriv(ib2,it2) - l2
c                                                 prodin contains analytical integrals of products of basis functions on [0,1]
c                                                 rlpow  contains powers of element length
                  cri = fac * (rlpow(id1 + id2 + 1) * prodin(ib1,l1,it1,ib2,l2,it2,ipro))  
                    do i = 1, iplbmu            ! Add each analytical contribution
                    totop(jb1,jb2,i) = totop(jb1,jb2,i) + cri * vmat(ir1,ir2,i)
                    end do
                  end do
                end do
                end do
              end do
            end do
            end do
	  end if ! (giplas = true) + + + + + + + + + + + + + + + + + + + + +

c	  End of Store results in totop (mu3tf2bl) --------------------------------

  
 719    continue ! End of other elements Gauss points loop
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c	==================================
      end do ! (ino) cERN ??????????????
c     ================================== 




cERN	cccccccccccc EVERYTHING READY cccccccccccccccccccccccc


c     Assembly into AEL:
      entry passem

      iplb = 0
      iplbs = nset0
        do 107 mav2 = 2*m1, 2*m2
c       ------------------------
          if(mav2.le.m1+m2)then
          kstop2 = min(klim, mav2-2*m1)
          else
          kstop2 = min(klim, 2*m2-mav2)
          end if
        if(mod(mav2-kstop2,2).ne.0)kstop2 = kstop2 - 1
        kstop1 = - kstop2
 
          do k = kstop1, kstop2, 2
          iplb = iplb + 1
          iplbs = iplbs + 1
          mpk = (mav2 + k) / 2
          mpkr = mpk - m1 + 1
          jnl = nmode(iel-1).gt.0 .and. mpk.ge.minf(iel-1) .and.
     ;          mpk.le.msup(iel-1)
          jnr = nmode(iel).gt.0 .and. mpk.ge.minf(iel) .and.
     ;          mpk.le.msup(iel)
          jndext = nmodc*ibulo
          if(jnl)jndext = jndext + (msup(iel-1) - mpk + 1) * icolo
          if(jnr)jndext = jndext + (mpk - minf(iel)) * icolo
            if(jnl)then
            jdeca = (mpk - minf(iel-1)) * icolo
            else
            jdeca = nmode(iel-1) * icolo
            end if

            do i = 1, icolo
            jndex(i+icpblo) = i + jndext
            end do

            if(eletyp.eq.'M23')then
            jndex(icolo+1) = (mpk - m1) * ibulo + 1
            if(jnl)jndex(icolo+1) = jndex(icolo+1)
     ;      + (msup(iel-1) - mpk + 1) * icolo
            end if

            if(jnl)then
            j1 = 1
	      else
            j1 = icolo + 1
	      end if
	      if(jnr)then
            j2 = idllo
            else
            j2 = icpblo
	      end if
 
c         kr=k-kstop1+1
          m = mpk - k
          mr = mpkr - k
          inl = nmode(iel-1).gt.0 .and. m.ge.minf(iel-1) .and.
     ;          m.le.msup(iel-1)
          inr = nmode(iel).gt.0 .and. m.ge.minf(iel) .and.
     ;          m.le.msup(iel)
          indext = nmodc * ibulo
          if(inl)indext = indext + (msup(iel-1) - m + 1) * icolo
          if(inr)indext = indext + (m - minf(iel)) * icolo

            do i = 1, icolo
            index(i+icpblo) = i + indext
            end do

            if(eletyp.eq.'M23')then
            index(icolo+1) = (m - m1) * ibulo + 1
            if(inl)index(icolo+1) = index(icolo+1)
     ;       + (msup(iel-1) - m + 1) * icolo
            end if

            if(inl)then
            i1 = 1
            ideca = (m - minf(iel-1)) * icolo
            else
            i1 = icolo + 1
            ideca = nmode(iel-1) * icolo
            end if

            if(inr)then
            i2 = idllo
            else
            i2 = icpblo
            end if
 
cERN            do j = j1, j2
cERN            jd = jdeca + jndex(j)
cERN              do i = i1, i2
cERN                 id = ideca + index(i)
cERN                 ael(id,jd) = ael(id,jd) + totop(i,j,iplb)
cERN              end do
cERN            end do
		  ael(ideca+index(i1:i2) , jdeca+jndex(j1:j2)) = 
     ;      ael(ideca+index(i1:i2) , jdeca+jndex(j1:j2)) + totop(i1:i2,j1:j2,iplb)


            if(.not.coldpl(ireg) .and. flrops(ireg))then
c Voir: Voir facteur yinv??? 
cERN              do j = j1, j2
cERN              jd = jdeca + jndex(j)
cERN                do i = i1, i2
cERN                   id = ideca + index(i)
cERN                   ael(id,jd) = ael(id,jd) + totop(i,j,iplbs)
cERN                end do
cERN              end do
		    ael(ideca+index(i1:i2) , jdeca+jndex(j1:j2)) = 
     ;        ael(ideca+index(i1:i2) , jdeca+jndex(j1:j2)) + totop(i1:i2,j1:j2,iplbs)
            end if
          end do
 107    continue
c       --------
      if(lefnod)return

      if(.not.giplas .and. iel.gt.1 .and. iel.ne.ilael(ireg))then
c     ===========================================================
c     Use results of right node to prepare assembly of left node in next
c     element. Assumes same number of modes in all elements.
c     See later for first element.
        do i = -2, 3
        rlpow(i) = FL(iel+1) ** i
        end do

      call zset(2*144*nset0, czero, totop, 1)
 
        do ic1 = 1, nficom
        do l1 = 0, 1
        ir1 = 2 * ic1 - 1 + l1
        it1 = ibafty(ic1,ielet)
          do ib1 = 1, nbafel(ic1,ielet)
          jb1 = jbfloc(ib1,ic1)
          id1 = ideriv(ib1,it1) - l1
            do ic2 = 1, nficom
            do l2 = 0, 1
            ir2 = 2 * ic2 - 1 + l2
            it2 = ibafty(ic2,ielet)
              do ib2 = 1, nbafel(ic2,ielet)
              jb2 = jbfloc(ib2,ic2)
              id2 = ideriv(ib2,it2) - l2
              cri = fac * 
     ;        (rlpow(id1 + id2 + 1) * prodin(ib1,l1,it1,ib2,l2,it2,0))
                do i = 1, iplbmu
                totop(jb1,jb2,i) = totop(jb1,jb2,i) + cri * VMAT(ir1,ir2,i)
                end do
              end do
            end do
            end do
          end do
        end do
        end do      
      end if 
c     ======

      return

cERN	Close cokpco files
	if(WRITE_OUTPUT) close(5555) ! Closing index_radii.dat file	  
	if(WRITE_REPORT) close(7777) ! Closing cokpco_report.txt file

      end
