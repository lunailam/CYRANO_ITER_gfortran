      subroutine asplas(onespe, icurspe)

c      SUBROUTINE asplas

	use dfport		! necessary for 'etime' function

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

	integer :: k_index, m_index
	logical :: resp
	real :: T1(2)	            !   Variables for evaluating
	real :: tempo, tempo0       !   the elapsed CPU time

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision 
     ;  cabs1, rle
     ;, second
     ;, rm

      complex*16 zdum

      complex*16 axl(3), axc(3), tem3, fac
     ;, v2(0:npfft,19)
     ;, to(12,12), toto(12*12), tos(12,12), totos(12*12)
     ;, big(12,12), bigs(12,12), big1(12*12), big1s(12*12)
     ;, cri

      save nset0, m1, m2, jndex, index, ielet, eletyp, idllo, icolo, ibulo
     ;, icpblo

      equivalence (toto,to), (totos,tos), (big1,big), (big1s,bigs)

      external second, nsum, i1st

      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))

c	tempo0 = etime (T1)
 
c     FAC: common factor to all matrix elements
c     remember field spectrum has dimensions:
c       - V/m in toroidal geometry
c       - V   in a cylinder
      fac = glofac * twopi ** 2
      if(.not.cyl)fac = fac * ra

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
        do  i = 1, icolo
        index(i) = i
        jndex(i) = i
        index(i+icpblo) = i + indext
        jndex(i+icpblo) = i + indext
        end do 
  
      if(iel.eq.1)then
c     ================
c     First element at magnetic axis or crown inner boudary
c     !nmode(0) = nmode(1) assumed!

      intab = 1 + (iel-1) * (ngauss+1) + ireg

      if(circ .and. eletyp.eq.'M23')call loctri(iel,.true.)
 
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
  
      do 717 ig = 1, ngauss
c     ---------------------
      intab = ig + (iel-1) * (ngauss+1) + ireg
      y = abscno(intab)
      iplb = 0
      iplbs = nset0
      ipath = 1

ccccccccccccccccccccccccccccccccccc ERNESTO ccccccccccccccccccccccccccccccccccc

! ------------------- First element (mag.axis) ------------------------------- !

      if(cokpco) then  !  k// = const. coordinates
	
	tempo0 = etime (T1)

        open (UNIT = 7, FILE = 'cokpco_report.txt', STATUS = "REPLACE", 
     ;        ACTION = "WRITE")

	  if (WRITE_REPORT) then
	    write(7,*)'=================================================='
		write(7,*)'Calculus of plasma response M12 for k// = constant'
		write(7,*)'=================================================='
		write(7,*)
	    write(7,*)'Number of poloidal mode diff.(m1-m2) =',2*klim+1
		write(7,*)'Number of A (or k//) per element =',2*(modva2-modva1)+1
	    if (READ_TABLES .eq. .FALSE.) then
		   write(7,*)'Method : Calculating M12 values explicitly'
	    else
		   write(7,*)'Method : Interpolating M12 values from table files'
	    end if
	  end if
	  
	  if (WRITE_SCREEN) then
	    print *,'=================================================='
	    print *,'Calculus of plasma response M12 for k// = constant'
	    print *,'=================================================='
	    print *
	    print *,'Number of poloidal mode diff.(m1-m2) =',2*klim+1
	    print *,'Number of A (or k//) per element =',2*(modva2-modva1)+1
	    if (READ_TABLES .eq. .FALSE.) then
		   print *, 'Method : Calculating M12 values explicitly'
	    else
		   print *, 'Method : Interpolating M12 values from table files'
	    end if
	  end if

c     1) Geometrical couplings for dielectric response in 'constant k//' coord.:
c        This call must be reviewed when non-Maxwellian behaviour will be 
c        implemented, because Max. and non-Max. need different geometric coefficients.
c        NEW: If reading M12 table files we don't need this call!!
	  if (READ_TABLES .eq. .FALSE.) then
	  call fougdr(1)
	  end if
c	   The geom. coef. are calculated using FFT with respect to variable Chi and 
c        are stored in the matrix gcdr(L,mdiff), where L is the Chi-harmonic index
c        and mdiff = m1-m2 is the poloidal mode difference.
c	      -> L range     : [-npfft/2 .... -1 0 +1 .... +npfft/2]
c	      -> mdiff range : [-klim .... -1 0 +1 .... +klim]	

c     2) Calculate M12 matrices (already summed for all species)
c        The results are stored in the common's (COKPCO.COPY):
c           m12left   (k_index, m_index) : for p = +1 polarization
c		  m12right  (k_index, m_index) : for p = -1 polarization
c           m12landau (k_index, m_index) : for p = 0  polarization
c	        -> k_index : index for k// (or A)
c             -> m_index : index for m1-m2 (poloidal couplings)

cPL@	  allow evaluation 1 species at a time:
	  call M12COKPCO(onespe, icurspe)
c	  NB: If onespe = TRUE, M12COKPCO will return the response of a single
c		  species (identified as icurspe) in m12left, m12right and m12landau.
c		  Otherwise, it returns the standard response summed over all species.

c	3) Write M12 elements to VMAT matrix --------------------------------
c       To be careful because data in m12 matrices and in VMAT have different 
c       structures:
c	    -> k// (or A) index : 'k_index' in m12 is 'mav2-2*m1+1' in VMAT
c	    -> mdiff index : 'm_index' in m12 is 'k+klim+1' in VMAT
	  if (WRITE_SCREEN) then
           print *
	     print *, 'Writing M12 elements to VMAT matrix'
	  end if

c       Loop over average pol. mode :
        do mav2 = 2*m1, 2*m2 ! mav2 = 2*(m_i+m_j) loop ++++++++++++++++
             
	     if(mav2 .le. m1+m2)then
              kstop2 = min(klim, mav2-2*m1)
           else
              kstop2 = min(klim, 2*m2-mav2)
           end if
           if(mod(mav2-kstop2,2).ne.0)kstop2 = kstop2 - 1
           kstop1 = - kstop2

           do k = kstop1, kstop2, 2 ! k = m_i-m_j loop ------------
              iplb = iplb + 1
		    k_index = mav2-2*m1+1
			m_index = k+klim+1

	        vmat(1,1,iplb) = m12left  (k_index, m_index) ! p = +1
	        vmat(3,3,iplb) = m12right (k_index, m_index) ! p = -1
	        vmat(5,5,iplb) = m12landau(k_index, m_index) ! p = 0 

           end do  ! end of k loop --------------------------------
          
	  end do ! end of mav2 loop ++++++++++++++++++++++++++++++++++++

	iplbmu = iplb
	end if  ! cokpco=true

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if(nspgdr.gt.0)then
          if(onespe)then
cPL@      to do: allow evaluation 1 species at a time:
          print *,'ASPLAS gdr: must implement calculation 1 species at a time'
	    stop
          end if
          
c       Species with a general dielectric response:
c       Compute x, v meshes at current radius and interpolate equil. distr. f0:
        call intf0
c       Add contribution to VMAT blocks. Currently 1 species allowed.
        isp = 1
c         do isp = 1, nspgdr
        call diresp(isp)
c         end do
        end if

c       Maxwellian species:
        if(.not.cokpco)then
c       Case of standard angle variables

c         Loop over average pol. mode (mav2 = 2 maverage):
          do mav2 = 2*minf(iel), 2*msup(iel)
            if(mav2 .le. minf(iel)+msup(iel))then
            kstop2 = min(klim, mav2-2*minf(iel))
            else
            kstop2 = min(klim, 2*msup(iel)-mav2)
            end if
          if(mod(mav2-kstop2,2).ne.0)kstop2 = kstop2 - 1
          kstop1 = - kstop2
c         Poloidal Fourier tf. of diel. coefficients:
c          call polfft(v2(0,1), npfft+1, .false., idum, .false., 1, trafo, ipath
c     ;  , onlyab_now)
          call polfft(v2(0,1), npfft+1, onespe, icurspe, .false., 1
     ;               ,trafo, ipath, onlyab_now)
c         NB ig loop sets ipath to 1. ipath=2 on exit bypasses some inits in polfft at subsequent calls.
cPL          call polfft(v2(0,1), npfft+1, .false., idum, .false., 1, trafo, ipath
c     ;  , nfft, onlyab)
 
            do k = kstop1, kstop2, 2
            mpk = (mav2 + k) / 2
            mpkr = mpk + 1 - m1
            iplb = iplb + 1
c           Considers all blocks sing:
            iplbs = iplbs + 1
c           kr = k - kstop1 + 1
            m = mpk - k
            mr = mpkr - k
c           Blocks are computed for current Gauss point and copied into VMAT:
            call plasmb(vmat(1,1,iplb), vmat(1,1,iplbs), k, k, v2(0,1), 0
     ;    , npfft+1, 1, .false.)
cPL     ;    , npfft+1, 1, nfft, .false.)
            end do
          end do
 
          if(.not.coldpl(ireg) .and. flrops(ireg))then
c         There are only 'singular' plasma blocks in this case (see PLASMB)
          iplbmu = iplbs
          else
          iplbmu = iplb
          end if

        end if ! cokpco=false

c         To be done above (polfft, plasmb) in noncircular:
          if(circ .and. eletyp.eq.'M23')then
c         nwkbl = iplbmu
c         Workspace to mul3bl: amount of free 6*6 blocks in totop(12,12,...):
          nwkbl = 4 * (bllen - 2 * nset0)
          if(nwkbl .le. 0)write(6,*)'error asplas: nwkbl=', nwkbl
          call mul3bl(bcih(1,1,ig), vmat, bci(1,1,ig), iplbmu
     ;             , totop(1,1,2*nset0+1), nwkbl)
          end if
c     no problem with cold plasma; quadratic elt. to extend for hot(erho')
 
      tem3 = rle * wga(ig) * fac
      call mu3tf2bl(totop, tem3, ldif(1,1,ig), vmat, iplbmu, idllo)

 717  continue
c     --------
 
      iplb = 0
      iplbs = nset0
      m1 = minf(iel)
      m2 = msup(iel)

        do mav2 = 2*minf(iel), 2*msup(iel)
c       ----------------------------------
          if(mav2 .le. minf(iel)+msup(iel))then
          kstop2 = min(klim, mav2-2*minf(iel))
          else
          kstop2 = min(klim, 2*msup(iel)-mav2)
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
              do j = 1, 12
                do i = 1, 12
                to(i,j) = to(i,j) + totop(i,j,iplbs)
                end do
              end do
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
C     conds. enforced by call to AXISBC in GENERA2
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
      end if
c     ======

c     Other elements:

      mr = 1
      ielm = iel
      if(nmode(iel-1).gt.nmode(iel))ielm = iel-1
      nmodc = nmode(ielm)

      m1 = minf(ielm)
      m2 = msup(ielm)
 
        if(klim.ge.nmodc-1)then
        nset0 = nmodc**2
        else
        nset0 = (klim+1)*(2*nmodc-klim) - nmodc
        end if

c     Terms will be summed into totop:
      call zset(144*2*nset0, czero, totop, 1)

        if(giplas)then
c       Gaussian quadrature over the element
        intab = 1 + (iel - 1) * (ngauss + 1) + ireg
        ig2 = ngauss
        nno = 1
        intast = 1
        else
c       Using analytical integrals of products of basis functions and linear
c       interpolation of plasma contribution over the element:
        ig2 = 1
          if(iel.eq.ifiel(ireg) .or. iel.eq.2)then
c         intab on left node; comput. for left and right nodes.
c         Assume abscissae in eqt etc unchanged.
          intab = (iel - 1) * (ngauss + 1) + ireg
          nno = 2
          intast = ngauss + 1
          else
c         intab on right node; assume abscissae in eqt etc unchanged:
          intab = iel * (ngauss + 1) + ireg
          ig2 = 1
          nno = 1
c         irrelevant value:
          intast = 1
          end if
        end if

        if(circ .and. eletyp.eq.'M23')call loctri(iel, giplas)
        intab = intab - intast

        do ino = 1, nno 
        ipro = 1
        if(nno.eq.2)ipro = ino - 1

        do 719 ig = 1, ig2
c       ------------------
        intab = intab + intast

        y = abscno(intab)
        iplb = 0
        iplbs = nset0
        ipath = 1

cccccccccccccccccccccccccccccc  ERNESTO  ccccccccccccccccccccccccccccccccccccccc

! --------------------------- Other elements ----------------------------------
      if(cokpco) then        

c     1) Geometrical couplings for dielectric response in 'constant k//' coord.:
c          This call must be reviewed when non-Maxwellian behaviour will be 
c          implemented, because Max. and non-Max. need different geometric coefficients.
c          NEW: If reading M12 table files avoid this call!
	  if (READ_TABLES .eq. .FALSE.) then
	    call fougdr(1)
	  end if
c	     The geom. coef. are calculated using FFT with respect to variable Chi and 
c          are stored in the matrix gcdr(L,mdiff), where L is the Chi-harmonic index
c          and mdiff = m1-m2 is the poloidal mode difference.
c	        -> L range     : [-npfft/2 .... -1 0 +1 .... +npfft/2]
c	        -> mdiff range : [-klim .... -1 0 +1 .... +klim]	

c     2) Calculate M12 matrices (already summed for all species)
c        The results are stored in the common's (COKPCO.COPY):
c           m12left   (k_index, m_index) : for p = +1 polarization
c		  m12right  (k_index, m_index) : for p = -1 polarization
c           m12landau (k_index, m_index) : for p = 0  polarization
c	        -> k_index : index for k// (or A)
c             -> m_index : index for m1-m2 (poloidal couplings)

cPL@    allow evaluation 1 species at a time:
	  call M12COKPCO(onespe, icurspe)
c	  NB: If onespe = TRUE, M12COKPCO will return the response of a single
c		  species (identified as icurspe) in m12left, m12right and m12landau.
c		  Otherwise, it returns the standard response summed over all species.


c	3) Writing M12 elements to VMAT matrix --------------------------------
c       To be careful because data in m12 matrices and in VMAT have different 
c       structures:
c	    -> k// (Acoef) index : 'i' in m12 is 'mav2-2*m1+1' in VMAT
c	    -> mdiff index : 'j' in m12 is 'k+klim+1' in VMAT
	    if (WRITE_SCREEN) then
             print *
	       print *, 'Writing M12 elements to VMAT matrix'
	    end if
ccc	    iplb=0

c       Loop over average pol. mode :
         do mav2 = 2*m1, 2*m2 ! mav2 = 2*(m_i+m_j) loop ++++++++++++++++
             
	      if(mav2 .le. m1+m2)then
               kstop2 = min(klim, mav2-2*m1)
            else
               kstop2 = min(klim, 2*m2-mav2)
            end if
            if(mod(mav2-kstop2,2).ne.0)kstop2 = kstop2 - 1
            kstop1 = - kstop2

            do k = kstop1, kstop2, 2 ! k = m_i-m_j loop ------------
               iplb = iplb + 1
		     k_index = mav2-2*m1+1
			 m_index = k+klim+1

	         vmat(1,1,iplb) = m12left  (k_index, m_index) ! p = +1
	         vmat(3,3,iplb) = m12right (k_index, m_index) ! p = -1
	         vmat(5,5,iplb) = m12landau(k_index, m_index) ! p = 0 

            end do  ! end of k loop --------------------------------
          
	   end do ! end of mav2 loop ++++++++++++++++++++++++++++++++++++

ccc	   iplbmu = iplb
	  
	  close(7) ! Closing cokpco_report.txt file
	
	end if ! cokpco = TRUE
c	print *, 'Stopping at ASPLAS.f line 626'
c	stop   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          if(nspgdr.gt.0)then
c       Species with a general dielectric response:
          if(onespe)then
cPL@      to do: allow evaluation 1 species at a time:
          print *,
     ;'ASPLAS gdr 2nd call: must implement calculation 1 species at a time'
	    stop
          end if
c       Compute x, v meshes at current radius and interpolate equil. distr. f0:
          call intf0
c       Add contribution to VMAT blocks. Currently 1 species allowed.
          isp = 1
c           do isp = 1, nspgdr
          call diresp(isp)
c           end do
          end if

c       Maxwellian species:

        if(.not.cokpco)then
c       Case of standard angle variables
c         Loop over average pol. mode (mav2 = 2 maverage):
          do mav2 = 2*m1, 2*m2
c         --------------------
            if(mav2 .le. m1+m2)then
            kstop2 = min(klim, mav2-2*m1)
            else
            kstop2 = min(klim, 2*m2-mav2)
            end if
          if(mod(mav2-kstop2,2).ne.0)kstop2 = kstop2 - 1
          kstop1 = - kstop2
 
c         Plasma terms poloidal integration (all species together):
          call polfft(v2(0,1), npfft+1, onespe, icurspe, .false., 1
     ;               ,trafo, ipath, onlyab_now)
c          call polfft(v2(0,1), npfft+1, .false., idum, .false., 1, trafo, ipath
c     ;    , onlyab_now)
cPL     ;    , nfft, onlyab)
 
            do k = kstop1, kstop2, 2
            iplb = iplb + 1
            iplbs = iplbs + 1
            mpk = (mav2 + k) / 2
            m = mpk - k
            mr = m + 1 - m1
            mpkr = mr + k
c           Blocks are computed for current Gauss point and copied into VMAT:
            call plasmb(vmat(1,1,iplb), vmat(1,1,iplbs), k, k, v2(0,1), 0
     ;      , npfft+1, 1, .false.)
cPL     ;      , npfft+1, 1, nfft, .false.)
            end do  ! k
          end do  ! mav2
c         ------

          if(.not.coldpl(ireg) .and. flrops(ireg))then
c         There are only 'singular' plasma blocks in this case (see PLASMB)
          iplbmu = iplbs
          else
          iplbmu = iplb
          end if
        end if     ! cokpco = false

c       To be done above (polfft, plasmb) in noncircular:
          if(circ .and. eletyp.eq.'M23')then
c         nwkbl = iplbmu
c         Workspace to mul3bl: amount of free 6*6 blocks in totop(12,12,...):
          nwkbl = 4 * (bllen - 2 * nset0)
          if(nwkbl .le. 0)write(6,*)'Error asplas: nwkbl=', nwkbl
          call mul3bl(bcih(1,1,ig), vmat, bci(1,1,ig), iplbmu,
     ;                totop(1,1,2*nset0+1), nwkbl)
          end if

          if(giplas)then
          tem3 = rle * wga(ig) * fac
          call mu3tf2bl(totop, tem3, ldif(1,1,ig), vmat, iplbmu, idllo)
          else  ! analytical integration, using linear interpolation of equilibrium (under development)
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
     ;            (rlpow(id1 + id2 + 1) * prodin(ib1,l1,it1,ib2,l2,it2,ipro))
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
  
 719    continue
        end do
C       --------  

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
 
            do j = j1, j2
            jd = jdeca + jndex(j)
              do i = i1, i2
              id = ideca + index(i)
              ael(id,jd) = ael(id,jd) + totop(i,j,iplb)
              end do
            end do

            if(.not.coldpl(ireg) .and. flrops(ireg))then
c Voir: Voir facteur yinv??? 
              do j = j1, j2
              jd = jdeca + jndex(j)
                do i = i1, i2
                id = ideca + index(i)
                ael(id,jd) = ael(id,jd) + totop(i,j,iplbs)
                end do
              end do
            end if
          end do
 107    continue
c       --------
      if(lefnod)return

      if(iel.gt.1 .and. .not.giplas .and. iel.ne.ilael(ireg))then
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

      end
