      subroutine foucog(ind)

ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)
      implicit none
      
      integer ind

c     ind = 1:
c     Computes Fourier coeffs. used in curl.curl with IMSL fft (formerly Cray CFFT2)
c
c     ind not = 1:
c     Computes Fourier coeffs. used in curl with IMSL fft (formerly Cray CFFT2)

c     Radial table index INTAB is input in COMSWE.
c     D-SHAPED or general cross section.

      include 'pardim.copy'
      include 'comfou.copy'
      include 'comgeo.copy'
      include 'commod.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comswe.copy'
      include 'comrot.copy'
      include 'comsub.copy'
      include 'comphy.copy'
      include 'comfft.copy'

      logical ax

      integer 
     ;  i, ii, ipo, ic, j, ka, l, syv(ncorot)
     ;, bsy(nficom,6), bssy(nficom,6)
     ;, ncoexp, wri(4)
c     ;, n2, ini

      double precision t, ts, b(nficom,6), bs(nficom,6), cr(npfft+1,ncorot)
     ;, ncrnor, nsrnor, rnllnr, aux1, aux2, aux3, aux4, ontpbt, nrnlpb, rnlsoc
     ;, crnor, srnor, ri, ntni, jacni, ds
     ;, rnor
      double precision sior, si1tn

      complex*16
     ;  cb(nficom,6), cbs(nficom,6), ccr(npfft+1,ncorot), ccri(npfft+1, ncorot)
     
      character*3 eletyp
      
      data wri/4*0/
	save wri
c     ;, ini/0/
c      save ini

c     This routine depends on the current finite element type:
      eletyp = styp(isubr)

      call iset(6*nficom, 0, bsy, 1)
      call iset(6*nficom, 0, bssy, 1)
      
      call zset(ncorot*(npfft+1), czero, cvec, 1)

      if(updsym)then
      npp = npft2 + 1
      else
      npp = npfft + 1
      end if

        if(.not.cokpco)then
c       ~~~~~~~~~~~~~~~~~~~        
c       Case of standard theta, phi Fourier expansion
       
            if(eletyp.eq.'M23')then
c           =======================
c           Curl is in rho, theta, phi components

      call dset(ncorot*(npfft+1), 0.d0, cr, 1)
      call dset(6*nficom, 0.d0, b, 1)
      call dset(6*nficom, 0.d0, bs, 1)

C-------------------------------------------------------------------------------
C     NNS: for each rf electric field component, number of nonsingular terms
C     NSI: for each rf electric field component, number of singular terms
C     NFACI: for each component, number of coefficients with a factor sqrt(-1)
C     resulting from toroidal Fourier expansion (written for 0 or 1 coefficient; 
C     modify dimensions of IFACI if new PDE problem with more than one
C     imaginary coefficient per component)
C-------------------------------------------------------------------------------
      nns(1) = 2
      nns(2) = 3
      nns(3) = 2
      nsi(1) = 1
      nsi(2) = 1
      nsi(3) = 4
      nfaci(1) = 1
      nfaci(2) = 1
      nfaci(3) = 0

C-------------------------------------------------------------------------------
C     INS:  localization of nonsingular terms
C     ISI:  localization of singular terms
C     INS, ISI: second index = field component;
C     first index of INS = (column) index of nonsingular term in B.
C     first index of ISI = (column) index of singular term in BS.
C     Contents of INS and ISI: 
C     give position of coefficient in vector of length 8
C       (Erho,   dErho/dtheta,
C        Etheta, dEtheta/dy, dEtheta/dtheta,
C        Ephi,   dEphi/dy,   dEphi/dtheta)
C
C     IFACI: gives position of imaginary factor in the same vector; works in
C            this particular case only.
C-------------------------------------------------------------------------------
      ins(1,1) = 3
      ins(2,1) = 6
      ins(1,2) = 1
      ins(2,2) = 6
      ins(3,2) = 7
      ins(1,3) = 1
      ins(2,3) = 4
      isi(1,1) = 8
      isi(1,2) = 8
      isi(1,3) = 1
      isi(2,3) = 2
      isi(3,3) = 3
      isi(4,3) = 5
      ifaci(1,1) = 3
      ifaci(2,1) = 1

C-------------------------------------------------------------------------------
C     CURL operator:
C     B:    for each field component, list of nonsingular coefficients
C     BS:   for each field component, list of (singular coefficients * y)
C     BSY:  up-down symmetry indices (=±1) of coefficients in B
C     BSSY: up-down symmetry indices (=±1) of coefficients in BS
C-------------------------------------------------------------------------------
      
      bsy(1,1) = uds(1)
      bsy(1,2) = uds(4) * uds(1) * uds(7)
      bssy(1,1) = uds(7)
      
      bsy(2,1) = bsy(1,1)
      bsy(2,2) = uds(6) * uds(1) * uds(7)
      bsy(2,3) = uds(7) * uds(9)
      bssy(2,1) = uds(8) * uds(7) * uds(9)

      bsy(3,1) = uds(11) * uds(7) * uds(9)
      bssy(3,1) = bssy(2,1) * uds(10)
      bssy(3,2) = bssy(1,1)
      bssy(3,3) = bssy(1,1) * uds(10)
      bsy(3,2) = bsy(2,3)
      bssy(3,4) = bssy(2,1)

c     Determine type of up-down symmetry for each coefficient:
      if(ind .eq. 1)then
c     ------------------
c     For curl.curl: take all pairs of coefficients for each component of curl.
      l = 1
c     For displacement current:
      syv(l) = 1
c     Nonsingular matrix elements:
      do ic = 1, nficom
        do j = 1, nns(ic)
          do ka = j, nns(ic)
          l = l + 1
          syv(l) = bsy(ic,j) * bsy(ic,ka)
          end do
          do ka = 1, nsi(ic)
          l = l + 1
          syv(l) = bsy(ic,j) * bssy(ic,ka)
          end do
        end do
      end do
c     Singular matrix elements * y:
      do ic = 1, nficom
        do j = 1, nsi(ic)
          do ka = j, nsi(ic)
          l = l + 1
          syv(l) = bssy(ic,j) * bssy(ic,ka)
          end do
        end do
      end do

      ncoexp = l
	  if(wri(1).eq.0)then
          write(6,*)'curl.curl: number of coefficients = ', ncoexp
	  wri(1) = 1
	  end if
      if(l.gt.ncorot)
     ;write(6,*)'ncorot has a wrong value: declare ncorot>=', l

      else
c     ----
c     For curl (times RNORM):
c     First coefficient is for surface element used in Poynting flux:
      l = 1
      syv(l) = 1
c     Nonsingular terms:
      do ic = 1, nficom
        do j = 1, nns(ic)
        l = l + 1
        syv(l) = bsy(ic,j)
        end do
      end do
c     Singular terms * y:
      do ic = 1, nficom
        do j = 1, nsi(ic)
        l = l + 1
        syv(l) = bssy(ic,j)
        end do
      end do
c     sin and cos of magnetic angle THETA, for coordinate transformations between modes:
      l = l + 1
      syv(l) = 1
      l = l + 1
      syv(l) = 1

      ncoexp = l
	  if(wri(2).eq.0)then
          write(6,*)'curl: number of coefficients = ', ncoexp
	  wri(2) = 1
	  end if
      end if
c     ------

      i = intab
      rho = abscis(i)
      do ipo = 1, npp
c     ###############

	r = eqt(i,ipo,1)
	drthn = eqt(i,ipo,4)
	dzthn = eqt(i,ipo,6)
      ntn = eqt(i,ipo,7)
      ntni = 1.d0 / ntn
	g12n = eqt(i,ipo,8)
      jacn = eqt(i,ipo,9)
      jacni = 1.d0 / jacn
	newmu = eqt(i,ipo,10)
	cn = eqt(i,ipo,11)
	co = eqt(i,ipo,14)
	si = eqt(i,ipo,15)
c      s7 = 1.d0 / eqt(i,ipo,7)
c      s9 = 1.d0 / eqt(i,ipo,9)
      
        if(cyl)then
        ri = 0.d0
c       factor of i=sqrt(-1) dealt with later:
        b(1,1) = - kprn
c       ds = Ntheta:
        ds = ntn * rho
        else
        ri = 1.d0 / r
c       factor of i dealt with later:
        b(1,1) = - n * ri * rnorm
c       ds = R * Ntheta:
        ds = r *  ntn * rho
        end if

      b(1,2) = drthn * ri * ntni * rnorm
      bs(1,1) = ntni

      b(2,1) = - b(1,1)
      b(2,2) = - ri * dzthn * ntni * rnorm
      b(2,3) = - ntn * jacni
      bs(2,1) = g12n * ntni * jacni

      b(3,1) = cn * ntni * jacni * rnorm
      bs(3,1) = - bs(2,1) * newmu
      bs(3,2) = - bs(1,1)
c 16/12/03 bug: following had a / instead of *
      bs(3,3) = bs(1,1) * newmu
      b(3,2) = - b(2,3)
      bs(3,4) = - bs(2,1)

      if(ind.eq.1)then
c     ----------------
c     Normalization for Cray RCFFT2 routine had a factor (R/RA*Jn) / (2*NPFFT):
c     Remains to (* y) all matrix elements and to /y all terms of BS.
c     Nonsingular times singular: *y/y cancel, will apply factor ts:
c      ts = jacn / dfloat(2 * npfft)
c     New: IMSL normalisation:
      ts = jacn / dfloat(npfft)
      if(.not.cyl)ts = ts * r * r0i

C     Nonsingular times nonsingular: multiplic. by y
      t = y * ts
      
      l = 1
c     For displacement current:
      cr(ipo,l) = t

c     Nonsingular matrix elements: 
      do ic = 1, nficom
        do j = 1, nns(ic)
c         All pairs of elements from each row of B, times t:
          do ka = j, nns(ic)
          l = l + 1
          cr(ipo,l) = b(ic,j) * b(ic,ka) * t
          end do
c      ...and  all pairs of elements from each row of B with the same row of BS, times ts:
          do ka = 1, nsi(ic)
          l = l + 1
          cr(ipo,l) = b(ic,j) * bs(ic,ka) * ts
          end do
        end do
      end do

c     Singular matrix elements * y: all pairs of elements from each row of BS, times ts.
c     There is an outstanding division by y, dealt with in routine vacmbg.
      do ic = 1, nficom
        do j = 1, nsi(ic)
          do ka = j, nsi(ic)
          l = l + 1
          cr(ipo,l) = bs(ic,j) * bs(ic,ka) * ts
          end do
        end do
      end do

      else
c     ----
c      Cray RCFFT2 normalization:
c      t = 1.d0 / dfloat(2 * npfft)
c     IMSL normalisation:
      t = 1.d0 / dfloat(npfft)
c     For surface element:
      l = 1
      cr(ipo,l) = ds * t
      do ic = 1, nficom
        do j = 1, nns(ic)
        l = l + 1
        cr(ipo,l) = b(ic,j) * t
        end do
      end do
      do ic = 1, nficom
        do j = 1, nsi(ic)
        l = l + 1
        cr(ipo,l) = bs(ic,j) * t
        end do
      end do
c     cos and sin, for coordinate transformations between modes:
      l = l + 1
      cr(ipo,l) = co * t
      l = l + 1
      cr(ipo,l) = si * t
      
      end if
c     ------
      end do
c     ######

c     Last correction remaining to be done is for 1/y terms
C     IND = 1:
C     (6,16-19,21-23,25,26,28,31).
C     This is taken into account in VACMBG. 
C     In file VACMBG2, the argument SEPAR of VACMBG allows two options:
C     SEPAR=.T.: coefficients of 1/y are stored apart without dividing by y;
C     SEPAR=.F.: division by y is performed and all terms are added.
  
C     Using symmetries:
      if(updsym)then
        do l = 1, ncoexp
        if(syv(l).eq.1)then
          do ipo = 1, npp - 1
          cr(npp+ipo,l) = cr(npp-ipo,l)
          end do
        else
          do ipo = 1, npp - 1
          cr(npp+ipo,l) = - cr(npp-ipo,l)
          end do
        end if
        end do
      end if

c     Call fft,n=2**m case: voir cas sym, antisym;
c     ------------------------------------------- ci-dessous general en reel.
cccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccccccc
c     WARNING: For FFT of REAL functions, the index of the Fourier coefs. 
c     in CRAY and IMSL routines are different:
      do l = 1, ncoexp
cERN      call rcfft2(0, 1, npfft, cr(1,l), workp, cvc(0,l))
	call df2trf(npfft, cr(1,l), cr(1,l), work2) ! -- IMSL --
      end do

c     Storage of IMSL transform in complex array:
c     NB we take the complex conjugate following convention on Fourier transform: we want
c     integral of exp(+i k theta)!
      do l = 1, ncoexp
	cvc(0,l) = dcmplx(cr(1,l), 0.d0)
	  do i = 1, npft2-1
	  cvc(i,l) = dcmplx(cr(2*i,l), -cr(2*i+1,l))
        end do
	cvc(npft2,l) = dcmplx(cr(npfft,l), 0.d0)
      end do

cc      do i = 0,5
cc         print *, cvc(i,1)
cc	end do
cc	print *, 'STOP in Foucog2 Line 364'
cc	stop     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	do l = 1, ncoexp
        do i = 1, npft2
        cvc(-i,l) = dconjg(cvc(i,l))   
        end do
      end do

c     Now CVC(K,L) contains harmonic K><0 of term #L.
      
      return

            else if(eletyp.eq.'HEC')then
C           ============================
C           Curl written in rho, eta, parallel components

C-------------------------------------------------------------------------------
C     NNS: for each rf electric field component, number of nonsingular terms
C     NSI: for each rf electric field component, number of singular terms
C-------------------------------------------------------------------------------
      nns(1) = 3
      nns(2) = 5
      nns(3) = 4
      nsi(1) = 1
      nsi(2) = 1
      nsi(3) = 4

C--------------------------------------------------------------------------------
C     INS:  localization of nonsingular terms
C     ISI:  localization of singular terms
C     INS, ISI: second index = field component;
C     first index of INS = (column) index of nonsingular term in CB.
C     first index of ISI = (column) index of singular term in CBS.
C     Contents of INS and ISI: 
C     give position of coefficient in vector of length 8
C       (Erho, dErho/dtheta,
C        Eeta, dEeta/drho,  dEeta/dtheta,
C        E//,  dE///drho,   dE///dtheta)
C--------------------------------------------------------------------------------
      ins(1,1) = 3
      ins(2,1) = 5
      ins(3,1) = 6
      ins(1,2) = 1
      ins(2,2) = 2
      ins(3,2) = 3
      ins(4,2) = 6
      ins(5,2) = 7
      ins(1,3) = 1
      ins(2,3) = 3
      ins(3,3) = 4
      ins(4,3) = 6
      isi(1,1) = 8
      isi(1,2) = 8
      isi(1,3) = 1
      isi(2,3) = 2
      isi(3,3) = 3
      isi(4,3) = 5

C-------------------------------------------------------------------------------
C     CURL operator:
C     CB:   for each field component, list of nonsingular coefficients (complex)
C     CBS:  for each field component, list of (singular coefficients * y)
C     BSY:  symmetry indices (=±1) of coefficients in CB
C     BSSY: symmetry indices (=±1) of coefficients in CBS
C     For complex entries, the parity of the real part is used;
C     the parity of the imaginary part is always opposite
C-------------------------------------------------------------------------------
      
      bsy(1,1) = -1
      bsy(1,2) =  1
      bsy(1,3) = -1
      bssy(1,1) = 1
      
      bsy(2,1) = -1
      bsy(2,2) =  1
      bsy(2,3) =  1
      bsy(2,4) =  1
      bsy(2,5) =  1
      bssy(2,1) = -1

      bsy(3,1) =  -1
      bsy(3,2) =   1
      bsy(3,3) =   1
      bsy(3,4) =   1
      bssy(3,1) = -1
      bssy(3,2) =  1
      bssy(3,3) =  1
      bssy(3,4) = -1

            else if(eletyp.eq.'CAR')then
C           ============================
C           Curl written in R, Y, phi components

C-------------------------------------------------------------------------------
C     NNS: for each rf electric field component, number of nonsingular terms
C     NSI: for each rf electric field component, number of singular terms.
C-------------------------------------------------------------------------------
      nns(1) = 2
      nsi(1) = 1
      nns(2) = 2
      nsi(2) = 1
      nns(3) = 2
      nsi(3) = 2

C--------------------------------------------------------------------------------
C     INS:  localization of nonsingular terms
C     ISI:  localization of singular terms
C     INS, ISI: second index = field component;
C     first index of INS = (column) index of nonsingular term in CB.
C     first index of ISI = (column) index of singular term in CBS.
C     Contents of INS and ISI: 
C     give position of coefficient in vector of length 9
C       (ER, dER/drho, dER/dtheta,
C        EY, dEY/drho, dEY/dtheta,
C        Ephi,  dEphi/drho,   dEphi/dtheta)
C--------------------------------------------------------------------------------
      ins(1,1) = 4
      ins(2,1) = 8
      isi(1,1) = 9
      ins(1,2) = 1
      ins(2,2) = 8
      isi(1,2) = 9
      ins(1,3) = 2
      isi(1,3) = 3
      ins(2,3) = 5
      isi(2,3) = 6

C-------------------------------------------------------------------------------
C     CURL operator:
C     CB:   for each field component, list of nonsingular coefficients (complex)
C     CBS:  for each field component, list of (singular coefficients * y)
C     BSY:  symmetry indices (=±1) of coefficients in CB
C     BSSY:  symmetry indices (=±1) of coefficients in CBS
C     For complex entries, the parity of the real part is used;
C     the parity of the imaginary part is always opposite
C-------------------------------------------------------------------------------
      
c      bsy(1,1) =  1
      bsy(1,1) = -1
c     (NB: remember convention: parity of real part, here zero)
      bsy(1,2) = -1
      bssy(1,1) =  1
      
c      bsy(2,1) =  1
      bsy(2,1) = -1
c     (NB: remember convention: parity of real part, here zero)
      bsy(2,2) =  1
      bssy(2,1) = -1

      bsy(3,1) = -1
      bssy(3,1) =  1
      bsy(3,2) =  1
      bssy(3,2) = -1

            end if
c           ======

        else
c       ~~~~      
c       Constant k// coordinates:

            if(eletyp.eq.'M23')then
c           =======================
c      ncoexp = 68

      nns(1) = 2
cERN      nns(2) = 3
cERN      nns(3) = 2
      nns(2) = 4		! theta: cERN 3 -> 4
      nns(3) = 4		! phi  : cERN 2 -> 4

      nsi(1) = 1
      nsi(2) = 1
      nsi(3) = 4

      ins(1,1) = 3
      ins(2,1) = 6
      ins(1,2) = 1
      ins(2,2) = 6
      ins(3,2) = 7
	ins(4,2) = 8	! cERN: a(theta,8) -> cb(2,4)

      ins(1,3) = 1
      ins(2,3) = 3
	ins(3,3) = 4	! cERN: a(phi,4) -> cb(3,3)
	ins(4,3) = 5	! cERN: a(phi,5) -> cb(3,4)


      isi(1,1) = 8
      isi(1,2) = 8
      isi(1,3) = 1
      isi(2,3) = 2
      isi(3,3) = 3
      isi(4,3) = 5
      
      bsy(1,1) = -1
c     (NB: remember convention: parity of real part, here zero)
      bsy(1,2) = -1
      bssy(1,1) = 1
      
      bsy(2,1) = bsy(1,1)
      bsy(2,2) = 1
      bsy(2,3) = 1
	bsy(2,4) = -1   ! cERN
      bssy(2,1) = -1

      bsy(3,1) = -1
      bssy(3,1) = -1
      bssy(3,2) = 1
      bssy(3,3) = 1
      bsy(3,2) = 1
	bsy(3,3) =  1	! cERN
	bsy(3,4) = -1	! cERN

c     (NB: remember convention: parity of real part, here zero)
      bssy(3,4) = -1

            else if(eletyp.eq.'HEC')then
C           ============================
C           Curl written in rho, eta, parallel components

c      ncoexp = 68

      nns(1) = 3
c      nns(2) = 5
c      nns(3) = 4
cPL23Nov04 these two lines were commented out:
      nns(2) = 6		! eta : cERN 5 -> 6
      nns(3) = 5		! //  : cERN 4 -> 5

      nsi(1) = 1
      nsi(2) = 1
      nsi(3) = 4

      ins(1,1) = 3
      ins(2,1) = 5
      ins(3,1) = 6
      ins(1,2) = 1
      ins(2,2) = 2
      ins(3,2) = 3
      ins(4,2) = 6
      ins(5,2) = 7
      ins(6,2) = 8	! cERN: a(eta,8) -> cb(2,6)
      ins(1,3) = 1
      ins(2,3) = 3
      ins(3,3) = 4
      ins(4,3) = 5	! cERN: a(//,5) -> cb(3,4)
c      ins(4,3) = 6
      ins(5,3) = 6	! a(//,6) -> cb(3,5) : OK

      isi(1,1) = 8
      isi(1,2) = 8
      isi(1,3) = 1
      isi(2,3) = 2
      isi(3,3) = 3
      isi(4,3) = 5

      bsy(1,1) = -1
      bsy(1,2) =  1
      bsy(1,3) = -1
      bssy(1,1) = 1
      
      bsy(2,1) = -1
      bsy(2,2) =  1
      bsy(2,3) =  1
      bsy(2,4) =  1
      bsy(2,5) =  1
	bsy(2,6) = -1	! cERN
      bssy(2,1) = -1

      bsy(3,1) =  -1
      bsy(3,2) =   1
      bsy(3,3) =   1
c	bsy(3,4) =   1
      bsy(3,4) = -1	! cERN
	bsy(3,5) =  1	! (3,4) -> (3,5)

      bssy(3,1) = -1
      bssy(3,2) =  1
      bssy(3,3) =  1
      bssy(3,4) = -1

            else if(eletyp.eq.'CAR')then
C           ============================
C           Curl written in R, Y, phi components
      print *, 'FOUCOG2: element type CAR in constant k// coord. to write'
	stop
            end if
c           ======

        end if
c       ~~~~~~       

      if(ind.eq.1)then
c     ----------------
c     Symmetries (only used when updsym=.T.):
      l = 1
c     For displacement current:
      syv(l) = 1
c     Nonsingular matrix elements:
      do ic = 1, nficom
        do j = 1, nns(ic)
          do ka = j, nns(ic)
          l = l + 1
          syv(l) = bsy(ic,j) * bsy(ic,ka)
          end do
          do ka = 1, nsi(ic)
          l = l + 1
          syv(l) = bsy(ic,j) * bssy(ic,ka)
          end do
        end do
      end do
C     Singular matrix elements * y:
      do ic = 1, nficom
        do j = 1, nsi(ic)
          do ka = j, nsi(ic)
          l = l + 1
          syv(l) = bssy(ic,j) * bssy(ic,ka)
          end do
        end do
      end do

c     Number of coefficients:
      ncoexp = l
	  if(wri(3).eq.0)then
          write(6,*)'curl.curl: number of coefficients = ', ncoexp
	  wri(3) = 1
	  end if
c       if(l.ne.ncoexp)then
c       write(6,*)'strange value found for number of curl coefficients: ', l
c       write(6,*)'expected ',ncoexp,' coefficients!'
c       write(6,*)'eletyp= ',eletyp,'; cokpco= ',cokpco
c       write(6,*)'nb: ncorot= ',ncorot,' must be >= above number'
c       stop
c       end if
      else
c     ----
      l = 1
      syv(l) = 1
c     Nonsingular terms:
      do ic = 1, nficom
        do j = 1, nns(ic)
        l = l + 1
        syv(l) = bsy(ic,j)
        end do
      end do
c     Singular terms * y:
      do ic = 1, nficom
        do j = 1, nsi(ic)
        l = l + 1
        syv(l) = bssy(ic,j)
        end do
      end do
c     sin and cos, for coordinate transormations between modes:
      l = l + 1
      syv(l) = 1
      l = l + 1
      syv(l) = 1
c@ CAR: some more needed!
      ncoexp = l
	  if(wri(4).eq.0)then
          write(6,*)'curl: number of coefficients = ', ncoexp
	  wri(4) = 1
	  end if
      end if
c     ------

c     Now build tables of coefficients to transform:
      call zset(ncorot*(npfft+1), czero, ccr, 1)
      call zset(6*nficom, czero, cb, 1)
      call zset(6*nficom, czero, cbs, 1)

        if(.not.cokpco)then
c       ~~~~~~~~~~~~~~~~~~~        
c       Case of theta, phi Fourier expansion
       
            if(eletyp.eq.'HEC')then
c           =======================
      i = intab
      rho = abscis(i)
      do ipo = 1, npp
c     ###############
      
	r = eqt(i,ipo,1)
	drthn = eqt(i,ipo,4)
	dzthn = eqt(i,ipo,6)
      ntn = eqt(i,ipo,7)
      ntni = 1.d0 / ntn
	g12n = eqt(i,ipo,8)
      jacn = eqt(i,ipo,9)
      jacni = 1.d0 / jacn
	newmu = eqt(i,ipo,10)
	cn = eqt(i,ipo,11)
	sior = eqt(i,ipo,12)
	tgth = eqt(i,ipo,13)
	co = eqt(i,ipo,14)
	si = eqt(i,ipo,15)
	si1 = eqt(i,ipo,16)
	si1tn = eqt(i,ipo,17)
      
        if(cyl)then
        crnor = 0.d0
        srnor = 0.d0
        ncrnor = kprn * co
        nsrnor = kprn * si
        rnllnr = 0.d0
        ds = ntn * rho
        else
        ri = 1.d0 / r
        crnor = ri * co * rnorm
        srnor = ri * si * rnorm
        ncrnor = n * crnor
        nsrnor = n * srnor
        rnllnr = - rnorm * ri * ntni * dzthn
        ds = r *  ntn * rho
        end if
        
      aux1 = si * rnllnr
     ;     + rnorm * sior * ntni * newmu
c     Rn*(sLc-cLs):
      aux2 = - rnorm * (g12n * si1tn * ntni - ntn * si1) / (jacn * co)
      
      cb(1,1) = dcmplx(-(si1tn * rnorm + srnor * drthn) * ntni, - ncrnor)
      cb(1,2) = - rnorm * sior  * ntni
      cb(1,3) = 
     ;dcmplx((-tgth*si1tn * rnorm + crnor * drthn) * ntni, - nsrnor )
      cbs(1,1) = co * ntni
      
      cb(2,1) = 
     ;dcmplx(- rnorm * (si*cn - g12n*sior*newmu) * ntni * jacni, ncrnor)
      cb(2,2) = - cb(1,2)
      cb(2,3) = - co * aux1 + aux2
      cb(2,4) = rnllnr - si * aux1
      cb(2,5) = - ntn  * jacni
      cbs(2,1) = g12n * ntni * jacni

      cb(3,1) = dcmplx(co * rnorm * cn * ntni * jacni, nsrnor)
      cb(3,2) = - si * aux1
      cb(3,3) = - cb(2,5)
      cb(3,4) = co * aux1 + aux2
      cbs(3,1) = - cbs(2,1) * co * newmu
      cbs(3,2) = - cbs(1,1)
      cbs(3,3) = ntni * newmu
      cbs(3,4) = - cbs(2,1)

      if(ind.eq.1)then
c     ----------------
C     Normalization for CFFT2 routine: fft over Theta
C     factor (R/RA*Jn) / (NPFFT); difference fac. 2 w.r. to RCFFT2!: 
C     Remains to (* y) all matrix elements and to /y all terms of cBS.
C     Singular: *y/y, do nothing:
      ts = jacn / dfloat(npfft)
      if(.not.cyl)ts = ts * r * r0i
c     Nonsingular: *y
      t = y * ts

      l = 1
c     For displacement current:
      ccr(ipo,l) = t

c     Nonsingular matrix elements: all pairs of elements from each row of cb
c     and  all pairs of elements from each row of cb with the same row of cbs.
c     NB: row indices correspond to test function, hence complex conjugate
c     the corresponding entries; column indices correspond to rf field.
      do ic = 1, nficom
c       Row sweep:
        do j = 1, nns(ic)
c         Column sweep:
          do ka = j, nns(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cb(ic,j)) * cb(ic,ka) * t
          end do
          do ka = 1, nsi(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cb(ic,j)) * cbs(ic,ka) * ts
          end do
        end do
      end do

c     Singular matrix elements * y: all pairs of elements from each row of BS.
      do ic = 1, nficom
c       Row sweep:
        do j = 1, nsi(ic)
c         Column sweep:
          do ka = j, nsi(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cbs(ic,j)) * cbs(ic,ka) * ts
          end do
        end do
      end do
      
      else
c     ----
      t = 1.d0 / dfloat(npfft)
c     For surface element:
      l = 1
      ccr(ipo,l) = ds * t
      
      do ic = 1, nficom
        do j = 1, nns(ic)
        l = l + 1
        ccr(ipo,l) = cb(ic,j) * t
        end do
      end do
      do ic = 1, nficom
        do j = 1, nsi(ic)
        l = l + 1
        ccr(ipo,l) = cbs(ic,j) * t
        end do
      end do
c     cos and sin, for coordinate transormations between modes:
      l = l + 1
      ccr(ipo,l) = co * t
      l = l + 1
      ccr(ipo,l) = si * t
      end if
c     ------

      end do
c     ######

            else if(eletyp.eq.'CAR')then
c           ============================
      print *, 'CAR element FOUCOG2: at work line 892'
	stop
      i = intab
      rho = abscis(i)
      do ipo = 1, npp
c     ###############
      
	r = eqt(i,ipo,1)
	drrho = eqt(i,ipo,3)
	drthn = eqt(i,ipo,4)
	dzrho = eqt(i,ipo,5)
	dzthn = eqt(i,ipo,6)
c      ntn = eqt(i,ipo,7)
c      ntni = 1.d0 / ntn
c	g12n = eqt(i,ipo,8)
      jacn = eqt(i,ipo,9)
      jacni = 1.d0 / jacn
c	newmu = eqt(i,ipo,10)
c	cn = eqt(i,ipo,11)
c	sior = eqt(i,ipo,12)
c	tgth = eqt(i,ipo,13)
	co = eqt(i,ipo,14)
	si = eqt(i,ipo,15)
c	si1 = eqt(i,ipo,16)
c	si1tn = eqt(i,ipo,17)
      
        if(cyl)then
        rnor = 0.d0
        rnllnr = 0.d0
        ds = ntn * rho
        else
        ri = 1.d0 / r
        rnor = ri * rnorm
        rnllnr = - rnorm * ri * ntni * dzthn
        ds = r *  ntn * rho
        end if
        
      cb(1,1) = dcmplx(0.d0, - n * rnor)
      cb(1,2) = - drthn * jacni
      cbs(1,1) = drrho * jacni
      
      cb(2,1) = - cb(1,1)
      cb(2,2) = - dzthn * jacni
      cbs(2,1) = dzrho * jacni

      cb(3,1) = - cb(1,2)
      cbs(3,1) = - cbs(1,1)
      cb(3,2) = - cb(2,2)
      cbs(3,2) = - cbs(2,1)

      if(ind.eq.1)then
c     ----------------
C     Normalization for CFFT2 routine: fft over Theta
C     factor (R/RA*Jn) / (NPFFT); difference fac. 2 w.r. to RCFFT2!: 
      ts = jacn / dfloat(npfft)
      if(.not.cyl)ts = ts * r * r0i
c     Nonsingular: *y
      t = y * ts

      l = 1
c     For displacement current:
      ccr(ipo,l) = t

c     Nonsingular matrix elements: all pairs of elements from each row of cb
c     and  all pairs of elements from each row of cb with the same row of cbs.
c     NB: row indices correspond to test function, hence complex conjugate
c     the corresponding entries; column indices correspond to rf field.
      do ic = 1, nficom
c       Row sweep:
        do j = 1, nns(ic)
c         Column sweep:
          do ka = j, nns(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cb(ic,j)) * cb(ic,ka) * t
          end do
          do ka = 1, nsi(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cb(ic,j)) * cbs(ic,ka) * ts
          end do
        end do
      end do

c     Singular matrix elements * y: all pairs of elements from each row of BS.
      do ic = 1, nficom
c       Row sweep:
        do j = 1, nsi(ic)
c         Column sweep:
          do ka = j, nsi(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cbs(ic,j)) * cbs(ic,ka) * ts
          end do
        end do
      end do
      
      else
c     ----
      t = 1.d0 / dfloat(npfft)
c     For surface element:
      l = 1
      ccr(ipo,l) = ds * t
      
      do ic = 1, nficom
        do j = 1, nns(ic)
        l = l + 1
        ccr(ipo,l) = cb(ic,j) * t
        end do
      end do
      do ic = 1, nficom
        do j = 1, nsi(ic)
        l = l + 1
        ccr(ipo,l) = cbs(ic,j) * t
        end do
      end do
c     cos and sin, for coordinate transormations between modes:
      l = l + 1
      ccr(ipo,l) = co * t
      l = l + 1
      ccr(ipo,l) = si * t
c@ CAR: and more!
      end if
c     ------

      end do
c     ######
            end if
c           ======

        else
c       ~~~~       
c       Case of Thbar, Phibar Fourier expansion

      i = intab
      rho = abscis(i)
      ax = .not.crown .and. i.eq.1
        if(ax)then
        ii = 2
        else
        ii = i
        end if

            if(eletyp.eq.'M23')then
c           =======================
      do ipo = 1, npp
c     ###############
	r = eqt(i,ipo,1)
	drthn = eqt(i,ipo,4)
	dzthn = eqt(i,ipo,6)
	drrho = eqt(i,ipo,3)  ! ERN
      ntn = eqt(i,ipo,7)
      ntni = 1.d0 / ntn
	g12n = eqt(i,ipo,8)
      jacn = eqt(i,ipo,9)
      jacni = 1.d0 / jacn
	newmu = eqt(i,ipo,10)
	cn = eqt(i,ipo,11)
	sior = eqt(i,ipo,12)
	tgth = eqt(i,ipo,13)
	co = eqt(i,ipo,14)
	si = eqt(i,ipo,15)
	si1 = eqt(i,ipo,16)
	si1tn = eqt(i,ipo,17)
c      ntni = 1.d0 / eqt(i,ipo,7)
c      jacni = 1.d0 / eqt(i,ipo,9)
c      s7 = 1.d0 / eqt(1,ipo,7)
c      s9 = 1.d0 / eqt(1,ipo,9)
      
        if(cyl)then
        ri = 0.d0
        cb(1,1) = dcmplx(0.d0, - kprn)
        ds = ntn * rho
        else
        ri = 1.d0 / r
        cb(1,1) = dcmplx(0.d0, - n * rnorm * ri)
        ds = r * ntn * rho
        end if
c     -------------------> BUG here: *abscis instead of /abscis
cERN      ontpbt = ckt(ii,ipo,2) * ntni * abscni(ii)
		ontpbt = ckt(ii,ipo,2) * ntni / abscni(ii)
c older, supposedly equivalent but division by 0 at axis:
cc       ontpbt = rnorm*(qfactor(i)*hachi(i)-co*ri) / si
      
      cb(1,2) = dcmplx(ri * rnorm * drthn * ntni, n * ontpbt)  ! OK

cc      cbs(1,1) = hachi(i) / sior
      cbs(1,1) = hachi(i) / sior ! ?????
	      
      cb(2,1) = dconjg(cb(1,1))

      nrnlpb = n * (ontpbt*g12n   - rnorm*ntn*ckt(i,ipo,7) ) * jacni
cERN	--------------------------------------------------------
      cb(2,2) = dcmplx(- ri * rnorm * dzthn * ntni, nrnlpb)

cERN	if(abscis(i) .eq. 0)abscis(i) = 1e-4
cc	cb(2,2) = 1.0d0* 
cc     ;   dcmplx(-ri*ntn*jacni*drrho*rnorm   + ri*g12n*jacni*ntni*drthn*rnorm, 
cc     ;           nrnlpb)
c	----------------------------------------------------
      cb(2,3) = - ntn * jacni
cERN	cb(2,3) = - ntn * jacni * rnorm

c      cb(2,4) = cb(2,3) * rnorm * dzthn
c 16/12/03
      cb(2,4) = cb(2,3) * rnorm * ckt(i,ipo,6)

      cbs(2,1) = cbs(1,1) * g12n * jacni

      cb(3,1) = dcmplx(rnorm * cn * ntni * jacni, - n * ontpbt)
c		CAGADA
c     cb(3,1) = dcmplx( rnorm * cn * ntni * jacni, 0)


      cb(3,2) = dcmplx(0.d0, -nrnlpb)


      cb(3,3) = - cb(2,3)
      cb(3,4) = - cb(2,4)
      cbs(3,2) = - cbs(1,1)
      cbs(3,3) = ntni * newmu
      cbs(3,1) = - cbs(3,3) * g12n * jacni
      cbs(3,4) = - cbs(2,1)

      if(ind.eq.1)then
c     ----------------
C     Normalization for CFFT2 routine: fft over Theta
C     factor (R/RA*Jn) / (NPFFT); difference fac. 2 w.r. to RCFFT2!: 
C     Remains to (* y) all matrix elements and to /y all terms of cBS.
C     Singular: *y/y, do nothing:
      ts = jacn / dfloat(npfft)
      if(.not.cyl)ts = ts * r * r0i
C     Nonsingular: *y
      t = y * ts

      l = 1
c     For displacement current:
      ccr(ipo,l) = t

c     Nonsingular matrix elements: all pairs of elements from each row of cB
c     and  all pairs of elements from each row of cB with the same row of cBS.
      do ic = 1, nficom
        do j = 1, nns(ic)
          do ka = j, nns(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cb(ic,j)) * cb(ic,ka) * t
          end do
          do ka = 1, nsi(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cb(ic,j)) * cbs(ic,ka) * ts
          end do
        end do
      end do

c     Singular matrix elements * y: all pairs of elements from each row of BS.
      do ic = 1, nficom
        do j = 1, nsi(ic)
          do ka = j, nsi(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cbs(ic,j)) * cbs(ic,ka) * ts
          end do
        end do
      end do
      
      else
c     ----
      t = 1.d0 / dfloat(npfft)
c     For surface element:
      l = 1
      ccr(ipo,l) = ds * t

      do ic = 1, nficom
        do j = 1, nns(ic)
        l = l + 1
        ccr(ipo,l) = cb(ic,j) * t
        end do
      end do
      do ic = 1, nficom
        do j = 1, nsi(ic)
        l = l + 1
        ccr(ipo,l) = cbs(ic,j) * t
        end do
      end do
c     cos and sin, for coordinate transormations between modes:
      l = l + 1
      ccr(ipo,l) = co * t
      l = l + 1
      ccr(ipo,l) = si * t
      end if
c     ------

      end do
c     ######
            else if(eletyp.eq.'HEC')then
c           ============================
      do ipo = 1, npp
c     ###############
	r = eqt(i,ipo,1)
	drthn = eqt(i,ipo,4)
	dzthn = eqt(i,ipo,6)
      ntn = eqt(i,ipo,7)
      ntni = 1.d0 / ntn
	g12n = eqt(i,ipo,8)
      jacn = eqt(i,ipo,9)
      jacni = 1.d0 / jacn
	newmu = eqt(i,ipo,10)
	cn = eqt(i,ipo,11)
	sior = eqt(i,ipo,12)
	tgth = eqt(i,ipo,13)
	co = eqt(i,ipo,14)
	si = eqt(i,ipo,15)
	si1 = eqt(i,ipo,16)
	si1tn = eqt(i,ipo,17)

c      ntni = 1.d0 / eqt(i,ipo,7)
c      jacni = 1.d0 / eqt(i,ipo,9)
c	co = eqt(ii,ipo,14)
c      s7 = 1.d0 / eqt(i,ipo,7)
c      s9 = 1.d0 / eqt(i,ipo,9)
      
cERN	q = eqta1d(i,4)
	q = qfactor(i)

        if(cyl)then
        ri = 0.d0
cERN       nb factor 1/r properly omitted in eqta1d( ,4):
c       nb factor 1/r properly omitted in qfactor:
cPL See why singular term contributes to cb: !
        aux4 =  kprn * (q*hachi(ii)*co - 1.d0) / dmax1(si, 1.d-10)
        ds = ntn * rho
        else
        ri = 1.d0 / r
c       aux3=(qhc-1/r)/s
cPL See why singular term contributes to cb: !
        aux3 = (q*hachi(ii)*co-1.d0/r) / dmax1(si, 1.d-10)
        aux4 = n * rnorm * aux3
        ds = r * ntn * rho
        end if
        
      cb(1,1) = dcmplx(
     ;  - rnorm * ntni * (ri * si * drthn + si1tn), - n * rnorm * q * hachi(i))
      cb(1,2) = - hachi(i) * rnorm
c     aux3=(qhc-1/r)/s
      ontpbt = ckt(ii,ipo,2) * ntni / abscno(ii)
cc  older:
cc      ontpbt = rnorm*(eqta1d(2,4)*hachi(2)-eqt(2,ipo,14)/eqt(2,ipo,1))
cc     ;         / eqt(2,ipo,15)
      nrnlpb = n * (ontpbt*g12n - rnorm*ntn*ckt(i,ipo,7)) * jacni
      cb(1,3) = dcmplx(rnorm * ntni * (ri * co * drthn - tgth * si1tn), aux4)

      cbs(1,1) = hachi(i) * co / sior
      
      cb(2,1) = dcmplx(
     ;  rnorm * ntni * jacni * (- si * cn + sior * g12n * newmu)
     ;, n * rnorm * q * hachi(i))

      cb(2,2) = - cb(1,2)

      rnlsoc = rnorm * jacni / co * (g12n * si1tn * ntni - ntn * si1) 

      aux1 = rnorm * ntni * (si * dzthn * ri - sior * newmu)

      cb(2,3) = - rnlsoc + co * aux1
      cb(2,4) = dcmplx(
     ;aux1 * si - rnorm * dzthn * ri * ntni, nrnlpb) 
      cb(2,5) = - ntn * jacni
      cb(2,6) = - rnorm * ntn * ckt(i,ipo,6) * jacni

      cbs(2,1) = g12n * hachi(i) * jacni / sior

      cb(3,1) = dcmplx(co * rnorm * cn * ntni * jacni, - aux4)
      cb(3,2) =  dcmplx(aux1 * si, - nrnlpb)
      cb(3,3) = - cb(2,5)
      cb(3,4) = - cb(2,6)
      cb(3,5) = - rnlsoc - co * aux1

      cbs(3,2) = - cbs(1,1)
      cbs(3,3) = ntni * newmu
      cbs(3,1) = - cbs(3,3) * co * g12n * jacni
      cbs(3,4) = - cbs(2,1)

      if(ind.eq.1)then
c     ----------------
c     Normalization for IMSL or Cray CFFT2 routines: fft over Theta
c     factor (R/RA*Jn) / (NPFFT) 
c     Remains to (* y) all matrix elements and to /y all terms of cBS.
c     Singular terms: *y/y, so do nothing:
      ts = jacn / dfloat(npfft)
      if(.not.cyl)ts = ts * r * r0i
c     Nonsingular: *y
      t = y * ts

      l = 1
c     For displacement current:
      ccr(ipo,l) = t

c     Nonsingular matrix elements: all pairs of elements from each row of cB
c     and  all pairs of elements from each row of cB with the same row of cBS.
      do ic = 1, nficom
        do j = 1, nns(ic)
          do ka = j, nns(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cb(ic,j)) * cb(ic,ka) * t
          end do
          do ka = 1, nsi(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cb(ic,j)) * cbs(ic,ka) * ts
          end do
        end do
      end do

c     Singular matrix elements * y: all pairs of elements from each row of BS.
      do ic = 1, nficom
        do j = 1, nsi(ic)
          do ka = j, nsi(ic)
          l = l + 1
          ccr(ipo,l) = dconjg(cbs(ic,j)) * cbs(ic,ka) * ts
          end do
        end do
      end do
      
      else
c     ----
      t = 1.d0 / dfloat(npfft)
c     For surface element:
      l = 1
      ccr(ipo,l) = ds * t

      do ic = 1, nficom
        do j = 1, nns(ic)
        l = l + 1
        ccr(ipo,l) = cb(ic,j) * t
        end do
      end do
      do ic = 1, nficom
        do j = 1, nsi(ic)
        l = l + 1
        ccr(ipo,l) = cbs(ic,j) * t
        end do
      end do
c     cos and sin, for coordinate transormations between modes:
      l = l + 1
      ccr(ipo,l) = co * t
      l = l + 1
      ccr(ipo,l) = si * t
      end if
c     ------

      end do
c     ######
            else if(eletyp.eq.'CAR')then
c           ============================
            print *, 'FOUCOG2: CAR element in const. k// to write, line 1298'
		stop 
            end if
c           ======
c     Interpolate to uniform Thbar grid:
      do l = 1, ncoexp
      call interp2c(ckt(i,1,3), nabplo, ccr(1,l), 1, npp, polang, ccri(1,l)
     ;, npp)
         do ipo = 1, npp
         ccr(ipo,l) = ccri(ipo,l)
         end do
      end do
      end if
c     ~~~~~~

c     Last correction remaining to be done is for 1/y terms.
c     This will be taken into account in VACMBG. 
c     In file VACMBG2, the argument SEPAR of VACMBG allows two options:
c     SEPAR=.T.: coefficients of 1/y are stored apart without dividing by y;
c     SEPAR=.F.: division by y is performed and all terms are added.
  
c     Using symmetries:
      if(updsym)then
        do l = 1, ncoexp
        if(syv(l).eq.1)then
          do ipo = 1, npp - 1
          ccr(npp+ipo,l) = dconjg(ccr(npp-ipo,l))
          end do
        else
          do ipo = 1, npp - 1
          ccr(npp+ipo,l) = - dconjg(ccr(npp-ipo,l))
          end do
        end if
        end do
      end if

c     Call fft,n=2**m case: voir cas sym, antisym;
c     --------------------------------------------ci-dessous general en cx.

cccccccccccccccccccccccc ERNESTO ccccccccccccccccccccccccccccc
      do l = 1, ncoexp
cERN      call cfft2(0, 1, npfft, ccr(1,l), cworkp, cvc(0,l))
	call df2tcb(npfft, ccr(1,l), cvc(0,l), cwork2, cpy) ! -- IMSL --
      end do
c      print *,  cvc(1:10,1)
c	print *, 'Stop at FOUCOG2'
c	stop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      n2 = npfft / 2
      do l = 1, ncoexp
        do i = 1, npft2
        cvc(-i,l) = cvc(npfft-i,l)
        end do
      end do

c     Now CVC(K,L) contains harmonic K><0 of term #L.
c     for k=-npft2 to npft2
      
c     Savings and only keep NCROT+1 terms: could sample abscissae      

      return
      end
