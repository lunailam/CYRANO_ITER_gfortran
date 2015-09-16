      subroutine foucoj(rtp)

ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)
      implicit none
      
      logical rtp

c     Computes Fourier coeffs. used in jump conditions at region boundaries
c     using IMSL fft (complex or/and real)
c     (formerly using Cray cfft2 or rcfft2)
c     Radial table index INTAB must be provided through COMSWE.
c     D-SHAPED or general cross section.
c     rtp: for HEC elements, choice of coordinate system to write boundary
c     conditions.
c     This version assumes all geom. coefficients radially continuous.

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

      integer 
     ;  ini, i, ipo, l, syv(ncorot)
     ;, bsy(nficom,6), bssy(nficom,6)
     ;, ncoexp

      double precision t, cr(npfft+1,ncorot)
     ;, ncrnor, nsrnor, rnllnr, aux1, aux2, aux4, ontpbt
     ;, jon
     ;, ri, ntni, sior, tr

      complex*16
     ;  ccr(npfft+1,ncorot), ccri(npfft+1, ncorot)
     ;, cdeid
     ;, facju
     
      external cdeid
      
      character*3 eletyp
      
      data ini/0/

      save ini, facju

c     This routine depends on the current finite element type; hence to be 
c     called on both sides of interface:
      eletyp = styp(isubr)

      call iset(6*nficom, 0, bsy, 1)
      call iset(6*nficom, 0, bssy, 1)
      
      call zset(ncorot*(npfft+1), czero, cvec, 1)

      if(ini.eq.0)then
c     -1/fac of old notes:
      facju = dcmplx(0.d0, - omegag * mu0 * rnorm)
      ini = 1
      end if

      yinv = abscni(intab)
      
        if(.not.cokpco)then
c       ~~~~~~~~~~~~~~~~~~~        
c       Case of theta, phi Fourier expansion
       
          if(eletyp.eq.'M23' .or. (eletyp.eq.'HEC'.and. rtp))then
c         =======================================================
c         Curl is in rho, theta, phi components

      call dset(ncorot*(npfft+1), 0.d0, cr, 1)

      ncoexp = 5
      syv(1) = 1
      syv(2) = -1
      syv(3) = 1
      syv(4) = 1
      syv(5) = 1
        if(rtp)then
        ncoexp = 7
        syv(6) = 1
        syv(7) = 1
        end if
      
      i = intab
c     Normalization for RCFFT2 routine: 
c      t = 1.d0 / dfloat(2 * npfft)
c     Normalization for IMSL: 
      t = 1.d0 / dfloat(npfft)
        if(cyl)then
          do ipo = 1, npp
	    ntn = eqt(i,ipo,7)
	    g12n = eqt(i,ipo,8)
	    jacn = eqt(i,ipo,9)
	    newmu = eqt(i,ipo,10)
	    cn = eqt(i,ipo,11)

c         Factor of i delayed:
          cr(ipo,1) = t * (- jacn * kprn) / ntn
          cr(ipo,2) = t * (cn * rnorm - g12n * newmu * yinv) / (ntn*ntn)
          cr(ipo,3) = - t * jacn * yinv / (ntn*ntn)

c         For current spectrum, resp. pol. and tor. (convol. by Morse to do):
ccc not in cyl:          cr(ipo,4) =  t / eqt(i,ipo,1)
c          cr(ipo,5) =  t * 0.d0
          end do
        else
          do ipo = 1, npp
	    r = eqt(i,ipo,1)
	    ri = 1.d0 / r
	    ntn = eqt(i,ipo,7)
	    ntni = 1.d0 / ntn
	    g12n = eqt(i,ipo,8)
	    jacn = eqt(i,ipo,9)
	    newmu = eqt(i,ipo,10)
	    cn = eqt(i,ipo,11)
	    co = eqt(i,ipo,14)
	    si = eqt(i,ipo,15)

c         Factor of i delayed:
          cr(ipo,1) = t * (- jacn * n * rnorm) * ri * ntni
          cr(ipo,2) = t * (cn * rnorm - g12n * newmu * yinv) * ntni * ntni
          cr(ipo,3) = - t * jacn * yinv * ntni * ntni

c         For current spectrum, resp. pol. and tor. (convol. by Morse pulse and
c         multiplication by facju to do):
          cr(ipo,4) = t * ri
          cr(ipo,5) = 0.d0

          cr(ipo,6) = t * co
          cr(ipo,7) = t * si
          end do
        end if
  
C     Using symmetries:
      if(updsym)then
        do l = 1, ncoexp
          if(syv(l).eq.1)then
            do ipo = 1, npp-1
            cr(npp+ipo,l) = cr(npp-ipo,l)
            end do
          else
            do ipo = 1, npp-1
            cr(npp+ipo,l) = - cr(npp-ipo,l)
            end do
          end if
        end do
      end if

c     call fft, n=2**m case: voir cas sym, antisym;
c     --------------------------------------------- ci-dessous general en reel.

cccccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccccc 
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

c      do i = 0,5
c         print *, cvc(i,1)
c	end do
c	print *, 'STOP in Foucoj Line 171'
c	stop     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l = 1, ncoexp
        do i = 1, npft2
        cvc(-i,l) = dconjg(cvc(i,l))   
        end do
      end do

c       Delayed factor of i on first tf, imaginary factor facju on currents:
        do i = - npft2, npft2
        cvc(i,1) = ci * cvc(i,1)
        cvc(i,4) = facju * cvc(i,4)
        cvc(i,5) = facju * cvc(i,5)
        end do
c     Now CVC(K,L) contains harmonic K><0 of term #L.
      
      return

            else if(eletyp.eq.'HEC')then
c           ============================
c             Curl written in rho, eta, parallel components
              ncoexp = 6
              syv(1) = -1
              syv(2) =  1
              syv(3) = -1
              syv(4) =  1
              syv(5) =  1
              syv(6) =  1

            else if(eletyp.eq.'CAR')then
c           ============================
c             Curl written in R, Y, phi components
            print *, 'FOUCOJ: CAR element type to write'
	      stop

            end if
c           ======

        else
c       ~~~~      
c      Constant k// coordinates:

            if(eletyp.eq.'M23' .or. (eletyp.eq.'HEC'.and. rtp))then
c           =======================================================
      ncoexp = 5
      syv(1) = -1
      syv(2) = -1
      syv(3) = 1
      syv(4) = 1
      syv(5) = 1
        if(rtp)then
        ncoexp = 7
        syv(6) =  1
        syv(7) =  1
        end if
            else if(eletyp.eq.'HEC')then
c           ============================
c             Curl written in rho, eta, parallel components
              ncoexp = 6
              syv(1) = -1
              syv(2) =  1
              syv(3) = -1
              syv(4) =  1
              syv(5) =  1
              syv(6) =  1

            else if(eletyp.eq.'CAR')then
c           ============================
c             Curl written in R, Y, phi components
            print *, 'FOUCOJ: const. k// coordinates to write for CAR element type'
	      stop
            end if
c           ======

        end if
c       ~~~~~~       

c     Now build tables of coefficients to transform:
      call zset(ncorot*(npfft+1), czero, ccr, 1)

        if(.not.cokpco)then
c       ~~~~~~~~~~~~~~~~~~~        
c       Case of theta, phi Fourier expansion
       
C     Normalization for Cray CFFT2 or IMSL routines: fft over Theta
      t = 1.d0 / dfloat(npfft)

c     (NB: case of 'M23' element type was dealt with above.)
            if(eletyp.eq.'HEC')then
c           =======================
      i = intab
      do ipo = 1, npp
c     ###############
	r = eqt(i,ipo,1)
	drrho = eqt(i,ipo,3)
	drthn = eqt(i,ipo,4)
	ntn = eqt(i,ipo,7)
	g12n = eqt(i,ipo,8)
	jacn = eqt(i,ipo,9)
	newmu = eqt(i,ipo,10)
	cn = eqt(i,ipo,11)
	sior = eqt(i,ipo,12)
	co = eqt(i,ipo,14)
	si = eqt(i,ipo,15)
	si1 = eqt(i,ipo,16)
      ntni = 1.d0 / ntn
        if(cyl)then
        ri = 0.d0
        aux1 = kprn
        else
        ri = 1.d0 / r
        aux1 = n * rnorm * ri
        end if
      ncrnor = aux1 * co
      nsrnor = aux1 * si
      rnllnr = rnorm * ri * jacn * (g12n * drthn * ntni - ntn * drrho)

      aux1 = si * rnllnr + rnorm * sior * ntni * newmu
C     Rn*(sLc-cLs):
      aux2 = - rnorm * (g12n * sior * ntni - ntn * si1) / (co * jacn)
            
      ccr(ipo,1) = - t * dcmplx(
     ;- rnorm * ntni*ntni * (si*cn - g12n*sior*newmu), ncrnor * jacn * ntni)
      ccr(ipo,2) = - t * rnorm * sior * jacn * ntni*ntni

c      cb33 = eqt(i,ipo,14) * aux1 + aux2

      ccr(ipo,3) = t * dcmplx(co * (rnorm * cn - g12n * newmu * yinv)
     ;                      , nsrnor * jacn * ntni)
      ccr(ipo,4) = - t * co * jacn * ntni * ntni * yinv

c     For currents:
      ccr(ipo,5) = facju * t * ri * co
      ccr(ipo,6) = facju * t * ri * si

      end do
c     ######

            end if
c           ======
        else
c       ~~~~       
c       Case of Thbar, Phibar Fourier expansion

      i = intab

            if(eletyp.eq.'M23' .or. rtp)then
c           ================================
      do ipo = 1, npp
c     ###############
            
c     Normalization for CFFT2 or IMSL routine: fft over Thbar
c     factor 1 / (NPFFT*dThbar/dtheta): 
      t = 1.d0 / (dfloat(npfft)*ckt(i,ipo,1))
	r = eqt(i,ipo,1)
	drrho = eqt(i,ipo,3)
	drthn = eqt(i,ipo,4)
	ntn = eqt(i,ipo,7)
	g12n = eqt(i,ipo,8)
	jacn = eqt(i,ipo,9)
	newmu = eqt(i,ipo,10)
	cn = eqt(i,ipo,11)
	sior = eqt(i,ipo,12)
	co = eqt(i,ipo,14)
	si = eqt(i,ipo,15)
	si1 = eqt(i,ipo,16)
      ntni = 1.d0 / ntn
        if(cyl)then
        ri = 0.d0
        else
        ri = 1.d0 / r
        end if
      
c      ontpbt = ckt(i,ipo,2) * ntni * yinv
      jon = jacn * ntni
      
        if(cyl)then
        ccr(ipo,1) = dcmplx(0.d0, - t * jon * kprn)
        else
        ccr(ipo,1) = dcmplx(0.d0, - t * jon * n * rnorm * ri)
        end if
      ccr(ipo,2) = t * ntni*ntni * 
     ;dcmplx(rnorm*cn - g12n * newmu * yinv, - n * ckt(i,ipo,2) * jacn * yinv)
      ccr(ipo,3) = - t * hachi(i) * yinv / sior * jon

c     For currents:
      ccr(ipo,4) =  facju * t * ri * cdeid(- n * ckt(i,ipo,4))
c      ccr(ipo,5) =  facju * t *

c     For eventual transfo to rtp:
      ccr(ipo,6) = t * co
      ccr(ipo,7) = t * si

      end do
c     ######
            else if(eletyp.eq.'HEC')then
C           ============================
      do ipo = 1, npp
c     ###############
            
c     Normalization for CFFT2 or IMSL routine: fft over Thbar
c     factor 1 / (NPFFT*dThbar/dtheta): 
      t = 1.d0 / (dfloat(npfft)*ckt(i,ipo,1))
	r = eqt(i,ipo,1)
	drrho = eqt(i,ipo,3)
	drthn = eqt(i,ipo,4)
	ntn = eqt(i,ipo,7)
	g12n = eqt(i,ipo,8)
	jacn = eqt(i,ipo,9)
	newmu = eqt(i,ipo,10)
	cn = eqt(i,ipo,11)
	sior = eqt(i,ipo,12)
	co = eqt(i,ipo,14)
	si = eqt(i,ipo,15)
	si1 = eqt(i,ipo,16)
      ntni = 1.d0 / ntn
        if(cyl)then
        ri = 0.d0
cERN       NB: 1/R properly omitted from eqta1d(, 4):
c       NB: 1/R properly omitted from qfactor:
cERN        tr = kprn * eqta1d(i,4) * hachi(i)
        tr = kprn * qfactor(i) * hachi(i)
        aux4 = (tr * co - kprn) / si
        else
        ri = 1.d0 / r
cERN        tr = n * rnorm * eqta1d(i,4) * hachi(i)
        tr = n * rnorm * qfactor(i) * hachi(i)
C       aux4=(qhc-1/R)/s * n * rnorm:
        aux4 = (tr * co - n * rnorm * ri) / si
        end if
      jon = jacn * ntni
      
      ccr(ipo,1) = - jon * t * dcmplx(
     ;rnorm * ntni / jacn * (- si * cn + sior * g12n * newmu), tr)
      ccr(ipo,2) = - jon * t * hachi(i) * rnorm

      ccr(ipo,3) = jon * t * dcmplx(
     ;co * (rnorm * cn - g12n*newmu * yinv) * ntni / jacn, - aux4)
      ccr(ipo,4) = - jon * t * hachi(i) * co * yinv / sior

c     For currents:
      ccr(ipo,5) =  facju * ri * co
      ccr(ipo,6) =  facju * ri * si

      end do
c     ######
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

c     Using symmetries:
      if(updsym)then
        do l = 1, ncoexp
          if(syv(l).eq.1)then
            do ipo = 1, npp-1
            ccr(npp+ipo,l) = dconjg(ccr(npp-ipo,l))
            end do
          else
            do ipo = 1, npp-1
            ccr(npp+ipo,l) = - dconjg(ccr(npp-ipo,l))
            end do
          end if
        end do
      end if

C     Call fft,n=2**m case: voir cas sym, antisym;
c     ------------------------------------------- ci-dessous general en cx.

cccccccccccccccccccccccc ERNESTO ccccccccccccccccccccccccccccc
      do l = 1, ncoexp
cERN     call cfft2(0, 1, npfft, ccr(1,l), cworkp, cvc(0,l))
	call df2tcb(npfft, ccr(1,l), cvc(0,l), cwork2, cpy) ! -- IMSL --
      end do

c     print *, cvc(1:10,1)
c	print *, 'Stop at FOUCOJ Line 424'
c	stop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l = 1, ncoexp
        do i = 1, npft2
        cvc(-i,l) = cvc(npfft-i,l)
        end do
      end do

c     Now cvc(k,l) contains harmonic K><0 of term #L.

      return
      end
