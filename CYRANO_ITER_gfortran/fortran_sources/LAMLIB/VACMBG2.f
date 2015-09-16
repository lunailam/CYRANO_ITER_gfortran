      subroutine vacmbg(v, v2, separ)

      implicit none
      logical separ
      complex*16 v(6,6), v2(6,6)

c     Wave operator in vacuum:
c     Computes submatrix v for harmonic k related to mode m with mpk = m + k
c     Components are rho,rho',theta,theta',phi,phi' for 'M23' element
c     +,+',-,-',//,//' for 'HEC' element
c     R,R',Y,Y',phi,phi' for 'CAR' element
c     General or D-shaped equilibrium; 
c     Fourier coefficients were evaluated in FOUCOG.
c     On output, V2 includes the coefficients of 1/y when SEPAR is true.

      include 'pardim.copy'
      include 'commod.copy'
      include 'comswe.copy'
      include 'comgeo.copy'
      include 'comfou.copy'
      include 'comrot.copy'
      include 'comfin.copy'
      include 'comsub.copy'
      include 'comin2.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comphy.copy'

      character*3 eletyp
      
      integer ini, k2, l, ic, j, ka, ii, ji, idest(9), id, jd, nold
      
      double precision tra, rno2r0, kpoln, a1
      
cPL4/5/04: idest, cimta, cimpk dimension 9 for CAR element type
      complex*16 
     ;  cita(nficom,8), citac(nficom,8)
     ;, cimta(9), cimpk(9)
     ;, rf, cf, rft, cft 
     ;, cmat(6,6), cmati(6,6)

cPL4/5/04: idest now initialized below
c     IDEST for sending coefficients stored into vector of length 8
C       (Erho,   dErho/dtheta,
C        Etheta, dEtheta/drho, dEtheta/dtheta,
C        Ephi,   dEphi/drho,   dEphi/dtheta)
C     into vector of length 6 (Erho,Erho',Etheta,Etheta',Ephi,Ephi'):      
c      data idest/1,1,3,4,3,5,6,5,0/

      data ini/0/
      save ini, idest, cita, citac, cimta, cimpk, cmat, cmati
      
      if(ini.eq.0)then

c     6 by 6 matrices for conversions to/from +,-,//:
      call zset(6*6, czero, cmat, 1)
      call zset(6*6, czero, cmati, 1)
      cmat(1,1) = dcmplx(sqrt2i,0.d0)
      cmat(1,3) = dcmplx(sqrt2i,0.d0)
      cmat(2,2) = dcmplx(sqrt2i,0.d0)
      cmat(2,4) = dcmplx(sqrt2i,0.d0)
      cmat(3,1) = dcmplx(0.d0,-sqrt2i)
      cmat(3,3) = dcmplx(0.d0, sqrt2i)
      cmat(4,2) = dcmplx(0.d0,-sqrt2i)
      cmat(4,4) = dcmplx(0.d0, sqrt2i)
      cmat(5,5) = cun
      cmat(6,6) = cun

      cmati(1,1) = dcmplx(sqrt2i,0.d0)
      cmati(3,1) = dcmplx(sqrt2i,0.d0)
      cmati(2,2) = dcmplx(sqrt2i,0.d0)
      cmati(4,2) = dcmplx(sqrt2i,0.d0)
      cmati(1,3) = dcmplx(0.d0, sqrt2i)
      cmati(3,3) = dcmplx(0.d0,-sqrt2i)
      cmati(2,4) = dcmplx(0.d0, sqrt2i)
      cmati(4,4) = dcmplx(0.d0,-sqrt2i)
      cmati(5,5) = cun
      cmati(6,6) = cun

c-------------------------------------------------------------------------------
C     CITA: factors of sqrt(-1) in vector of length 8 given above,
C     resulting from toroidal angle derivatives.
C     CITAC: complex conjugate of CITA.
C     Both only used for M23 elements and COKPCO=.F., 'for historical reasons'
c-------------------------------------------------------------------------------
        do ic = 1, nficom
cERN          do j = 1, 8
        cita(ic,1:8) = cun
        citac(ic,1:8) = cun
cERN          end do

          if(nfaci(ic).gt.0)then
          j = ifaci(ic,1)
          cita(ic,j) = ci
          citac(ic,j) = - ci
          end if
        end do
cPL4/5/04 now done below:
c        do j = 1, 9
c        cimta(j) = cun
c        cimpk(j) = cun
c        end do
      ini = 1
      end if

      eletyp = styp(isubr)
      call zset(36, czero, v, 1)
      
        if(separ)then
        call zset(36, czero, v2, 1)
        else
        yinv = abscni(intab)
        end if
 
      if(geneq.or.dshape)then
c     #######################

	if(eletyp .eq. 'CAR')then
C     IDEST for sending coefficients stored into vector of length 9
C       (ER,   dER/drho,   dER/dtheta,
C        EY,   dEY/drho,   dEY/dtheta,
C        Ephi, dEphi/drho, dEphi/dtheta)
C     into vector of length 6 (ER,ER',EY,EY',Ephi,Ephi'):      
      idest(1) = 1
      idest(2) = 2
      idest(3) = 1
      idest(4) = 3
      idest(5) = 4
      idest(6) = 3
      idest(7) = 5
      idest(8) = 6
      idest(9) = 5
c-------------------------------------------------------------------------------
C     CIMPK: factors of i*(m+k) in vector of length 9 given above.
C     CIMTA: factors of -i*m in vector of length 9.
C     They result from poloidal Fourier expansion of rf field (CIMPK)
C     and test function (CIMTA).
c-------------------------------------------------------------------------------
      cimpk(1) = cun
      cimpk(2) = cun
      cimpk(3) = dcmplx(0.d0,dfloat(mpk))
      cimpk(4) = cun
      cimpk(5) = cun
      cimpk(6) = cimpk(3)
      cimpk(7) = cun
      cimpk(8) = cun
      cimpk(9) = cimpk(3)
      cimta(1) = cun
      cimta(2) = cun
      cimta(3) = dcmplx(0.d0,-dfloat(m))
      cimta(4) = cun
      cimta(5) = cun
      cimta(6) = cimta(3)
      cimta(7) = cun
      cimta(8) = cun
      cimta(9) = cimta(3)

	else
c     IDEST for sending coefficients stored into vector of length 8
C       (Erho,   dErho/dtheta,
C        Etheta, dEtheta/drho, dEtheta/dtheta,
C        Ephi,   dEphi/drho,   dEphi/dtheta)
C     into vector of length 6 (Erho,Erho',Etheta,Etheta',Ephi,Ephi'):      
      idest(1) = 1
      idest(2) = 1
      idest(3) = 3
      idest(4) = 4
      idest(5) = 3
      idest(6) = 5
      idest(7) = 6
      idest(8) = 5
	
c-------------------------------------------------------------------------------
C     CIMPK: factors of i*(m+k) in vector of length 8 given above.
C     CIMTA: factors of -i*m in vector of length 8.
C     They result from poloidal Fourier expansion of rf field (CIMPK)
C     and test function (CIMTA).
c-------------------------------------------------------------------------------
      cimpk(1) = cun
      cimpk(2) = dcmplx(0.d0,dfloat(mpk))
      cimpk(3) = cun
      cimpk(4) = cun
      cimpk(5) = cimpk(2)
      cimpk(6) = cun
      cimpk(7) = cun
      cimpk(8) = cimpk(2)
      cimta(1) = cun
      cimta(2) = dcmplx(0.d0,-dfloat(m))
      cimta(3) = cun
      cimta(4) = cun
      cimta(5) = cimta(2)
      cimta(6) = cun
      cimta(7) = cun
      cimta(8) = cimta(2)
	end if

      l = 1
c     Displacement current:
      if(.not.quasis)then
      tra = - k0rn2 * cvc(k,1)
      v(1,1) = v(1,1) + tra
      v(3,3) = v(3,3) + tra
      v(5,5) = v(5,5) + tra
      end if

      if(eletyp.eq.'M23' .and. .not.cokpco)then
c     ========================================c     This case was programmed using only fft of real coefficients.
c     Nonsingular matrix elements:
      do ic = 1, nficom
c       Row sweep:
        do j = 1, nns(ic)
        ii = ins(j,ic)
        rf = citac(ic,ii) * cimta(ii)
        rft = cita(ic,ii) * cimpk(ii)
        id = idest(ii)
c         Column sweep:
          do ka = j, nns(ic)
          l = l + 1
          ji = ins(ka,ic)
          cf = cita(ic,ji) * cimpk(ji)
          jd = idest(ji)
          v(id,jd) = v(id,jd) + cvc(k,l) * rf * cf
            if(ka .gt. j)then
c           Using symmetry of curl.curl:
            cft = citac(ic,ji) * cimta(ji)
            v(jd,id) = v(jd,id) + cvc(k,l) * rft * cft
            end if
          end do
          do ka = 1, nsi(ic)
          l = l + 1
          ji = isi(ka,ic)
          cf = cita(ic,ji) * cimpk(ji)
          cft = citac(ic,ji) * cimta(ji)
          jd = idest(ji)
          v(id,jd) = v(id,jd) + cvc(k,l) * rf * cf
          v(jd,id) = v(jd,id) + cvc(k,l) * rft * cft
          end do
        end do
      end do

      if(separ)then
c     Singular matrix elements * y are stored in V2: final division by y
c     or appropriate treatment is deferred.
      do ic = 1, nficom
        do j = 1, nsi(ic)
        ii = isi(j,ic)
        rf = citac(ic,ii) * cimta(ii)
        rft = cita(ic,ii) * cimpk(ii)
        id = idest(ii)
          do ka = j, nsi(ic)
          l = l + 1
          ji = isi(ka,ic)
          cf = cita(ic,ji) * cimpk(ji)
          jd = idest(ji)
          v2(id,jd) = v2(id,jd) + cvc(k,l) * rf * cf
            if(ka.gt.j)then
c           Using symmetry of curl.curl:
            cft = citac(ic,ji) * cimta(ji)
            v2(jd,id) = v2(jd,id) + cvc(k,l) * rft * cft
            end if
          end do
        end do
      end do
      
      else
c     All matrix elements are stored in V; division by y is performed here:
      do ic = 1, nficom
        do j = 1, nsi(ic)
        ii = isi(j,ic)
        rf = citac(ic,ii) * cimta(ii) * yinv
        rft = cita(ic,ii) * cimpk(ii) * yinv
        id = idest(ii)
          do ka = j, nsi(ic)
          l = l + 1
          ji = isi(ka,ic)
          cf = cita(ic,ji) * cimpk(ji)
          jd = idest(ji)
          v(id,jd) = v(id,jd) + cvc(k,l) * rf * cf
            if(ka.gt.j)then
c           Using symmetry of curl.curl:
            cft = citac(ic,ji) * cimta(ji)
            v(jd,id) = v(jd,id) + cvc(k,l) * rft * cft
            end if
          end do
        end do
      end do
      
      end if
      
      else
c     ===c     Cases: 'HEC' elements with COKPCO=.FALSE. or .TRUE.; or 'M23' with COKPCO=.T.
c     Written using FFT of complex coefficients.

c     Nonsingular matrix elements:
      do ic = 1, nficom
        do j = 1, nns(ic)
        ii = ins(j,ic)
        rf = cimta(ii)
        rft = cimpk(ii)
        id = idest(ii)
        
          do ka = j, nns(ic)
          l = l + 1
          ji = ins(ka,ic)
          cf = cimpk(ji)
          jd = idest(ji)
          v(id,jd) = v(id,jd) + cvc(k,l) * rf * cf
            if(ka.gt.j)then
            cft = cimta(ji)
c  corr:
            v(jd,id) = v(jd,id) + dconjg(cvc(-k,l)) * rft * cft
c            v(jd,id) = v(jd,id) + cvc(k,l) * rft * cft
            end if
          end do
          do ka = 1, nsi(ic)
          l = l + 1
          ji = isi(ka,ic)
          cf = cimpk(ji)
c  corr:
          cft = cimta(ji)
          jd = idest(ji)
          v(id,jd) = v(id,jd) + cvc(k,l) * rf * cf
c  corr:
          v(jd,id) = v(jd,id) + dconjg(cvc(-k,l)) * rft * cft
c          v(jd,id) = v(jd,id) + cvc(k,l) * rft * cft
          end do
        end do
      end do

      if(separ)then
c     Singular matrix elements * y are stored in V2: final division by y
c     or appropriate treatment is deferred.
      do ic = 1, nficom
        do j = 1, nsi(ic)
        ii = isi(j,ic)
        rf = cimta(ii)
        rft = cimpk(ii)
        id = idest(ii)
          do ka = j, nsi(ic)
          l = l + 1
          ji = isi(ka,ic)
          cf = cimpk(ji)
          jd = idest(ji)
          v2(id,jd) = v2(id,jd) + cvc(k,l) * rf * cf
            if(ka.gt.j)then
            cft = cimta(ji)
            v2(jd,id) = v2(jd,id) + dconjg(cvc(-k,l)) * rft * cft
c            v2(jd,id) = v2(jd,id) + cvc(k,l) * rft * cft
            end if
          end do
        end do
      end do
      
      else
c     All matrix elements are stored in V; division by y is performed here:
      do ic = 1, nficom
        do j = 1, nsi(ic)
        ii = isi(j,ic)
        rf = cimta(ii) * yinv
        rft = cimpk(ii) * yinv
        id = idest(ii)
          do ka = j, nsi(ic)
          l = l + 1
          ji = isi(ka,ic)
          cf = cimpk(ji)
          jd = idest(ji)
          v(id,jd) = v(id,jd) + cvc(k,l) * rf * cf
            if(ka.gt.j)then
            cft = cimta(ji)
            v(jd,id) = v(jd,id) + dconjg(cvc(-k,l)) * rft * cft
c            v(jd,id) = v(jd,id) + cvc(k,l) * rft * cft
            end if
          end do
        end do
      end do
      
      end if
      
      end if
c     ====c     Eventual transformation from rho, eta, // to +,-,// components:
        if(eletyp.eq.'HEC')then
cERN        call mul3(v, cmati, v, cmat, 6, 6)
		  call mul3_fast(v, cmati, v, cmat, 6, 6)
            if(separ)then 
cERN		     call mul3(v2, cmati, v2, cmat, 6, 6)
			 call mul3_fast(v2, cmati, v2, cmat, 6, 6)
		  end if
	  end if

      else  ! Circular concentric case: 
c     ####

      rno2r0 = rnor0 * 0.5d0
      k2 = iabs(k)
 
      if(k.eq.0)then
c     Uses AFOU(0)=1.
      v(1,1) = dcmplx(bfou(0) * y * kprn2, 0.d0)
      v(1,4) = dcmplx(0.d0, dfloat(m))
      v(1,5) = dcmplx(0.d0, y*cfou(0)*kprn)
      v(1,6) = dcmplx(0.d0, y*kprn)
      v(4,1) = - v(1,4)
      v(5,1) = - v(1,5)
      v(6,1) = - v(1,6)
      v(3,3) = v(1,1)
      v(3,4) = cun
      v(4,3) = cun
      v(4,4) = dcmplx(y, 0.d0)
      v(3,5) = dcmplx(- m * kprn, 0.d0)
      v(5,3) = v(3,5)
      v(5,5) = dcmplx(bfou(0) * rnor0**2 * y, 0.d0)
      v(6,6) = v(4,4)
 
        if(separ)then
        v2(1,1) = dcmplx(dfloat(m*m), 0.d0)
        v2(1,3) = dcmplx(0.d0, dfloat(m))
        v2(3,1) = - v2(1,3)
        v2(3,3) = cun
        v2(5,5) = v2(1,1)
        else
        kpoln = m * yinv
        v(1,1) = v(1,1) + dcmplx(kpoln * m, 0.d0)
        v(1,3) = v(1,3) + dcmplx(0.d0, kpoln)
        v(3,1) = v(3,1) - dcmplx(0.d0, kpoln)
        v(3,3) = v(3,3) + dcmplx(yinv, 0.d0)
        v(5,5) = v(5,5) + dcmplx(kpoln * m, 0.d0)
        end if
 
      else if(k2.eq.1)then
c     Cases k=±1:
      a1 = rno2r0 * y
      v(1,1) = dcmplx(bfou(1) * y * kprn2, 0.d0)
      v(1,5) = dcmplx(0.d0, y * cfou(1) * kprn)
      v(5,1) = - v(1,5)
      v(3,3) = v(1,1)
cPL-ERN25Nov04 flip sign
cOLD      v(3,5) = dcmplx(- dfou(1) * kprn, 0.d0)
cPL28/1/05 use negative harmonics
c      v(3,5) = dcmplx(dfou(1) * kprn, 0.d0)
      v(3,5) = dcmplx(dfou(k) * kprn, 0.d0)
cPL28/1/05 use negative harmonics
c      if(k.lt.0)v(3,5) = - v(3,5)
      v(5,3) = - v(3,5)
cOLD      v(5,3) = v(3,5)
      v(5,5) = dcmplx(bfou(1) * rnor0**2 * y, 0.d0)
      v(5,6) = dcmplx(y * gfou(1), 0.d0)
      v(6,5) = v(5,6)
      
      v(1,1) = v(1,1) +  dcmplx(rno2r0 * mpk * m, 0.d0)
      v(1,3) = dcmplx(0.d0, rno2r0*m)
      v(1,4) = dcmplx(0.d0, m*a1)
      v(3,1) = dcmplx(0.d0, -rno2r0*mpk)
      v(4,1) = dcmplx(0.d0, -mpk*a1)
      v(3,3) = v(3,3) +  dcmplx(rno2r0, 0.d0)
      v(3,4) = dcmplx(a1, 0.d0)
      v(4,3) = v(3,4)
      v(4,4) = dcmplx(a1*y, 0.d0)
      v(5,5) = v(5,5) + dcmplx(rno2r0*(m*mpk+1), 0.d0)
      v(6,6) = v(4,4)
 
      else
      v(1,1) = dcmplx(bfou(k2)*y*kprn2, 0.d0)
      v(1,5) = dcmplx(0.d0, y*cfou(k2)*kprn)
      v(5,1) = - v(1,5)
      v(3,3) = v(1,1)
cPL-ERN25Nov04 flip sign
cOLD      v(3,5) = dcmplx(-dfou(k2)*kprn, 0.d0)
cPL28/1/05 use negative harmonics
c      v(3,5) = dcmplx(dfou(k2)*kprn, 0.d0)
      v(3,5) = dcmplx(dfou(k)*kprn, 0.d0)
cPL28/1/05 use negative harmonics
c      if(k.lt.0)v(3,5) = - v(3,5)
      v(5,3) = - v(3,5)
      v(5,5) = dcmplx(bfou(k2) * rnor0**2 * y, 0.d0)
      end if
 


cERN	-----------------------------------------------------------------------------------

cERN	NEW: Add k//=const. corrections to Vij elements if cokpco=TRUE 
c	     Obs: Fourier coef. already computed in FOUCO.f (common COMFOU)
c		 Storage: acokfou[k=-npft/2,...,-1,0,+1,....+npft/2] (complex)
c	nold=n
c	n=1
	if(cokpco) then
c	................

c	1) Radial contributions:
c	if(k.gt.0) then
	v(3,5) = v(3,5) - n * kprn * acokfou(k)		 ! ~ FFT(  dPhibar/dtheta) - even, real
	v(5,3) = v(5,3) - n * kprn * acokfou(k)		 ! ~ FFT(  dPhibar/dtheta) - even, real
c	else
c	v(3,5) = v(3,5) + n * kprn * acokfou(k)		 ! ~ FFT(  dPhibar/dtheta)
c	v(5,3) = v(5,3) + n * kprn * acokfou(k)		 ! ~ FFT(  dPhibar/dtheta)
c	end if

	v(5,5) = v(5,5) + n*n*yinv * ccokfou(k)		 ! ~ FFT(R/Ro.(dPhibar/dtheta)^2) - even, real
     ;                + (m+mpk) * n * yinv * bcokfou(k)
cPL6Dec04 Added 1 line above

c	2) Poloidal contributions:

	v(1,5) = v(1,5) - n * kprn * y * dcokfou(k)	 ! ~ FFT(a.dPhibar/drho) 
	v(5,1) = v(5,1) - n * kprn * y * dcokfou(k)	 ! ~ FFT(a.dPhibar/drho) 

	v(5,6) = v(5,6) - ci * n * y * ecokfou(k)	       ! ~ FFT(a.R/Ro.dPhibar/drho) 
	v(6,5) = v(6,5) + ci * n * y * ecokfou(k)	       ! ~ FFT(a.R/Ro.dPhibar/drho) 

	v(5,5) = v(5,5)  + n*n * y * fcokfou(k)		 ! ~ FFT(a2.R/Ro.(dPhibar/drho)^2) - even, real


c	3) Toroidal contributions:
	v(1,1) = v(1,1) + n*n*yinv * ccokfou(k)		 ! ~ FFT(R/Ro.(dPhibar/dtheta)^2)  - even, real
     ;                + (m+mpk) * n * yinv * bcokfou(k)
cPL6Dec04 Added 1 line above
	v(3,3) = v(3,3) + n*n * y * fcokfou(k)		 ! ~ FFT(a2.R/Ro.(dPhibar/drho)^2) - even, real

c	if(k.gt.0)then  
	v(1,3) = v(1,3) + ci * n * yinv * bcokfou(k)     ! ~ FFT(R. dPhibar/dtheta))
     ;                - n * m * ecokfou(k)		 ! ~ FFT(a.R/Ro. dPhibar/drho))
     ;                - n * n * gcokfou(k)             ! ~ FFT(a.R/Ro.dPhibar/dtheta.dPhibar/drho))

	v(3,1) = v(3,1) - ci * n * yinv * bcokfou(k)     ! ~ FFT(R. dPhibar/dtheta))
     ;                - n * mpk * ecokfou(k)		 ! ~ FFT(a.R/Ro. dPhibar/drho))
     ;                - n * n * gcokfou(k)             ! ~ FFT(a.R/Ro.dPhibar/dtheta.dPhibar/drho))

c	else
c	v(1,3) = v(1,3) + ci * n * yinv * bcokfou(k) ! ~ FFT(R. dPhibar/dtheta))
c     ;                + n * m * ecokfou(k)		 ! ~ FFT(a.R/Ro. dPhibar/drho))
c     ;				+ n * n * gcokfou(k)		 ! ~ FFT(a.R/Ro.dPhibar/dtheta.dPhibar/drho))

c	v(3,1) = v(3,1) - ci * n * yinv * bcokfou(k) ! ~ FFT(R. dPhibar/dtheta))
c     ;                + n * mpk * ecokfou(k)		 ! ~ FFT(a.R/Ro. dPhibar/drho))
c     ;				+ n * n * gcokfou(k)		 ! ~ FFT(a.R/Ro.dPhibar/dtheta.dPhibar/drho))

c	end if 
 
      v(1,4) = v(1,4) + ci * n * bcokfou(k)		 ! ~ FFT(R/Ro.dPhibar/dtheta) - even, real
      v(4,1) = v(4,1) - ci * n * bcokfou(k)		 ! ~ FFT(R/Ro.dPhibar/dtheta) - even, real

      v(3,4) = v(3,4) - ci * n * y * ecokfou(k)	 ! ~ FFT(a.R/Ro.dPhibar/drho)
      v(4,3) = v(4,3) + ci * n * y * ecokfou(k)	 ! ~ FFT(a.R/Ro.dPhibar/drho)

c	  if(abscis(intab)>0.1) then
c	     print *, 'tempo....'
c	  end if

	end if !(cokpco)
c	................	
c	n=nold
c	--------------------------------------------------------------------------------------
c     Displacement current:

c	if(cokpco) then
c        if(.not.quasis)then
c        tra = - k0rn2 * y * afou(k2)
c        v(1,1) = v(1,1) + dcmplx(tra, 0.d0) - n*n*yinv * ccokfou(k)	
c        v(3,3) = v(3,3) + dcmplx(tra, 0.d0) - n*n *   rnorm**2 * y * fcokfou(k)
c        v(5,5) = v(5,5) + dcmplx(tra, 0.d0) - n*n*yinv * ccokfou(k) - n*n*rnorm**2 * y * fcokfou(k)	
c        end if
c	else
        if(.not.quasis)then
        tra = - k0rn2 * y * afou(k2)
        v(1,1) = v(1,1) + dcmplx(tra, 0.d0)
        v(3,3) = v(3,3) + dcmplx(tra, 0.d0)
        v(5,5) = v(5,5) + dcmplx(tra, 0.d0)
        end if
c	end if

c     Eventual transformation to +,-,// components:
         if(eletyp.eq.'HEC')then
cERN         call mul3(v, bch(1,1,ig), v, bc(1,1,ig), 6, 6)
		   call mul3_fast(v, bch(1,1,ig), v, bc(1,1,ig), 6, 6)
             if(separ)then
cERN	          call mul3(v2, bch(1,1,ig), v2, bc(1,1,ig), 6, 6)
			  call mul3_fast(v2, bch(1,1,ig), v2, bc(1,1,ig), 6, 6)	
		   end if
	   end if

      end if	! (dshape or circular)
c     ######
 
      return
      end
