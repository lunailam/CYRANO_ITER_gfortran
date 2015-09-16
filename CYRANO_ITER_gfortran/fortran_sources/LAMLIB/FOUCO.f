      subroutine fouco

      implicit none

c     Computes analytic Fourier coeffs. used in curl.curl term
c     (circular concentric cross section)
c     The notations for the coeffs(k) are as in thesis Appendix 1 for
c     coeffs(-k), i.e. the forward-backward FT convention is opposite; 
c     fsyfou is f/y.

      include 'pardim.copy'
      include 'comfou.copy'
      include 'comgeo.copy'
      include 'commod.copy'
      include 'commag.copy'
      include 'compla.copy'
	include 'comfft.copy'
	include 'comswe.copy'

      integer i
      double precision t, t2, yp


c        do i = 0, ncrot
c        afou(i) = 0.d0
c        bfou(i) = 0.d0
c        cfou(i) = 0.d0
c        dfou(i) = 0.d0
c        efou(i) = 0.d0
c        fsyfou(i) = 0.d0
c        gfou(i) = 0.d0
c        hfou(i) = 0.d0
c        jfou(i) = 0.d0
c        end do
cPL28/1/05 Store negative harmonics as well, and be happier:
        afou(-ncrot:ncrot) = 0.d0
        bfou(-ncrot:ncrot) = 0.d0
        cfou(-ncrot:ncrot) = 0.d0
        dfou(-ncrot:ncrot) = 0.d0
        efou(-ncrot:ncrot) = 0.d0
        fsyfou(-ncrot:ncrot) = 0.d0
        gfou(-ncrot:ncrot) = 0.d0
        hfou(-ncrot:ncrot) = 0.d0
        jfou(-ncrot:ncrot) = 0.d0


        afou(0) = 1.d0
        bfou(0) = 1.d0

        if(cyl)return
c       rho / ra:
        yp = y * rnor0

        afou(1) = yp * 0.5d0
c Thesis:
c        efou(1) =  yp * 0.5d0
c Opposite to thesis:
        efou(1) = - yp * 0.5d0
        gfou(1) = 0.5d0

        if(yp .lt. 1.d-2)then
        t = - yp * 0.5d0
        else
        t = - (1.d0 - dsqrt(1.d0-yp*yp)) / yp
        end if
        bfou(0) = 1.d0 / dsqrt(1.d0-yp*yp)

c         Two extra terms for bfou, allow recursions for some other terms:
          do i = 1, ncrot+2
          bfou(i) = bfou(0) * t ** i
          end do

        cfou(0) = bfou(1)
	  dfou(0) = 0.d0
        bfou(-1) = bfou(1)
        fsyfou(0) = (bfou(0)-bfou(2))*yp*0.5d0
        hfou(0) = 0.5d0*(bfou(0)+bfou(2))

        do i = 1, ncrot
        cfou(i) = 0.5d0 * (bfou(i+1) + bfou(i-1))
c Thesis:
c        DFOU(I) = YP * 0.5D0 * (BFOU(I+1) - BFOU(I-1))
c Opposite to thesis:
      dfou(i) = yp * 0.5d0 * (bfou(i-1) - bfou(i+1))
c ERN 17Nov04 Restore thesis convention
cPL 23/11/04 undo:
c	  dfou(i) = yp * 0.5d0 * (bfou(i+1) - bfou(i-1))

        t2 = (bfou(i+2) + bfou(i-2)) * 0.5d0
        fsyfou(i) = yp * 0.5d0 * (bfou(i) - t2)
        hfou(i) = (bfou(i) + t2) * 0.5d0
c JFOU not used???
        jfou(i) = yp * 0.25d0 * (bfou(i-2) - bfou(i+2))
        end do

c       Normalization:
        do i = 0, ncrot
        cfou(i) = cfou(i) * rnor0
        gfou(i) = gfou(i) * rnor0
        fsyfou(i) = fsyfou(i) * rnor0
        hfou(i) = hfou(i) * rnor02
        end do

cPL28/1/05 Store negative harmonics:
          do i = 1, ncrot
          afou(-i) = afou(i)
          bfou(-i) = bfou(i)
          cfou(-i) = cfou(i)
          dfou(-i) = -dfou(i)
          efou(-i) = -efou(i)
          fsyfou(-i) = fsyfou(i)
          gfou(-i) = gfou(i)
          hfou(-i) = hfou(i)
          jfou(-i) = jfou(i)
	    end do

cERN	NEW: coefficients for cokpco case ---------------------------------------------------
c	    Use this FFT: (same convention as circular, opposite to thesis)
c		FFT(k) = int ( f(theta) exp(+i.k.theta) dtheta ) -> df2tcb (IMSL)
c

	if(cokpco) then
	
		acokfou = 0
		bcokfou = 0
		ccokfou = 0
		dcokfou = 0
		ecokfou = 0
		fcokfou = 0
		gcokfou = 0

c	 1) acokfou = FFT(dPhibar/dtheta)
		call df2tcb(npfft, dcmplx(ckt(intab,1:npfft,2),0), 
     ;                       acokfou(0:npfft-1), cwork2, cpy) ! -- IMSL --
      acokfou(0) = 0.d0

c	 2) bcokfou = FFT(R/Ro.dPhibar/dtheta)
		call df2tcb(npfft, dcmplx(eqt(intab,1:npfft,1)/r0*ckt(intab,1:npfft,2),0), 
     ;                       bcokfou(0:npfft-1), cwork2, cpy)
c	call df2tcf(npfft, ernaux(1:npfft), bcokfou(0:npfft-1), cwork2, cpy)

c	 3) ccokfou = FFT(R/Ro.(dPhibar/dtheta)^2)
		call df2tcb(npfft, dcmplx(eqt(intab,1:npfft,1)/r0*ckt(intab,1:npfft,2)**2,0), 
     ;                       ccokfou(0:npfft-1), cwork2, cpy) 

c	 4) dcokfou = FFT(a.dPhibar/drho)
		call df2tcb(npfft, dcmplx(rnorm*ckt(intab,1:npfft,7),0), 
     ;                       dcokfou(0:npfft-1), cwork2, cpy) 
     
c	 5) ecokfou = FFT(a.R/Ro.dPhibar/drho)
		call df2tcb(npfft, dcmplx(rnorm/r0*eqt(intab,1:npfft,1)*ckt(intab,1:npfft,7),0), 
     ;                       ecokfou(0:npfft-1), cwork2, cpy) 

c	 6) fcokfou = FFT(a^2.R/Ro.(dPhibar/drho)^2)
		call df2tcb(npfft, dcmplx(rnorm*rnorm/r0*eqt(intab,1:npfft,1)*ckt(intab,1:npfft,7)**2,0), 
     ;                       fcokfou(0:npfft-1), cwork2, cpy) ! -- IMSL --


c	 7) gcokfou = FFT(a.R/Ro.(dPhibar/dtheta)*(dPhibar/drho))
		call df2tcb(npfft, dcmplx(rnorm/r0*eqt(intab,1:npfft,1)*ckt(intab,1:npfft,2)*ckt(intab,1:npfft,7),0), 
     ;		           gcokfou(0:npfft-1), cwork2, cpy) ! -- IMSL --




c	Apply normalization (1/Npoints)

	t=1/dfloat(npfft)

	acokfou(0:npfft) = acokfou(0:npfft) * t
	bcokfou(0:npfft) = bcokfou(0:npfft) * t
	ccokfou(0:npfft) = ccokfou(0:npfft) * t
	dcokfou(0:npfft) = dcokfou(0:npfft) * t
	ecokfou(0:npfft) = ecokfou(0:npfft) * t
	fcokfou(0:npfft) = fcokfou(0:npfft) * t
	gcokfou(0:npfft) = gcokfou(0:npfft) * t

c	Shift values for for k<0 
        do i = 1, npft2
		 acokfou(-i) = acokfou(npfft-i)
		 bcokfou(-i) = bcokfou(npfft-i)
		 ccokfou(-i) = ccokfou(npfft-i)
		 dcokfou(-i) = dcokfou(npfft-i)
		 ecokfou(-i) = ecokfou(npfft-i)
		 fcokfou(-i) = fcokfou(npfft-i)
		 gcokfou(-i) = gcokfou(npfft-i)
        end do

	end if

c	Write Fourier coefficients to files (only DEBUG)
c	(Should be commented OUT in usual runs)
c	if (abscis(intab)>0.1)then
c	do i=1,npfft
c	  write(111,"(f20.8, f20.8)"), polang(i), ckt(intab,i,2)
c	  write(222,"(f20.8, f20.8)"), polang(i), eqt(intab,i,1)/r0*ckt(intab,i,2)
c	  write(333,"(f20.8, f20.8)"), polang(i), eqt(intab,i,1)/r0*ckt(intab,i,2)**2
c	  write(444,"(f20.8, f20.8)"), polang(i), rnorm*ckt(intab,i,7)
c	  write(555,"(f20.8, f20.8)"), polang(i), rnorm/r0*eqt(intab,i,1)*ckt(intab,i,7)
c	  write(666,"(f20.8, f20.8)"), polang(i), rnorm*rnorm/r0*eqt(intab,i,1)*ckt(intab,i,7)**2
c	  write(777,"(f20.8, f20.8)"), polang(i), rnorm/r0*eqt(intab,i,1)*ckt(intab,i,2)*ckt(intab,i,7)
c	end do
c	
c	do i=-npft2,npft2
c	  write(100,"(i5, f20.8, f20.8)"), i, dreal(acokfou(i)), dimag(acokfou(i))
c	  write(200,"(i5, f20.8, f20.8)"), i, dreal(bcokfou(i)), dimag(bcokfou(i))
c	  write(300,"(i5, f20.8, f20.8)"), i, dreal(ccokfou(i)), dimag(ccokfou(i))
c	  write(400,"(i5, f20.8, f20.8)"), i, dreal(dcokfou(i)), dimag(dcokfou(i))
c	  write(500,"(i5, f20.8, f20.8)"), i, dreal(ecokfou(i)), dimag(ecokfou(i))
c	  write(600,"(i5, f20.8, f20.8)"), i, dreal(fcokfou(i)), dimag(fcokfou(i))
c	  write(700,"(i5, f20.8, f20.8)"), i, dreal(gcokfou(i)), dimag(gcokfou(i))
c	end do
c
c	print *, 'Stop in FOUCO'
c	stop
c	end if

c ---------------------------------------------------------------------------------


        return
        end
C
C************************************************
C
CCC        TEST FOUCO
CCC
CC      include 'COMFOU.COPY'
CC      include 'COMMOD.COPY'
CC      include 'NAMPOM.COPY'
CC      DOUBLE PRECISION Y
CCC
CC      Y=0.2
CC      CALL FOUCO(Y)
CC      WRITE(5,100)(AFOU(I),BFOU(I),CFOU(I),DFOU(I),EFOU(I),FSYFOU(I
CC     ;    0,NCROT)
CC      STOP
CC  100   FORMAT(1H ,8(1X,G13.5))
CC      END

