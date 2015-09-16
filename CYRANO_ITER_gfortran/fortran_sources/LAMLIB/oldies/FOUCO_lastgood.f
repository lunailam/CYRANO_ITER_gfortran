      subroutine fouco

      implicit none

c     Computes analytic Fourier coeffs. used in curl.curl term
c     (circular concentric cross section)
c     The notations for the coeffs(k) are as in thesis Appendix 1 for
c     coeffs(-k); fsyfou is f/y.

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
	complex*16 :: acokfou(0:npfft), bcokfou(0:npfft)

        do 4 i = 0, ncrot
        afou(i) = 0.d0
        bfou(i) = 0.d0
        cfou(i) = 0.d0
        dfou(i) = 0.d0
        efou(i) = 0.d0
        fsyfou(i) = 0.d0
        gfou(i) = 0.d0
        hfou(i) = 0.d0
        jfou(i) = 0.d0
   4    continue

        afou(0) = 1.d0
        bfou(0) = 1.d0

        if(cyl)return
c       rho / Ra:
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

          do i = 1, ncrot+2
          bfou(i) = bfou(0) * t ** i
          end do

        cfou(0) = bfou(1)
        bfou(-1) = bfou(1)
        fsyfou(0) = (bfou(0)-bfou(2))*yp*0.5d0
        hfou(0) = 0.5d0*(bfou(0)+bfou(2))

        do 2 i = 1, ncrot
        cfou(i) = 0.5d0 * (bfou(i+1) + bfou(i-1))

c Thesis:
c        DFOU(I) = YP * 0.5D0 * (BFOU(I+1) - BFOU(I-1))
c Opposite to thesis:
cERN_17Nov04       dfou(i) = yp * 0.5d0 * (bfou(i-1) - bfou(i+1))
c	Restore thesis convention
	  dfou(i) = yp * 0.5d0 * (bfou(i+1) - bfou(i-1))

        t2 = (bfou(i+2) + bfou(i-2)) * 0.5d0
        fsyfou(i) = yp * 0.5d0 * (bfou(i) - t2)
        hfou(i) = (bfou(i) + t2) * 0.5d0
c JFOU not used???
        jfou(i) = yp * 0.25d0 * (bfou(i-2) - bfou(i+2))
  2     continue

c       Normalization:
        do i = 0, ncrot
        cfou(i) = cfou(i) * rnor0
        gfou(i) = gfou(i) * rnor0
        fsyfou(i) = fsyfou(i) * rnor0
        hfou(i) = hfou(i) * rnor02
      end do

cERN	New coefficients for cokpco case

c	if(cokpco) then

c	 1) acokfou = FFT(dPhibar/dtheta)
c	call df2tcb(npfft, dcmplx(ckt(intab,1:2*npp-1,2),0), acokfou, cwork2, cpy) ! -- IMSL --
c	 2) bcokfou = FFT(dPhibar/drho)
c	call df2tcb(npfft, dcmplx(ckt(intab,1:2*npp-1,7),0), bcokfou, cwork2, cpy) ! -- IMSL --

c        do i = 1, npft2
c        bcokfou(-i) = bcokfou(npfft-i)
c        end do

c	end if

c	if (abscis(intab)>0.1)then
c	write(888,"(f15.5)"),ckt(intab,1:2*npp-1,7)
c	write(999,"(f15.5, f15.5)"),real(bcokfou(:)),imag(bcokfou(:))
c	print *, npfft
c	stop
c	end if


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
