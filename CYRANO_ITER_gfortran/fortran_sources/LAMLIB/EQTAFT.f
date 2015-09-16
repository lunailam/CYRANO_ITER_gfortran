      subroutine eqtaft(index,foutra)

ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)
      implicit none
      
c     Computes poloidal Fourier transform of a term (1<=index<=17)
c     tabulated in 3D array eqt using cray rcfft2
c     Radial table index intab in common comswe.
c     D-shaped or general cross section.

      include 'pardim.copy'
      include 'comfou.copy'
      include 'comgeo.copy'
      include 'commod.copy'
      include 'commag.copy'
      include 'compla.copy'
      include 'comswe.copy'
      include 'comrot.copy'
      include 'comsub.copy'
      include 'comphy.copy'
	include 'comfft.copy'

      integer index
      complex*16 foutra(-npft2:npft2)


c      integer ini
	integer i, ipo

      double precision t1, cr(npfft+1)

      complex*16 work(3*npft2+2)
     
c      data ini/0/

c      save ini
c      if(ini.eq.0)then
c     Init FFT once, OK if always same number of points:
c     Convention: direct transform with exp[+i k theta] !
cERN  call rcfft2(1, 1, npfft, cr, work, cvc)
c      ini = 1
c      end if

      i = intab
C     Normalization for RCFFT2: 
c      t1 = abscno(i) * r0i / dfloat(2 * npfft)
C     Normalization for IMSL: 
      t1 = abscno(i) * r0i / dfloat(npfft)

        if(.not.cokpco)then
c       ~~~~~~~~~~~~~~~~~~~        
c       Case of theta Fourier expansion
c       RCFFT2: factor y*(R/Ra*Jn) / (2*npfft): 
c       IMSL: factor y*(R/Ra*Jn) / npfft: 
          do ipo = 1, npp
          cr(ipo) = t1 * eqt(i,ipo,1) * eqt(i,ipo,9) * eqt(i,ipo,index)
          end do

        else
c       ~~~~      
c       Constant k// coordinates: thbar expansion
c       Normalization for RCFFT2 routine: fft over Thbar
c       factor y*(R/RA*Jn) / (2*npfft*dThbar/dtheta): 
c       Normalization for IMSL routine: fft over Thbar
c       factor y*(R/RA*Jn) / (npfft*dThbar/dtheta): 
          do ipo = 1, npp
          cr(ipo) = t1 * eqt(i,ipo,1) * eqt(i,ipo,9) * eqt(i,ipo,index)
     ;     /  ckt(i,ipo,1)
          end do

        end if
c       ~~~~~~       
  
c     Using symmetries:
      if(updsym)then
        if(uds(index).eq.1)then
          do ipo = 1, npp-1
          cr(npp+ipo) = cr(npp-ipo)
          end do
        else
          do ipo = 1, npp-1
          cr(npp+ipo) = - cr(npp-ipo)
          end do
        end if
      end if

c     Call fft,n=2**m case: voir cas sym, antisym;
c     --------------------------------------------- ci-dessous general en reel.

cccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccccccc
c     WARNING: For FFT of REAL functions, the index of the Fourier coefs. 
c     in CRAY and IMSL routines are different:
cERN      call rcfft2(0, 1, npfft, cr, work, foutra(0))
	call df2trf(npfft, cr, cr, work2) ! -- IMSL --

c     Storage of IMSL transform in complex array:
c     NB we take the complex conjugate following convention on Fourier transform: we want
c     integral of exp(+i k theta)!
	foutra(0) = dcmplx(cr(1), 0.d0)
	  do i = 1, npft2-1
	  foutra(i) = dcmplx(cr(2*i), -cr(2*i+1))
        end do
	foutra(npft2) = dcmplx(cr(npfft), 0.d0)

c      do i = 0,5
c         print *, foutra(i)
c	end do
c	print *, 'STOP in Eqtaft Line 364'
c	stop     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      do i = 1, npft2
      foutra(-i) = dconjg(foutra(i))
      end do

c     Now foutra(K) contains harmonic K><0 of term.
      
      return
      end
