      subroutine fougdr_ern(icase)
      
ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)
      implicit none
      
	integer icase

c     generates Fourier coefficients giving geometrical couplings in 
c     dielectric response at current radial abscissa (intab).
c     icase: input switch selecting what is transformed:
c            icase = 1 for geom. coefficients of Maxwellians
c            icase = 2 for non-Maxwellian.
c  
      include 'pardim.copy'
      include 'comfou.copy'
      include 'comgeo.copy'
      include 'commod.copy'
      include 'commag.copy'
      include 'compla.copy'
      include 'comswe.copy'
      include 'comrot.copy'
      include 'comwrr.copy'
      include 'comphy.copy'
      include 'comfft.copy'
      include 'cokpco.copy'

      integer i, ka, l, ipo, method

	parameter (method=1)

      double precision thbmkh(npfft+1), bvskhi(npfft+1), rnfi, khothb(npfft+1)
     ;, bvsthb(npfft+1), dthdkh(npfft+1), dtdkvk(npfft+1)
     ;, dkhdth(npfft+1)

      complex*16 
     ;  cr(npfft+1,0:maxcou), ei(npfft+1)
     ;, cdeid

      external cdeid
      
c     9/11/03: FFT sine and cosine tables now initialized once for all in Cyrano
c     FFT:
c     Convention: direct transform over l-m with exp[-i(l-m)khi]
      
      call zset((maxcou+1)*(npfft+1), czero, cr, 1)
      call zset((2*maxcou+1)*(2*npfft+1), czero, gcdr1, 1)

        if(updsym)then
c       Using up-down symmetry:
        npp = npft2 + 1
        else
        npp = npfft + 1
        end if

      if(method.eq.1)then
c     -------------------
c     Using FFT over khi.
c
c     In thbmkh, store (Thbar-khi) on equidistant khi mesh: 
c     (contents of existing vector polang are used for equidistant abscissae)
c     polang, khi and Thbar all range from 0 to 2pi; 
c     NB: Thbar was tabulated vs theta in ckt(,,3), khi in ckt(,,5).
      call interp2(ckt(intab,1,5), nabplo, ckt(intab,1,3), nabplo, npp
     ;             ,polang, thbmkh, npp)

      thbmkh(1:npp) = thbmkh(1:npp) - polang(1:npp) 
      rnfi = 1.d0 / dfloat(npfft)

c     Get dkhi/dtheta by Fourier tf:

	cr(1:npfft,0) = (ckt(intab,1:npfft,5) - polang(1:npfft)) * rnfi

cERN		call cfft2(0, -1, npfft, cr(1,0), cworkm, cr(1,0))
          call df2tcf(npfft, cr(1,0), cr(1,0), cwork2, cpy) ! -- IMSL --
	    do i = 1, npft2
	       cr(i,0) = dcmplx(0.d0,dfloat(i-1)) * cr(i,0)
	       cr(npfft+1-i,0) = dcmplx(0.d0,dfloat(-i)) * cr(npfft+1-i,0)
          end do
	    
cERN		call cfft2(0, 1, npfft, cr(1,0), cworkm, cr(1,0))
	    call df2tcb(npfft, cr(1,0), cr(1,0), cwork2, cpy) ! -- IMSL --
            dkhdth(1:npfft) = dreal(cr(1:npfft,0)) + 1.d0
c	Impose last point (2*pi) is equal to first point (0)
		  dkhdth(npfft+1) = dkhdth(1) 


c     Table of dThetabar/dkhi on standard theta mesh:
cERN      dthdkh(1) = (ckt(intab,2,3)-ckt(intab,1,3)) /
cERN     ;            (ckt(intab,2,5)-ckt(intab,1,5))
cERN        do ipo = 2, npp-1
cERN        dthdkh(ipo) =
cERN     ;     (ckt(intab,ipo+1,3)-ckt(intab,ipo-1,3)) /
cERN     ;     (ckt(intab,ipo+1,5)-ckt(intab,ipo-1,5))
cERN        end do
cERN      dthdkh(npp) = 
cERN     ;     (ckt(intab,npp,3)-ckt(intab,npp-1,3)) /
cERN     ;     (ckt(intab,npp,5)-ckt(intab,npp-1,5))

	dthdkh(1:npp) = ckt(intab,1:npp,1) / dkhdth(1:npp)

c     Additional factor for Maxwellian case:
      if(icase.eq.1)then
	   dthdkh(1:npp) = dthdkh(1:npp)*b0/bmotab(intab,1:npp)
      end if

c     Interpolate to uniform khi mesh:
      call interp2(ckt(intab,1,5), nabplo, dthdkh, 1, npp
     ;, polang, dtdkvk, npp)

c     In cr(*,k), store exp{ik(Thbar-khi)} dThbar/dkhi / npfft
c     on equidistant khi mesh:
	   cr(1:npp,0) = dcmplx(dtdkvk(1:npp)*rnfi, 0.d0)

c       Exponential of i*(thbar-khi):
 	   ei(1:npp) = dcmplx(dcos(thbmkh(1:npp)), dsin(thbmkh(1:npp)))
         cr(1:npp,1) = ei(1:npp) * cr(1:npp,0)   

c       Recurrence using properties of trigonometric functions:
        do ka = 2, klim
	     cr(1:npp,ka) = (2.d0*dreal(ei(1:npp)))*cr(1:npp,ka-1) - cr(1:npp,ka-2)
        end do
            
      if(updsym)then
c     Using up-down symmetry:
          do ipo = 1, npp-1
             cr(npp+ipo,0:klim) = dconjg(cr(npp-ipo,0:klim))
          end do
      end if

      if(wrifou)then
	open(unit=27, file='..\..\integrand_ka0', status='unknown')
        do L = 1,npfft+1
        write(27,*)L, dreal(cr(L,0)), dimag(cr(L,0))
        end do
      end if

c     Call fft, n=2**m case: voir cas sym, antisym
c     --------------------------------------------
c     klim vs npfft: ?
      do ka = 0, klim

ccccccccccccccccccccccccc ERNESTO ccccccccccccccccccccccccccccc

cERN 		 call cfft2(0, -1, npfft, cr(1,ka), cworkm, cr(1,ka))
	     call df2tcf(npfft, cr(1,ka), cr(1,ka), cwork2, cpy) ! -- IMSL --
c          To compare with Num.Recipes routine:
c          call four1(cr(1,ka), npfft, -1)
c          On output: npfft coefficients in cr(1,ka). 
c                Indices 1 to npfft/2+1 correspond to harmonics 0 to npfft/2
c                Indices npfft/2+1 to npfft correspond to harmonics -npfft/2 to -1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      if(wrifou)then
      write(6,*)'Geometrical Fourier coefficients Bmj: ka=',ka
        do L = 1,npfft
        write(6,*)L, dreal(cr(L,ka)), dimag(cr(L,ka))
        end do
      end if

c     Sort harmonics in conventional order:
      if(.not.updsym)then
	l = 1 - ka
	if(l.le.0)l = l + npfft
      gcdr(0,ka) = cr(l,ka)
        do i = 1, npfft/2-ka
	  l = i + 1 - ka
	  if(l.le.0)l = l + npfft
        gcdr(i,ka) = cr(l,ka) + cr(npfft+1-i-ka,ka)
        gcdr(-i,ka) = gcdr(i,ka)
        end do
      else
c     use up-down symmetry: the spectrum is real.
	l = 1 - ka
	if(l.le.0)l = l + npfft
      gcdr(0,ka) = dreal(cr(l,ka))
        do i = 1, npfft/2-ka
	  l = i + 1 - ka
	  if(l.le.0)l = l + npfft
        gcdr(i,ka) = dreal(cr(l,ka)) + dreal(cr(npfft+1-i-ka,ka))
        gcdr(-i,ka) = gcdr(i,ka)
        end do
      end if

c     Pad the remaining of gcdr with zeros:
        do i = npfft/2-ka+1, npfft
	  gcdr(i,ka) = czero
	  gcdr(-i,ka) = czero
        end do
      end do

      else if(method.eq.2)then
c     -------------------
c     Using FFT over Thetabar.
c
c     In khothb, store khi on equidistant Thetabar mesh: use contents
c     of existing vector polang:
c     polang, khi and Thbar all range from 0 to 2pi; 
c     NB: Thbar was tabulated in ckt(,,3), khi in ckt(,,5).
      call interp2(ckt(intab,1,3), nabplo, ckt(intab,1,5), nabplo, npp
     ;, polang, khothb, npp)

c     For Maxwellian case: different geom. coeffs, B0 required:
      if(icase.eq.1)
     ;call interp2(ckt(intab,1,3), nabplo, bmotab(intab,1), nabplo, npp
     ;, polang, bvsthb, npp)

c     In cr, first store cos(l*khi) / (npfft)
c     on equidistant thetabar mesh:
      rnfi = 1.d0 / dfloat(npfft)
        do ipo = 1, npp
        cr(ipo,0) = rnfi
        end do

c       Additional factor for Maxwellian case:
        if(icase.eq.1)then
	    do ipo = 1, npp
	    cr(ipo,0) = cr(ipo,0) * b0 / bvsthb(ipo)
	    end do
	  end if

        do ipo = 1, npp
c       cos(khi):
        ei(ipo) = dcmplx(dcos(khothb(ipo)),0.d0)
        cr(ipo,1) = ei(ipo) * cr(ipo,0)   
        end do
c       Recurrence using properties of trigonometric functions:
        do ka = 2, klim
            do ipo = 1, npp
            cr(ipo,ka) = (2.d0*dreal(ei(ipo)))*cr(ipo,ka-1) - cr(ipo,ka-2)
            end do
        end do
            
      if(updsym)then
c     Using up-down symmetry:
        do ka = 0, klim
          do ipo = 1, npp-1
          cr(npp+ipo,ka) = dconjg(cr(npp-ipo,ka))
          end do
        end do
      end if

cccccccccccccccccccc  ERNESTO cccccccccccccccccc
	print *, 'method=2 not working yet!'
	print *, 'STOP at FOUGDR Line 303'
	stop 
cccccccccccccccccccccccccccccccccccccccccccccc

c  AT WORK HERE: method=2 not working yet!

c     CALL CRAY FAST FOURIER TRANSFORM cfft2, N=2**M CASE: VOIR CAS SYM, ANTISYM
c     --------------------------------------------------------------------------
c     klim vs npfft: ?
      do ka = 0, klim
cERN      call cfft2(0, -1, npfft, cr(1,ka), cworkm, cr(1,ka))
	call df2tcf(npfft, cr(1,ka), cr(1,ka), cwork2, cpy) ! -- IMSL --
c     On output: npfft coefficients in cr(1,ka). 
c                Indices 1 to npfft/2+1 correspond to harmonics 0 to npfft/2
c                Indices npfft/2+1 to npfft correspond to harmonics -npfft/2 to -1
c     Sort harmonics in conventional order:
      if(.not.updsym)then
	l = 1 - ka
	if(l.le.0)l = l + npfft
      gcdr(0,ka) = cr(l,ka)
        do i = 1, npfft/2-ka
	  l = i + 1 - ka
	  if(l.le.0)l = l + npfft
        gcdr(i,ka) = cr(l,ka) + cr(npfft+1-i-ka,ka)
        gcdr(-i,ka) = gcdr(i,ka)
        end do
      else
c     use up-down symmetry: the spectrum is real.
	l = 1 - ka
	if(l.le.0)l = l + npfft
      gcdr(0,ka) = dreal(cr(l,ka))
        do i = 1, npfft/2-ka
	  l = i + 1 - ka
	  if(l.le.0)l = l + npfft
        gcdr(i,ka) = dreal(cr(l,ka)) + dreal(cr(npfft+1-i-ka,ka))
        gcdr(-i,ka) = gcdr(i,ka)
        end do
      end if

c     Pad the remaining of gcdr with zeros:
        do i = npfft/2-ka+1, npfft
	  gcdr(i,ka) = czero
	  gcdr(-i,ka) = czero
        end do
      end do

      end if
c     ------

      do ka = 1, klim
c     Build the Pml coefficients for m<0
        do i = -npfft, npfft
        gcdr(i,-ka) = dconjg(gcdr(-i,ka))
        end do
      end do
c     Now gcdr(L,m) contains PmL
c     for -npfft/2 <= L-k <= npfft/2 (zeros outside this range),
c     and -klim <= m <= klim

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cERN 14Feb2005: Correction FACTOR !!!!!!
c	We know that gcdr(L,m=0) must be 1 (ONE).
c	CORR = 1.0d0 / gcdr(0,0)
c	gcdr = 1.0d0 / dreal(gcdr(0,0)) * gcdr


cERN  ccccccccccccccc Writing OUTPUT file ccccccccccccccccccccccccc
      
	if (WRITE_TABLES) then

		open (UNIT = 777, FILE = TABFOLDER // '/Plm.dat', 
     ;		  STATUS = "REPLACE", ACTION = "WRITE")
		write(777,*)'Geometrical Fourier coef. for M12 response: rho(m) = ', 
     ;                 abscis(intab)
		do ka = 0, klim
		   write(777,*)'ka=',ka
		   do L = -npfft/2,npfft/2
			  write(777,*)L, dreal(gcdr(L,ka)), dimag(gcdr(L,ka))
		   end do
		end do
c		write(6,1000)
	
	else
	
		if (WRITE_OUTPUT) then
			open (UNIT = 777, FILE = COKFOLDER // '/Plm.dat', 
     ;			  STATUS = "REPLACE", ACTION = "WRITE")
			write(777,*)'Geometrical Fourier coef. for M12 response: rho(m) =' , 
     ;                 abscis(intab)
			do ka = 0, klim
			   write(777,*)'ka=',ka
			   do L = -npfft/2,npfft/2
				  write(777,*)L, dreal(gcdr(L,ka)), dimag(gcdr(L,ka))
			   end do
			end do
c			write(6,1000)
		end if
		close(777)

	end if ! WRITE_TABLES

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
	return

 1000 format(1h ,i2,2(2x,g12.4))
      end
