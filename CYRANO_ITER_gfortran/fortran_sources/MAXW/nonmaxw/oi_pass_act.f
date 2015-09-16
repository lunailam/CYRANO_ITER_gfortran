      SUBROUTINE oi_pass_act(xin, vin, sigm, mav, Nmdiff, qovm, xtg, xo,
     ;                                     thein, residue, FLAG)

      IMPLICIT NONE


c     Compute orbit integral for given (x,v,k//,mdiff)
																	
c	PASSING branch (resonant + non-resonant)

c	######### VERSION for STANDARD coordinates ###############

																		


      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'cokpco.copy'
      include 'comgeo.copy'
      include 'compla.copy'
      include 'comma2.copy'
      include 'commod.copy'
	 
!	Input	                                                                 
	integer, intent(in) :: Nmdiff
	real*8, intent(in)  :: xin, vin, mav, qovm , sigm, xtg, xo, thein
	real*8, intent(out) :: residue(1:Nmdiff) 
	logical, intent(in) :: FLAG
c      real*8 sigm
     
!     Variables
      integer j, ipo, itmax
	real*8 :: xmax, xsep, vpar2, vpar, vperp2


c	real*8 f_reso_std
c	external f_reso_std

	real*8 f_ic_std
	external f_ic_std

	real*8 f_vp_std
	external f_vp_std

c	integer itmax2, ilast_R, ilast_L, Nins

      real*8 :: teste, thaux, sintt, btt, ini, fim
	real*8 :: sai1, sai2, limit, delti, xnout
	real*8 :: gammadd, theres, mdiff, kp, kpp, deltwa 


c	integer ISRCHFGE
c	external ISRCHFGE



c     Set global variables (for f_reso.f routine)
      VNOW = vin
	XNOW = xin
	sigNOW = sigm
	qomNOW = qovm
c	kpNOW = kp    !!!!!!!!!!!!!!!!!!!!!!
	mavNOW = mav  !!!!!!!!!!!!!!!!!!!!!!


      xmax = B0 / bmin(intab)
	xnout = -sigm*(1.d0-xin/xmax)

c	k//(theta) passed through COMMON kptab(1:npfft+1)  
c	          (DIRESP_standard_1 line 489)



	residue(1:Nmdiff)=0.d0
                  
c           itmax  = npfft/2+1
c           itmax2 = itmax
                  



ccccccccccccc Residue at theres cccccccccccccccccccccccccccccccccccccccc



	    theres = thein

c     d/dt(gammadot) = |v//|.sin(THETA)/r * d/dtheta(gammadot)
c     where gammadot = w-wc-k//v//
	if(xin>xtg)then


			vpar   = f_vp_std(theres)
			vperp2 = vin*vin - vpar**2



c	Find sin(THETA) form tables at polang=theres	
c      call spline0(3, polang(ipo-1:ipo+1),eqt(intab,ipo-1:ipo+1,12),
c     ;             1, theres, sintt)

	sintt= eqt(intab,1,12)

c	gamma_dot_dot
	deltwa=0.3d-3
      gammadd = dabs(vpar)*sintt * (f_ic_std(theres+deltwa) - f_ic_std(theres-deltwa)) 
     ;                     / (2*deltwa)

c	(Nmdiff = klim + 1)  
	do j = 1, Nmdiff
	   mdiff = dfloat(j - 1)  
	   residue(j) = 0.5d0 * pi * dcos(mdiff*theres) * vperp2 / (vin*vin) / gammadd
	end do

      end if



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c      if(FLAG .and. vin.gt.0.9d6 .and. vin<1.1d6)then
      if(FLAG)then
c     Write data to file (only for debugging) ----------------------

c                write(1552,"(200G16.8)"), xin, xnout, npfft/2+1, polang(itmax), (aux(1:npfft/2+1))
c                write(1552,"(200G16.8)"), xin, xnout, npfft/2+1, polang(itmax), polang(1:npfft/2+1)
c                write(1552,"(200G16.8)"), xin, xnout, itmax2, th(itmax2), (aux2(1:itmax2))
c                write(1552,"(200G16.8)"), xin, xnout, sai2, sai1, th(1:itmax2)

c               write(1553,"(400G16.8)")
               write(1553,"(400G16.8)"), xin, xnout, xtg, pi, theres, 0.d0
c			 write(1553,"(400G16.8)"), theres, 0, 0, 0, 0, residue 
c		     write(1553,"(400G16.8)"), 0, 0, 0, 0, 0, 0


      end if


      return
 
      END SUBROUTINE oi_pass_act