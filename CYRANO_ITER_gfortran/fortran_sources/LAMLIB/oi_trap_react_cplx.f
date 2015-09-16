      SUBROUTINE oi_trap_react_cplx(xin, vin, sigm, mav, Nmdiff, qovm, 
     ;                              residue,  FLAG)

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
	real*8, intent(in)  :: xin, vin, mav, qovm , sigm
	real*8, intent(out) :: residue(1:Nmdiff) 
	logical, intent(in) :: FLAG
c      real*8 sigm
     
!     Variables
      integer j, ipo, itmax
	real*8 :: xmax, xsep, vpar2, vpar, vperp2

c	real*8 f_reso_std
c	external f_reso_std

	real*8 f_ic_anal
	external f_ic_anal

	real*8 f_vp_std
	external f_vp_std

c	integer itmax2, ilast_R, ilast_L, Nins

      real*8 :: teste, thaux, sintt, btt, ini, fim
	real*8 :: sai1, sai2, limit, delti, xnout
	real*8 :: gammadd, theres, mdiff, kp, kpp, deltwa 
	real*8 :: aux(2000), soma2(20),th(2000), themax
	integer NTT, range,auxi, iL, iR, ind2

c	For analytical integral
c	real*8 :: thref, thtilde, coefA(1:2), coefB(1:2), soma_anal(1:Nmdiff)
c	real*8 :: meuteste(1000), fff

	integer isrchfge
	external isrchfge


c     Set global variables (for f_reso.f routine)
      VNOW = vin
	XNOW = xin
	sigNOW = sigm
	qomNOW = qovm
	mavNOW = mav  !!!!!!!!!!!!!!!!!!!!!!


      xmax = B0 / bmin(intab)
	xnout = -sigm*(1.d0-xin/xmax)
	sintt= eqt(intab,1,12)



	residue(1:Nmdiff)=0.d0




	theres = 0.d0

	co = eqt(intab,1,14)
	themax = dacos(r0/abscis(intab)*(xin/co-1.d0))
	
c	Construct theta grid for Principal Value integration
c     [0, ... , theres-delta , theres+delta , ... , theta_max]	
	NTT=100;
	deltwa=1.d-3
	do ipo=1, NTT
	   th(ipo) = dfloat(ipo-1)/dfloat(NTT-1) * themax  !! NOT pi as for pass
	   if( dabs(th(ipo) - themax)  < deltwa)then
            th(ipo) = themax - deltwa
	   end if
	end do




c	refine grid near themax
c	only LEFT
c	only LEFT
	range = 4
	   ind2 = min(NTT,NTT-range)
	   call insert_points_new(NTT, th(1:NTT), ind2, NTT, 20, th(1:NTT+20))
         NTT = NTT + 20
	   call insert_points_new(NTT, th(1:NTT), NTT-1, NTT, 8, th(1:NTT+8))
         NTT = NTT + 8






	do ipo=1,NTT

		vpar   = f_vp_std(th(ipo))
		vperp2 = vin*vin - vpar**2
		aux(ipo) = 0.5d0 * vperp2 / (vin*vin) / abs(vpar) / sintt
     ;							  / f_ic_anal(th(ipo))

	end do


	soma2(1:Nmdiff) = 0.d0
c	soma_anal(1:Nmdiff) = 0.d0
c	meuteste(1)=0.d0

	do j = 1, Nmdiff

	   mdiff = dfloat(j - 1)
         
	   do ipo = 2, NTT
              	
			if( ipo>1 .and. aux(ipo)*aux(ipo-1)>0 )then  ! exclude singularity
				soma2(j) = soma2(j)
     ;            + 0.5d0 * ( aux(ipo)*cos(mdiff*th(ipo)) + aux(ipo-1)*cos(mdiff*th(ipo-1)) )
     ;                    * ( th(ipo) - th(ipo-1) )
			end if			 

  
         end do  

	end do


              residue(1:Nmdiff) = soma2(1:Nmdiff)
c			testanal(1:Nmdiff) = soma_anal(1:Nmdiff)
c	        meuteste=fff*meuteste

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c	meuteste(1)=meuteste(2)

      if(FLAG .and. vin.gt.0.9d6 .and. vin<1.1d6)then
c     if(FLAG)then
         write(1552,"(1000G16.8)"), xin, xnout, NTT, 0.d0, 0.d0, aux(1:NTT)
         write(1552,"(1000G16.8)"), xin, xnout, NTT, 0.d0, 0.d0, th(1:NTT)

      end if


      return
 
      END SUBROUTINE oi_trap_react_cplx