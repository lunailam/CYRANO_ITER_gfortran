      SUBROUTINE oi_trap_react_cplx_anal(xin, vin, sigm, mav, Nmdiff, qovm, 
     ;                                   testanal, FLAG)

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
	real*8, intent(out) :: testanal(1:Nmdiff) 
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
	real*8 :: thref, thtilde, coefA(1:2), coefB(1:2), soma_anal(1:Nmdiff)
	real*8 :: meuteste(1000), fff

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
	co = eqt(intab,1,14)

 	testanal(1:Nmdiff) = 0.d0
                        

	theres = 0.d0
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




	soma_anal(1:Nmdiff) = 0.d0
	meuteste(1)=0.d0

	do j = 1, Nmdiff

	   mdiff = dfloat(j - 1)
         
	   do ipo = 2, NTT
              	

cNEW  Analytical (expansion around thref) ----------------------------
	if(th(ipo).ne.th(ipo-1))then 
	           thref = 0.5d0*(th(ipo)+th(ipo-1))
			   call expand_A(thref, mdiff, coefA)  ! numerator
			   call expand_B(thref, coefB)  ! denominator

				meuteste(ipo) = ( coefA(1) + coefA(2) * (th(ipo)-thref) ) 
     ;                          / ( coefB(1) + coefB(2) * (th(ipo)-thref) )   
				

c			   Apply correction factor
			   fff = 0.5d0/sintt

   
		if(coefB(2).ne.0.d0)then	
		        thtilde = thref - coefB(1)/coefB(2);	
				soma_anal(j) = soma_anal(j)
     ;            +  (coefA(1)/coefB(2)-coefA(2)*coefB(1)/coefB(2)**2) 
     ;			    * fff * ( dlog(dabs(th(ipo)-thtilde)) - dlog(dabs(th(ipo-1)-thtilde)) )
     ;			+   coefA(2)/coefB(2) * fff *(th(ipo)-th(ipo-1)) 
		end if
	
	else

				meuteste(ipo) = meuteste(ipo-1) 

	endif
c			-------------------------------------------------------
  
         end do  

	end do

		  testanal(1:Nmdiff) = soma_anal(1:Nmdiff)
	        meuteste=fff*meuteste

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	meuteste(1)=meuteste(2)

      if(FLAG .and. vin.gt.0.9d6 .and. vin<1.1d6)then
c      if(FLAG)then
         write(1559,"(1000G16.8)"), xin, xnout, NTT, 0.d0, 0.d0, meuteste(1:NTT)
         write(1559,"(1000G16.8)"), xin, xnout, NTT, 0.d0, 0.d0, th(1:NTT)
      end if


      return
 
      END SUBROUTINE oi_trap_react_cplx_anal