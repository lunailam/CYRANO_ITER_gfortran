      SUBROUTINE oi_pass_react(xin, vin, sigm, mav, Nmdiff, qovm, xtg, xo,
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

	real*8 f_ic_anal
	external f_ic_anal

	real*8 f_vp_std
	external f_vp_std

c	integer itmax2, ilast_R, ilast_L, Nins

      real*8 :: teste, thaux, sintt, btt, ini, fim
	real*8 :: sai1, sai2, limit, delti, xnout
	real*8 :: gammadd, theres, mdiff, kp, kpp, deltwa 
	real*8 :: aux(2000), soma2(20),th(2000)
	integer NTT, range, auxi, iL, iR, ind2

c	For analytical integral
c	real*8 :: thref, thtilde, coefA(1:2), coefB(1:2), soma_anal(1:Nmdiff), fff
c	real*8 :: meuteste(1000)

	integer isrchfge
	external isrchfge



c     Set global variables (for f_reso.f routine)
      VNOW = vin
	XNOW = xin
	sigNOW = sigm
	qomNOW = qovm
c	kpNOW = kp    !!!!!!!!!!!!!!!!!!!!!!
	mavNOW = mav  !!!!!!!!!!!!!!!!!!!!!!


      xmax = B0 / bmin(intab)
	xnout = -sigm*(1.d0-xin/xmax)
	sintt= eqt(intab,1,12)


	residue(1:Nmdiff)=0.d0




ccccccccccccc Residue at theres cccccccccccccccccccccccccccccccccccccccc


	theres = thein

	
	
c	Construct theta grid for Principal Value integration
c     [0, ... , theres-delta , theres+delta , ... , pi]	
	NTT=100;
c	deltwa=1.d-3
	do ipo=1, NTT
	   th(ipo) = dfloat(ipo-1)/dfloat(NTT-1) * pi
c	   if( dabs(th(ipo) - theres)  < deltwa)then
c            th(ipo) = theres + dsign(1.d0,th(ipo)-theres) * deltwa
c	   end if
	end do

c	Correct points near theres
	deltwa=2.d-3
	auxi = isrchfge(NTT, th, 1, theres)

	iL = auxi-1
	iR = auxi

	if(auxi<2 .or. auxi>NTT)goto 1234

	   if(auxi.eq.2)then  ! special case for theres near 0

c		insert 1 point between 0 and theres (ipo=2)
c		th(iL=1)=0 ! don't touch
		if(deltwa<theres)then
			th(3:NTT+1)=th(2:NTT)  ! shift original values
			th(2)= theres-deltwa   ! insert point
			NTT=NTT+1
			iL=2
			iR=3
		end if
		
c		correct right side
		th(iR)  = theres+deltwa

	   elseif(auxi.eq.NTT)then ! special case for theres near pi


c		insert 1 point before last one 'pi'  (ipo=NTT)
c		th(iL=1)=0 ! don't touch
		if(deltwa<(pi-theres))then
			th(NTT+1)=th(NTT)  ! shift original value (last point only)
			th(NTT)= theres+deltwa   ! insert point
			NTT=NTT+1
			iL=NTT-2
			iR=NTT-1
		end if

c		USE same delta:
c		correct left side
		th(iL)  = theres-dmin1(deltwa,th(NTT)-theres)

	   else  ! (standard procedure - shift closest points)

c			iL = auxi-1
c			iR = auxi			
			th(iL)= theres-deltwa
	            th(iR)  = theres+deltwa


	   end if



c	refine grid near theres


	range=min(4,min(4,iL),min(4,NTT-iR))

c	LEFT
c	  range = 4
	 ind2 = max(iL-range, 1)
c	  range=iL-ind2
	  call insert_points_new(NTT, th(1:NTT), ind2, iL, 20, th(1:NTT+20))
        NTT = NTT + 20
	  iL=iL+20; iR=iR+20 !!!!!!!!

c	RIGHT
	range = min(4,range)  ! use same range as before if needed
	  ind2 = min(iR+range, NTT)
	  call insert_points_new(NTT, th(1:NTT), iR, ind2, 20, th(1:NTT+20))
        NTT = NTT + 20




1234	continue



c	refine grid near pi (if not refined by theres yet)
c	only LEFT
c	range = 4
c	if(NTT-iR>range)
c	  ind2 = NTT-range
c	  call insert_points_new(NTT, th(1:NTT), ind2, NTT, 20, th(1:NTT+20))
c       NTT = NTT + 20
c	end if

      if(iR.eq.0)then
	range = 4
         range=min(abs(NTT-iR),4)
	   ind2 = NTT-range
	   call insert_points_new(NTT, th(1:NTT), ind2, NTT, 20, th(1:NTT+20))
         NTT = NTT + 20
	call insert_points_new(NTT, th(1:NTT), NTT-1, NTT, 8, th(1:NTT+8))
         NTT = NTT + 8
      end if


      if(iL.eq.0)then
	range = 4
c         range=min(iL,4)
	   ind2 = 1+range
	   call insert_points_new(NTT, th(1:NTT), 1, ind2, 20, th(1:NTT+20))
         NTT = NTT + 20
	call insert_points_new(NTT, th(1:NTT), 1, 2, 8, th(1:NTT+8))
         NTT = NTT + 8
      end if


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
     ;            + 0.5d0 * ( aux(ipo)*dcos(mdiff*th(ipo)) + aux(ipo-1)*dcos(mdiff*th(ipo-1)) )
     ;                    * ( th(ipo) - th(ipo-1) )
			end if			 

c      goto 8888  
cNEW  Analytical (expansion around thref) ----------------------------
c	if(th(ipo).ne.th(ipo-1))then 
c	           thref = 0.5d0*(th(ipo)+th(ipo-1))
c			   call expand_A(thref, mdiff, coefA)  ! numerator
c			   call expand_B(thref, coefB)  ! denominator

c				meuteste(ipo) = ( coefA(1) + coefA(2) * (th(ipo)-thref) ) 
c     ;                          / ( coefB(1) + coefB(2) * (th(ipo)-thref) )   
				

cc			   Apply correction factor
c			   fff = 0.5d0/sintt

   
c		if(coefB(2).ne.0.d0)then	
c		        thtilde = thref - coefB(1)/coefB(2);	
c				soma_anal(j) = soma_anal(j)
c     ;            +  (coefA(1)/coefB(2)-coefA(2)*coefB(1)/coefB(2)**2) 
c     ;			    * fff * ( dlog(dabs(th(ipo)-thtilde)) - dlog(dabs(th(ipo-1)-thtilde)) )
c     ;			+   coefA(2)/coefB(2) * fff *(th(ipo)-th(ipo-1)) 
c		end if
	
c	else

c				meuteste(ipo) = meuteste(ipo-1) 

c	endif
c			-------------------------------------------------------
c8888  continue

         end do  

	end do


              residue(1:Nmdiff) = soma2(1:Nmdiff)
c			testanal(1:Nmdiff) = soma_anal(1:Nmdiff)
c	        meuteste=fff*meuteste

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	meuteste(1)=meuteste(2)


      if(FLAG .and. vin.gt.0.9d6 .and. vin<1.1d6)then
c     if(FLAG)then
         write(1552,"(1000G16.8)"), xin, xnout, NTT, theres, 0.d0, aux(1:NTT)
         write(1552,"(1000G16.8)"), xin, xnout, NTT, theres, 0.d0, th(1:NTT)
      end if


      return
 
      END SUBROUTINE oi_pass_react