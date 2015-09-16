      SUBROUTINE oi_trap_react_double(xin, vin, sigm, mav, Nmdiff, qovm, xtg, xtg2,xo,
     ;                                     thein, thein2, residue, FLAG)

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
	real*8, intent(in)  :: xin, vin, mav, qovm , sigm, xtg, xtg2,xo, thein, thein2
	real*8, intent(out) :: residue(1:Nmdiff)
	logical, intent(in) :: FLAG
c      real*8 sigm
     
!     Variables
      integer j, ipo, itmax
	real*8 :: xmax, xsep, vpar2, vpar, vperp2

c	For analytical integral
	real*8 :: thref, thtilde, coefA(1:2), coefB(1:2), soma_anal(1:Nmdiff)
	real*8 :: meuteste(1000), fff

c	real*8 f_reso_std
c	external f_reso_std

	real*8 f_ic_anal
	external f_ic_anal

	real*8 f_vp_std
	external f_vp_std

c	integer itmax2, ilast_R, ilast_L, Nins

      real*8 :: teste, thaux, sintt, btt, ini, fim
	real*8 :: sai1, sai2, limit, delti, xnout
	real*8 :: gammadd, theres, mdiff, kp, kpp, deltwa, theres2 
	real*8 :: aux(2000), soma2(20),th(2000), themax
	integer NTT, range,auxi, iL, iR, ind2, auxi2

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

c	k//(theta) passed through COMMON kptab(1:npfft+1)  
c	          (DIRESP_standard_1 line 489)



	residue(1:Nmdiff)=0.d0
c 	testanal(1:Nmdiff) = 0.d0


	sintt= eqt(intab,1,12)
	co = eqt(intab,1,14)
	themax = dacos(r0/abscis(intab)*(xin/co-1.d0))

	theres = dmin1(thein,thein2)  ! sort the roots
	theres2= dmax1(thein,thein2)		






ccccccccccccc Residue at theres cccccccccccccccccccccccccccccccccccccccc

	
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


	auxi = isrchfge(NTT, th, 1, theres)
	auxi2 = isrchfge(NTT, th, 1, theres2)

ccccccccccccccccccccccccccccccccccccc
c	if(auxi.eq.auxi2)then  ! too close roots
c		residue(1:Nmdiff)=0.d0
c		goto 9999
c	endif
ccccccccccccccccccccccccccccccccccccc




c	Correct points near theres
	iL = auxi-1
	iR = auxi
	deltwa=2.d-3

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

	   elseif(auxi.eq.NTT)then ! special case for theres near themax



c		insert 1 point before last one 'themax'  (ipo=NTT)
c		th(iL=1)=0 ! don't touch
		if(deltwa<(th(NTT)-theres))then
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


	range=min(4,min(4,iL),min(4,NTT-iR+1))

c	LEFT
c	  range = 4
	 ind2 = max(iL-range, 1)
c	  range=iL-ind2
	  call insert_points_new(NTT, th(1:NTT), ind2, iL, 20, th(1:NTT+20))
        NTT = NTT + 20
	  iL=iL+20; iR=iR+20 !!!!!!!!

c	RIGHT
c	range = min(4,range)  ! use same range as before if needed
	  ind2 = min(iR+range, NTT)
	  call insert_points_new(NTT, th(1:NTT), iR, ind2, 20, th(1:NTT+20))
        NTT = NTT + 20




1234	continue




c	Correct points near theres 2

	if(theres2<themax)then ! -----------------------------------

	deltwa=2.d-3
	
	auxi = isrchfge(NTT, th, 1, theres2)
	iL = auxi-1
	iR = auxi

	if(auxi<2 .or. auxi>NTT)goto 1235

	   if(auxi.eq.2)then  ! special case for theres2 near 0

c		insert 1 point between 0 and theres (ipo=2)
c		th(iL=1)=0 ! don't touch
		if(deltwa<theres2)then
			th(3:NTT+1)=th(2:NTT)  ! shift original values
			th(2)= theres2-deltwa   ! insert point
			NTT=NTT+1
			iL=2
			iR=3
		end if
		
c		correct right side
		th(iR)  = theres2+deltwa

	   elseif(auxi.eq.NTT)then ! special case for theres near themax

c		insert 1 point before last one 'themax'  (ipo=NTT)
c		th(iL=1)=0 ! don't touch
		if(deltwa<(th(NTT)-theres2))then
			th(NTT+1)=th(NTT)  ! shift original value (last point only)
			th(NTT)= theres2+deltwa   ! insert point
			NTT=NTT+1
			iL=NTT-2
			iR=NTT-1
		end if

c		USE same delta:
c		correct left side
		th(iL)  = theres2-dmin1(deltwa,th(NTT)-theres2)


	   else  ! (standard procedure - shift closest points)

c			iL = auxi-1
c			iR = auxi			
			th(iL)= theres2-deltwa
	        th(iR)  = theres2+deltwa


	   end if


c	refine grid near theres2


	range=min(4,min(4,iL),min(4,NTT-iR+1))

c	LEFT
c	  range = 4
	 ind2 = max(iL-range, 1)
c	  range=iL-ind2
	  call insert_points_new(NTT, th(1:NTT), ind2, iL, 20, th(1:NTT+20))
        NTT = NTT + 20
	  iL=iL+20; iR=iR+20 !!!!!!!!

c	RIGHT
c	range = min(4,range)  ! use same range as before if needed
	  ind2 = min(iR+range, NTT)
	  call insert_points_new(NTT, th(1:NTT), iR, ind2, 20, th(1:NTT+20))
        NTT = NTT + 20




1235	continue


	endif  ! (theres2<themax) -----------------------------------



c	refine grid near themax
c	only LEFT
	range = 4
c	if(iR<NTT)then
	range=min(abs(NTT-iR),4)
	   ind2 = NTT-range
	   call insert_points_new(NTT, th(1:NTT), ind2, NTT, 20, th(1:NTT+20))
         NTT = NTT + 20

	   call insert_points_new(NTT, th(1:NTT), NTT-1, NTT, 8, th(1:NTT+8))
         NTT = NTT + 8

c	end if





	do ipo=1,NTT

		vpar   = f_vp_std(th(ipo))
		vperp2 = vin*vin - vpar**2
		aux(ipo) = 0.5d0 * vperp2 / (vin*vin) / dabs(vpar) / sintt
     ;							  / f_ic_anal(th(ipo))

	end do


	soma2(1:Nmdiff) = 0.d0

	do j = 1, Nmdiff

	   mdiff = dfloat(j - 1)
         
	   do ipo = 2, NTT
              	
			if( ipo>1 .and. aux(ipo)*aux(ipo-1)>0 )then  ! exclude singularity
				soma2(j) = soma2(j)
     ;            + 0.5d0 * ( aux(ipo)*dcos(mdiff*th(ipo)) + aux(ipo-1)*dcos(mdiff*th(ipo-1)) )
     ;                    * ( th(ipo) - th(ipo-1) )
			end if			 

         end do  

	end do


              residue(1:Nmdiff) = soma2(1:Nmdiff)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

9999	continue

	meuteste(1)=meuteste(2)

      if(FLAG .and. vin.gt.0.9d6 .and. vin<1.1d6)then
c      if(FLAG)then
         write(1552,"(1000G16.8)"), xin, xnout, NTT, theres, theres2, aux(1:NTT)
         write(1552,"(1000G16.8)"), xin, xnout, NTT, theres, theres2, th(1:NTT)

      end if


      return
 
      END SUBROUTINE oi_trap_react_double