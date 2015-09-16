      SUBROUTINE orbit_integ_trap_std(xin, vin, sigm, mav, Nmdiff, qovm,  
     ;           xtg, xo, resultat, res2, resid, FLAG)

      IMPLICIT NONE

c     Compute orbit integral for given (x,v,k//,mdiff)
																	
c	TRAPPED real branch (resonant + non-resonant)

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
	real*8, intent(in)  :: xin, vin, mav, qovm , xo, sigm, xtg
	real*8, intent(out) :: resultat, res2(1:Nmdiff), resid(1:Nmdiff)
	logical, intent(in) :: FLAG
c      real*8 sigm
     
!     Variables
      integer j, ipo, itmax
   	real*8 :: xmax, xsep, vpar2, vpar, vperp2
	real*8 :: wc, aux(npfft+1), aux2(npfft+1+400), soma, th(npfft+1+400), 
     ;          soma2(1:Nmdiff)

	real*8 f_reso_std
	external f_reso_std

	real*8 f_ic_std
	external f_ic_std

	integer  ind2, itmax2, ilast_R, ilast_L
     ;        , ilast_R2, ilast_L2
c	real*8 epsabs, epsrel, abserr, alist(200), blist(200), vlast_R, vlast_L
c     ;     , rlist(200), elist(200)  
      integer Nins

      real*8 :: teste, thaux, sintt, btt, ini, fim, kp, kpp
	real*8 :: sai1, limit, urtimo, delti, xnout,  
     ;          themax, theres, theins, theres2, deltwa, mdiff

      real*8 :: gammadd, residue(1:Nmdiff), residue2(1:Nmdiff)


	real*8 :: sai1L, sai1R, sai2L, sai2R, pro1, pro2, proi1, proi2 
	integer :: achei1, achei2
	real*8 :: qua1, qua2, fac, g1, g2

      logical insertpts, insertpts2, tricky
	character(6) :: direct, tip1, tip2

	integer ISRCHFGE
	external ISRCHFGE

	integer find_nearest
	external find_nearest


c     Set global variables (for f_reso.f routine)
      VNOW = vin
	XNOW = xin
	sigNOW = sigm
	qomNOW = qovm
c	kpNOW = kp    !!!!!!!!!!!!!!!!!!!!!!
	mavNOW = mav  !!!!!!!!!!!!!!!!!!!!!!

c	k//(theta) passed through COMMON kptab(1:npfft+1)  
c	          (DIRESP_standard_1 line 489)

c	Radius dependent quantities
      xnsep = - 2.d0 * delb(intab) / bmax(intab) 
      xmax = B0 / bmin(intab)
      xsep = (1.d0-dabs(xnsep))*xmax
	xnout = -sigm*(1.d0-xin/xmax)

	th(1:npfft/2+1) = polang(1:npfft/2+1)
      th(npfft/2+2:npfft) = pi

	aux(1:npfft/2+1)=0
	aux2(1:npfft/2+1+400)=0
	ilast_L=0
	ilast_R=0
	ilast_L2=0
	ilast_R2=0
	sai1=0
c	sai2=0
	themax=0.d0
	theres=0.d0
	theres2=0.d0
      theins=0.d0

	sai1L=0.d0
	sai1R=0.d0
	sai2L=0.d0
	sai2R=0.d0
	pro1=0.d0
	pro2=0.d0
	proi1=0.d0
	proi2=0.d0

      insertpts = .false.
      insertpts2 = .false.
      residue(1:Nmdiff) = 0.d0
      residue2(1:Nmdiff) = 0.d0

      tricky = .false.


	achei1=0
	achei2=0
	qua1=0
	qua2=0
      tip1='XXX'
	tip2='XXX'
	g1=0
	g2=0


c          Find itmax in fine table
           ind2 = isrchfge(1001, xtabchi, 1, xin)
	     !ind2 = find_nearest(1001, xtabchi, xin)
		 itmax = ithemaxchi(ind2)
	     itmax2 = itmax


c      if(FLAG .and. vin.gt.4000000 .and.abscis(intab)>0.0685)then
c		print*
c	end if

cccccccccccccccccccccccc Theta maximum (v// -> 0) cccccccccccccccccccccccccc

c	if(circ)then

c		themax = dacos ( r0/abscis(intab) * (xin/eqt(intab,1,14)-1.d0) );
c		 theins = themax
c	else

      if(itmax2>1)then  !!!! .and. itmax2.gt.ilast_R
                              
	       limit = 1.d-8 !vmin*vmin / (vin*vin)			  
             ini = th(itmax2-1)
             fim = dmin1(th(itmax2+9),pi)  ! +4
                            
             do k = 1,50

                thaux = 0.5d0*(ini+fim)
	          call spline0(6,polang(itmax2-1:itmax2+4),bmotab(intab,itmax2-1:itmax2+4),
     ;                       1, thaux, btt)
                teste  = (1.d0 - xin*btt/B0)
                ! Always descending
	          if(teste.gt.limit)then
	             ini = thaux
	          else
	             fim = thaux
	          end if
                if(abs(teste-limit)/limit<1.d-2)goto 1223

	       end do
	                      
1223  continue

c     ------ Theta maximum ---------
	      themax = thaux
	      theins = thaux			
c     ------------------------------			
					
      end if	


c	end if  ! if circ




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c	Loop over poloidal angle

	do ipo = 1, npfft/2+1  ! =================================================================
		
	   if(polang(ipo) > themax)then
            itmax = ipo-1
	      goto 7777
	   end if
			 
	   vpar2  = vin*vin*(1.d0 - xin*bmotab(intab,ipo)/B0)
c	   vpar2  = dmax1(vpar2,1.d0)
	   vperp2 = vin*vin - vpar2
	   vpar   = sigm*dsqrt(abs(vpar2))
	   wc     = qovm * bmotab(intab,ipo)
cccccccccccccc
	   kp = kptab(ipo)  
cccccccccccc

	   aux(ipo) = 0.5d0 * vperp2 / (vin*vin) / abs(vpar) / eqt(intab,ipo,12)
     ;              / (omegag - wc - kp*vpar)
   	   aux2(ipo)= aux(ipo)


c	   Detect IC resonance point	

         if(ipo>1 .and. aux2(ipo)*aux2(ipo-1)<0.d0 .and. ilast_R.eq.0) then ! 1st resonance ----------------------------

		  if(aux2(ipo)<aux2(ipo-1))then
		     direct='descen'
		  else
	           direct='ascend'
		  end if

c           FIRST resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	      ilast_R = ipo      ! first point to the RIGHT of resonance
	      ilast_L = ipo-1    ! last  point to the LEFT  of resonance
            insertpts = .true.

c           Accuracy
            limit = 1.d-8
            
c	      Search for resonance (w-wc-k//v// = limit ~ 0)
		ini = th(ipo-1)
		fim = dmin1(th(ipo+1),themax)

		do k=1,80
               thaux = 0.5d0*(ini+fim)
	         call spline0(3, polang(ipo-1:ipo+1), bmotab(intab,ipo-1:ipo+1),
     ;                      1, thaux, btt)
               call spline0(3,polang(ipo-1:ipo+1),kptab(ipo-1:ipo+1),
     ;                      1, thaux, kpp)

		   vpar2  = vin*vin*(1.d0 - xin*btt/B0)
	         vpar   = sigm*dsqrt(abs(vpar2))
		   wc     = qovm * btt
               teste  = omegag - wc - kpp*vpar

		   if(teste.gt.limit)then
                  if(direct.eq.'descen')ini = thaux
			if(direct.eq.'ascend')fim = thaux
	         else
	            if(direct.eq.'descen')fim = thaux
			if(direct.eq.'ascend')ini = thaux
	         end if

               if(abs(teste-limit)/limit<1.d-2)goto 1233

	      end do

1233		continue
	           
c           FIRST resonance angle ----- 
	            theres = thaux
c           ---------------------------
	if(themax-theres<1.d-6)goto 7777

		achei1 = 1
	qua1 = teste
	tip1=direct
	     
c           Correct two points near FIRST resonance
            deltwa = 0.002d0  ! wanted delta

            if(th(ilast_R).lt.themax)then ! OK - - - - - - - - - - - - - 	     

c                Correct ipo = ilast_R2 (RIGHT)  +++++++++
                 delti     = dmin1(deltwa, 0.5*(themax-theres))
	           th(ipo)   = theres + delti
		       aux2(ipo) = f_reso_std(theres + delti)
c                Correct ipo-1 = ilast_L2 (LEFT) +++++++++
                 delti       = dmin1(delti, theres)
	           th(ipo-1)   = theres - delti
	           aux2(ipo-1) = f_reso_std(theres - delti)

		  else  !! CAREFUL (ilast_R == itmax2) - - - - - - - - - - - - - - 

                  print*, 'WARNING: Resonance 1 at BANANA tip !!!!!!!'
c			stop
                  insertpts = .false.
			goto 7777
	
              end if  !! CAREFUL (ilast_R == itmax2) - - - - - - - - - - - - -


ccccccccccccc (1) Residue at theres cccccccccccccccccccccccccccccccccccccccc

c     d/dt(gammadot) = |v//|.sin(THETA)/r * d/dtheta(gammadot)
c     where gammadot = w-wc-k//v//
	if(xin>xtg .or. direct.eq.'ascend')then

      deltwa = dmin1(0.0003d0,themax-theres)
c		deltwa = 0.002d0



      call spline0(3, polang(ipo-1:ipo+1),eqt(intab,ipo-1:ipo+1,12),
     ;             1, theres, sintt)
      vperp2 = vin*vin - vpar2
      
      gammadd = dabs(vpar)*sintt * (f_ic_std(theres+deltwa) - f_ic_std(theres-deltwa)) 
     ;                     / (2.d0*deltwa)

c	Residue  
      if(direct.eq.'descen')fac = 1.d0
	if(direct.eq.'ascend')fac = -1.d0
	do j = 1, Nmdiff  ! (Nmdiff = klim + 1)
	   mdiff = dfloat(j - 1)  
	   residue(j) = fac * dcos(mdiff*theres) * vperp2 / (vin*vin) / gammadd

	end do



	end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	----------------------
	sai1L = th(ipo-1) 
	sai1R = th(ipo) 
		pro1 = polang(ipo+1)
	proi1 = polang(ipo-1)
	g1 = gammadd
c	-----------------------




c           END FIRST resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        elseif(ipo>1 .and. aux2(ipo)*aux2(ipo-1)<0 .and. ilast_R.ne.0) then ! 2nd resonance ----------------------------

		  if(aux2(ipo)<aux2(ipo-1))then
		     direct='descen'
		  else
	           direct='ascend'
		  end if

c         ! SECOND resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	      ilast_R2 = ipo      ! first point to the RIGHT of resonance
	      ilast_L2 = ipo-1    ! last  point to the LEFT  of resonance
            insertpts2 = .true.

c           Accuracy
            limit = 1.d-8

c	      Search for resonance (w-wc-k//v// = limit ~ 0)
		ini = th(ipo-1)
		fim = dmin1(th(ipo+1),themax)

		do k=1,80
               thaux = 0.5d0*(ini+fim)
	         call spline0(3, polang(ipo-1:ipo+1), bmotab(intab,ipo-1:ipo+1),
     ;                      1, thaux, btt)
               call spline0(3,polang(ipo-1:ipo+1),kptab(ipo-1:ipo+1),
     ;                      1, thaux, kpp)

		   vpar2  = vin*vin*(1.d0 - xin*btt/B0)
	         vpar   = sigm*dsqrt(abs(vpar2))
		   wc     = qovm * btt
               teste  = omegag - wc - kpp*vpar

		   if(teste.gt.limit)then
                  if(direct.eq.'descen')ini = thaux
			if(direct.eq.'ascend')fim = thaux
	         else
	            if(direct.eq.'descen')fim = thaux
			if(direct.eq.'ascend')ini = thaux
	         end if
               if(abs(teste-limit)/limit<1.d-2)goto 1244

	      end do

1244		continue
	          
c           SECOND resonance angle -------- 
	             theres2 = thaux
c           -------------------------------
	if(themax-theres2<1.d-6)goto 7777

	achei2 = 1
	qua2=teste
		tip2=direct

c           Correct 2 points near SECOND resonance
            deltwa = 0.001d0  ! wanted delta
		
              if(th(ilast_R2).lt. themax)then ! OK  - - - - - - - - - - - - - 	     

c                Correct ipo = ilast_R2 (RIGHT)  +++++++++
                 delti     = dmin1(deltwa, 0.5*(themax-theres2))
	           th(ipo)   = theres2 + delti
		     aux2(ipo) = f_reso_std(theres2 + delti)

c                Correct ipo-1 = ilast_L2 (LEFT) +++++++++
                 delti       = dmin1(delti,theres2)
	           th(ipo-1)   = theres2 - delti
	           aux2(ipo-1) = f_reso_std(theres2 - delti)


		  else  !! CAREFUL (ilast_R2 == itmax2) - - - - - - - - - - - - - - 

                        print*, 'WARNING: Resonance 2 at BANANA tip !!!!!!!'
c	                  STOP
				insertpts2 = .false.
				goto 7777
	
              end if  !! CAREFUL (ilast_R == itmax2) - - - - - - - - - - - - -

ccccccccccccc (2) Residue at theres2 cccccccccccccccccccccccccccccccccccccccc

c     d/dt(gammadot) = |v//|.sin(THETA)/r * d/dtheta(gammadot)
c     where gammadot = w-wc-k//v//

	if(xin>xtg)then

      deltwa = dmin1(0.0003d0,themax-theres2)
c		deltwa = 0.002d0
      call spline0(3, polang(ipo-1:ipo+1),eqt(intab,ipo-1:ipo+1,12),
     ;             1, theres2, sintt)
      vperp2 = vin*vin - vpar2

      gammadd = dabs(vpar)*sintt * (f_ic_std(theres2+deltwa) - f_ic_std(theres2-deltwa)) 
     ;                     / (2.d0*deltwa)
	  
c	Residue
      if(direct.eq.'descen')fac = 1.d0
	if(direct.eq.'ascend')fac = -1.d0
	do j = 1, Nmdiff   ! (Nmdiff = klim + 1)
	   mdiff = dfloat(j - 1)  
	   residue2(j) = fac * dcos(mdiff*theres) * vperp2 / (vin*vin) / gammadd

	end do

	end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	----------------------
	sai2L = th(ipo-1) 
	sai2R = th(ipo) 

	pro2 = polang(ipo+1)
	proi2 = polang(ipo-1)
	g2=gammadd
c	-----------------------

c           END second resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


         end if  ! ---------------------------------------- aux(ipo)*aux(ipo-1)<0 (resonance)

	end do ! ipo ==============================================================================


7777	continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     theta max -----------------------------------------
c     polang(itmax) is last pol. angle lower than themax     
c     Insert additional (themax) value      	
	itmax2 = itmax+1
      th(itmax2)   = themax
	aux2(itmax2) = f_reso_std(themax)
c     ---------------------------------------------------      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






cccccccccc  SPECIAL: CHECK if some resonance was not detected near themax ccccccccccccccc

	if(itmax2>2 .and. aux2(itmax2)*aux2(itmax2-1)<0)then 
      ! THERE IS A RESONANCE BETWEEN th(itmax2-1) and th(itmax2) !!!!


		  if(aux2(itmax2)<aux2(itmax2-1))then
		     direct='descen'
		  else
	         direct='ascend'
		  end if


	if(ilast_R.ne.0)then ! RES2 - - - - - - - - - - - -    

c      ! SECOND resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c        Accuracy
         limit = 1.d-8

c	   Search for resonance (w-wc-k//v// = limit ~ 0)
	   ini = th(itmax2-1)
	   fim = th(itmax2)

	   do k=1,80
            thaux = 0.5d0*(ini+fim)
	      call spline0(3, polang(itmax2-2:itmax2), bmotab(intab,itmax2-2:itmax2),
     ;                   1, thaux, btt)
            call spline0(3,polang(itmax2-2:itmax2),kptab(itmax2-2:itmax2),
     ;                   1, thaux, kpp)

		vpar2  = vin*vin*(1.d0 - xin*btt/B0)
	      vpar   = sigm*dsqrt(abs(vpar2))
		wc     = qovm * btt
            teste  = omegag - wc - kpp*vpar

		   if(teste.gt.limit)then
                  if(direct.eq.'descen')ini = thaux
		      if(direct.eq.'ascend')fim = thaux
	         else
	            if(direct.eq.'descen')fim = thaux
			if(direct.eq.'ascend')ini = thaux
	         end if
            if(abs(teste-limit)/limit<1.d-2)goto 1245

	   end do

1245	   continue
	           
c        SECOND resonance angle -------- 
	          theres2 = thaux
c        -------------------------------
	if(themax-theres2<1.d-6)goto 7879

	achei2 = 1000
	qua2=teste
		tip2=direct

c        Correct points near SECOND resonance in 3 steps:
c        1) correct th(itmax2-1) to theres2-delta  (ilast_L2=itmax2-1)
c        2) create new th(itmax2+1) = themax       (itmax2=itmax2+1)
c        3) correct th(itmax2)   to theres2+delta  (ilast_R2=itmax2)

	   ilast_L2 = itmax2-1    ! last  point to the LEFT 
	   ilast_R2 = itmax2      ! first point to the RIGHT (themax updated later)            
         insertpts2 = .true.

c        1) correct th(itmax2-1)       
            deltwa = 0.001d0  ! wanted delta
		delti = dmin1(deltwa, 0.5d0*(themax - theres2))
            delti = dmin1(delti, theres2)  ! Avoid th<0
	      th(ilast_L2)   = theres2 - delti
	      aux2(ilast_L2) = f_reso_std(theres2 - delti)

c        2) create new th(itmax2+1) = themax 
		th(itmax2+1) = th(itmax2)
	      aux2(itmax2+1)=aux2(itmax2)
		itmax2=itmax2+1

c        3) correct th(itmax2)
		delti = dmin1(deltwa,0.5d0*(themax - theres2))
		th(ilast_R2)   = theres2 + delti
		aux2(ilast_R2) = f_reso_std(theres2 + delti)

ccccccccccccc (3) Residue at theres2 cccccccccccccccccccccccccccccccccccccccc

c     d/dt(gammadot) = |v//|.sin(THETA)/r * d/dtheta(gammadot)
c     where gammadot = w-wc-k//v//

	if(xin>xtg)then

      deltwa = dmin1(0.0003d0,themax-theres2)
c		deltwa = 0.002d0
      call spline0(3, polang(itmax2-2:itmax2),eqt(intab,itmax2-2:itmax2,12),
     ;             1, theres2, sintt)
      vperp2 = vin*vin - vpar2

      gammadd = dabs(vpar)*sintt * (f_ic_std(theres2+deltwa) - f_ic_std(theres2-deltwa)) 
     ;                     / (2.d0*deltwa)
	
c	Residue	  
      if(direct.eq.'descen')fac = 1.d0
	if(direct.eq.'ascend')fac = -1.d0
	do j = 1, Nmdiff ! (Nmdiff = klim + 1)
	   mdiff = dfloat(j - 1)  
	   residue2(j) = fac*dcos(mdiff*theres) * vperp2 / (vin*vin) / gammadd

	end do

	end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	----------------------
	sai2L = th(ilast_L2) 
	sai2R = th(ilast_R2) 

	pro2 = polang(itmax2)
	proi2 = polang(itmax2-2)
	g2=gammadd
c	-----------------------

c      ! END SECOND resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	elseif(ilast_R.eq.0)then ! RES1 - - - - - - - - - - - - - - -  

c      ! FIRST resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c        Accuracy
         limit = 1.d-8

c	   Search for resonance (w-wc-k//v// = limit ~ 0)
	   ini = th(itmax2-1)
	   fim = th(itmax2)

	   do k=1,80
            thaux = 0.5d0*(ini+fim)
	      call spline0(3, polang(itmax2-2:itmax2), bmotab(intab,itmax2-2:itmax2),
     ;                   1, thaux, btt)
	      call spline0(3,polang(itmax2-2:itmax2),kptab(itmax2-2:itmax2),
     ;                   1, thaux, kpp)

		vpar2  = vin*vin*(1.d0 - xin*btt/B0)
	      vpar   = sigm*dsqrt(abs(vpar2))
		wc     = qovm * btt
            teste  = omegag - wc - kpp*vpar

		   if(teste.gt.limit)then
                  if(direct.eq.'descen')ini = thaux
			if(direct.eq.'ascend')fim = thaux
	         else
	            if(direct.eq.'descen')fim = thaux
			if(direct.eq.'ascend')ini = thaux
	         end if
            if(abs(teste-limit)/limit<1.d-2)goto 1246

	   end do

1246	   continue
	           
c        FIRST resonance angle -------- 
	          theres = thaux
c        -------------------------------
	if(themax-theres<1.d-6)goto 7879
	achei1 = 1000
	qua1=teste
		tip1=direct

c        Correct points near FIRST resonance in 3 steps:
c        1) correct th(itmax2-1) to theres-delta  (ilast_L=itmax2-1)
c        2) create new th(itmax2+1) = themax      (itmax2=itmax2+1)
c        3) correct th(itmax2)   to theres+delta  (ilast_R=itmax2)

	   ilast_L = itmax2-1    ! last  point to the LEFT  
	   ilast_R = itmax2      ! first point to the RIGHT (themax updated later)           
         insertpts2 = .true.

c        1) correct th(itmax2-1)       
            deltwa = 0.001d0  ! wanted delta
		delti = dmin1(deltwa, 0.5d0*(themax - theres))
            delti = dmin1(delti, theres)  ! Avoid th<0
	      th(ilast_L)   = theres - delti
	      aux2(ilast_L) = f_reso_std(theres - delti)

c        2) create new th(itmax2+1) = themax 
		th(itmax2+1)   = th(itmax2)
	      aux2(itmax2+1) = aux2(itmax2)
		itmax2 = itmax2+1

c        3) correct th(itmax2)
		delti = dmin1(deltwa,0.5d0*(themax - theres))
		th(ilast_R)   = theres + delti
		aux2(ilast_R) = f_reso_std(theres + delti)

ccccccccccccc (4) Residue at theres cccccccccccccccccccccccccccccccccccccccc

c     d/dt(gammadot) = |v//|.sin(THETA)/r * d/dtheta(gammadot)
c     where gammadot = w-wc-k//v//

	if(xin>xtg .or. direct.eq.'ascend')then

      deltwa = dmin1(0.0003d0,themax-theres)
c	deltwa = 0.002d0
      call spline0(3, polang(itmax2-2:itmax2),eqt(intab,itmax2-2:itmax2,12),
     ;             1, theres, sintt)
      vperp2 = vin*vin - vpar2

      gammadd = dabs(vpar)*sintt * (f_ic_std(theres+deltwa) - f_ic_std(theres-deltwa)) 
     ;                     / (2.d0*deltwa)

c	Residue	  
      if(direct.eq.'descen')fac = 1.d0
	if(direct.eq.'ascend')fac = -1.d0
	do j = 1, Nmdiff  ! (Nmdiff = klim + 1)
	   mdiff = dfloat(j - 1)  
	   residue(j) = fac*dcos(mdiff*theres) * vperp2 / (vin*vin) / gammadd

	end do

	end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	----------------------
	sai1L = th(ilast_L) 
	sai1R = th(ilast_R) 

	pro1 = polang(itmax2)
	proi1 = polang(itmax2-2)
	g1=gammadd
c	-----------------------

c      ! END FIRST resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	end if  ! RES2 or RES1 - - - - - - - - - - - - - - - - 

      tricky = .true.

	end if   ! RESONANCE DETECTED near itmax2 !!!!

cccccccccc  END SPECIAL: CHECK if some resonance was not detected near themax ccccccccccc


7879	continue



c		ALL roots determined !!!!!	
      resid(1:Nmdiff) = 0.5*pi*(residue(1:Nmdiff)+residue2(1:Nmdiff))

		if(only_active)goto 9876


cccccccccccccccccccccccccccccccc Refine theta mesh cccccccccccccccccccccccccc
    
c       1) Insert points near turning point (themax) ------------------------

        Nins = 30
        delti = 0.d0

         if(itmax2>1)then

            ! original point (themax)
            urtimo = aux2(itmax2)

            if(ilast_R2>0)then
	         delti = dmin1(2.d-2,0.5*(th(itmax2)-th(ilast_R2)))
            elseif(ilast_R>0)then
               delti = dmin1(2.d-2,0.5*(th(itmax2)-th(ilast_R)))	  
		else
               delti = dmin1(3.d-2,th(itmax2))
	      end if

c        Insert points (always LEFT)         
         call insert_points_deltax(itmax2, th(1:itmax2), th(itmax2), 'L', 
     ;                             delti, Nins, 0.d0, th(1:itmax2+Nins))
 
      ! update new points
	
	  if(itmax2>2)then
           do ipo = itmax2-2, itmax2+Nins-1
              aux2(ipo)=f_reso_std(th(ipo))
           end do
	  else
           do ipo = itmax2-1, itmax2+Nins-1
              aux2(ipo)=f_reso_std(th(ipo))
           end do
	  end if

	! update original point
        aux2(itmax2+Nins)=urtimo
	  itmax2 = itmax2 + Nins

c     ! Store left-limit of point insertion
        theins = th(itmax2)-delti

        end if

c     -----------------------------------------------------------------------


c       2) Insert points near SECOND resonance (theres2) --------------------

        ! RIGHT side of resonance
        if(ilast_R2>1 .and. ilast_R2<itmax2)then

           Nins = 30
	     delti = dmin1(2.d-2,theins-th(ilast_R2))
              if(delti>0.d0)then
		     ! Shift original values
                 aux2(ilast_R2+1+Nins:itmax2+Nins)=aux2(ilast_R2+1:itmax2)
                 ! Insert new points
                 call insert_points_deltax(itmax2, th(1:itmax2), th(ilast_R2), 'R', 
     ;                                     delti, Nins, 0.d0, th(1:itmax2+Nins))
                 ! Update new points
                 do ipo = ilast_R2+1, ilast_R2+1+Nins
                    aux2(ipo) = f_reso_std(th(ipo))
                 end do
                 itmax2 = itmax2+Nins
              end if
        end if

        ! LEFT side of resonance
        if(ilast_L2>2)then

           Nins = 30
           delti = dmin1(delti,th(ilast_L2)-th(ilast_R))
           if(delti>0.d0)then
		     ! Shift original values
	           aux2(ilast_L2+Nins:itmax2+Nins)=aux2(ilast_L2:itmax2)
c	           ! Insert new points
                 call insert_points_deltax(itmax2, th(1:itmax2), th(ilast_L2), 'L', 
     ;                                     delti, Nins, 0.d0, th(1:itmax2+Nins))
                 ! update new points
                 do ipo = ilast_L2-1, ilast_L2+Nins-1
                    aux2(ipo)=f_reso_std(th(ipo))
                 end do
                 itmax2 = itmax2+Nins
c	           ilast_R = ilast_R+Nins  ! Don't need it any more (from RIGHT to LEFT)
              end if

c     ------------------------------------------------------------------------------------
c            put some extra points to cover original delta
              if(delti<2.d-2 .and. delti>0.d0)then
	           
                 if(xin > 0.995d0*xo)then
                    Nins=60
	              delti = dmin1(2.0d-2-delti,th(ilast_L2)-th(ilast_R) )
	           else
                    Nins=30
	              delti = dmin1(2.d-2-delti,th(ilast_L2)-th(ilast_R))
		     end if
                 
                 ! Shift original values
	           aux2(ilast_L2+Nins:itmax2+Nins)=aux2(ilast_L2:itmax2)
c	           ! Insert new points
                 call insert_points_deltax(itmax2, th(1:itmax2), th(ilast_L2), 'L', 
     ;                                     delti, Nins, 0.d0, th(1:itmax2+Nins))
                 ! update new points
                 do ipo = ilast_L2-1, ilast_L2+Nins-1
                    aux2(ipo)=f_reso_std(th(ipo))
                 end do
                 itmax2 = itmax2+Nins
	        end if
c     --------------------------------------------------------------------------------

        end if


c     -----------------------------------------------------------------------


c       2) Insert points near FIRST resonance (theres) --------------------


        ! RIGHT side of resonance
        if(ilast_R>1 .and. ilast_R<itmax2)then

           Nins = 30
	     delti = dmin1(2.d-2,theins-th(ilast_R))
              
              if(delti>0.d0)then
		     ! Shift original values
                 aux2(ilast_R+1+Nins:itmax2+Nins)=aux2(ilast_R+1:itmax2)
                 ! Insert new points
                 call insert_points_deltax(itmax2, th(1:itmax2), th(ilast_R), 'R', 
     ;                                     delti, Nins, 0.d0, th(1:itmax2+Nins))
                 ! Update new points
                 do ipo = ilast_R+1, ilast_R+1+Nins
                    aux2(ipo) = f_reso_std(th(ipo))
                 end do
                 itmax2 = itmax2+Nins
              end if
        end if

        ! LEFT side of resonance
        if(ilast_L>2)then

           Nins = 30
	!    same as befoire: delti = dmin1(3.d-2,theins-th(ilast_R))
           delti = dmin1(delti, th(ilast_L)-0.d0)

              if(delti>0.d0)then
		     ! Shift original values
	           aux2(ilast_L+Nins:itmax2+Nins)=aux2(ilast_L:itmax2)
c	           ! Insert new points
                 call insert_points_deltax(itmax2, th(1:itmax2), th(ilast_L), 'L', 
     ;                                     delti, Nins, 0.d0, th(1:itmax2+Nins))
                 ! update new points
                 do ipo = ilast_L-1, ilast_L+Nins-1
                    aux2(ipo)=f_reso_std(th(ipo))
                 end do
                 itmax2 = itmax2+Nins
c	           ilast_R = ilast_R+Nins  ! Don't need it any more (from RIGHT to LEFT)
              end if

c     ------------------------------------------------------------------------------------
c            put some extra points to cover original delta
              if(delti<2.d-2 .and. delti>0.d0)then
                 if(xin > 0.995d0*xo)then
                    Nins=60
	              delti = dmin1(2.0d-2-delti,th(ilast_L))
	           else
	              Nins=30
                    delti = dmin1(2.d-2-delti,th(ilast_L))
		     end if
                 ! Shift original values
	           aux2(ilast_L+Nins:itmax2+Nins)=aux2(ilast_L:itmax2)
c	           ! Insert new points
                 call insert_points_deltax(itmax2, th(1:itmax2), th(ilast_L), 'L', 
     ;                                     delti, Nins, 0.d0, th(1:itmax2+Nins))
                 ! update new points
                 do ipo = ilast_L-1, ilast_L+Nins-1
                    aux2(ipo)=f_reso_std(th(ipo))
                 end do
                 itmax2 = itmax2+Nins
	        end if
c     --------------------------------------------------------------------------------


        end if


ccccccccccccccccccccccccccc

cccccccc Testing xtg>xtir cccccccc

c      if(xtg<xo)then
c	if(xin>xo)then

c	resid(1:Nmdiff) = 0.5*pi*residue(1:Nmdiff)

c	else


      resid(1:Nmdiff) = 0.5*pi*(residue(1:Nmdiff)+residue2(1:Nmdiff))
c       else
c      resid(1:Nmdiff) = -0.5*pi*(residue(1:Nmdiff)+residue2(1:Nmdiff))
c	end if
c	end if






c            Integration

	       
c             do ipo = 2, itmax2

c                  if( ipo>1 .and. aux2(ipo)*aux2(ipo-1)>0 )then  ! exclude singularity
c				   soma2 = soma2 + 0.5d0*(aux2(ipo)+aux2(ipo-1))*(th(ipo)-th(ipo-1))
c			end if			 

c		 end do  

c            	 resultat = soma
c                   res2 = soma2


	soma2(1:Nmdiff) = 0.d0

	do j = 1, Nmdiff

	   mdiff = dfloat(j - 1)
         
	   do ipo = 2, itmax2
              	
			if( ipo>1 .and. aux2(ipo)*aux2(ipo-1)>0 )then  ! exclude singularity
				soma2(j) = soma2(j)
     ;            + 0.5d0 * ( aux2(ipo)*cos(mdiff*th(ipo)) + aux2(ipo-1)*cos(mdiff*th(ipo-1)) )
     ;                    * ( th(ipo) - th(ipo-1) )
			end if			 
  
         end do  

	end do

              res2(1:Nmdiff) = soma2(1:Nmdiff)



c		 epsabs = 1.d-10
c		 epsrel = 1.d-4

c      call dqagse(f_reso,0.d0,polang(itmax), 
c     ;            epsabs,epsrel,200,res2,abserr,neval,
c     ;            ier,alist,blist,rlist,elist,iord,last)




c     =========================================================================

9876	continue

	if(FLAG)then

      if(FLAG .and. vin.gt.4.5d6)then

c     Write data to file (only for debugging) ----------------------

c                write(1552,"(400G16.8)"), xin, xnout, npfft/2+1, polang(itmax), (aux(1:npfft/2+1))
c                write(1552,"(400G16.8)"), xin, xnout, npfft/2+1, polang(itmax), polang(1:npfft/2+1)
c                write(1552,"(400G16.8)"), xin, xnout, itmax2, th(itmax2), (aux2(1:itmax2))
c                write(1552,"(400G16.8)"), xin, xnout, theres, theres2, th(1:itmax2)

c               write(1553,"(400G16.8)")
c               write(1553,"(400G16.8)"), xin, xnout, xtg, th(itmax2)
c			 write(1553,"(400G16.8)"), theres, achei1, tip1, qua1, g1, residue(1) 
c		     write(1553,"(400G16.8)"), theres2, achei2, tip2, qua2, g2, residue2(1)


      end if


	end if



      return
 
      END SUBROUTINE orbit_integ_trap_std
