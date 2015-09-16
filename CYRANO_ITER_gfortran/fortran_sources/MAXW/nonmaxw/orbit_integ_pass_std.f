      SUBROUTINE orbit_integ_pass_std(xin, vin, sigm, mav, Nmdiff, qovm, xtg, xo,
     ;                                     resultat, res2, residue, FLAG)

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
	real*8, intent(in)  :: xin, vin, mav, qovm , sigm, xtg, xo
	real*8, intent(out) :: resultat, res2(1:Nmdiff), residue(1:Nmdiff) 
	logical, intent(in) :: FLAG
c      real*8 sigm
     
!     Variables
      integer j, ipo, itmax
	real*8 :: xmax, xsep, vpar2, vpar, vperp2
	real*8 :: wc, aux(npfft+1), aux2(npfft+1+100), soma, th(npfft+1+100), soma2(1:Nmdiff)

	real*8 f_reso_std
	external f_reso_std

	real*8 f_ic_std
	external f_ic_std

	integer itmax2, ilast_R, ilast_L, Nins

      real*8 :: teste, thaux, sintt, btt, ini, fim
	real*8 :: sai1, sai2, limit, delti, xnout
	real*8 :: gammadd, theres, mdiff, kp, kpp 


	integer ISRCHFGE
	external ISRCHFGE



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
	aux2(1:npfft/2+1+100)=0
	ilast_L=0
	ilast_R=0
	sai1=0
c	sai2=0

	residue(1:Nmdiff)=0.d0
                  
           itmax  = npfft/2+1
           itmax2 = itmax
                  




c	Loop over poloidal angle

	do ipo = 1, itmax  ! =================================================================
			 
	   vpar2  = vin*vin*(1.d0 - xin*bmotab(intab,ipo)/B0)
	   vpar2  = dmax1(vpar2,1.d0)
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

         if(ipo>1 .and. aux(ipo)*aux(ipo-1)<0) then ! resonance ----------------------------

	      ilast_R = ipo      ! first point to the RIGHT of resonance
	      ilast_L = ipo-1    ! last  point to the LEFT  of resonance
           
 

          if(ipo>1)then

            limit = 1.d-8
              
c	      Search for w-wc-k//v//=limit~0

			ini = polang(ipo-1)
			fim = polang(ipo+1)
			
			do k=1,100
                     thaux = 0.5d0*(ini+fim)
	               call spline0(3,polang(ipo-1:ipo+1),bmotab(intab,ipo-1:ipo+1),
     ;                            1, thaux, btt)
	               call spline0(3,polang(ipo-1:ipo+1),kptab(ipo-1:ipo+1),
     ;                            1, thaux, kpp)

			   vpar2  = vin*vin*(1.d0 - xin*btt/B0)
	             !vpar2  = dmax1(vpar2,1.d0)
	           vperp2 = vin*vin - vpar2
	           vpar   = sigm*dsqrt(dabs(vpar2))
			   wc     = qovm * btt

                 teste  = omegag - wc - kpp*vpar
c			   teste2 = (omegag - wc)/(kp*vpar)

			  if(teste.gt.limit)then
                       ini = thaux
	              else
	                 fim = thaux
	              end if
                    if(dabs(teste-limit)/limit<1.d-2)goto 1233

	         end do

1233		  continue
	           
              theres = thaux

		end if

 
cccccccccccccccccccccccccccccccccccc

c           Correct ipo-1 ++++++++++++++++++++++++++++++++++++++



	           th(ipo-1) = thaux-0.002d0
			   aux2(ipo-1)= f_reso_std(thaux-0.002d0)



c           Correct ipo ++++++++++++++++++++++++++++++++++++++

	           th(ipo) = thaux+0.002d0
		     aux2(ipo)= f_reso_std(thaux+0.002d0)




ccccccccccccccccccccccccccccccccccccccc


ccccccccccccc Residue at theres cccccccccccccccccccccccccccccccccccccccc

c     d/dt(gammadot) = |v//|.sin(THETA)/r * d/dtheta(gammadot)
c     where gammadot = w-wc-k//v//
	if(xin>xtg)then

      call spline0(3, polang(ipo-1:ipo+1),eqt(intab,ipo-1:ipo+1,12),
     ;             1, theres, sintt)
      gammadd = dabs(vpar)*sintt * (f_ic_std(theres+2.d-3) - f_ic_std(theres-2.d-3)) 
     ;                     / (2*2.d-3)

c	(Nmdiff = klim + 1)  
	do j = 1, Nmdiff
	   mdiff = dfloat(j - 1)  
	   residue(j) = 0.5d0 * pi * cos(mdiff*theres) * vperp2 / (vin*vin) / gammadd
	end do

      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



         end if  !  ---------------------------------------- ipo*(ipo-1)<0


2323     continue


	end do ! ipo ==============================================================================


		if(only_active)goto 9876



	sai1 = 0.d0



        if(ilast_L>1)then
        ! LEFT side of resonance
c        delti = (th(ilast_L)-th(ilast_L-1))/10.d0
        Nins = 40
	  aux2(ilast_L+Nins:itmax2+Nins)=aux2(ilast_L:itmax2)
c	  call insert_points_simple(itmax2, th(1:itmax2), ilast_L-2, ilast_L, Nins, th(1:itmax2+Nins))

        delti=dmin1(4.d-2,th(ilast_L))	
        call insert_points_deltax(itmax2, th(1:itmax2), th(ilast_L), 'L', delti, Nins, 0.d0, th(1:itmax2+Nins))
      ! update new points
        do ipo = ilast_L-1, ilast_L+Nins-1
            aux2(ipo)=f_reso_std(th(ipo))
        end do

        itmax2 = itmax2+Nins
	  ilast_R = ilast_R+Nins

        end if

        if(ilast_R>1 .and. ilast_R<itmax2)then
        ! RIGHT side of resonance
c        delti = (th(ilast_R+1)-th(ilast_R))/10.d0
        Nins = 40
	  aux2(ilast_R+1+Nins:itmax2+Nins)=aux2(ilast_R+1:itmax2)
c	  call insert_points_simple(itmax2, th(1:itmax2), ilast_R, ilast_R+2, Nins, th(1:itmax2+Nins))

        delti=dmin1(4.d-2,pi-th(ilast_R))	
        call insert_points_deltax(itmax2, th(1:itmax2), th(ilast_R), 'R', delti, Nins, 0.d0, th(1:itmax2+Nins))

      ! update new points
        do ipo = ilast_R+1, ilast_R+1+Nins
            aux2(ipo)=f_reso_std(th(ipo))
        end do

      itmax2 = itmax2+Nins
      end if




c            Integration

c			 soma = 0.d0



c	       do ipo = 2, itmax

c                  if( ipo>1 .and. aux(ipo)*aux(ipo-1)>0 )then  ! exclude singularity
c				   soma = soma + 0.5d0*(aux(ipo)+aux(ipo-1))*(polang(ipo)-polang(ipo-1))
c			end if			 

c		 end do    

c             itmax2 = itmax+1
	       




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

      if(FLAG .and. vin.gt.4.5d6)then

c     Write data to file (only for debugging) ----------------------

c                write(1552,"(200G16.8)"), xin, xnout, npfft/2+1, polang(itmax), (aux(1:npfft/2+1))
c                write(1552,"(200G16.8)"), xin, xnout, npfft/2+1, polang(itmax), polang(1:npfft/2+1)
c                write(1552,"(200G16.8)"), xin, xnout, itmax2, th(itmax2), (aux2(1:itmax2))
c                write(1552,"(200G16.8)"), xin, xnout, sai2, sai1, th(1:itmax2)

c               write(1553,"(400G16.8)")
c               write(1553,"(400G16.8)"), xin, xnout, xtg, pi
c			 write(1553,"(400G16.8)"), theres, 0, 0, 0, 0, residue 
c		     write(1553,"(400G16.8)"), 0, 0, 0, 0, 0, 0


      end if


      return
 
      END SUBROUTINE orbit_integ_pass_std
