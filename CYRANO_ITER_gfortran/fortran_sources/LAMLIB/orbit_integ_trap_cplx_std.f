      SUBROUTINE orbit_integ_trap_cplx_std(xin, vin, sigm, mav, Nmdiff, qovm, 
     ;                                 resultat, res2, FLAG)

      IMPLICIT NONE

c     Compute orbit integral for given (x,v,k//,mdiff)
																	
c	TRAPPED - complex branch (non-resonant)

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
	real*8, intent(out) :: resultat, res2(1:Nmdiff) 
	logical, intent(in) :: FLAG
     
!     Variables
      integer j, ipo, itmax
	real*8 :: xmax, vpar2, vpar, vperp2
	real*8 :: wc, aux(npfft+1), aux2(npfft+1+100), soma, th(npfft+1+100), 
     ;          soma2(1:Nmdiff), mdiff

	real*8 f_reso_std
	external f_reso_std

	integer last, ind2, itmax2, isea

      integer Nins

      real*8 :: teste, thaux, btt, ini, fim, kp, kpp
	real*8 :: limit, urtimo, delti, xnout, themax


	integer ISRCHFGE
	external ISRCHFGE

c     Set global variables (for f_reso.f routine)
      VNOW = vin
	XNOW = xin
	sigNOW = sigm
	qomNOW = qovm
c	kpNOW = kp
	mavNOW = mav  !!!!!!!!!!!!!!!!!!!!!!

c	k//(theta) passed through COMMON kptab(1:npfft+1)  
c	          (DIRESP_standard_1 line 489)

c	Radius dependent quantities
c      xnsep = - 2.d0 * delb(intab) / bmax(intab) 
      xmax = B0 / bmin(intab)
c      xsep = (1.d0-dabs(xnsep))*xmax
	xnout = -sigm*(1.d0-xin/xmax)

	th(1:npfft/2+1) = polang(1:npfft/2+1)
      th(npfft/2+2:npfft) = pi

	aux(1:npfft/2+1)=0
	aux2(1:npfft/2+1+100)=0
c
	themax=0.d0



		if(only_active)RETURN




c          Find itmax in fine table
           ind2 = isrchfge(1001, xtabchi, 1, xin)
	     itmax = ithemaxchi(ind2)
	     itmax2 = itmax


cccccccccccccccccccccccccccccccccccccccccccccccccc

c	First find theta maximum

	   ! Find theta max where v// -> 0 and put in th(itmax+1)

         if(itmax2>0)then  !!!! .and. itmax2.gt.ilast_R
                              

	       limit = 1.d-6 !vmin*vmin / (vin*vin)
		   isea=4			  
	       ! Near stagnation orbit
		   if(dabs(xnout).le.2.d-3)isea=20
             
		   ini = th(itmax2)
             fim = th(itmax2+isea)  ! +4
                         
             do k=1,60

                thaux = 0.5d0*(ini+fim)
	          call spline0(isea+1,polang(itmax2:itmax2+isea),bmotab(intab,itmax2:itmax2+isea),
     ;                            1, thaux, btt)

                teste  = (1.d0 - xin*btt/B0)
	  
	          if(teste.gt.limit)then
	             ini = thaux
	          else
	             fim = thaux
	          end if
                if(abs(teste-limit)/limit<1.d-2)goto 1223

             end do
	                      
1223  continue
             
c		   ----- theta max -----			       
		       themax = thaux
c		   ---------------------

                    th(itmax2+1) = thaux
				aux2(itmax2+1) = f_reso_std(thaux)
	                  
                        ! in case original value was already higher then thaux
                        if(th(itmax2)>thaux)then
                           th(itmax2) = th(itmax2+1)
	                   aux2(itmax2) = aux2(itmax2+1)
                        end if  
                      

                      itmax2 = itmax2+1

					
                      end if	





ccccccccccccccccccccccccccccccccccccccccccccccccccc



c	Loop over poloidal angle

	do ipo = 1, itmax+1  ! =================================================================
			 
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


c	   Last point (itmax+1 == itmax2) -----------------
	   if(ipo.eq.itmax+1)then
	          aux(ipo)=0.d0  ! dsign(1.d0,aux(ipo))*aux(ipo-1)
c		--> aux2(ipo) already computed above == theta max
	   else
            aux2(ipo)= aux(ipo)
	   end if
c	   --------------------------------------------------

c	   Detect IC resonance point	
         if(ipo>1 .and. aux2(ipo)*aux2(ipo-1)<0) then ! resonance -------
            print*, 'There is a problem in orbit_integ_trap_cplx.f !!!'
	      print*, 'IC resonance found where it should NOT exist  !!!'
c	      stop
         end if  ! ---------------------------------------- ipo*(ipo-1)<0  


	end do ! ipo ==============================================================================







cccccccccccccccccccccccccccccccccc last points
c     c       Insert points near v//=0 point (itmax+1)

        Nins = 30  ! MAX 100
        delti = 0.d0


         if(itmax2>1)then

            urtimo = aux2(itmax2)

            delti = dmin1(2.d-2,th(itmax2))
	       ! Near stagnation orbit
		   if(dabs(xnout).le.2.d-3)delti = dmin1(1.d-1,th(itmax2))

        call insert_points_deltax(itmax2, th(1:itmax2), th(itmax2), 'L', delti, Nins, 0.d0, th(1:itmax2+20))
 
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


        end if






c            Integration


c	             soma2 = 0.d0	       
c             do ipo = 2, itmax2

c                  if( ipo>1 .and. aux2(ipo)*aux2(ipo-1)>0 )then  ! exclude singularity
c				   soma2 = soma2 + 0.5d0*(aux2(ipo)+aux2(ipo-1))*(th(ipo)-th(ipo-1))
c			end if			 

c		 end do  



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

c     =========================================================================

      if(FLAG .and. vin.gt.0.9d6 .and. vin<1.1d6)then

c     Write data to file (only for debugging) ----------------------

c                write(1552,"(200G16.8)"), xin, xnout, npfft/2+1, polang(itmax), (aux(1:npfft/2+1))
c                write(1552,"(200G16.8)"), xin, xnout, npfft/2+1, polang(itmax), polang(1:npfft/2+1)
c                write(1552,"(200G16.8)"), xin, xnout, itmax2, th(itmax2), (aux2(1:itmax2))
c                write(1552,"(200G16.8)"), xin, xnout, delti, 0.d0, th(1:itmax2)

      end if


      return
 
      END SUBROUTINE orbit_integ_trap_cplx_std
