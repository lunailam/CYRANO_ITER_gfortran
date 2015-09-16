      SUBROUTINE thetamax_CHI(LX,xx)

      IMPLICIT NONE

c     Compute maximum poloidal angle of the orbits for 
c     given x-mesh : xx = [0...Xmax]

c     INPUT  : LX (length of xx), xx (vector with x values)

c     OUPUT (COMMONS):
c       - ithemaxchi(x)    : index of max. theta value
																		
      include 'pardim.copy'
      include 'commag.copy'
      include 'comphy.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'cokpco.copy'
      include 'comgeo.copy'
      include 'compow.copy'

		 
!	Input	                                                                 
	integer, intent(in) :: LX
	real*8, intent(in)  :: xx(LX)
 
     
!     Variables
      integer j, k, OpenStat, ipo
      character(100) :: FILE_NAME
	real*8 :: xmax, xsep, Btot(1:npfft+1), aux(1:npfft+1)
     ;          , ini, fim, limit, btt, teste, thaux,  auxs(2000)



c	Radius dependent quantities
      xnsep = - 2.d0 * delb(intab) / bmax(intab) 
      xmax = B0 / bmin(intab)
      xsep = (1.d0-dabs(xnsep))*xmax
	Btot(1:npfft+1) = bmotab(intab,1:npfft+1)



c     x - loop ===============================================

	do j = 1, LX   

         xtabchi(j) = xx(j)         

         if(xtabchi(j) .eq. 0.0d0)then     ! Special case: x = 0
c        ----------------------------------------------------        

		   ithemaxchi(j) = npfft/2+1
			theta_max(j) = pi

         elseif(xtabchi(j) .eq. xmax)then  ! Special case: x = Xmax
c        -------------------------------------------------------   

		   ithemaxchi(j) = 1
			theta_max(j) = 0.d0

         elseif(xtabchi(j) .eq. xsep)then  ! Special case: x = xsep
c        -------------------------------------------------------  

		   ithemaxchi(j) = npfft/2+1
		   theta_max(j) = pi

         else   ! General case : 0 < x < xsep(pass) and xsep < x < Xmax(trap)
c        --------------------------------------------------------------------               


            if(xtabchi(j) < xsep)then ! >>>>>>>>>>>>>>>>>>>  passing orbits
               
               ithemaxchi(j) = npfft/2+1
		     theta_max(j) = pi
                                                
            else                ! >>>>>>>>>>>>>>>>>>>>>>>  trapped orbits    
               
			 ! Find roots of 1-x.Btot(rho,theta)/Bax=0
			 
			 
			 aux(1) = 1.d0 - xtabchi(j)*Btot(1)/B0

			 do ipo = 2, npfft/2+1
			 
				aux(ipo) = 1.d0 - xtabchi(j)*Btot(ipo)/B0

				if(aux(ipo)*aux(ipo-1) < 0.d0)then

				   ithemaxchi(j) = ipo-1

					! Refine search for theta_max
                        ini = polang(ipo-1)
	                  fim = polang(ipo)
					  limit = 1.d-12
                            
                        do k=1,30
                           thaux = 0.5d0*(ini+fim)
	                      call spline0(npfft/2+1,polang(1:npfft/2+1),bmotab(intab,1:npfft/2+1),
     ;                                   1, thaux, btt)

                            teste  = (1.d0 - xtabchi(j)*btt/B0)

                            if(abs(teste-limit)/limit<1.d-2)goto 1234

                               if(teste.gt.limit)then
	                            ini = thaux
	                         else
	                            fim = thaux
	                         end if

	                  end do

1234  continue	

					 theta_max(j) = thaux
					 auxs(j) = teste

				     goto 111

				end if
			 
			 end do     



111		  continue				        
						   
            end if  ! >>>>>>>>>>>>>>>>>>>>>>> 


         end if ! Special cases (x=0, x=Xmax or x=xsep) ----------------------
             
		      
	end do    ! x-loop


c     =========================================================================


c     Write data to file (only for debugging) ----------------------


		if(cokpco)then
          FILE_NAME = trim(COKFOLDER)  // "/thetamax_chi.dat" 
	    else
	    FILE_NAME = trim(STDFOLDER)  // "/thetamax_chi.dat"
		end if

	    open (UNIT = 999, FILE = FILE_NAME, STATUS = "REPLACE",
     ;          IOSTAT = OpenStat, ACTION = "WRITE")
	          if (OpenStat > 0) then
		        print *, 'Error writing file: ', FILE_NAME
	              stop
		    end if	  
		    write(999,*),"Table of thetamax_chi"
		    write(999,*),"rho=", abscis(intab) 
		    write(999,*),"xsep=", xsep 
                write(999,*),"xmax=", xmax		      
		    write(999,*),"Nx "
                write(999,*), LX
                do j = 1,LX 
                   write(999,"(G16.8,I8,G16.8,G16.8,G16.8)"), xtabchi(j), ithemaxchi(j),
     ;                                            polang(ithemaxchi(j)), theta_max(j), auxs(j) 
	          end do
           close(999)




      return
 
      END SUBROUTINE thetamax_CHI
