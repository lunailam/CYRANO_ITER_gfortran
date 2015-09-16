      subroutine ypro_read(inscr)

      implicit none
      logical inscr

c     computes density, temperature, thermal velocity and plasma freq.
c     for all species at local abscissa y and angle phi.
c     y,ireg passed by comswe
 
      include 'pardim.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comreg.copy'
      include 'comant.copy'
      include 'compla.copy'
      include 'comswe.copy'
      include 'comphy.copy'
      include 'coequi.copy'

      integer i, file_stat
      double precision y2, y3

      do 3 i = 1, nspec
      vt(i) = 0.d0
      omp(i) = 0.d0
      densit(i) = 0.d0
      temper(i) = 0.d0
  3   continue
 
      if(vacuum(ireg))return
 
      if(inscr)then
 
        if(.not.coldpl(ireg))then
	      do i = 1, nspec
	      vt(i) = dsqrt(2.d0 * tfac * tinscr(i) / (mh*amass(i)))
	      temper(i) = tinscr(i)
	      end do
        end if

	    do 2 i = 1, nspec
	    densit(i) = ninscr(i)
	    omp(i) = dsqrt(dfac/(amass(i)*mh*eps0)
     ;         * ninscr(i)) * abs(zch(i)*eel)
  2     continue
      else
 
c a revoir: parametris.

      y3 = (y - rx0(ireg-1)) / rl(ireg)
	y3 = dmax1(y3, 0.d0)
	y3 = dmin1(y3, 1.d0)

        if(ireg.eq.1)then
        y2 = 1.d0 - y3 * y3
        else
        y2 = 1.d0 - y3
        end if
 
        if( abs(y3-1.d0) .lt. 1.d-4 )then
          do 1 i = 1,nspec
 
          omp(i) = dsqrt( dfac/(amass(i)*mh*eps0)
     ;                    * nb(i,ireg) ) * dabs(zch(i)*eel)
          densit(i) = nb(i,ireg)
            if(.not.coldpl(ireg) .or. sccol(ireg))then
            vt(i) = dsqrt( 2.d0 * tfac/(mh*amass(i)) * tb(i,ireg) )
            temper(i) = tb(i,ireg)
              if(sccol(ireg))then
c               Electron-ion collisions (formula is for deuterium)
c               coulomb log = 15
                damp(1,ireg) = 1.55d9 * densit(1) / temper(1) ** 1.5
              end if
            end if
  1       continue
 
        else
          
		do 4 i = 1, nspec

! Particle density
      file_stat = 0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc		
		  call read_prof('density', paname(i), y, file_stat, 
     ;                     densit(i) )

c	Correction for density values (Cyrano uses 1e+20/m3)		
			densit(i)=densit(i)*1.0d-1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	
		  if (file_stat .eq. 0) then 
                densit(i) =
     ;          ( n0(i,ireg) - nb(i,ireg) ) * y2 ** idexp(ireg)
     ;          + nb(i,ireg)
	      end if

! Plasma frequency
          
		  omp(i) = dsqrt( dfac * densit(i) / (amass(i)*mh*eps0) ) 
     ;               * dabs(zch(i)*eel)
 
! Temperature

      file_stat = 0          
			
cERN		if(.not.coldpl(ireg) .or. sccol(ireg))then	

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		    call read_prof('temperature', paname(i), y, file_stat,
     ;                        temper(i) )
c			print *, rho, y
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		  

	  	
		     if (file_stat .eq. 0) then 
                   temper(i) =
     ;             ( t0(i,ireg)-tb(i,ireg) ) * y2 ** itexp(i,ireg)
     ;             + tb(i,ireg)
		     end if

! Thermal velocity
          
		     vt(i) = dsqrt( 2.d0 * tfac /(mh*amass(i)) * temper(i) )

               if(sccol(ireg))then
c                 electron-ion collisions (formula is for deuterium)
c                 coulomb log = 15
                  damp(1,ireg) = 1.55d9 * densit(1) / temper(1) ** 1.5
               end if
            
cERN		 end if ! NOT coldpl OR sccol

4         continue
          end if
 
      end if
 
 	if(SAME_DENS)then
	  do i = 2,nspec
	     densit(i) = densit(1) * spefra(i,1)
	  end do
	end if

 	if(SAME_TEMP)then
	temper(2:nspec) = temper(1)
	end if

      return
      end
