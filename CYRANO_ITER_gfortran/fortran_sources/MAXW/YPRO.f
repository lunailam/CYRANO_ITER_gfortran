      subroutine ypro(inscr)

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

      integer i
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
c             Electron-ion collisions (formula is for deuterium)
c             coulomb log = 15
              damp(1,ireg) = 1.55d9 * densit(1) / temper(1) ** 1.5
              end if
            end if
  1       continue
 
        else
          do 4 i = 1, nspec
          densit(i) =
     ;      ( n0(i,ireg) - nb(i,ireg) ) * y2 ** idexp(ireg)
     ;     +  nb(i,ireg)
cc   ;        n0(i,ireg) * y2 **idexp(ireg)
cc   ;     +  nb(i,ireg) * y3 **2
          omp(i) = dsqrt(
     ;        dfac / (amass(i)*mh*eps0)
     ;        * densit(i)) * dabs(zch(i)*eel)
 
            if(.not.coldpl(ireg) .or. sccol(ireg))then
            temper(i) =
     ;        ( t0(i,ireg)-tb(i,ireg) ) * y2 ** itexp(i,ireg)
     ;       +  tb(i,ireg)
cc   ;          t0(i,ireg) * y2 **itexp(i,ireg)
cc   ;       +  tb(i,ireg) * y3 **2
            vt(i) = dsqrt( 2.d0 *
     ;         tfac /(mh*amass(i)) * temper(i) )
              if(sccol(ireg))then
c             electron-ion collisions (formula is for deuterium)
c             coulomb log = 15
              damp(1,ireg) = 1.55d9 * densit(1) / temper(1) ** 1.5
              end if
            end if
 
  4       continue
          end if
 
      end if
 
      return
      end
